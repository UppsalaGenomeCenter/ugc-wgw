"""Failure kinds: what a finished attempt's error documents and task logs say about why it failed.

Pure functions over a run directory (docs/DESIGN.md §8.1, decision log 2026-09-11). The kind
decides what the driver does next: `transient` failures are re-attempted with backoff, every
other kind blocks the subject until `ugc-wgw retry`.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

from . import slurm
from .config import Config
from .engine import RunResult
from .layout import run_files

KINDS = ("transient", "resource", "input", "tool", "version", "cancelled", "unknown")
BACKOFF_CAP_SECONDS = 3600.0
TAIL_BYTES = 64 * 1024
MESSAGE_CHARS = 200

TRANSIENT_CLASSES = frozenset({"driver_lost", "Interrupted", "Terminated", "launch_error", "killed", "NoResult"})
INPUT_CLASSES = frozenset({"InputError", "DownloadFailed"})
COMMAND_CLASSES = frozenset({"CommandFailed", "OutputError"})
RESOURCE_EXIT = frozenset({137, 253})  # 137 = killed by SIGKILL (cgroup OOM), 253 = sbatch --wait out-of-memory

RESOURCE_TEXT = re.compile(r"oom_kill|OUT_OF_MEMORY|Out of memory|DUE TO TIME LIMIT|TIMEOUT|DUE TO MEMORY"
                           # sbatch refused the request outright (node shape, partition limits): docs/guide/12-resources.md
                           r"|Requested node configuration is not available|CPU count per node can not be satisfied"
                           r"|Memory specification can not be satisfied|Requested time limit is invalid"
                           r"|exceeds partition limit")
TRANSIENT_TEXT = re.compile(r"DUE TO PREEMPTION|PREEMPTED|DUE TO NODE FAILURE|NODE_FAIL|Socket timed out"
                            r"|Unable to contact slurm controller|slurm_load_jobs error")
INPUT_TEXT = re.compile(r"No such file or directory|does not exist|Permission denied")


@dataclass
class Classification:
    kind: str
    error_class: str | None
    message: str
    evidence: dict[str, object] = field(default_factory=dict)


def _tail_lines(path: Path) -> list[str]:
    try:
        size = path.stat().st_size
        with open(path, "rb") as fh:
            if size > TAIL_BYTES:
                fh.seek(size - TAIL_BYTES)
            data = fh.read()
    except OSError:
        return []
    return data.decode("utf-8", errors="replace").splitlines()


def evidence_lines(task_dir: str | None, attempt_path: Path) -> list[tuple[str, str]]:
    """(source file, line) pairs to scan: the failed task's SLURM log and stderr, else miniwdl's stderr."""
    out: list[tuple[str, str]] = []
    if task_dir:
        d = Path(task_dir)
        if not d.is_absolute():
            d = attempt_path / d
        for name in (slurm.LOG_NAME, "stderr.txt"):
            p = d / name
            if p.exists():
                out.extend((str(p), line) for line in _tail_lines(p))
    if not out:
        p = run_files(attempt_path).stderr
        if p.exists():
            out.extend((str(p), line) for line in _tail_lines(p))
    return out


def _match(lines: list[tuple[str, str]], pattern: re.Pattern[str]) -> tuple[str, str] | None:
    for source, line in lines:
        if pattern.search(line):
            return source, line.strip()
    return None


def classify(result: RunResult, attempt_path: Path) -> Classification:
    """First matching rule wins; see the table in docs/guide/09-troubleshooting.md."""
    cls = result.error_class or None
    if cls == "Unknown":
        cls = None
    evidence: dict[str, object] = {"task_dir": result.task_dir, "node": result.node, "exit_status": result.exit_status}
    kind: str
    hit: tuple[str, str] | None = None
    if cls == "version_mismatch":
        kind = "version"
    elif cls in TRANSIENT_CLASSES:
        kind = "transient"
    elif cls in INPUT_CLASSES:
        kind = "input"
    elif result.exit_status in RESOURCE_EXIT or result.exit_code in RESOURCE_EXIT:
        kind = "resource"
    else:
        lines = evidence_lines(result.task_dir, Path(attempt_path))
        if (hit := _match(lines, RESOURCE_TEXT)) is not None:
            kind = "resource"
        elif (hit := _match(lines, TRANSIENT_TEXT)) is not None:
            kind = "transient"
        elif cls in COMMAND_CLASSES and (hit := _match(lines, INPUT_TEXT)) is not None:
            kind = "input"
        elif cls:
            kind = "tool"
        else:
            kind = "unknown"
    if hit is not None:
        evidence["source"], evidence["line"] = hit
    parts = [cls or "no class"]
    if result.exit_status is not None:
        parts.append(f"exit_status={result.exit_status}")
    if result.node:
        parts.append(f"node={result.node}")
    message = " ".join(parts)
    detail = hit[1] if hit is not None else result.error_message
    if detail:
        message += ": " + detail
    if len(message) > MESSAGE_CHARS:
        message = message[:MESSAGE_CHARS - 1] + "…"
    return Classification(kind, cls, message, evidence)


def backoff_delay(cfg: Config, attempt: int) -> float:
    """Seconds to wait before re-attempting after transient failure number `attempt` (doubling, capped)."""
    base = max(0.0, float(cfg.backoff_seconds))
    return min(base * (2 ** max(0, attempt - 1)), BACKOFF_CAP_SECONDS)
