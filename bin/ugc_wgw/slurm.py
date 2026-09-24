"""SLURM job ids of an attempt, read from miniwdl-slurm's per-task submission logs, and orphan cancellation.

miniwdl-slurm 0.4.0 runs `sbatch --wait --parsable`; sbatch's stdout (the job id, `<id>` or `<id>;<cluster>`)
and slurmd's messages land in `<task dir>/slurm_singularity.log.txt`. Nothing else records the id.
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
from pathlib import Path
from typing import Iterator

LOG_NAME = "slurm_singularity.log.txt"
_ID_RE = re.compile(r"^\s*(\d+)(?:;|\s|$)")


def task_dirs(attempt_path: Path) -> Iterator[Path]:
    """Task directories holding a submission log, walking only call-*/failed* subtrees (never work/ or out/)."""
    for root, dirs, files in os.walk(attempt_path):
        dirs[:] = sorted(d for d in dirs if d.startswith("call-") or d.startswith("failed"))
        if LOG_NAME in files:
            yield Path(root)


def job_ids(attempt_path: Path) -> list[str]:
    ids: list[str] = []
    for d in task_dirs(Path(attempt_path)):
        try:
            with open(d / LOG_NAME, errors="replace") as fh:
                first = fh.readline()
        except OSError:
            continue
        m = _ID_RE.match(first)
        if m and m.group(1) not in ids:
            ids.append(m.group(1))
    return ids


def scancel(ids: list[str], *, timeout: float = 60.0) -> tuple[str, str]:
    """('ok' | 'skipped' | 'failed', detail)."""
    if not ids:
        return "skipped", "no job ids recorded"
    exe = shutil.which("scancel")
    if exe is None:
        return "skipped", "scancel not on PATH"
    try:
        res = subprocess.run([exe, *ids], capture_output=True, text=True, timeout=timeout, check=False)
    except (OSError, subprocess.SubprocessError) as exc:
        return "failed", str(exc)
    if res.returncode != 0:
        return "failed", (res.stderr or res.stdout).strip() or f"exit {res.returncode}"
    return "ok", f"cancelled {len(ids)} job(s)"
