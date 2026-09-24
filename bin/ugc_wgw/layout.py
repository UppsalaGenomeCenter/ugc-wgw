"""The results layout (docs/DESIGN.md §9.1, attempt directories per §18 2026-09-10). Nothing else knows these paths."""
from __future__ import annotations

import os
import re
from dataclasses import dataclass
from pathlib import Path

from .util import UgcError

_ATTEMPT_RE = re.compile(r"^attempt-(\d+)$")


def subject_root(results_dir: Path, subject_type: str, subject_id: str) -> Path:
    sub = "samples" if subject_type == "sample" else "cohorts"
    return Path(results_dir) / sub / subject_id


def stage_dir(results_dir: Path, subject_type: str, subject_id: str, ugc_wgw_version: str, stage: str) -> Path:
    return subject_root(results_dir, subject_type, subject_id) / ugc_wgw_version / stage


def attempt_dir(results_dir: Path, subject_type: str, subject_id: str, ugc_wgw_version: str, stage: str, attempt: int) -> Path:
    return stage_dir(results_dir, subject_type, subject_id, ugc_wgw_version, stage) / f"attempt-{attempt}"


def existing_attempts(stage_path: Path) -> list[int]:
    if not stage_path.is_dir():
        return []
    out = []
    for child in stage_path.iterdir():
        m = _ATTEMPT_RE.match(child.name)
        if m and child.is_dir():
            out.append(int(m.group(1)))
    return sorted(out)


def make_attempt_dir(path: Path) -> Path:
    """Create the attempt directory; it must be empty (miniwdl refuses a directory holding an old out/)."""
    path.mkdir(parents=True, exist_ok=True)
    if any(path.iterdir()):
        raise UgcError(f"attempt directory is not empty: {path}")
    return path


def point_current(stage_path: Path, attempt: int) -> None:
    """Atomically (re)point `<stage>/current` at `attempt-<n>` with a relative symlink."""
    tmp = stage_path / ".current.tmp"
    if tmp.is_symlink() or tmp.exists():
        tmp.unlink()
    os.symlink(f"attempt-{attempt}", tmp)
    os.replace(tmp, stage_path / "current")


@dataclass(frozen=True)
class RunFiles:
    inputs: Path
    run_json: Path
    outputs: Path
    error: Path
    workflow_log: Path
    workflow_log_json: Path
    stdout: Path
    stderr: Path
    manifest: Path


def run_files(attempt_path: Path) -> RunFiles:
    p = Path(attempt_path)
    return RunFiles(
        inputs=p / "inputs.json",
        run_json=p / "run.json",
        outputs=p / "outputs.json",
        error=p / "error.json",
        workflow_log=p / "workflow.log",
        workflow_log_json=p / "workflow.log.json",
        stdout=p / "miniwdl.stdout",
        stderr=p / "miniwdl.stderr",
        manifest=p / "run_manifest.json",
    )
