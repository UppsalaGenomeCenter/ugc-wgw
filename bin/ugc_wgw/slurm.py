"""SLURM job ids of an attempt, read from miniwdl-slurm's per-task submission logs; `sacct` for their accounting;
orphan cancellation.

miniwdl-slurm 0.4.0 runs `sbatch --wait --parsable --job-name <call id>`; sbatch's stdout (the job id, `<id>` or
`<id>;<cluster>`) and slurmd's messages land in `<task dir>/slurm_singularity.log.txt`. Nothing else records the id.
A retried task keeps its final attempt in `call-x/` and the earlier ones in `call-x/failedN/`.
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

LOG_NAME = "slurm_singularity.log.txt"
_ID_RE = re.compile(r"^\s*(\d+)(?:;|\s|$)")
_FAILED_RE = re.compile(r"^failed(\d*)$")
_SHARD_RE = re.compile(r"(-\d+(-[^-]+)?)+$")   # miniwdl's scatter suffixes: `-3`, `-03-chr3`; task names carry no hyphen

# What `sacct` is asked for, per job: the parent line carries the allocation and the times, the `.batch`/`.extern`
# step lines the peak memory and the disk volumes (accounting.fold). SACCT_OPTIONAL are dropped on an older sacct
# that rejects them.
SACCT_FIELDS = ("JobID", "JobName%100", "State", "ExitCode", "Submit", "Start", "End", "ElapsedRaw", "ReqCPUS",
                "AllocCPUS", "CPUTimeRAW", "TotalCPU", "ReqMem", "MaxRSS", "AllocTRES", "Partition", "NodeList",
                "TimelimitRaw", "MaxDiskRead", "MaxDiskWrite")
SACCT_OPTIONAL = frozenset({"ReqCPUS", "TimelimitRaw", "MaxDiskRead", "MaxDiskWrite"})


@dataclass(frozen=True)
class JobRef:
    """A job id and the task directory it was read from."""
    job_id: str
    rel_dir: str      # task dir relative to the attempt, as on disk: `call-sub-1/call-b-0`, `call-a/failed1`
    call_path: str    # rel_dir without a trailing failedN: unique per call within the attempt
    call_id: str      # the last component of call_path, miniwdl's run id and sbatch's --job-name
    retry: int        # 0 for the final attempt of the call, n for failedN


def task_dirs(attempt_path: Path) -> Iterator[Path]:
    """Task directories holding a submission log, walking only call-*/failed* subtrees (never work/ or out/)."""
    for root, dirs, files in os.walk(attempt_path):
        dirs[:] = sorted(d for d in dirs if d.startswith("call-") or d.startswith("failed"))
        if LOG_NAME in files:
            yield Path(root)


def job_refs(attempt_path: Path) -> list[JobRef]:
    attempt_path = Path(attempt_path)
    refs: list[JobRef] = []
    for d in task_dirs(attempt_path):
        try:
            with open(d / LOG_NAME, errors="replace") as fh:
                first = fh.readline()
        except OSError:
            continue
        m = _ID_RE.match(first)
        if not m:
            continue
        rel = d.relative_to(attempt_path).as_posix()
        parts = rel.split("/")
        retry = 0
        fm = _FAILED_RE.match(parts[-1]) if parts else None
        if fm and len(parts) > 1:
            retry = int(fm.group(1) or 1)
            parts = parts[:-1]
        refs.append(JobRef(m.group(1), rel, "/".join(parts), parts[-1], retry))
    return refs


def job_ids(attempt_path: Path) -> list[str]:
    ids: list[str] = []
    for ref in job_refs(attempt_path):
        if ref.job_id not in ids:
            ids.append(ref.job_id)
    return ids


def call_base(call_id: str) -> str:
    """`call-pbmm2_align_wgs-3` -> `pbmm2_align_wgs`, `call-deepvariant-03-chr3` -> `deepvariant`: the call name
    without `call-` and miniwdl's scatter suffixes. An aliased
    call (`call x as y`) gives the alias, not the task; the workflow log names the task."""
    name = call_id[len("call-"):] if call_id.startswith("call-") else call_id
    return _SHARD_RE.sub("", name)


def sacct(ids: list[str], *, timeout: float = 120.0, chunk: int = 400) -> tuple[str, list[dict[str, str]], str]:
    """('ok' | 'skipped' | 'failed', rows keyed by field name, detail). One `sacct -j` per chunk of ids; a sacct
    that rejects one of SACCT_OPTIONAL is asked again without them."""
    if not ids:
        return "skipped", [], "no job ids"
    exe = shutil.which("sacct")
    if exe is None:
        return "skipped", [], "sacct not on PATH"
    env = dict(os.environ, SLURM_TIME_FORMAT="standard")
    fields = list(SACCT_FIELDS)
    rows: list[dict[str, str]] = []
    i = 0
    while i < len(ids):
        part = ids[i:i + chunk]
        cmd = [exe, "-j", ",".join(part), "-P", "--noheader", "--units=M", "--format=" + ",".join(fields)]
        try:
            res = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, check=False, env=env)
        except (OSError, subprocess.SubprocessError) as exc:
            return "failed", rows, str(exc)
        if res.returncode != 0:
            err = (res.stderr or res.stdout).strip()
            optional = [f for f in fields if f.split("%", 1)[0] in SACCT_OPTIONAL]
            if "Invalid field" in err and optional:
                fields = [f for f in fields if f.split("%", 1)[0] not in SACCT_OPTIONAL]
                continue   # same chunk, fewer fields
            return "failed", rows, err.splitlines()[-1] if err else f"sacct exited {res.returncode}"
        names = [f.split("%", 1)[0] for f in fields]
        for line in res.stdout.splitlines():
            if not line.strip():
                continue
            vals = line.split("|")
            if len(vals) < len(names):
                vals += [""] * (len(names) - len(vals))
            rows.append(dict(zip(names, vals)))
        i += chunk
    names = [f.split("%", 1)[0] for f in fields]
    return "ok", rows, f"{len(rows)} rows for {len(ids)} jobs; fields {','.join(names)}"


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
