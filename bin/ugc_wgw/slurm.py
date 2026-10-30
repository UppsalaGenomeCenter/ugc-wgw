"""SLURM job ids of an attempt, read from miniwdl-slurm's per-task submission logs; `sacct` for their accounting;
orphan cancellation.

miniwdl-slurm 0.4.0 runs `sbatch --wait --parsable --job-name <call id>` in the task directory; sbatch's stdout
(the job id, `<id>` or `<id>;<cluster>`) lands in `<task dir>/slurm_singularity.log.txt`, which miniwdl truncates
for every try, and SLURM writes each job's own output to `<task dir>/slurm-<id>.out`, which stays. A retried task
(miniwdl 1.15) keeps every try in the same directory: `work`, `stdout.txt`, `stderr.txt` for the first,
`work2`, `stdout2.txt`, `stderr2.txt` for the second; older layouts put earlier tries under `call-x/failedN/`.
Job ids therefore come from the sbatch log and from the `slurm-<id>.out` names, so a retried task's earlier jobs
reach the accounting too (the 864-sample run's first cohort_call had 52 such jobs, 2026-10-30).
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
    rel_dir: str      # task dir relative to the attempt, as on disk: `call-sub-1/call-b-0` (`call-a/failed1` in the older layout)
    call_path: str    # rel_dir without a trailing failedN: unique per call within the attempt
    call_id: str      # the last component of call_path, miniwdl's run id and sbatch's --job-name
    retry: int        # 0 for the call's final job, n for an earlier one (an earlier `slurm-<id>.out` of the task dir, or failedN)


_OUT_RE = re.compile(r"^slurm-(\d+)\.out$")


def task_dirs(attempt_path: Path) -> Iterator[Path]:
    """Task directories holding a submission log or a job output file, walking only call-*/failed* subtrees
    (never work/ or out/)."""
    for root, dirs, files in os.walk(attempt_path):
        dirs[:] = sorted(d for d in dirs if d.startswith("call-") or d.startswith("failed"))
        if LOG_NAME in files or any(_OUT_RE.match(f) for f in files):
            yield Path(root)


def job_refs(attempt_path: Path) -> list[JobRef]:
    attempt_path = Path(attempt_path)
    refs: list[JobRef] = []
    for d in task_dirs(attempt_path):
        ids: list[str] = []
        try:
            with open(d / LOG_NAME, errors="replace") as fh:
                m = _ID_RE.match(fh.readline())
            if m:
                ids.append(m.group(1))
        except OSError:
            pass
        final = ids[0] if ids else None
        # the job output files name every job the task submitted; the sbatch log only the last try's
        for f in sorted(os.listdir(d), key=lambda n: int(_OUT_RE.match(n).group(1)) if _OUT_RE.match(n) else -1):
            om = _OUT_RE.match(f)
            if om and om.group(1) not in ids:
                ids.append(om.group(1))
        if not ids:
            continue
        rel = d.relative_to(attempt_path).as_posix()
        parts = rel.split("/")
        retry = 0
        fm = _FAILED_RE.match(parts[-1]) if parts else None
        if fm and len(parts) > 1:
            retry = int(fm.group(1) or 1)
            parts = parts[:-1]
        # earlier tries of the same directory get retry numbers below the final try's 0: 1 for the one before it
        earlier = [i for i in ids if i != final]   # ascending job ids: the oldest try first
        for idx, job_id in enumerate(earlier):
            refs.append(JobRef(job_id, rel, "/".join(parts), parts[-1], retry + len(earlier) - idx))
        if final is not None:
            refs.append(JobRef(final, rel, "/".join(parts), parts[-1], retry))
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
