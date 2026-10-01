"""SLURM accounting of one run attempt: `accounting.json` next to the manifest (guide chapter 07, DESIGN §9.2).

Every task of an attempt ran as one SLURM job whose id miniwdl-slurm left in the task directory (slurm.job_refs).
`sacct` is asked for those ids once per finished attempt (submit.finalize_run, or `ugc-wgw usage --collect`
later) and its rows are folded into one record per job: the allocation and the times from the parent line, the
peak memory and the disk volumes from the `.batch`/`.extern` step lines. Wall time from the workflow log includes
the queue wait; `Elapsed` here does not (End − Start), and `Start − Submit` is the wait itself.

Nothing in this module may fail a run: `capture()` reports through an event and returns.
"""
from __future__ import annotations

import datetime as dt
import os
import re
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Optional

from . import slurm
from .config import Config
from .db import RunRecord
from .layout import run_files
from .log import LOGGER, Events
from .util import UgcError, read_json, utc_now, write_json
from .wdllog import TaskStat, parse_workflow_log

ACCOUNTING_SCHEMA = 1
TERMINAL_STATES = frozenset({"COMPLETED", "FAILED", "CANCELLED", "TIMEOUT", "OUT_OF_MEMORY", "NODE_FAIL", "PREEMPTED",
                             "BOOT_FAIL", "DEADLINE", "REVOKED"})
_NO_VALUE = frozenset({"", "UNLIMITED", "Partition_Limit", "Unknown", "None", "INVALID", "N/A"})
_MEM_RE = re.compile(r"^\s*([0-9]*\.?[0-9]+)\s*([KMGTPkmgtp]?)([nc]?)\s*$")
_UNIT_MB = {"": 1.0, "K": 1 / 1024, "M": 1.0, "G": 1024.0, "T": 1024.0 ** 2, "P": 1024.0 ** 3}


# ---- value parsing ---------------------------------------------------------------------

def parse_seconds(text: object) -> Optional[float]:
    """sacct durations: `D-HH:MM:SS`, `HH:MM:SS`, `MM:SS`, `MM:SS.mmm`; None for empty or unlimited."""
    s = str(text or "").strip()
    if s in _NO_VALUE:
        return None
    days = 0
    if "-" in s:
        d, s = s.split("-", 1)
        if not d.isdigit():
            return None
        days = int(d)
    parts = s.split(":")
    try:
        nums = [float(p) for p in parts]
    except ValueError:
        return None
    if len(nums) == 3:
        h, m, sec = nums
    elif len(nums) == 2:
        h, (m, sec) = 0.0, nums
    elif len(nums) == 1:
        h, m, sec = 0.0, 0.0, nums[0]
    else:
        return None
    return days * 86400 + h * 3600 + m * 60 + sec


def parse_mem(text: object) -> tuple[Optional[float], bool]:
    """sacct memory under `--units=M`: (MB, per_cpu). `25600M`, `1.5G`, `1024K`; old releases append `n` (per node)
    or `c` (per CPU) to ReqMem. A bare number is taken as MB (the unit asked for)."""
    s = str(text or "").strip()
    if s in _NO_VALUE:
        return None, False
    m = _MEM_RE.match(s)
    if not m:
        return None, False
    value, unit, scope = float(m.group(1)), m.group(2).upper(), m.group(3)
    return value * _UNIT_MB[unit], scope == "c"


def parse_mb(text: object) -> Optional[float]:
    return parse_mem(text)[0]


def parse_int(text: object) -> Optional[int]:
    s = str(text or "").strip()
    try:
        return int(float(s))
    except ValueError:
        return None


def parse_tres(text: object) -> dict[str, str]:
    """`billing=48,cpu=48,gres/gpu=4,gres/gpu:a100=4,mem=192000M,node=1` -> dict."""
    out: dict[str, str] = {}
    for item in str(text or "").split(","):
        if "=" in item:
            k, v = item.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def gpu_of_tres(tres: dict[str, str]) -> tuple[int, Optional[str]]:
    """(GPU count, type): `gres/gpu=4` gives the count, `gres/gpu:a100=4` the type."""
    count = parse_int(tres.get("gres/gpu", "")) or 0
    gpu_type = None
    for k, v in tres.items():
        if k.startswith("gres/gpu:"):
            gpu_type = k.split(":", 1)[1]
            if not count:
                count = parse_int(v) or 0
    return count, gpu_type


def parse_time(text: object) -> Optional[dt.datetime]:
    """sacct's `YYYY-MM-DDTHH:MM:SS` (SLURM_TIME_FORMAT=standard, cluster-local, no zone)."""
    s = str(text or "").strip()
    if s in _NO_VALUE:
        return None
    try:
        return dt.datetime.fromisoformat(s)
    except ValueError:
        return None


# ---- folding sacct rows ------------------------------------------------------------------

def fold(rows: list[dict[str, str]]) -> dict[str, dict[str, object]]:
    """One record per base job id: `parent` (the job line's fields, or None when sacct returned only steps),
    `steps` (the step rows) and the folded maxima `max_rss_mb`, `disk_read_mb`, `disk_write_mb`."""
    out: dict[str, dict[str, object]] = {}
    for row in rows:
        jid = str(row.get("JobID", "")).strip()
        if not jid:
            continue
        base = jid.split(".", 1)[0]
        rec = out.setdefault(base, {"parent": None, "steps": [], "max_rss_mb": None, "disk_read_mb": None,
                                    "disk_write_mb": None})
        if "." in jid:
            rec["steps"].append(row)  # type: ignore[union-attr]
        else:
            rec["parent"] = row
        for key, col in (("max_rss_mb", "MaxRSS"), ("disk_read_mb", "MaxDiskRead"), ("disk_write_mb", "MaxDiskWrite")):
            v = parse_mb(row.get(col, ""))
            if v is not None and (rec[key] is None or v > rec[key]):  # type: ignore[operator]
                rec[key] = v
    return out


# ---- the record ---------------------------------------------------------------------------

@dataclass
class JobAcct:
    job_id: str
    call_path: str                   # `call-sub-1/call-b-0`; "" when the task directory is gone
    call_id: str
    retry: int
    task: str                        # inventory task name (from the workflow log), else the call's base name
    task_source: str                 # "workflow.log" | "call id"
    job_name: str
    state: str
    exit_code: str
    partition: str
    nodes: str
    submit: Optional[str]
    start: Optional[str]
    end: Optional[str]
    queue_seconds: Optional[float]
    elapsed_seconds: Optional[float]
    timelimit_minutes: Optional[int]
    req_cpus: Optional[int]
    alloc_cpus: Optional[int]
    cpu_seconds_alloc: Optional[float]   # AllocCPUS × Elapsed (CPUTimeRAW)
    cpu_seconds_used: Optional[float]    # TotalCPU
    req_mem_mb: Optional[float]
    alloc_mem_mb: Optional[float]
    max_rss_mb: Optional[float]
    gpus: int
    gpu_type: Optional[str]
    disk_read_mb: Optional[float]
    disk_write_mb: Optional[float]
    final: bool
    raw: dict[str, str] = field(default_factory=dict)


@dataclass
class Accounting:
    run_id: str
    attempt: int
    stage: str
    subject_type: str
    subject_id: str
    collected_at: str
    status: str                       # "ok" | "partial" (jobs sacct did not return, or not yet final)
    sacct_format: str
    jobs: list[JobAcct] = field(default_factory=list)
    missing: list[str] = field(default_factory=list)      # ids sacct returned no row for
    unmatched: list[str] = field(default_factory=list)    # ids whose task directory is gone (manifest only)
    notes: list[str] = field(default_factory=list)
    schema: int = ACCOUNTING_SCHEMA
    source: str = "sacct"

    def to_json(self) -> dict[str, object]:
        d = asdict(self)
        d["subject"] = {"type": d.pop("subject_type"), "id": d.pop("subject_id")}
        return d

    @classmethod
    def from_json(cls, doc: dict[str, object]) -> "Accounting":
        if doc.get("schema") != ACCOUNTING_SCHEMA:
            raise UgcError(f"unsupported accounting schema {doc.get('schema')!r}")
        subject = doc.get("subject") or {}
        jobs = [JobAcct(**j) for j in doc.get("jobs", [])]  # type: ignore[arg-type]
        return cls(run_id=str(doc.get("run_id", "")), attempt=int(doc.get("attempt", 0)), stage=str(doc.get("stage", "")),  # type: ignore[arg-type]
                   subject_type=str(subject.get("type", "")), subject_id=str(subject.get("id", "")),  # type: ignore[union-attr]
                   collected_at=str(doc.get("collected_at", "")), status=str(doc.get("status", "")),
                   sacct_format=str(doc.get("sacct_format", "")), jobs=jobs, missing=list(doc.get("missing", [])),  # type: ignore[arg-type]
                   unmatched=list(doc.get("unmatched", [])), notes=list(doc.get("notes", [])),  # type: ignore[arg-type]
                   source=str(doc.get("source", "sacct")))


def _job(ref: Optional[slurm.JobRef], job_id: str, rec: dict[str, object], task: Optional[TaskStat]) -> JobAcct:
    parent = rec.get("parent")
    steps: list[dict[str, str]] = rec.get("steps") or []  # type: ignore[assignment]
    row: dict[str, str] = parent if isinstance(parent, dict) else (steps[0] if steps else {})
    call_path = ref.call_path if ref else ""
    call_id = ref.call_id if ref else str(row.get("JobName", ""))
    if task is not None:
        name, src = task.name, "workflow.log"
    else:
        name, src = slurm.call_base(call_id) if call_id else "", "call id"
    submit, start, end = parse_time(row.get("Submit")), parse_time(row.get("Start")), parse_time(row.get("End"))
    elapsed = parse_seconds(row.get("ElapsedRaw"))
    if elapsed is None and start and end:
        elapsed = (end - start).total_seconds()
    alloc_cpus = parse_int(row.get("AllocCPUS"))
    cpu_alloc = parse_seconds(row.get("CPUTimeRAW"))
    if cpu_alloc is None and alloc_cpus is not None and elapsed is not None:
        cpu_alloc = alloc_cpus * elapsed
    tres = parse_tres(row.get("AllocTRES"))
    gpus, gpu_type = gpu_of_tres(tres)
    req_mem, per_cpu = parse_mem(row.get("ReqMem"))
    if req_mem is not None and per_cpu:
        req_mem *= (parse_int(row.get("ReqCPUS")) or alloc_cpus or 1)
    alloc_mem = parse_mb(tres.get("mem", ""))
    state = str(row.get("State", "")).split()[0] if str(row.get("State", "")).strip() else ""
    return JobAcct(
        job_id=job_id, call_path=call_path, call_id=call_id, retry=ref.retry if ref else 0, task=name, task_source=src,
        job_name=str(row.get("JobName", "")), state=state, exit_code=str(row.get("ExitCode", "")),
        partition=str(row.get("Partition", "")), nodes=str(row.get("NodeList", "")),
        submit=submit.isoformat() if submit else None, start=start.isoformat() if start else None,
        end=end.isoformat() if end else None,
        queue_seconds=(start - submit).total_seconds() if submit and start else None,
        elapsed_seconds=elapsed, timelimit_minutes=parse_int(row.get("TimelimitRaw")),
        req_cpus=parse_int(row.get("ReqCPUS")), alloc_cpus=alloc_cpus, cpu_seconds_alloc=cpu_alloc,
        cpu_seconds_used=parse_seconds(row.get("TotalCPU")), req_mem_mb=req_mem,
        alloc_mem_mb=alloc_mem if alloc_mem is not None else req_mem,
        max_rss_mb=rec.get("max_rss_mb"), gpus=gpus, gpu_type=gpu_type,  # type: ignore[arg-type]
        disk_read_mb=rec.get("disk_read_mb"), disk_write_mb=rec.get("disk_write_mb"),  # type: ignore[arg-type]
        final=state in TERMINAL_STATES, raw={k: v for k, v in row.items() if v != ""})


def build(run: RunRecord, refs: list[slurm.JobRef], rows: list[dict[str, str]], tasks: list[TaskStat],
          manifest_ids: list[str], *, sacct_format: str = "", collected_at: Optional[str] = None) -> Accounting:
    """Join the task directories' job ids with sacct's rows and the workflow log's task names."""
    folded = fold(rows)
    attempt = run.run_path
    by_dir: dict[str, TaskStat] = {}
    by_path: dict[str, TaskStat] = {}
    for t in tasks:
        if t.dir:
            try:
                by_dir[os.path.relpath(t.dir, attempt)] = t
            except ValueError:
                pass
        if t.call_path:
            by_path.setdefault(t.call_path, t)
    acct = Accounting(run_id=run.run_id, attempt=run.attempt, stage=run.stage, subject_type=run.subject_type,
                      subject_id=run.subject_id, collected_at=collected_at or utc_now(), status="ok",
                      sacct_format=sacct_format)
    seen: set[str] = set()
    order: list[tuple[str, Optional[slurm.JobRef]]] = [(r.job_id, r) for r in refs]
    order += [(i, None) for i in manifest_ids if i not in {r.job_id for r in refs}]
    for job_id, ref in order:
        if job_id in seen:
            continue
        seen.add(job_id)
        rec = folded.get(job_id)
        if rec is None:
            acct.missing.append(job_id)
            continue
        if ref is None:
            acct.unmatched.append(job_id)
        task = (by_dir.get(ref.call_path) or by_path.get(ref.call_path)) if ref else None
        acct.jobs.append(_job(ref, job_id, rec, task))
    not_final = [j.job_id for j in acct.jobs if not j.final]
    if acct.missing:
        acct.notes.append(f"sacct returned no row for {len(acct.missing)} job(s)")
    if not_final:
        acct.notes.append(f"{len(not_final)} job(s) not in a final state when collected")
    if acct.unmatched:
        acct.notes.append(f"{len(acct.unmatched)} job(s) without a task directory (from the manifest only)")
    if any(j.task_source == "call id" for j in acct.jobs):
        acct.notes.append("some task names come from the call id (no task setup line in the workflow log)")
    acct.status = "ok" if not acct.missing and not not_final else "partial"
    return acct


def read(attempt_path: Path) -> Optional[Accounting]:
    path = run_files(attempt_path).accounting
    if not path.exists():
        return None
    try:
        doc = read_json(path)
    except (OSError, ValueError):
        return None
    if not isinstance(doc, dict):
        return None
    try:
        return Accounting.from_json(doc)
    except (UgcError, TypeError, KeyError, ValueError):
        return None


def write(attempt_path: Path, acct: Accounting) -> Path:
    path = run_files(attempt_path).accounting
    write_json(path, acct.to_json())
    return path


def collect(cfg: Config, run: RunRecord) -> tuple[str, Optional[Accounting], str]:
    """('ok' | 'partial' | 'skipped' | 'failed', the record written, detail). Skipped without job ids or sacct."""
    refs = slurm.job_refs(run.run_path)
    manifest_ids = [str(i) for i in (run.meta.get("slurm_job_ids") or [])]
    ids: list[str] = []
    for i in [r.job_id for r in refs] + manifest_ids:
        if i not in ids:
            ids.append(i)
    if not ids:
        return "skipped", None, "no job ids recorded"
    status, rows, detail = slurm.sacct(ids, timeout=cfg.accounting_timeout)
    if status != "ok":
        return status, None, detail
    fmt = detail.split("fields ", 1)[1] if "fields " in detail else ""
    tasks, _ = parse_workflow_log(run_files(run.run_path).workflow_log_json)
    acct = build(run, refs, rows, tasks, manifest_ids, sacct_format=fmt)
    write(run.run_path, acct)
    return acct.status, acct, f"{len(acct.jobs)} job(s), {len(acct.missing)} missing"


def capture(cfg: Config, events: Events, run: RunRecord) -> str:
    """collect() for a finished run, reported as a `run.accounting` event when sacct was consulted (nothing is
    recorded for a run without job ids or a machine without sacct); never raises."""
    t0 = time.monotonic()
    try:
        status, acct, detail = collect(cfg, run)
    except Exception as exc:  # noqa: BLE001 - accounting must never take a run down
        LOGGER.warning("accounting of run %s failed: %s", run.run_id, exc)
        status, acct, detail = "failed", None, str(exc)
    if status == "skipped":
        LOGGER.info("accounting of run %s skipped: %s", run.run_id, detail)
        return status
    events.emit("run.accounting", run_id=run.run_id, level="warning" if status == "failed" else "info",
                subject=run.subject_id, stage=run.stage, attempt=run.attempt, status=status,
                jobs=len(acct.jobs) if acct else 0, missing=len(acct.missing) if acct else 0,
                seconds=round(time.monotonic() - t0, 1), detail=detail)
    return status
