"""`ugc-wgw usage` and the report's "Resource usage" section: core-hours, GPU-hours, memory, queue wait, disk and a
cost estimate per stage, task and subject, with a per-sample figure (guide chapter 07).

Each attempt contributes either its accounting.json (sacct: compute time without the queue wait, allocated and
requested CPUs, peak memory, GPUs, disk volumes) or, without one, an estimate: the workflow log's wall time per
task (which includes the queue wait) times the cpu the task was launched with (policy row, miniwdl's cap, else the
inventory's effective value). Prices are the user's (`config.json` `prices`, `--price`); the raw figures are always
shown, so a cost is a transparent multiplication.
"""
from __future__ import annotations

import os
import shutil
import stat
import statistics
from collections import Counter, defaultdict
from dataclasses import dataclass, field, fields
from pathlib import Path
from typing import Optional

from . import accounting, resources
from .config import Config
from .db import DB, RunRecord
from .layout import run_files
from .log import Events
from .resources import GPU_TASKS, Effective, GpuChoice
from .stages import mode_stages
from .util import UgcError
from .wdllog import TaskStat, parse_workflow_log

PRICE_KEYS = ("cpu_hour", "gpu_hour", "mem_gb_hour", "storage_gb_month")
BASES = ("allocated", "requested")
FINISHED = ("success", "failed", "cancelled")


# ---- prices ------------------------------------------------------------------------------

@dataclass
class Prices:
    cpu_hour: Optional[float] = None
    gpu_hour: Optional[float] = None
    mem_gb_hour: Optional[float] = None
    storage_gb_month: Optional[float] = None
    currency: str = ""

    def any(self) -> bool:
        return any(getattr(self, k) is not None for k in PRICE_KEYS)

    def describe(self) -> str:
        parts = [f"{k} {getattr(self, k):g}" for k in PRICE_KEYS if getattr(self, k) is not None]
        text = ", ".join(parts) if parts else "none"
        return f"{text} {self.currency}".rstrip() if parts and self.currency else text

    def money(self, x: Optional[float]) -> str:
        if x is None:
            return "–"
        return f"{x:,.2f} {self.currency}".rstrip()


def parse_prices(cfg: Config, overrides: Optional[list[str]] = None) -> Prices:
    """config.json `prices`, then `--price KEY=VALUE`; keys cpu_hour, gpu_hour, mem_gb_hour, storage_gb_month, currency."""
    p = Prices()
    known = ", ".join((*PRICE_KEYS, "currency"))
    for source, items in (("config.json prices", list(cfg.prices.items())),
                          ("--price", [tuple(o.split("=", 1)) if "=" in o else (o, "") for o in (overrides or [])])):
        for key, value in items:
            key = str(key).strip()
            if key == "currency":
                p.currency = str(value).strip()
                continue
            if key not in PRICE_KEYS:
                raise UgcError(f"{source}: unknown price {key!r} (known: {known})")
            try:
                v = float(value)  # type: ignore[arg-type]
            except (TypeError, ValueError):
                raise UgcError(f"{source}: {key} must be a number, not {value!r}") from None
            if v < 0:
                raise UgcError(f"{source}: {key} must not be negative")
            setattr(p, key, v)
    return p


# ---- the additive figures ---------------------------------------------------------------------

@dataclass
class Usage:
    jobs: int = 0
    failed_jobs: int = 0             # jobs whose state is not COMPLETED (sacct), or failed tasks (estimate)
    attempts: int = 0
    estimated_attempts: int = 0
    elapsed_s: float = 0.0           # sum of the jobs' elapsed (compute) time
    queue_s: list[float] = field(default_factory=list)
    core_h_alloc: float = 0.0        # AllocCPUS × elapsed
    core_h_req: float = 0.0          # ReqCPUS × elapsed (what a cloud would be provisioned with)
    cpu_h_used: float = 0.0          # TotalCPU
    core_h_with_used: float = 0.0    # allocated core-hours of the jobs that reported TotalCPU (efficiency denominator)
    gpu_h: float = 0.0
    mem_gb_h: float = 0.0            # allocated memory × elapsed
    peak_rss_mb: Optional[float] = None
    req_mem_mb_max: Optional[float] = None
    disk_read_mb: float = 0.0
    disk_write_mb: float = 0.0
    wasted_core_h: float = 0.0       # allocated core-hours of jobs that did not complete
    out_bytes: Optional[int] = None
    attempt_bytes: Optional[int] = None

    _SUMS = ("jobs", "failed_jobs", "attempts", "estimated_attempts", "elapsed_s", "core_h_alloc", "core_h_req",
             "cpu_h_used", "core_h_with_used", "gpu_h", "mem_gb_h", "disk_read_mb", "disk_write_mb", "wasted_core_h")

    def add(self, o: "Usage") -> None:
        for k in self._SUMS:
            setattr(self, k, getattr(self, k) + getattr(o, k))
        self.queue_s += o.queue_s
        for k in ("peak_rss_mb", "req_mem_mb_max"):
            a, b = getattr(self, k), getattr(o, k)
            setattr(self, k, b if a is None else a if b is None else max(a, b))
        for k in ("out_bytes", "attempt_bytes"):
            a, b = getattr(self, k), getattr(o, k)
            setattr(self, k, b if a is None else a if b is None else a + b)

    def scaled(self, f: float) -> "Usage":
        """The additive figures times f (a cohort's share per member, a mean per sample); counts, peaks and the
        queue list are not meaningful per sample and are left empty."""
        u = Usage()
        for k in ("elapsed_s", "core_h_alloc", "core_h_req", "cpu_h_used", "core_h_with_used", "gpu_h", "mem_gb_h",
                  "disk_read_mb", "disk_write_mb", "wasted_core_h"):
            setattr(u, k, getattr(self, k) * f)
        for k in ("out_bytes", "attempt_bytes"):
            v = getattr(self, k)
            setattr(u, k, None if v is None else int(round(v * f)))
        return u

    @property
    def efficiency(self) -> Optional[float]:
        return 100.0 * self.cpu_h_used / self.core_h_with_used if self.core_h_with_used > 0 else None

    def core_h(self, basis: str) -> float:
        return self.core_h_alloc if basis == "allocated" else self.core_h_req

    @property
    def queue_median(self) -> Optional[float]:
        return statistics.median(self.queue_s) if self.queue_s else None

    @property
    def queue_max(self) -> Optional[float]:
        return max(self.queue_s) if self.queue_s else None

    @property
    def mem_use(self) -> Optional[float]:
        if self.peak_rss_mb is None or not self.req_mem_mb_max:
            return None
        return 100.0 * self.peak_rss_mb / self.req_mem_mb_max

    def cost(self, prices: Prices, basis: str) -> Optional[float]:
        if not prices.any():
            return None
        total = 0.0
        if prices.cpu_hour is not None:
            total += self.core_h(basis) * prices.cpu_hour
        if prices.gpu_hour is not None:
            total += self.gpu_h * prices.gpu_hour
        if prices.mem_gb_hour is not None:
            total += self.mem_gb_h * prices.mem_gb_hour
        return total

    def storage_per_month(self, prices: Prices) -> Optional[float]:
        if prices.storage_gb_month is None or self.out_bytes is None:
            return None
        return self.out_bytes / 1024 ** 3 * prices.storage_gb_month

    def to_dict(self) -> dict[str, object]:
        d: dict[str, object] = {}
        for f in fields(self):
            v = getattr(self, f.name)
            if f.name == "queue_s":
                d["queue_median_s"], d["queue_max_s"] = self.queue_median, self.queue_max
                continue
            d[f.name] = round(v, 4) if isinstance(v, float) else v
        d["efficiency_pct"] = None if self.efficiency is None else round(self.efficiency, 1)
        return d


# ---- per attempt ----------------------------------------------------------------------------

@dataclass
class AttemptUsage:
    run: RunRecord
    source: str                      # "sacct" | "estimate" | "none"
    status: str                      # accounting.json status (ok | partial), "" otherwise
    usage: Usage
    per_task: dict[str, Usage]
    members: list[str] = field(default_factory=list)   # a cohort subject's members


def from_accounting(acct: accounting.Accounting) -> tuple[Usage, dict[str, Usage]]:
    total, per_task = Usage(), defaultdict(Usage)
    for j in acct.jobs:
        u = Usage(jobs=1, failed_jobs=0 if j.state == "COMPLETED" else 1)
        elapsed = j.elapsed_seconds or 0.0
        hours = elapsed / 3600.0
        u.elapsed_s = elapsed
        if j.queue_seconds is not None:
            u.queue_s = [max(0.0, j.queue_seconds)]
        alloc = j.alloc_cpus or 0
        u.core_h_alloc = (j.cpu_seconds_alloc / 3600.0) if j.cpu_seconds_alloc is not None else alloc * hours
        u.core_h_req = (j.req_cpus if j.req_cpus is not None else alloc) * hours
        if j.cpu_seconds_used is not None:
            u.cpu_h_used = j.cpu_seconds_used / 3600.0
            u.core_h_with_used = u.core_h_alloc
        u.gpu_h = j.gpus * hours
        mem = j.alloc_mem_mb if j.alloc_mem_mb is not None else j.req_mem_mb
        if mem is not None:
            u.mem_gb_h = mem / 1024.0 * hours
        u.peak_rss_mb, u.req_mem_mb_max = j.max_rss_mb, j.req_mem_mb
        u.disk_read_mb, u.disk_write_mb = j.disk_read_mb or 0.0, j.disk_write_mb or 0.0
        if j.state != "COMPLETED":
            u.wasted_core_h = u.core_h_alloc
        total.add(u)
        per_task[j.task or j.call_id].add(u)
    return total, dict(per_task)


def _gres_count(gres: Optional[str]) -> Optional[int]:
    """`gpu:a100:4` / `gpu:4` -> 4."""
    if not gres:
        return None
    try:
        return int(gres.rsplit(":", 1)[-1])
    except ValueError:
        return None


def estimate(tasks: list[TaskStat], inv: dict[str, Effective], gpu: GpuChoice) -> tuple[Usage, dict[str, Usage]]:
    """Without accounting: each finished, uncached task's wall time (queue wait included) times the cpu it was
    launched with, and the inventory's memory; GPUs from the plugin's NOTICE, else the project's GPU choice."""
    total, per_task = Usage(), defaultdict(Usage)
    for t in tasks:
        if t.duration is None or t.cached:
            continue
        eff = inv.get(t.name)
        hours = t.duration / 3600.0
        cpu = t.cpu_launched if t.cpu_launched is not None else (eff.cpu if eff else None)
        u = Usage(jobs=1, failed_jobs=1 if t.failed else 0, elapsed_s=t.duration)
        if cpu is not None:
            u.core_h_alloc = u.core_h_req = cpu * hours
        gpus = _gres_count(t.gres)
        if gpus is None and t.name in GPU_TASKS:
            gpus = gpu.count(t.name)
        u.gpu_h = (gpus or 0) * hours
        if eff and eff.memory:
            u.mem_gb_h = eff.memory / 1024 ** 3 * hours
            u.req_mem_mb_max = eff.memory / 1024 ** 2
        if t.failed:
            u.wasted_core_h = u.core_h_alloc
        total.add(u)
        per_task[t.name].add(u)
    return total, dict(per_task)


def dir_bytes(path: Path) -> int:
    """Bytes under a directory, counting a hardlinked file once (miniwdl's `output_hardlinks` shares out/ and work/)."""
    total, seen = 0, set()
    for root, _dirs, files in os.walk(path):
        for f in files:
            try:
                st = os.lstat(os.path.join(root, f))
            except OSError:
                continue
            if stat.S_ISREG(st.st_mode) and st.st_nlink > 1:
                key = (st.st_dev, st.st_ino)
                if key in seen:
                    continue
                seen.add(key)
            total += st.st_size
    return total


def inventory(cfg: Config) -> tuple[dict[str, Effective], GpuChoice, list[str]]:
    """The effective requests per task (resources.load), or an empty inventory with a note when it cannot be read."""
    try:
        res = resources.load(cfg)
    except UgcError as exc:
        return {}, GpuChoice(cfg.deepvariant, cfg.gpu_type, cfg.parabricks_gpus), [f"inventory unavailable, estimates carry no cpu or memory: {exc}"]
    return {e.task: e for e in res.rows()}, res.gpu, []


def attempt_usage(cfg: Config, run: RunRecord, *, inv: dict[str, Effective], gpu: GpuChoice, sizes: bool,
                  members: Optional[list[str]] = None) -> AttemptUsage:
    acct = accounting.read(run.run_path)
    if acct is not None and acct.jobs:
        usage, per_task = from_accounting(acct)
        source, status = "sacct", acct.status
    else:
        tasks, _ = parse_workflow_log(run_files(run.run_path).workflow_log_json)
        if tasks:
            usage, per_task = estimate(tasks, inv, gpu)
            usage.estimated_attempts = 1
            source, status = "estimate", ""
        else:
            usage, per_task, source, status = Usage(), {}, "none", ""
    usage.attempts = 1
    if sizes:
        out = run.run_path / "out"
        usage.out_bytes = dir_bytes(out) if out.exists() else 0
        usage.attempt_bytes = dir_bytes(run.run_path) if run.run_path.exists() else 0
    return AttemptUsage(run, source, status, usage, per_task, list(members or []))


# ---- aggregation -------------------------------------------------------------------------------

@dataclass
class Summary:
    totals: Usage
    by_stage: dict[str, Usage]
    by_task: dict[str, Usage]
    by_subject: dict[tuple[str, str], Usage]
    clean: Usage                          # the latest successful attempt of every (subject, stage)
    n_samples: int                        # samples with at least one successful run in the selection
    sources: Counter
    notes: list[str] = field(default_factory=list)

    @property
    def per_sample_clean(self) -> Optional[Usage]:
        return self.clean.scaled(1.0 / self.n_samples) if self.n_samples else None

    @property
    def per_sample_as_run(self) -> Optional[Usage]:
        return self.totals.scaled(1.0 / self.n_samples) if self.n_samples else None


def aggregate(attempts: list[AttemptUsage]) -> Summary:
    totals, clean = Usage(), Usage()
    by_stage: dict[str, Usage] = defaultdict(Usage)
    by_task: dict[str, Usage] = defaultdict(Usage)
    by_subject: dict[tuple[str, str], Usage] = defaultdict(Usage)
    latest: dict[tuple[str, str, str], AttemptUsage] = {}
    for a in attempts:
        totals.add(a.usage)
        by_stage[a.run.stage].add(a.usage)
        by_subject[(a.run.subject_type, a.run.subject_id)].add(a.usage)
        for task, u in a.per_task.items():
            by_task[task].add(u)
        if a.run.status == "success":
            key = (a.run.subject_type, a.run.subject_id, a.run.stage)
            if key not in latest or a.run.attempt > latest[key].run.attempt:
                latest[key] = a
    for a in latest.values():
        clean.add(a.usage)
    samples = {a.run.subject_id for a in latest.values() if a.run.subject_type == "sample"}
    sources = Counter(a.source for a in attempts)
    notes = []
    if sources["estimate"]:
        notes.append(f"{sources['estimate']} attempt(s) estimated from workflow.log wall time (queue wait included) × "
                     "requested cpu; no sacct accounting")
    if sources["none"]:
        notes.append(f"{sources['none']} attempt(s) without task logs contribute nothing")
    partial = sum(1 for a in attempts if a.status == "partial")
    if partial:
        notes.append(f"{partial} attempt(s) with partial accounting (jobs missing or not final when collected): "
                     "`ugc-wgw usage --collect` re-reads them")
    return Summary(totals, dict(by_stage), dict(by_task), dict(by_subject), clean, len(samples), sources, notes)


def select_runs(db: DB, *, mode: Optional[str], cohort: Optional[str], samples: Optional[list[str]],
                ugc_wgw_version: Optional[str]) -> list[RunRecord]:
    members: set[str] = set()
    if cohort:
        c = db.get_cohort(cohort)
        if c is None:
            raise UgcError(f"unknown cohort {cohort}")
        members = set(c.members)
    out = []
    for run in db.all_runs(ugc_wgw_version):
        if mode and run.mode != mode:
            continue
        if cohort and not (run.cohort_id == cohort or run.subject_id == cohort or run.subject_id in members):
            continue
        if samples is not None and not (run.subject_type == "sample" and run.subject_id in samples):
            continue
        out.append(run)
    return out


def build(cfg: Config, db: DB, *, mode: Optional[str], cohort: Optional[str], samples: Optional[list[str]],
          ugc_wgw_version: Optional[str], sizes: bool) -> tuple[Summary, list[AttemptUsage]]:
    inv, gpu, notes = inventory(cfg)
    cohorts = {c.cohort_id: list(c.members) for c in db.list_cohorts()}
    attempts = [attempt_usage(cfg, run, inv=inv, gpu=gpu, sizes=sizes,
                              members=cohorts.get(run.subject_id) if run.subject_type == "cohort" else None)
                for run in select_runs(db, mode=mode, cohort=cohort, samples=samples, ugc_wgw_version=ugc_wgw_version)]
    summary = aggregate(attempts)
    summary.notes = notes + summary.notes
    return summary, attempts


def collect_missing(cfg: Config, db: DB, events: Events, runs: list[RunRecord], *, refresh: bool = False) -> Counter:
    """accounting.capture for finished runs without accounting.json (or with a partial one; all of them with
    refresh); needs `sacct`. Returns the statuses."""
    if shutil.which("sacct") is None:
        raise UgcError("--collect reads SLURM accounting with sacct, which is not on PATH (run it on the cluster)")
    statuses: Counter = Counter()
    for run in runs:
        if run.status not in FINISHED:
            continue
        acct = accounting.read(run.run_path)
        if not refresh and acct is not None and acct.status == "ok":
            statuses["kept"] += 1
            continue
        statuses[accounting.capture(cfg, events, run)] += 1
    return statuses


# ---- rows and text ------------------------------------------------------------------------------

def fmt_dur(seconds: Optional[float]) -> str:
    if seconds is None:
        return "–"
    s = int(round(seconds))
    if s < 60:
        return f"{s}s"
    if s < 3600:
        return f"{s // 60}m {s % 60:02d}s"
    return f"{s // 3600}h {(s % 3600) // 60:02d}m"


def fmt_h(hours: Optional[float]) -> str:
    if hours is None:
        return "–"
    return f"{hours:.2f}" if hours < 10 else f"{hours:,.1f}"


def fmt_pct(x: Optional[float]) -> str:
    return "–" if x is None else f"{x:.0f}%"


def fmt_mb(mb: Optional[float]) -> str:
    if mb is None:
        return "–"
    return f"{mb / 1024:.1f} GB" if mb >= 1024 else f"{mb:.0f} MB"


def fmt_bytes(n: Optional[int]) -> str:
    if n is None:
        return "–"
    x = float(n)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if x < 1024 or unit == "TB":
            return f"{x:.0f} {unit}" if unit == "B" else f"{x:.1f} {unit}"
        x /= 1024
    return str(n)


def _row(u: Usage, prices: Prices, basis: str, *, sizes: bool) -> dict[str, object]:
    d: dict[str, object] = {
        "jobs": u.jobs, "failed": u.failed_jobs, "elapsed": fmt_dur(u.elapsed_s), "queue median": fmt_dur(u.queue_median),
        "queue max": fmt_dur(u.queue_max), "core-h alloc": fmt_h(u.core_h_alloc), "core-h req": fmt_h(u.core_h_req),
        "cpu-h used": fmt_h(u.cpu_h_used) if u.core_h_with_used else "–", "eff %": fmt_pct(u.efficiency),
        "gpu-h": fmt_h(u.gpu_h), "mem GB-h": fmt_h(u.mem_gb_h), "peak RSS": fmt_mb(u.peak_rss_mb),
        "disk written": fmt_mb(u.disk_write_mb) if u.disk_write_mb else "–",
    }
    if sizes:
        d["out"] = fmt_bytes(u.out_bytes)
        d["attempt dir"] = fmt_bytes(u.attempt_bytes)
    if prices.any():
        d["cost"] = prices.money(u.cost(prices, basis))
    return d


def rows_stage(s: Summary, prices: Prices, basis: str, *, mode: Optional[str], sizes: bool) -> list[dict[str, object]]:
    order = [st for st in (mode_stages(mode) if mode else []) if st in s.by_stage]
    order += sorted((st for st in s.by_stage if st not in order), key=lambda st: -s.by_stage[st].core_h(basis))
    rows = []
    for st in order:
        u = s.by_stage[st]
        rows.append({"stage": st, "attempts": u.attempts, **_row(u, prices, basis, sizes=sizes)})
    if rows:
        rows.append({"stage": "total", "attempts": s.totals.attempts, **_row(s.totals, prices, basis, sizes=sizes)})
    return rows


def rows_task(s: Summary, prices: Prices, basis: str, *, top: int) -> list[dict[str, object]]:
    rows = []
    for task, u in sorted(s.by_task.items(), key=lambda kv: -kv[1].core_h(basis))[:top]:
        d: dict[str, object] = {"task": task, "jobs": u.jobs, "failed": u.failed_jobs,
                                "mean elapsed": fmt_dur(u.elapsed_s / u.jobs if u.jobs else None),
                                "queue median": fmt_dur(u.queue_median), "core-h": fmt_h(u.core_h(basis)),
                                "cpu-h used": fmt_h(u.cpu_h_used) if u.core_h_with_used else "–", "eff %": fmt_pct(u.efficiency),
                                "gpu-h": fmt_h(u.gpu_h), "peak RSS": fmt_mb(u.peak_rss_mb), "req mem": fmt_mb(u.req_mem_mb_max),
                                "mem use %": fmt_pct(u.mem_use), "disk written": fmt_mb(u.disk_write_mb) if u.disk_write_mb else "–"}
        if prices.any():
            d["cost"] = prices.money(u.cost(prices, basis))
        rows.append(d)
    return rows


def _subject_row(label: str, kind: str, u: Usage, prices: Prices, basis: str, *, sizes: bool) -> dict[str, object]:
    d: dict[str, object] = {"subject": label, "type": kind, "core-h": fmt_h(u.core_h(basis)), "gpu-h": fmt_h(u.gpu_h),
                            "mem GB-h": fmt_h(u.mem_gb_h), "elapsed": fmt_dur(u.elapsed_s)}
    if sizes:
        d["out"] = fmt_bytes(u.out_bytes)
    if prices.any():
        d["cost"] = prices.money(u.cost(prices, basis))
        if sizes and prices.storage_gb_month is not None:
            d["storage/month"] = prices.money(u.storage_per_month(prices))
    return d


def rows_subject(s: Summary, prices: Prices, basis: str, *, sizes: bool) -> list[dict[str, object]]:
    rows = []
    for (kind, sid), u in sorted(s.by_subject.items(), key=lambda kv: (kv[0][0] != "sample", kv[0][1])):
        rows.append({**_subject_row(sid, kind, u, prices, basis, sizes=sizes), "attempts": u.attempts, "jobs": u.jobs})
    if s.n_samples:
        rows.append({**_subject_row("mean per sample, clean", f"{s.n_samples} samples", s.per_sample_clean, prices, basis, sizes=sizes),
                     "attempts": "", "jobs": ""})
        rows.append({**_subject_row("mean per sample, as run", f"{s.n_samples} samples", s.per_sample_as_run, prices, basis, sizes=sizes),
                     "attempts": "", "jobs": ""})
    return rows


def rows_attempts(attempts: list[AttemptUsage], prices: Prices, basis: str, *, sizes: bool) -> list[dict[str, object]]:
    rows = []
    for a in sorted(attempts, key=lambda a: (a.run.subject_id, a.run.stage, a.run.attempt)):
        u = a.usage
        d: dict[str, object] = {"subject": a.run.subject_id, "stage": a.run.stage, "attempt": a.run.attempt,
                                "status": a.run.status, "source": a.source + (f" ({a.status})" if a.status else ""),
                                "jobs": u.jobs, "elapsed": fmt_dur(u.elapsed_s), "queue median": fmt_dur(u.queue_median),
                                "core-h": fmt_h(u.core_h(basis)), "gpu-h": fmt_h(u.gpu_h), "peak RSS": fmt_mb(u.peak_rss_mb)}
        if sizes:
            d["out"] = fmt_bytes(u.out_bytes)
            d["attempt dir"] = fmt_bytes(u.attempt_bytes)
        if prices.any():
            d["cost"] = prices.money(u.cost(prices, basis))
        rows.append(d)
    return rows


def summary_line(s: Summary, prices: Prices, basis: str) -> str:
    src = ", ".join(f"{k} {v}" for k, v in sorted(s.sources.items()))
    t = s.totals
    text = (f"{t.attempts} attempt(s) ({src}), {s.n_samples} sample(s) with results; core-hours allocated {fmt_h(t.core_h_alloc)}"
            f" (requested {fmt_h(t.core_h_req)}), cpu-hours used {fmt_h(t.cpu_h_used) if t.core_h_with_used else '–'}"
            f" ({fmt_pct(t.efficiency)}), gpu-hours {fmt_h(t.gpu_h)}, memory {fmt_h(t.mem_gb_h)} GB-h")
    if t.queue_s:
        text += f", queue wait median {fmt_dur(t.queue_median)} max {fmt_dur(t.queue_max)}"
    if t.out_bytes is not None:
        text += f", results on disk {fmt_bytes(t.out_bytes)}"
    if s.n_samples:
        pc, pr = s.per_sample_clean, s.per_sample_as_run
        text += f"; per sample {fmt_h(pc.core_h(basis))} core-h clean, {fmt_h(pr.core_h(basis))} as run"  # type: ignore[union-attr]
    text += f"; basis {basis}"
    if prices.any():
        text += f"; cost {prices.money(t.cost(prices, basis))} at {prices.describe()}"
        if s.n_samples:
            text += f" ({prices.money(s.per_sample_clean.cost(prices, basis))} per sample, clean)"  # type: ignore[union-attr]
        if prices.storage_gb_month is not None and t.out_bytes is not None:
            text += f"; storage {prices.money(t.storage_per_month(prices))} per month"
    else:
        text += "; no prices (config.json `prices` or --price KEY=VALUE: cpu_hour, gpu_hour, mem_gb_hour, storage_gb_month, currency)"
    return text


def as_json(s: Summary, attempts: list[AttemptUsage], prices: Prices, basis: str) -> dict[str, object]:
    def with_cost(u: Usage) -> dict[str, object]:
        d = u.to_dict()
        d["cost"] = u.cost(prices, basis)
        d["storage_per_month"] = u.storage_per_month(prices)
        return d
    return {
        "basis": basis,
        "prices": {k: getattr(prices, k) for k in (*PRICE_KEYS, "currency")},
        "sources": dict(s.sources),
        "n_samples": s.n_samples,
        "notes": s.notes,
        "totals": with_cost(s.totals),
        "clean": with_cost(s.clean),
        "per_sample_clean": with_cost(s.per_sample_clean) if s.per_sample_clean else None,
        "per_sample_as_run": with_cost(s.per_sample_as_run) if s.per_sample_as_run else None,
        "by_stage": {k: with_cost(v) for k, v in s.by_stage.items()},
        "by_task": {k: with_cost(v) for k, v in s.by_task.items()},
        "by_subject": [{"type": k[0], "id": k[1], **with_cost(v)} for k, v in s.by_subject.items()],
        "attempts": [{"run_id": a.run.run_id, "subject": a.run.subject_id, "subject_type": a.run.subject_type,
                      "stage": a.run.stage, "attempt": a.run.attempt, "status": a.run.status, "source": a.source,
                      "accounting_status": a.status, **with_cost(a.usage)} for a in attempts],
    }
