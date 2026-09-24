"""Progress of a selection: done/total per stage, what the active attempts are doing, and an ETA.

Counts are exact (every candidate of the plan is done, active, blocked or runnable). Times are estimates:
the median wall time of finished runs of the same stage at this version, scaled by what is left and by
the driver's concurrency; active attempts contribute the fraction of tasks they have finished according to
miniwdl's workflow.log.json, against the task count of a completed run of the same stage.
"""
from __future__ import annotations

import math
import statistics
from dataclasses import dataclass, field

from .config import Config
from .db import DB, RunRecord
from .layout import run_files
from .plan import Plan
from .report import parse_workflow_log
from .stages import mode_stages, stage_spec
from .util import parse_utc


@dataclass
class ActiveInfo:
    subject_id: str
    tasks_done: int
    tasks_expected: int | None

    @property
    def fraction(self) -> float | None:
        if not self.tasks_expected:
            return None
        return min(1.0, self.tasks_done / self.tasks_expected)

    def label(self) -> str:
        if self.tasks_expected:
            return f"{self.subject_id} {self.tasks_done}/{self.tasks_expected} tasks"
        return f"{self.subject_id} {self.tasks_done} tasks"


@dataclass
class StageProgress:
    stage: str
    total: int = 0
    done: int = 0
    runnable: int = 0
    waiting: int = 0
    backing_off: int = 0
    failed: int = 0
    active: list[ActiveInfo] = field(default_factory=list)
    median: float | None = None
    eta: float | None = None

    @property
    def remaining(self) -> int:
        """Runs still to complete that the driver can do without operator action."""
        return self.total - self.done - self.failed


@dataclass
class Progress:
    mode: str
    stages: list[StageProgress]

    @property
    def total(self) -> int:
        return sum(s.total for s in self.stages)

    @property
    def done(self) -> int:
        return sum(s.done for s in self.stages)

    @property
    def eta(self) -> float | None:
        etas = [s.eta for s in self.stages if s.eta is not None]
        return sum(etas) if etas else None

    def signature(self) -> tuple:
        return tuple((s.stage, s.done, s.failed, s.backing_off, len(s.active), s.waiting) for s in self.stages)

    def format(self) -> str:
        pct = f"{100 * self.done // self.total}%" if self.total else "–"
        lines = [f"progress mode={self.mode} {self.done}/{self.total} ({pct}) eta≈{fmt_seconds(self.eta)}"]
        width = max((len(s.stage) for s in self.stages), default=8)
        for s in self.stages:
            bits = [f"done {s.done}/{s.total}"]
            if s.active:
                bits.append(f"active {len(s.active)} (" + ", ".join(a.label() for a in s.active) + ")")
            if s.runnable:
                bits.append(f"runnable {s.runnable}")
            if s.waiting:
                bits.append(f"waiting {s.waiting}")
            if s.backing_off:
                bits.append(f"backing off {s.backing_off}")
            if s.failed:
                bits.append(f"failed {s.failed} (needs retry)")
            bits.append(f"median {fmt_seconds(s.median)}")
            bits.append(f"eta≈{fmt_seconds(s.eta)}")
            lines.append(f"  {s.stage.ljust(width)}  " + "  ".join(bits))
        return "\n".join(lines)

    def as_detail(self) -> dict[str, object]:
        return {
            "done": self.done, "total": self.total, "eta_seconds": None if self.eta is None else int(self.eta),
            "stages": {s.stage: {"done": s.done, "total": s.total, "active": len(s.active), "runnable": s.runnable,
                                 "waiting": s.waiting, "backing_off": s.backing_off, "failed": s.failed,
                                 "median_seconds": None if s.median is None else int(s.median),
                                 "eta_seconds": None if s.eta is None else int(s.eta)} for s in self.stages},
        }


def fmt_seconds(seconds: float | None) -> str:
    if seconds is None:
        return "–"
    s = int(round(seconds))
    if s < 60:
        return f"{s}s"
    if s < 3600:
        return f"{s // 60}m{s % 60:02d}s"
    if s < 86400:
        return f"{s // 3600}h{(s % 3600) // 60:02d}m"
    return f"{s // 86400}d{(s % 86400) // 3600:02d}h"


def _duration(run: RunRecord) -> float | None:
    if not run.started_at or not run.finished_at:
        return None
    return max(0.0, (parse_utc(run.finished_at) - parse_utc(run.started_at)).total_seconds())


class TaskExpectations:
    """Task count of a completed run per stage, learnt lazily and cached for the driver's lifetime."""

    def __init__(self) -> None:
        self.expected: dict[str, int | None] = {}

    def get(self, db: DB, stage: str, ugc_wgw_version: str) -> int | None:
        if stage not in self.expected:
            self.expected[stage] = None
            for run in db.finished_runs(stage, ugc_wgw_version, "success", limit=3):
                tasks, _ = parse_workflow_log(run_files(run.run_path).workflow_log_json)
                if tasks:
                    self.expected[stage] = len(tasks)
                    break
        return self.expected[stage]


def summarize(db: DB, cfg: Config, plan: Plan, ugc_wgw_version: str, *, mode: str, max_inflight: int,
              expectations: TaskExpectations | None = None) -> Progress:
    expectations = expectations or TaskExpectations()
    stages = {s: StageProgress(s) for s in mode_stages(mode)}
    for cand in plan.done:
        stages.setdefault(cand.stage, StageProgress(cand.stage))
        stages[cand.stage].total += 1
        stages[cand.stage].done += 1
    for cand in plan.runnable:
        stages.setdefault(cand.stage, StageProgress(cand.stage))
        stages[cand.stage].total += 1
        stages[cand.stage].runnable += 1
    for b in plan.blocked:
        sp = stages.setdefault(b.candidate.stage, StageProgress(b.candidate.stage))
        sp.total += 1
        if b.not_before:
            sp.backing_off += 1
        elif b.reason.startswith(("failed", "cancelled")):
            sp.failed += 1
        else:
            sp.waiting += 1
    for run in plan.active:
        sp = stages.setdefault(run.stage, StageProgress(run.stage))
        sp.total += 1
        tasks, _ = parse_workflow_log(run_files(run.run_path).workflow_log_json)
        sp.active.append(ActiveInfo(run.subject_id, sum(1 for t in tasks if t.end is not None),
                                    expectations.get(db, run.stage, ugc_wgw_version)))
    out = []
    for stage, sp in stages.items():
        if sp.total == 0:
            continue
        durations = [d for d in (_duration(r) for r in db.finished_runs(stage, ugc_wgw_version, "success", limit=20))
                     if d is not None]
        sp.median = statistics.median(durations) if durations else None
        if sp.remaining == 0:
            sp.eta = 0.0
        elif sp.median is not None:
            parallel = 1 if stage_spec(stage).subject_type == "cohort" else max(1, max_inflight)
            credit = sum(a.fraction or 0.0 for a in sp.active)  # work already done inside active attempts
            sp.eta = max(0.0, (math.ceil(sp.remaining / parallel) - credit / parallel) * sp.median)
        out.append(sp)
    return Progress(mode, out)
