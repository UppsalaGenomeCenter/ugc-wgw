"""What is runnable now: mode-derived gating over the state DB (docs/DESIGN.md §5.1, §8.1)."""
from __future__ import annotations

from dataclasses import dataclass, field

from .config import Config
from .db import DB, Cohort, RunRecord
from .stages import mode_stages, stage_spec
from .util import UgcError, utc_now


@dataclass
class Selection:
    mode: str
    stage: str | None = None
    samples: list[str] | None = None
    cohort: str | None = None


@dataclass
class Candidate:
    stage: str
    subject_type: str
    subject_id: str
    cohort_id: str | None = None
    note: str | None = None   # why a previously failed subject is runnable again (automatic retry)

    def label(self) -> str:
        return f"{self.subject_type} {self.subject_id} / {self.stage}"


@dataclass
class Blocked:
    candidate: Candidate
    reason: str
    not_before: str | None = None   # set when the subject is only waiting out a backoff


@dataclass
class Plan:
    runnable: list[Candidate] = field(default_factory=list)
    blocked: list[Blocked] = field(default_factory=list)
    active: list[RunRecord] = field(default_factory=list)
    done: list[Candidate] = field(default_factory=list)

    @property
    def backing_off(self) -> list[Blocked]:
        return [b for b in self.blocked if b.not_before]

    def next_wake(self) -> str | None:
        """Earliest not_before among backing-off subjects, or None."""
        times = [b.not_before for b in self.blocked if b.not_before]
        return min(times) if times else None


def _cohort(db: DB, sel: Selection) -> Cohort | None:
    if not sel.cohort:
        return None
    cohort = db.get_cohort(sel.cohort)
    if cohort is None:
        raise UgcError(f"unknown cohort {sel.cohort}")
    return cohort


def _subjects(db: DB, sel: Selection, cohort: Cohort | None) -> list[str]:
    if sel.samples:
        unknown = [s for s in sel.samples if not db.sample_exists(s)]
        if unknown:
            raise UgcError("unknown samples: " + ", ".join(unknown))
        if cohort is not None:
            outside = [s for s in sel.samples if s not in cohort.members]
            if outside:
                raise UgcError(f"not members of cohort {cohort.cohort_id}: " + ", ".join(outside))
        return list(sel.samples)
    if cohort is not None:
        return list(cohort.members)
    return db.list_sample_ids()


def prerequisites(db: DB, mode: str, stage: str, subject_type: str, subject_id: str,
                  cohort: Cohort | None, version: str | None) -> str | None:
    """None when every earlier stage of the mode is satisfied, else the reason it is not."""
    for earlier in mode_stages(mode):
        if earlier == stage:
            break
        espec = stage_spec(earlier)
        if espec.subject_type == "sample":
            members = [subject_id] if subject_type == "sample" else (cohort.members if cohort else [])
            missing = [m for m in members if db.latest_success("sample", m, earlier, version) is None]
            if missing:
                shown = ", ".join(missing[:5]) + (f" (+{len(missing) - 5} more)" if len(missing) > 5 else "")
                return f"waiting for {earlier} of {shown}"
        else:
            if cohort is None:
                return f"{earlier} needs --cohort"
            if db.latest_success("cohort", cohort.cohort_id, earlier, version) is None:
                return f"waiting for {earlier} of cohort {cohort.cohort_id}"
    return None


def compute(db: DB, cfg: Config, sel: Selection, ugc_wgw_version: str, *, any_version: bool = False,
            retry_failed: bool = False, max_attempts: int | None = None, now: str | None = None) -> Plan:
    """`now` (UTC_FMT) decides whether a backoff has elapsed; tests inject it."""
    now = now or utc_now()
    stages = mode_stages(sel.mode)
    if sel.stage:
        if sel.stage not in stages:
            raise UgcError(f"stage {sel.stage} is not part of mode {sel.mode} ({' -> '.join(stages)})")
        selected = (sel.stage,)
    else:
        selected = stages
    cohort = _cohort(db, sel)
    if sel.mode == "joint" and cohort is None:
        raise UgcError("joint mode needs --cohort (cohort_call and downstream run in a cohort context)")
    samples = _subjects(db, sel, cohort)
    lookup_version = None if any_version else ugc_wgw_version
    plan = Plan()

    for stage in selected:
        spec = stage_spec(stage)
        if spec.subject_type == "sample":
            candidates = [Candidate(stage, "sample", s, cohort.cohort_id if cohort else None) for s in samples]
        elif cohort is None:
            plan.blocked.append(Blocked(Candidate(stage, "cohort", "?", None), f"{stage} needs --cohort <id>"))
            continue
        else:
            candidates = [Candidate(stage, "cohort", cohort.cohort_id, cohort.cohort_id)]
        for cand in candidates:
            latest = db.latest_run(cand.subject_type, cand.subject_id, stage, ugc_wgw_version)
            if latest is not None and latest.status == "success":
                plan.done.append(cand)
                continue
            if latest is not None and latest.is_active:
                plan.active.append(latest)
                continue
            if latest is not None:  # failed or cancelled
                if retry_failed:  # `ugc-wgw retry`: the manual path ignores the kind and any backoff
                    if max_attempts is not None and latest.attempt >= max_attempts:
                        plan.blocked.append(Blocked(cand, f"{latest.status} after {latest.attempt} attempts (max {max_attempts})"))
                        continue
                else:
                    auto = (latest.status == "failed" and latest.error_kind == "transient"
                            and latest.attempt <= cfg.auto_retry_max)
                    if not auto:
                        extra = ", auto-retries exhausted" if latest.error_kind == "transient" else ""
                        plan.blocked.append(Blocked(cand, f"{latest.status} (attempt {latest.attempt}, "
                                                    f"{latest.error_class or 'no class'}, {latest.error_kind or 'no kind'}{extra}); "
                                                    f"run `ugc-wgw retry --stage {stage}`"))
                        continue
                    if latest.not_before and latest.not_before > now:
                        plan.blocked.append(Blocked(cand, f"backing off until {latest.not_before} after transient "
                                                    f"{latest.error_class} (attempt {latest.attempt})", not_before=latest.not_before))
                        continue
                    cand.note = (f"auto-retry after transient {latest.error_class} "
                                 f"(attempt {latest.attempt} of {cfg.auto_retry_max + 1})")
            reason = prerequisites(db, sel.mode, stage, cand.subject_type, cand.subject_id, cohort, lookup_version)
            if reason:
                plan.blocked.append(Blocked(cand, reason))
            else:
                plan.runnable.append(cand)
    return plan
