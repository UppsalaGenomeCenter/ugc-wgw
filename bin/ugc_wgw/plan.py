"""What is runnable now: mode-derived gating over the state DB (docs/DESIGN.md §5.1, §8.1)."""
from __future__ import annotations

from dataclasses import dataclass, field

from .config import Config
from .db import DB, Cohort, RunRecord
from .stages import ENTRY_KINDS, member_sequence, mode_stages, sample_stages, stage_spec
from .util import UgcError, utc_now


@dataclass
class Selection:
    mode: str
    stage: str | None = None
    samples: list[str] | None = None
    cohort: str | None = None
    allow_standalone_members: bool = False   # joint mode: a vcf member without a gVCF phases its own VCFs and skips cohort_call
    rerun: bool = False   # `submit --rerun`: a new attempt of the one named stage for subjects whose latest run succeeded


@dataclass
class Candidate:
    stage: str
    subject_type: str
    subject_id: str
    cohort_id: str | None = None
    note: str | None = None   # why a previously failed subject is runnable again (automatic retry)
    rerun_of: int | None = None   # `submit --rerun`: the successful attempt this candidate runs again

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


def _kind(db: DB, sid: str) -> tuple[str, bool]:
    """(entry kind, has a gVCF) of a registered sample; reads for an unknown one."""
    rec = db.get_sample(sid)
    return (rec.entry_kind or "reads", bool(rec.gvcf) or rec.entry_kind != "vcf") if rec else ("reads", True)


def member_stages(db: DB, mode: str, sid: str, *, allow_standalone: bool = False) -> tuple[str, ...]:
    """The stages member `sid` takes part in, by its entry kind (stages.member_sequence)."""
    kind, has_gvcf = _kind(db, sid)
    return member_sequence(mode, kind, has_gvcf=has_gvcf, allow_standalone=allow_standalone)


def planned_stages(db: DB, mode: str, samples: list[str]) -> tuple[str, ...]:
    """The mode's stages with the per-sample stages of every entry kind among `samples` in the place of the mode's
    (reads: singleton or upstream; bam: upstream and downstream; vcf: intake and downstream), in order."""
    kinds = sorted({_kind(db, sid)[0] for sid in samples} or {"reads"}, key=ENTRY_KINDS.index)
    out: list[str] = []
    for stage in mode_stages(mode):
        if stage_spec(stage).subject_type == "sample" and stage in ("singleton", "upstream", "assembly"):
            for kind in kinds:
                own = sample_stages(mode, kind)
                head = own[:1] if mode == "joint" else own
                out.extend(s for s in head if s not in out)
            continue
        if stage not in out:
            out.append(stage)
    return tuple(out)


def gvcf_blockers(db: DB, mode: str, cohort: Cohort | None) -> list[str]:
    """Joint mode: the vcf members without a gVCF (cohort_call needs one from every member)."""
    if mode != "joint" or cohort is None:
        return []
    return [sid for sid in cohort.members if _kind(db, sid) == ("vcf", False)]


def prerequisites(db: DB, mode: str, stage: str, subject_type: str, subject_id: str,
                  cohort: Cohort | None, version: str | None, *, allow_standalone: bool = False) -> str | None:
    """None when every earlier stage is satisfied, else the reason it is not. A sample's earlier stages are those of
    its own sequence (its entry kind's per-sample stages, the mode's cohort stages); a cohort stage waits for every
    member's per-sample stages before it and for the cohort stages before it (the optional ones, cohort_qc and
    cohort_tr, never gate)."""
    if subject_type == "sample":
        seq = member_stages(db, mode, subject_id, allow_standalone=allow_standalone)
        if stage not in seq:
            kind, _ = _kind(db, subject_id)
            return f"{stage} is not a stage of a {kind} sample in {mode} mode ({' -> '.join(seq) or 'none'})"
        for earlier in seq:
            if earlier == stage:
                break
            espec = stage_spec(earlier)
            if espec.optional:
                continue
            if espec.subject_type == "sample":
                if db.latest_success("sample", subject_id, earlier, version) is None:
                    return f"waiting for {earlier} of {subject_id}"
            else:
                if cohort is None:
                    return f"{earlier} needs --cohort"
                if db.latest_success("cohort", cohort.cohort_id, earlier, version) is None:
                    return f"waiting for {earlier} of cohort {cohort.cohort_id}"
        return None
    if cohort is None:
        return f"{stage} needs --cohort"
    if stage == "cohort_call" and not allow_standalone:
        blockers = gvcf_blockers(db, mode, cohort)
        if blockers:
            shown = ", ".join(blockers[:5]) + (f" (+{len(blockers) - 5} more)" if len(blockers) > 5 else "")
            return (f"joint mode needs a gVCF from every member; vcf member(s) without one: {shown} (register the gvcf column, "
                    "or submit with --allow-standalone-members to phase them from their own VCFs outside the joint call)")
    missing: list[str] = []
    for sid in cohort.members:
        seq = member_stages(db, mode, sid, allow_standalone=allow_standalone)
        if stage not in seq:
            continue   # a member with no place in this cohort stage (a vcf member left out of cohort_call)
        for earlier in seq:
            if earlier == stage:
                break
            if stage_spec(earlier).subject_type == "sample" and db.latest_success("sample", sid, earlier, version) is None:
                missing.append(f"{earlier} of {sid}")
                break
    if missing:
        shown = ", ".join(missing[:5]) + (f" (+{len(missing) - 5} more)" if len(missing) > 5 else "")
        return f"waiting for {shown}"
    for earlier in mode_stages(mode):
        if earlier == stage:
            break
        espec = stage_spec(earlier)
        if espec.subject_type != "cohort" or espec.optional:
            continue
        if db.latest_success("cohort", cohort.cohort_id, earlier, version) is None:
            return f"waiting for {earlier} of cohort {cohort.cohort_id}"
    return None


def check_rerun(sel: Selection) -> None:
    """`--rerun` names what it repeats: one stage, and for a per-sample stage the samples. Never a whole mode, never
    every member of a cohort by default."""
    if not sel.rerun:
        return
    if not sel.stage:
        raise UgcError("--rerun needs --stage: it runs one stage again, never a whole mode")
    if stage_spec(sel.stage).subject_type == "sample" and not sel.samples:
        raise UgcError(f"--rerun of the per-sample stage {sel.stage} needs --samples or --samples-file: name the samples to run again")


def rerun_followups(db: DB, sel: Selection, ugc_wgw_version: str, *, any_version: bool = False) -> list[str]:
    """The stages after `sel.stage` that already succeeded for the rerun's subjects and so keep the outputs of the
    attempt being replaced until they are run again themselves: `stage (n samples)` or `stage (cohort id)`, in the
    mode's order. Empty for an optional stage, which no stage reads."""
    if not sel.rerun or not sel.stage or stage_spec(sel.stage).optional:
        return []
    cohort = _cohort(db, sel)
    version = None if any_version else ugc_wgw_version
    samples = _subjects(db, sel, cohort)
    stages = planned_stages(db, sel.mode, samples)
    if sel.stage not in stages:
        return []
    out: list[str] = []
    for later in stages[stages.index(sel.stage) + 1:]:
        if stage_spec(later).subject_type == "sample":
            n = sum(1 for sid in samples
                    if later in member_stages(db, sel.mode, sid, allow_standalone=sel.allow_standalone_members)
                    and db.latest_success("sample", sid, later, version) is not None)
            if n:
                out.append(f"{later} ({n} sample{'s' if n != 1 else ''})")
        elif cohort is not None and db.latest_success("cohort", cohort.cohort_id, later, version) is not None:
            out.append(f"{later} (cohort {cohort.cohort_id})")
    return out


def compute(db: DB, cfg: Config, sel: Selection, ugc_wgw_version: str, *, any_version: bool = False,
            retry_failed: bool = False, max_attempts: int | None = None, now: str | None = None,
            rerun_started: set[tuple[str, str, str]] | None = None) -> Plan:
    """`now` (UTC_FMT) decides whether a backoff has elapsed; tests inject it. With `sel.rerun`, a subject whose
    latest run of the stage succeeded is a candidate again (once: `rerun_started` holds the subjects this invocation
    has already given their new attempt, as (subject_type, subject_id, stage)); a failed or active latest run is
    handled as always."""
    now = now or utc_now()
    check_rerun(sel)
    cohort = _cohort(db, sel)
    if sel.mode == "joint" and cohort is None:
        raise UgcError("joint mode needs --cohort (cohort_call and downstream run in a cohort context)")
    samples = _subjects(db, sel, cohort)
    stages = planned_stages(db, sel.mode, samples)
    if sel.stage:
        if sel.stage not in stages:
            raise UgcError(f"stage {sel.stage} is not part of mode {sel.mode} for these samples ({' -> '.join(stages)})")
        selected = (sel.stage,)
    else:
        selected = stages
    lookup_version = None if any_version else ugc_wgw_version
    plan = Plan()

    allow = sel.allow_standalone_members
    for stage in selected:
        spec = stage_spec(stage)
        if spec.subject_type == "sample":
            candidates = [Candidate(stage, "sample", s, cohort.cohort_id if cohort else None) for s in samples
                          if stage in member_stages(db, sel.mode, s, allow_standalone=allow)]
            for s_id in samples:
                if stage not in member_stages(db, sel.mode, s_id, allow_standalone=allow) and sel.stage:
                    kind, _ = _kind(db, s_id)
                    plan.blocked.append(Blocked(Candidate(stage, "sample", s_id, cohort.cohort_id if cohort else None),
                                                f"not a stage of a {kind} sample in {sel.mode} mode"))
        elif cohort is None:
            plan.blocked.append(Blocked(Candidate(stage, "cohort", "?", None), f"{stage} needs --cohort <id>"))
            continue
        else:
            candidates = [Candidate(stage, "cohort", cohort.cohort_id, cohort.cohort_id)]
        for cand in candidates:
            latest = db.latest_run(cand.subject_type, cand.subject_id, stage, ugc_wgw_version)
            if latest is not None and latest.status == "success":
                if not sel.rerun or (cand.subject_type, cand.subject_id, stage) in (rerun_started or ()):
                    plan.done.append(cand)
                    continue
                cand.rerun_of = latest.attempt
                cand.note = f"rerun of the successful attempt {latest.attempt}"
                reason = prerequisites(db, sel.mode, stage, cand.subject_type, cand.subject_id, cohort, lookup_version, allow_standalone=allow)
                if reason:
                    plan.blocked.append(Blocked(cand, reason))
                else:
                    plan.runnable.append(cand)
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
                    esc = latest.meta.get("escalation") if latest.status == "failed" and latest.error_kind == "resource" else None
                    if isinstance(esc, dict) and esc.get("task"):
                        # the task ran out of time or memory and failures/escalate planned a bigger request: no backoff
                        cand.note = (f"escalation: {esc['task']} {esc['resource']} ×{float(esc['factor']):g} "
                                     f"(escalation {esc['n']} of {esc['max']})")
                        reason = prerequisites(db, sel.mode, stage, cand.subject_type, cand.subject_id, cohort, lookup_version, allow_standalone=allow)
                        if reason:
                            plan.blocked.append(Blocked(cand, reason))
                        else:
                            plan.runnable.append(cand)
                        continue
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
            reason = prerequisites(db, sel.mode, stage, cand.subject_type, cand.subject_id, cohort, lookup_version, allow_standalone=allow)
            if reason:
                plan.blocked.append(Blocked(cand, reason))
            else:
                plan.runnable.append(cand)
    return plan
