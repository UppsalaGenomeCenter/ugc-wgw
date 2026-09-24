"""The foreground submit loop: create attempts, launch `miniwdl run`, poll, finalize, advance stages."""
from __future__ import annotations

import fcntl
import json
import os
import signal
import sys
import time
from dataclasses import dataclass
from typing import IO

from . import engine, failures, inputs, layout, manifest, progress, resources, slurm
from . import lease as lease_mod
from .config import Config
from .db import DB, RunRecord
from .log import LOGGER, Events
from .outputs import MissingOutputError
from .plan import Candidate, Plan, Selection, compute
from .stages import namespace, stage_spec
from .util import UgcError, diag, hostname, seconds_until, utc_now, utc_plus, utc_stamp

GRACE_SECONDS = 120.0


class Lock:
    """The project lock: a flock on `.ugc-wgw/lock` plus the lease line other hosts can read (lease.py)."""

    def __init__(self, fh: IO[str], path: "os.PathLike[str]", previous: lease_mod.Lease | None, since: str):
        self.fh, self.path, self.previous, self.since = fh, path, previous, since

    def _write(self, last_seen: str) -> None:
        self.fh.seek(0)
        self.fh.truncate()
        self.fh.write(lease_mod.render(os.getpid(), hostname(), self.since, last_seen))
        self.fh.flush()
        os.fsync(self.fh.fileno())

    def refresh(self) -> bool:
        """Heartbeat: rewrite last_seen. False when the file no longer names this driver (taken over)."""
        self.fh.seek(0)
        current = lease_mod.parse(self.fh.read())
        if current is not None and not current.is_mine():
            return False
        self._write(utc_now())
        return True

    @property
    def lease(self) -> lease_mod.Lease | None:
        return lease_mod.read(self.path)  # type: ignore[arg-type]

    def close(self) -> None:
        self.fh.close()


def acquire_lock(cfg: Config, *, takeover: bool = False) -> Lock:
    fh = open(cfg.lock_path, "a+")
    try:
        fcntl.flock(fh, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        fh.seek(0)
        holder = fh.read().strip()
        fh.close()
        raise UgcError(f"another ugc-wgw submit/retry holds {cfg.lock_path} ({holder or 'unknown holder'})") from None
    fh.seek(0)
    previous = lease_mod.parse(fh.read())
    if previous is not None and previous.host and previous.host != hostname():
        if not previous.expired(cfg.lease_seconds):
            if not takeover:
                fh.close()
                raise UgcError(f"another driver holds the lease on {previous.host} ({previous.describe()}; lease "
                               f"{int(cfg.lease_seconds)}s); wait for it to expire or use --takeover only if that "
                               "driver is known to be dead")
            LOGGER.warning("taking over a live lease: %s", previous.describe())
        else:
            LOGGER.warning("lease expired, taking over: %s", previous.describe())
    lock = Lock(fh, cfg.lock_path, previous, utc_now())
    lock._write(lock.since)
    return lock


def cancel_orphans(cfg: Config, events: Events, run: RunRecord) -> tuple[str, str]:
    """scancel the SLURM jobs a dead or killed miniwdl left behind; recorded as a run.scancel event."""
    ids = [str(i) for i in (run.meta.get("slurm_job_ids") or [])] or slurm.job_ids(run.run_path)  # type: ignore[union-attr]
    if not cfg.cancel_orphans:
        status, detail = "skipped", "cancel_orphans is false"
    else:
        status, detail = slurm.scancel(ids)
    events.emit("run.scancel", run_id=run.run_id, level="warning" if status == "failed" else "info",
                subject=run.subject_id, stage=run.stage, job_ids=ids, status=status, detail=detail)
    return status, detail


def finalize_run(cfg: Config, db: DB, events: Events, code: manifest.CodeInfo, engine_info: dict[str, str],
                 run: RunRecord, result: engine.RunResult, status_override: str | None = None,
                 finished_at: str | None = None) -> str:
    """Single finalization path (submit loop, cancellation, reconcile): manifest first, then the guarded DB update."""
    status = status_override or result.status
    error_class = result.error_class
    if status == "success" and result.outputs is not None:
        echoed = result.outputs.get(f"{namespace(run.stage)}.ugc_wgw_workflow_version")
        if echoed != run.ugc_wgw_version:
            status, error_class = "failed", "version_mismatch"
            result = engine.RunResult("failed", result.exit_code, "version_mismatch",
                                      f"workflow echoed ugc_wgw_version {echoed!r}, expected {run.ugc_wgw_version!r}")
    run.finished_at = finished_at or utc_now()
    run.error_class, run.error_kind, run.error_message, run.not_before = error_class, None, None, None
    run.meta.pop("error_evidence", None)
    if status != "success":
        cls = failures.classify(result, run.run_path)
        run.error_kind = "cancelled" if status == "cancelled" else cls.kind
        run.error_message = cls.message
        run.meta["error_evidence"] = {k: v for k, v in cls.evidence.items() if v is not None}
        if status == "failed" and run.error_kind == "transient":
            run.not_before = utc_plus(failures.backoff_delay(cfg, run.attempt), run.finished_at)
    run.meta["slurm_job_ids"] = slurm.job_ids(run.run_path)
    members = list(run.meta.get("members", []))  # type: ignore[arg-type]
    doc = manifest.build(run, result, code, engine_info, members, status)
    manifest.write(run, doc)
    changed = db.finalize_run(run.run_id, status, run.finished_at, error_class=error_class, exit_code=result.exit_code,
                              error_kind=run.error_kind, error_message=run.error_message, not_before=run.not_before,
                              meta=run.meta)
    if not changed:
        LOGGER.info("run %s was already finalized; manifest rewritten", run.run_id)
        return status
    run.status = status
    stage_path = run.run_path.parent
    try:
        layout.point_current(stage_path, run.attempt)
    except OSError as exc:
        LOGGER.warning("could not update %s/current: %s", stage_path, exc)
    events.emit(f"run.{status}", run_id=run.run_id, subject=run.subject_id, stage=run.stage, attempt=run.attempt,
                exit_code=result.exit_code, error_class=error_class, kind=run.error_kind, message=run.error_message,
                engine_message=result.error_message, task_dir=result.task_dir, node=result.node,
                slurm_job_ids=run.meta.get("slurm_job_ids", []), not_before=run.not_before)
    return status


@dataclass
class Child:
    run: RunRecord
    proc: "object"
    signalled: bool = False


class Submitter:
    def __init__(self, cfg: Config, db: DB, events: Events, code: manifest.CodeInfo, engine_info: dict[str, str],
                 ugc_wgw_version: str, sel: Selection, *, max_inflight: int, poll: float, any_version: bool = False,
                 retry_failed: bool = False, max_attempts: int | None = None, grace: float = GRACE_SECONDS,
                 lock: Lock | None = None):
        self.cfg, self.db, self.events, self.code, self.engine_info = cfg, db, events, code, engine_info
        self.lock = lock
        self.ugc_wgw_version, self.sel = ugc_wgw_version, sel
        self.max_inflight, self.poll, self.grace = max(1, max_inflight), max(0.05, poll), grace
        self.any_version, self.retry_failed, self.max_attempts = any_version, retry_failed, max_attempts
        self.children: dict[str, Child] = {}
        self.expectations = progress.TaskExpectations()
        self.progress_signature: tuple | None = None
        self.progress_at = 0.0
        self.stop = False
        self.failures = 0
        self.failed_subjects: set[tuple[str, str, str]] = set()  # subjects whose latest attempt this session failed
        self.launched = 0
        self.blocked_reported: set[str] = set()

    # ---- run creation -----------------------------------------------------

    def _context(self, cand: Candidate) -> inputs.BuildContext:
        cohort = self.db.get_cohort(cand.cohort_id) if cand.cohort_id else None
        return inputs.BuildContext(self.cfg, self.db, self.sel.mode, self.ugc_wgw_version, cand.subject_type,
                                   cand.subject_id, cohort=cohort, any_version=self.any_version)

    def create_run(self, cand: Candidate) -> RunRecord:
        ctx = self._context(cand)
        doc, members = inputs.generate(ctx, cand.stage)   # may raise MissingOutputError
        stage_path = layout.stage_dir(self.cfg.results_dir, cand.subject_type, cand.subject_id, self.ugc_wgw_version, cand.stage)
        attempt = max(self.db.next_attempt(cand.subject_type, cand.subject_id, cand.stage, self.ugc_wgw_version),
                      (layout.existing_attempts(stage_path) or [0])[-1] + 1)
        attempt_path = layout.make_attempt_dir(stage_path / f"attempt-{attempt}")
        files = layout.run_files(attempt_path)
        sha = inputs.write_inputs(files.inputs, doc)
        run = RunRecord(
            run_id=f"{cand.subject_id}-{cand.stage}-a{attempt}-{utc_stamp()}",
            subject_type=cand.subject_type, subject_id=cand.subject_id, stage=cand.stage, mode=self.sel.mode,
            ugc_wgw_version=self.ugc_wgw_version, run_dir=str(attempt_path), status="pending", attempt=attempt,
            cohort_id=cand.cohort_id, inputs_path=str(files.inputs), inputs_sha256=sha, host=hostname(),
            meta={"members": members, "resource_policy": resources.policy_ref(self.cfg)},
        )
        self.db.insert_run(run)
        self.events.emit("run.created", run_id=run.run_id, subject=run.subject_id, stage=run.stage, attempt=attempt,
                         run_dir=run.run_dir)
        return run

    def start(self, run: RunRecord) -> Child:
        cmd = engine.miniwdl_cmd(self.cfg, run.stage, run.run_path)
        try:
            proc = engine.launch(cmd, run.run_path)
        except OSError as exc:
            result = engine.RunResult("failed", None, "launch_error", str(exc))
            finalize_run(self.cfg, self.db, self.events, self.code, self.engine_info, run, result)
            self.failures += 1
            raise UgcError(f"could not launch miniwdl for {run.run_id}: {exc}") from None
        run.pid, run.started_at, run.status = proc.pid, utc_now(), "submitted"
        self.db.set_status(run.run_id, "submitted", pid=proc.pid, host=hostname(), started_at=run.started_at)
        self.events.emit("run.submitted", run_id=run.run_id, subject=run.subject_id, stage=run.stage, pid=proc.pid,
                         cmd=" ".join(cmd))
        self.launched += 1
        return Child(run, proc)

    # ---- polling ----------------------------------------------------------

    def _launch_wave(self) -> Plan:
        plan = compute(self.db, self.cfg, self.sel, self.ugc_wgw_version, any_version=self.any_version,
                       retry_failed=self.retry_failed, max_attempts=self.max_attempts)
        for b in plan.blocked:
            key = f"{b.candidate.label()}: {b.reason}"
            if key not in self.blocked_reported:
                self.blocked_reported.add(key)
                self.events.emit("run.blocked", subject=b.candidate.subject_id, stage=b.candidate.stage, reason=b.reason)
        while len(self.children) < self.max_inflight and plan.runnable and not self.stop:
            cand = plan.runnable.pop(0)
            if cand.note:
                self.events.emit("run.auto_retry", subject=cand.subject_id, stage=cand.stage, reason=cand.note)
            try:
                run = self.create_run(cand)
                child = self.start(run)
            except MissingOutputError as exc:
                self.events.emit("run.blocked", level="warning", subject=cand.subject_id, stage=cand.stage, reason=str(exc))
                continue
            except UgcError as exc:
                self.events.emit("run.blocked", level="error", subject=cand.subject_id, stage=cand.stage, reason=str(exc))
                continue
            self.children[run.run_id] = child
            plan.active.append(run)  # keep the plan's accounting right for report_progress
        return plan

    def report_progress(self, plan: Plan, *, force: bool = False) -> None:
        """Log the progress summary when the counts changed or the interval elapsed; emit submit.progress."""
        interval = float(self.cfg.progress_interval)
        if interval <= 0 and not force:
            return
        summary = progress.summarize(self.db, self.cfg, plan, self.ugc_wgw_version, mode=self.sel.mode,
                                     max_inflight=self.max_inflight, expectations=self.expectations)
        sig = summary.signature()
        now = time.monotonic()
        if not force and sig == self.progress_signature and now - self.progress_at < interval:
            return
        self.progress_signature, self.progress_at = sig, now
        for line in summary.format().splitlines():
            LOGGER.info(line)
        self.events.emit("submit.progress", **summary.as_detail())

    def poll_once(self) -> list[RunRecord]:
        finished: list[RunRecord] = []
        for run_id, child in list(self.children.items()):
            run = child.run
            files = layout.run_files(run.run_path)
            if run.status == "submitted" and files.workflow_log.exists():
                run.status = "running"
                self.db.set_status(run.run_id, "running")
                self.events.emit("run.running", run_id=run.run_id, subject=run.subject_id, stage=run.stage)
            rc = child.proc.poll()  # type: ignore[attr-defined]
            if rc is None:
                continue
            result = engine.read_result(run.run_path, rc) or engine.RunResult(
                "failed", rc, "NoResult", f"miniwdl exited {rc} without a result")
            status = finalize_run(self.cfg, self.db, self.events, self.code, self.engine_info, run, result,
                                  status_override="cancelled" if child.signalled else None)
            key = (run.subject_type, run.subject_id, run.stage)
            if status != "success":
                self.failures += 1
                self.failed_subjects.add(key)
            else:
                self.failed_subjects.discard(key)
            finished.append(run)
            del self.children[run_id]
        return finished

    def _wait(self, limit: float | None = None) -> None:
        """Sleep up to poll seconds (or `limit`, when shorter), waking early on a stop request or a child exit."""
        deadline = time.monotonic() + (self.poll if limit is None else max(0.0, min(self.poll, limit)))
        while time.monotonic() < deadline and not self.stop:
            if any(c.proc.poll() is not None for c in self.children.values()):  # type: ignore[attr-defined]
                return
            time.sleep(min(0.2, max(0.0, deadline - time.monotonic())))

    def stop_all(self) -> None:
        for child in self.children.values():
            if not child.signalled:
                try:
                    os.killpg(os.getpgid(child.proc.pid), signal.SIGTERM)  # type: ignore[attr-defined]
                except (ProcessLookupError, PermissionError):
                    pass
                child.signalled = True
                self.events.emit("run.terminating", run_id=child.run.run_id, subject=child.run.subject_id, stage=child.run.stage)
        deadline = time.monotonic() + self.grace
        while self.children and time.monotonic() < deadline:
            self.poll_once()
            if self.children:
                time.sleep(1.0)
        for child in list(self.children.values()):
            try:
                os.killpg(os.getpgid(child.proc.pid), signal.SIGKILL)  # type: ignore[attr-defined]
            except (ProcessLookupError, PermissionError):
                pass
            child.proc.wait(timeout=10)  # type: ignore[attr-defined]
            result = engine.RunResult("failed", child.proc.returncode, "killed", "SIGKILL after grace period")  # type: ignore[attr-defined]
            finalize_run(self.cfg, self.db, self.events, self.code, self.engine_info, child.run, result, status_override="cancelled")
            cancel_orphans(self.cfg, self.events, child.run)  # miniwdl could not scancel its own jobs
            del self.children[child.run.run_id]

    # ---- the loop ---------------------------------------------------------

    def run_loop(self) -> int:
        def on_signal(signum: int, _frame: object) -> None:
            if not self.stop:
                diag(f"signal {signum}: stopping (terminating {len(self.children)} run(s), grace {int(self.grace)}s)")
            self.stop = True

        previous = {s: signal.signal(s, on_signal) for s in (signal.SIGINT, signal.SIGTERM)}
        self.events.emit("submit.start", mode=self.sel.mode, stage=self.sel.stage, cohort=self.sel.cohort,
                         max_inflight=self.max_inflight, ugc_wgw_version=self.ugc_wgw_version, pid=os.getpid())
        try:
            while True:
                plan = self._launch_wave() if not self.stop else None
                if plan is not None:
                    self.report_progress(plan)
                if self.lock is not None and not self.stop and not self.lock.refresh():
                    diag("the project lease was taken over by another driver; stopping")
                    self.stop = True
                if not self.children:
                    if self.stop or plan is None:
                        break
                    wake = plan.next_wake()
                    if wake is not None:
                        key = f"backing off until {wake}"
                        if key not in self.blocked_reported:
                            self.blocked_reported.add(key)
                            diag(f"nothing runnable; {len(plan.backing_off)} subject(s) backing off, next at {wake}")
                        self._wait(limit=seconds_until(wake))
                        continue
                    self.report_progress(plan, force=True)
                    self._summary(plan)
                    break
                self._wait()
                self.poll_once()
                if self.stop and self.children:
                    self.stop_all()
        finally:
            for s, h in previous.items():
                signal.signal(s, h)
            self.events.emit("submit.stop", launched=self.launched, failures=self.failures,
                             failed_subjects=len(self.failed_subjects), interrupted=self.stop)
        if self.stop:
            return 130
        return 1 if self.failed_subjects else 0

    def _summary(self, plan: Plan) -> None:
        diag(f"nothing runnable: {len(plan.done)} done, {len(plan.blocked)} blocked, "
             f"{len(plan.active)} active elsewhere; launched {self.launched}, failures {self.failures}")
        for b in plan.blocked:
            diag(f"  {'backoff' if b.not_before else 'blocked'}: {b.candidate.label()}: {b.reason}")
        for r in plan.active:
            diag(f"  active:  {r.subject_type} {r.subject_id} / {r.stage} attempt {r.attempt} ({r.status}, pid {r.pid} on {r.host})")


def dry_run(cfg: Config, db: DB, ugc_wgw_version: str, sel: Selection, *, any_version: bool = False,
            retry_failed: bool = False, max_attempts: int | None = None, out: IO[str] | None = None) -> int:
    """Print the first wave (inputs and command lines) plus blocked/active lists. No side effects."""
    out = out or sys.stdout
    plan = compute(db, cfg, sel, ugc_wgw_version, any_version=any_version, retry_failed=retry_failed, max_attempts=max_attempts)
    print(f"# ugc-wgw dry-run: mode={sel.mode} version={ugc_wgw_version} runnable={len(plan.runnable)} "
          f"blocked={len(plan.blocked)} active={len(plan.active)} done={len(plan.done)}", file=out)
    for cand in plan.runnable:
        cohort = db.get_cohort(cand.cohort_id) if cand.cohort_id else None
        ctx = inputs.BuildContext(cfg, db, sel.mode, ugc_wgw_version, cand.subject_type, cand.subject_id, cohort=cohort,
                                  any_version=any_version)
        stage_path = layout.stage_dir(cfg.results_dir, cand.subject_type, cand.subject_id, ugc_wgw_version, cand.stage)
        attempt = max(db.next_attempt(cand.subject_type, cand.subject_id, cand.stage, ugc_wgw_version),
                      (layout.existing_attempts(stage_path) or [0])[-1] + 1)
        attempt_path = stage_path / f"attempt-{attempt}"
        print(f"\n## {cand.label()} -> {attempt_path}", file=out)
        try:
            doc, _ = inputs.generate(ctx, cand.stage)
        except UgcError as exc:
            print(f"# blocked: {exc}", file=out)
            continue
        print(json.dumps(doc, indent=2, sort_keys=True), file=out)
        print(" ".join(engine.miniwdl_cmd(cfg, cand.stage, attempt_path)), file=out)
    for b in plan.blocked:
        print(f"\n# blocked: {b.candidate.label()}: {b.reason}", file=out)
    for r in plan.active:
        print(f"\n# active: {r.subject_type} {r.subject_id} / {r.stage} attempt {r.attempt} ({r.status})", file=out)
    return 0
