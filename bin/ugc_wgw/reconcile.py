"""After a driver restart or crash: settle runs the DB still calls active from what their run dirs say."""
from __future__ import annotations

from . import engine, manifest
from . import lease as lease_mod
from .config import Config
from .db import DB, RunRecord
from .layout import run_files
from .log import LOGGER, Events
from .submit import cancel_orphans, finalize_run
from .util import hostname, iso_from_mtime, pid_alive


def reconcile(cfg: Config, db: DB, events: Events, code: manifest.CodeInfo, engine_info: dict[str, str],
              dry: bool = False) -> list[RunRecord]:
    """Returns the runs that were (or, when dry, would be) finalized."""
    settled: list[RunRecord] = []
    lease = lease_mod.read(cfg.lock_path)
    holder_is_me = lease is not None and lease.is_mine()

    def lost(run: RunRecord, message: str, source: str) -> None:
        result = engine.RunResult("failed", None, "driver_lost", message)
        finalize_run(cfg, db, events, code, engine_info, run, result)
        events.emit("run.reconciled", run_id=run.run_id, subject=run.subject_id, stage=run.stage,
                    status="failed", source=source)
        cancel_orphans(cfg, events, run)

    for run in db.active_runs():
        files = run_files(run.run_path)
        result = engine.read_result(run.run_path, None)
        if result is not None:
            finished_at = iso_from_mtime(files.run_json if files.run_json.exists() else
                                         files.outputs if files.outputs.exists() else files.error)
            if not dry:
                finalize_run(cfg, db, events, code, engine_info, run, result, finished_at=finished_at)
                events.emit("run.reconciled", run_id=run.run_id, subject=run.subject_id, stage=run.stage,
                            status=run.status, source="run directory")
            settled.append(run)
            continue
        same_host = (run.host or hostname()) == hostname()
        if same_host and not pid_alive(run.pid):
            if not dry:
                lost(run, "driver process disappeared before the run finished; no result in the run directory",
                     "dead pid")
            settled.append(run)
        elif same_host:
            LOGGER.info("run %s still has a live driver pid %s", run.run_id, run.pid)
        else:
            # Another host: trust the lease. The driver that recorded the run held the lease; if nobody holds a
            # fresh one for that host any more (or this process holds it now), that driver is gone.
            stale = lease is None or holder_is_me or (lease.host == run.host and lease.expired(cfg.lease_seconds))
            if stale:
                if not dry:
                    lost(run, f"driver on host {run.host} is gone (lease expired or taken over); "
                              "no result in the run directory", "lease expired")
                settled.append(run)
            else:
                LOGGER.warning("run %s is recorded as %s on host %s; its lease is still fresh (%s); left alone",
                               run.run_id, run.status, run.host, lease.describe() if lease else "none")
    return settled
