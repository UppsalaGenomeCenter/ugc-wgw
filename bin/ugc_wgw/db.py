"""SQLite state of a ugc-wgw project (docs/DESIGN.md §8.3), with a few extra columns on `runs`.

Rules: every multi-statement write goes through `DB.tx()` (BEGIN IMMEDIATE); terminal
runs are never mutated (`finalize_run` is guarded and idempotent); state is never
inferred from the filesystem alone.
"""
from __future__ import annotations

import json
import sqlite3
from contextlib import contextmanager
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator

from .util import UgcError

SCHEMA_VERSION = 2

# Forward-only migrations applied in order on open; each entry brings an older file to that version.
MIGRATIONS: dict[int, tuple[str, ...]] = {
    2: (
        "ALTER TABLE runs ADD COLUMN error_kind TEXT",
        "ALTER TABLE runs ADD COLUMN error_message TEXT",
        "ALTER TABLE runs ADD COLUMN not_before TEXT",
    ),
}

ACTIVE = ("pending", "submitted", "running")
TERMINAL = ("success", "failed", "cancelled")
STATUSES = ACTIVE + TERMINAL

SCHEMA = """
CREATE TABLE IF NOT EXISTS schema_version (version INTEGER NOT NULL);
CREATE TABLE IF NOT EXISTS samples (
  sample_id TEXT PRIMARY KEY,
  sex TEXT,
  father_id TEXT,
  mother_id TEXT,
  added_at TEXT NOT NULL,
  meta_json TEXT NOT NULL DEFAULT '{}'
);
CREATE TABLE IF NOT EXISTS sample_inputs (
  sample_id TEXT NOT NULL REFERENCES samples(sample_id) ON DELETE CASCADE,
  kind TEXT NOT NULL CHECK (kind IN ('hifi_reads', 'fail_reads')),
  position INTEGER NOT NULL,
  path TEXT NOT NULL,
  checksum TEXT,
  PRIMARY KEY (sample_id, kind, position)
);
CREATE TABLE IF NOT EXISTS cohorts (
  cohort_id TEXT PRIMARY KEY,
  frozen_at TEXT NOT NULL,
  sample_list_sha256 TEXT NOT NULL,
  meta_json TEXT NOT NULL DEFAULT '{}'
);
CREATE TABLE IF NOT EXISTS cohort_members (
  cohort_id TEXT NOT NULL REFERENCES cohorts(cohort_id),
  sample_id TEXT NOT NULL REFERENCES samples(sample_id),
  position INTEGER NOT NULL,
  PRIMARY KEY (cohort_id, sample_id)
);
CREATE TABLE IF NOT EXISTS runs (
  run_id TEXT PRIMARY KEY,
  subject_type TEXT NOT NULL CHECK (subject_type IN ('sample', 'cohort')),
  subject_id TEXT NOT NULL,
  stage TEXT NOT NULL,
  mode TEXT NOT NULL,
  ugc_wgw_version TEXT NOT NULL,
  cohort_id TEXT,
  inputs_path TEXT,
  inputs_sha256 TEXT,
  run_dir TEXT NOT NULL,
  status TEXT NOT NULL CHECK (status IN ('pending', 'submitted', 'running', 'success', 'failed', 'cancelled')),
  attempt INTEGER NOT NULL,
  error_class TEXT,
  exit_code INTEGER,
  pid INTEGER,
  host TEXT,
  started_at TEXT,
  finished_at TEXT,
  meta_json TEXT NOT NULL DEFAULT '{}',
  error_kind TEXT,
  error_message TEXT,
  not_before TEXT
);
CREATE INDEX IF NOT EXISTS runs_subject ON runs (subject_type, subject_id, stage, ugc_wgw_version);
CREATE INDEX IF NOT EXISTS runs_status ON runs (status);
CREATE TABLE IF NOT EXISTS events (
  id INTEGER PRIMARY KEY AUTOINCREMENT,
  ts TEXT NOT NULL,
  run_id TEXT,
  level TEXT NOT NULL,
  event TEXT NOT NULL,
  detail_json TEXT NOT NULL DEFAULT '{}'
);
CREATE INDEX IF NOT EXISTS events_run ON events (run_id);
"""


@dataclass
class SampleRecord:
    sample_id: str
    sex: str | None
    father_id: str | None
    mother_id: str | None
    added_at: str
    hifi_reads: list[str] = field(default_factory=list)
    fail_reads: list[str] = field(default_factory=list)
    meta: dict[str, object] = field(default_factory=dict)


@dataclass
class Cohort:
    cohort_id: str
    frozen_at: str
    sample_list_sha256: str
    members: list[str] = field(default_factory=list)


@dataclass
class RunRecord:
    run_id: str
    subject_type: str
    subject_id: str
    stage: str
    mode: str
    ugc_wgw_version: str
    run_dir: str
    status: str
    attempt: int
    cohort_id: str | None = None
    inputs_path: str | None = None
    inputs_sha256: str | None = None
    error_class: str | None = None
    exit_code: int | None = None
    pid: int | None = None
    host: str | None = None
    started_at: str | None = None
    finished_at: str | None = None
    meta: dict[str, object] = field(default_factory=dict)
    error_kind: str | None = None      # failures.KINDS
    error_message: str | None = None   # short classified message
    not_before: str | None = None      # earliest automatic re-attempt (transient failures)

    @property
    def run_path(self) -> Path:
        return Path(self.run_dir)

    @property
    def is_active(self) -> bool:
        return self.status in ACTIVE


def _run_from_row(row: sqlite3.Row) -> RunRecord:
    return RunRecord(
        run_id=row["run_id"], subject_type=row["subject_type"], subject_id=row["subject_id"],
        stage=row["stage"], mode=row["mode"], ugc_wgw_version=row["ugc_wgw_version"], run_dir=row["run_dir"],
        status=row["status"], attempt=row["attempt"], cohort_id=row["cohort_id"],
        inputs_path=row["inputs_path"], inputs_sha256=row["inputs_sha256"], error_class=row["error_class"],
        exit_code=row["exit_code"], pid=row["pid"], host=row["host"], started_at=row["started_at"],
        finished_at=row["finished_at"], meta=json.loads(row["meta_json"] or "{}"),
        error_kind=row["error_kind"], error_message=row["error_message"], not_before=row["not_before"],
    )


class DB:
    def __init__(self, path: Path):
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.conn = sqlite3.connect(str(self.path), isolation_level=None, timeout=30)
        self.conn.row_factory = sqlite3.Row
        self.conn.execute("PRAGMA journal_mode=WAL")
        self.conn.execute("PRAGMA foreign_keys=ON")
        self.conn.execute("PRAGMA busy_timeout=5000")
        # executescript commits any open transaction, so the DDL runs in autocommit mode
        self.conn.executescript(SCHEMA)
        self.migrated: list[int] = []
        try:
            self._migrate()
        except BaseException:
            self.conn.close()
            raise

    def _migrate(self) -> None:
        """Bring an older file to SCHEMA_VERSION (single transaction); refuse a newer one."""
        with self.tx() as c:
            row = c.execute("SELECT version FROM schema_version ORDER BY version DESC LIMIT 1").fetchone()
            if row is None:
                c.execute("INSERT INTO schema_version (version) VALUES (?)", (SCHEMA_VERSION,))
                return
            stored = int(row[0])
            if stored > SCHEMA_VERSION:
                raise UgcError(f"database schema {stored} is newer than this driver ({SCHEMA_VERSION}): {self.path}")
            for version in range(stored + 1, SCHEMA_VERSION + 1):
                for stmt in MIGRATIONS[version]:
                    c.execute(stmt)
                c.execute("DELETE FROM schema_version")
                c.execute("INSERT INTO schema_version (version) VALUES (?)", (version,))
                self.migrated.append(version)

    def close(self) -> None:
        self.conn.close()

    @contextmanager
    def tx(self) -> Iterator[sqlite3.Connection]:
        self.conn.execute("BEGIN IMMEDIATE")
        try:
            yield self.conn
        except BaseException:
            self.conn.execute("ROLLBACK")
            raise
        self.conn.execute("COMMIT")

    # ---- samples ---------------------------------------------------------

    def upsert_sample(self, rec: SampleRecord, replace: bool = False) -> None:
        with self.tx() as c:
            exists = c.execute("SELECT 1 FROM samples WHERE sample_id = ?", (rec.sample_id,)).fetchone()
            if exists and not replace:
                raise UgcError(f"sample {rec.sample_id} already registered (use --replace to overwrite)")
            if exists:
                c.execute("DELETE FROM sample_inputs WHERE sample_id = ?", (rec.sample_id,))
                c.execute("DELETE FROM samples WHERE sample_id = ?", (rec.sample_id,))
            c.execute(
                "INSERT INTO samples (sample_id, sex, father_id, mother_id, added_at, meta_json) VALUES (?, ?, ?, ?, ?, ?)",
                (rec.sample_id, rec.sex, rec.father_id, rec.mother_id, rec.added_at, json.dumps(rec.meta, sort_keys=True)),
            )
            for kind, paths in (("hifi_reads", rec.hifi_reads), ("fail_reads", rec.fail_reads)):
                for i, path in enumerate(paths):
                    c.execute(
                        "INSERT INTO sample_inputs (sample_id, kind, position, path) VALUES (?, ?, ?, ?)",
                        (rec.sample_id, kind, i, path),
                    )

    def get_sample(self, sample_id: str) -> SampleRecord | None:
        row = self.conn.execute("SELECT * FROM samples WHERE sample_id = ?", (sample_id,)).fetchone()
        if row is None:
            return None
        rec = SampleRecord(
            sample_id=row["sample_id"], sex=row["sex"], father_id=row["father_id"], mother_id=row["mother_id"],
            added_at=row["added_at"], meta=json.loads(row["meta_json"] or "{}"),
        )
        for kind in ("hifi_reads", "fail_reads"):
            rows = self.conn.execute(
                "SELECT path FROM sample_inputs WHERE sample_id = ? AND kind = ? ORDER BY position", (sample_id, kind)
            ).fetchall()
            setattr(rec, kind, [r["path"] for r in rows])
        return rec

    def sample_exists(self, sample_id: str) -> bool:
        return self.conn.execute("SELECT 1 FROM samples WHERE sample_id = ?", (sample_id,)).fetchone() is not None

    def list_sample_ids(self) -> list[str]:
        return [r["sample_id"] for r in self.conn.execute("SELECT sample_id FROM samples ORDER BY sample_id")]

    def list_samples(self) -> list[SampleRecord]:
        out = []
        for sid in self.list_sample_ids():
            rec = self.get_sample(sid)
            if rec is not None:
                out.append(rec)
        return out

    def sample_cohorts(self, sample_id: str) -> list[str]:
        rows = self.conn.execute("SELECT cohort_id FROM cohort_members WHERE sample_id = ? ORDER BY cohort_id", (sample_id,))
        return [r["cohort_id"] for r in rows]

    def delete_sample(self, sample_id: str, *, with_runs: bool) -> dict[str, int]:
        """Unregister a sample. Refused for cohort members and for samples with active runs; recorded runs
        need with_runs (their rows go, results on disk stay)."""
        with self.tx() as c:
            if not c.execute("SELECT 1 FROM samples WHERE sample_id = ?", (sample_id,)).fetchone():
                raise UgcError(f"unknown sample {sample_id}")
            cohorts = [r["cohort_id"] for r in c.execute(
                "SELECT cohort_id FROM cohort_members WHERE sample_id = ? ORDER BY cohort_id", (sample_id,))]
            if cohorts:
                raise UgcError(f"sample {sample_id} is a member of cohort(s) {', '.join(cohorts)}; cohorts are immutable")
            active = c.execute(
                "SELECT COUNT(*) FROM runs WHERE subject_type = 'sample' AND subject_id = ? "
                "AND status IN ('pending', 'submitted', 'running')", (sample_id,)).fetchone()[0]
            if active:
                raise UgcError(f"sample {sample_id} has {active} active run(s); stop or reconcile them first")
            n_runs = c.execute("SELECT COUNT(*) FROM runs WHERE subject_type = 'sample' AND subject_id = ?",
                               (sample_id,)).fetchone()[0]
            if n_runs and not with_runs:
                raise UgcError(f"sample {sample_id} has {n_runs} recorded run(s); use --force to delete the rows "
                               "(results on disk are kept)")
            runs_deleted = 0
            if with_runs:
                runs_deleted = c.execute("DELETE FROM runs WHERE subject_type = 'sample' AND subject_id = ?",
                                         (sample_id,)).rowcount
            inputs_deleted = c.execute("DELETE FROM sample_inputs WHERE sample_id = ?", (sample_id,)).rowcount
            c.execute("DELETE FROM samples WHERE sample_id = ?", (sample_id,))
            return {"samples": 1, "sample_inputs": inputs_deleted, "runs": runs_deleted}

    # ---- cohorts ---------------------------------------------------------

    def freeze_cohort(self, cohort_id: str, sample_ids: list[str], sha256: str, frozen_at: str) -> None:
        with self.tx() as c:
            if c.execute("SELECT 1 FROM cohorts WHERE cohort_id = ?", (cohort_id,)).fetchone():
                raise UgcError(f"cohort {cohort_id} is already frozen (cohorts are immutable)")
            unknown = [s for s in sample_ids if not c.execute("SELECT 1 FROM samples WHERE sample_id = ?", (s,)).fetchone()]
            if unknown:
                raise UgcError("unknown samples: " + ", ".join(unknown))
            c.execute(
                "INSERT INTO cohorts (cohort_id, frozen_at, sample_list_sha256) VALUES (?, ?, ?)",
                (cohort_id, frozen_at, sha256),
            )
            for i, sid in enumerate(sample_ids):
                c.execute(
                    "INSERT INTO cohort_members (cohort_id, sample_id, position) VALUES (?, ?, ?)", (cohort_id, sid, i)
                )

    def get_cohort(self, cohort_id: str) -> Cohort | None:
        row = self.conn.execute("SELECT * FROM cohorts WHERE cohort_id = ?", (cohort_id,)).fetchone()
        if row is None:
            return None
        members = [
            r["sample_id"]
            for r in self.conn.execute(
                "SELECT sample_id FROM cohort_members WHERE cohort_id = ? ORDER BY position", (cohort_id,)
            )
        ]
        return Cohort(cohort_id=row["cohort_id"], frozen_at=row["frozen_at"], sample_list_sha256=row["sample_list_sha256"], members=members)

    def list_cohorts(self) -> list[Cohort]:
        ids = [r["cohort_id"] for r in self.conn.execute("SELECT cohort_id FROM cohorts ORDER BY cohort_id")]
        return [c for c in (self.get_cohort(i) for i in ids) if c is not None]

    # ---- runs ------------------------------------------------------------

    def insert_run(self, run: RunRecord) -> None:
        with self.tx() as c:
            c.execute(
                """INSERT INTO runs (run_id, subject_type, subject_id, stage, mode, ugc_wgw_version, cohort_id, inputs_path,
                   inputs_sha256, run_dir, status, attempt, error_class, exit_code, pid, host, started_at, finished_at,
                   meta_json, error_kind, error_message, not_before)
                   VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
                (run.run_id, run.subject_type, run.subject_id, run.stage, run.mode, run.ugc_wgw_version, run.cohort_id,
                 run.inputs_path, run.inputs_sha256, run.run_dir, run.status, run.attempt, run.error_class, run.exit_code,
                 run.pid, run.host, run.started_at, run.finished_at, json.dumps(run.meta, sort_keys=True),
                 run.error_kind, run.error_message, run.not_before),
            )

    def set_status(self, run_id: str, status: str, *, pid: int | None = None, host: str | None = None,
                   started_at: str | None = None) -> None:
        if status not in ACTIVE:
            raise ValueError(f"set_status only moves between active states, not to {status}")
        with self.tx() as c:
            c.execute(
                """UPDATE runs SET status = ?, pid = COALESCE(?, pid), host = COALESCE(?, host),
                   started_at = COALESCE(?, started_at) WHERE run_id = ? AND status IN ('pending', 'submitted', 'running')""",
                (status, pid, host, started_at, run_id),
            )

    def finalize_run(self, run_id: str, status: str, finished_at: str, *, error_class: str | None = None,
                     exit_code: int | None = None, error_kind: str | None = None, error_message: str | None = None,
                     not_before: str | None = None, meta: dict[str, object] | None = None) -> bool:
        """Move an active run to a terminal state in one guarded write. Returns False if it was already
        terminal (idempotent). `meta`, when given, replaces meta_json in the same transition."""
        if status not in TERMINAL:
            raise ValueError(f"finalize_run needs a terminal status, not {status}")
        with self.tx() as c:
            cur = c.execute(
                """UPDATE runs SET status = ?, finished_at = ?, error_class = ?, exit_code = ?, error_kind = ?,
                   error_message = ?, not_before = ?, meta_json = COALESCE(?, meta_json)
                   WHERE run_id = ? AND status IN ('pending', 'submitted', 'running')""",
                (status, finished_at, error_class, exit_code, error_kind, error_message, not_before,
                 json.dumps(meta, sort_keys=True) if meta is not None else None, run_id),
            )
            return cur.rowcount == 1

    def get_run(self, run_id: str) -> RunRecord | None:
        row = self.conn.execute("SELECT * FROM runs WHERE run_id = ?", (run_id,)).fetchone()
        return _run_from_row(row) if row else None

    def _latest(self, subject_type: str, subject_id: str, stage: str, ugc_wgw_version: str | None,
                where_extra: str = "") -> RunRecord | None:
        sql = "SELECT * FROM runs WHERE subject_type = ? AND subject_id = ? AND stage = ?"
        params: list[object] = [subject_type, subject_id, stage]
        if ugc_wgw_version is not None:
            sql += " AND ugc_wgw_version = ?"
            params.append(ugc_wgw_version)
        sql += where_extra + " ORDER BY attempt DESC, rowid DESC LIMIT 1"
        row = self.conn.execute(sql, params).fetchone()
        return _run_from_row(row) if row else None

    def latest_run(self, subject_type: str, subject_id: str, stage: str, ugc_wgw_version: str | None) -> RunRecord | None:
        return self._latest(subject_type, subject_id, stage, ugc_wgw_version)

    def latest_success(self, subject_type: str, subject_id: str, stage: str, ugc_wgw_version: str | None) -> RunRecord | None:
        return self._latest(subject_type, subject_id, stage, ugc_wgw_version, " AND status = 'success'")

    def next_attempt(self, subject_type: str, subject_id: str, stage: str, ugc_wgw_version: str) -> int:
        row = self.conn.execute(
            "SELECT COALESCE(MAX(attempt), 0) FROM runs WHERE subject_type = ? AND subject_id = ? AND stage = ? AND ugc_wgw_version = ?",
            (subject_type, subject_id, stage, ugc_wgw_version),
        ).fetchone()
        return int(row[0]) + 1

    def finished_runs(self, stage: str, ugc_wgw_version: str | None, status: str, limit: int = 20) -> list[RunRecord]:
        """Most recent finished runs of a stage with the given status (newest first)."""
        sql = "SELECT * FROM runs WHERE stage = ? AND status = ? AND finished_at IS NOT NULL"
        params: list[object] = [stage, status]
        if ugc_wgw_version is not None:
            sql += " AND ugc_wgw_version = ?"
            params.append(ugc_wgw_version)
        sql += " ORDER BY finished_at DESC, rowid DESC LIMIT ?"
        params.append(limit)
        return [_run_from_row(r) for r in self.conn.execute(sql, params)]

    def all_runs(self, ugc_wgw_version: str | None = None) -> list[RunRecord]:
        sql, params = "SELECT * FROM runs", []
        if ugc_wgw_version is not None:
            sql, params = sql + " WHERE ugc_wgw_version = ?", [ugc_wgw_version]
        return [_run_from_row(r) for r in self.conn.execute(sql + " ORDER BY rowid", params)]

    def active_runs(self) -> list[RunRecord]:
        rows = self.conn.execute(
            "SELECT * FROM runs WHERE status IN ('pending', 'submitted', 'running') ORDER BY rowid"
        ).fetchall()
        return [_run_from_row(r) for r in rows]

    def runs_for(self, subject_type: str, subject_id: str, stage: str | None = None) -> list[RunRecord]:
        sql = "SELECT * FROM runs WHERE subject_type = ? AND subject_id = ?"
        params: list[object] = [subject_type, subject_id]
        if stage:
            sql += " AND stage = ?"
            params.append(stage)
        sql += " ORDER BY ugc_wgw_version, stage, attempt, rowid"
        return [_run_from_row(r) for r in self.conn.execute(sql, params)]

    def latest_per_stage(self, ugc_wgw_version: str | None) -> dict[tuple[str, str, str], RunRecord]:
        """Latest run (highest attempt) for every (subject_type, subject_id, stage), optionally at one version."""
        sql = "SELECT * FROM runs"
        params: list[object] = []
        if ugc_wgw_version is not None:
            sql += " WHERE ugc_wgw_version = ?"
            params.append(ugc_wgw_version)
        sql += " ORDER BY attempt, rowid"
        out: dict[tuple[str, str, str], RunRecord] = {}
        for row in self.conn.execute(sql, params):
            run = _run_from_row(row)
            out[(run.subject_type, run.subject_id, run.stage)] = run
        return out

    # ---- events ----------------------------------------------------------

    def add_event(self, ts: str, run_id: str | None, level: str, event: str, detail: dict[str, object]) -> None:
        with self.tx() as c:
            c.execute(
                "INSERT INTO events (ts, run_id, level, event, detail_json) VALUES (?, ?, ?, ?, ?)",
                (ts, run_id, level, event, json.dumps(detail, sort_keys=True, default=str)),
            )

    def all_events(self) -> list[dict[str, object]]:
        rows = self.conn.execute("SELECT * FROM events ORDER BY id").fetchall()
        return [{"ts": r["ts"], "run_id": r["run_id"], "level": r["level"], "event": r["event"],
                 "detail": json.loads(r["detail_json"] or "{}")} for r in rows]

    def events_for(self, run_id: str) -> list[dict[str, object]]:
        rows = self.conn.execute("SELECT * FROM events WHERE run_id = ? ORDER BY id", (run_id,)).fetchall()
        return [
            {"ts": r["ts"], "run_id": r["run_id"], "level": r["level"], "event": r["event"],
             "detail": json.loads(r["detail_json"] or "{}")}
            for r in rows
        ]
