import sqlite3
import tempfile
import unittest
from pathlib import Path

from ugc_wgw import db as db_mod
from ugc_wgw.db import DB
from ugc_wgw.util import UgcError, utc_now

V1_RUNS = """
CREATE TABLE schema_version (version INTEGER NOT NULL);
CREATE TABLE samples (sample_id TEXT PRIMARY KEY, sex TEXT, father_id TEXT, mother_id TEXT, added_at TEXT NOT NULL,
  meta_json TEXT NOT NULL DEFAULT '{}');
CREATE TABLE runs (
  run_id TEXT PRIMARY KEY, subject_type TEXT NOT NULL, subject_id TEXT NOT NULL, stage TEXT NOT NULL, mode TEXT NOT NULL,
  ugc_wgw_version TEXT NOT NULL, cohort_id TEXT, inputs_path TEXT, inputs_sha256 TEXT, run_dir TEXT NOT NULL,
  status TEXT NOT NULL, attempt INTEGER NOT NULL, error_class TEXT, exit_code INTEGER, pid INTEGER, host TEXT,
  started_at TEXT, finished_at TEXT, meta_json TEXT NOT NULL DEFAULT '{}');
INSERT INTO schema_version (version) VALUES (1);
INSERT INTO runs (run_id, subject_type, subject_id, stage, mode, ugc_wgw_version, run_dir, status, attempt)
  VALUES ('S1-singleton-a1-old', 'sample', 'S1', 'singleton', 'standalone', '0.1.0', '/r/a1', 'running', 1);
"""


class MigrationTest(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.path = self.tmp / "state.sqlite"

    def columns(self, conn, table):
        return [r[1] for r in conn.execute(f"PRAGMA table_info({table})")]

    def test_v1_database_is_migrated(self):
        conn = sqlite3.connect(str(self.path))
        conn.executescript(V1_RUNS)
        conn.close()
        db = DB(self.path)
        self.assertEqual(db.migrated, [2])
        self.assertEqual([r[0] for r in db.conn.execute("SELECT version FROM schema_version")], [db_mod.SCHEMA_VERSION])
        for col in ("error_kind", "error_message", "not_before"):
            self.assertIn(col, self.columns(db.conn, "runs"))
        run = db.get_run("S1-singleton-a1-old")
        self.assertIsNone(run.error_kind)
        self.assertTrue(db.finalize_run(run.run_id, "failed", utc_now(), error_class="CommandFailed", error_kind="tool",
                                        error_message="m", meta={"slurm_job_ids": ["1"]}))
        run = db.get_run(run.run_id)
        self.assertEqual((run.error_kind, run.error_message, run.meta["slurm_job_ids"]), ("tool", "m", ["1"]))
        db.close()
        self.assertEqual(DB(self.path).migrated, [])  # second open: nothing to do

    def test_newer_schema_refused(self):
        db = DB(self.path)
        db.conn.execute("UPDATE schema_version SET version = 99")
        db.close()
        with self.assertRaisesRegex(UgcError, "newer than this driver"):
            DB(self.path)

    def test_fresh_database_is_current(self):
        db = DB(self.path)
        rows = [r[0] for r in db.conn.execute("SELECT version FROM schema_version")]
        self.assertEqual(rows, [db_mod.SCHEMA_VERSION])
        self.assertEqual(db.migrated, [])
        db.close()


if __name__ == "__main__":
    unittest.main()
