import tempfile
import unittest
from pathlib import Path

from ugc_wgw import layout, manifest
from ugc_wgw import lease as lease_mod
from ugc_wgw.db import DB, RunRecord
from ugc_wgw.log import Events
from ugc_wgw.reconcile import reconcile
from ugc_wgw.submit import acquire_lock
from ugc_wgw.util import UgcError, hostname, utc_now

from .helpers import REPO, VERSION, make_project, run_cli, write_lease, write_tsv


class LeaseTest(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.cfg = make_project(self.tmp, lease_seconds=100)
        self.proj = str(self.cfg.project_dir)

    def test_lock_line_and_refresh(self):
        lock = acquire_lock(self.cfg)
        try:
            lease = lease_mod.read(self.cfg.lock_path)
            self.assertTrue(lease.is_mine())
            self.assertEqual(lease.last_seen, lease.since)
            self.assertTrue(lock.refresh())
            self.assertGreaterEqual(lease_mod.read(self.cfg.lock_path).last_seen, lease.since)
            self.assertFalse(lease.expired(100))
        finally:
            lock.close()
        self.assertIsNone(lease_mod.parse(""))
        old = lease_mod.parse("pid=1 host=h since=2026-01-01T00:00:00Z\n")  # phase-1 line without last_seen
        self.assertTrue(old.expired(100))

    def test_fresh_lease_on_other_host_refused_unless_takeover(self):
        write_lease(self.cfg, host="other-node")
        with self.assertRaisesRegex(UgcError, "holds the lease on other-node"):
            acquire_lock(self.cfg)
        lock = acquire_lock(self.cfg, takeover=True)
        try:
            self.assertTrue(lease_mod.read(self.cfg.lock_path).is_mine())
        finally:
            lock.close()

    def test_expired_lease_is_taken(self):
        write_lease(self.cfg, host="other-node", age_seconds=1000)
        lock = acquire_lock(self.cfg)
        try:
            self.assertTrue(lease_mod.read(self.cfg.lock_path).is_mine())
        finally:
            lock.close()

    def test_refresh_detects_takeover(self):
        lock = acquire_lock(self.cfg)
        try:
            write_lease(self.cfg, host="other-node", pid=42)
            self.assertFalse(lock.refresh())
        finally:
            lock.close()

    def test_submit_takeover_via_cli(self):
        run_cli(["--project", self.proj, "samples", "add", str(write_tsv(self.tmp, [{"sample_id": "S1", "hifi_reads": "a.bam"}]))])
        write_lease(self.cfg, host="other-node")
        code, _, err = run_cli(["--project", self.proj, "submit", "--mode", "standalone", "--samples", "S1"])
        self.assertEqual(code, 1)
        self.assertIn("holds the lease", err)
        code, _, err = run_cli(["--project", self.proj, "submit", "--mode", "standalone", "--samples", "S1", "--takeover"])
        self.assertEqual(code, 0, err)

    def active_run(self, db, sid, host):
        stage_path = layout.stage_dir(self.cfg.results_dir, "sample", sid, VERSION, "singleton")
        path = layout.make_attempt_dir(stage_path / "attempt-1")
        run = RunRecord(run_id=f"{sid}-singleton-a1-t", subject_type="sample", subject_id=sid, stage="singleton",
                        mode="standalone", ugc_wgw_version=VERSION, run_dir=str(path), status="pending", attempt=1,
                        host=host, started_at=utc_now())
        db.insert_run(run)
        db.set_status(run.run_id, "running", pid=12345)
        return run

    def test_reconcile_other_host_lease_rules(self):
        db = DB(self.cfg.db_path)
        events = Events(db, self.cfg.logs_dir / "events.jsonl")
        code = manifest.load_code_info(REPO)
        try:
            run = self.active_run(db, "S1", "other-node")
            write_lease(self.cfg, host="other-node")  # fresh: left alone
            self.assertEqual(reconcile(self.cfg, db, events, code, {}), [])
            write_lease(self.cfg, host="other-node", age_seconds=1000)  # expired: settled
            settled = reconcile(self.cfg, db, events, code, {})
            self.assertEqual([r.run_id for r in settled], [run.run_id])
            r = db.get_run(run.run_id)
            self.assertEqual((r.status, r.error_class), ("failed", "driver_lost"))
            ev = [e for e in db.events_for(run.run_id) if e["event"] == "run.reconciled"]
            self.assertEqual(ev[0]["detail"]["source"], "lease expired")
            run2 = self.active_run(db, "S2", "other-node")
            write_lease(self.cfg, host="other-node")
            lock = acquire_lock(self.cfg, takeover=True)  # holding the lease myself proves the other driver is gone
            try:
                settled = reconcile(self.cfg, db, events, code, {})
            finally:
                lock.close()
            self.assertEqual([r.run_id for r in settled], [run2.run_id])
            self.assertEqual(db.get_run(run2.run_id).status, "failed")
        finally:
            db.close()


if __name__ == "__main__":
    unittest.main()
