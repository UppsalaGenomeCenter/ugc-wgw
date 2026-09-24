import json
import tempfile
import unittest
from pathlib import Path

from ugc_wgw import layout, manifest
from ugc_wgw.db import DB, RunRecord
from ugc_wgw.log import Events
from ugc_wgw.reconcile import reconcile
from ugc_wgw.util import hostname, utc_now

from .helpers import REPO, VERSION, make_project, write_lease


class ReconcileTest(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.cfg = make_project(self.tmp)
        self.db = DB(self.cfg.db_path)
        self.events = Events(self.db, self.cfg.logs_dir / "events.jsonl")
        self.code = manifest.load_code_info(REPO)

    def tearDown(self):
        self.db.close()

    def active_run(self, sid: str, *, host: str | None = None, pid: int = 999999999, status: str = "running") -> RunRecord:
        stage_path = layout.stage_dir(self.cfg.results_dir, "sample", sid, VERSION, "singleton")
        path = layout.make_attempt_dir(stage_path / "attempt-1")
        run = RunRecord(run_id=f"{sid}-singleton-a1-t", subject_type="sample", subject_id=sid, stage="singleton",
                        mode="standalone", ugc_wgw_version=VERSION, run_dir=str(path), status="pending", attempt=1,
                        host=host or hostname(), started_at=utc_now())
        self.db.insert_run(run)
        self.db.set_status(run.run_id, status, pid=pid)
        (path / "inputs.json").write_text("{}\n")
        return run

    def test_success_and_failure_from_run_dir(self):
        ok = self.active_run("S1")
        outputs = {"ugc_wgw_singleton.ugc_wgw_workflow_version": VERSION, "ugc_wgw_singleton.stats_file": "/x"}
        (ok.run_path / "outputs.json").write_text(json.dumps(outputs))
        (ok.run_path / "run.json").write_text(json.dumps({"dir": ok.run_dir, "outputs": outputs}))
        bad = self.active_run("S2")
        (bad.run_path / "error.json").write_text(json.dumps({"error": "RunFailed", "cause": {"error": "CommandFailed", "exit_status": 137}}))
        settled = reconcile(self.cfg, self.db, self.events, self.code, {})
        self.assertEqual({r.run_id for r in settled}, {ok.run_id, bad.run_id})
        self.assertEqual(self.db.get_run(ok.run_id).status, "success")
        self.assertTrue((ok.run_path / "run_manifest.json").exists())
        b = self.db.get_run(bad.run_id)
        self.assertEqual((b.status, b.error_class), ("failed", "CommandFailed"))
        self.assertEqual(b.error_kind, "resource")  # exit_status 137
        self.assertIsNone(b.not_before)
        # idempotent: a second pass touches nothing
        self.assertEqual(reconcile(self.cfg, self.db, self.events, self.code, {}), [])

    def test_dead_pid_same_host_is_driver_lost(self):
        run = self.active_run("S3", pid=2**22 + 12345)
        reconcile(self.cfg, self.db, self.events, self.code, {})
        r = self.db.get_run(run.run_id)
        self.assertEqual((r.status, r.error_class), ("failed", "driver_lost"))
        m = json.loads((r.run_path / "run_manifest.json").read_text())
        self.assertEqual(m["error"]["class"], "driver_lost")

    def test_other_host_left_alone_while_its_lease_is_fresh(self):
        run = self.active_run("S4", host="some-other-node")
        write_lease(self.cfg, host="some-other-node")
        self.assertEqual(reconcile(self.cfg, self.db, self.events, self.code, {}), [])
        self.assertEqual(self.db.get_run(run.run_id).status, "running")

    def test_dry_reports_without_changes(self):
        run = self.active_run("S5", pid=2**22 + 1)
        settled = reconcile(self.cfg, self.db, self.events, self.code, {}, dry=True)
        self.assertEqual([r.run_id for r in settled], [run.run_id])
        self.assertEqual(self.db.get_run(run.run_id).status, "running")


if __name__ == "__main__":
    unittest.main()
