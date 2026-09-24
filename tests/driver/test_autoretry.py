import json
import tempfile
import time
import unittest
from pathlib import Path

from ugc_wgw import plan
from ugc_wgw.db import DB, RunRecord
from ugc_wgw.util import utc_now, utc_plus

from .helpers import VERSION, make_project, run_cli, write_tsv


class AutoRetryTest(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.cfg = make_project(self.tmp, backoff_seconds=1, auto_retry_max=2)
        self.proj = str(self.cfg.project_dir)
        code, _, err = run_cli(["--project", self.proj, "samples", "add",
                                str(write_tsv(self.tmp, [{"sample_id": "S1", "hifi_reads": "s1.bam"}]))])
        self.assertEqual(code, 0, err)

    def events(self):
        return [json.loads(l) for l in (self.cfg.logs_dir / "events.jsonl").read_text().splitlines()]

    def failed_row(self, attempt, kind, not_before=None, status="failed", cls="Interrupted"):
        db = DB(self.cfg.db_path)
        run = RunRecord(run_id=f"S1-singleton-a{attempt}-x", subject_type="sample", subject_id="S1", stage="singleton",
                        mode="standalone", ugc_wgw_version=VERSION, run_dir=str(self.tmp / f"a{attempt}"), status=status,
                        attempt=attempt, error_class=cls, error_kind=kind, not_before=not_before)
        db.insert_run(run)
        db.close()

    def test_transient_failure_is_retried_after_backoff(self):
        t0 = time.monotonic()
        code, _, err = run_cli(["--project", self.proj, "submit", "--mode", "standalone", "--samples", "S1"],
                               env={"UGC_WGW_FAKE_FAIL": "S1", "UGC_WGW_FAKE_FAIL_CLASS": "Interrupted:143",
                                    "UGC_WGW_FAKE_FAIL_ATTEMPTS": "1"})
        self.assertEqual(code, 0, err)  # the subject recovered, so the session is a success
        self.assertLess(time.monotonic() - t0, 30.0)
        db = DB(self.cfg.db_path)
        runs = db.runs_for("sample", "S1", "singleton")
        db.close()
        self.assertEqual([r.attempt for r in runs], [1, 2])
        a1, a2 = runs
        self.assertEqual((a1.status, a1.error_class, a1.error_kind), ("failed", "Interrupted", "transient"))
        self.assertIsNotNone(a1.not_before)
        self.assertGreaterEqual(a1.not_before, a1.finished_at)
        self.assertGreaterEqual(a2.started_at, a1.not_before)  # the second attempt waited out the backoff
        self.assertEqual(a2.status, "success")
        names = [e["event"] for e in self.events()]
        self.assertIn("run.auto_retry", names)
        blocked = [e for e in self.events() if e["event"] == "run.blocked"]
        self.assertTrue(any(str(e["detail"]["reason"]).startswith("backing off until") for e in blocked), blocked)
        code, out, _ = run_cli(["--project", self.proj, "status", "--json"])
        self.assertEqual(json.loads(out)[0]["singleton"], "success(2)")

    def test_tool_failure_is_not_auto_retried(self):
        code, _, err = run_cli(["--project", self.proj, "submit", "--mode", "standalone", "--samples", "S1"],
                               env={"UGC_WGW_FAKE_FAIL": "S1"})
        self.assertEqual(code, 1)
        db = DB(self.cfg.db_path)
        runs = db.runs_for("sample", "S1", "singleton")
        db.close()
        self.assertEqual(len(runs), 1)
        self.assertEqual((runs[0].error_kind, runs[0].not_before), ("tool", None))
        self.assertIn("CommandFailed exit_status=1 node=call-fake", runs[0].error_message)
        code, out, _ = run_cli(["--project", self.proj, "status", "--failed"])
        self.assertIn("kind", out.splitlines()[0])
        self.assertIn("tool", out)
        code, out, _ = run_cli(["--project", self.proj, "status", "--json"])
        self.assertEqual(json.loads(out)[0]["kind"], "tool")
        code, out, _ = run_cli(["--project", self.proj, "status"])
        self.assertNotIn("kind", out.splitlines()[0])
        m = json.loads((runs[0].run_path / "run_manifest.json").read_text())
        self.assertEqual(m["error"]["kind"], "tool")
        self.assertEqual(m["error"]["exit_status"], 1)
        self.assertTrue(m["error"]["task_dir"].endswith("call-fake"))

    def test_compute_backoff_states(self):
        sel = plan.Selection(mode="standalone", stage="singleton", samples=["S1"])
        db = DB(self.cfg.db_path)
        try:
            future = utc_plus(3600)
            self.failed_row(1, "transient", not_before=future)
            p = plan.compute(db, self.cfg, sel, VERSION)
            self.assertEqual(len(p.blocked), 1)
            self.assertTrue(p.blocked[0].reason.startswith("backing off until"))
            self.assertEqual(p.blocked[0].not_before, future)
            self.assertEqual(p.next_wake(), future)
            p = plan.compute(db, self.cfg, sel, VERSION, now=utc_plus(7200))
            self.assertEqual(len(p.runnable), 1)
            self.assertIn("auto-retry after transient Interrupted", p.runnable[0].note)
            p = plan.compute(db, self.cfg, sel, VERSION, retry_failed=True)  # manual retry ignores the backoff
            self.assertEqual(len(p.runnable), 1)
            self.assertIsNone(p.runnable[0].note)
            self.failed_row(3, "transient", not_before=utc_now())
            p = plan.compute(db, self.cfg, sel, VERSION)
            self.assertIn("auto-retries exhausted", p.blocked[0].reason)
            self.failed_row(4, "tool", cls="CommandFailed")
            p = plan.compute(db, self.cfg, sel, VERSION)
            self.assertIn("run `ugc-wgw retry --stage singleton`", p.blocked[0].reason)
            self.assertIn("tool", p.blocked[0].reason)
            self.assertIsNone(p.blocked[0].not_before)
        finally:
            db.close()

    def test_retry_ignores_backoff(self):
        self.failed_row(1, "transient", not_before=utc_plus(3600))
        code, out, err = run_cli(["--project", self.proj, "submit", "--mode", "standalone", "--samples", "S1", "--dry-run"])
        self.assertEqual(code, 0, err)
        self.assertIn("backing off until", out)  # a plain submit would wait out the hour; retry does not
        code, _, err = run_cli(["--project", self.proj, "retry", "--mode", "standalone", "--stage", "singleton",
                                "--samples", "S1"])
        self.assertEqual(code, 0, err)
        db = DB(self.cfg.db_path)
        self.assertEqual(db.latest_run("sample", "S1", "singleton", VERSION).status, "success")
        db.close()


if __name__ == "__main__":
    unittest.main()
