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


class EscalationRetryTest(unittest.TestCase):
    """A task that ran out of memory or time is re-attempted at once with a bigger request, and the project learns."""

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.cfg = make_project(self.tmp, backoff_seconds=1, auto_retry_max=2)
        self.proj = str(self.cfg.project_dir)
        code, _, err = run_cli(["--project", self.proj, "samples", "add",
                                str(write_tsv(self.tmp, [{"sample_id": "S1", "hifi_reads": "s1.bam"},
                                                          {"sample_id": "S2", "hifi_reads": "s2.bam"}]))])
        self.assertEqual(code, 0, err)

    def events(self):
        return [json.loads(l) for l in (self.cfg.logs_dir / "events.jsonl").read_text().splitlines()]

    def test_memory_failure_is_escalated_and_learned(self):
        code, _, err = run_cli(["--project", self.proj, "submit", "--mode", "standalone", "--samples", "S1"],
                               env={"UGC_WGW_FAKE_FAIL": "S1", "UGC_WGW_FAKE_FAIL_RESOURCE": "memory"})
        self.assertEqual(code, 0, err)
        db = DB(self.cfg.db_path)
        a1, a2 = db.runs_for("sample", "S1", "singleton")
        db.close()
        self.assertEqual((a1.status, a1.error_kind, a1.exit_code), ("failed", "resource", 137))
        self.assertIsNone(a1.not_before)   # no backoff: the cause is known
        self.assertEqual(a1.meta["error_evidence"]["resource"], "memory")
        self.assertEqual(a1.meta["escalation"], {"task": "fake_task", "resource": "memory", "factor": 2.0, "n": 1, "max": 2})
        self.assertIn("next attempt: fake_task memory ×2 (escalation 1 of 2)", a1.error_message)
        self.assertEqual((a2.status, a2.attempt), ("success", 2))
        self.assertEqual(a2.meta["escalations"], {"fake_task": {"memory_factor": 2.0}})
        m1 = json.loads((a1.run_path / "run_manifest.json").read_text())
        m2 = json.loads((a2.run_path / "run_manifest.json").read_text())
        self.assertEqual(m1["escalation_next"]["task"], "fake_task")
        self.assertEqual(m2["escalations"], {"fake_task": {"memory_factor": 2.0}})
        ev = self.events()
        submitted = [e for e in ev if e["event"] == "run.submitted" and e["detail"].get("stage") == "singleton"]
        self.assertEqual(len(submitted), 2)
        self.assertNotIn("env", submitted[0]["detail"])
        self.assertEqual(json.loads(submitted[1]["detail"]["env"]["MINIWDL__UGC_WGW__ESCALATIONS"]), {"fake_task": {"memory_factor": 2.0}})
        retries = [e for e in ev if e["event"] == "run.auto_retry"]
        self.assertEqual(len(retries), 1)
        self.assertTrue(retries[0]["detail"]["reason"].startswith("escalation: fake_task memory ×2"))
        # the project learned what the task needed
        learned = json.loads(self.cfg.learned_path.read_text())["tasks"]
        self.assertEqual(learned["fake_task"]["memory"], 16 * 1024 ** 3)
        self.assertEqual(learned["fake_task"]["run_id"], a2.run_id)
        self.assertTrue(any(e["event"] == "resources.learned" and e["detail"]["lessons"] == ["fake_task memory 17179869184"] for e in ev))
        # the next sample starts at the learned request: the fake would fail at the base, it succeeds at once
        code, _, err = run_cli(["--project", self.proj, "submit", "--mode", "standalone", "--samples", "S2"],
                               env={"UGC_WGW_FAKE_FAIL": "S2", "UGC_WGW_FAKE_FAIL_RESOURCE": "memory"})
        self.assertEqual(code, 0, err)
        db = DB(self.cfg.db_path)
        s2 = db.runs_for("sample", "S2", "singleton")
        db.close()
        self.assertEqual([(r.attempt, r.status) for r in s2], [(1, "success")])
        env = [e for e in self.events() if e["event"] == "run.submitted" and e["detail"].get("subject") == "S2"][0]["detail"]["env"]
        self.assertEqual(json.loads(env["MINIWDL__UGC_WGW__ESCALATIONS"]), {"fake_task": {"memory_min": 16 * 1024 ** 3}})
        code, out, _ = run_cli(["--project", self.proj, "status", "--json"])
        self.assertEqual({r["subject"]: r["singleton"] for r in json.loads(out)}, {"S1": "success(2)", "S2": "success"})
        code, out, _ = run_cli(["--project", self.proj, "resources", "--learned"])
        self.assertEqual(code, 0)
        self.assertIn("fake_task\t-\t16G\t-\t-\t-", out)
        code, out, _ = run_cli(["--project", self.proj, "report", "--out", str(self.tmp / "r.html")])
        self.assertEqual(code, 0)
        html = (self.tmp / "r.html").read_text()
        self.assertIn("escalation", html)
        self.assertIn("resource escalations", html)

    def test_time_failure_is_escalated(self):
        code, _, err = run_cli(["--project", self.proj, "submit", "--mode", "standalone", "--samples", "S1"],
                               env={"UGC_WGW_FAKE_FAIL": "S1", "UGC_WGW_FAKE_FAIL_RESOURCE": "time"})
        self.assertEqual(code, 0, err)
        db = DB(self.cfg.db_path)
        a1, a2 = db.runs_for("sample", "S1", "singleton")
        db.close()
        self.assertEqual((a1.error_kind, a1.meta["error_evidence"]["resource"], a1.exit_code), ("resource", "time", 143))
        self.assertEqual(a1.meta["escalation"]["resource"], "time")
        self.assertEqual(a2.status, "success")
        self.assertEqual(json.loads(self.cfg.learned_path.read_text())["tasks"]["fake_task"]["time_minutes"], 240)

    def test_escalation_exhausts_then_blocks(self):
        code, _, err = run_cli(["--project", self.proj, "submit", "--mode", "standalone", "--samples", "S1"],
                               env={"UGC_WGW_FAKE_FAIL": "S1", "UGC_WGW_FAKE_FAIL_RESOURCE": "memory",
                                    "UGC_WGW_FAKE_FAIL_RESOURCE_INCURABLE": "1"})
        self.assertEqual(code, 1)
        db = DB(self.cfg.db_path)
        runs = db.runs_for("sample", "S1", "singleton")
        db.close()
        self.assertEqual([(r.attempt, r.status, r.error_kind) for r in runs], [(1, "failed", "resource"), (2, "failed", "resource"), (3, "failed", "resource")])
        self.assertEqual(runs[1].meta["escalations"], {"fake_task": {"memory_factor": 2.0}})
        self.assertEqual(runs[2].meta["escalations"], {"fake_task": {"memory_factor": 4.0}})
        self.assertNotIn("escalation", runs[2].meta)
        self.assertIn("escalation exhausted: fake_task memory already ×4 (2 of 2)", runs[2].error_message)
        self.assertIn("ugc-wgw retry", err)
        self.assertFalse(self.cfg.learned_path.exists())   # nothing finished with an escalation
        # a manual retry carries the last escalation on (×4), no further doubling
        code, _, err = run_cli(["--project", self.proj, "retry", "--mode", "standalone", "--stage", "singleton", "--samples", "S1",
                                "--max-attempts", "5"], env={"UGC_WGW_FAKE_FAIL": "S1", "UGC_WGW_FAKE_FAIL_RESOURCE": "memory"})
        self.assertEqual(code, 0, err)
        db = DB(self.cfg.db_path)
        a4 = db.runs_for("sample", "S1", "singleton")[-1]
        db.close()
        self.assertEqual((a4.attempt, a4.status, a4.meta["escalations"]), (4, "success", {"fake_task": {"memory_factor": 4.0}}))

    def test_resource_retry_max_zero_blocks_as_before(self):
        self.cfg.resource_retry_max = 0
        from ugc_wgw import config as config_mod
        config_mod.save(self.cfg)
        code, _, err = run_cli(["--project", self.proj, "submit", "--mode", "standalone", "--samples", "S1"],
                               env={"UGC_WGW_FAKE_FAIL": "S1", "UGC_WGW_FAKE_FAIL_RESOURCE": "memory"})
        self.assertEqual(code, 1)
        db = DB(self.cfg.db_path)
        runs = db.runs_for("sample", "S1", "singleton")
        db.close()
        self.assertEqual(len(runs), 1)
        self.assertIn("not escalated: resource_retry_max is 0", runs[0].error_message)
        self.assertIn("ugc-wgw retry", err)

