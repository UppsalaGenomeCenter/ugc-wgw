import json
import tempfile
import unittest
from pathlib import Path

from ugc_wgw import layout, plan, progress
from ugc_wgw.db import DB, RunRecord
from ugc_wgw.util import utc_now, utc_plus

from .helpers import VERSION, make_project, run_cli, seed_success, write_ids, write_tsv


def task_log(path: Path, n_done: int, n_total: int) -> None:
    lines = []
    for i in range(n_total):
        src = f"wdl.w:ugc_wgw_singleton.t:call-t{i}"
        lines.append({"message": "task setup", "name": f"t{i}", "source": src, "timestamp": 100.0 + i, "level": "NOTICE"})
        if i < n_done:
            lines.append({"message": "done", "source": src, "timestamp": 200.0 + i, "level": "NOTICE"})
    path.write_text("\n".join(json.dumps(l) for l in lines) + "\n")


class SummarizeTest(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.cfg = make_project(self.tmp, max_inflight=4)
        self.proj = str(self.cfg.project_dir)
        run_cli(["--project", self.proj, "samples", "add", str(write_tsv(self.tmp, [
            {"sample_id": f"S{i}", "hifi_reads": f"s{i}.bam"} for i in range(1, 6)]))])
        run_cli(["--project", self.proj, "cohort", "freeze", "C1", "--samples", str(write_ids(self.tmp, [f"S{i}" for i in range(1, 6)]))])

    def test_counts_fractions_and_eta(self):
        db = DB(self.cfg.db_path)
        try:
            for sid in ("S1", "S2"):  # two successes of 100 s each, with 4-task logs
                run = seed_success(self.cfg, db, "singleton", "sample", sid, {"sample_id": sid, "hifi_reads": ["/x"]})
                db.conn.execute("UPDATE runs SET started_at = ?, finished_at = ? WHERE run_id = ?",
                                (utc_plus(-100), utc_now(), run.run_id))
                task_log(layout.run_files(run.run_path).workflow_log_json, 4, 4)
            stage_path = layout.stage_dir(self.cfg.results_dir, "sample", "S3", VERSION, "singleton")
            path = layout.make_attempt_dir(stage_path / "attempt-1")
            task_log(layout.run_files(path).workflow_log_json, 2, 4)
            active = RunRecord(run_id="S3-singleton-a1-x", subject_type="sample", subject_id="S3", stage="singleton",
                               mode="standalone", ugc_wgw_version=VERSION, run_dir=str(path), status="running", attempt=1,
                               started_at=utc_now())
            db.insert_run(active)
            db.insert_run(RunRecord(run_id="S4-singleton-a1-x", subject_type="sample", subject_id="S4", stage="singleton",
                                    mode="standalone", ugc_wgw_version=VERSION, run_dir=str(self.tmp / "s4"), status="failed",
                                    attempt=1, error_class="CommandFailed", error_kind="tool"))
            sel = plan.Selection(mode="standalone", cohort="C1")
            p = plan.compute(db, self.cfg, sel, VERSION)
            summary = progress.summarize(db, self.cfg, p, VERSION, mode="standalone", max_inflight=4)
            by = {s.stage: s for s in summary.stages}
            s = by["singleton"]
            self.assertEqual((s.total, s.done, s.runnable, s.failed, s.waiting, len(s.active)), (5, 2, 1, 1, 0, 1))
            self.assertEqual((s.active[0].tasks_done, s.active[0].tasks_expected), (2, 4))
            self.assertAlmostEqual(s.median, 100.0)
            # remaining 2 of which one is half done: ceil(2/4)=1 wave minus 0.5/4 credit -> 87.5 s
            self.assertAlmostEqual(s.eta, 87.5)
            c = by["cohort_merge"]
            self.assertEqual((c.total, c.waiting, c.median, c.eta), (1, 1, None, None))
            self.assertEqual((summary.done, summary.total), (2, 7))
            text = summary.format()
            self.assertIn("progress mode=standalone 2/7 (28%) eta≈1m28s", text)
            self.assertIn("active 1 (S3 2/4 tasks)", text)
            self.assertIn("failed 1 (needs retry)", text)
            self.assertIn("waiting 1", text)
            detail = summary.as_detail()
            self.assertEqual(detail["stages"]["singleton"]["eta_seconds"], 87)
        finally:
            db.close()

    def test_progress_command_and_events(self):
        self.cfg.progress_interval = 0.01
        from ugc_wgw import config as config_mod
        config_mod.save(self.cfg)
        code, _, err = run_cli(["--project", self.proj, "submit", "--mode", "standalone", "--samples", "S1", "S2"],
                               env={"UGC_WGW_FAKE_SLEEP": "0.3"})
        self.assertEqual(code, 0, err)
        events = [json.loads(l) for l in (self.cfg.logs_dir / "events.jsonl").read_text().splitlines()]
        prog = [e for e in events if e["event"] == "submit.progress"]
        self.assertGreaterEqual(len(prog), 2)
        last = prog[-1]["detail"]
        self.assertEqual((last["done"], last["total"]), (2, 4))  # cohort_merge and cohort_freq count as waiting for --cohort
        self.assertEqual((last["stages"]["singleton"]["done"], last["stages"]["singleton"]["total"]), (2, 2))
        self.assertEqual(last["stages"]["cohort_merge"]["waiting"], 1)
        self.assertIn("progress mode=standalone", (self.cfg.logs_dir / "ugc-wgw.log").read_text())
        code, out, err = run_cli(["--project", self.proj, "progress", "--mode", "standalone", "--cohort", "C1"])
        self.assertEqual(code, 0, err)
        self.assertIn("progress mode=standalone 2/7", out)
        code, out, err = run_cli(["--project", self.proj, "progress", "--mode", "standalone", "--cohort", "C1", "--json"])
        self.assertEqual(json.loads(out)["stages"]["singleton"]["done"], 2)


if __name__ == "__main__":
    unittest.main()
