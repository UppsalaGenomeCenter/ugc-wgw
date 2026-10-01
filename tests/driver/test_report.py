import json
import tempfile
import unittest
from pathlib import Path

from ugc_wgw import report

from .helpers import HtmlChecker as _Checker
from .helpers import make_project, run_cli, sacct_env, write_ids, write_tsv


class WorkflowLogParseTest(unittest.TestCase):
    def test_pairs_setup_with_done(self):
        tmp = Path(tempfile.mkdtemp())
        log = tmp / "workflow.log.json"
        src = "wdl.w:ugc_wgw_x.w:call-sub.t:call-pbmm2-0"
        src2 = "wdl.w:ugc_wgw_x.t:call-mosdepth"
        lines = [
            {"message": "workflow start", "source": "wdl.w:ugc_wgw_x", "timestamp": 100.0, "level": "NOTICE"},
            {"message": "task setup", "name": "pbmm2_align_wgs", "source": src, "timestamp": 101.0, "level": "NOTICE"},
            {"message": "ignored runtime settings", "keys": ["disk", "zones"], "source": src, "timestamp": 101.5, "level": "WARNING"},
            {"message": "runtime.cpu adjusted to host limit", "original": 32, "adjusted": 24, "source": src, "timestamp": 102.0, "level": "WARNING"},
            {"message": "ugc-wgw resource policy applied", "task": "pbmm2_align_wgs", "summary": "cpu 24→16, partition -→fat",
             "changes": {"cpu": [24, 16], "slurm_partition": [None, "fat"]}, "source": src, "timestamp": 102.5, "level": "NOTICE"},
            {"message": "ugc-wgw gpu request", "task": "pbmm2_align_wgs", "gres": "gpu:a100:1", "partition": "gpu", "source": src,
             "timestamp": 102.6, "level": "NOTICE"},
            {"message": "failed task will be retried", "source": src, "timestamp": 150.0, "level": "ERROR"},
            {"message": "done", "source": src, "timestamp": 161.0, "level": "NOTICE"},
            {"message": "task setup", "name": "mosdepth", "source": src2, "dir": "/x/attempt-1/call-mosdepth", "timestamp": 162.0, "level": "NOTICE"},
            {"message": "done (cached)", "source": src2, "timestamp": 162.5, "level": "NOTICE"},
            {"message": "done", "source": "wdl.w:ugc_wgw_x", "timestamp": 170.0, "level": "NOTICE"},
            {"message": "not json"},
        ]
        log.write_text("\n".join(json.dumps(l) for l in lines[:-1]) + "\nnot json\n")
        tasks, ignored = report.parse_workflow_log(log)
        self.assertEqual([(t.name, t.call_id, t.duration, t.cached, t.retries) for t in tasks],
                         [("pbmm2_align_wgs", "call-pbmm2-0", 60.0, False, 1), ("mosdepth", "call-mosdepth", 0.5, True, 0)])
        self.assertEqual((tasks[0].cpu_requested, tasks[0].cpu_granted), (32, 24))
        self.assertEqual((tasks[0].policy, tasks[1].policy), ("cpu 24→16, partition -→fat", ""))
        self.assertEqual((tasks[0].call_path, tasks[1].call_path), ("call-sub/call-pbmm2-0", "call-mosdepth"))
        self.assertEqual((tasks[0].cpu_policy, tasks[0].cpu_launched, tasks[0].gres), (16, 16, "gpu:a100:1"))
        self.assertEqual((tasks[1].dir, tasks[1].cpu_launched, tasks[1].gres), ("/x/attempt-1/call-mosdepth", None, None))
        self.assertEqual(ignored, {"disk": 1, "zones": 1})
        self.assertEqual(report.parse_workflow_log(tmp / "missing.json"), ([], {}))


class ReportTest(unittest.TestCase):
    def test_report_end_to_end(self):
        tmp = Path(tempfile.mkdtemp())
        cfg = make_project(tmp)
        proj = str(cfg.project_dir)
        run_cli(["--project", proj, "samples", "add", str(write_tsv(tmp, [
            {"sample_id": "S1", "hifi_reads": "s1.bam"}, {"sample_id": "S2", "hifi_reads": "s2.bam"}]))])
        run_cli(["--project", proj, "cohort", "freeze", "C1", "--samples", str(write_ids(tmp, ["S1", "S2"]))])
        env, _ = sacct_env(tmp, "31337|call-fake|COMPLETED|0:0|2026-09-25T10:00:00|2026-09-25T10:00:30|2026-09-25T10:05:30|300|8|8|2400|"
                                "30:00.000|16000M||billing=8,cpu=8,mem=16000M,node=1|core|m1|120||\n")
        env.update({"UGC_WGW_FAKE_FAIL": "S2", "UGC_WGW_FAKE_SLURM_LOG": "31337", "UGC_WGW_FAKE_TASK_LOG": "1"})
        code, _, err = run_cli(["--project", proj, "submit", "--mode", "standalone", "--cohort", "C1"], env=env)
        self.assertEqual(code, 1, err)
        out = tmp / "report.html"
        code, stdout, err = run_cli(["--project", proj, "report", "--mode", "standalone", "--cohort", "C1", "--sizes",
                                     "--price", "cpu_hour=0.05", "--price", "currency=EUR", "--out", str(out)])
        self.assertEqual(code, 0, err)
        self.assertEqual(stdout.strip(), str(out))
        text = out.read_text()
        checker = _Checker()
        checker.feed(text)
        self.assertEqual(checker.errors, [])
        self.assertEqual(checker.stack, [])
        self.assertEqual(checker.titles, 1)
        for needle in ("proj run report", "Provenance", "Stages", "Resource usage", "Timeline", "Failures", "S2", "tool",
                       "CommandFailed", "singleton", "Driver sessions", "31337", "hifi-human-wgs-wdl", "<svg",
                       "SLURM accounting (sacct) for 2 of 2 attempt(s)", "core-hours allocated", "prices: cpu_hour 0.05 EUR",
                       "estimated compute cost", "Core-hours per task", "fake_task", "mean per sample, clean",
                       "results on disk, all attempts", "sacct (ok), 1 jobs"):
            self.assertIn(needle, text, needle)
        self.assertIn("EUR", text)
        self.assertNotIn("src=\"http", text)
        self.assertNotIn("href=\"http", text)
        # default location and any-version
        code, stdout, err = run_cli(["--project", proj, "report", "--any-version"])
        self.assertEqual(code, 0, err)
        self.assertTrue(Path(stdout.strip()).is_relative_to(cfg.results_dir / "reports"))
        code, _, err = run_cli(["--project", proj, "report", "--cohort", "nope"])
        self.assertEqual(code, 1)
        self.assertIn("unknown cohort", err)
        code, _, err = run_cli(["--project", proj, "report", "--price", "nope=1"])
        self.assertEqual(code, 1)
        self.assertIn("unknown price", err)

    def test_report_without_accounting_says_estimate(self):
        tmp = Path(tempfile.mkdtemp())
        cfg = make_project(tmp)
        proj = str(cfg.project_dir)
        run_cli(["--project", proj, "samples", "add", str(write_tsv(tmp, [{"sample_id": "S1", "hifi_reads": "s1.bam"}]))])
        code, _, err = run_cli(["--project", proj, "submit", "--mode", "standalone", "--samples", "S1"],
                               env={"UGC_WGW_FAKE_TASK_LOG": "1", "PATH": "/usr/bin:/bin"})
        self.assertEqual(code, 0, err)
        out = tmp / "report.html"
        code, _, err = run_cli(["--project", proj, "report", "--out", str(out)])
        self.assertEqual(code, 0, err)
        text = out.read_text()
        self.assertIn("SLURM accounting (sacct) for 0 of 1 attempt(s); 1 estimated from workflow.log wall time", text)
        self.assertIn("no prices given", text)
        self.assertIn("estimate", text)
        checker = _Checker()
        checker.feed(text)
        self.assertEqual((checker.errors, checker.stack, checker.titles), ([], [], 1))


if __name__ == "__main__":
    unittest.main()
