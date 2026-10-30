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
            {"message": "ugc-wgw resource escalation applied", "task": "pbmm2_align_wgs", "summary": "memory 32G→64G",
             "changes": {"memory_reservation": [32 * 1024 ** 3, 64 * 1024 ** 3]}, "source": src, "timestamp": 102.7, "level": "NOTICE"},
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
        self.assertEqual((tasks[0].escalation, tasks[0].escalation_changes, tasks[1].escalation),
                         ("memory 32G→64G", {"memory_reservation": [32 * 1024 ** 3, 64 * 1024 ** 3]}, ""))
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
        cfg.miniwdl_cfg.write_text("[scheduler]\ncontainer_backend = ugc_slurm_singularity\n[ugc_wgw]\nscratch_env = TMPDIR\n"
                                   "scratch_mount = /mnt/scratch\nresources = /site/resources.tsv\n")
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
                       "results on disk, all attempts", "sacct (ok), 1 jobs",
                       "container backend", "ugc_slurm_singularity", "node-local scratch",
                       "TMPDIR (bound at /mnt/scratch, TMPDIR inside every task)", "resource policy", "/site/resources.tsv"):
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
        self.assertIn("node-local scratch", text)
        self.assertIn("off (temporaries in the task", text)   # the plain test cfg has no scratch_env
        self.assertIn("(miniwdl default)", text)              # and no container_backend line


if __name__ == "__main__":
    unittest.main()


class ScaleTest(unittest.TestCase):
    """The report at hundreds of runs: lanes instead of rows, per-pixel downsampling, folded tables."""

    def test_per_pixel_keeps_the_shape(self):
        pts = [(float(t), (t * 7919) % 13) for t in range(0, 20000)]
        out = report._per_pixel(pts, 0.0, 20000.0, 500.0)
        self.assertLess(len(out), 4 * 500 + 10)
        self.assertEqual(out[0], pts[0]); self.assertEqual(out[-1], pts[-1])
        self.assertEqual(max(v for _, v in out), 12); self.assertEqual(min(v for _, v in out), 0)
        self.assertEqual(report._per_pixel([], 0.0, 1.0, 10.0), [])

    def test_lanes_above_the_row_limit(self):
        from types import SimpleNamespace
        now = 1_700_000_000
        runs = []
        for i in range(report.TIMELINE_ROWS + 50):
            stage = "upstream" if i % 3 else "downstream"
            status = "failed" if i % 97 == 0 else "success"
            s = now + i * 60; e = s + 3600
            run = SimpleNamespace(subject_id=f"S{i:04d}", stage=stage, attempt=1, status=status,
                                  started_at=report._iso(s) if hasattr(report, "_iso") else __import__("datetime").datetime.fromtimestamp(s, __import__("datetime").timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
                                  finished_at=__import__("datetime").datetime.fromtimestamp(e, __import__("datetime").timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"))
            runs.append(SimpleNamespace(run=run, duration=3600.0, tasks=[]))
        svg = report.timeline_svg(runs)
        self.assertIn("<polygon", svg)
        self.assertEqual(svg.count("<polygon"), 2)                 # one lane per stage
        self.assertEqual(svg.count("<line") - svg.count("stroke=\"#ddd\""), 3)   # the three failures as ticks (axis lines excluded)
        self.assertLess(len(svg), 60_000)
        height = int(__import__("re").search(r'height="(\d+)"', svg).group(1))
        self.assertLess(height, 400)
        small = report.timeline_svg(runs[:10])
        self.assertNotIn("<polygon", small); self.assertEqual(small.count("<rect"), 10)

    def test_big_table_folds_and_filters(self):
        rows = [[f"S{i}", "x"] for i in range(5)]
        self.assertNotIn("<details", report.big_table(["a", "b"], rows, what="rows"))
        rows = [[f"S{i}", "x"] for i in range(report.BIG_TABLE_FROM)]
        html = report.big_table(["a", "b"], rows, what="100 rows (3 bad)")
        self.assertIn("<details><summary>100 rows (3 bad): click to open", html)
        self.assertIn('class="filter"', html)
        self.assertIn("input.filter", report.JS)

