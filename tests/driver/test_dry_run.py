import json
import tempfile
import unittest
from pathlib import Path

from ugc_wgw.db import DB

from .helpers import VERSION, make_project, run_cli, write_bam, write_ids, write_tsv


class DryRunTest(unittest.TestCase):
    def test_dry_run_prints_commands_and_changes_nothing(self):
        tmp = Path(tempfile.mkdtemp())
        cfg = make_project(tmp)
        proj = str(cfg.project_dir)
        run_cli(["--project", proj, "samples", "add", str(write_tsv(tmp, [{"sample_id": "S1", "hifi_reads": "a.bam"},
                                                                          {"sample_id": "S2", "hifi_reads": "b.bam"}]))])
        run_cli(["--project", proj, "cohort", "freeze", "C1", "--samples", str(write_ids(tmp, ["S1", "S2"]))])
        code, out, err = run_cli(["--project", proj, "submit", "--mode", "standalone", "--cohort", "C1", "--dry-run"])
        self.assertEqual(code, 0, err)
        self.assertIn("runnable=2 blocked=2", out)
        self.assertEqual(out.count("miniwdl run"), 0)  # the fake is not called miniwdl; check the command shape instead
        cmds = [l for l in out.splitlines() if " run " in l and "ugc_wgw_singleton.wdl" in l]
        self.assertEqual(len(cmds), 2)
        self.assertIn("--cfg", cmds[0])
        self.assertIn(f"/{VERSION}/singleton/attempt-1/.", cmds[0])
        self.assertIn("-o", cmds[0])
        self.assertIn('"ugc_wgw_singleton.sample_id": "S1"', out)
        self.assertEqual(out.count("# inputs: 1 read file(s), 0.0 GB, present and complete"), 2)
        self.assertIn("blocked: cohort C1 / cohort_merge: waiting for singleton", out)
        db = DB(cfg.db_path)
        self.assertEqual(db.active_runs(), [])
        self.assertIsNone(db.latest_run("sample", "S1", "singleton", VERSION))
        db.close()
        self.assertFalse((cfg.results_dir / "samples").exists())
        # inputs command
        code, out, err = run_cli(["--project", proj, "inputs", "S1", "--stage", "singleton"])
        self.assertEqual(code, 0, err)
        self.assertIn('"ugc_wgw_singleton.hifi_reads"', out)
        code, out, err = run_cli(["--project", proj, "inputs", "C1", "--stage", "cohort_merge"])
        self.assertEqual(code, 1)
        self.assertIn("no output", err)
        code, out, err = run_cli(["--project", proj, "submit", "--mode", "assembly", "--dry-run"])
        self.assertEqual(code, 0, err)
        self.assertIn("runnable=2 blocked=0", out)
        self.assertEqual(len([l for l in out.splitlines() if " run " in l and "ugc_wgw_assembly.wdl" in l]), 2)

    def test_preflight_fails_the_run_without_a_launch(self):
        tmp = Path(tempfile.mkdtemp())
        cfg = make_project(tmp)
        proj = str(cfg.project_dir)
        run_cli(["--project", proj, "samples", "add", str(write_tsv(tmp, [{"sample_id": "S1", "hifi_reads": "a.bam,b.bam"}]))])
        db = DB(cfg.db_path)
        a, b = db.get_sample("S1").hifi_reads
        good = Path(b).read_bytes()
        write_bam(Path(b), reads=300, truncated=True)
        code, out, err = run_cli(["--project", proj, "submit", "--mode", "standalone", "--samples", "S1", "--dry-run"])
        self.assertEqual(code, 0, err)
        self.assertIn(f"# input problem: {b}: size changed since registration", out)
        self.assertNotIn("# inputs:", out)
        trace = tmp / "trace.txt"
        code, out, err = run_cli(["--project", proj, "submit", "--mode", "standalone", "--samples", "S1"],
                                 env={"UGC_WGW_FAKE_TRACE": str(trace)})
        self.assertEqual(code, 1)
        self.assertFalse(trace.exists())   # miniwdl was never launched
        run = db.latest_run("sample", "S1", "singleton", VERSION)
        self.assertEqual((run.status, run.error_class, run.error_kind, run.attempt, run.pid), ("failed", "InputError", "input", 1, None))
        self.assertIn("1 input file problem(s) before launch", run.error_message)
        self.assertIn("size changed since registration", run.error_message)
        self.assertEqual(run.meta["slurm_job_ids"], [])
        self.assertTrue((run.run_path / "run_manifest.json").exists())
        events = [json.loads(l) for l in (cfg.logs_dir / "events.jsonl").read_text().splitlines() if '"run_id": "' + run.run_id in l]
        self.assertEqual([e["event"] for e in events], ["run.created", "run.failed"])
        code, out, _ = run_cli(["--project", proj, "status", "--failed", "--json"])
        self.assertEqual(json.loads(out)[0]["singleton"], "failed")
        # a second submit does not retry an input failure on its own
        code, out, err = run_cli(["--project", proj, "submit", "--mode", "standalone", "--samples", "S1", "--dry-run"])
        self.assertIn("blocked: sample S1 / singleton: failed (attempt 1, InputError, input)", out)
        # the data is restored (same size): retry runs it
        Path(b).write_bytes(good)
        code, out, err = run_cli(["--project", proj, "retry", "--mode", "standalone", "--stage", "singleton", "--samples", "S1"])
        self.assertEqual(code, 0, err)
        run = db.latest_run("sample", "S1", "singleton", VERSION)
        self.assertEqual((run.status, run.attempt), ("success", 2))
        events = [json.loads(l) for l in (cfg.logs_dir / "events.jsonl").read_text().splitlines() if '"run_id": "' + run.run_id in l]
        self.assertEqual([e["event"] for e in events][:3], ["run.created", "run.preflight", "run.submitted"])
        self.assertEqual(events[1]["detail"]["files"], 2)
        db.close()


if __name__ == "__main__":
    unittest.main()
