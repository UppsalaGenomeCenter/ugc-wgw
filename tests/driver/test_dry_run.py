import tempfile
import unittest
from pathlib import Path

from ugc_wgw.db import DB

from .helpers import VERSION, make_project, run_cli, write_ids, write_tsv


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


if __name__ == "__main__":
    unittest.main()
