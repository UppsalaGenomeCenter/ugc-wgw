import os
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path

from .helpers import REPO, VERSION, make_project, run_cli, write_tsv


class FollowTest(unittest.TestCase):
    def test_follow_until_terminal(self):
        tmp = Path(tempfile.mkdtemp())
        cfg = make_project(tmp)
        proj = str(cfg.project_dir)
        run_cli(["--project", proj, "samples", "add", str(write_tsv(tmp, [{"sample_id": "S1", "hifi_reads": "a.bam"}]))])
        env = dict(os.environ, UGC_WGW_FAKE_SLEEP="2")
        proc = subprocess.Popen([sys.executable, str(REPO / "bin" / "ugc-wgw"), "--project", proj, "submit", "--mode", "standalone",
                                 "--samples", "S1", "--poll-interval", "0.1"], env=env,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        log = cfg.results_dir / "samples" / "S1" / VERSION / "singleton" / "attempt-1" / "workflow.log"
        deadline = time.time() + 20
        while not log.exists() and time.time() < deadline:
            time.sleep(0.1)
        self.assertTrue(log.exists(), "fake miniwdl never started")
        code, out, err = run_cli(["--project", proj, "logs", "S1", "--stage", "singleton", "--follow"])
        proc.communicate(timeout=30)
        self.assertEqual(code, 0, err)
        self.assertIn("fake miniwdl run singleton S1", out)
        self.assertIn("finished: status=success", out)
        self.assertEqual(proc.returncode, 0)
        code, out, err = run_cli(["--project", proj, "logs", "S1", "--stage", "singleton", "--follow", "--json"])
        self.assertEqual(code, 0, err)
        self.assertIn('"message": "fake miniwdl run singleton S1"', out)
        self.assertIn("finished: status=success", out)
        code, out, err = run_cli(["--project", proj, "logs", "S1", "--stage", "singleton", "--json", "--tail", "1"])
        self.assertEqual(code, 0, err)
        self.assertIn("workflow.log.json", out)


if __name__ == "__main__":
    unittest.main()
