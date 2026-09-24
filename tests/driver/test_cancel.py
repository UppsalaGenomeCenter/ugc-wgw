import os
import signal
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path

from ugc_wgw.db import DB

from .helpers import REPO, VERSION, make_project, run_cli, write_tsv


class CancelTest(unittest.TestCase):
    def test_sigterm_cancels_children(self):
        tmp = Path(tempfile.mkdtemp())
        cfg = make_project(tmp)
        proj = str(cfg.project_dir)
        run_cli(["--project", proj, "samples", "add", str(write_tsv(tmp, [{"sample_id": "S1", "hifi_reads": "a.bam"}]))])
        env = dict(os.environ, UGC_WGW_FAKE_HANG="1")
        proc = subprocess.Popen([sys.executable, str(REPO / "bin" / "ugc-wgw"), "--project", proj, "submit", "--mode", "standalone",
                                 "--samples", "S1", "--poll-interval", "0.1"], env=env,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        log = cfg.results_dir / "samples" / "S1" / VERSION / "singleton" / "attempt-1" / "workflow.log"
        deadline = time.time() + 20
        while not log.exists() and time.time() < deadline:
            time.sleep(0.1)
        self.assertTrue(log.exists(), "fake miniwdl never started")
        time.sleep(0.3)
        proc.send_signal(signal.SIGTERM)
        try:
            _, err = proc.communicate(timeout=20)
        except subprocess.TimeoutExpired:
            proc.kill()
            self.fail("driver did not stop after SIGTERM")
        self.assertEqual(proc.returncode, 130, err)
        db = DB(cfg.db_path)
        run = db.latest_run("sample", "S1", "singleton", VERSION)
        db.close()
        self.assertEqual(run.status, "cancelled")
        self.assertTrue((run.run_path / "run_manifest.json").exists())


if __name__ == "__main__":
    unittest.main()
