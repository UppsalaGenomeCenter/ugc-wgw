import json
import os
import signal
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path

from ugc_wgw import layout, manifest, slurm
from ugc_wgw.db import DB, RunRecord
from ugc_wgw.log import Events
from ugc_wgw.reconcile import reconcile
from ugc_wgw.util import hostname, utc_now

from .helpers import REPO, VERSION, make_project, run_cli, scancel_env, write_tsv


class JobIdsTest(unittest.TestCase):
    def test_job_ids_pruned_walk(self):
        tmp = Path(tempfile.mkdtemp())
        for rel, text in (("call-a", "123;cluster\n"), ("call-sub/call-b", "456\n"), ("call-a/failed1", "789\n"),
                          ("call-a/work/deep", "999\n"), ("out/x", "111\n")):
            d = tmp / rel
            d.mkdir(parents=True)
            (d / slurm.LOG_NAME).write_text(text)
        self.assertEqual(slurm.job_ids(tmp), ["123", "789", "456"])

    def test_scancel_skipped_without_binary(self):
        env = dict(os.environ, PATH="/nonexistent")
        old = os.environ.copy()
        os.environ.clear()
        os.environ.update(env)
        try:
            self.assertEqual(slurm.scancel(["1"])[0], "skipped")
            self.assertEqual(slurm.scancel([])[0], "skipped")
        finally:
            os.environ.clear()
            os.environ.update(old)


class OrphanTest(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.cfg = make_project(self.tmp)
        self.proj = str(self.cfg.project_dir)
        run_cli(["--project", self.proj, "samples", "add", str(write_tsv(self.tmp, [{"sample_id": "S1", "hifi_reads": "a.bam"}]))])

    def test_manifest_and_meta_record_job_ids(self):
        code, _, err = run_cli(["--project", self.proj, "submit", "--mode", "standalone", "--samples", "S1"],
                               env={"UGC_WGW_FAKE_SLURM_LOG": "777"})
        self.assertEqual(code, 0, err)
        db = DB(self.cfg.db_path)
        run = db.latest_run("sample", "S1", "singleton", VERSION)
        events = db.events_for(run.run_id)
        db.close()
        self.assertEqual(run.meta["slurm_job_ids"], ["777"])
        m = json.loads((run.run_path / "run_manifest.json").read_text())
        self.assertEqual(m["slurm_job_ids"], ["777"])
        done = [e for e in events if e["event"] == "run.success"]
        self.assertEqual(done[0]["detail"]["slurm_job_ids"], ["777"])

    def test_reconcile_driver_lost_cancels_orphans(self):
        db = DB(self.cfg.db_path)
        events = Events(db, self.cfg.logs_dir / "events.jsonl")
        code = manifest.load_code_info(REPO)
        stage_path = layout.stage_dir(self.cfg.results_dir, "sample", "S1", VERSION, "singleton")
        path = layout.make_attempt_dir(stage_path / "attempt-1")
        (path / "call-fake").mkdir()
        (path / "call-fake" / slurm.LOG_NAME).write_text("888\n")
        run = RunRecord(run_id="S1-singleton-a1-t", subject_type="sample", subject_id="S1", stage="singleton",
                        mode="standalone", ugc_wgw_version=VERSION, run_dir=str(path), status="pending", attempt=1,
                        host=hostname(), started_at=utc_now())
        db.insert_run(run)
        db.set_status(run.run_id, "running", pid=2**22 + 999)
        env, log = scancel_env(self.tmp)
        old = os.environ.copy()
        os.environ.update(env)
        try:
            reconcile(self.cfg, db, events, code, {})
        finally:
            os.environ.clear()
            os.environ.update(old)
        r = db.get_run(run.run_id)
        self.assertEqual((r.status, r.error_class, r.error_kind), ("failed", "driver_lost", "transient"))
        self.assertEqual(r.meta["slurm_job_ids"], ["888"])
        self.assertEqual(log.read_text().split(), ["888"])
        ev = [e for e in db.events_for(run.run_id) if e["event"] == "run.scancel"]
        self.assertEqual((ev[0]["detail"]["status"], ev[0]["detail"]["job_ids"]), ("ok", ["888"]))
        db.close()

    def test_cancel_orphans_disabled(self):
        self.cfg.cancel_orphans = False
        from ugc_wgw import config as config_mod
        config_mod.save(self.cfg)
        db = DB(self.cfg.db_path)
        events = Events(db, self.cfg.logs_dir / "events.jsonl")
        code = manifest.load_code_info(REPO)
        stage_path = layout.stage_dir(self.cfg.results_dir, "sample", "S1", VERSION, "singleton")
        path = layout.make_attempt_dir(stage_path / "attempt-1")
        (path / "call-fake").mkdir()
        (path / "call-fake" / slurm.LOG_NAME).write_text("889\n")
        run = RunRecord(run_id="S1-singleton-a1-u", subject_type="sample", subject_id="S1", stage="singleton",
                        mode="standalone", ugc_wgw_version=VERSION, run_dir=str(path), status="pending", attempt=1,
                        host=hostname(), started_at=utc_now())
        db.insert_run(run)
        db.set_status(run.run_id, "running", pid=2**22 + 998)
        cfg = config_mod.load(self.cfg.project_dir)
        reconcile(cfg, db, events, code, {})
        ev = [e for e in db.events_for(run.run_id) if e["event"] == "run.scancel"]
        self.assertEqual(ev[0]["detail"]["status"], "skipped")
        db.close()

    def test_sigkill_cancels_orphans(self):
        env, log = scancel_env(self.tmp)
        env = dict(os.environ, **env, UGC_WGW_FAKE_HANG="1", UGC_WGW_FAKE_IGNORE_TERM="1", UGC_WGW_FAKE_SLURM_LOG="4242")
        proc = subprocess.Popen([sys.executable, str(REPO / "bin" / "ugc-wgw"), "--project", self.proj, "submit",
                                 "--mode", "standalone", "--samples", "S1", "--poll-interval", "0.1", "--grace", "1"],
                                env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        wlog = self.cfg.results_dir / "samples" / "S1" / VERSION / "singleton" / "attempt-1" / "workflow.log"
        deadline = time.time() + 20
        while not wlog.exists() and time.time() < deadline:
            time.sleep(0.1)
        self.assertTrue(wlog.exists(), "fake miniwdl never started")
        time.sleep(0.3)
        proc.send_signal(signal.SIGTERM)
        try:
            _, err = proc.communicate(timeout=30)
        except subprocess.TimeoutExpired:
            proc.kill()
            self.fail("driver did not stop")
        self.assertEqual(proc.returncode, 130, err)
        db = DB(self.cfg.db_path)
        run = db.latest_run("sample", "S1", "singleton", VERSION)
        ev = [e for e in db.events_for(run.run_id) if e["event"] == "run.scancel"]
        db.close()
        self.assertEqual((run.status, run.error_class, run.error_kind), ("cancelled", "killed", "cancelled"))
        self.assertEqual(run.meta["slurm_job_ids"], ["4242"])
        self.assertEqual(log.read_text().split(), ["4242"])
        self.assertEqual(ev[0]["detail"]["status"], "ok")


if __name__ == "__main__":
    unittest.main()
