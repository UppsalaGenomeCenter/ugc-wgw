import json
import os
import tempfile
import unittest
from pathlib import Path

from ugc_wgw import layout
from ugc_wgw.db import DB

from .helpers import REPO, VERSION, make_project, run_cli, write_ids, write_tsv


class SubmitTest(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.cfg = make_project(self.tmp)
        self.proj = str(self.cfg.project_dir)
        code, _, err = run_cli(["--project", self.proj, "samples", "add", str(write_tsv(self.tmp, [
            {"sample_id": "S1", "sex": "MALE", "hifi_reads": "s1.bam"},
            {"sample_id": "S2", "hifi_reads": "s2.bam"},
            {"sample_id": "S3", "hifi_reads": "s3.bam"},
        ]))])
        self.assertEqual(code, 0, err)
        code, _, err = run_cli(["--project", self.proj, "cohort", "freeze", "C1", "--samples", str(write_ids(self.tmp, ["S1", "S2", "S3"]))])
        self.assertEqual(code, 0, err)

    def db(self) -> DB:
        return DB(self.cfg.db_path)

    def events(self) -> list[dict]:
        return [json.loads(l) for l in (self.cfg.logs_dir / "events.jsonl").read_text().splitlines()]

    def test_standalone_end_to_end(self):
        trace = self.tmp / "trace.txt"
        code, out, err = run_cli(["--project", self.proj, "submit", "--mode", "standalone", "--cohort", "C1", "--max-inflight", "2"],
                                 env={"UGC_WGW_FAKE_SLEEP": "0.3", "UGC_WGW_FAKE_TRACE": str(trace)})
        self.assertEqual(code, 0, err)
        db = self.db()
        for sid in ("S1", "S2", "S3"):
            run = db.latest_run("sample", sid, "singleton", VERSION)
            self.assertEqual(run.status, "success")
            self.assertEqual(run.attempt, 1)
            self.assertTrue((run.run_path / "run_manifest.json").exists())
            self.assertEqual(os.readlink(run.run_path.parent / "current"), "attempt-1")
        merge = db.latest_run("cohort", "C1", "cohort_merge", VERSION)
        self.assertEqual(merge.status, "success")
        freq = db.latest_run("cohort", "C1", "cohort_freq", VERSION)
        self.assertEqual(freq.status, "success")
        self.assertGreaterEqual(freq.started_at, merge.finished_at)
        fm = json.loads((freq.run_path / "run_manifest.json").read_text())
        self.assertEqual({m["stage"] for m in fm["cohort_members"]} >= {"cohort_merge", "singleton"}, True)
        self.assertTrue((freq.run_path / "out" / "sv_freq_vcf").exists())
        # transitions in order for one run
        names = [e["event"] for e in self.events() if e["run_id"] == merge.run_id]
        self.assertEqual(names, ["run.created", "run.submitted", "run.running", "run.success"])
        # cohort_merge started only after all singletons ended
        lines = [l.split() for l in trace.read_text().splitlines()]
        ends = [float(l[2]) for l in lines if l[0] == "end" and l[1] != "C1"]
        start_c1 = [float(l[2]) for l in lines if l[0] == "start" and l[1] == "C1"][0]
        self.assertGreaterEqual(start_c1, max(ends))
        # at most 2 singletons in flight
        overlap, running = 0, 0
        for l in sorted(lines, key=lambda x: float(x[2])):
            if l[1] == "C1":
                continue
            running += 1 if l[0] == "start" else -1
            overlap = max(overlap, running)
        self.assertLessEqual(overlap, 2)
        # manifest content
        m = json.loads((merge.run_path / "run_manifest.json").read_text())
        self.assertEqual(m["schema"], 1)
        self.assertEqual(m["ugc_pacbio_wgw"]["version"], VERSION)
        self.assertEqual(m["stage"], "cohort_merge")
        self.assertEqual(m["mode"], "standalone")
        self.assertEqual(m["subject"], {"type": "cohort", "id": "C1"})
        self.assertEqual([x["sample_id"] for x in m["cohort_members"]], ["S1", "S2", "S3"])
        self.assertIn("hifi-human-wgs-wdl", m["upstream"])
        self.assertIn("hifi-wdl-resources", m["references"])
        self.assertTrue(any("svx@sha256" in c["image"] for c in m["containers"]))
        self.assertEqual(m["slurm_job_ids"], [])
        self.assertEqual(m["engine"]["miniwdl"], "miniwdl v0.0.0-fake")
        self.assertEqual(m["ugc_wgw_manifest"]["subject"], {"type": "cohort", "id": "C1"})
        from ugc_wgw.util import sha256_file
        self.assertEqual(m["inputs_sha256"], sha256_file(merge.run_path / "inputs.json"))
        self.assertEqual(m["status"], "success")
        # idempotent re-run
        n_runs = len(db.runs_for("cohort", "C1"))
        db.close()
        code, out, err = run_cli(["--project", self.proj, "submit", "--mode", "standalone", "--cohort", "C1"])
        self.assertEqual(code, 0, err)
        self.assertIn("nothing runnable", err)
        self.assertEqual(len(self.db().runs_for("cohort", "C1")), n_runs)
        # status
        code, out, _ = run_cli(["--project", self.proj, "status", "--json"])
        rows = json.loads(out)
        self.assertEqual({r["subject"]: r.get("cohort_merge", r.get("singleton")) for r in rows},
                         {"S1": "success", "S2": "success", "S3": "success", "C1": "success"})

    def test_failure_blocks_cohort_then_retry(self):
        code, _, err = run_cli(["--project", self.proj, "submit", "--mode", "standalone", "--cohort", "C1"],
                               env={"UGC_WGW_FAKE_FAIL": "S2"})
        self.assertEqual(code, 1)
        db = self.db()
        s2 = db.latest_run("sample", "S2", "singleton", VERSION)
        self.assertEqual(s2.status, "failed")
        self.assertEqual(s2.error_class, "CommandFailed")
        self.assertEqual(s2.exit_code, 1)
        self.assertIsNone(db.latest_run("cohort", "C1", "cohort_merge", VERSION))
        m = json.loads((s2.run_path / "run_manifest.json").read_text())
        self.assertEqual(m["status"], "failed")
        self.assertEqual(m["error"]["class"], "CommandFailed")
        self.assertTrue(any(e["event"] == "run.blocked" and e["detail"]["stage"] == "cohort_merge" for e in self.events()))
        db.close()
        # submit again: still blocked, no new attempt
        code, _, err = run_cli(["--project", self.proj, "submit", "--mode", "standalone", "--cohort", "C1"])
        self.assertEqual(code, 0)
        self.assertIn("ugc-wgw retry", err)
        self.assertEqual(self.db().latest_run("sample", "S2", "singleton", VERSION).attempt, 1)
        # retry creates attempt-2, then cohort_merge can run
        code, _, err = run_cli(["--project", self.proj, "retry", "--mode", "standalone", "--stage", "singleton", "--cohort", "C1"])
        self.assertEqual(code, 0, err)
        db = self.db()
        s2 = db.latest_run("sample", "S2", "singleton", VERSION)
        self.assertEqual((s2.status, s2.attempt), ("success", 2))
        self.assertEqual(os.readlink(s2.run_path.parent / "current"), "attempt-2")
        self.assertTrue((s2.run_path.parent / "attempt-1" / "run_manifest.json").exists())
        db.close()
        code, _, err = run_cli(["--project", self.proj, "submit", "--mode", "standalone", "--cohort", "C1"])
        self.assertEqual(code, 0, err)
        self.assertEqual(self.db().latest_run("cohort", "C1", "cohort_merge", VERSION).status, "success")
        code, out, _ = run_cli(["--project", self.proj, "status", "--json"])
        self.assertEqual([r for r in json.loads(out) if r["subject"] == "S2"][0]["singleton"], "success(2)")

    def test_assembly_end_to_end(self):
        code, _, err = run_cli(["--project", self.proj, "samples", "add", str(write_tsv(self.tmp, [
            {"sample_id": "K1", "hifi_reads": "k1.bam", "father_id": "S1", "mother_id": "S2"}], "kids.tsv"))])
        self.assertEqual(code, 0, err)
        code, out, err = run_cli(["--project", self.proj, "submit", "--mode", "assembly", "--samples", "S1", "S2", "K1"])
        self.assertEqual(code, 0, err)
        db = self.db()
        for sid in ("S1", "S2", "K1"):
            self.assertEqual(db.latest_run("sample", sid, "assembly", VERSION).status, "success")
        k1 = db.latest_run("sample", "K1", "assembly", VERSION)
        s1 = db.latest_run("sample", "S1", "assembly", VERSION)
        k1_in = json.loads((k1.run_path / "inputs.json").read_text())
        self.assertIn("ugc_wgw_assembly.father_hifi_reads", k1_in)
        self.assertNotIn("ugc_wgw_assembly.father_hifi_reads", json.loads((s1.run_path / "inputs.json").read_text()))
        k1_out = json.loads((k1.run_path / "outputs.json").read_text())
        self.assertEqual(k1_out["ugc_wgw_assembly.hap1_parent"], "S1")
        self.assertTrue(k1_out["ugc_wgw_assembly.trio"])
        self.assertEqual(len(k1_out["ugc_wgw_assembly.zipped_assembly_fastas"]), 2)
        m = json.loads((k1.run_path / "run_manifest.json").read_text())
        self.assertEqual(m["mode"], "assembly")
        self.assertEqual(m["ugc_wgw_manifest"]["cohort_members"], ["S1", "S2"])
        self.assertFalse(any(e["event"] == "run.blocked" for e in self.events()))
        db.close()
        code, out, _ = run_cli(["--project", self.proj, "status", "--mode", "assembly", "--json"])
        rows = json.loads(out)
        self.assertTrue(all(r["type"] == "sample" for r in rows))
        self.assertEqual({r["subject"]: r["assembly"] for r in rows}, {"S1": "success", "S2": "success", "S3": "-", "K1": "success"})
        code, _, err = run_cli(["--project", self.proj, "retry", "--mode", "assembly", "--stage", "assembly", "--samples", "K1"])
        self.assertEqual(code, 0, err)
        self.assertIn("nothing runnable", err)

    def test_version_mismatch(self):
        code, _, err = run_cli(["--project", self.proj, "submit", "--mode", "standalone", "--samples", "S1"],
                               env={"UGC_WGW_FAKE_VERSION": "bogus"})
        self.assertEqual(code, 1)
        run = self.db().latest_run("sample", "S1", "singleton", VERSION)
        self.assertEqual((run.status, run.error_class), ("failed", "version_mismatch"))

    def test_lock_and_logs(self):
        lock = self.cfg.lock_path
        from ugc_wgw.submit import acquire_lock
        fh = acquire_lock(self.cfg)
        try:
            code, _, err = run_cli(["--project", self.proj, "submit", "--mode", "standalone", "--samples", "S1"])
            self.assertEqual(code, 1)
            self.assertIn("holds", err)
        finally:
            fh.close()
        code, _, err = run_cli(["--project", self.proj, "submit", "--mode", "standalone", "--samples", "S1"])
        self.assertEqual(code, 0, err)
        code, out, _ = run_cli(["--project", self.proj, "logs", "S1", "--stage", "singleton", "--tail", "5"])
        self.assertEqual(code, 0)
        self.assertIn("fake miniwdl run singleton S1", out)
        self.assertTrue(lock.exists())


if __name__ == "__main__":
    unittest.main()
