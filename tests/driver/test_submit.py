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

    def test_cohort_qc_never_gates(self):
        """cohort_qc fails (or is simply not there): cohort_merge and cohort_freq run anyway; a later `--stage cohort_qc`
        runs it alone, and the cohort's status shows it."""
        code, out, err = run_cli(["--project", self.proj, "submit", "--mode", "standalone", "--cohort", "C1"], env={"UGC_WGW_FAKE_FAIL_STAGES": "cohort_qc"})
        self.assertEqual(code, 1, err)                       # a failure is still reported ...
        self.assertIn("7 done, 1 blocked", err)              # ... but the stages after cohort_qc ran
        db = self.db()
        try:
            self.assertEqual(db.latest_run("cohort", "C1", "cohort_qc", VERSION).status, "failed")
            self.assertEqual(db.latest_run("cohort", "C1", "cohort_merge", VERSION).status, "success")
            self.assertEqual(db.latest_run("cohort", "C1", "cohort_freq", VERSION).status, "success")
        finally:
            db.close()
        code, out, err = run_cli(["--project", self.proj, "retry", "--mode", "standalone", "--cohort", "C1", "--stage", "cohort_qc"])
        self.assertEqual(code, 0, err)
        db = self.db()
        try:
            qc = db.latest_run("cohort", "C1", "cohort_qc", VERSION)
            self.assertEqual((qc.status, qc.attempt), ("success", 2))
            self.assertEqual(len(db.runs_for("cohort", "C1", "cohort_merge")), 1)   # nothing else re-ran
        finally:
            db.close()
        code, out, _ = run_cli(["--project", self.proj, "status", "--json"])
        self.assertEqual([r.get("cohort_qc") for r in json.loads(out) if r["subject"] == "C1"], ["success(2)"])

    def test_rerun(self):
        """`submit --stage S --rerun`: one new attempt for a subject whose latest run succeeded, once per invocation;
        nothing else runs; the later stages that keep the old outputs are named; the refusals."""
        base = ["--project", self.proj, "submit", "--mode", "standalone", "--cohort", "C1"]
        code, _, err = run_cli(base)
        self.assertEqual(code, 0, err)
        # the refusals: no stage; a per-sample stage without named samples; retry has no --rerun
        code, _, err = run_cli(base + ["--rerun"])
        self.assertNotEqual(code, 0)
        self.assertIn("--rerun needs --stage", err)
        code, _, err = run_cli(base + ["--stage", "singleton", "--rerun"])
        self.assertNotEqual(code, 0)
        self.assertIn("needs --samples or --samples-file", err)
        db = self.db()
        try:
            self.assertEqual(len(db.runs_for("sample", "S1", "singleton")), 1)   # a refusal launches nothing
        finally:
            db.close()
        # a plain submit of a finished stage still does nothing
        code, _, err = run_cli(base + ["--stage", "cohort_merge"])
        self.assertEqual(code, 0, err)
        self.assertIn("launched 0", err)
        # the dry run shows the new attempt and the later stages, and changes nothing
        code, out, err = run_cli(base + ["--stage", "cohort_merge", "--rerun", "--dry-run"])
        self.assertEqual(code, 0, err)
        self.assertIn("runnable=1", out)
        self.assertIn("attempt-2  (rerun of the successful attempt 1)", out)
        self.assertIn("cohort_freq (cohort C1), cohort_tr (cohort C1), cohort_methyl (cohort C1)", out)
        # the rerun: attempt 2 of cohort_merge alone, exactly once, `current` moved, the reason recorded
        code, _, err = run_cli(base + ["--stage", "cohort_merge", "--rerun"])
        self.assertEqual(code, 0, err)
        self.assertIn("launched 1", err)
        self.assertIn("keep the outputs of the attempt being replaced", err)
        db = self.db()
        try:
            merge = db.latest_run("cohort", "C1", "cohort_merge", VERSION)
            self.assertEqual((merge.status, merge.attempt, merge.meta.get("rerun_of")), ("success", 2, 1))
            self.assertEqual(len(db.runs_for("cohort", "C1", "cohort_merge")), 2)
            self.assertEqual(os.readlink(merge.run_path.parent / "current"), "attempt-2")
            self.assertTrue((merge.run_path.parent / "attempt-1" / "run_manifest.json").exists())
            self.assertEqual(len(db.runs_for("cohort", "C1", "cohort_freq")), 1)     # later stages are not rerun for the user
            self.assertEqual(len(db.runs_for("sample", "S1", "singleton")), 1)
        finally:
            db.close()
        reruns = [e for e in self.events() if e["event"] == "run.rerun"]
        self.assertEqual([(e["detail"]["stage"], e["detail"]["reason"]) for e in reruns],
                         [("cohort_merge", "rerun of the successful attempt 1")])
        # a later stage rerun picks the new attempt up; an optional stage names no later stage
        code, _, err = run_cli(base + ["--stage", "cohort_tr", "--rerun"])
        self.assertEqual(code, 0, err)
        self.assertNotIn("keep the outputs", err)
        db = self.db()
        try:
            tr = db.latest_run("cohort", "C1", "cohort_tr", VERSION)
            self.assertEqual((tr.status, tr.attempt), ("success", 2))
            used = {m["stage"]: m["run_id"] for m in json.loads((tr.run_path / "run_manifest.json").read_text())["cohort_members"]}
            self.assertEqual(used["cohort_merge"], db.latest_run("cohort", "C1", "cohort_merge", VERSION).run_id)
        finally:
            db.close()
        # per-sample: only the named sample, and the later stages are counted for it
        code, _, err = run_cli(base + ["--stage", "singleton", "--samples", "S2", "--rerun"])
        self.assertEqual(code, 0, err)
        self.assertIn("cohort_merge (cohort C1)", err)
        db = self.db()
        try:
            self.assertEqual([len(db.runs_for("sample", s, "singleton")) for s in ("S1", "S2", "S3")], [1, 2, 1])
        finally:
            db.close()
        code, out, _ = run_cli(["--project", self.proj, "status", "--json"])
        self.assertEqual([r for r in json.loads(out) if r["subject"] == "S2"][0]["singleton"], "success(2)")

    def test_rerun_leaves_failed_and_missing_alone(self):
        """--rerun changes nothing for a subject whose latest run failed (retry's business) or that never ran."""
        base = ["--project", self.proj, "submit", "--mode", "standalone", "--cohort", "C1"]
        code, _, err = run_cli(base, env={"UGC_WGW_FAKE_FAIL": "S2"})
        self.assertEqual(code, 1)
        code, _, err = run_cli(base + ["--stage", "singleton", "--samples", "S1", "S2", "--rerun"])
        self.assertEqual(code, 0, err)
        self.assertIn("ugc-wgw retry", err)                  # S2 is still blocked on its failure
        db = self.db()
        try:
            self.assertEqual([(r.status, r.attempt) for r in (db.latest_run("sample", s, "singleton", VERSION) for s in ("S1", "S2"))],
                             [("success", 2), ("failed", 1)])
        finally:
            db.close()
        # a cohort stage that never ran and cannot yet: blocked on its prerequisite as without --rerun
        code, _, err = run_cli(base + ["--stage", "cohort_merge", "--rerun"])
        self.assertIn("waiting for singleton of S2", err)

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
        self.assertEqual(len(k1_out["ugc_wgw_assembly.hap_fastas"]), 2)
        self.assertEqual(k1_out["ugc_wgw_assembly.haplotype_labels"], ["pat", "mat"])
        # the mode continues into assembly_call, fed by the contract; the smoke-sized fake infers karyotype unknown
        for sid in ("S1", "S2", "K1"):
            call = db.latest_run("sample", sid, "assembly_call", VERSION)
            self.assertEqual(call.status, "success", sid)
            call_in = json.loads((call.run_path / "inputs.json").read_text())
            self.assertEqual(call_in["ugc_wgw_assembly_call.hap_fastas"], k1_out["ugc_wgw_assembly.hap_fastas"] if sid == "K1" else call_in["ugc_wgw_assembly_call.hap_fastas"])
            self.assertEqual(call_in["ugc_wgw_assembly_call.trio"], sid == "K1")
            call_out = json.loads((call.run_path / "outputs.json").read_text())
            self.assertTrue(call_out["ugc_wgw_assembly_call.dipcall_vcf"].endswith(f"{sid}.dipcall.GRCh38_GIABv3.vcf.gz"))
            self.assertEqual(call_out["ugc_wgw_assembly_call.inferred_karyotype"], "unknown")
        m = json.loads((k1.run_path / "run_manifest.json").read_text())
        self.assertEqual(m["mode"], "assembly")
        self.assertEqual(m["ugc_wgw_manifest"]["cohort_members"], ["S1", "S2"])
        # only assembly_call waits (for assembly) and cohort_pangenome (for --cohort); the first stage is never blocked
        self.assertTrue(all("assembly_call" in json.dumps(e) or "cohort_pangenome" in json.dumps(e)
                            for e in self.events() if e["event"] == "run.blocked"))
        db.close()
        code, out, _ = run_cli(["--project", self.proj, "status", "--mode", "assembly", "--json"])
        rows = json.loads(out)
        self.assertEqual({r["subject"]: r["assembly"] for r in rows if r["type"] == "sample"}, {"S1": "success", "S2": "success", "S3": "-", "K1": "success"})
        self.assertEqual({r["subject"]: r["cohort_pangenome"] for r in rows if r["type"] == "cohort"}, {"C1": "-"})   # the mode's cohort stage, not run
        code, _, err = run_cli(["--project", self.proj, "retry", "--mode", "assembly", "--stage", "assembly", "--samples", "K1"])
        self.assertEqual(code, 0, err)
        self.assertIn("nothing runnable", err)
        # with --cohort the mode continues into cohort_pangenome over the members' assemblies (S3 is assembled on the way)
        code, _, err = run_cli(["--project", self.proj, "submit", "--mode", "assembly", "--cohort", "C1"])
        self.assertEqual(code, 0, err)
        db = self.db()
        self.assertEqual(db.latest_run("sample", "S3", "assembly_call", VERSION).status, "success")
        pg = db.latest_run("cohort", "C1", "cohort_pangenome", VERSION)
        self.assertEqual(pg.status, "success")
        pg_in = json.loads((pg.run_path / "inputs.json").read_text())
        self.assertEqual(pg_in["ugc_wgw_cohort_pangenome.sample_ids"], ["S1", "S2", "S3"])
        self.assertEqual([p.split("/")[-1] for p in pg_in["ugc_wgw_cohort_pangenome.hap1_fastas"]],
                         [f"{s}.asm.bp.hap1.p_ctg.fasta.gz" for s in ("S1", "S2", "S3")])
        self.assertEqual(pg_in["ugc_wgw_cohort_pangenome.trios"], [False, False, False])
        self.assertEqual(pg_in["ugc_wgw_cohort_pangenome.sample_sexes"], ["", "", ""])
        self.assertNotIn("ugc_wgw_cohort_pangenome.sv_vcf", pg_in)                  # no cohort_merge of C1: truvari skipped
        pg_out = json.loads((pg.run_path / "outputs.json").read_text())
        self.assertEqual(pg_out["ugc_wgw_cohort_pangenome.chromosomes_built"], ["chr1", "chr20", "chrX"])
        self.assertTrue(pg_out["ugc_wgw_cohort_pangenome.pangenome_wave_vcf"].endswith("C1.pangenome.GRCh38_GIABv3.wave.vcf.gz"))
        pm = json.loads((pg.run_path / "run_manifest.json").read_text())
        self.assertEqual(pm["ugc_wgw_manifest"]["cohort_members"], ["S1", "S2", "S3"])
        db.close()
        code, out, _ = run_cli(["--project", self.proj, "status", "--mode", "assembly", "--json"])
        rows = json.loads(out)
        self.assertEqual({r["subject"]: r.get("cohort_pangenome") for r in rows if r["type"] == "cohort"}, {"C1": "success"})

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
