import tempfile
import unittest
from pathlib import Path

from ugc_wgw import cohorts, config, inputs, samples
from ugc_wgw.db import DB
from ugc_wgw.log import Events
from ugc_wgw.outputs import MissingOutputError
from ugc_wgw.stages import declared_inputs, wdl_path

from .helpers import FAKE, REPO, VERSION, make_project, run_cli, seed_success, write_ids, write_tsv


class InputsTest(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.cfg = make_project(self.tmp)
        self.db = DB(self.cfg.db_path)
        self.events = Events(self.db, self.cfg.logs_dir / "events.jsonl")
        samples.add_samples(self.db, self.events, write_tsv(self.tmp, [
            {"sample_id": "S1", "sex": "MALE", "hifi_reads": "s1.bam", "fail_reads": "s1f.bam"},
            {"sample_id": "S2", "hifi_reads": "s2a.bam,s2b.bam"},
            {"sample_id": "S3", "sex": "FEMALE", "hifi_reads": "s3.bam"},
            {"sample_id": "K1", "hifi_reads": "k1.bam", "father_id": "S1", "mother_id": "S3"},
            {"sample_id": "K2", "hifi_reads": "k2.bam", "father_id": "S1"},
            {"sample_id": "K3", "hifi_reads": "k3.bam", "father_id": "NOPE", "mother_id": "S3"},
        ]))
        self.cohort = cohorts.freeze(self.db, self.events, "C1", write_ids(self.tmp, ["S2", "S1", "S3"]))

    def tearDown(self):
        self.db.close()

    def ctx(self, stage_subject: str, subject_id: str, mode: str = "standalone", **kw) -> inputs.BuildContext:
        return inputs.BuildContext(self.cfg, self.db, mode, VERSION, stage_subject, subject_id, cohort=self.cohort, **kw)

    def check_doc(self, stage: str, doc: dict) -> dict:
        declared = declared_inputs(wdl_path(REPO, stage))
        ns = f"ugc_wgw_{stage}."
        names = {k[len(ns):] for k in doc}
        plain = {n for n in names if "." not in n}   # dotted names are nested call inputs, checked by miniwdl
        self.assertTrue(plain <= set(declared), plain - set(declared))
        self.assertTrue({n for n, r in declared.items() if r} <= names)
        self.assertFalse(any(v is None for v in doc.values()))
        return {k[len(ns):]: v for k, v in doc.items()}

    def test_singleton_and_upstream(self):
        doc, members = inputs.generate(self.ctx("sample", "S1"), "singleton")
        d = self.check_doc("singleton", doc)
        self.assertNotIn("sex", d)  # v4 infers sex from coverage; the sheet's sex only strata cohort_freq
        self.assertEqual(len(d["fail_reads"]), 1)
        self.assertEqual(d["ugc_wgw_version"], VERSION)
        self.assertEqual(d["backend"], "HPC")
        self.assertTrue(d["ref_map_file"].endswith("ugc_wgw_ref_map.GRCh38_GIABv3.tsv"))
        self.assertEqual(members, [])
        d2 = self.check_doc("singleton", inputs.generate(self.ctx("sample", "S2"), "singleton")[0])
        self.assertNotIn("fail_reads", d2)
        self.assertEqual(len(d2["hifi_reads"]), 2)
        d3 = self.check_doc("upstream", inputs.generate(self.ctx("sample", "S1", "joint"), "upstream")[0])
        self.assertNotIn("run_sawfish_call", d3)

    def test_deepvariant_switch(self):
        d = self.check_doc("singleton", inputs.generate(self.ctx("sample", "S1"), "singleton")[0])
        self.assertEqual((d["use_gpu"], d["use_parabricks_deepvariant"]), (False, False))
        self.assertNotIn("gpuType", d)
        self.cfg.deepvariant = "gpu"
        d = self.check_doc("upstream", inputs.generate(self.ctx("sample", "S1", "joint"), "upstream")[0])
        self.assertEqual((d["use_gpu"], d["use_parabricks_deepvariant"]), (True, False))
        self.assertNotIn("gpuType", d)                      # empty gpu_type: the plugin turns "" into gpu:N
        self.assertNotIn(inputs.PARABRICKS_GPUS_INPUT, d)
        self.cfg.deepvariant, self.cfg.gpu_type, self.cfg.parabricks_gpus = "parabricks", "a100", 2
        doc = inputs.generate(self.ctx("sample", "S1"), "singleton")[0]
        d = self.check_doc("singleton", doc)
        self.assertEqual((d["use_gpu"], d["use_parabricks_deepvariant"], d["gpuType"]), (True, True, "a100"))
        self.assertEqual(d[inputs.PARABRICKS_GPUS_INPUT], 2)
        self.assertIn("ugc_wgw_singleton." + inputs.PARABRICKS_GPUS_INPUT, doc)
        self.cfg.stage_inputs = {"ugc_wgw_singleton": {"use_gpu": False}}   # overrides still win
        d = self.check_doc("singleton", inputs.generate(self.ctx("sample", "S1"), "singleton")[0])
        self.assertFalse(d["use_gpu"])
        # a nested key must start with a call of the entrypoint
        self.cfg.stage_inputs = {"ugc_wgw_singleton": {"nope.task.input": 1}}
        with self.assertRaises(Exception):
            inputs.generate(self.ctx("sample", "S1"), "singleton")
        self.cfg.stage_inputs = {"ugc_wgw_singleton": {"downstream.some_task.threads": 8}}
        d = self.check_doc("singleton", inputs.generate(self.ctx("sample", "S1"), "singleton")[0])
        self.assertEqual(d["downstream.some_task.threads"], 8)

    def test_config_keys_round_trip_and_validation(self):
        self.cfg.deepvariant, self.cfg.gpu_type, self.cfg.parabricks_gpus = "gpu", "l40s", 6
        config.save(self.cfg)
        c = config.load(self.cfg.project_dir)
        self.assertEqual((c.deepvariant, c.gpu_type, c.parabricks_gpus), ("gpu", "l40s", 6))
        doc = c.to_json()
        doc["deepvariant"] = "tpu"
        with self.assertRaises(Exception):
            config.Config.from_json(self.cfg.project_dir, doc)
        doc["deepvariant"], doc["parabricks_gpus"] = "cpu", 0
        with self.assertRaises(Exception):
            config.Config.from_json(self.cfg.project_dir, doc)
        del doc["deepvariant"], doc["gpu_type"], doc["parabricks_gpus"]   # a 0.2.0 config.json
        c = config.Config.from_json(self.cfg.project_dir, doc)
        self.assertEqual((c.deepvariant, c.gpu_type, c.parabricks_gpus), ("cpu", "", 4))
        proj = self.tmp / "p2"
        code, out, err = run_cli(["init", str(proj), "--code", str(REPO), "--miniwdl", str(FAKE), "--cfg", str(self.cfg.miniwdl_cfg),
                                  "--ref-map", str(self.cfg.ref_map_file), "--deepvariant", "parabricks", "--gpu-type", "a100",
                                  "--parabricks-gpus", "2"])
        self.assertEqual(code, 0, err)
        c = config.load(proj)
        self.assertEqual((c.deepvariant, c.gpu_type, c.parabricks_gpus), ("parabricks", "a100", 2))
        code, out, err = run_cli(["init", str(self.tmp / "p3"), "--code", str(REPO), "--miniwdl", str(FAKE), "--cfg", str(self.cfg.miniwdl_cfg),
                                  "--ref-map", str(self.cfg.ref_map_file), "--parabricks-gpus", "0"])
        self.assertEqual(code, 1)
        self.assertIn("parabricks_gpus", err)

    def test_overrides_win_and_are_checked(self):
        self.cfg.stage_inputs = {"ugc_wgw_singleton": {"use_gpu": True, "use_alignment_chunking": False}}
        d = self.check_doc("singleton", inputs.generate(self.ctx("sample", "S1"), "singleton")[0])
        self.assertTrue(d["use_gpu"])
        self.assertFalse(d["use_alignment_chunking"])
        self.cfg.stage_inputs = {"ugc_wgw_singleton": {"not_an_input": 1}}
        with self.assertRaises(Exception):
            inputs.generate(self.ctx("sample", "S1"), "singleton")

    def test_cohort_merge_standalone(self):
        for sid in ("S1", "S2", "S3"):
            seed_success(self.cfg, self.db, "singleton", "sample", sid, {"sample_id": sid, "inferred_sex": "MALE"})
        doc, members = inputs.generate(self.ctx("cohort", "C1"), "cohort_merge")
        d = self.check_doc("cohort_merge", doc)
        self.assertEqual(d["sample_ids"], ["S2", "S1", "S3"])
        self.assertEqual([Path(p).name.split(".")[0] for p in d["sv_vcfs"]], ["S2", "S1", "S3"])
        self.assertTrue(all(p.endswith("structural_variants.phased.vcf.gz") for p in d["sv_vcfs"]))
        self.assertTrue(all(p.endswith("trgt.sorted.vcf.gz") for p in d["trgt_vcfs"]))
        self.assertTrue(d["run_glnexus"])
        self.assertEqual(len(d["gvcfs"]), 3)
        self.assertNotIn("phased_small_variant_vcfs", d)
        self.assertEqual(d["ugc_wgw_container_registry"], config.DEFAULT_REGISTRY)
        self.assertEqual([m["sample_id"] for m in members], ["S2", "S1", "S3"])
        self.assertTrue(all(m["ugc_wgw_version"] == VERSION and m["stage"] == "singleton" for m in members))
        # toggles change the key set
        self.cfg.stage_inputs = {"ugc_wgw_cohort_merge": {"run_glnexus": False, "merge_phased_small_variants": True}}
        d = self.check_doc("cohort_merge", inputs.generate(self.ctx("cohort", "C1"), "cohort_merge")[0])
        self.assertNotIn("gvcfs", d)
        self.assertEqual(len(d["phased_small_variant_vcfs"]), 3)

    def test_cohort_freq_standalone_and_joint(self):
        # sheet: S1 MALE, S2 (no sex), S3 FEMALE (see setUp); S2's sex comes from the singleton's inferred_sex
        for sid, sex in (("S1", "FEMALE"), ("S2", "FEMALE"), ("S3", "")):
            seed_success(self.cfg, self.db, "singleton", "sample", sid, {"sample_id": sid, "inferred_sex": sex})
        seed_success(self.cfg, self.db, "cohort_merge", "cohort", "C1",
                     {"cohort_id": "C1", "sample_ids": ["S2", "S1", "S3"], "run_glnexus": True})
        doc, members = inputs.generate(self.ctx("cohort", "C1"), "cohort_freq")
        d = self.check_doc("cohort_freq", doc)
        self.assertEqual(d["sample_ids"], ["S2", "S1", "S3"])
        self.assertEqual(d["sample_sexes"], ["XX", "XY", "XX"])          # S1: sheet wins over the inferred FEMALE
        self.assertEqual(d["sample_sex_sources"], ["inferred", "sheet", "sheet"])
        self.assertTrue(d["sv_vcf"].endswith("C1.merged.GRCh38_GIABv3.structural_variants.vcf.gz"))
        self.assertTrue(d["small_variant_vcf"].endswith("C1.merged.GRCh38_GIABv3.small_variants.vcf.gz"))
        self.assertEqual({(m.get("sample_id") or m.get("cohort_id"), m["stage"]) for m in members},
                         {("C1", "cohort_merge"), ("S1", "singleton"), ("S2", "singleton"), ("S3", "singleton")})
        # no joint VCF (run_glnexus off): SV frequencies only, the optional inputs are simply absent
        seed_success(self.cfg, self.db, "cohort_merge", "cohort", "C1",
                     {"cohort_id": "C1", "sample_ids": ["S2", "S1", "S3"], "run_glnexus": False}, attempt=2)
        d = self.check_doc("cohort_freq", inputs.generate(self.ctx("cohort", "C1"), "cohort_freq")[0])
        self.assertNotIn("small_variant_vcf", d)
        self.assertNotIn("small_variant_vcf_index", d)
        # joint mode: the joint VCF comes from cohort_call, the inferred sex from upstream
        for sid in ("S1", "S2", "S3"):
            seed_success(self.cfg, self.db, "upstream", "sample", sid, {"sample_id": sid, "inferred_sex": "MALE"}, mode="joint")
        seed_success(self.cfg, self.db, "cohort_call", "cohort", "C1",
                     {"cohort_id": "C1", "sample_ids": ["S2", "S1", "S3"]}, mode="joint")
        seed_success(self.cfg, self.db, "cohort_merge", "cohort", "C1",
                     {"cohort_id": "C1", "sample_ids": ["S2", "S1", "S3"], "run_glnexus": False}, mode="joint", attempt=3)
        d = self.check_doc("cohort_freq", inputs.generate(self.ctx("cohort", "C1", "joint"), "cohort_freq")[0])
        self.assertTrue(d["small_variant_vcf"].endswith("C1.joint.GRCh38_GIABv3.small_variants.vcf.gz"))
        self.assertEqual(d["sample_sexes"], ["XY", "XY", "XX"])
        self.assertEqual(d["sample_sex_sources"], ["inferred", "sheet", "sheet"])

    def test_missing_prior_output_names_the_member(self):
        seed_success(self.cfg, self.db, "singleton", "sample", "S1", {"sample_id": "S1"})
        with self.assertRaises(MissingOutputError) as cm:
            inputs.generate(self.ctx("cohort", "C1"), "cohort_merge")
        self.assertEqual(cm.exception.subject_id, "S2")
        self.assertEqual(cm.exception.stage, "singleton")

    def test_joint_chain(self):
        for sid in ("S1", "S2", "S3"):
            seed_success(self.cfg, self.db, "upstream", "sample", sid,
                         {"sample_id": sid, "inferred_sex": "" if sid == "S2" else "MALE", "fail_reads": ["/x"] if sid == "S1" else None},
                         mode="joint")
        doc, _ = inputs.generate(self.ctx("cohort", "C1", "joint"), "cohort_call")
        d = self.check_doc("cohort_call", doc)
        self.assertEqual(len(d["gvcfs"]), 3)
        self.assertNotIn("discover_tars", d)
        self.cfg.stage_inputs = {"ugc_wgw_cohort_call": {"run_sawfish_joint_call": True}}
        d = self.check_doc("cohort_call", inputs.generate(self.ctx("cohort", "C1", "joint"), "cohort_call")[0])
        self.assertEqual(len(d["discover_tars"]), 3)
        self.assertTrue(d["run_sawfish_joint_call"])
        self.cfg.stage_inputs = {}
        seed_success(self.cfg, self.db, "cohort_call", "cohort", "C1",
                     {"cohort_id": "C1", "sample_ids": ["S2", "S1", "S3"]}, mode="joint")
        doc, _ = inputs.generate(self.ctx("sample", "S2", "joint"), "downstream")
        d = self.check_doc("downstream", doc)
        self.assertEqual(d["sex"], "")  # empty inferred sex is kept
        self.assertTrue(d["small_variant_vcf"].endswith("S2.C1.joint.GRCh38_GIABv3.small_variants.vcf.gz"))
        self.assertTrue(d["sv_vcf"].endswith("S2.GRCh38_GIABv3.structural_variants.vcf.gz"))  # per-sample sawfish VCF
        self.assertNotIn("aligned_fail_reads", d)
        self.assertEqual(d["upstream_msg"], ["S2: fake upstream message"])
        d1 = self.check_doc("downstream", inputs.generate(self.ctx("sample", "S1", "joint"), "downstream")[0])
        self.assertIn("aligned_fail_reads", d1)
        self.assertTrue(d1["small_variant_vcf"].endswith("S1.C1.joint.GRCh38_GIABv3.small_variants.vcf.gz"))
        for sid in ("S1", "S2", "S3"):
            seed_success(self.cfg, self.db, "downstream", "sample", sid, {"sample_id": sid, "inferred_sex": "MALE"}, mode="joint")
        d = self.check_doc("cohort_merge", inputs.generate(self.ctx("cohort", "C1", "joint"), "cohort_merge")[0])
        self.assertTrue(all(p.endswith("trgt.sorted.vcf.gz") for p in d["trgt_vcfs"]))
        self.assertFalse(d["run_glnexus"])  # joint mode: cohort_call already made the GLnexus VCF
        self.assertNotIn("gvcfs", d)
        self.cfg.stage_inputs = {"cohort_merge": {"run_glnexus": True}}
        d = self.check_doc("cohort_merge", inputs.generate(self.ctx("cohort", "C1", "joint"), "cohort_merge")[0])
        self.assertTrue(all("small_variants.g.vcf.gz" in p for p in d["gvcfs"]))  # from the upstream runs
        self.cfg.stage_inputs = {}

    def test_assembly_sample_and_trio(self):
        d = self.check_doc("assembly", inputs.generate(self.ctx("sample", "K1", "assembly"), "assembly")[0])
        self.assertEqual((d["father_id"], d["mother_id"]), ("S1", "S3"))
        self.assertEqual(len(d["father_hifi_reads"]), 1)
        self.assertTrue(d["mother_hifi_reads"][0].endswith("s3.bam"))
        for k in ("sex", "fail_reads"):
            self.assertNotIn(k, d)
        self.assertEqual(d["ugc_wgw_container_registry"], self.cfg.ugc_wgw_container_registry)  # our hifiasm image
        for sid in ("S1", "K2", "K3"):
            d = self.check_doc("assembly", inputs.generate(self.ctx("sample", sid, "assembly"), "assembly")[0])
            self.assertNotIn("father_id", d)
            self.assertNotIn("father_hifi_reads", d)
            self.assertNotIn("mother_hifi_reads", d)
        self.cfg.assembly_use_parents = False
        d = self.check_doc("assembly", inputs.generate(self.ctx("sample", "K1", "assembly"), "assembly")[0])
        self.assertNotIn("father_hifi_reads", d)
        self.cfg.assembly_use_parents = True
        self.cfg.stage_inputs = {"ugc_wgw_assembly": {"hifiasm_extra_params": "--telo-m CCCTAA", "hifiasm_mem_gb": 200}}
        d = self.check_doc("assembly", inputs.generate(self.ctx("sample", "K1", "assembly"), "assembly")[0])
        self.assertEqual(d["hifiasm_extra_params"], "--telo-m CCCTAA")
        self.assertEqual(d["hifiasm_mem_gb"], 200)

    def test_any_version(self):
        seed_success(self.cfg, self.db, "singleton", "sample", "S1", {"sample_id": "S1"}, version="0.0.9")
        ctx = inputs.BuildContext(self.cfg, self.db, "standalone", VERSION, "cohort", "C1", cohort=self.cohort)
        with self.assertRaises(MissingOutputError):
            inputs.generate(ctx, "cohort_merge")
        for sid in ("S2", "S3"):
            seed_success(self.cfg, self.db, "singleton", "sample", sid, {"sample_id": sid})
        ctx = inputs.BuildContext(self.cfg, self.db, "standalone", VERSION, "cohort", "C1", cohort=self.cohort, any_version=True)
        _, members = inputs.generate(ctx, "cohort_merge")
        self.assertEqual({m["sample_id"]: m["ugc_wgw_version"] for m in members}, {"S1": "0.0.9", "S2": VERSION, "S3": VERSION})


if __name__ == "__main__":
    unittest.main()
