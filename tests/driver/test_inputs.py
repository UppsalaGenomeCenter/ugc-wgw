import tempfile
import unittest
from pathlib import Path

from ugc_wgw import cohorts, config, inputs, samples
from ugc_wgw.db import DB
from ugc_wgw.log import Events
from ugc_wgw.outputs import MissingOutputError
from ugc_wgw.stages import declared_inputs, wdl_path
from ugc_wgw.util import UgcError

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
        self.assertEqual((d["use_gpu"], d["use_parabricks_deepvariant"]), (True, True))
        self.assertNotIn("gpuType", d)   # the type reaches the plugin through the run's environment, not the inputs
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

    def test_pblab_choices(self):
        """aligner and the pbdv flavours: each alone and together, what the run records as its tools, the refusals."""
        tools = lambda stage, doc, kind="reads": inputs.tools_record(stage, doc, kind)
        doc = inputs.generate(self.ctx("sample", "S1"), "singleton")[0]
        d = self.check_doc("singleton", doc)
        self.assertNotIn("aligner", d)                       # a default project keeps the inputs document it always had
        self.assertNotIn("use_pbdv", d)
        self.assertEqual(tools("singleton", doc), {"aligner": "pbmm2", "small_variant_caller": "deepvariant", "small_variant_device": "cpu"})
        self.cfg.aligner = "pbweave"                          # the aligner alone
        doc = inputs.generate(self.ctx("sample", "S1"), "singleton")[0]
        d = self.check_doc("singleton", doc)
        self.assertEqual((d["aligner"], d["use_gpu"], d.get("use_pbdv")), ("pbweave", False, None))
        self.assertEqual(tools("singleton", doc)["aligner"], "pbweave")
        self.cfg.aligner, self.cfg.deepvariant = "pbmm2", "pbdv"   # the caller alone, on a GPU
        doc = inputs.generate(self.ctx("sample", "S1", "joint"), "upstream")[0]
        d = self.check_doc("upstream", doc)
        self.assertEqual((d.get("aligner"), d["use_gpu"], d["use_pbdv"], d["use_parabricks_deepvariant"]), (None, True, True, False))
        self.assertNotIn(inputs.PARABRICKS_GPUS_INPUT, d)
        self.assertEqual(tools("upstream", doc), {"aligner": "pbmm2", "small_variant_caller": "pbdv", "small_variant_device": "gpu"})
        self.assertEqual(tools("upstream", doc, "bam")["aligner"], "external")   # an aligned BAM came with the sample
        self.cfg.aligner, self.cfg.deepvariant = "pbweave", "pbdv-cpu"   # both, without a GPU
        doc = inputs.generate(self.ctx("sample", "S1"), "singleton")[0]
        d = self.check_doc("singleton", doc)
        self.assertEqual((d["aligner"], d["use_gpu"], d["use_pbdv"]), ("pbweave", False, True))
        self.assertEqual(tools("singleton", doc), {"aligner": "pbweave", "small_variant_caller": "pbdv", "small_variant_device": "cpu"})
        self.cfg.stage_inputs = {"ugc_wgw_singleton": {"aligner": "pbmm2"}}   # an override is told as it ran
        doc = inputs.generate(self.ctx("sample", "S1"), "singleton")[0]
        self.assertEqual(tools("singleton", doc)["aligner"], "pbmm2")
        self.cfg.stage_inputs = {}
        self.cfg.deepvariant = "parabricks"
        doc = inputs.generate(self.ctx("sample", "S1"), "singleton")[0]
        self.assertEqual(tools("singleton", doc)["small_variant_caller"], "parabricks_deepvariant")
        self.assertEqual(tools("intake", {}), {"aligner": "external", "small_variant_caller": "external"})
        self.assertIsNone(tools("downstream", {}))
        # the engine environment: a GPU type only for the flavours that ask for a GPU
        from ugc_wgw import engine
        for mode, expect in (("pbdv", True), ("pbdv-cpu", False), ("cpu", False), ("gpu", True)):
            self.cfg.deepvariant = mode
            self.assertEqual(bool(engine.run_env(self.cfg)), expect, mode)
        # validation, in the config and in a profile
        self.assertEqual(config.check_aligner(" PBweave "), "pbweave")
        for bad in ("bwa", "", "minimap2"):
            with self.assertRaises(UgcError):
                config.check_aligner(bad)
        self.assertEqual(config.check_deepvariant("pbdv-cpu"), "pbdv-cpu")
        self.cfg.aligner, self.cfg.deepvariant = "pbweave", "pbdv"
        config.save(self.cfg)
        again = config.load(self.cfg.project_dir)
        self.assertEqual((again.aligner, again.deepvariant), ("pbweave", "pbdv"))

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
        # pbrun's partition mode wants an even GPU count: refused at init, in config.json and in a profile
        code, _, err = run_cli(["init", str(self.tmp / "odd"), "--code", str(REPO), "--miniwdl", str(FAKE), "--cfg", str(self.cfg.miniwdl_cfg),
                                "--ref-map", str(self.cfg.ref_map_file), "--deepvariant", "parabricks", "--parabricks-gpus", "3"])
        self.assertEqual(code, 1)
        self.assertIn("parabricks_gpus must be even", err)
        for odd in (1, 3, 5):
            with self.assertRaisesRegex(UgcError, "must be even"):
                config.check_gpus(odd)
        self.assertEqual([config.check_gpus(n) for n in (2, 4, 6)], [2, 4, 6])
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
        # the mitochondrial sets: every member has one (0.15 runs), so no mitorsaw inputs travel
        self.assertTrue(d["merge_mitochondria"])
        self.assertEqual([Path(p).name for p in d["mt_small_variant_vcfs"]], [f"{s}.GRCh38_GIABv3.mt.small_variants.vcf.gz" for s in ("S2", "S1", "S3")])
        self.assertEqual(len(d["mt_sv_vcf_indices"]), 3)
        self.assertTrue(all(Path(p).name.endswith(".mt.haplogroups.tsv") for p in d["mt_haplogroups"]))
        self.assertNotIn("mitorsaw_vcfs", d)
        # toggles change the key set
        self.cfg.stage_inputs = {"ugc_wgw_cohort_merge": {"run_glnexus": False, "merge_phased_small_variants": True}}
        d = self.check_doc("cohort_merge", inputs.generate(self.ctx("cohort", "C1"), "cohort_merge")[0])
        self.assertNotIn("gvcfs", d)
        self.assertEqual(len(d["phased_small_variant_vcfs"]), 3)
        # a member whose run predates the collapse: null in the mt arrays, its mitorsaw outputs travel, a warning names it
        seed_success(self.cfg, self.db, "singleton", "sample", "S1", {"sample_id": "S1", "inferred_sex": "MALE", "legacy": True}, attempt=2)
        self.cfg.stage_inputs = {}
        with self.assertLogs("ugc-wgw", level="WARNING") as logs:
            d = self.check_doc("cohort_merge", inputs.generate(self.ctx("cohort", "C1"), "cohort_merge")[0])
        self.assertTrue(any("no collapsed mitochondrial set" in m and "S1" in m for m in logs.output), logs.output)
        self.assertEqual([p is None for p in d["mt_small_variant_vcfs"]], [False, True, False])
        self.assertEqual([p is None for p in d["mt_haplogroups"]], [False, True, False])
        self.assertEqual([None if p is None else Path(p).name for p in d["mitorsaw_vcfs"]], [None, "S1.GRCh38_GIABv3.mitorsaw.vcf.gz", None])
        self.assertEqual([p is None for p in d["mitorsaw_hap_stats"]], [True, False, True])
        # the off switch leaves the branch out
        self.cfg.stage_inputs = {"ugc_wgw_cohort_merge": {"merge_mitochondria": False}}
        d = self.check_doc("cohort_merge", inputs.generate(self.ctx("cohort", "C1"), "cohort_merge")[0])
        self.assertFalse(d["merge_mitochondria"])
        self.assertNotIn("mt_small_variant_vcfs", d)

    def test_reference_map_keys(self):
        """A map rendered by an older version lacks the keys the entrypoints read: refused before any job, naming them
        (the 864's cohort_methyl at 0.16.0 failed in the engine with `Map key not found` on a 0.15.0 map)."""
        from ugc_wgw.stages import ref_map_keys, wdl_path
        keys = ref_map_keys(wdl_path(self.cfg.code_dir, "cohort_methyl"))
        self.assertTrue({"methbat_map", "methbat_report_regions_tsv", "methbat_atlas_tsv", "name"} <= keys, keys)
        self.assertIn("fasta", ref_map_keys(wdl_path(self.cfg.code_dir, "singleton")))   # through the override's imports
        for sid in ("S1", "S2", "S3"):
            seed_success(self.cfg, self.db, "singleton", "sample", sid, {"sample_id": sid, "inferred_sex": "MALE"})
        ref = Path(self.cfg.ref_map_file)
        full = ref.read_text()
        ref.write_text("".join(line for line in full.splitlines(True) if not line.startswith(("methbat_map\t", "methbat_atlas_tsv\t"))))
        with self.assertRaises(UgcError) as cm:
            inputs.generate(self.ctx("cohort", "C1"), "cohort_methyl")
        self.assertIn("lacks methbat_atlas_tsv, methbat_map", str(cm.exception))
        self.assertIn("upgrade --install", str(cm.exception))
        ref.write_text(full)
        inputs.generate(self.ctx("cohort", "C1"), "cohort_methyl")

    def test_cohort_methyl(self):
        """The methylation layer: pileup, profile and report per member (null report + a warning for a run from before the
        task), the labels (sex group, then the project's methyl_label_columns as column=value), a member without a pileup
        left out, an unknown column refused; joint mode reads downstream."""
        for sid in ("S1", "S2", "S3"):
            seed_success(self.cfg, self.db, "singleton", "sample", sid, {"sample_id": sid, "inferred_sex": "MALE" if sid != "S3" else "FEMALE"})
        d = self.check_doc("cohort_methyl", inputs.generate(self.ctx("cohort", "C1"), "cohort_methyl")[0])
        self.assertEqual(d["sample_ids"], ["S2", "S1", "S3"])
        self.assertEqual(d["sample_labels"], ["XY", "XY", "XX"])   # S2 unsexed in the sheet: coverage says MALE
        self.assertEqual(len(d["cpg_pileup_beds"]), 3)
        self.assertTrue(all(p.endswith(".cpg.bed.gz") for p in d["cpg_pileup_beds"]))
        self.assertTrue(all(p.endswith(".methbat.profile.tsv") for p in d["methbat_profiles"]))
        self.assertTrue(all(p and p.endswith(".methbat.report.tsv") for p in d["methbat_reports"]))
        self.assertNotIn("compare_pairs", d)   # the WDL default (XX:XY) unless the project sets it
        # a legacy member: null report, the warning names it
        seed_success(self.cfg, self.db, "singleton", "sample", "S1", {"sample_id": "S1", "inferred_sex": "MALE", "legacy": True}, attempt=2)
        with self.assertLogs("ugc-wgw", level="WARNING") as logs:
            d = self.check_doc("cohort_methyl", inputs.generate(self.ctx("cohort", "C1"), "cohort_methyl")[0])
        self.assertEqual(d["methbat_reports"][1], None)
        self.assertTrue(any("1 member(s) have no methbat_report" in m and "S1" in m for m in logs.output), logs.output)
        # sheet columns as labels: the sheet's extra columns live in meta
        self.db.conn.execute("UPDATE samples SET meta_json = ? WHERE sample_id = ?", ('{"project": "plate 3/A", "phenotype": "0"}', "S1"))
        self.db.conn.execute("UPDATE samples SET meta_json = ? WHERE sample_id = ?", ('{"project": "plate 3/A"}', "S2"))
        self.db.conn.commit()
        self.cfg.methyl_label_columns = ["project"]
        d = self.check_doc("cohort_methyl", inputs.generate(self.ctx("cohort", "C1"), "cohort_methyl")[0])
        self.assertEqual(d["sample_labels"], ["XY;project=plate_3_A", "XY;project=plate_3_A", "XX"])
        self.cfg.methyl_label_columns = ["project", "batch"]
        with self.assertRaises(UgcError) as cm:
            inputs.generate(self.ctx("cohort", "C1"), "cohort_methyl")
        self.assertIn("'batch'", str(cm.exception))
        self.cfg.methyl_label_columns = []
        # joint mode: downstream's pileups, upstream's sex
        for sid in ("S1", "S2", "S3"):
            seed_success(self.cfg, self.db, "upstream", "sample", sid, {"sample_id": sid, "inferred_sex": "MALE"}, mode="joint")
            seed_success(self.cfg, self.db, "downstream", "sample", sid, {"sample_id": sid}, mode="joint")
        d = self.check_doc("cohort_methyl", inputs.generate(self.ctx("cohort", "C1", "joint"), "cohort_methyl")[0])
        self.assertEqual(d["sample_labels"], ["XY", "XY", "XX"])   # S3's sheet sex wins over the coverage call
        self.assertTrue(all("/downstream/" in p for p in d["cpg_pileup_beds"]))

    def test_cohort_mqtl(self):
        """The association stage reads cohort_merge's SV and TRGT VCFs, cohort_methyl's annotated profiles and cell-type
        table, the members' sex groups, and cohort_qc's ancestry table only when a run has one. It needs the first two."""
        for sid, sex in (("S1", "FEMALE"), ("S2", "FEMALE"), ("S3", "")):
            seed_success(self.cfg, self.db, "singleton", "sample", sid, {"sample_id": sid, "inferred_sex": sex})
        seed_success(self.cfg, self.db, "cohort_merge", "cohort", "C1", {"cohort_id": "C1", "sample_ids": ["S2", "S1", "S3"], "run_glnexus": True})
        with self.assertRaises((UgcError, MissingOutputError)):
            inputs.generate(self.ctx("cohort", "C1"), "cohort_mqtl")              # no cohort_methyl run yet
        seed_success(self.cfg, self.db, "cohort_methyl", "cohort", "C1", {"cohort_id": "C1", "sample_ids": ["S2", "S1"]})   # S3 has no methylation
        doc, members = inputs.generate(self.ctx("cohort", "C1"), "cohort_mqtl")
        d = self.check_doc("cohort_mqtl", doc)
        self.assertEqual((d["sample_ids"], d["sample_sexes"]), (["S2", "S1", "S3"], ["XX", "XY", "XX"]))
        self.assertEqual([p.rsplit("/", 1)[-1] for p in d["scored_profiles"]],
                         ["S2.GRCh38_GIABv3.methbat.scored.annotated.tsv", "S1.GRCh38_GIABv3.methbat.scored.annotated.tsv"])
        self.assertTrue(d["methyl_deconvolution"].endswith("C1.methyl.GRCh38_GIABv3.deconvolution.tsv"))
        self.assertTrue(d["sv_vcf"].endswith("C1.merged.GRCh38_GIABv3.structural_variants.vcf.gz"))
        self.assertTrue(d["sv_vcf_index"].endswith(".structural_variants.vcf.gz.tbi"))
        self.assertTrue(d["trgt_vcf"].endswith("C1.merged.GRCh38_GIABv3.trgt.vcf.gz"))
        self.assertNotIn("ancestry", d)                                           # no cohort_qc run: the model goes without the components
        self.assertNotIn("window_bp", d)
        stages_used = {(m.get("sample_id") or m.get("cohort_id"), m["stage"]) for m in members}
        self.assertTrue({("C1", "cohort_merge"), ("C1", "cohort_methyl")} <= stages_used)
        seed_success(self.cfg, self.db, "cohort_qc", "cohort", "C1", {"cohort_id": "C1", "sample_ids": ["S2", "S1", "S3"], "run_ancestry": True})
        self.cfg.stage_inputs = {"cohort_mqtl": {"min_members": 4, "exclude_ids": ["S3"]}}
        d = self.check_doc("cohort_mqtl", inputs.generate(self.ctx("cohort", "C1"), "cohort_mqtl")[0])
        self.assertTrue(d["ancestry"].endswith("C1.qc.GRCh38_GIABv3.somalier-ancestry.tsv"))
        self.assertEqual((d["min_members"], d["exclude_ids"]), (4, ["S3"]))
        self.cfg.stage_inputs = {}

    def test_cohort_tr(self):
        """The tandem-repeat stage reads cohort_merge's TRGT VCF, its index and the LPS table, and the members' sex groups;
        no per-sample file. Without a cohort_merge run it has nothing to read."""
        for sid, sex in (("S1", "FEMALE"), ("S2", "FEMALE"), ("S3", "")):
            seed_success(self.cfg, self.db, "singleton", "sample", sid, {"sample_id": sid, "inferred_sex": sex})
        with self.assertRaises((UgcError, MissingOutputError)):
            inputs.generate(self.ctx("cohort", "C1"), "cohort_tr")
        seed_success(self.cfg, self.db, "cohort_merge", "cohort", "C1", {"cohort_id": "C1", "sample_ids": ["S2", "S1", "S3"], "run_glnexus": True})
        doc, members = inputs.generate(self.ctx("cohort", "C1"), "cohort_tr")
        d = self.check_doc("cohort_tr", doc)
        self.assertEqual(d["sample_ids"], ["S2", "S1", "S3"])
        self.assertEqual(d["sample_sexes"], ["XX", "XY", "XX"])          # S1: the sheet's MALE wins over the inferred FEMALE
        self.assertTrue(d["trgt_vcf"].endswith("C1.merged.GRCh38_GIABv3.trgt.vcf.gz"))
        self.assertTrue(d["trgt_vcf_index"].endswith("C1.merged.GRCh38_GIABv3.trgt.vcf.gz.tbi"))
        self.assertTrue(d["trgt_lps"].endswith(".trgt.lps.tsv"))
        self.assertNotIn("outlier_min_called", d)                        # the WDL defaults unless the project's stage_inputs set them
        self.assertIn(("C1", "cohort_merge"), {(m.get("sample_id") or m.get("cohort_id"), m["stage"]) for m in members})
        self.cfg.stage_inputs = {"cohort_tr": {"outlier_min_called": 4}}
        d = self.check_doc("cohort_tr", inputs.generate(self.ctx("cohort", "C1"), "cohort_tr")[0])
        self.assertEqual(d["outlier_min_called"], 4)
        self.cfg.stage_inputs = {}

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
        self.assertEqual(d["sample_sex_sources"], ["coverage", "declared", "declared"])
        self.assertTrue(d["sv_vcf"].endswith("C1.merged.GRCh38_GIABv3.structural_variants.vcf.gz"))
        self.assertTrue(d["small_variant_vcf"].endswith("C1.merged.GRCh38_GIABv3.small_variants.vcf.gz"))
        self.assertNotIn("mt_small_variant_vcf", d)   # that cohort_merge seed merged no mitochondrial set (0.14 shape)
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
        self.assertEqual(d["sample_sex_sources"], ["coverage", "declared", "declared"])
        # a cohort_merge that merged the mitochondrial sets (0.15): the mt VCFs travel to the mt branch
        seed_success(self.cfg, self.db, "cohort_merge", "cohort", "C1", {"cohort_id": "C1", "sample_ids": ["S2", "S1", "S3"], "run_glnexus": True,
                                                                         "merge_mitochondria": True}, attempt=7)
        d = self.check_doc("cohort_freq", inputs.generate(self.ctx("cohort", "C1"), "cohort_freq")[0])
        self.assertTrue(d["mt_small_variant_vcf"].endswith("C1.merged.GRCh38_GIABv3.mt_small_variants.vcf.gz"))
        self.assertTrue(d["mt_sv_vcf_index"].endswith("C1.merged.GRCh38_GIABv3.mt_svs.vcf.gz.tbi"))

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

    def test_assembly_call(self):
        seed_success(self.cfg, self.db, "assembly", "sample", "K1", {"sample_id": "K1", "hifi_reads": ["/k"], "father_id": "S1", "mother_id": "S3",
                                                                      "father_hifi_reads": ["/f"], "mother_hifi_reads": ["/m"]}, mode="assembly")
        seed_success(self.cfg, self.db, "assembly", "sample", "S1", {"sample_id": "S1", "hifi_reads": ["/s"]}, mode="assembly")
        # trio-binned, no read-based run: the contract travels, nothing optional
        ctx = self.ctx("sample", "K1", "assembly")
        doc, members = inputs.generate(ctx, "assembly_call")
        d = self.check_doc("assembly_call", doc)
        self.assertTrue(d["trio"])
        self.assertEqual(d["haplotype_labels"], ["pat", "mat"])
        self.assertEqual(len(d["hap_fastas"]), 2)
        self.assertTrue(d["hap_fastas"][0].endswith("K1.asm.dip.hap1.p_ctg.fasta.gz"))
        self.assertNotIn("read_small_variant_vcf", d)
        self.assertTrue(d["flagger_bed"].endswith("K1.asm.flagger.bed"))
        self.assertEqual([m["stage"] for m in members], ["assembly"])
        # not trio-binned, with a singleton of an older version: its phased VCF is passed and recorded, never required
        seed_success(self.cfg, self.db, "singleton", "sample", "S1", {"sample_id": "S1", "hifi_reads": ["/s"]}, version="0.0.9")
        doc, members = inputs.generate(self.ctx("sample", "S1", "assembly"), "assembly_call")
        d = self.check_doc("assembly_call", doc)
        self.assertFalse(d["trio"])
        self.assertEqual(d["haplotype_labels"], ["hap1", "hap2"])
        self.assertTrue(d["read_small_variant_vcf"].endswith("phased.vcf.gz") or "small_variants" in d["read_small_variant_vcf"], d["read_small_variant_vcf"])
        self.assertTrue(d["read_small_variant_vcf_index"].endswith(".tbi"))
        self.assertEqual([(m["stage"], m["ugc_wgw_version"]) for m in members], [("assembly", VERSION), ("singleton", "0.0.9")])
        # downstream (joint) serves when there is no singleton
        seed_success(self.cfg, self.db, "assembly", "sample", "S2", {"sample_id": "S2", "hifi_reads": ["/s2"]}, mode="assembly")
        seed_success(self.cfg, self.db, "downstream", "sample", "S2", {"sample_id": "S2", "sex": "MALE"}, mode="joint")
        d = self.check_doc("assembly_call", inputs.generate(self.ctx("sample", "S2", "assembly"), "assembly_call")[0])
        self.assertIn("read_small_variant_vcf", d)
        # no assembly yet: a prerequisite, not optional
        with self.assertRaises(MissingOutputError):
            inputs.generate(self.ctx("sample", "S3", "assembly"), "assembly_call")

    def test_cohort_pangenome(self):
        for sid in ("S2", "S1", "S3"):
            seed_success(self.cfg, self.db, "assembly", "sample", sid, {"sample_id": sid, "hifi_reads": ["/s"]}, mode="assembly")
        # assembly_call of every member is a prerequisite (the karyotype strata), like the assembly itself
        with self.assertRaises(MissingOutputError):
            inputs.generate(self.ctx("cohort", "C1", "assembly"), "cohort_pangenome")
        for sid in ("S2", "S1", "S3"):
            seed_success(self.cfg, self.db, "assembly_call", "sample", sid, {"sample_id": sid, "hap_fastas": ["/a", "/b"], "haplotype_labels": ["hap1", "hap2"], "trio": False}, mode="assembly")
        doc, members = inputs.generate(self.ctx("cohort", "C1", "assembly"), "cohort_pangenome")
        d = self.check_doc("cohort_pangenome", doc)
        self.assertEqual(d["sample_ids"], ["S2", "S1", "S3"])                       # the frozen order
        self.assertEqual([p.split("/")[-1] for p in d["hap1_fastas"]], [f"{s}.asm.bp.hap1.p_ctg.fasta.gz" for s in ("S2", "S1", "S3")])
        self.assertEqual([p.split("/")[-1] for p in d["hap2_fastas"]], [f"{s}.asm.bp.hap2.p_ctg.fasta.gz" for s in ("S2", "S1", "S3")])
        self.assertEqual(d["trios"], [False, False, False])
        self.assertEqual(d["sample_sexes"], ["", "", ""])                           # the fake infers karyotype unknown
        self.assertEqual(d["sample_sex_sources"], ["", "", ""])
        self.assertNotIn("sv_vcf", d)                                               # no cohort_merge: truvari is skipped, not required
        self.assertNotIn("chromosomes", d)                                          # the map's list applies
        self.assertEqual([(m["sample_id"], m["stage"]) for m in members],
                         [(s, st) for s in ("S2", "S1", "S3") for st in ("assembly", "assembly_call")])
        # the cohort's own SV merge (any mode, any version) is compared with; pangenome_sv_cohort points elsewhere
        seed_success(self.cfg, self.db, "cohort_merge", "cohort", "C1", {"cohort_id": "C1", "sample_ids": ["S2", "S1", "S3"]}, version="0.0.9")
        doc, members = inputs.generate(self.ctx("cohort", "C1", "assembly"), "cohort_pangenome")
        d = self.check_doc("cohort_pangenome", doc)
        self.assertTrue(d["sv_vcf"].endswith("C1.merged.GRCh38_GIABv3.structural_variants.vcf.gz"))
        self.assertTrue(d["sv_vcf_index"].endswith(".tbi"))
        self.assertEqual(members[-1], {"cohort_id": "C1", "stage": "cohort_merge", "ugc_wgw_version": "0.0.9", "run_id": "C1-cohort_merge-a1-seed"})
        self.cfg.pangenome_sv_cohort = "ELSEWHERE"
        d = self.check_doc("cohort_pangenome", inputs.generate(self.ctx("cohort", "C1", "assembly"), "cohort_pangenome")[0])
        self.assertNotIn("sv_vcf", d)
        self.cfg.pangenome_sv_cohort = ""
        # a trio-binned member keeps its phase; the chromosomes override travels
        seed_success(self.cfg, self.db, "assembly", "sample", "S1", {"sample_id": "S1", "hifi_reads": ["/s"], "father_id": "S2", "mother_id": "S3",
                                                                      "father_hifi_reads": ["/f"], "mother_hifi_reads": ["/m"]}, mode="assembly", attempt=2)
        self.cfg.stage_inputs = {"cohort_pangenome": {"chromosomes": ["chr20"]}}
        d = self.check_doc("cohort_pangenome", inputs.generate(self.ctx("cohort", "C1", "assembly"), "cohort_pangenome")[0])
        self.assertEqual(d["trios"], [False, True, False])
        self.assertEqual(d["chromosomes"], ["chr20"])
        self.cfg.stage_inputs = {}

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
