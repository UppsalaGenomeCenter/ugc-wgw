"""`ugc-wgw stage-inputs`: the who-fills-what table against the real builders, the listing, the override checks,
and the nested listing through the engine's Python when one with miniwdl is at hand."""
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from ugc_wgw import cohorts, config, inputs, samples, stage_inputs
from ugc_wgw.db import DB
from ugc_wgw.engine import engine_python
from ugc_wgw.log import Events
from ugc_wgw.stages import STAGES, declared_input_specs, input_descriptions, wdl_path

from .helpers import REPO, VERSION, make_project, run_cli, seed_success, write_ids, write_tsv

MEMBERS = ["S2", "S1", "S3"]


class SourcesTest(unittest.TestCase):
    """Every key a builder can produce is in SOURCES/NESTED_SOURCES, and every key there is produced by some
    configuration: the listing's `set_by` column never drifts from inputs.py."""

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
        ]))
        self.cohort = cohorts.freeze(self.db, self.events, "C1", write_ids(self.tmp, MEMBERS))
        self.produced: dict[str, set] = {s: set() for s in STAGES}

    def tearDown(self):
        self.db.close()

    def gen(self, stage: str, subject_type: str, subject_id: str, mode: str = "standalone", **overrides) -> None:
        self.cfg.stage_inputs = {stage: overrides} if overrides else {}
        ctx = inputs.BuildContext(self.cfg, self.db, mode, VERSION, subject_type, subject_id, cohort=self.cohort)
        doc, _ = inputs.generate(ctx, stage)
        ns = f"ugc_wgw_{stage}."
        self.produced[stage] |= {k[len(ns):] for k in doc}
        self.cfg.stage_inputs = {}

    def seed(self, stage: str, subject_type: str, subject_id: str, mode: str = "standalone", attempt: int = 1, **extra) -> None:
        key = "sample_id" if subject_type == "sample" else "cohort_id"
        doc = {key: subject_id, **extra}
        if subject_type == "cohort":
            doc["sample_ids"] = MEMBERS
        seed_success(self.cfg, self.db, stage, subject_type, subject_id, doc, mode=mode, attempt=attempt)

    def test_sources_pin_the_builders(self):
        # per-sample stages, both DeepVariant flavours
        self.gen("singleton", "sample", "S1")
        self.gen("upstream", "sample", "S1", "joint")
        self.cfg.deepvariant, self.cfg.gpu_type = "parabricks", "a100"
        self.gen("singleton", "sample", "S1")
        self.gen("upstream", "sample", "S1", "joint")
        self.cfg.deepvariant, self.cfg.gpu_type = "cpu", ""
        # joint chain
        for sid in MEMBERS:
            self.seed("upstream", "sample", sid, "joint", inferred_sex="MALE", fail_reads=["/x"] if sid == "S1" else None)
        self.gen("cohort_call", "cohort", "C1", "joint")
        self.gen("cohort_call", "cohort", "C1", "joint", run_sawfish_joint_call=True)
        self.seed("cohort_call", "cohort", "C1", "joint")
        self.gen("downstream", "sample", "S1", "joint")   # with fail reads
        self.gen("downstream", "sample", "S2", "joint")
        for sid in MEMBERS:
            self.seed("downstream", "sample", sid, "joint", inferred_sex="MALE")
        self.gen("cohort_merge", "cohort", "C1", "joint")
        self.seed("cohort_merge", "cohort", "C1", "joint", run_glnexus=False)
        self.gen("cohort_freq", "cohort", "C1", "joint")
        # standalone chain
        for sid in MEMBERS:
            self.seed("singleton", "sample", sid, inferred_sex="MALE")
        self.gen("cohort_merge", "cohort", "C1")
        self.gen("cohort_merge", "cohort", "C1", merge_phased_small_variants=True)
        self.seed("cohort_merge", "cohort", "C1", attempt=2, run_glnexus=True)
        self.gen("cohort_freq", "cohort", "C1")
        # assembly, trio and sample mode
        self.gen("assembly", "sample", "K1", "assembly")
        self.gen("assembly", "sample", "S1", "assembly")
        for stage in STAGES:
            plain = {k for k in self.produced[stage] if "." not in k}
            dotted = {k for k in self.produced[stage] if "." in k}
            self.assertEqual(plain, set(stage_inputs.SOURCES[stage]), stage)
            self.assertEqual(dotted, set(stage_inputs.NESTED_SOURCES.get(stage, {})), stage)
        # the labels themselves: every driver/config key of SOURCES really is one of the three kinds
        for stage, table in stage_inputs.SOURCES.items():
            for name, (who, note) in table.items():
                self.assertIn(who, (stage_inputs.DRIVER, stage_inputs.CONFIG, stage_inputs.FREE), f"{stage}.{name}")
                if who == stage_inputs.FREE:
                    self.assertTrue(note.startswith("switch:"), f"{stage}.{name}: a FREE key in SOURCES is a switch the driver reads")


class ParserTest(unittest.TestCase):
    def test_specs_and_descriptions_of_the_real_entrypoints(self):
        specs = {s.name: s for s in declared_input_specs(wdl_path(REPO, "singleton"))}
        self.assertEqual((specs["use_alignment_chunking"].type, specs["use_alignment_chunking"].default,
                          specs["use_alignment_chunking"].required), ("Boolean", "true", False))
        self.assertEqual((specs["fail_reads"].type, specs["fail_reads"].default, specs["fail_reads"].required), ("Array[File]?", None, False))
        self.assertEqual((specs["hifi_reads"].type, specs["hifi_reads"].required), ("Array[File]", True))
        merge = {s.name: s for s in declared_input_specs(wdl_path(REPO, "cohort_merge"))}
        self.assertEqual(merge["sv_merge_method"].default, '"svx"')
        self.assertEqual(merge["backend"].default, '"HPC"')
        desc = input_descriptions(wdl_path(REPO, "cohort_merge"))
        self.assertEqual(desc["sv_merge_method"], "SV merge method (choices: svx, bcftools)")
        self.assertIn("GLnexus", desc["run_glnexus"])
        for stage in STAGES:   # every input of every entrypoint carries a parameter_meta description
            described = input_descriptions(wdl_path(REPO, stage))
            missing = [s.name for s in declared_input_specs(wdl_path(REPO, stage)) if not described.get(s.name)]
            self.assertEqual(missing, [], f"{stage}: inputs without parameter_meta description")


class ListingTest(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.cfg = make_project(self.tmp)
        self.proj = str(self.cfg.project_dir)

    @staticmethod
    def header(out: str) -> str:
        return next(line for line in out.splitlines() if line.startswith("input "))

    def rows(self, *args) -> tuple[list, str]:
        code, out, err = run_cli(["--project", self.proj, "stage-inputs", "--json", *args])
        self.assertEqual(code, 0, err)
        return json.loads(out), err

    def test_every_stage_is_described_and_required_inputs_are_filled(self):
        rows, err = self.rows()
        self.assertEqual(err, "")
        self.assertEqual([r["stage"] for r in rows][:1], ["singleton"])
        self.assertEqual(set(r["stage"] for r in rows), set(STAGES))
        self.assertEqual(sum(1 for r in rows if r["nested"]), 0)
        for r in rows:
            self.assertTrue(r["description"], r)
            self.assertIn(r["set_by"], (stage_inputs.DRIVER, stage_inputs.CONFIG, stage_inputs.FREE), r)
            if r["required"]:
                self.assertNotEqual(r["set_by"], stage_inputs.FREE, r)
        by = {(r["stage"], r["input"]): r for r in rows}
        self.assertEqual(by[("singleton", "use_gpu")]["set_by"], "config.json")
        self.assertEqual(by[("singleton", "use_gpu")]["from"], "deepvariant")
        self.assertEqual((by[("singleton", "use_alignment_chunking")]["set_by"], by[("singleton", "use_alignment_chunking")]["default"]),
                         ("stage_inputs", "true"))
        self.assertEqual(by[("singleton", "hifi_reads")]["from"], "sample sheet")
        self.assertTrue(by[("cohort_merge", "run_glnexus")]["from"].startswith("switch:"))
        self.assertEqual(by[("cohort_merge", "sv_merge_method")]["default"], '"svx"')
        self.assertIn("choices: svx, bcftools", by[("cohort_merge", "sv_merge_method")]["description"])
        self.assertNotIn("override", by[("cohort_merge", "sv_merge_method")])

    def test_stage_and_mode_filters_and_the_table(self):
        rows, _ = self.rows("--mode", "standalone")
        self.assertEqual(sorted(set(r["stage"] for r in rows)), ["cohort_freq", "cohort_merge", "singleton"])
        rows, _ = self.rows("--stage", "assembly")
        self.assertEqual(set(r["stage"] for r in rows), {"assembly"})
        self.assertIn("hifiasm_threads", [r["input"] for r in rows])
        code, out, err = run_cli(["--project", self.proj, "stage-inputs", "--stage", "assembly"])
        self.assertEqual(code, 0, err)
        self.assertIn("# 1 stage, 23 workflow inputs:", out)
        self.assertIn("--nested adds", out)
        self.assertIn("# assembly (ugc_wgw_assembly.wdl): 23 inputs", out)
        self.assertRegex(out, r"hifiasm_threads +Int +48 +stage_inputs")
        self.assertRegex(out, r"ugc_wgw_container_registry +String +- +config.json +ugc_wgw_container_registry")
        self.assertNotIn("override", self.header(out))
        with self.assertRaises(SystemExit) as cm:   # argparse: --stage and --mode are mutually exclusive
            run_cli(["--project", self.proj, "stage-inputs", "--stage", "assembly", "--mode", "joint"])
        self.assertEqual(cm.exception.code, 2)

    def test_overrides_are_shown_and_checked(self):
        self.cfg.stage_inputs = {
            "cohort_merge": {"sv_merge_method": "bcftools", "sample_ids": ["X"]},
            "ugc_wgw_cohort_call": {"glnexus_mem_gb": 200},
            "singleton": {"upstream.parabricks_deepvariant.run_parabricks_deepvariant.gpuCount": 3,
                          "nosuch.call.x": 1, "typo_input": 2},
            "trio": {"a": 1},
            "assembly": 7,
        }
        config.save(self.cfg)
        rows, err = self.rows()
        by = {(r["stage"], r["input"]): r for r in rows}
        self.assertEqual(by[("cohort_merge", "sv_merge_method")]["override"], "bcftools")
        self.assertEqual(by[("cohort_call", "glnexus_mem_gb")]["override"], 200)
        self.assertEqual(by[("cohort_merge", "sample_ids")]["override"], ["X"])
        self.assertNotIn("override", by[("cohort_merge", "svx_threads")])
        self.assertIn("stage_inputs key 'trio' names no stage", err)
        self.assertIn("stage_inputs.assembly must be an object", err)
        self.assertIn("stage_inputs.singleton.typo_input is not an input of ugc_wgw_singleton.wdl; submit would refuse it", err)
        self.assertIn("stage_inputs.singleton.nosuch.call.x: nosuch is not a call of ugc_wgw_singleton.wdl", err)
        self.assertIn("stage_inputs.cohort_merge.sample_ids overrides a value the driver fills (frozen cohort, in its order)", err)
        self.assertNotIn("gpuCount", err)   # a call-qualified key of a real call: miniwdl checks the rest at start
        code, out, err = run_cli(["--project", self.proj, "stage-inputs", "--stage", "cohort_merge"])
        self.assertEqual(code, 0)
        self.assertRegex(out, r'sv_merge_method +String +"svx" +stage_inputs .*"bcftools"')
        self.assertIn("override", self.header(out))   # the column appears when a row has one


def python_with_wdl() -> Path | None:
    for cand in (REPO / ".venv" / "bin" / "python3", Path(sys.executable)):
        if cand.exists() and subprocess.run([str(cand), "-c", "import WDL"], capture_output=True).returncode == 0:
            return cand
    return None


class NestedTest(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def test_needs_the_engine_python(self):
        cfg = make_project(self.tmp)   # the fake miniwdl has no python next to it
        code, out, err = run_cli(["--project", str(cfg.project_dir), "stage-inputs", "--nested", "--stage", "singleton"])
        self.assertEqual(code, 1)
        self.assertIn("engine's Python", err)
        self.assertIn(str(engine_python(cfg)), err)

    def test_nested_listing_from_miniwdl(self):
        python = python_with_wdl()
        if python is None:
            self.skipTest("no Python with miniwdl's WDL package (.venv or sys.executable)")
        cfg = make_project(self.tmp, venv=python.parents[1])
        self.assertTrue(engine_python(cfg).exists(), engine_python(cfg))
        proj = str(cfg.project_dir)
        code, out, err = run_cli(["--project", proj, "stage-inputs", "--nested", "--stage", "singleton", "--json"])
        self.assertEqual(code, 0, err)
        self.assertEqual(err, "", err)   # the driver's parser and miniwdl agree on the entrypoint's own inputs
        rows = json.loads(out)
        nested = {r["input"]: r for r in rows if r["nested"]}
        self.assertGreater(len(nested), 50)
        threads = nested["upstream.pbmm2.pbmm2_align_wgs.threads"]
        self.assertEqual((threads["type"], threads["default"], threads["required"], threads["set_by"], threads["description"]),
                         ("Int", "32", False, "stage_inputs", "Number of threads to use"))
        gpus = nested["upstream.parabricks_deepvariant.run_parabricks_deepvariant.gpuCount"]
        self.assertEqual((gpus["default"], gpus["set_by"], gpus["from"]), ("4", "config.json", "parabricks_gpus, parabricks flavour only"))
        self.assertTrue(all(r["description"] for r in rows if not r["nested"]))
        # the project's call-qualified overrides are attached, or refused when miniwdl knows no such input
        cfg.stage_inputs = {"singleton": {"upstream.pbmm2.pbmm2_align_wgs.threads": 24, "upstream.nosuch": 1}}
        config.save(cfg)
        code, out, err = run_cli(["--project", proj, "stage-inputs", "--nested", "--stage", "singleton", "--json"])
        self.assertEqual(code, 0, err)
        rows = {r["input"]: r for r in json.loads(out)}
        self.assertEqual(rows["upstream.pbmm2.pbmm2_align_wgs.threads"]["override"], 24)
        self.assertIn("stage_inputs.singleton.upstream.nosuch: miniwdl lists no such input of ugc_wgw_singleton.wdl", err)
        code, out, err = run_cli(["--project", proj, "stage-inputs", "--nested", "--mode", "assembly"])
        self.assertEqual(code, 0, err)
        self.assertIn("call-qualified inputs of the tasks inside", out)
        self.assertRegex(out, r"# assembly \(ugc_wgw_assembly.wdl\): 23 inputs, \d+ call-qualified")


if __name__ == "__main__":
    unittest.main()
