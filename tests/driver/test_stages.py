import json
import subprocess
import unittest

from ugc_wgw.stages import MODES, STAGES, declared_calls, declared_inputs, mode_stages, wdl_path
from ugc_wgw.util import UgcError

from .helpers import REPO


class StagesTest(unittest.TestCase):
    def test_modes_reference_known_stages(self):
        for mode, stages in MODES.items():
            self.assertEqual(mode_stages(mode), stages)
            for s in stages:
                self.assertIn(s, STAGES)
        self.assertEqual(mode_stages("assembly"), ("assembly",))
        with self.assertRaises(UgcError):
            mode_stages("trio")

    def test_declared_inputs_of_real_entrypoints(self):
        expect = {
            "singleton": {"ugc_wgw_version": True, "sample_id": True, "hifi_reads": True, "fail_reads": False,
                          "ref_map_file": True, "use_alignment_chunking": False, "backend": False, "preemptible": False},
            "upstream": {"run_sawfish_call": False, "use_alignment_chunking": False, "hifi_reads": True},
            "cohort_call": {"cohort_id": True, "sample_ids": True, "gvcfs": True, "ref_map_file": True,
                            "split_keep_homref": False, "discover_tars": False},
            "downstream": {"sex": True, "small_variant_vcf": True, "sv_vcf": True, "upstream_msg": False,
                           "stat_depth_mean": False, "aligned_fail_reads": False},
            "cohort_merge": {"sv_vcfs": True, "trgt_vcfs": True, "gvcfs": False, "ugc_wgw_container_registry": True,
                             "run_glnexus": False, "svx_mem_gb": False},
            "cohort_freq": {"cohort_id": True, "sample_ids": True, "sample_sexes": True, "sample_sex_sources": False,
                            "small_variant_vcf": False, "sv_vcf": True, "sv_vcf_index": True, "freq_mem_gb": False,
                            "ref_map_file": True},
            "assembly": {"ugc_wgw_version": True, "sample_id": True, "hifi_reads": True, "ref_map_file": True,
                         "ugc_wgw_container_registry": True,
                         "father_id": False, "mother_id": False, "father_hifi_reads": False, "mother_hifi_reads": False,
                         "hifiasm_extra_params": False, "hifiasm_mem_gb": False, "trio_low_depth_gb": False},
        }
        for stage, subset in expect.items():
            declared = declared_inputs(wdl_path(REPO, stage))
            for name, required in subset.items():
                self.assertIn(name, declared, f"{stage}: {name}")
                self.assertEqual(declared[name], required, f"{stage}: {name} required={required}")

    def test_declared_calls_of_real_entrypoints(self):
        calls = declared_calls(wdl_path(REPO, "singleton"))
        self.assertTrue({"upstream", "downstream", "filter_messages", "backend_configuration"} <= calls, calls)
        self.assertTrue({"upstream", "ugc_wgw_manifest_write"} <= declared_calls(wdl_path(REPO, "upstream")))
        for stage in STAGES:
            self.assertTrue(declared_calls(wdl_path(REPO, stage)), stage)

    def test_against_miniwdl_input_template(self):
        miniwdl = REPO / ".venv" / "bin" / "miniwdl"
        if not miniwdl.exists():
            self.skipTest("miniwdl not installed in .venv")
        for stage in STAGES:
            wdl = wdl_path(REPO, stage)
            res = subprocess.run([str(miniwdl), "input_template", str(wdl)], capture_output=True, text=True, check=True)
            template = json.loads(res.stdout)
            ns = f"ugc_wgw_{stage}."
            declared = declared_inputs(wdl)
            tpl_names = {k[len(ns):] for k in template if k.startswith(ns)}
            required = {n for n, r in declared.items() if r}
            self.assertTrue(required <= tpl_names, f"{stage}: required inputs missing from input_template: {required - tpl_names}")
            self.assertTrue(tpl_names <= set(declared), f"{stage}: template has undeclared inputs: {tpl_names - set(declared)}")


if __name__ == "__main__":
    unittest.main()
