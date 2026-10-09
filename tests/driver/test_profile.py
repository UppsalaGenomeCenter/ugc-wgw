"""Site profiles: `ugc-wgw init --profile` from the shipped examples, an install prefix or a file; flags win over
the profile; the record in config.json; the listing names the profile as the source of its overrides."""
import json
import os
import re
import tempfile
import unittest
from pathlib import Path

from ugc_wgw import config
from ugc_wgw.util import UgcError

from .helpers import FAKE, REPO, run_cli


class ProfileTest(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.ref = self.tmp / "ugc_wgw_ref_map.GRCh38_GIABv3.tsv"
        self.ref.write_text("name\tGRCh38_GIABv3\nscatter_regions\t/dev/null\n")
        self.cfg_file = self.tmp / "miniwdl.cfg"
        self.cfg_file.write_text("[scheduler]\n")
        self.n = 0

    def init(self, *extra: str) -> tuple[int, str, config.Config | None]:
        self.n += 1
        proj = self.tmp / f"proj{self.n}"
        code, _, err = run_cli(["init", str(proj), "--code", str(REPO), "--miniwdl", str(FAKE), "--cfg", str(self.cfg_file),
                                "--ref-map", str(self.ref), *extra])
        return code, err, config.load(proj) if code == 0 else None

    def test_shipped_profiles_apply_and_their_keys_are_real_inputs(self):
        for name in ("cpu", "parabricks"):
            code, err, cfg = self.init("--profile", name)
            self.assertEqual(code, 0, err)
            doc = json.loads((REPO / "backends" / "hpc" / "profiles" / f"{name}.json").read_text())
            self.assertEqual(cfg.stage_inputs, doc["stage_inputs"])
            self.assertEqual((cfg.max_inflight, cfg.deepvariant), (20, doc["deepvariant"]))
            self.assertEqual(cfg.gpu_type, doc.get("gpu_type", ""))
            self.assertEqual(cfg.profile["name"], name)
            self.assertTrue(cfg.profile["path"].endswith(f"backends/hpc/profiles/{name}.json"))
            self.assertRegex(cfg.profile["sha256"], r"^[0-9a-f]{64}$")
            self.assertIn("stage_inputs.singleton.use_alignment_chunking", cfg.profile["applied"])
            self.assertIn("max_inflight", cfg.profile["applied"])
            self.assertNotIn("description", cfg.profile["applied"])
            self.assertIn(f"profile {name} (", err)
            # every key the profile sets is an input submit would accept: the listing has nothing to warn about
            code, out, err = run_cli(["--project", str(cfg.project_dir), "stage-inputs", "--json"])
            self.assertEqual(code, 0, err)
            self.assertEqual(err, "", err)
            rows = {(r["stage"], r["input"]): r for r in json.loads(out)}
            row = rows[("singleton", "use_alignment_chunking")]
            self.assertEqual((row["override"], row["override_source"]), (False, f"profile {name}"))
            code, out, _ = run_cli(["--project", str(cfg.project_dir), "stage-inputs", "--stage", "upstream"])
            self.assertRegex(out, r"use_alignment_chunking +Boolean +true +stage_inputs .*false \(profile\)")
            # the report's provenance names it
            report = self.tmp / f"report-{name}.html"
            code, _, err = run_cli(["--project", str(cfg.project_dir), "report", "--out", str(report)])
            self.assertEqual(code, 0, err)
            self.assertIn(f"{name} ({cfg.profile['sha256'][:12]}) from", report.read_text())

    def test_flags_win_over_the_profile(self):
        code, err, cfg = self.init("--profile", "parabricks", "--deepvariant", "cpu", "--max-inflight", "3", "--gpu-type", "a100")
        self.assertEqual(code, 0, err)
        self.assertEqual((cfg.deepvariant, cfg.max_inflight, cfg.gpu_type, cfg.parabricks_gpus), ("cpu", 3, "a100", 4))
        self.assertEqual(cfg.stage_inputs["singleton"]["upstream.pbmm2.pbmm2_align_wgs.threads"], 24)
        code, err, cfg = self.init()
        self.assertEqual(code, 0, err)
        self.assertEqual((cfg.deepvariant, cfg.max_inflight, cfg.gpu_type, cfg.parabricks_gpus, cfg.stage_inputs, cfg.profile),
                         ("cpu", 4, "", 4, {}, {}))

    def test_profile_by_file_and_install_root(self):
        mine = self.tmp / "mine.json"
        mine.write_text(json.dumps({"description": "x", "prices": {"cpu_hour": 0.14, "currency": "SEK"}, "poll_interval": 15,
                                    "summary_thresholds": {"depth_mean_min": 18}, "assembly_use_parents": False,
                                    "stage_inputs": {"ugc_wgw_cohort_merge": {"sv_merge_method": "bcftools"}}}))
        code, err, cfg = self.init("--profile", str(mine))
        self.assertEqual(code, 0, err)
        self.assertEqual((cfg.prices, cfg.poll_interval, cfg.summary_thresholds, cfg.assembly_use_parents),
                         ({"cpu_hour": 0.14, "currency": "SEK"}, 15.0, {"depth_mean_min": 18}, False))
        self.assertEqual(cfg.profile["applied"], ["assembly_use_parents", "poll_interval", "prices",
                                                  "stage_inputs.ugc_wgw_cohort_merge.sv_merge_method", "summary_thresholds"])
        code, out, _ = run_cli(["--project", str(cfg.project_dir), "stage-inputs", "--stage", "cohort_merge", "--json"])
        row = {r["input"]: r for r in json.loads(out)}["sv_merge_method"]
        self.assertEqual((row["override"], row["override_source"]), ("bcftools", "profile mine"))
        # an install prefix: <root>/profiles/<name>.json wins over the code's examples
        root = self.tmp / "prefix"
        (root / "versions" / "0.0.0").mkdir(parents=True)
        os.symlink("versions/0.0.0", root / "current")
        (root / "profiles").mkdir()
        (root / "profiles" / "cpu.json").write_text(json.dumps({"max_inflight": 7}))
        self.assertEqual(config.install_root(root / "current"), root)
        self.assertEqual(config.install_root(root / "versions" / "0.0.0"), root)
        self.assertIsNone(config.install_root(self.tmp))
        self.assertEqual(config.find_profile("cpu", install=root / "current", code=REPO), (root / "profiles" / "cpu.json").resolve())
        self.assertEqual(config.find_profile("parabricks", install=root / "current", code=REPO),
                         (REPO / "backends" / "hpc" / "profiles" / "parabricks.json").resolve())
        with self.assertRaises(UgcError) as cm:
            config.find_profile("nope", install=root / "current", code=REPO)
        self.assertIn(str(root / "profiles" / "nope.json"), str(cm.exception))
        self.assertIn("backends/hpc/profiles/nope.json", str(cm.exception))
        with self.assertRaises(UgcError):
            config.find_profile(str(self.tmp / "missing.json"), install=None, code=REPO)

    def test_bad_profiles_are_refused(self):
        cases = {
            "unknown.json": ({"bogus": 1}, "unknown key 'bogus'"),
            "inputs.json": ({"stage_inputs": {"singleton": 7}}, "stage_inputs.singleton must be an object"),
            "dv.json": ({"deepvariant": "tpu"}, "deepvariant must be one of"),
            "inflight.json": ({"max_inflight": 0}, "max_inflight must be a positive integer"),
            "prices.json": ({"prices": [1]}, "prices must be an object"),
            "list.json": ([1, 2], "not a JSON object"),
        }
        for name, (doc, msg) in cases.items():
            path = self.tmp / name
            path.write_text(json.dumps(doc))
            with self.assertRaises(UgcError) as cm:
                config.load_profile(path)
            self.assertIn(msg, str(cm.exception), name)
        (self.tmp / "broken.json").write_text("{")
        with self.assertRaises(UgcError):
            config.load_profile(self.tmp / "broken.json")
        code, err, _ = self.init("--profile", str(self.tmp / "unknown.json"))
        self.assertEqual(code, 1)
        self.assertIn("unknown key", err)


    def test_fs_warnings(self):
        code, err, cfg = self.init()
        self.assertEqual(code, 0, err)
        same_fs = config.device_of(cfg.results_dir) == config.device_of(cfg.code_dir)   # depends on the machine
        self.assertEqual("different file systems" in err, not same_fs, err)
        self.assertEqual(len(config.fs_warnings(cfg)), 0 if same_fs else 1)
        devs = {cfg.results_dir: 1, cfg.code_dir: 2}
        warnings = config.fs_warnings(cfg, device=lambda p: devs.get(p))
        self.assertEqual(len(warnings), 1)
        self.assertIn(f"results {cfg.results_dir} and the install {cfg.code_dir} are on different file systems", warnings[0])
        self.assertIn("[call_cache] get = false", warnings[0])
        self.assertEqual(config.fs_warnings(cfg, device=lambda p: None), [])   # unknown devices: no warning


if __name__ == "__main__":
    unittest.main()
