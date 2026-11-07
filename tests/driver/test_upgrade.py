"""`ugc-wgw upgrade --install`: a project moves to another installed version; the driver refuses a project whose
engine paths are another version's (the 864 cohort's project at 0.12.6 under the 0.15.0 driver, 2026-11-07)."""
from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from bin.ugc_wgw import config
from bin.ugc_wgw.util import UgcError

from .helpers import REPO, run_cli

REPO_VERSION = (REPO / "VERSION").read_text().strip()


def fake_install(root: Path, version: str) -> Path:
    """`<root>/versions/<version>` with what derive_from_install and installed_version_dir look for."""
    v = root / "versions" / version
    (v / "code" / "references").mkdir(parents=True)
    (v / "code" / "VERSION").write_text(version + "\n")
    (v / "code" / "references" / "ugc_wgw_ref_map.GRCh38_GIABv3.tsv").write_text("name\tGRCh38_GIABv3\n")
    (v / "venv" / "bin").mkdir(parents=True)
    (v / "venv" / "bin" / "miniwdl").write_text("#!/bin/sh\n")
    (v / "miniwdl.cfg").write_text("[singularity]\nimage_cache = " + str(v / "sif") + "\n")
    (v / "manifest.json").write_text(json.dumps({"version": version}))
    return v


class UpgradeTest(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.v1 = fake_install(self.tmp, "0.1.0")
        self.v2 = fake_install(self.tmp, REPO_VERSION)
        (self.tmp / "current").symlink_to(self.v2, target_is_directory=True)
        self.cfg = config.init_project(self.tmp / "proj", code=self.v1 / "code", miniwdl=self.v1 / "venv" / "bin" / "miniwdl",
                                       cfg=self.v1 / "miniwdl.cfg", venv=self.v1 / "venv",
                                       ref_map=self.v1 / "code" / "references" / "ugc_wgw_ref_map.GRCh38_GIABv3.tsv")
        self.proj = str(self.cfg.project_dir)

    def test_refused_under_another_driver(self):
        # this driver (the repository) is not 0.1.0: every project command stops with the remedy, nothing is opened
        code, _, err = run_cli(["--project", self.proj, "samples", "list"])
        self.assertEqual(code, 1)
        self.assertIn(f"this driver is ugc-wgw {REPO_VERSION}", err)
        self.assertIn("code_dir is version 0.1.0", err)
        self.assertIn("upgrade --install", err)

    def test_upgrade_moves_every_engine_path(self):
        code, _, err = run_cli(["--project", self.proj, "upgrade", "--install", str(self.tmp / "current")])
        self.assertEqual(code, 0, err)
        self.assertIn(f"ugc-wgw 0.1.0 -> {REPO_VERSION}", err)
        cfg = config.load(Path(self.proj))
        self.assertEqual(cfg.code_dir, self.v2 / "code")                       # `current` resolved to the version dir
        self.assertEqual(cfg.miniwdl, self.v2 / "venv" / "bin" / "miniwdl")
        self.assertEqual(cfg.miniwdl_cfg, self.v2 / "miniwdl.cfg")
        self.assertEqual(cfg.venv_dir, self.v2 / "venv")
        self.assertEqual(cfg.ref_map_file, self.v2 / "code" / "references" / "ugc_wgw_ref_map.GRCh38_GIABv3.tsv")
        self.assertTrue((cfg.config_path.with_name("config.json.0.1.0")).exists())
        self.assertEqual(json.loads(cfg.config_path.with_name("config.json.0.1.0").read_text())["code_dir"], str(self.v1 / "code"))
        # the project now opens under this driver; the upgrade is recorded
        code, out, err = run_cli(["--project", self.proj, "samples", "list"])
        self.assertEqual(code, 0, err)
        events = (cfg.logs_dir / "events.jsonl").read_text()
        self.assertIn('"project.upgraded"', events)
        self.assertIn('"old_version": "0.1.0"', events)
        # a second upgrade to the same version changes nothing
        code, _, err = run_cli(["--project", self.proj, "upgrade", "--install", str(self.v2)])
        self.assertEqual(code, 0, err)
        self.assertIn("nothing changed", err)

    def test_reference_map_outside_the_code_dir_is_kept(self):
        own = self.tmp / "site_ref_map.tsv"
        own.write_text("name\tGRCh38_GIABv3\n")
        cfg = config.load(Path(self.proj))
        cfg.ref_map_file = own
        config.save(cfg)
        old, new, changes = config.upgrade_project(config.load(Path(self.proj)), self.v2)
        self.assertEqual((old, new), ("0.1.0", REPO_VERSION))
        self.assertNotIn("ref_map_file", changes)
        self.assertEqual(config.load(Path(self.proj)).ref_map_file, own)

    def test_missing_reference_map_at_the_new_version(self):
        (self.v2 / "code" / "references" / "ugc_wgw_ref_map.GRCh38_GIABv3.tsv").unlink()
        with self.assertRaises(UgcError) as cm:
            config.upgrade_project(config.load(Path(self.proj)), self.v2)
        self.assertIn("does not exist", str(cm.exception))
        self.assertEqual(config.load(Path(self.proj)).code_dir, self.v1 / "code")   # nothing written

    def test_mixed_engine_paths_are_refused(self):
        # code_dir at this driver's version, miniwdl.cfg still another version's: miniwdl would resolve images from the
        # wrong SIF directory and try the registry (the cluster has no network)
        cfg = config.load(Path(self.proj))
        cfg.code_dir, cfg.miniwdl, cfg.venv_dir = self.v2 / "code", self.v2 / "venv" / "bin" / "miniwdl", self.v2 / "venv"
        config.save(cfg)
        code, _, err = run_cli(["--project", self.proj, "samples", "list"])
        self.assertEqual(code, 1)
        self.assertIn("miniwdl_cfg = " + str(self.v1 / "miniwdl.cfg"), err)
        self.assertIn("every engine path must belong to one installed version", err)
        # a repository checkout as code_dir (no manifest.json next to it) is never checked for stray paths
        cfg = config.load(Path(self.proj))
        cfg.code_dir = REPO
        config.save(cfg)
        code, _, err = run_cli(["--project", self.proj, "samples", "list"])
        self.assertEqual(code, 0, err)


if __name__ == "__main__":
    unittest.main()
