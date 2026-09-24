"""Tests for `ugc-wgw resources` and the resource-policy preflight of submit (bin/ugc_wgw/resources.py)."""
from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from ugc_wgw import resources
from ugc_wgw.db import DB

from ugc_wgw import config

from .helpers import VERSION, make_project, run_cli, write_tsv

G = 1024 ** 3
HEADER = "task\tcpu\tmemory\ttime\tpartition\tconstraint\n"
DEFAULTS = {"docker": "ubuntu:20.04", "maxRetries": 1, "time_minutes": 4320, "slurm_partition": "core"}


def write_cfg(cfg, *, policy_path=None, cpu_max=48, memory_max="386547056640", extra_args="", defaults=DEFAULTS,
              ugc_wgw_section=True, backend="singularity", run_options=None):
    text = f"[scheduler]\ncontainer_backend = {backend}\n\n"
    if run_options is not None:
        text += f"[singularity]\nrun_options = {json.dumps(run_options)}\n\n"
    text += "[task_runtime]\n"
    text += f"cpu_max = {cpu_max}\nmemory_max = {memory_max}\ndefaults = {json.dumps(defaults)}\n\n"
    text += f'[slurm]\nextra_args = "{extra_args}"\n'
    if ugc_wgw_section:
        text += f"\n[ugc_wgw]\nresources = {policy_path or ''}\n"
    cfg.miniwdl_cfg.write_text(text)


class ResourcesCommandTest(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.cfg = make_project(self.tmp)
        self.proj = str(self.cfg.project_dir)
        self.policy = self.tmp / "resources.tsv"
        self.policy.write_text(HEADER + "ugc_wgw_hifiasm_assemble\t-\t-\t2-00:00:00\tfat\t-\ndeepvariant_*\t-\t-\t12:00:00\t-\t-\n")
        write_cfg(self.cfg, policy_path=self.policy)

    def rows(self, *args):
        code, out, err = run_cli(["--project", self.proj, "resources", "--json", *args])
        self.assertEqual(code, 0, err)
        return {r["task"]: r for r in json.loads(out)}, err

    def test_declared_clamp_and_policy(self):
        rows, err = self.rows()
        self.assertEqual(err, "")
        dv = rows["deepvariant_call_variants_cpu"]
        self.assertEqual((dv["declared"], dv["cpu"], dv["memory"]), ("64 / 28G", "48 (clamp)", "28G"))
        self.assertEqual((dv["time"], dv["partition"], dv["stages"]), ("12:00:00 (policy)", "core", "singleton,upstream"))
        self.assertEqual(dv["rules"], f"{self.policy}:3")
        asm = rows["ugc_wgw_hifiasm_assemble"]
        self.assertEqual((asm["cpu"], asm["memory"], asm["time"], asm["partition"]),
                         ("48", "288G", "2-00:00:00 (policy)", "fat (policy)"))
        mos = rows["mosdepth"]
        self.assertEqual((mos["cpu"], mos["memory"], mos["time"], mos["partition"], mos["rules"]),
                         ("4", "8G", "3-00:00:00", "core", ""))  # mosdepth 0.3.14 task since v4.0.0
        self.assertEqual(rows["glnexus"]["memory"], "-")
        self.assertIn("scales with N", rows["glnexus"]["notes"])

    def test_filters_and_text_output(self):
        changed, _ = self.rows("--changed")
        self.assertIn("deepvariant_call_variants_cpu", changed)
        self.assertIn("deepvariant_make_examples", changed)   # the glob row sets its time
        self.assertIn("ugc_wgw_hifiasm_assemble", changed)
        self.assertNotIn("mosdepth", changed)
        asm, _ = self.rows("--stage", "assembly")
        self.assertIn("ugc_wgw_yak_count", asm)
        self.assertNotIn("mosdepth", asm)
        code, out, err = run_cli(["--project", self.proj, "resources", "--changed"])
        self.assertEqual(code, 0, err)
        self.assertIn(f"# resources: cpu_max=48 memory_max=360G policy={self.policy} (2 rows)", out)
        self.assertIn("task", out.splitlines()[2])
        self.assertIn("48 (clamp)", out)

    def test_warnings(self):
        self.policy.write_text(HEADER + "deepvarient_*\t2\t-\t-\t-\t-\nugc_wgw_hifiasm_assemble\t-\t-\t-\tfat\t-\n")
        write_cfg(self.cfg, policy_path=self.policy, extra_args="--partition core --account a1")
        _, err = self.rows()
        self.assertIn("matches no task in the inventory", err)
        self.assertIn("sbatch keeps the last one", err)
        write_cfg(self.cfg, policy_path=self.tmp / "nope.tsv")
        _, err = self.rows()
        self.assertIn("resource policy file missing", err)
        write_cfg(self.cfg, ugc_wgw_section=False)
        rows, err = self.rows()
        self.assertIn("no [ugc_wgw] resources entry", err)
        self.assertEqual(rows["deepvariant_call_variants_cpu"]["cpu"], "48 (clamp)")
        write_cfg(self.cfg, policy_path=self.policy, cpu_max=0, memory_max="0")
        rows, _ = self.rows()
        self.assertEqual(rows["deepvariant_call_variants_cpu"]["cpu"], "64")

    def test_gpu_rows_and_warnings(self):
        rows, err = self.rows()
        self.assertEqual((rows["deepvariant_call_variants_gpu"]["gpus"], rows["run_parabricks_deepvariant"]["gpus"],
                          rows["mosdepth"]["gpus"]), ("1", "4", "-"))
        self.assertEqual(err, "")   # deepvariant=cpu: nothing to warn about
        code, out, _ = run_cli(["--project", self.proj, "resources"])
        self.assertIn("deepvariant=cpu", out.splitlines()[0])
        self.assertIn("gpus", out.splitlines()[2])
        self.cfg.deepvariant, self.cfg.gpu_type, self.cfg.parabricks_gpus = "parabricks", "a100", 2
        config.save(self.cfg)
        rows, err = self.rows()
        self.assertEqual((rows["deepvariant_call_variants_gpu"]["gpus"], rows["run_parabricks_deepvariant"]["gpus"]), ("a100:1", "a100:2"))
        self.assertIn("run_options has no --nv", err)
        self.assertNotIn("no GPU partition", err)          # plain singularity backend: no partitions at all
        code, out, _ = run_cli(["--project", self.proj, "resources"])
        self.assertIn("deepvariant=parabricks gpu_type=a100", out.splitlines()[0])
        write_cfg(self.cfg, policy_path=self.policy, backend="slurm_singularity", run_options=["--containall", "--nv"])
        rows, err = self.rows()
        self.assertNotIn("--nv", err)
        self.assertIn("no GPU partition", err)
        self.assertIn("run_parabricks_deepvariant", err)
        self.assertEqual(rows["run_parabricks_deepvariant"]["partition"], "core")
        write_cfg(self.cfg, policy_path=self.policy, backend="slurm_singularity", run_options=["--nv"],
                  defaults={**DEFAULTS, "slurm_partition_gpu": "gpu"})
        rows, err = self.rows()
        self.assertEqual(err, "")
        self.assertEqual((rows["run_parabricks_deepvariant"]["partition"], rows["deepvariant_call_variants_gpu"]["partition"],
                          rows["mosdepth"]["partition"]), ("gpu", "gpu", "core"))
        self.policy.write_text(HEADER + "run_parabricks_deepvariant\t-\t-\t-\tgpuq\t-\n")
        write_cfg(self.cfg, policy_path=self.policy, backend="slurm_singularity", run_options=["--nv"])
        rows, err = self.rows()
        self.assertEqual(err, "")   # the row names the partition
        self.assertEqual(rows["run_parabricks_deepvariant"]["partition"], "gpuq (policy)")
        self.cfg.deepvariant = "gpu"
        config.save(self.cfg)
        rows, err = self.rows()
        self.assertIn("deepvariant_call_variants_gpu", err)   # the selected task has no partition

    def test_malformed_policy_is_refused_by_resources_and_submit(self):
        self.policy.write_text(HEADER + "mosdepth\tlots\t-\t-\t-\t-\n")
        code, out, err = run_cli(["--project", self.proj, "resources"])
        self.assertEqual(code, 1)
        self.assertIn(f"{self.policy}:2: column cpu", err)
        run_cli(["--project", self.proj, "samples", "add", str(write_tsv(self.tmp, [{"sample_id": "S1", "hifi_reads": "a.bam"}]))])
        code, out, err = run_cli(["--project", self.proj, "submit", "--mode", "standalone", "--samples", "S1", "--dry-run"])
        self.assertEqual(code, 1)
        self.assertIn(f"{self.policy}:2: column cpu", err)
        self.assertNotIn("ugc_wgw_singleton.wdl", out)
        code, out, err = run_cli(["--project", self.proj, "submit", "--mode", "standalone", "--samples", "S1"])
        self.assertEqual(code, 1)
        self.assertIn("column cpu", err)

    def test_dry_run_and_submit_record_the_policy(self):
        run_cli(["--project", self.proj, "samples", "add", str(write_tsv(self.tmp, [{"sample_id": "S1", "hifi_reads": "a.bam"}]))])
        code, out, err = run_cli(["--project", self.proj, "submit", "--mode", "standalone", "--samples", "S1", "--dry-run"])
        self.assertEqual(code, 0, err)
        self.assertIn(f"# resources: cpu_max=48 memory_max=360G policy={self.policy} (2 rows)", out)
        code, out, err = run_cli(["--project", self.proj, "submit", "--mode", "standalone", "--samples", "S1"])
        self.assertEqual(code, 0, err)
        self.assertIn("no ugc_wgw_resources task plugin", err)   # the fake engine has none
        db = DB(self.cfg.db_path)
        run = db.latest_run("sample", "S1", "singleton", VERSION)
        db.close()
        self.assertEqual(run.status, "success")
        m = json.loads((run.run_path / "run_manifest.json").read_text())
        self.assertEqual(m["resource_policy"]["path"], str(self.policy))
        self.assertEqual(len(m["resource_policy"]["sha256"]), 64)
        ref = resources.policy_ref(self.cfg)
        self.assertEqual(ref, m["resource_policy"])


class EngineCfgTest(unittest.TestCase):
    def test_reads_rendered_and_hand_edited_values(self):
        tmp = Path(tempfile.mkdtemp())
        p = tmp / "miniwdl.cfg"
        p.write_text('[task_runtime]\ncpu_max = 24\nmemory_max = 370G\ndefaults = {"docker": "x", "time_minutes": 60}\n'
                     '[slurm]\nextra_args = "--qos a"\n[ugc_wgw]\nresources = "/x/y.tsv"\n')
        eng = resources.read_engine_cfg(p)
        self.assertEqual((eng.cpu_max, eng.memory_max, eng.policy_path, eng.extra_args), (24, 370 * 10**9, "/x/y.tsv", "--qos a"))
        self.assertEqual(eng.defaults["time_minutes"], 60)
        p.write_text("[task_runtime]\nmemory_max = -1\n")
        self.assertEqual(resources.read_engine_cfg(p).memory_max, -1)
        self.assertEqual(resources.read_engine_cfg(tmp / "missing.cfg").cpu_max, 0)
        p.write_text("[task_runtime]\nmemory_max = lots\n")
        with self.assertRaises(Exception):
            resources.read_engine_cfg(p)


if __name__ == "__main__":
    unittest.main()
