"""The ugc_slurm_singularity backend and its exec script: node-local scratch bound into every task's container."""
from __future__ import annotations

import logging
import os
import pathlib
import stat
import subprocess
import tempfile
import unittest

from ugc_wgw_miniwdl import scratch

SCRIPT = pathlib.Path(scratch.EXEC_SCRIPT)


def fake_apptainer(tmp: pathlib.Path) -> pathlib.Path:
    """An `apptainer` that prints its arguments one per line and exits 0."""
    exe = tmp / "apptainer"
    exe.write_text('#!/usr/bin/env bash\nfor a in "$@"; do printf "%s\\n" "$a"; done\n')
    exe.chmod(exe.stat().st_mode | stat.S_IXUSR)
    return exe


def run_script(args: list[str], env: dict[str, str]) -> subprocess.CompletedProcess:
    return subprocess.run(["bash", str(SCRIPT), *args], env={"PATH": os.environ["PATH"], **env},
                          capture_output=True, text=True)


class ExecScriptTest(unittest.TestCase):
    def setUp(self):
        self.tmp = pathlib.Path(tempfile.mkdtemp())
        self.app = fake_apptainer(self.tmp)
        self.scratch = self.tmp / "job-scratch"
        self.scratch.mkdir()
        self.cmd = [str(self.app), "exec", "--pwd", "/mnt/miniwdl_task_container/work", "--containall", "--no-mount",
                    "hostfs", "--bind", "/data/a.bam:/mnt/miniwdl_task_container/work/_miniwdl_inputs/0/a.bam:ro",
                    "/sif/pbmm2.sif", "/bin/bash", "-c", "exec pbmm2 align"]

    def test_script_is_shipped_and_executable(self):
        self.assertTrue(SCRIPT.is_file(), SCRIPT)
        self.assertTrue(SCRIPT.stat().st_mode & stat.S_IXUSR)
        self.assertTrue(SCRIPT.read_text().startswith("#!/usr/bin/env bash\n"))

    def test_bind_and_env_inserted_after_exec(self):
        res = run_script(["TMPDIR", "/mnt/scratch", "--", *self.cmd], {"TMPDIR": str(self.scratch)})
        self.assertEqual(res.returncode, 0, res.stderr)
        got = res.stdout.splitlines()                    # the fake prints its arguments, without its own name
        self.assertEqual(got[:5], ["exec", "--bind", f"{self.scratch}:/mnt/scratch", "--env", "TMPDIR=/mnt/scratch"])
        self.assertEqual(got[5:], self.cmd[2:])          # the rest untouched, including the literal "exec" in the command
        self.assertEqual(res.stderr.strip(), f"[ugc-wgw] scratch: {self.scratch} -> /mnt/scratch (TMPDIR)")

    def test_verbose_before_exec_is_kept(self):
        cmd = [self.cmd[0], "--verbose", *self.cmd[1:]]
        res = run_script(["SNIC_TMP", "/scratch", "--", *cmd], {"SNIC_TMP": str(self.scratch)})
        self.assertEqual(res.returncode, 0, res.stderr)
        got = res.stdout.splitlines()
        self.assertEqual(got[:4], ["--verbose", "exec", "--bind", f"{self.scratch}:/scratch"])

    def test_unset_or_unusable_variable_runs_the_command_unchanged(self):
        for env, reason in (({}, "TMPDIR is unset"),
                            ({"TMPDIR": str(self.tmp / "missing")}, "is not a writable directory"),
                            ({"TMPDIR": str(self.tmp / "file")}, "is not a writable directory")):
            (self.tmp / "file").write_text("x")
            res = run_script(["TMPDIR", "/mnt/scratch", "--", *self.cmd], env)
            self.assertEqual(res.returncode, 0, res.stderr)
            self.assertEqual(res.stdout.splitlines(), self.cmd[1:], env)
            self.assertIn("[ugc-wgw] scratch: none", res.stderr)
            self.assertIn(reason, res.stderr)

    def test_tmpfs_is_refused(self):
        shm = pathlib.Path("/dev/shm")
        if not shm.is_dir() or not os.access(shm, os.W_OK):
            self.skipTest("no writable /dev/shm")
        fstype = subprocess.run(["stat", "-f", "-c", "%T", str(shm)], capture_output=True, text=True).stdout.strip()
        if fstype != "tmpfs":
            self.skipTest(f"/dev/shm is {fstype}, not tmpfs")
        res = run_script(["TMPDIR", "/mnt/scratch", "--", *self.cmd], {"TMPDIR": str(shm)})
        self.assertEqual(res.returncode, 0, res.stderr)
        self.assertEqual(res.stdout.splitlines(), self.cmd[1:])
        self.assertIn("is tmpfs, memory", res.stderr)

    def test_no_exec_subcommand_runs_unchanged(self):
        res = run_script(["TMPDIR", "/mnt/scratch", "--", str(self.app), "run", "/sif/x.sif"], {"TMPDIR": str(self.scratch)})
        self.assertEqual(res.returncode, 0, res.stderr)
        self.assertEqual(res.stdout.splitlines(), ["run", "/sif/x.sif"])
        self.assertIn("no exec subcommand", res.stderr)

    def test_usage_errors(self):
        for args in ([], ["TMPDIR"], ["TMPDIR", "/mnt/scratch", str(self.app), "exec"], ["TMPDIR", "/mnt/scratch", "-x", "a", "b"]):
            res = run_script(args, {})
            self.assertEqual(res.returncode, 2, args)
            self.assertIn("usage:", res.stderr)


class BackendTest(unittest.TestCase):
    def test_with_exec_script(self):
        stock = ["sbatch", "--wait", "--parsable", "--job-name", "call-x", "--ntasks", "1", "--cpus-per-task", "4",
                 "/venv/lib/python3.12/site-packages/miniwdl_slurm/scripts/exec_script.sh"]
        self.assertEqual(scratch.with_exec_script(stock, "", "/mnt/scratch"), stock)
        self.assertIsNot(scratch.with_exec_script(stock, "", "/mnt/scratch"), stock)
        out = scratch.with_exec_script(stock, "TMPDIR", "/mnt/scratch")
        self.assertEqual(out[:-4], stock[:-1])
        self.assertEqual(out[-4:], [scratch.EXEC_SCRIPT, "TMPDIR", "/mnt/scratch", "--"])
        with self.assertRaisesRegex(RuntimeError, "do not end with its exec script"):
            scratch.with_exec_script(stock[:-1], "TMPDIR", "/mnt/scratch")
        with self.assertRaises(RuntimeError):
            scratch.with_exec_script([], "TMPDIR", "/mnt/scratch")

    def test_settings_from_cfg(self):
        class Cfg:
            def __init__(self, values):
                self.values = values

            def get(self, section, key, default=""):
                return self.values.get((section, key), default)

        self.assertEqual(scratch.scratch_settings(Cfg({})), ("", "/mnt/scratch"))
        self.assertEqual(scratch.scratch_settings(Cfg({("ugc_wgw", "scratch_env"): " TMPDIR ",
                                                       ("ugc_wgw", "scratch_mount"): "/scratch"})), ("TMPDIR", "/scratch"))
        self.assertEqual(scratch.scratch_settings(Cfg({("ugc_wgw", "scratch_env"): "SNIC_TMP",
                                                       ("ugc_wgw", "scratch_mount"): ""})), ("SNIC_TMP", "/mnt/scratch"))

    def test_backend_subclass_with_miniwdl_slurm(self):
        if scratch.SlurmSingularity is None:
            self.skipTest("miniwdl-slurm not importable")
        from WDL.runtime.config import Loader  # type: ignore
        tmp = pathlib.Path(tempfile.mkdtemp())
        for env, expect_ours in (("TMPDIR", True), ("", False)):
            cfg_file = tmp / f"miniwdl-{env or 'off'}.cfg"
            cfg_file.write_text("[scheduler]\ncontainer_backend = ugc_slurm_singularity\n[ugc_wgw]\n"
                                f"scratch_env = {env}\nscratch_mount = /mnt/scratch\n[slurm]\nextra_args = \"--qos short\"\n")
            cfg = Loader(logging.getLogger("test"), filenames=[str(cfg_file)])
            host_dir = tmp / f"host-{env or 'off'}"
            host_dir.mkdir()
            backend = scratch.SlurmSingularityScratch(cfg, "call-pbmm2_align_wgs-3", str(host_dir))
            self.assertEqual(backend.cli_name, "slurm_singularity")   # the log file the driver reads keeps its name
            backend.runtime_values = {"cpu": 4, "memory_reservation": 8 * 1024 ** 3}
            args = backend._slurm_invocation()
            self.assertEqual(args[:3], ["sbatch", "--wait", "--parsable"])
            self.assertIn("--qos", args)
            if expect_ours:
                self.assertEqual(args[-4:], [scratch.EXEC_SCRIPT, "TMPDIR", "/mnt/scratch", "--"])
            else:
                self.assertTrue(args[-1].endswith("miniwdl_slurm/scripts/exec_script.sh"), args[-1])

    def test_entry_point_registered(self):
        import importlib.metadata as m
        eps = m.entry_points()
        group = eps.select(group="miniwdl.plugin.container_backend") if hasattr(eps, "select") else eps.get("miniwdl.plugin.container_backend", [])
        hit = [e for e in group if e.name == "ugc_slurm_singularity"]
        if not hit:
            self.skipTest("plugin not installed in this interpreter")
        self.assertEqual(hit[0].value, "ugc_wgw_miniwdl.scratch:SlurmSingularityScratch")
        self.assertIs(hit[0].load(), scratch.SlurmSingularityScratch)


if __name__ == "__main__":
    unittest.main()
