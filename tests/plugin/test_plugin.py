"""Tests for ugc_wgw_miniwdl.resources: the miniwdl task-plugin coroutine driven by hand, no container backend."""
from __future__ import annotations

import logging
import os
import sys
import tempfile
import unittest
from types import SimpleNamespace

from ugc_wgw_miniwdl import policy, resources

G = 1024 ** 3
HEADER = "task\tcpu\tmemory\ttime\tpartition\tconstraint\n"


class FakeCfg:
    def __init__(self, **sections):
        self.sections = sections

    def get(self, section, key, default=None):
        value = self.sections.get(section, {}).get(key)
        if value is None:
            if default is not None:
                return default
            raise KeyError(f"[{section}] {key}")
        return value


class FakeContainer:
    limits = {"cpu": sys.maxsize, "mem_bytes": sys.maxsize}

    def __init__(self, **runtime_values):
        self.runtime_values = dict(runtime_values)

    @classmethod
    def detect_resource_limits(cls, cfg, logger):
        return cls.limits


class Capture(logging.Handler):
    def __init__(self):
        super().__init__(level=1)
        self.records = []

    def emit(self, record):
        self.records.append(record)


def drive(cfg, task_name, container, logger, run_dir="/tmp/run"):
    """Run the three stages the way miniwdl's compose_coroutines does; returns the stage-2 dict."""
    gen = resources.task(cfg, logger, ["wf", "call-x"], run_dir, SimpleNamespace(name=task_name), inputs="IN")
    first = next(gen)
    assert first == {"inputs": "IN"}
    out = gen.send({"command": "cmd", "container": container})
    third = gen.send({"outputs": "OUT"})
    assert third == {"outputs": "OUT"}
    gen.close()
    return out


class PluginTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.path = os.path.join(self.tmp.name, "resources.tsv")
        with open(self.path, "w") as fh:
            fh.write(HEADER + "deepvariant_call_variants_cpu\t48\t-\t-\t-\t-\nugc_wgw_hifiasm_assemble\t-\t-\t2-00\tfat\t-\n")
        self.logger = logging.getLogger(f"test.{id(self)}")
        self.logger.setLevel(1)
        self.capture = Capture()
        self.logger.addHandler(self.capture)
        policy._CACHE.clear()

    def tearDown(self):
        self.logger.removeHandler(self.capture)
        self.tmp.cleanup()

    def messages(self, level=None):
        return [r.msg for r in self.capture.records if level is None or r.levelno == level]

    def test_applies_and_logs_once(self):
        cfg = FakeCfg(ugc_wgw={"resources": self.path})
        c = FakeContainer(cpu=64, memory_reservation=256 * G, time_minutes=4320)
        out = drive(cfg, "deepvariant_call_variants_cpu", c, self.logger)
        self.assertIs(out["container"], c)
        self.assertEqual(out["command"], "cmd")
        self.assertEqual(c.runtime_values["cpu"], 48)
        notices = self.messages(resources.NOTICE_LEVEL)
        self.assertEqual(len(notices), 1)
        msg = notices[0]
        self.assertEqual(msg.message, resources.MESSAGE)
        self.assertEqual(msg.kwargs["task"], "deepvariant_call_variants_cpu")
        self.assertEqual(msg.kwargs["changes"], {"cpu": [64, 48]})
        self.assertEqual(msg.kwargs["summary"], "cpu 64→48")
        self.assertEqual(msg.kwargs["rules"], [f"{self.path}:2"])
        self.assertIn('changes: {"cpu": [64, 48]}', str(msg))

    def test_partition_and_time(self):
        cfg = FakeCfg(ugc_wgw={"resources": self.path})
        c = FakeContainer(cpu=48, memory_reservation=288 * G, time_minutes=4320)
        drive(cfg, "ugc_wgw_hifiasm_assemble", c, self.logger)
        self.assertEqual(c.runtime_values["slurm_partition"], "fat")
        self.assertEqual(c.runtime_values["time_minutes"], 2880)
        self.assertEqual(self.messages(logging.WARNING), [])

    def test_no_match_no_notice(self):
        cfg = FakeCfg(ugc_wgw={"resources": self.path})
        c = FakeContainer(cpu=8)
        drive(cfg, "mosdepth", c, self.logger)
        self.assertEqual(c.runtime_values, {"cpu": 8})
        self.assertEqual(self.messages(resources.NOTICE_LEVEL), [])

    def test_unconfigured_missing_and_empty_policy_pass_through(self):
        for cfg in (FakeCfg(), FakeCfg(ugc_wgw={"resources": ""}), FakeCfg(ugc_wgw={"resources": self.path + ".none"})):
            c = FakeContainer(cpu=64)
            drive(cfg, "deepvariant_call_variants_cpu", c, self.logger)
            self.assertEqual(c.runtime_values, {"cpu": 64})
        with open(self.path, "w") as fh:
            fh.write("# empty\n")
        os.utime(self.path, (1_600_000_000, 1_600_000_000))
        c = FakeContainer(cpu=64)
        drive(FakeCfg(ugc_wgw={"resources": self.path}), "deepvariant_call_variants_cpu", c, self.logger)
        self.assertEqual(c.runtime_values, {"cpu": 64})
        self.assertEqual([m.message for m in self.messages(resources.VERBOSE_LEVEL)],
                         ["ugc-wgw resource policy missing", "ugc-wgw resource policy empty"])
        self.assertEqual(self.messages(resources.NOTICE_LEVEL), [])

    def test_backend_limits_bound_the_policy(self):
        cfg = FakeCfg(ugc_wgw={"resources": self.path})

        class Local(FakeContainer):
            limits = {"cpu": 24, "mem_bytes": 64 * G}

        c = Local(cpu=24)
        drive(cfg, "deepvariant_call_variants_cpu", c, self.logger)
        self.assertEqual(c.runtime_values["cpu"], 24)
        self.assertEqual(self.messages(resources.NOTICE_LEVEL), [])

    def test_malformed_policy_raises_with_line(self):
        with open(self.path, "w") as fh:
            fh.write(HEADER + "x\tlots\t-\t-\t-\t-\n")
        gen = resources.task(FakeCfg(ugc_wgw={"resources": self.path}), self.logger, [], "/tmp", SimpleNamespace(name="x"),
                             inputs=None)
        next(gen)
        with self.assertRaises(policy.PolicyError) as cm:
            gen.send({"command": "c", "container": FakeContainer(cpu=1)})
        self.assertIn(f"{self.path}:2: column cpu", str(cm.exception))

    def test_warns_when_extra_args_carry_a_partition(self):
        cfg = FakeCfg(ugc_wgw={"resources": self.path}, slurm={"extra_args": "--partition core --account a1"})
        c = FakeContainer(cpu=48)
        drive(cfg, "ugc_wgw_hifiasm_assemble", c, self.logger)
        warnings = self.messages(logging.WARNING)
        self.assertEqual(len(warnings), 1)
        self.assertIn("partition overridden by [slurm] extra_args", warnings[0].message)
        self.assertIn("SLURM_PARTITION", warnings[0].kwargs["hint"])
        for harmless in ("--qos short", "-p", "--partitioned x"):
            self.capture.records.clear()
            cfg = FakeCfg(ugc_wgw={"resources": self.path}, slurm={"extra_args": harmless})
            drive(cfg, "ugc_wgw_hifiasm_assemble", FakeContainer(cpu=48), self.logger)
            expect = 1 if harmless == "-p" else 0
            self.assertEqual(len(self.messages(logging.WARNING)), expect, harmless)

    def test_empty_gpu_type_is_dropped_and_the_request_logged(self):
        # upstream's HPC backend sets gpuType "" when the workflow input is unset: gres must be gpu:N, not gpu::N
        c = FakeContainer(cpu=8, gpuCount=1, gpuType="")
        drive(FakeCfg(), "deepvariant_call_variants_gpu", c, self.logger)  # no policy configured at all
        self.assertEqual(c.runtime_values, {"cpu": 8, "gpuCount": 1})
        notices = [m for m in self.messages(resources.NOTICE_LEVEL) if getattr(m, "message", "") == "ugc-wgw gpu request"]
        self.assertEqual(len(notices), 1)
        self.assertEqual(notices[0].kwargs["gres"], "gpu:1")
        c = FakeContainer(cpu=48, gpuCount=4, gpuType="a100")
        drive(FakeCfg(ugc_wgw={"resources": self.path}), "run_parabricks_deepvariant", c, self.logger)
        self.assertEqual(c.runtime_values["gpuType"], "a100")
        self.assertEqual(self.messages(resources.NOTICE_LEVEL)[-1].kwargs["gres"], "gpu:a100:4")
        before = len(self.messages())
        c = FakeContainer(cpu=8)
        drive(FakeCfg(ugc_wgw={"resources": self.path}), "mosdepth", c, self.logger)
        self.assertEqual(len(self.messages()), before)  # a task without gpuCount logs no GPU request
        self.assertIsNone(resources.gpu_request({"cpu": 8}))
        self.assertEqual(resources.gpu_request({"gpuCount": 2, "gpuType": " "}), "gpu:2")

    def test_partition_row_reaches_the_gpu_partition_key(self):
        with open(self.path, "w") as fh:
            fh.write(HEADER + "run_parabricks_deepvariant\t-\t-\t-\tgpu\t-\nugc_wgw_hifiasm_assemble\t-\t-\t-\tfat\t-\n")
        cfg = FakeCfg(ugc_wgw={"resources": self.path})
        c = FakeContainer(cpu=48, gpuCount=4, gpuType="", slurm_partition="core", slurm_partition_gpu="gpu-default")
        drive(cfg, "run_parabricks_deepvariant", c, self.logger)
        self.assertEqual(c.runtime_values["slurm_partition"], "gpu")
        self.assertEqual(c.runtime_values["slurm_partition_gpu"], "gpu")  # the row wins over the site GPU default
        self.assertNotIn("gpuType", c.runtime_values)
        self.assertEqual(self.messages(resources.NOTICE_LEVEL)[-1].kwargs["partition"], "gpu")
        c = FakeContainer(cpu=48, slurm_partition="core")
        drive(cfg, "ugc_wgw_hifiasm_assemble", c, self.logger)
        self.assertEqual(c.runtime_values["slurm_partition"], "fat")
        self.assertNotIn("slurm_partition_gpu", c.runtime_values)  # not a GPU task

    def test_container_missing_is_tolerated(self):
        gen = resources.task(FakeCfg(ugc_wgw={"resources": self.path}), self.logger, [], "/tmp", SimpleNamespace(name="x"),
                             inputs=None)
        next(gen)
        out = gen.send({"command": "c"})
        self.assertEqual(out, {"command": "c"})
        gen.close()


class RealLoaderTests(unittest.TestCase):
    def setUp(self):
        try:
            from WDL.runtime.config import Loader  # noqa: F401
        except ImportError:
            raise unittest.SkipTest("miniwdl not importable")

    def test_loader_custom_section_and_default(self):
        from WDL.runtime.config import Loader
        cfg = Loader(logging.getLogger("t"), filenames=[], overrides={"ugc_wgw": {"resources": "/x/y.tsv"}})
        self.assertEqual(resources.policy_path(cfg), "/x/y.tsv")
        self.assertEqual(resources.policy_path(Loader(logging.getLogger("t"), filenames=[])), "")
        self.assertEqual(resources.cfg_get(cfg, "slurm", "extra_args", ""), "")


if __name__ == "__main__":
    unittest.main()
