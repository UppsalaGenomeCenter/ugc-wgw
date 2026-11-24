"""Escalation inside a run (ugc_wgw_miniwdl.inrun): a try killed for wall time or memory is retried with a
larger request, and never submitted again with the same one. The pure parts by hand; then miniwdl itself, run
through the ugc_slurm_singularity backend against stand-ins for sbatch, scontrol and singularity."""
from __future__ import annotations

import json
import logging
import os
import pathlib
import stat
import subprocess
import sys
import tempfile
import unittest
from types import SimpleNamespace

from ugc_wgw_miniwdl import inrun, resources, scratch

G = 1024 ** 3
TIME_LINE = "slurmstepd: error: *** JOB 77 ON n1 CANCELLED AT 2026-11-24T10:00:00 DUE TO TIME LIMIT ***"
OOM_LINE = "slurmstepd: error: Detected 1 oom_kill event in StepId=77.batch. Some of the step tasks have been OOM Killed."


class Cfg:
    def __init__(self, **values):
        self.values = values

    def get(self, section, key, default=None):
        return self.values.get(f"{section}.{key}", default)


class Container:
    """What inrun needs of a miniwdl TaskContainer."""
    cli_name = "slurm_singularity"
    limits = {"cpu": sys.maxsize, "mem_bytes": sys.maxsize}

    def __init__(self, host_dir, cfg=None, **runtime_values):
        self.host_dir = str(host_dir)
        self.cfg = cfg or Cfg()
        self.runtime_values = dict(runtime_values)
        self.try_counter = 1
        self.run_id = "call-t"

    @classmethod
    def detect_resource_limits(cls, cfg, logger):
        return cls.limits

    def _stream(self, stem):
        return os.path.join(self.host_dir, f"{stem}{self.try_counter if self.try_counter > 1 else ''}.txt")

    def host_stderr_txt(self):
        return self._stream("stderr")

    def host_stdout_txt(self):
        return self._stream("stdout")

    def success_exit_code(self, code):
        return code == 0


class Capture(logging.Handler):
    def __init__(self):
        super().__init__(level=1)
        self.records = []

    def emit(self, record):
        self.records.append(record)


def no_scontrol(*args, **kwargs):
    raise FileNotFoundError("scontrol")


def scontrol_says(state):
    def runner(cmd, **kwargs):
        return SimpleNamespace(returncode=0, stdout=f"JobId={cmd[4]} JobName=x JobState={state} Reason=None\n")
    return runner


class Base(unittest.TestCase):
    def setUp(self):
        self.tmp = pathlib.Path(tempfile.mkdtemp())
        self.capture = Capture()
        self.logger = logging.getLogger(f"inrun.{self.id()}")
        self.logger.setLevel(1)
        self.logger.propagate = False
        self.logger.addHandler(self.capture)
        self.slept = []

    def messages(self, text):
        return [r.msg for r in self.capture.records if getattr(r.msg, "message", str(r.msg)) == text]

    def fields(self, msg):
        return getattr(msg, "kwargs", None) or getattr(msg, "_kwargs", None) or {}

    def kill(self, exit_code=1, job="77", runner=no_scontrol, streams=(), since=0.0):
        return inrun.kill_reason(str(self.tmp), exit_code, job, "", list(streams), since, runner, self.slept.append)


class KillReasonTest(Base):
    def test_job_id_from_the_sbatch_log(self):
        log = self.tmp / "slurm_singularity.log.txt"
        for text, want in (("4711\n", ("4711", "")), ("4711;rackham\nsbatch: note\n", ("4711", "rackham")),
                           ("sbatch: error: Requested time limit is invalid\n", ("", "")), ("", ("", ""))):
            log.write_text(text)
            self.assertEqual(inrun.job_of(str(log)), want, text)
        self.assertEqual(inrun.job_of(str(self.tmp / "missing")), ("", ""))

    def test_time_and_memory_from_the_jobs_own_output(self):
        (self.tmp / "slurm-77.out").write_text("[ugc-wgw] scratch: /scratch/77 -> /mnt/scratch (TMPDIR)\n" + TIME_LINE + "\n")
        self.assertEqual(self.kill(), ("time", TIME_LINE))
        self.assertEqual(self.kill(exit_code=137), ("time", TIME_LINE))      # a job out of time is killed too
        (self.tmp / "slurm-77.out").write_text(OOM_LINE + "\n")
        self.assertEqual(self.kill(), ("memory", OOM_LINE))
        # an earlier try's output says nothing about this one
        (self.tmp / "slurm-76.out").write_text(TIME_LINE + "\n")
        (self.tmp / "slurm-77.out").write_text("")
        self.assertIsNone(self.kill(runner=scontrol_says("FAILED")))
        self.assertEqual(self.slept, [])

    def test_state_exit_status_and_pinned_memory(self):
        self.assertEqual(self.kill(runner=scontrol_says("TIMEOUT")), ("time", "job 77 state TIMEOUT"))
        self.assertEqual(self.kill(runner=scontrol_says("OUT_OF_MEMORY")), ("memory", "job 77 state OUT_OF_MEMORY"))
        self.assertEqual(self.kill(exit_code=253, runner=scontrol_says("FAILED")), ("memory", "exit status 253"))
        self.assertEqual(self.kill(exit_code=137), ("memory", "exit status 137"))
        out = self.tmp / "stdout.txt"
        out.write_text("[PB Error] cudaErrorMemoryAllocation: out of memory\n")
        self.assertEqual(self.kill(exit_code=255, runner=scontrol_says("FAILED"), streams=[str(self.tmp / "stderr.txt"), str(out)])[0], "memory")
        # a tool's own failure is not a kill, whatever its stderr says about time
        err = self.tmp / "stderr.txt"
        err.write_text("TIMEOUT waiting for the license server\n")
        self.assertIsNone(self.kill(exit_code=1, runner=scontrol_says("FAILED"), streams=[str(err)]))
        self.assertEqual(self.slept, [])

    def test_a_silent_controller_gives_the_output_file_time(self):
        self.assertIsNone(self.kill())
        self.assertEqual(self.slept, [inrun.RECHECK_SECONDS] * inrun.RECHECKS)
        self.slept.clear()

        def arrive(seconds):
            self.slept.append(seconds)
            (self.tmp / "slurm-77.out").write_text(TIME_LINE + "\n")
        self.assertEqual(inrun.kill_reason(str(self.tmp), 1, "77", "", [], 0.0, no_scontrol, arrive), ("time", TIME_LINE))
        self.assertEqual(len(self.slept), 1)

    def test_without_a_job_id_only_output_written_since_the_try_began(self):
        old = self.tmp / "slurm-70.out"
        old.write_text(TIME_LINE + "\n")
        os.utime(old, (1000, 1000))
        self.assertIsNone(self.kill(job="", since=2000.0))
        new = self.tmp / "slurm-71.out"
        new.write_text(OOM_LINE + "\n")
        os.utime(new, (2005, 2005))
        self.assertEqual(self.kill(job="", since=2000.0), ("memory", OOM_LINE))


class EscalateTest(Base):
    def container(self, cfg=None, **values):
        c = Container(self.tmp, cfg, **values)
        (self.tmp / "slurm_singularity.log.txt").write_text("77\n")
        return c

    def fail(self, c, line=TIME_LINE, exit_code=1):
        (self.tmp / "slurm-77.out").write_text(line + "\n")
        return inrun.after_failure(c, self.logger, exit_code, 0.0, no_scontrol, self.slept.append)

    def test_time_doubles_twice_then_stops(self):
        c = self.container(time_minutes=120, memory_reservation=64 * G)
        inrun.prepare(c.cfg, c, "ugc_wgw_cactus_chrom")
        self.assertEqual(self.fail(c), "time")
        self.assertEqual(self.fail(c), "time")
        self.assertEqual((c.runtime_values["time_minutes"], c.runtime_values["memory_reservation"]), (480, 64 * G))
        notices = self.messages(resources.ESCALATION_MESSAGE)
        self.assertEqual([self.fields(m)["changes"] for m in notices], [{"time_minutes": [120, 240]}, {"time_minutes": [240, 480]}])
        self.assertEqual([(self.fields(m)["n"], self.fields(m)["run_factor"], self.fields(m)["in_run"]) for m in notices],
                         [(1, 2.0, True), (2, 4.0, True)])
        self.assertEqual(self.fields(notices[0])["task"], "ugc_wgw_cactus_chrom")
        self.assertIsNone(inrun.skip(c, self.logger))
        # the third kill: nothing left, and the tries miniwdl still has are answered without a job
        self.assertEqual(self.fail(c), "time")
        self.assertEqual(c.runtime_values["time_minutes"], 480)
        stop = self.fields(self.messages(inrun.KILL_STOP_MESSAGE)[0])
        self.assertEqual((stop["code"], stop["reason"], stop["resource"], stop["n"]), ("exhausted", "escalation exhausted (2 of 2)", "time", 2))
        c.try_counter = 4
        self.assertEqual(inrun.skip(c, self.logger), 1)
        self.assertIn("not submitted again", pathlib.Path(c.host_stderr_txt()).read_text())
        self.assertTrue(os.path.exists(c.host_stdout_txt()))
        self.assertEqual(len(self.messages(inrun.SKIPPED_MESSAGE)), 1)

    def test_the_drivers_factor_counts_and_memory_has_its_own_ladder(self):
        cfg = Cfg(**{"ugc_wgw.escalations": json.dumps({"ugc_wgw_pav": {"time_factor": 2.0, "memory_min": 80 * G}})})
        c = self.container(cfg, time_minutes=240, memory_reservation=80 * G, memory_limit=80 * G)
        state = inrun.prepare(cfg, c, "ugc_wgw_pav")
        self.assertEqual(state.launched, {"time": 1})
        self.fail(c)
        self.assertEqual(c.runtime_values["time_minutes"], 480)
        self.fail(c)                                   # one launched with, one made here: exhausted
        self.assertEqual(c.runtime_values["time_minutes"], 480)
        self.assertIsNotNone(state.stop)
        state.stop = None
        self.fail(c, OOM_LINE, 137)
        self.fail(c, OOM_LINE, 137)
        self.assertEqual((c.runtime_values["memory_reservation"], c.runtime_values["memory_limit"]), (320 * G, 320 * G))

    def test_settings_from_the_driver(self):
        cfg = Cfg(**{"ugc_wgw.resource_retry_factor": "1.5", "ugc_wgw.resource_retry_max": "1"})
        self.assertEqual(inrun.retry_settings(cfg), (1.5, 1))
        self.assertEqual(inrun.retry_settings(Cfg()), (2.0, 2))
        self.assertEqual(inrun.retry_settings(Cfg(**{"ugc_wgw.resource_retry_factor": "x", "ugc_wgw.resource_retry_max": ""})), (2.0, 2))
        c = self.container(cfg, time_minutes=100)
        inrun.prepare(cfg, c, "t")
        self.fail(c)
        self.assertEqual(c.runtime_values["time_minutes"], 150)
        self.fail(c)
        self.assertEqual(c.runtime_values["time_minutes"], 150)
        # resource_retry_max 0: no escalation, and still no second job with the request that was killed
        zero = Cfg(**{"ugc_wgw.resource_retry_max": "0"})
        c = self.container(zero, time_minutes=100)
        inrun.prepare(zero, c, "t")
        self.fail(c)
        self.assertEqual((c.runtime_values["time_minutes"], inrun.state_of(c).stop), (100, (1, "resource_retry_max is 0")))

    def test_memory_stops_at_the_site_cap(self):
        cfg = Cfg(**{"task_runtime.memory_max": str(375 * G)})
        c = self.container(cfg, memory_reservation=200 * G)
        inrun.prepare(cfg, c, "ugc_wgw_pav")
        self.fail(c, OOM_LINE, 137)
        self.assertEqual(c.runtime_values["memory_reservation"], 375 * G)
        self.fail(c, OOM_LINE, 137)
        self.assertEqual(c.runtime_values["memory_reservation"], 375 * G)
        stop = self.fields(self.messages(inrun.KILL_STOP_MESSAGE)[0])
        self.assertEqual((stop["code"], stop["reason"]), ("capped", "memory is at the site's cap (375G)"))
        # a policy row put the task above the cap, on a partition of its own: that is not capped
        c = self.container(cfg, memory_reservation=500 * G)
        inrun.prepare(cfg, c, "ugc_wgw_cactus_chrom")
        self.fail(c, OOM_LINE, 137)
        self.assertEqual(c.runtime_values["memory_reservation"], 1000 * G)

    def test_no_request_to_raise_and_other_failures(self):
        c = self.container(memory_reservation=8 * G)       # the site sets no wall time
        inrun.prepare(c.cfg, c, "t")
        self.fail(c)
        self.assertEqual(self.fields(self.messages(inrun.KILL_STOP_MESSAGE)[0])["code"], "no_request")
        c = self.container(time_minutes=60)
        (self.tmp / "slurm-77.out").write_text("")
        self.assertIsNone(inrun.after_failure(c, self.logger, 3, 0.0, scontrol_says("FAILED"), self.slept.append))
        self.assertEqual(c.runtime_values, {"time_minutes": 60})
        self.assertIsNone(inrun.skip(c, self.logger))

    def test_a_request_sbatch_refused_is_not_offered_again(self):
        c = self.container(time_minutes=120)
        inrun.prepare(c.cfg, c, "ugc_wgw_cactus_chrom")
        self.fail(c)
        self.assertEqual(c.runtime_values["time_minutes"], 240)
        refusal = "sbatch: error: Batch job submission failed: Requested time limit is invalid (missing or exceeds some limit)"
        (self.tmp / "slurm_singularity.log.txt").write_text(refusal + "\n")
        self.assertIsNone(inrun.after_failure(c, self.logger, 1, 0.0, no_scontrol, self.slept.append))
        stop = self.fields(self.messages(inrun.KILL_STOP_MESSAGE)[0])
        self.assertEqual(stop["code"], "refused")
        self.assertIn("Requested time limit is invalid", stop["reason"])
        self.assertEqual(inrun.skip(c, self.logger), 1)
        self.assertEqual(self.slept, [])
        # a controller that cannot be reached is not a refusal: miniwdl retries
        c = self.container(time_minutes=120)
        (self.tmp / "slurm_singularity.log.txt").write_text("sbatch: error: Batch job submission failed: Socket timed out on send/recv operation\n")
        (self.tmp / "slurm-77.out").unlink()
        self.assertIsNone(inrun.after_failure(c, self.logger, 1, 1e12, no_scontrol, self.slept.append))
        self.assertIsNone(inrun.skip(c, self.logger))

    def test_the_task_plugin_attaches_the_state_and_run_try_uses_it(self):
        cfg = Cfg()
        c = self.container(cfg, time_minutes=120)
        gen = resources.task(cfg, self.logger, ["wf", "call-t"], "/tmp/run", SimpleNamespace(name="ugc_wgw_pav"), inputs="IN")
        next(gen)
        gen.send({"command": "cmd", "container": c})
        self.assertEqual(inrun.state_of(c).task, "ugc_wgw_pav")
        submitted = []

        def submit():
            submitted.append(c.runtime_values["time_minutes"])
            (self.tmp / "slurm-77.out").write_text(TIME_LINE + "\n")
            return 1
        real_run, real_sleep = inrun.subprocess.run, inrun.time.sleep
        inrun.subprocess.run, inrun.time.sleep = no_scontrol, self.slept.append
        try:
            for attempt in range(1, 5):
                c.try_counter = attempt
                self.assertEqual(scratch.run_try(c, self.logger, lambda: False, submit), 1)
            # a run that is being cancelled is left alone, and so is a try that succeeded
            c2 = self.container(cfg, time_minutes=120)
            self.assertEqual(scratch.run_try(c2, self.logger, lambda: True, submit), 1)
            self.assertEqual(scratch.run_try(c2, self.logger, lambda: False, lambda: 0), 0)
            self.assertEqual(c2.runtime_values["time_minutes"], 120)
        finally:
            inrun.subprocess.run, inrun.time.sleep = real_run, real_sleep
        self.assertEqual(submitted[:3], [120, 240, 480])   # the fourth try had no job
        self.assertEqual(len(submitted), 4)                # (the fourth entry is c2's cancelled try)


SBATCH = r'''#!/usr/bin/env python3
"""A stand-in for `sbatch --wait --parsable`: a job that needs FAKE_NEED_MINUTES of wall time and FAKE_NEED_MEM_MB
of memory. Asked for less, it ends the way SLURM ends such a job; asked for enough, it runs the task's command
where the container's binds point."""
import os, subprocess, sys
fake, args = os.environ["FAKE_SLURM"], sys.argv[1:]


def opt(name):
    return args[args.index(name) + 1] if name in args else ""


counter = os.path.join(fake, "next")
job = int(open(counter).read()) if os.path.exists(counter) else 100
open(counter, "w").write(str(job + 1))
minutes, mem = int(opt("--time") or 0), int((opt("--mem") or "0M").rstrip("M"))
with open(os.path.join(fake, "calls"), "a") as fh:
    fh.write(f"{job}\t{minutes}\t{mem}\n")
print(job, flush=True)


def end(state, line, code):
    open(os.path.join(fake, f"state.{job}"), "w").write(state)
    if line:
        open(f"slurm-{job}.out", "w").write(line + "\n")
    sys.exit(code)


if minutes < int(os.environ.get("FAKE_NEED_MINUTES", "0")):
    end("TIMEOUT", f"slurmstepd: error: *** JOB {job} ON n1 CANCELLED AT 2026-11-24T10:00:00 DUE TO TIME LIMIT ***", 1)
if mem < int(os.environ.get("FAKE_NEED_MEM_MB", "0")):
    end("OUT_OF_MEMORY", "" if os.environ.get("FAKE_QUIET_OOM") else
        f"slurmstepd: error: Detected 1 oom_kill event in StepId={job}.batch. Some of the step tasks have been OOM Killed.", 253)
binds = {}
for i, a in enumerate(args):
    if a == "--bind":
        host, inside = args[i + 1].split(":")[:2]
        binds[inside] = host
base = "/mnt/miniwdl_task_container"
cmd = f"bash {binds[base + '/command']} >> {binds[base + '/stdout.txt']} 2>> {binds[base + '/stderr.txt']}"
rc = subprocess.run(["/bin/sh", "-c", cmd], cwd=binds[base + "/work"]).returncode
end("COMPLETED" if rc == 0 else "FAILED", "", rc)
'''

SCONTROL = r'''#!/usr/bin/env python3
import os, sys
path = os.path.join(os.environ["FAKE_SLURM"], f"state.{sys.argv[-1]}")
if os.environ.get("FAKE_SCONTROL") == "off" or not os.path.exists(path):
    sys.exit(1)
print(f"JobId={sys.argv[-1]} JobName=x JobState={open(path).read()} Reason=None")
'''

SINGULARITY = r'''#!/usr/bin/env python3
import sys
if sys.argv[1:2] == ["--version"]:
    print("stand-in 1.0")
elif sys.argv[1:2] == ["pull"]:
    open(sys.argv[2], "w").close()
'''

WDL = '''version 1.0
workflow w {
  input { Int code = 0 }
  call t { input: code = code }
  output { String said = t.said }
}
task t {
  input { Int code }
  command <<<
    echo hello
    exit ~{code}
  >>>
  output { String said = read_string(stdout()) }
  runtime {
    docker: "ubuntu:20.04"
    cpu: 1
    memory: "1 GiB"
    maxRetries: 3
  }
}
'''


class MiniwdlTest(unittest.TestCase):
    """miniwdl's own retry loop around our backend: what each try asked SLURM for."""

    @classmethod
    def setUpClass(cls):
        try:
            import WDL  # noqa: F401
            import miniwdl_slurm  # noqa: F401
            import importlib.metadata as m
        except ImportError:
            raise unittest.SkipTest("miniwdl and miniwdl-slurm not importable")
        eps = m.entry_points()
        group = eps.select(group="miniwdl.plugin.container_backend") if hasattr(eps, "select") else eps.get("miniwdl.plugin.container_backend", [])
        if not any(e.name == scratch.BACKEND_NAME for e in group):
            raise unittest.SkipTest("plugin not installed in this interpreter")

    def setUp(self):
        self.tmp = pathlib.Path(tempfile.mkdtemp())
        self.bin = self.tmp / "bin"
        self.fake = self.tmp / "slurm"
        for d in (self.bin, self.fake, self.tmp / "sif"):
            d.mkdir()
        for name, text in (("sbatch", SBATCH), ("scontrol", SCONTROL), ("singularity", SINGULARITY)):
            exe = self.bin / name
            exe.write_text(text.replace("#!/usr/bin/env python3", f"#!{sys.executable}", 1))
            exe.chmod(exe.stat().st_mode | stat.S_IXUSR)
        (self.tmp / "w.wdl").write_text(WDL)
        defaults = json.dumps({"docker": "ubuntu:20.04", "maxRetries": 1, "time_minutes": 120})
        (self.tmp / "miniwdl.cfg").write_text(
            "[scheduler]\ncontainer_backend = ugc_slurm_singularity\n"
            f"[singularity]\nexe = [\"singularity\"]\nimage_cache = {self.tmp / 'sif'}\nrun_options = [\"--containall\"]\n"
            f"[task_runtime]\ndefaults = {defaults}\n[ugc_wgw]\nresources =\nscratch_env =\n")

    def run_wdl(self, code=0, **env):
        run_dir = self.tmp / "run"
        cmd = [sys.executable, "-m", "WDL", "run", str(self.tmp / "w.wdl"), f"code={code}", "--dir", f"{run_dir}/.",
               "--cfg", str(self.tmp / "miniwdl.cfg"), "--no-color", "--log-json", "-o", str(run_dir / "run.json")]
        run_dir.mkdir()
        environ = {k: v for k, v in os.environ.items() if not k.startswith("MINIWDL")}
        environ.update({"PATH": f"{self.bin}{os.pathsep}{os.environ['PATH']}", "FAKE_SLURM": str(self.fake), **env})
        res = subprocess.run(cmd, env=environ, capture_output=True, text=True, cwd=str(self.tmp))
        calls = [tuple(int(x) for x in line.split("\t")[1:]) for line in (self.fake / "calls").read_text().splitlines()]
        log = [json.loads(line) for line in (run_dir / "workflow.log.json").read_text().splitlines() if line.startswith("{")] \
            if (run_dir / "workflow.log.json").exists() else []
        return res, calls, log, run_dir

    def named(self, log, message):
        return [d for d in log if d.get("message") == message]

    def test_a_timeout_is_retried_with_twice_the_time(self):
        res, calls, log, run_dir = self.run_wdl(FAKE_NEED_MINUTES="300")
        self.assertEqual(res.returncode, 0, res.stderr[-3000:])
        self.assertEqual(calls, [(120, 1024), (240, 1024), (480, 1024)])
        self.assertEqual(json.loads((run_dir / "run.json").read_text())["outputs"]["w.said"], "hello")
        notices = self.named(log, resources.ESCALATION_MESSAGE)
        self.assertEqual([(d["changes"], d["run_factor"], d["task"], d["in_run"]) for d in notices],
                         [({"time_minutes": [120, 240]}, 2.0, "t", True), ({"time_minutes": [240, 480]}, 4.0, "t", True)])
        self.assertIn("DUE TO TIME LIMIT", notices[0]["evidence"])
        self.assertEqual(len(self.named(log, "failed task will be retried")), 2)

    def test_exhausted_escalation_submits_no_fourth_job(self):
        res, calls, log, run_dir = self.run_wdl(FAKE_NEED_MINUTES="10000")
        self.assertNotEqual(res.returncode, 0)
        self.assertEqual(calls, [(120, 1024), (240, 1024), (480, 1024)])
        stop = self.named(log, inrun.KILL_STOP_MESSAGE)
        self.assertEqual([(d["code"], d["resource"], d["n"], d["max"]) for d in stop], [("exhausted", "time", 2, 2)])
        self.assertEqual(len(self.named(log, inrun.SKIPPED_MESSAGE)), 1)
        error = json.loads((run_dir / "run.json").read_text())
        self.assertEqual((error["cause"]["error"], error["cause"]["exit_status"]), ("CommandFailed", 1))
        # the driver finds the kill where it always has: the newest job output in the failed task's directory
        task_dir = pathlib.Path(error["cause"]["dir"])
        self.assertEqual(len(list(task_dir.glob("slurm-*.out"))), 3)
        self.assertIn("not submitted again", pathlib.Path(error["cause"]["stderr_file"]).read_text())

    def test_memory_by_exit_status_alone(self):
        res, calls, log, _ = self.run_wdl(FAKE_NEED_MEM_MB="3000", FAKE_QUIET_OOM="1", FAKE_SCONTROL="off")
        self.assertEqual(res.returncode, 0, res.stderr[-3000:])
        self.assertEqual(calls, [(120, 1024), (120, 2048), (120, 4096)])
        self.assertEqual(self.named(log, resources.ESCALATION_MESSAGE)[0]["evidence"], "exit status 253")

    def test_the_drivers_settings_and_factors_reach_the_backend(self):
        res, calls, log, _ = self.run_wdl(FAKE_NEED_MINUTES="300", MINIWDL__UGC_WGW__RESOURCE_RETRY_MAX="0",
                                          MINIWDL__UGC_WGW__RESOURCE_RETRY_FACTOR="2")
        self.assertNotEqual(res.returncode, 0)
        self.assertEqual(calls, [(120, 1024)])
        self.assertEqual(len(self.named(log, inrun.SKIPPED_MESSAGE)), 3)
        self.assertNotIn("unused", " ".join(d.get("message", "") for d in log if "RESOURCE_RETRY" in json.dumps(d)))

    def test_launched_with_a_factor_one_escalation_is_left(self):
        (self.fake / "calls").unlink(missing_ok=True)
        res, calls, log, _ = self.run_wdl(FAKE_NEED_MINUTES="10000",
                                          MINIWDL__UGC_WGW__ESCALATIONS=json.dumps({"t": {"time_factor": 2.0}}))
        self.assertNotEqual(res.returncode, 0)
        self.assertEqual(calls, [(240, 1024), (480, 1024)])
        self.assertEqual(self.named(log, inrun.KILL_STOP_MESSAGE)[0]["n"], 2)

    def test_a_failing_command_is_retried_as_before(self):
        res, calls, log, _ = self.run_wdl(code=3)
        self.assertNotEqual(res.returncode, 0)
        self.assertEqual(calls, [(120, 1024)] * 4)
        self.assertEqual(self.named(log, resources.ESCALATION_MESSAGE) + self.named(log, inrun.KILL_STOP_MESSAGE), [])


if __name__ == "__main__":
    unittest.main()
