"""Escalation inside a run: a task killed for wall time or memory is retried with a larger request, never the same.

miniwdl retries a failed command `runtime.maxRetries` times with the request it failed with, and upstream's HPC
backend sets that to 3: a task that hits its wall time ran four times at the same limit before the driver saw the
failure and doubled it (the 864-sample cohort, 2026-11-24: eight hours for one chromosome graph). The number of
retries is read by miniwdl before any plugin sees the container, so it cannot be lowered from here. What the
container backend can do is change what a retry asks for: the request is rebuilt from `runtime_values` at every try.

So after each failed try (scratch.SlurmSingularityScratch._run) this module asks why the job ended:

- killed for time or memory, and an escalation is left: the request is multiplied by `[ugc_wgw]
  resource_retry_factor` (2) before miniwdl's own retry, up to `[ugc_wgw] resource_retry_max` (2) escalations
  per task and resource, counting those the driver launched the run with. The change is logged with the message
  the driver already reads (resources.ESCALATION_MESSAGE), so the request the task finishes with is learned.
- killed, and no escalation is left (the maximum, the site's memory cap, no request to raise): the tries miniwdl
  still has are answered with the same exit status at once, without a job. The run then fails as it would have,
  and the driver reads from the log why (KILL_STOP_MESSAGE).
- not submitted at all because sbatch refused the request (a wall time above the partition's limit, a node
  shape that does not exist): the same request is not offered again either.
- any other failure: nothing changes; miniwdl retries as before.

The evidence of a kill is this try's own: the job's `slurm-<jobid>.out` in the task directory (slurmstepd's
"CANCELLED ... DUE TO TIME LIMIT", "oom_kill event"), the job's state from `scontrol`, sbatch's exit status (253
out of memory, 137 SIGKILL), and for the one tool that fails cleanly when pinned host memory does not fit
(Parabricks), its message in the try's stdout or stderr. The same patterns as bin/ugc_wgw/failures.py.
"""
from __future__ import annotations

import glob
import json
import logging
import math
import os
import re
import subprocess
import time
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional, Tuple

from . import policy, resources

KILL_STOP_MESSAGE = "ugc-wgw resource kill not resubmitted"
SKIPPED_MESSAGE = "ugc-wgw try skipped after a resource kill"
RUNTIME_KEY = {"time": "time_minutes", "memory": "memory_reservation"}
MEMORY_EXIT = frozenset({137, 253})   # 137 = SIGKILL (cgroup OOM), 253 = sbatch --wait out-of-memory
TAIL_BYTES = 64 * 1024

TIME_TEXT = re.compile(r"DUE TO TIME LIMIT")
OOM_TEXT = re.compile(r"oom_kill|oom-kill|OUT_OF_MEMORY|Out of memory|DUE TO MEMORY")
PINNED_TEXT = re.compile(r"cudaErrorMemoryAllocation|CUDA_ERROR_OUT_OF_MEMORY|cudaErrorHostMemoryAlreadyRegistered")
REFUSED_TEXT = re.compile(r"Requested node configuration is not available|CPU count per node can not be satisfied"
                          r"|Memory specification can not be satisfied|Requested time limit is invalid"
                          r"|exceeds partition limit")
JOB_STATE = re.compile(r"\bJobState=(\w+)")
STATE_RESOURCE = {"TIMEOUT": "time", "OUT_OF_MEMORY": "memory"}
RECHECKS = 2            # the job's output file may trail the job on a shared file system: look again, twice
RECHECK_SECONDS = 5.0


@dataclass
class State:
    """What one task's container carries between its tries."""
    task: str
    factor: float = 2.0
    max: int = 2
    launched: Dict[str, int] = field(default_factory=dict)   # resource -> escalations the run started with
    steps: Dict[str, int] = field(default_factory=dict)      # resource -> escalations made here
    stop: Optional[Tuple[int, str]] = None                   # (exit status, why): tries are no longer submitted
    skipped: int = 0


def retry_settings(cfg: Any) -> Tuple[float, int]:
    """([ugc_wgw] resource_retry_factor, resource_retry_max): the driver passes its own per run; 2 and 2 without."""
    factor, most = 2.0, 2
    try:
        factor = float(resources.cfg_get(cfg, resources.CFG_SECTION, "resource_retry_factor", "2") or 2)
    except ValueError:
        pass
    try:
        most = int(float(resources.cfg_get(cfg, resources.CFG_SECTION, "resource_retry_max", "2") or 2))
    except ValueError:
        pass
    return (factor if factor > 1 else 2.0), max(0, most)


def launched_steps(cfg: Any, task_name: str, factor: float) -> Dict[str, int]:
    """How many escalations the driver's factors for this task amount to (resource -> n)."""
    raw = resources.cfg_get(cfg, resources.CFG_SECTION, resources.ESCALATIONS_KEY, "").strip()
    if not raw:
        return {}
    try:
        doc = json.loads(raw)
    except ValueError:
        return {}
    spec = doc.get(task_name) if isinstance(doc, dict) else None
    out: Dict[str, int] = {}
    if isinstance(spec, dict):
        for resource in RUNTIME_KEY:
            given = spec.get(f"{resource}_factor")
            if isinstance(given, (int, float)) and given > 1:
                out[resource] = int(round(math.log(float(given), factor)))
    return out


def prepare(cfg: Any, container: Any, task_name: str) -> State:
    """Attach the state to a container (the task plugin calls this for every task it launches, so the driver's
    per-run settings count as read)."""
    factor, most = retry_settings(cfg)
    state = State(task_name, factor, most, launched_steps(cfg, task_name, factor))
    try:
        container._ugc_wgw_inrun = state
    except AttributeError:  # pragma: no cover - a container class with __slots__
        pass
    return state


def state_of(container: Any) -> State:
    state = getattr(container, "_ugc_wgw_inrun", None)
    if not isinstance(state, State):   # the task plugin did not run for this task: no driver factors, the call's id
        state = prepare(getattr(container, "cfg", None), container, str(getattr(container, "run_id", "") or "task"))
    return state


def _tail(path: str) -> List[str]:
    try:
        size = os.path.getsize(path)
        with open(path, "rb") as fh:
            if size > TAIL_BYTES:
                fh.seek(size - TAIL_BYTES)
            return fh.read().decode("utf-8", errors="replace").splitlines()
    except OSError:
        return []


def job_of(cli_log: str) -> Tuple[str, str]:
    """(job id, cluster) from the try's sbatch log (`--parsable`: "id" or "id;cluster" on the first line); ('', '')
    when sbatch printed no id."""
    try:
        with open(cli_log, "rt", errors="replace") as fh:
            first = fh.readline().strip()
    except OSError:
        return "", ""
    job, _, cluster = first.partition(";")
    return (job, cluster.strip()) if job.isdigit() else ("", "")


def job_outputs(host_dir: str, job: str, since: float) -> List[str]:
    """The job's own output file; without a job id, the ones written since the try began."""
    if job:
        return [os.path.join(host_dir, f"slurm-{job}.out")]
    out = []
    for path in sorted(glob.glob(os.path.join(host_dir, "slurm-*.out"))):
        try:
            if os.path.getmtime(path) >= since - 2:
                out.append(path)
        except OSError:
            continue
    return out


def job_state(job: str, cluster: str, runner: Callable[..., Any] = subprocess.run) -> str:
    """The job's state from the controller, '' when it cannot be asked."""
    if not job:
        return ""
    cmd = ["scontrol", "show", "job", "-o", job] + (["-M", cluster] if cluster else [])
    try:
        res = runner(cmd, capture_output=True, text=True, timeout=30)
    except (OSError, subprocess.SubprocessError):
        return ""
    if getattr(res, "returncode", 1) != 0:
        return ""
    m = JOB_STATE.search(str(getattr(res, "stdout", "")))
    return m.group(1) if m else ""


def _in_job_output(paths: List[str]) -> Optional[Tuple[str, str]]:
    for path in paths:
        lines = _tail(path)
        for pattern, resource in ((TIME_TEXT, "time"), (OOM_TEXT, "memory")):   # a job out of time is also killed
            for line in lines:
                if pattern.search(line):
                    return resource, line.strip()[:200]
    return None


def kill_reason(host_dir: str, exit_code: int, job: str, cluster: str, streams: List[str], since: float,
                runner: Callable[..., Any] = subprocess.run, sleep: Callable[[float], None] = time.sleep,
                ) -> Optional[Tuple[str, str]]:
    """("time" | "memory", the evidence) when this try's job was killed for a resource, else None."""
    paths = job_outputs(host_dir, job, since)
    hit = _in_job_output(paths)
    if hit:
        return hit
    state = job_state(job, cluster, runner)
    if state in STATE_RESOURCE:
        return STATE_RESOURCE[state], f"job {job} state {state}"
    if exit_code in MEMORY_EXIT:
        return "memory", f"exit status {exit_code}"
    for path in streams:
        for line in _tail(path):
            if PINNED_TEXT.search(line):
                return "memory", line.strip()[:200]
    if not state:   # the controller did not answer: give the job's output file time to arrive
        for _ in range(RECHECKS):
            sleep(RECHECK_SECONDS)
            hit = _in_job_output(job_outputs(host_dir, job, since))
            if hit:
                return hit
    return None


def _stop(state: State, logger: logging.Logger, exit_code: int, resource: str, code: str, why: str, evidence: str) -> None:
    state.stop = (exit_code, why)
    logger.warning(resources._structured(
        KILL_STOP_MESSAGE, task=state.task, resource=resource, code=code, reason=why, evidence=evidence,
        n=state.launched.get(resource, 0) + state.steps.get(resource, 0), max=state.max))


def escalate(container: Any, logger: logging.Logger, state: State, exit_code: int, resource: str, evidence: str,
             job: str = "") -> bool:
    """Raise the container's request for `resource` for miniwdl's next try; False (and no further job) when no
    escalation is left."""
    done = state.launched.get(resource, 0) + state.steps.get(resource, 0)
    if done >= state.max:
        _stop(state, logger, exit_code, resource, "exhausted",
              f"escalation exhausted ({done} of {state.max})" if state.max else "resource_retry_max is 0", evidence)
        return False
    cfg = getattr(container, "cfg", None)
    values = container.runtime_values
    cap = resources.memory_cap(cfg, resources._limits(cfg, logger, container)) if resource == "memory" else None
    change, capped = resources.raise_request(values, resource, cap, factor=state.factor)
    if change is None:
        if capped:
            _stop(state, logger, exit_code, resource, "capped",
                  f"memory is at the site's cap ({policy.format_memory(int(values[RUNTIME_KEY[resource]]))})", evidence)
        else:
            _stop(state, logger, exit_code, resource, "no_request",
                  "the task has no wall time request to raise" if resource == "time" else "the task has no memory request to raise",
                  evidence)
        return False
    state.steps[resource] = state.steps.get(resource, 0) + 1
    changes = {RUNTIME_KEY[resource]: change}
    logger.log(resources.NOTICE_LEVEL, resources._structured(
        resources.ESCALATION_MESSAGE, task=state.task, summary=policy.describe(changes),
        changes={k: list(v) for k, v in changes.items()}, in_run=True, resource=resource, evidence=evidence, job=job,
        n=done + 1, max=state.max, run_factor=state.factor ** state.steps[resource]))
    return True


def after_failure(container: Any, logger: logging.Logger, exit_code: int, since: float,
                  runner: Callable[..., Any] = subprocess.run, sleep: Callable[[float], None] = time.sleep) -> Optional[str]:
    """Called with a try's failing exit status, before miniwdl decides about a retry. Returns the resource the try
    was killed for (None for any other failure, a request sbatch refused among them)."""
    state = state_of(container)
    host_dir = str(container.host_dir)
    cli_log = os.path.join(host_dir, f"{container.cli_name}.log.txt")
    job, cluster = job_of(cli_log)
    if not job:
        for line in _tail(cli_log):
            if REFUSED_TEXT.search(line):
                _stop(state, logger, exit_code, "request", "refused", f"sbatch refused the request ({line.strip()[:160]})",
                      line.strip()[:200])
                return None
    streams = [container.host_stderr_txt(), container.host_stdout_txt()]
    hit = kill_reason(host_dir, exit_code, job, cluster, streams, since, runner, sleep)
    if hit is None:
        return None
    resource, evidence = hit
    escalate(container, logger, state, exit_code, resource, evidence, job)
    return resource


def skip(container: Any, logger: logging.Logger) -> Optional[int]:
    """The exit status to answer a try with instead of submitting it: set once a kill could not be escalated."""
    state = getattr(container, "_ugc_wgw_inrun", None)
    if not isinstance(state, State) or state.stop is None:
        return None
    exit_code, why = state.stop
    state.skipped += 1
    note = f"[ugc-wgw] not submitted again: an earlier try of this task was killed and {why}\n"
    for path, text in ((container.host_stderr_txt(), note), (container.host_stdout_txt(), "")):
        try:
            with open(path, "a") as fh:
                fh.write(text)
        except OSError:
            pass
    logger.warning(resources._structured(SKIPPED_MESSAGE, task=state.task, reason=why, exit_status=exit_code))
    return exit_code
