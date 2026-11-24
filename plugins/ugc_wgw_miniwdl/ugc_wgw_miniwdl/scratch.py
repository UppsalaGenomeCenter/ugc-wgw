"""Node-local scratch for tasks: miniwdl-slurm's backend with the job's scratch directory bound in.

Registered as `miniwdl.plugin.container_backend` entry point `ugc_slurm_singularity` and selected by
`[scheduler] container_backend = ugc_slurm_singularity` in the rendered miniwdl.cfg. It is miniwdl-slurm's
`slurm_singularity` (sbatch --wait per task, Apptainer on the node) with one change: when `[ugc_wgw]
scratch_env` names a job environment variable (`TMPDIR`, `SNIC_TMP`), the sbatch script that runs the
container is `scripts/exec_scratch.sh` instead of miniwdl-slurm's `exec "$@"`. On the node, where the
variable is known, that script binds the directory it names at `[ugc_wgw] scratch_mount` (default
/mnt/scratch) and sets TMPDIR to it inside the container, so tools that honour TMPDIR keep their
temporaries off the shared file system. Everything else (the sbatch arguments, the log file the driver
reads, cancellation) is miniwdl-slurm's. With `scratch_env` empty the exec script is miniwdl-slurm's own.

The backend is also where a try that was killed for wall time or memory gets a larger request before
miniwdl retries it (inrun.py): `_run` is the one place that sees every try's exit status.

Why here and not in the configuration: Apptainer's `--containall` (our run_options) clears the environment,
so TMPDIR and SLURM_JOB_ID are unknown inside the container, and miniwdl-slurm passes the container command
to sbatch as literal script arguments, so a `$TMPDIR` in `run_options` would never be expanded. Only a
script running on the node can resolve the per-job directory. See docs/guide/04-install.md.
"""
from __future__ import annotations

import logging
import os
import time
from typing import Any, Callable, List, Tuple

from . import inrun
from .resources import cfg_get

CFG_SECTION = "ugc_wgw"
ENV_KEY = "scratch_env"
MOUNT_KEY = "scratch_mount"
DEFAULT_MOUNT = "/mnt/scratch"
BACKEND_NAME = "ugc_slurm_singularity"
EXEC_SCRIPT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "scripts", "exec_scratch.sh")
STOCK_SCRIPT = "exec_script.sh"   # miniwdl-slurm's, always the last sbatch argument


def scratch_settings(cfg: Any) -> Tuple[str, str]:
    """(variable name, mount point) from the miniwdl configuration; ('', mount) when scratch is off."""
    env = cfg_get(cfg, CFG_SECTION, ENV_KEY, "").strip()
    mount = cfg_get(cfg, CFG_SECTION, MOUNT_KEY, DEFAULT_MOUNT).strip() or DEFAULT_MOUNT
    return env, mount


def with_exec_script(sbatch_args: List[str], env: str, mount: str) -> List[str]:
    """miniwdl-slurm's sbatch arguments with its trailing exec script replaced by ours and its parameters;
    unchanged when `env` is empty."""
    if not env:
        return list(sbatch_args)
    if not sbatch_args or not sbatch_args[-1].endswith(STOCK_SCRIPT):
        raise RuntimeError("miniwdl-slurm's sbatch arguments do not end with its exec script; "
                           f"{BACKEND_NAME} expects miniwdl-slurm 0.4")
    return list(sbatch_args[:-1]) + [EXEC_SCRIPT, env, mount, "--"]


def run_try(container: Any, logger: logging.Logger, terminating: Callable[[], bool], submit: Callable[[], int]) -> int:
    """One try of a task: `submit` runs the job and returns its exit status. A try that follows a kill no
    escalation was left for is answered without a job; a failed try is looked at before miniwdl decides about
    the next one (inrun.after_failure), and nothing found there can change the try's own result."""
    skipped = inrun.skip(container, logger)
    if skipped is not None:
        return skipped
    since = time.time()
    exit_code = submit()
    if not terminating() and not container.success_exit_code(exit_code):
        try:
            inrun.after_failure(container, logger, exit_code, since)
        except Exception as exc:  # the try failed either way; miniwdl's own retry goes ahead
            logger.warning(f"ugc-wgw: the failed try could not be examined for a resource kill: {exc!r}")
    return exit_code


try:
    from miniwdl_slurm import SlurmSingularity  # type: ignore
except ImportError:  # pragma: no cover - a venv without miniwdl-slurm never selects this backend
    SlurmSingularity = None  # type: ignore[assignment,misc]


if SlurmSingularity is not None:

    class SlurmSingularityScratch(SlurmSingularity):  # type: ignore[misc,valid-type]
        """miniwdl-slurm's backend; the exec script is ours when [ugc_wgw] scratch_env is set."""

        def _slurm_invocation(self) -> List[str]:
            env, mount = scratch_settings(self.cfg)
            return with_exec_script(super()._slurm_invocation(), env, mount)

        def _run(self, logger: logging.Logger, terminating: Callable[[], bool], command: str) -> int:
            return run_try(self, logger, terminating, lambda: super(SlurmSingularityScratch, self)._run(logger, terminating, command))

else:

    class SlurmSingularityScratch:  # type: ignore[no-redef]
        """Placeholder so the entry point loads in a venv without miniwdl-slurm; selecting it is an error."""

        def __init__(self, *args: Any, **kwargs: Any) -> None:
            raise RuntimeError(f"container backend {BACKEND_NAME} needs miniwdl-slurm, which is not installed")
