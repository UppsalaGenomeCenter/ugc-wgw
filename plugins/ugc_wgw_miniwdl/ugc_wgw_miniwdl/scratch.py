"""Node-local scratch for tasks: miniwdl-slurm's backend with the job's scratch directory bound in.

Registered as `miniwdl.plugin.container_backend` entry point `ugc_slurm_singularity` and selected by
`[scheduler] container_backend = ugc_slurm_singularity` in the rendered miniwdl.cfg. It is miniwdl-slurm's
`slurm_singularity` (sbatch --wait per task, Apptainer on the node) with one change: when `[ugc_wgw]
scratch_env` names a job environment variable (`TMPDIR`, `SNIC_TMP`), the sbatch script that runs the
container is `scripts/exec_scratch.sh` instead of miniwdl-slurm's `exec "$@"`. On the node, where the
variable is known, that script binds the directory it names at `[ugc_wgw] scratch_mount` (default
/mnt/scratch) and sets TMPDIR to it inside the container, so tools that honour TMPDIR keep their
temporaries off the shared file system. Everything else (the sbatch arguments, the log file the driver
reads, cancellation) is miniwdl-slurm's. With `scratch_env` empty the backend behaves exactly like
`slurm_singularity`.

Why here and not in the configuration: Apptainer's `--containall` (our run_options) clears the environment,
so TMPDIR and SLURM_JOB_ID are unknown inside the container, and miniwdl-slurm passes the container command
to sbatch as literal script arguments, so a `$TMPDIR` in `run_options` would never be expanded. Only a
script running on the node can resolve the per-job directory. See docs/guide/04-install.md.
"""
from __future__ import annotations

import os
from typing import Any, List, Tuple

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

else:

    class SlurmSingularityScratch:  # type: ignore[no-redef]
        """Placeholder so the entry point loads in a venv without miniwdl-slurm; selecting it is an error."""

        def __init__(self, *args: Any, **kwargs: Any) -> None:
            raise RuntimeError(f"container backend {BACKEND_NAME} needs miniwdl-slurm, which is not installed")
