"""miniwdl task plugin: apply the site resource policy to a task's evaluated runtime values.

Registered as `miniwdl.plugin.task` entry point `ugc_wgw_resources`. miniwdl calls the generator once
per task attempt, after the call-cache lookup: stage 1 sees the inputs, stage 2 sees the command
and the TaskContainer with `runtime_values` already evaluated and clamped (cpu_max, memory_max)
but not yet submitted, stage 3 sees the outputs. Only stage 2 does anything here.

Configuration: `[ugc_wgw] resources = <path>` in miniwdl.cfg (rendered by scripts/install-bundle.sh).
Unset or empty means no policy; a missing file is an empty policy. See docs/guide/12-resources.md.

GPU tasks (runtime `gpuCount`, `gpuType`): upstream's HPC backend fills `gpuType` with "" when no type is
given, which miniwdl-slurm would render as `--gres gpu::N`; the empty type is dropped here so the request
becomes `gpu:N` (or `gpu:<type>:N` when a type is set). A policy row that moves a GPU task to a partition
also sets `slurm_partition_gpu`, the key miniwdl-slurm prefers for tasks with `gpuCount` (it is otherwise
the site default `SLURM_PARTITION_GPU`, which would win over the row).
"""
from __future__ import annotations

import json
import logging
import os
import re
from typing import Any, Dict, Iterator, List, Optional

from . import policy

CFG_SECTION = "ugc_wgw"
CFG_KEY = "resources"
NOTICE_LEVEL = 25
VERBOSE_LEVEL = 15
MESSAGE = "ugc-wgw resource policy applied"
_PARTITION_ARG = re.compile(r"(^|\s)(-p|--partition)(\s|=|$)")
_CONSTRAINT_ARG = re.compile(r"(^|\s)(-C|--constraint)(\s|=|$)")


class _Message:
    """Stand-in for WDL._util.StructuredLogMessage when miniwdl is not importable (unit tests)."""

    def __init__(self, _message: str, **kwargs: Any) -> None:
        self.message = _message
        self.kwargs = kwargs

    def __str__(self) -> str:
        return f"{self.message} :: {', '.join(k + ': ' + json.dumps(v) for k, v in self.kwargs.items())}"


def _structured(message: str, **kwargs: Any) -> object:
    try:
        from WDL._util import StructuredLogMessage  # type: ignore
    except ImportError:  # pragma: no cover - only outside miniwdl
        return _Message(message, **kwargs)
    return StructuredLogMessage(message, **kwargs)


def cfg_get(cfg: Any, section: str, key: str, default: str = "") -> str:
    try:
        value = cfg.get(section, key, default)
    except Exception:  # a section-less fake or ConfigMissing
        return default
    return default if value is None else str(value)


def policy_path(cfg: Any) -> str:
    return cfg_get(cfg, CFG_SECTION, CFG_KEY, "").strip()


def gpu_request(runtime_values: Dict[str, Any]) -> Optional[str]:
    """The gres miniwdl-slurm will ask for (`gpu:N` or `gpu:<type>:N`), None for a task without gpuCount.
    Drops an empty `gpuType` in place: upstream's HPC backend sets "" when no type is configured, and
    `--gres gpu::N` is not a request SLURM accepts."""
    gpu_type = runtime_values.get("gpuType")
    if gpu_type is not None and str(gpu_type).strip() == "":
        runtime_values.pop("gpuType")
        gpu_type = None
    count = runtime_values.get("gpuCount")
    if count is None:
        return None
    return f"gpu:{gpu_type}:{count}" if gpu_type else f"gpu:{count}"


def apply_policy(cfg: Any, logger: logging.Logger, task_name: str, container: Any) -> Dict[str, Any]:
    """Apply the configured policy to `container.runtime_values`; returns the changes (empty when none).
    Also normalises the GPU request (gpu_request) and logs it, policy or not."""
    gres = gpu_request(container.runtime_values)
    changes = _apply_rows(cfg, logger, task_name, container, gres)
    if gres is not None:
        logger.log(NOTICE_LEVEL, _structured("ugc-wgw gpu request", task=task_name, gres=gres,
                                             partition=container.runtime_values.get("slurm_partition_gpu",
                                                                                    container.runtime_values.get("slurm_partition"))))
    return changes


def _apply_rows(cfg: Any, logger: logging.Logger, task_name: str, container: Any, gres: Optional[str]) -> Dict[str, Any]:
    path = policy_path(cfg)
    if not path:
        return {}
    if not os.path.exists(path):
        logger.log(VERBOSE_LEVEL, _structured("ugc-wgw resource policy missing", policy=path))
        return {}
    rules = policy.load_policy(path)  # PolicyError propagates: the task fails with file:line
    if not rules:
        logger.log(VERBOSE_LEVEL, _structured("ugc-wgw resource policy empty", policy=path))
        return {}
    limits = None
    detect = getattr(container, "detect_resource_limits", None)
    if callable(detect):
        try:
            limits = detect(cfg, logger)
        except Exception:  # pragma: no cover - a backend without limits
            limits = None
    changes, sources = policy.apply(rules, task_name, container.runtime_values, limits)
    if not changes:
        return {}
    if gres is not None and "slurm_partition" in changes:
        # miniwdl-slurm submits a task with gpuCount to slurm_partition_gpu when that key exists
        container.runtime_values["slurm_partition_gpu"] = changes["slurm_partition"][1]
    extra_args = cfg_get(cfg, "slurm", "extra_args", "")
    if "slurm_partition" in changes and _PARTITION_ARG.search(extra_args):
        logger.warning(_structured("ugc-wgw resource policy: partition overridden by [slurm] extra_args",
                                   task=task_name, extra_args=extra_args,
                                   hint="move --partition from SLURM_EXTRA_ARGS to SLURM_PARTITION in site.cfg"))
    if "slurm_constraint" in changes and _CONSTRAINT_ARG.search(extra_args):
        logger.warning(_structured("ugc-wgw resource policy: constraint overridden by [slurm] extra_args",
                                   task=task_name, extra_args=extra_args))
    logger.log(NOTICE_LEVEL, _structured(MESSAGE, task=task_name, summary=policy.describe(changes),
                                         changes={k: list(v) for k, v in changes.items()},
                                         rules=sources, policy=path))
    return changes


def task(cfg: Any, logger: logging.Logger, run_id_stack: List[str], run_dir: str, task: Any,
         **recv: Any) -> Iterator[Dict[str, Any]]:
    """The miniwdl task plugin coroutine (see the module docstring for the three stages)."""
    recv = yield recv                      # stage 1: inputs, untouched
    container = recv.get("container")
    if container is not None:
        apply_policy(cfg, logger, str(task.name), container)
    recv = yield recv                      # stage 2: command and container
    yield recv                             # stage 3: outputs, untouched
