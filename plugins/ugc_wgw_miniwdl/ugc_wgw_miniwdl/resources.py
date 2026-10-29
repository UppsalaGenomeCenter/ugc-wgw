"""miniwdl task plugin: apply the site resource policy to a task's evaluated runtime values.

Registered as `miniwdl.plugin.task` entry point `ugc_wgw_resources`. miniwdl calls the generator once
per task attempt, after the call-cache lookup: stage 1 sees the inputs, stage 2 sees the command
and the TaskContainer with `runtime_values` already evaluated and clamped (cpu_max, memory_max)
but not yet submitted, stage 3 sees the outputs. Only stage 2 does anything here.

Configuration: `[ugc_wgw] resources = <path>` in miniwdl.cfg (rendered by scripts/install-bundle.sh).
Unset or empty means no policy; a missing file is an empty policy. See docs/guide/12-resources.md.

GPU tasks (runtime `gpuCount`, `gpuType`): the gres type comes from `[ugc_wgw] gpu_type`, which the driver
passes per run through miniwdl's environment override (MINIWDL__UGC_WGW__GPU_TYPE, the project's config.json
`gpu_type`), not through a workflow input: upstream's `RuntimeAttributes.gpuType` is an input of every task,
so a type given that way changes every task's call-cache key and switching GPU nodes re-aligns everything.
A non-empty `gpuType` the WDL set itself (an explicit stage input) wins. Upstream's HPC backend fills
`gpuType` with "" when no type is given, which miniwdl-slurm would render as `--gres gpu::N`; an empty type
is dropped so the request becomes `gpu:N` (any GPU of the partition) or `gpu:<type>:N`. A policy row that
moves a GPU task to a partition also sets `slurm_partition_gpu`, the key miniwdl-slurm prefers for tasks with
`gpuCount` (it is otherwise the site default `SLURM_PARTITION_GPU`, which would win over the row).
"""
from __future__ import annotations

import json
import math
import logging
import os
import re
from typing import Any, Dict, Iterator, List, Optional

from . import crossdev, policy  # noqa: F401  (importing crossdev installs the hardlink-else-symlink fallback)

CFG_SECTION = "ugc_wgw"
CFG_KEY = "resources"
GPU_TYPE_KEY = "gpu_type"
NOTICE_LEVEL = 25
VERBOSE_LEVEL = 15
MESSAGE = "ugc-wgw resource policy applied"
ESCALATION_MESSAGE = "ugc-wgw resource escalation applied"
ESCALATIONS_KEY = "escalations"   # [ugc_wgw] escalations, set per run by the driver as MINIWDL__UGC_WGW__ESCALATIONS
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


def gpu_type(cfg: Any) -> str:
    """`[ugc_wgw] gpu_type`: the SLURM gres type for tasks that ask for GPUs ('' = any type in the partition)."""
    return cfg_get(cfg, CFG_SECTION, GPU_TYPE_KEY, "").strip()


def apply_gpu_type(cfg: Any, runtime_values: Dict[str, Any]) -> Optional[str]:
    """Set `gpuType` from the configuration on a task with `gpuCount` whose WDL gave none (or ""); returns the
    type applied, None when nothing changed. Read for every task so the driver's per-run override counts as used."""
    configured = gpu_type(cfg)
    if "gpuCount" not in runtime_values or not configured:
        return None
    current = runtime_values.get("gpuType")
    if current is not None and str(current).strip():
        return None   # the workflow (an explicit stage input) set a type: it wins
    runtime_values["gpuType"] = configured
    return configured


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
    Also applies the configured GPU type (apply_gpu_type), normalises the GPU request (gpu_request) and logs
    it, policy or not."""
    applied = apply_gpu_type(cfg, container.runtime_values)
    gres = gpu_request(container.runtime_values)
    changes = _apply_rows(cfg, logger, task_name, container, gres)
    changes.update(apply_escalations(cfg, logger, task_name, container))
    if gres is not None:
        logger.log(NOTICE_LEVEL, _structured("ugc-wgw gpu request", task=task_name, gres=gres,
                                             gpu_type_from="config" if applied else "workflow",
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
    changes, sources = policy.apply(rules, task_name, container.runtime_values, _limits(cfg, logger, container))
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


def _limits(cfg: Any, logger: logging.Logger, container: Any) -> Optional[Dict[str, int]]:
    """The backend's detect_resource_limits() answer ({cpu, mem_bytes}), None without one."""
    detect = getattr(container, "detect_resource_limits", None)
    if not callable(detect):
        return None
    try:
        return detect(cfg, logger)
    except Exception:  # pragma: no cover - a backend without limits
        return None


def apply_escalations(cfg: Any, logger: logging.Logger, task_name: str, container: Any) -> Dict[str, Any]:
    """The driver's escalations for this task (bin/ugc_wgw/escalate.py), after the policy rows: `{task:
    {time_factor, memory_factor, time_minutes_min, memory_min}}` as JSON in [ugc_wgw] escalations (the environment
    override MINIWDL__UGC_WGW__ESCALATIONS). A factor multiplies the task's current wall time or memory, a minimum
    lifts it to at least that value (a learned request); memory never exceeds the backend's limit, and a request
    already at that limit is logged as capped. Returns the changes."""
    raw = cfg_get(cfg, CFG_SECTION, ESCALATIONS_KEY, "").strip()
    if not raw:
        return {}
    try:
        doc = json.loads(raw)
    except ValueError as exc:
        logger.warning(_structured("ugc-wgw resource escalations unreadable", error=str(exc)))
        return {}
    spec = doc.get(task_name) if isinstance(doc, dict) else None
    if not isinstance(spec, dict) or not spec:
        return {}
    limits = _limits(cfg, logger, container)
    values = container.runtime_values
    changes: Dict[str, Any] = {}
    capped: List[str] = []
    for resource, key in (("time", "time_minutes"), ("memory", "memory_reservation")):
        old = values.get(key)
        if not isinstance(old, int) or old <= 0:
            continue
        new = old
        factor = spec.get(f"{resource}_factor")
        if isinstance(factor, (int, float)) and factor > 1:
            new = int(math.ceil(old * float(factor)))
        minimum = spec.get(f"{resource}_min")
        if isinstance(minimum, (int, float)) and int(minimum) > new:
            new = int(minimum)
        if resource == "memory" and limits and limits.get("mem_bytes", 0) > 0 and new > int(limits["mem_bytes"]):
            new = int(limits["mem_bytes"])
            if new <= old:
                capped.append(policy.format_memory(old))
        if new == old:
            continue
        if resource == "memory" and isinstance(values.get("memory_limit"), int):
            values["memory_limit"] = int(values["memory_limit"] * new / old)
        values[key] = new
        changes[key] = (old, new)
    if capped:
        logger.warning(_structured("ugc-wgw resource escalation capped at the node's memory", task=task_name,
                                   memory=capped[0], hint="a policy row with a larger partition (docs/guide/12-resources.md)"))
    if changes:
        logger.log(NOTICE_LEVEL, _structured(ESCALATION_MESSAGE, task=task_name, summary=policy.describe(changes),
                                             changes={k: list(v) for k, v in changes.items()}, spec=spec))
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
