"""Resource escalation and learned resources (docs/DESIGN.md §8.1, decision log 2026-10-29; guide chapters 06
and 12).

A run that failed because one task ran out of wall time or memory (failures.classify: kind `resource`,
resource `time` or `memory`) is re-attempted by the submit loop without backoff, with that task's limit
multiplied by `resource_retry_factor` (2), up to `resource_retry_max` (2) escalations per task and resource.
The multiplier is not a workflow input: the driver passes it to the ugc_wgw_resources plugin through miniwdl's
environment override MINIWDL__UGC_WGW__ESCALATIONS, and the plugin applies it after the policy rows, so the
task's call-cache key is unchanged and every other task of the run replays from the cache. Memory stops at the
site cap (the plugin clamps to the backend's limit); a failure at the cap blocks the subject with a message
that names the policy row it needs.

What an escalated task needed is remembered: when a task that ran with an escalation finished, the request it
finished with goes into the project's `.ugc-wgw/learned_resources.json` (task -> time_minutes, memory), and
every later run of the project starts that task at the learned value (passed as a minimum in the same
environment variable). `ugc-wgw resources` shows learned values as `(learned)` and prints them as policy rows
with `--learned`, for promotion to the site policy.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

from .config import Config
from .db import RunRecord
from .failures import Classification
from .layout import run_files
from .log import LOGGER, Events
from .util import read_json, utc_now, write_json
from . import wdllog

ESCALATION_ENV = "MINIWDL__UGC_WGW__ESCALATIONS"
ESCALATION_MESSAGE = "ugc-wgw resource escalation applied"   # the plugin's NOTICE (plugins/ugc_wgw_miniwdl)
RUNTIME_KEY = {"time": "time_minutes", "memory": "memory_reservation"}   # resource -> miniwdl runtime key
LEARNED_KEY = {"time_minutes": "time_minutes", "memory_reservation": "memory"}   # runtime key -> learned store key


def task_name_for(run_path: Path, node: str | None) -> str | None:
    """The WDL task name behind a failed call id (`call-ugc_wgw_pav`, `call-glnexus-00`, `call-concat_svx`), from
    the run's workflow.log.json (`task setup` names the task); None when the log does not have it."""
    if not node:
        return None
    tasks, _ = wdllog.parse_workflow_log(run_files(run_path).workflow_log_json)
    for t in tasks:
        if t.call_id == node or t.call_path.endswith(node):
            return t.name
    return None


def escalations_of(run: RunRecord) -> dict[str, dict[str, float]]:
    """What the run was launched with: task -> {time_factor, memory_factor}."""
    raw = run.meta.get("escalations")
    return {str(k): dict(v) for k, v in raw.items()} if isinstance(raw, dict) else {}


def plan(cfg: Config, run: RunRecord, cls: Classification) -> tuple[dict[str, object] | None, str]:
    """The escalation the next attempt should carry after this run's resource failure, and a short reason when
    there is none ("" when there is one): {task, resource, factor, n, max}."""
    if cls.kind != "resource":
        return None, ""
    if cls.resource not in ("time", "memory"):
        return None, "not escalated: the request was refused or the resource is unknown"
    if cfg.resource_retry_max <= 0:
        return None, "not escalated: resource_retry_max is 0"
    node = cls.evidence.get("node")
    task = task_name_for(run.run_path, str(node) if node else None)
    if not task:
        return None, "not escalated: the failed task is unknown"
    so_far = float(escalations_of(run).get(task, {}).get(f"{cls.resource}_factor", 1.0))
    base = max(1.000001, float(cfg.resource_retry_factor))
    n_so_far = int(round(math.log(so_far, base))) if so_far > 1.0 else 0
    if n_so_far >= cfg.resource_retry_max:
        return None, (f"escalation exhausted: {task} {cls.resource} already ×{so_far:g} ({n_so_far} of "
                      f"{cfg.resource_retry_max}); a policy row or a stage input is needed (guide chapter 12)")
    return {"task": task, "resource": cls.resource, "factor": so_far * base, "n": n_so_far + 1,
            "max": cfg.resource_retry_max}, ""


def describe(esc: dict[str, object]) -> str:
    return f"{esc['task']} {esc['resource']} ×{float(esc['factor']):g} (escalation {esc['n']} of {esc['max']})"


def carry(previous: RunRecord | None) -> dict[str, dict[str, float]]:
    """The escalations for the attempt after `previous`: what it ran with, plus the one its failure asked for."""
    if previous is None:
        return {}
    out = escalations_of(previous)
    nxt = previous.meta.get("escalation")
    if isinstance(nxt, dict) and nxt.get("task") and nxt.get("resource") in RUNTIME_KEY:
        out.setdefault(str(nxt["task"]), {})[f"{nxt['resource']}_factor"] = float(nxt["factor"])  # type: ignore[arg-type]
    return out


def load_learned(cfg: Config) -> dict[str, dict[str, object]]:
    path = cfg.learned_path
    if not path.exists():
        return {}
    try:
        doc = read_json(path)
    except (OSError, ValueError) as exc:
        LOGGER.warning("learned resources unreadable (%s): %s", path, exc)
        return {}
    return {str(k): dict(v) for k, v in (doc.get("tasks") or {}).items()} if isinstance(doc, dict) else {}


def save_learned(cfg: Config, tasks: dict[str, dict[str, object]]) -> None:
    write_json(cfg.learned_path, {"schema": 1, "updated": utc_now(), "tasks": dict(sorted(tasks.items()))})


def env_for(cfg: Config, run: RunRecord) -> dict[str, str]:
    """The plugin's environment for this run: the project's learned values as minimums and the run's escalation
    factors; {} when there is neither."""
    spec: dict[str, dict[str, object]] = {}
    for task, learned in load_learned(cfg).items():
        for key, value in learned.items():
            if key in ("time_minutes", "memory") and isinstance(value, int) and value > 0:
                spec.setdefault(task, {})[f"{'time' if key == 'time_minutes' else 'memory'}_min"] = value
    for task, factors in escalations_of(run).items():
        for key, value in factors.items():
            if key in ("time_factor", "memory_factor") and float(value) > 1.0:
                spec.setdefault(task, {})[key] = float(value)
    return {ESCALATION_ENV: json.dumps(spec, sort_keys=True)} if spec else {}


def learn(cfg: Config, events: Events, run: RunRecord) -> list[str]:
    """Record what the run's escalated tasks finished with. Every task that the plugin escalated and that
    completed (whatever the run's own outcome) teaches the project its request; returns the lessons."""
    if not cfg.learn_resources or not escalations_of(run):
        return []
    tasks, _ = wdllog.parse_workflow_log(run_files(run.run_path).workflow_log_json)
    learned = load_learned(cfg)
    lessons: list[str] = []
    for t in tasks:
        if not t.escalation_changes or t.end is None or t.failed:
            continue
        for runtime_key, pair in t.escalation_changes.items():
            key = LEARNED_KEY.get(runtime_key)
            if key is None:
                continue
            try:
                new = int(pair[1])
            except (TypeError, ValueError):
                continue
            entry = learned.setdefault(t.name, {})
            if isinstance(entry.get(key), int) and int(entry[key]) >= new:  # type: ignore[call-overload]
                continue
            entry[key] = new
            entry["run_id"] = run.run_id
            entry["learned_at"] = utc_now()
            lessons.append(f"{t.name} {key} {new}")
    if lessons:
        save_learned(cfg, learned)
        events.emit("resources.learned", run_id=run.run_id, subject=run.subject_id, stage=run.stage,
                    lessons=lessons, path=str(cfg.learned_path))
    return lessons


def policy_rows(cfg: Config, policy: object) -> list[str]:
    """The learned values as rows of the site policy (`ugc-wgw resources --learned`), for promotion."""
    rows = []
    for task, learned in load_learned(cfg).items():
        mem = learned.get("memory")
        time = learned.get("time_minutes")
        if not isinstance(mem, int) and not isinstance(time, int):
            continue
        rows.append("\t".join([task, "-", policy.format_memory(mem) if isinstance(mem, int) else "-",  # type: ignore[attr-defined]
                               policy.format_time(time) if isinstance(time, int) else "-", "-", "-"]))  # type: ignore[attr-defined]
    return rows


__all__ = ["ESCALATION_ENV", "ESCALATION_MESSAGE", "plan", "describe", "carry", "env_for", "learn", "load_learned",
           "policy_rows", "task_name_for"]
