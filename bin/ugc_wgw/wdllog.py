"""miniwdl's `workflow.log.json`, read per task: timings, cache hits, retries, and the resource lines miniwdl and
the ugc_wgw_resources plugin write. Shared by report.py, progress.py, accounting.py and usage.py (guide chapter 07).

NOTICE lines pair `task setup` with `done` / `done (cached)` per task logger (`source`); the `task setup` line
also names the task (`name`) and its directory (`dir`), which is where miniwdl-slurm leaves the job id.
"""
from __future__ import annotations

import json
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


@dataclass
class TaskStat:
    name: str
    call_id: str
    start: float
    end: Optional[float] = None
    cached: bool = False
    cpu_requested: Optional[int] = None
    cpu_granted: Optional[int] = None
    policy: str = ""              # the ugc_wgw_resources plugin's summary, e.g. "cpu 64→48, partition -→fat"
    retries: int = 0
    failed: bool = False
    call_path: str = ""           # call ids from the entrypoint down, `call-sub-1/call-b-0` (unique within an attempt)
    dir: str = ""                 # the task directory miniwdl logged, when it did
    gres: Optional[str] = None    # the plugin's `ugc-wgw gpu request`, e.g. `gpu:a100:4`
    cpu_policy: Optional[int] = None   # cpu after the policy row (`changes.cpu[1]`)

    @property
    def duration(self) -> Optional[float]:
        return None if self.end is None else max(0.0, self.end - self.start)

    @property
    def cpu_launched(self) -> Optional[int]:
        """The cpu the task was launched with, as far as the log tells: the policy's value, else miniwdl's cap."""
        return self.cpu_policy if self.cpu_policy is not None else self.cpu_granted


def call_path_of(source: str) -> str:
    """`wdl.w:ugc_wgw_x.w:call-sub-1.t:call-b-0` -> `call-sub-1/call-b-0`: the call components of a task logger."""
    parts = []
    for seg in source.split("."):
        name = seg.split(":", 1)[1] if ":" in seg else seg
        if name.startswith("call-"):
            parts.append(name)
    return "/".join(parts)


def _int(v: object) -> Optional[int]:
    try:
        return int(v)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None


def parse_workflow_log(path: Path) -> tuple[list[TaskStat], Counter]:
    """Per-task timings from miniwdl's JSON-lines workflow log; the task logger name is the key."""
    tasks: dict[str, TaskStat] = {}
    ignored: Counter = Counter()
    try:
        fh = open(path, errors="replace")
    except OSError:
        return [], ignored
    with fh:
        for line in fh:
            try:
                d = json.loads(line)
            except ValueError:
                continue
            src = str(d.get("source", ""))
            msg = str(d.get("message", ""))
            ts = d.get("timestamp")
            if not isinstance(ts, (int, float)) or ".t:" not in src:
                continue
            if msg == "task setup":
                call_id = src.rsplit(".t:", 1)[-1]
                tasks[src] = TaskStat(str(d.get("name", call_id)), call_id, float(ts), call_path=call_path_of(src),
                                      dir=str(d.get("dir") or ""))
            elif src in tasks:
                t = tasks[src]
                if msg.startswith("done"):
                    t.end = float(ts)
                    t.cached = "cached" in msg
                elif msg == "runtime.cpu adjusted to host limit":
                    t.cpu_requested, t.cpu_granted = _int(d.get("original")), _int(d.get("adjusted"))
                elif msg == "ugc-wgw resource policy applied":
                    t.policy = str(d.get("summary") or "")
                    changes = d.get("changes")
                    if isinstance(changes, dict) and isinstance(changes.get("cpu"), list) and len(changes["cpu"]) == 2:
                        t.cpu_policy = _int(changes["cpu"][1])
                elif msg == "ugc-wgw gpu request":
                    t.gres = str(d.get("gres") or "") or None
                elif msg == "failed task will be retried":
                    t.retries += 1
                elif msg == "ignored runtime settings":
                    for k in d.get("keys") or []:
                        ignored[str(k)] += 1
                elif d.get("level") == "ERROR" and "failed" in msg:
                    t.failed = True
    return list(tasks.values()), ignored
