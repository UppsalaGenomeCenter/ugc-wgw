"""The per-task resource policy: a TSV the operator edits, applied to a task's evaluated runtime values.

Format (tab-separated, header required, `#` lines and blank lines ignored):

    task	cpu	memory	time	partition	constraint
    deepvariant_call_variants_cpu	48	-	-	-	-
    ugc_wgw_hifiasm_assemble	-	-	2-00:00:00	fat	-
    glnexus	-	200G	-	-	-
    deepvariant_*	-	-	12:00:00	-	-

`task` is a task name or a glob (fnmatch). `-` or an empty cell leaves that value as the task
declared it. Every row whose pattern matches a task applies, in file order; a later row's cells
override an earlier row's. Memory needs a unit; units are binary whatever the spelling (`G`,
`GB` and `GiB` all mean 2^30 bytes, as SLURM's --mem counts), so `256G` equals the WDL's
`"256 GiB"`. Time is SLURM's syntax: minutes, `MM:SS`, `HH:MM:SS`, `D-HH`, `D-HH:MM` or
`D-HH:MM:SS`, rounded up to whole minutes. Extra informational columns `stages`, `source`
and `notes` (as in the generated inventory backends/hpc/resources.declared.tsv) are accepted
and ignored, so inventory rows can be copied as they are.

Stdlib only: the driver (bin/ugc_wgw/resources.py) imports this module too, outside the venv.
"""
from __future__ import annotations

import csv
import fnmatch
import math
import os
import re
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Tuple

COLUMNS = ("task", "cpu", "memory", "time", "partition", "constraint")
INFO_COLUMNS = ("stages", "source", "notes")
VALUE_COLUMNS = COLUMNS[1:]
UNSET = frozenset({"", "-"})

# policy column -> key in miniwdl's container.runtime_values (miniwdl-slurm reads the last three)
RUNTIME_KEYS = {
    "cpu": "cpu",
    "memory": "memory_reservation",
    "time": "time_minutes",
    "partition": "slurm_partition",
    "constraint": "slurm_constraint",
}

KIB = 1024
BYTE_UNITS = {
    "B": 1,
    "K": KIB, "KB": KIB, "KI": KIB, "KIB": KIB,
    "M": KIB ** 2, "MB": KIB ** 2, "MI": KIB ** 2, "MIB": KIB ** 2,
    "G": KIB ** 3, "GB": KIB ** 3, "GI": KIB ** 3, "GIB": KIB ** 3,
    "T": KIB ** 4, "TB": KIB ** 4, "TI": KIB ** 4, "TIB": KIB ** 4,
}
_MEMORY_RE = re.compile(r"^(\d+(?:\.\d+)?)\s*([A-Za-z]+)$")
_TIME_DAYS_RE = re.compile(r"^(\d+)-(\d+)(?::(\d+))?(?::(\d+))?$")
_TIME_COLON_RE = re.compile(r"^(\d+):(\d+)(?::(\d+))?$")
_NAME_RE = re.compile(r"^[A-Za-z0-9_.*?\[\]!-]+$")


class PolicyError(ValueError):
    """A malformed policy file; the message names file and line."""


@dataclass(frozen=True)
class Rule:
    pattern: str
    cpu: Optional[int] = None
    memory: Optional[int] = None      # bytes
    time: Optional[int] = None        # minutes
    partition: Optional[str] = None
    constraint: Optional[str] = None
    source: str = ""                  # "file:line"

    def matches(self, task_name: str) -> bool:
        return fnmatch.fnmatchcase(task_name, self.pattern)

    def values(self) -> Dict[str, object]:
        """The cells this row sets, keyed by policy column."""
        out: Dict[str, object] = {}
        for col in VALUE_COLUMNS:
            v = getattr(self, col)
            if v is not None:
                out[col] = v
        return out


# ---- units -----------------------------------------------------------------

def parse_cpu(text: str) -> int:
    if not text.isdigit() or int(text) < 1:
        raise ValueError(f"cpu must be a positive integer, got {text!r}")
    return int(text)


def parse_memory(text: str) -> int:
    """'256G', '1.5T', '300GiB', '512M' -> bytes (binary units). A bare number is an error."""
    m = _MEMORY_RE.match(text.strip())
    if not m:
        raise ValueError(f"memory needs a number and a unit such as 256G or 1.5T, got {text!r}")
    number, unit = float(m.group(1)), m.group(2).upper()
    if unit not in BYTE_UNITS:
        raise ValueError(f"unknown memory unit {m.group(2)!r} in {text!r} (use K, M, G or T)")
    nbytes = int(number * BYTE_UNITS[unit])
    if nbytes < 1:
        raise ValueError(f"memory must be positive, got {text!r}")
    return nbytes


def parse_time(text: str) -> int:
    """SLURM time syntax -> whole minutes, rounded up."""
    t = text.strip()
    if t.isdigit():
        seconds = int(t) * 60
    else:
        m = _TIME_DAYS_RE.match(t)
        if m:
            days, hours, minutes, secs = (int(g) if g is not None else 0 for g in m.groups())
            seconds = ((days * 24 + hours) * 60 + minutes) * 60 + secs
        else:
            m = _TIME_COLON_RE.match(t)
            if not m:
                raise ValueError(f"time must be minutes, MM:SS, HH:MM:SS, D-HH, D-HH:MM or D-HH:MM:SS, got {text!r}")
            a, b, c = m.groups()
            if c is None:  # MM:SS, as sbatch reads it
                seconds = int(a) * 60 + int(b)
            else:  # HH:MM:SS
                seconds = (int(a) * 60 + int(b)) * 60 + int(c)
    minutes = math.ceil(seconds / 60)
    if minutes < 1:
        raise ValueError(f"time must be at least one minute, got {text!r}")
    return minutes


def format_memory(nbytes: int) -> str:
    for unit, size in (("T", KIB ** 4), ("G", KIB ** 3), ("M", KIB ** 2)):
        if nbytes % size == 0:
            return f"{nbytes // size}{unit}"
    return f"{nbytes / KIB ** 3:.2f}G"


def format_time(minutes: int) -> str:
    days, rem = divmod(minutes, 24 * 60)
    hours, mins = divmod(rem, 60)
    return f"{days}-{hours:02d}:{mins:02d}:00" if days else f"{hours:02d}:{mins:02d}:00"


_PARSERS = {"cpu": parse_cpu, "memory": parse_memory, "time": parse_time,
            "partition": lambda s: s, "constraint": lambda s: s}


# ---- parsing ---------------------------------------------------------------

def parse_policy_text(text: str, name: str = "<policy>") -> List[Rule]:
    lines = text.splitlines()
    header: Optional[List[str]] = None
    header_line = 0
    rules: List[Rule] = []
    for lineno, raw in enumerate(lines, 1):
        line = raw.rstrip("\r\n")
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        cells = next(csv.reader([line], delimiter="\t"))
        cells = [c.strip() for c in cells]
        if header is None:
            header = [c.lower() for c in cells]
            header_line = lineno
            if not header or header[0] != "task":
                raise PolicyError(f"{name}:{lineno}: the first column of the header must be 'task'")
            unknown = [c for c in header if c not in COLUMNS and c not in INFO_COLUMNS]
            if unknown:
                raise PolicyError(f"{name}:{lineno}: unknown column(s) {', '.join(unknown)}; "
                                  f"known: {', '.join(COLUMNS + INFO_COLUMNS)}")
            if len(set(header)) != len(header):
                raise PolicyError(f"{name}:{lineno}: duplicate column in the header")
            continue
        if len(cells) > len(header):
            raise PolicyError(f"{name}:{lineno}: {len(cells)} cells for {len(header)} columns "
                              f"(header at line {header_line})")
        cells += [""] * (len(header) - len(cells))
        row = dict(zip(header, cells))
        pattern = row["task"]
        if not pattern or not _NAME_RE.match(pattern):
            raise PolicyError(f"{name}:{lineno}: bad task name or glob {pattern!r}")
        values: Dict[str, object] = {}
        for col in VALUE_COLUMNS:
            cell = row.get(col, "")
            if cell in UNSET:
                continue
            try:
                values[col] = _PARSERS[col](cell)
            except ValueError as exc:
                raise PolicyError(f"{name}:{lineno}: column {col}: {exc}") from None
        rules.append(Rule(pattern=pattern, source=f"{name}:{lineno}", **values))  # type: ignore[arg-type]
    if header is None and any(l.strip() and not l.lstrip().startswith("#") for l in lines):
        raise PolicyError(f"{name}: no header line")
    return rules


def parse_policy(path: str) -> List[Rule]:
    with open(path, encoding="utf-8") as fh:
        return parse_policy_text(fh.read(), path)


_CACHE: Dict[str, Tuple[Tuple[float, int], List[Rule]]] = {}


def load_policy(path: str) -> List[Rule]:
    """parse_policy with a per-process cache keyed on mtime and size; a missing file is an empty policy."""
    try:
        st = os.stat(path)
    except FileNotFoundError:
        return []
    key = (st.st_mtime, st.st_size)
    hit = _CACHE.get(path)
    if hit is not None and hit[0] == key:
        return hit[1]
    rules = parse_policy(path)
    _CACHE[path] = (key, rules)
    return rules


# ---- resolution ------------------------------------------------------------

def matching(rules: Iterable[Rule], task_name: str) -> List[Rule]:
    return [r for r in rules if r.matches(task_name)]


def resolve(rules: Iterable[Rule], task_name: str) -> Dict[str, Tuple[object, Rule]]:
    """Merged cells for one task: policy column -> (value, the row it came from); later rows win."""
    out: Dict[str, Tuple[object, Rule]] = {}
    for rule in matching(rules, task_name):
        for col, value in rule.values().items():
            out[col] = (value, rule)
    return out


def apply(rules: Iterable[Rule], task_name: str, runtime_values: Dict[str, object],
          limits: Optional[Dict[str, int]] = None) -> Tuple[Dict[str, Tuple[object, object]], List[str]]:
    """Write the resolved cells into miniwdl's runtime_values dict (in place).

    Returns ({runtime key: (old, new)} for the keys that changed, [sources of the rows used]).
    `limits` is the backend's detect_resource_limits() answer ({"cpu", "mem_bytes"}); cpu and memory
    never exceed it (miniwdl-slurm reports sys.maxsize, the local backends the host).
    """
    resolved = resolve(rules, task_name)
    changes: Dict[str, Tuple[object, object]] = {}
    sources: List[str] = []
    for col, (value, rule) in resolved.items():
        key = RUNTIME_KEYS[col]
        new: object = value
        if col == "cpu" and limits and limits.get("cpu", 0) > 0:
            new = max(1, min(int(value), int(limits["cpu"])))  # type: ignore[arg-type]
        if col == "memory" and limits and limits.get("mem_bytes", 0) > 0:
            new = min(int(value), int(limits["mem_bytes"]))  # type: ignore[arg-type]
        old = runtime_values.get(key)
        if old == new:
            continue
        if col == "memory" and isinstance(old, int) and old > 0 and isinstance(runtime_values.get("memory_limit"), int):
            limit = runtime_values["memory_limit"]
            runtime_values["memory_limit"] = int(limit * int(new) / old)  # type: ignore[arg-type]
        runtime_values[key] = new
        changes[key] = (old, new)
        if rule.source not in sources:
            sources.append(rule.source)
    return changes, sources


def describe(changes: Dict[str, Tuple[object, object]]) -> str:
    """'cpu 64→48, memory 256G→200G, partition -→fat' for a log line."""
    parts = []
    for key, (old, new) in changes.items():
        col = next((c for c, k in RUNTIME_KEYS.items() if k == key), key)
        if col == "memory":
            fmt = lambda v: format_memory(int(v)) if isinstance(v, int) else "-"  # noqa: E731
        elif col == "time":
            fmt = lambda v: format_time(int(v)) if isinstance(v, int) else "-"  # noqa: E731
        else:
            fmt = lambda v: str(v) if v is not None else "-"  # noqa: E731
        parts.append(f"{col} {fmt(old)}→{fmt(new)}")
    return ", ".join(parts)
