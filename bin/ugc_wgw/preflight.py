"""The input check a run gets just before it starts (docs/DESIGN.md §8.1).

Registration (`ugc-wgw samples add`) reads every BAM once; this is the cheap repeat for the moment a
run is created: the raw read files named in the generated inputs must still be there, non-empty, the
size registration recorded, and end with the BGZF marker. A failure becomes a failed attempt of kind
`input` (error class `InputError`) without a SLURM job, blocked until `ugc-wgw retry`.
"""
from __future__ import annotations

from pathlib import Path

from . import bam
from .db import DB

READ_KEYS = ("hifi_reads", "fail_reads", "father_hifi_reads", "mother_hifi_reads")


def read_paths(doc: dict[str, object]) -> list[str]:
    """The raw read files an inputs document names (keys may be namespaced `ugc_wgw_<stage>.hifi_reads`)."""
    out: list[str] = []
    for key, value in doc.items():
        if key.rsplit(".", 1)[-1] in READ_KEYS and isinstance(value, list):
            out.extend(str(v) for v in value)
    return out


def check(doc: dict[str, object], db: DB) -> tuple[int, int, list[str]]:
    """(files checked, their bytes, problems) for the raw read files of an inputs document."""
    paths = read_paths(doc)
    sizes = db.input_sizes(paths) if paths else {}
    problems: list[str] = []
    total = 0
    for p in paths:
        problem = bam.quick_check(Path(p), sizes.get(p))
        if problem:
            problems.append(f"{p}: {problem}")
        else:
            total += Path(p).stat().st_size
    return len(paths), total, problems
