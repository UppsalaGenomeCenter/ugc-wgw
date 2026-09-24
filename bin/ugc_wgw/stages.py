"""The stage table: which entrypoint each stage runs, its subject type, and the mode sequences (docs/DESIGN.md §5)."""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

from .util import UgcError


@dataclass(frozen=True)
class StageSpec:
    name: str
    subject_type: str  # "sample" | "cohort"
    wdl: str           # file under <code>/workflows/
    builder: str       # function name in inputs.py


STAGES: dict[str, StageSpec] = {
    "singleton": StageSpec("singleton", "sample", "ugc_wgw_singleton.wdl", "build_singleton"),
    "upstream": StageSpec("upstream", "sample", "ugc_wgw_upstream.wdl", "build_upstream"),
    "cohort_call": StageSpec("cohort_call", "cohort", "ugc_wgw_cohort_call.wdl", "build_cohort_call"),
    "downstream": StageSpec("downstream", "sample", "ugc_wgw_downstream.wdl", "build_downstream"),
    "cohort_merge": StageSpec("cohort_merge", "cohort", "ugc_wgw_cohort_merge.wdl", "build_cohort_merge"),
    "cohort_freq": StageSpec("cohort_freq", "cohort", "ugc_wgw_cohort_freq.wdl", "build_cohort_freq"),
    "assembly": StageSpec("assembly", "sample", "ugc_wgw_assembly.wdl", "build_assembly"),
}

MODES: dict[str, tuple[str, ...]] = {
    "standalone": ("singleton", "cohort_merge", "cohort_freq"),
    "joint": ("upstream", "cohort_call", "downstream", "cohort_merge", "cohort_freq"),
    "assembly": ("assembly",),
}


def stage_spec(name: str) -> StageSpec:
    try:
        return STAGES[name]
    except KeyError:
        raise UgcError(f"unknown stage {name!r}; known: {', '.join(STAGES)}") from None


def mode_stages(mode: str) -> tuple[str, ...]:
    try:
        return MODES[mode]
    except KeyError:
        raise UgcError(f"unknown mode {mode!r}; known: {', '.join(MODES)}") from None


def namespace(stage: str) -> str:
    return f"ugc_wgw_{stage}"


def wdl_path(code_dir: Path, stage: str) -> Path:
    return Path(code_dir) / "workflows" / stage_spec(stage).wdl


_DECL_RE = re.compile(r"^\s*([A-Za-z]+(?:\[[^\]]+\])?\??)\s+([A-Za-z_]\w*)\s*(=.*)?$")
_CALL_RE = re.compile(r"^\s*call\s+([A-Za-z_][\w.]*)(?:\s+as\s+([A-Za-z_]\w*))?", re.M)


def declared_calls(wdl: Path) -> set[str]:
    """Names of the calls in the entrypoint (the alias, else the last segment of the callee): the first
    segment of a nested call input such as `upstream.parabricks_deepvariant.run_parabricks_deepvariant.gpuCount`."""
    text = "\n".join(line.split("#", 1)[0] for line in Path(wdl).read_text().splitlines())
    return {alias or callee.rsplit(".", 1)[-1] for callee, alias in _CALL_RE.findall(text)}


def declared_inputs(wdl: Path) -> dict[str, bool]:
    """Parse the workflow's `input {` block: {name: required}. Required = no default and not optional (`?`)."""
    text = Path(wdl).read_text()
    m = re.search(r"^\s*workflow\s+\w+\s*\{", text, re.M)
    if not m:
        raise UgcError(f"{wdl}: no workflow block")
    start = text.find("input {", m.end())
    if start < 0:
        raise UgcError(f"{wdl}: no input block")
    depth = 0
    i = text.find("{", start)
    body_start = i + 1
    while i < len(text):
        ch = text[i]
        if ch == "{":
            depth += 1
        elif ch == "}":
            depth -= 1
            if depth == 0:
                break
        i += 1
    block = text[body_start:i]
    out: dict[str, bool] = {}
    for line in block.splitlines():
        code = line.split("#", 1)[0].rstrip()
        if not code.strip():
            continue
        dm = _DECL_RE.match(code)
        if not dm:
            continue
        typ, name, default = dm.group(1), dm.group(2), dm.group(3)
        out[name] = default is None and not typ.endswith("?")
    if not out:
        raise UgcError(f"{wdl}: input block is empty")
    return out
