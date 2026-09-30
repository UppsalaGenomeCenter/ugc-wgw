"""The stage table: which entrypoint each stage runs, its subject type, and the mode sequences (docs/DESIGN.md §5)."""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

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


@dataclass(frozen=True)
class InputSpec:
    """One declaration of an entrypoint's `input {` block, as written."""
    name: str
    type: str                    # `Array[File]?`, `Int`, ...
    default: Optional[str]       # the default expression's text, None when there is none
    required: bool               # no default and not optional (`?`)


def _block_after(text: str, opener: str, start: int) -> tuple[int, int]:
    """(body start, body end) of the first `opener {` at or after `start`; (-1, -1) when absent."""
    at = text.find(opener, start)
    if at < 0:
        return -1, -1
    i = text.find("{", at)
    depth, body_start = 0, i + 1
    while i < len(text):
        ch = text[i]
        if ch == "{":
            depth += 1
        elif ch == "}":
            depth -= 1
            if depth == 0:
                return body_start, i
        i += 1
    return body_start, len(text)


def _workflow_start(wdl: Path, text: str) -> int:
    m = re.search(r"^\s*workflow\s+\w+\s*\{", text, re.M)
    if not m:
        raise UgcError(f"{wdl}: no workflow block")
    return m.end()


def declared_input_specs(wdl: Path) -> list[InputSpec]:
    """Parse the workflow's `input {` block, in declaration order."""
    text = Path(wdl).read_text()
    start, end = _block_after(text, "input {", _workflow_start(wdl, text))
    if start < 0:
        raise UgcError(f"{wdl}: no input block")
    out: list[InputSpec] = []
    for line in text[start:end].splitlines():
        code = line.split("#", 1)[0].rstrip()
        if not code.strip():
            continue
        dm = _DECL_RE.match(code)
        if not dm:
            continue
        typ, name, default = dm.group(1), dm.group(2), dm.group(3)
        expr = default.lstrip("=").strip() if default is not None else None
        out.append(InputSpec(name, typ, expr, default is None and not typ.endswith("?")))
    if not out:
        raise UgcError(f"{wdl}: input block is empty")
    return out


def declared_inputs(wdl: Path) -> dict[str, bool]:
    """{name: required} of the workflow's input block. Required = no default and not optional (`?`)."""
    return {spec.name: spec.required for spec in declared_input_specs(wdl)}


_META_KEY_RE = re.compile(r"([A-Za-z_]\w*)\s*:\s*")
_STRING_RE = re.compile(r'"((?:[^"\\]|\\.)*)"')


def input_descriptions(wdl: Path) -> dict[str, str]:
    """The workflow's `parameter_meta`: {input: description}, with the `choices` list appended when there is one.
    An entry may be a bare string or an object with `description` (and `choices`), as upstream writes them."""
    text = Path(wdl).read_text()
    start, end = _block_after(text, "parameter_meta {", _workflow_start(wdl, text))
    if start < 0:
        return {}
    block, out, i = text[start:end], {}, 0
    while True:
        km = _META_KEY_RE.search(block, i)
        if not km:
            break
        name, i = km.group(1), km.end()
        if i < len(block) and block[i] == '"':
            sm = _STRING_RE.match(block, i)
            if not sm:
                break
            out[name], i = sm.group(1), sm.end()
            continue
        if i < len(block) and block[i] == "{":
            body_start, body_end = _block_after(block, "{", i)
            body = block[body_start:body_end]
            dm = re.search(r"description\s*:\s*" + _STRING_RE.pattern, body)
            cm = re.search(r"choices\s*:\s*\[([^\]]*)\]", body)
            desc = dm.group(1) if dm else ""
            if cm:
                choices = ", ".join(_STRING_RE.findall(cm.group(1)))
                desc = f"{desc} (choices: {choices})" if desc else f"choices: {choices}"
            out[name], i = desc, body_end + 1
            continue
        i += 1
    return out
