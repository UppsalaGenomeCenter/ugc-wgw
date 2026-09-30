#!/usr/bin/env python3
"""Engine-side half of `ugc-wgw stage-inputs --nested`: every input miniwdl accepts for a workflow, the
call-qualified ones included (`upstream.pbmm2.pbmm2_align_wgs.threads`), as JSON.

Runs under the engine venv's Python because it imports miniwdl's `WDL` package; the driver stays
stdlib-only and only reads the output (bin/ugc_wgw/stage_inputs.py). Nothing is executed: the documents are
loaded and their `available_inputs` walked, the same list `miniwdl run <wdl>` prints without inputs.

usage: python wdl_available_inputs.py <wdl>...
stdout: {"<wdl>": [{"name", "type", "default", "required", "description"}, ...], ...}
"""
from __future__ import annotations

import json
import sys

try:
    import WDL
except ImportError:  # pragma: no cover - the driver reports it
    sys.stderr.write("wdl_available_inputs: this interpreter has no WDL package (miniwdl); run it with the engine venv's python\n")
    sys.exit(3)


def find_call(node, name):
    """The call named `name` in a workflow body, through scatter and conditional sections."""
    for el in getattr(node, "body", []) or []:
        if isinstance(el, WDL.Tree.Call) and el.name == name:
            return el
        if isinstance(el, WDL.Tree.WorkflowSection):
            found = find_call(el, name)
            if found is not None:
                return found
    return None


def describe(meta) -> str:
    """A parameter_meta entry as one line: a bare string, or an object's description plus its choices."""
    if meta is None:
        return ""
    if isinstance(meta, str):
        return meta
    if isinstance(meta, dict):
        desc = str(meta.get("description") or "")
        choices = meta.get("choices")
        if isinstance(choices, list) and choices:
            joined = ", ".join(str(c) for c in choices)
            return f"{desc} (choices: {joined})" if desc else f"choices: {joined}"
        return desc
    return str(meta)


def description(wf, name: str) -> str:
    """parameter_meta of the declaration behind a (possibly call-qualified) input name."""
    parts = name.split(".")
    node = wf
    for part in parts[:-1]:
        call = find_call(node, part)
        if call is None:
            return ""
        node = call.callee
    meta = getattr(node, "parameter_meta", None) or {}
    return describe(meta.get(parts[-1]))


def available(path: str) -> list:
    doc = WDL.load(path)
    wf = doc.workflow
    if wf is None:
        raise SystemExit(f"{path}: no workflow")
    required = {b.name for b in wf.required_inputs}
    rows = []
    for b in wf.available_inputs:
        decl = b.value
        rows.append({
            "name": b.name,
            "type": str(decl.type),
            "default": str(decl.expr) if decl.expr is not None else None,
            "required": b.name in required,
            "description": description(wf, b.name),
        })
    return rows


def main(argv: list) -> int:
    if not argv:
        sys.stderr.write(__doc__)
        return 2
    json.dump({path: available(path) for path in argv}, sys.stdout)
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
