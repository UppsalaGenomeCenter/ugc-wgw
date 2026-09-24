#!/usr/bin/env python3
"""Inventory of the resource requests of every task reachable from the entrypoints.

Two generated files, never edited by hand (tests/check.sh fails when they are stale):

  backends/hpc/resources.declared.tsv   one row per task: cpu, memory, informational columns;
                                        the same columns as the site policy, so rows can be
                                        copied into <prefix>/resources.tsv and edited
  docs/guide/task-resources.md          the same data as tables, one per entrypoint

Values are folded statically from the WDL: literal declarations and arithmetic on them,
task input defaults, and call-site values that fold from the calling workflow's own literal
declarations (DeepVariant's `tasks_per_shard`). Anything else (input sizes via size(),
workflow inputs, `n_samples`) is left as the expression in the notes column. A note also says
when the task's command interpolates the same declaration (`--threads ~{threads}`): lowering
the request then oversubscribes, lowering memory below `--mem-gbytes` can kill the tool.

Usage:
  python3 scripts/wdl-resources.py --write     # regenerate both files
  python3 scripts/wdl-resources.py --check     # exit 1 with a diff when either is stale
  python3 scripts/wdl-resources.py             # the TSV on stdout (--format md for the Markdown)

Dev-time only: needs miniwdl importable (activate .venv). See docs/guide/12-resources.md.
"""
from __future__ import annotations

import argparse
import difflib
import importlib.util
import math
import pathlib
import sys
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Set, Tuple

try:
    import WDL
    from WDL import Env, Expr, StdLib, Tree, Type, Value
    from WDL._util import parse_byte_size
except ImportError:  # pragma: no cover
    sys.exit("scripts/wdl-resources.py needs miniwdl (dev-time only): activate .venv or pip install miniwdl")

REPO = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "plugins" / "ugc_wgw_miniwdl"))
from ugc_wgw_miniwdl import policy  # noqa: E402  (format_memory; the policy columns)


def _load_graph_script():
    spec = importlib.util.spec_from_file_location("wdl_graph", REPO / "scripts" / "wdl-graph.py")
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules.setdefault("wdl_graph", mod)
    spec.loader.exec_module(mod)
    return mod


graph = _load_graph_script()
graph.EXPR_CHARS = 120  # the diagrams truncate expressions early; the notes need the whole formula
ENTRYPOINTS: List[str] = graph.ENTRYPOINTS
TSV_OUTPUT = REPO / "backends" / "hpc" / "resources.declared.tsv"
MD_OUTPUT = REPO / "docs" / "guide" / "task-resources.md"
DYNAMIC_FUNCTIONS = {"size", "glob", "stdout", "stderr"}
TSV_COLUMNS = policy.COLUMNS + policy.INFO_COLUMNS
LARGE_CPU = 16
LARGE_MEMORY = 64 * 1024 ** 3


class Dynamic(Exception):
    """The expression depends on something not known statically."""


@dataclass
class Folded:
    value: Optional[Value.Base]
    notes: List[str] = field(default_factory=list)


@dataclass
class TaskInfo:
    name: str
    task: Tree.Task
    doc: Tree.Document
    stages: List[str] = field(default_factory=list)
    call_sites: List[Dict[str, Optional[Value.Base]]] = field(default_factory=list)
    cpu: Optional[int] = None
    cpu_text: str = ""
    memory: Optional[int] = None
    memory_text: str = ""
    notes: List[str] = field(default_factory=list)

    @property
    def source(self) -> str:
        return f"{graph.relpath(self.task.pos.abspath)}:{self.task.pos.line}"


# ---- static folding ----------------------------------------------------------

def idents(expr: Expr.Base) -> Iterable[Expr.Ident]:
    if isinstance(expr, Expr.Ident):
        yield expr
    for ch in expr.children:
        if isinstance(ch, Expr.Base):
            yield from idents(ch)


def applies(expr: Expr.Base) -> Iterable[Expr.Apply]:
    if isinstance(expr, Expr.Apply):
        yield expr
    for ch in expr.children:
        if isinstance(ch, Expr.Base):
            yield from applies(ch)


def dynamic_call(expr: Expr.Base) -> bool:
    return any(a.function_name in DYNAMIC_FUNCTIONS or a.function_name.startswith("read_") for a in applies(expr))


def fold(expr: Expr.Base, static: Dict[int, Value.Base], std: StdLib.Base) -> Value.Base:
    """Evaluate `expr` from the static values of the declarations it references (by identity)."""
    if dynamic_call(expr):
        raise Dynamic()
    env: Env.Bindings[Value.Base] = Env.Bindings()
    for ident in idents(expr):
        ref = getattr(ident, "referee", None)
        if not isinstance(ref, Tree.Decl) or id(ref) not in static:
            raise Dynamic()
        env = env.bind(str(ident.name), static[id(ref)])
    try:
        return expr.eval(env, std)
    except Exception:
        raise Dynamic()


def workflow_statics(wf: Tree.Workflow, std: StdLib.Base) -> Dict[int, Value.Base]:
    """Literal-foldable declarations of a workflow body, and the defaults of its inputs."""
    static: Dict[int, Value.Base] = {}
    for decl in wf.inputs or []:
        if decl.expr is not None:
            try:
                static[id(decl)] = fold(decl.expr, static, std)
            except Dynamic:
                pass

    def visit(body: Iterable[Tree.WorkflowNode]) -> None:
        for node in body:
            if isinstance(node, Tree.Decl) and node.expr is not None:
                try:
                    static[id(node)] = fold(node.expr, static, std)
                except Dynamic:
                    pass
            elif isinstance(node, Tree.WorkflowSection):
                visit(node.body)

    visit(wf.body)
    return static


# ---- reachability --------------------------------------------------------------

def site_text(expr: Expr.Base, doc: Tree.Document, inherited: Dict[str, str]) -> str:
    """A call input as written; when it is just a name, also what that name is: the workflow declaration's
    expression, or (for a subworkflow input) what the calling workflow passed."""
    text = graph.expr_text(expr, doc)
    if isinstance(expr, Expr.Get) and expr.member is None and isinstance(expr.expr, Expr.Ident):
        ref = getattr(expr.expr, "referee", None)
        if isinstance(ref, Tree.Decl):
            if ref.expr is not None:
                text += f" = {graph.expr_text(ref.expr, doc)}"
            elif ref.name in inherited:
                text += f" = {inherited[ref.name]}"
    return text


def site_key(entry: Dict[str, object]) -> str:
    """Identity of a call site for de-duplication (a subworkflow called twice with the same values)."""
    parts = []
    for k in sorted(entry):
        v = entry[k]
        parts.append(f"{k}={v.json if isinstance(v, Value.Base) else v}")
    return "|".join(str(p) for p in parts)


def wdl_version(doc: Tree.Document) -> str:
    return getattr(doc, "effective_wdl_version", None) or getattr(doc, "wdl_version", None) or "1.0"

def walk(doc: Tree.Document, stage: str, tasks: Dict[str, TaskInfo], docs: Dict[str, Tree.Document],
         inherited: Optional[Dict[str, str]] = None) -> None:
    """Record every task call under this workflow (recursing into called workflows) for one entrypoint."""
    wf = doc.workflow
    assert wf is not None
    inherited = {} if inherited is None else inherited
    std = StdLib.Base(wdl_version(doc))
    static = workflow_statics(wf, std)
    wf_inputs = {id(d) for d in (wf.inputs or [])}

    def visit(body: Iterable[Tree.WorkflowNode]) -> None:
        for node in body:
            if isinstance(node, Tree.Call):
                callee = node.callee
                if isinstance(callee, Tree.Task):
                    site: Dict[str, Optional[Value.Base]] = {}
                    notes: Dict[str, str] = {}
                    texts: Dict[str, str] = {}
                    for name, expr in node.inputs.items():
                        texts[name] = site_text(expr, doc, inherited)
                        try:
                            site[name] = fold(expr, static, std)
                            refs = [getattr(i, "referee", None) for i in idents(expr)]
                            notes[name] = ("workflow input default" if any(id(r) in wf_inputs for r in refs)
                                           else "call site")
                        except Dynamic:
                            site[name] = None
                    key = callee.name
                    info = tasks.get(key)
                    if info is None:
                        info = TaskInfo(key, callee, docs[callee.pos.abspath])
                        tasks[key] = info
                    elif (info.task.pos.abspath, info.task.pos.line) != (callee.pos.abspath, callee.pos.line):
                        raise SystemExit(f"task name {key} is defined twice: {info.source} and "
                                         f"{graph.relpath(callee.pos.abspath)}:{callee.pos.line}; the policy keys on "
                                         "bare task names, rename one")
                    if stage not in info.stages:
                        info.stages.append(stage)
                    entry = {"__notes__": notes, "__texts__": texts, **site}
                    if site_key(entry) not in {site_key(e) for e in info.call_sites}:
                        info.call_sites.append(entry)  # type: ignore[arg-type]
                elif isinstance(callee, Tree.Workflow):
                    passed = {name: site_text(expr, doc, inherited) for name, expr in node.inputs.items()}
                    walk(docs[callee.pos.abspath], stage, tasks, docs, passed)
            elif isinstance(node, Tree.WorkflowSection):
                visit(node.body)

    visit(wf.body)


# ---- per-task evaluation ----------------------------------------------------

def decl_names_behind(expr: Optional[Expr.Base], by_name: Dict[str, Tree.Decl], seen: Optional[Set[str]] = None) -> Set[str]:
    """Names of the task declarations an expression depends on, transitively through their own expressions."""
    seen = set() if seen is None else seen
    if expr is None:
        return seen
    for ident in idents(expr):
        name = str(ident.name)
        if name in by_name and name not in seen:
            seen.add(name)
            decl_names_behind(by_name[name].expr, by_name, seen)
    return seen


def evaluate(info: TaskInfo) -> None:
    task = info.task
    std = StdLib.Base(task.effective_wdl_version)
    inputs = list(task.inputs or [])
    by_name = {d.name: d for d in inputs + list(task.postinputs)}
    results: List[Tuple[Optional[int], Optional[int], List[str]]] = []
    for site in info.call_sites:
        site_notes: Dict[str, str] = site.get("__notes__", {})  # type: ignore[assignment]
        static: Dict[int, Value.Base] = {}
        notes: List[str] = []
        for decl in inputs:
            given = site.get(decl.name)
            if given is not None:
                static[id(decl)] = given
            elif decl.name not in site and decl.expr is not None:
                try:
                    static[id(decl)] = fold(decl.expr, static, std)
                    site_notes[decl.name] = "input default"
                except Dynamic:
                    pass
        for decl in task.postinputs:
            if decl.expr is None:
                continue
            try:
                static[id(decl)] = fold(decl.expr, static, std)
            except Dynamic:
                pass
        cpu = memory = None
        for key in ("cpu", "memory"):
            expr = task.runtime.get(key)
            if expr is None:
                continue
            try:
                v = fold(expr, static, std)
            except Dynamic:
                continue
            if key == "cpu":
                cpu = max(1, math.ceil(v.coerce(Type.Float()).value))
            else:
                memory = parse_byte_size(v.coerce(Type.String()).value)
            used = decl_names_behind(expr, by_name)
            origins = sorted({site_notes[n] for n in used if n in site_notes})
            for o in origins:
                note = f"{key} from {o}"
                if note not in notes:
                    notes.append(note)
        results.append((cpu, memory, notes))
    cpus = {r[0] for r in results}
    mems = {r[1] for r in results}
    info.notes = []
    if len(cpus) == 1:
        info.cpu = cpus.pop()
    else:
        info.notes.append("cpu varies by call site: " + ", ".join(str(r[0] or "?") for r in results))
    if len(mems) == 1:
        info.memory = mems.pop()
    else:
        info.notes.append("memory varies by call site: " + ", ".join(policy.format_memory(r[1]) if r[1] else "?" for r in results))
    for _, _, notes in results:
        for n in notes:
            if n not in info.notes:
                info.notes.append(n)
    input_names = {d.name for d in inputs}
    for key in ("cpu", "memory"):
        expr = task.runtime.get(key)
        text = graph.expr_text(expr, info.doc) if expr is not None else "-"
        setattr(info, f"{key}_text", text)
        if expr is None or (key == "cpu" and info.cpu is not None) or (key == "memory" and info.memory is not None):
            continue
        used = decl_names_behind(expr, by_name)
        texts = [text]
        shown: List[str] = []
        for n in sorted(used):
            if n in input_names:
                for site in info.call_sites:
                    t = site.get("__texts__", {}).get(n)  # type: ignore[union-attr]
                    if t and f"{n} = {t}" not in shown:
                        shown.append(f"{n} = {t}")
            elif by_name[n].expr is not None:
                shown.append(f"{n} = {graph.expr_text(by_name[n].expr, info.doc)}")
        texts += shown
        info.notes.append(f"{key}: {text}" + (f" where {'; '.join(shown)}" if shown else ""))
        if any("n_samples" in t for t in texts):
            info.notes.append(f"{key} scales with N")
    # the command's closure (a `sort -S ~{sort_mem}` derived from mem_gb counts) against the declarations
    # the runtime keys name directly (mem_gb = threads * 6 does not make the command "use mem")
    cmd_refs = decl_names_behind(task.command, by_name)
    cpu_names = {str(i.name) for e in [task.runtime.get("cpu")] if e is not None for i in idents(e)}
    mem_names = {str(i.name) for e in [task.runtime.get("memory")] if e is not None for i in idents(e)}
    if cmd_refs & cpu_names:
        info.notes.append("command uses threads")
    if cmd_refs & mem_names:
        info.notes.append("command uses mem")
    if "gpuCount" in task.runtime or "gpu" in task.runtime:
        info.notes.append("gpu")


def inventory() -> Dict[str, TaskInfo]:
    tasks: Dict[str, TaskInfo] = {}
    for ep in ENTRYPOINTS:
        doc = WDL.load(str(REPO / "workflows" / f"{ep}.wdl"))
        docs = graph.collect_docs(doc)
        walk(doc, ep, tasks, docs)
    for info in tasks.values():
        evaluate(info)
    return tasks


# ---- rendering ---------------------------------------------------------------

def stage_label(ep: str) -> str:
    return ep[len("ugc_wgw_"):] if ep.startswith("ugc_wgw_") else ep


def tsv_text(tasks: Dict[str, TaskInfo]) -> str:
    out = [
        "# Generated by scripts/wdl-resources.py from workflows/ugc_wgw_*.wdl; do not edit.",
        "# Regenerate: python3 scripts/wdl-resources.py --write (tests/check.sh fails on drift).",
        "# Declared cpu and memory of every task the entrypoints can call. Copy rows into the site",
        "# policy (<prefix>/resources.tsv) and edit the values; '-' in a cell means 'as declared'.",
        "# stages, source and notes are informational and ignored by the policy.",
        "\t".join(TSV_COLUMNS),
    ]
    for name in sorted(tasks):
        t = tasks[name]
        out.append("\t".join([
            name,
            str(t.cpu) if t.cpu is not None else "-",
            policy.format_memory(t.memory) if t.memory is not None else "-",
            "-", "-", "-",
            ",".join(stage_label(s) for s in t.stages),
            t.source,
            "; ".join(t.notes),
        ]))
    return "\n".join(out) + "\n"


def md_cell(s: str) -> str:
    return s.replace("|", "\\|")


def md_text(tasks: Dict[str, TaskInfo]) -> str:
    out = [
        "<!-- Generated by scripts/wdl-resources.py from workflows/ugc_wgw_*.wdl; do not edit.",
        "     Regenerate: python3 scripts/wdl-resources.py --write (tests/check.sh fails on drift). -->",
        "",
        "# Task resources",
        "",
        "What every task the entrypoints can call asks SLURM for, folded from the WDL",
        "by `scripts/wdl-resources.py`. How to cap or change these requests without",
        "editing WDL is chapter 12; the same rows in policy-file form are",
        "`backends/hpc/resources.declared.tsv`.",
        "",
        "Reading the tables:",
        "",
        "- `cpu` and `memory` are the declared requests (`runtime.cpu`, `runtime.memory`).",
        "  `-` means the value depends on the inputs; the note shows the expression.",
        "- `command uses threads` (`mem`): the command interpolates the same",
        "  declaration, so lowering the request oversubscribes the cores (slower, still",
        "  correct) or lowers the tool's own memory budget (may be killed).",
        "- `input default`, `workflow input default`: the value is an input default,",
        "  overridable through `stage_inputs` (chapter 05); `call site`: fixed by the",
        "  calling workflow.",
        "- `source` is the task definition; `stages` lists the entrypoints that call it.",
        "",
    ]
    large = sorted((t for t in tasks.values() if (t.cpu or 0) >= LARGE_CPU or (t.memory or 0) >= LARGE_MEMORY),
                   key=lambda t: (-(t.cpu or 0), -(t.memory or 0), t.name))
    out += [f"## Largest requests (cpu ≥ {LARGE_CPU} or memory ≥ {policy.format_memory(LARGE_MEMORY)})", "",
            "| task | cpu | memory | stages | notes |", "|---|---|---|---|---|"]
    for t in large:
        out.append(f"| `{t.name}` | {t.cpu if t.cpu is not None else '-'} | "
                   f"{policy.format_memory(t.memory) if t.memory is not None else '-'} | "
                   f"{', '.join(stage_label(s) for s in t.stages)} | {md_cell('; '.join(t.notes)) or '-'} |")
    out.append("")
    seen_stages = {s for t in tasks.values() for s in t.stages}
    stages = [ep for ep in ENTRYPOINTS if ep in seen_stages] + sorted(seen_stages - set(ENTRYPOINTS))
    for ep in stages:
        rows = [t for t in tasks.values() if ep in t.stages]
        out += [f"## {ep}", "", f"{len(rows)} tasks, in first-call order.", "",
                "| task | cpu | memory | notes | source |", "|---|---|---|---|---|"]
        for t in rows:
            out.append(f"| `{t.name}` | {t.cpu if t.cpu is not None else '-'} | "
                       f"{policy.format_memory(t.memory) if t.memory is not None else '-'} | "
                       f"{md_cell('; '.join(t.notes)) or '-'} | `{t.source}` |")
        out.append("")
    return "\n".join(out)


def check(text: str, target: pathlib.Path) -> int:
    rel = graph.relpath(str(target))
    try:
        current = target.read_text()
    except FileNotFoundError:
        print(f"{rel}: missing", file=sys.stderr)
        return 1
    if current == text:
        return 0
    sys.stdout.writelines(difflib.unified_diff(current.splitlines(True), text.splitlines(True),
                                               f"{rel} (committed)", "generated"))
    print(f"{rel} is stale; run: python3 scripts/wdl-resources.py --write", file=sys.stderr)
    return 1


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--write", action="store_true", help="write both generated files")
    ap.add_argument("--check", action="store_true", help="exit 1 with a diff when a committed file is stale")
    ap.add_argument("--format", choices=["tsv", "md"], default="tsv", help="what to print on stdout")
    a = ap.parse_args(argv)
    tasks = inventory()
    tsv, md = tsv_text(tasks), md_text(tasks)
    if a.check:
        return max(check(tsv, TSV_OUTPUT), check(md, MD_OUTPUT))
    if a.write:
        for text, target in ((tsv, TSV_OUTPUT), (md, MD_OUTPUT)):
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text(text)
            print(f"wrote {graph.relpath(str(target))}", file=sys.stderr)
        return 0
    sys.stdout.write(tsv if a.format == "tsv" else md)
    return 0


if __name__ == "__main__":
    sys.exit(main())
