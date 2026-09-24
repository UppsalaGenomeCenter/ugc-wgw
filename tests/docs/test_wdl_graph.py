"""Tests for scripts/wdl-graph.py: extraction rules on an inline WDL fixture, size budget on the real entrypoints."""
from __future__ import annotations

import importlib.util
import pathlib
import re
import sys
import tempfile
import unittest

REPO = pathlib.Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "wdl-graph.py"

SUB_WDL = """version 1.1
task t3 { input { File a } command <<< echo ~{a} >>> output { File o = stdout() } }
workflow sub {
  input {
    File x
    Array[File] xs
  }
  scatter (i in xs) { call t3 { input: a = x } }
  output { File out = t3.o[0] }
}
"""

MAIN_WDL = """version 1.1
import "sub.wdl" as S
task t1 { command <<< echo 1 >>> output { File o = stdout() } }
task t2 { input { File a } command <<< echo ~{a} >>> output { File o = stdout() } }
task t4 { input { Array[File] xs } command <<< echo ~{sep(" ", xs)} >>> output { File o = stdout() } }
workflow main {
  call t1
  scatter (i in range(2)) { call t2 { input: a = t1.o } }
  if (length(t2.o) > 1) { call t4 as merge { input: xs = t2.o } }
  File hub = select_first([merge.o, t2.o[0]])
  scatter (j in range(2)) { String only_name = "~{j}" }
  call S.sub { input: x = hub, xs = t2.o }
  call t4 as after_one after t1 { input: xs = [hub] }
  String passthru = basename(sub.out)
  call t4 as last { input: xs = [passthru] }
  if (defined(hub) && basename(hub) == "q\\"x" && length(t2.o) > 0) { call t1 as quoted }
  output { File h = hub }
}
"""


def setUpModule():
    try:
        import WDL  # noqa: F401
    except ImportError:
        raise unittest.SkipTest("miniwdl not importable (activate .venv)")


def load_script():
    spec = importlib.util.spec_from_file_location("wdl_graph", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules["wdl_graph"] = mod  # dataclasses resolve annotations through sys.modules
    spec.loader.exec_module(mod)
    return mod


NODE_RE = re.compile(r'^\s*(n\d+)(?:\[\[|\(\[|\[)"([^"]*)"')
SECTION_RE = re.compile(r'^\s*subgraph (s\d+)\["([^"]*)"\]')
EDGE_RE = re.compile(r'^\s*([ns]\d+) (-->|-\.->) ([ns]\d+)$')


def parse(text: str):
    """label -> id for nodes (last wins), title -> id for sections, the edge set, and every declared id."""
    nodes, sections, edges, ids = {}, {}, set(), []
    for line in text.splitlines():
        m = NODE_RE.match(line)
        if m:
            nodes[m.group(2)] = m.group(1)
            ids.append(m.group(1))
            continue
        m = SECTION_RE.match(line)
        if m:
            sections[m.group(2)] = m.group(1)
            ids.append(m.group(1))
            continue
        m = EDGE_RE.match(line)
        if m:
            edges.add((m.group(1), m.group(2), m.group(3)))
    return nodes, sections, edges, ids


def unescape(label: str) -> str:
    for ent, ch in (("#quot;", '"'), ("#lt;", "<"), ("#gt;", ">"), ("#38;", "&"), ("#35;", "#")):
        label = label.replace(ent, ch)
    return label


class FixtureTest(unittest.TestCase):
    def setUp(self):
        self.mod = load_script()
        self.tmp = tempfile.TemporaryDirectory()
        d = pathlib.Path(self.tmp.name)
        (d / "sub.wdl").write_text(SUB_WDL)
        self.main = d / "main.wdl"
        self.main.write_text(MAIN_WDL)

    def tearDown(self):
        self.tmp.cleanup()

    def render(self, depth, hide=frozenset()):
        text, graph = self.mod.graph_for(self.main, depth, hide)
        return text, graph, parse(text)[:3]

    def test_depth0_nodes_sections_edges(self):
        text, graph, (nodes, sections, edges) = self.render(0)
        for label in ("t1", "t2", "merge", "hub", "sub", "after_one", "last", "quoted"):
            self.assertIn(label, nodes, text)
        self.assertIn(f'{nodes["sub"]}[["sub"]]', text)  # workflow call drawn as one node
        self.assertIn(f'{nodes["hub"]}(["hub"])', text)  # select_first over two calls is a hub
        self.assertIn("scatter i in range(2)", sections)
        self.assertIn("if length(t2.o) #gt; 1", sections)
        self.assertFalse(any(t.startswith("scatter j") for t in sections), "section without calls must be elided")
        n, s = nodes, sections
        expected = {
            (n["t1"], "-->", n["t2"]),
            (n["t2"], "-.->", s["if length(t2.o) #gt; 1"]),
            (n["t2"], "-->", n["merge"]),
            (n["t2"], "-->", n["hub"]),
            (n["merge"], "-->", n["hub"]),
            (n["hub"], "-->", n["sub"]),
            (n["hub"], "-->", n["after_one"]),
            (n["t1"], "-->", n["after_one"]),  # `after t1`
            (n["sub"], "-->", n["last"]),  # through the pass-through decl `passthru`
        }
        self.assertTrue(expected <= edges, f"missing: {expected - edges}\n{text}")
        quoted_if = [t for t in sections if t.startswith("if defined(hub)")]
        self.assertEqual(len(quoted_if), 1, sections)
        self.assertIn((n["hub"], "-.->", s[quoted_if[0]]), edges)
        self.assertEqual(graph.linked, [("sub", str(pathlib.Path(self.tmp.name) / "sub.wdl"))])

    def test_escape_and_truncate(self):
        _, _, (_, sections, _) = self.render(0)
        title = next(t for t in sections if t.startswith("if defined(hub)"))
        self.assertIn("#quot;q\\#quot;x", title)
        self.assertIn("#38;#38;", title)
        self.assertTrue(title.endswith("…"), title)
        self.assertLessEqual(len(unescape(title)), len("if ") + self.mod.EXPR_CHARS)
        self.assertNotIn('"', title)

    def test_depth1_inlines_subworkflow(self):
        text, graph, (nodes, sections, edges) = self.render(1)
        self.assertIn("sub (sub.wdl)", sections)
        self.assertNotIn("sub", nodes)
        self.assertIn("t3", nodes)
        self.assertIn("scatter i in xs", sections)
        self.assertIn((nodes["hub"], "-->", nodes["t3"]), edges)  # input x maps to t3's consumer
        self.assertIn((nodes["t2"], "-.->", sections["scatter i in xs"]), edges)  # input xs controls the scatter
        self.assertIn((sections["sub (sub.wdl)"], "-->", nodes["last"]), edges)
        self.assertEqual(graph.inlined, [("sub", str(pathlib.Path(self.tmp.name) / "sub.wdl"))])
        self.assertEqual(graph.linked, [])

    def test_hide(self):
        text, _, (nodes, _, edges) = self.render(0, frozenset({"t1"}))
        self.assertNotIn("t1", nodes)
        self.assertIn("quoted", nodes)  # hidden by call name, not task name
        self.assertFalse(any(dst == nodes["t2"] for _, _, dst in edges))

    def test_size_budget(self):
        with self.assertRaises(self.mod.SizeError):
            self.mod.graph_for(self.main, 1, frozenset(), max_chars=50)


class EntrypointsTest(unittest.TestCase):
    def setUp(self):
        self.mod = load_script()

    def test_all_entrypoints_render_within_budget(self):
        for ep in self.mod.ENTRYPOINTS:
            text, _ = self.mod.graph_for(REPO / "workflows" / f"{ep}.wdl")
            self.assertLessEqual(len(text), self.mod.MAX_CHARS, ep)
            nodes, sections, edges, ids = parse(text)
            declared = set(ids)
            self.assertEqual(len(ids), len(declared), f"{ep}: duplicate ids")
            for src, _, dst in edges:
                self.assertIn(src, declared, f"{ep}: {src}")
                self.assertIn(dst, declared, f"{ep}: {dst}")
            node_ids = [i for i in ids if i.startswith("n")]
            self.assertEqual(node_ids, [f"n{i}" for i in range(1, len(node_ids) + 1)], ep)
            self.assertIn("ugc_wgw_manifest_write", nodes, ep)
            self.assertNotIn("backend_configuration", nodes, ep)

    def test_markdown_is_deterministic(self):
        a = self.mod.markdown()
        b = self.mod.markdown()
        self.assertEqual(a, b)
        self.assertEqual(a.count("```mermaid"), len(self.mod.ENTRYPOINTS))
        self.assertTrue(a.startswith("<!-- Generated by scripts/wdl-graph.py"))


if __name__ == "__main__":
    unittest.main()
