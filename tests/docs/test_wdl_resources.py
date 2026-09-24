"""Tests for scripts/wdl-resources.py: static folding on an inline WDL fixture and facts about the real closure."""
from __future__ import annotations

import importlib.util
import pathlib
import sys
import tempfile
import unittest

REPO = pathlib.Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "wdl-resources.py"
G = 1024 ** 3

FIXTURE = """version 1.0
task lit {
  Int threads = 4
  Int mem_gb = threads * 2
  command <<< echo ~{threads} >>>
  output { String o = read_string(stdout()) }
  runtime { docker: "x" cpu: threads memory: "~{mem_gb} GiB" }
}
task withinput {
  input { Int n  Int t = 2 }
  Int mem_gb = 8 + n
  command <<< echo hi >>>
  output { String o = read_string(stdout()) }
  runtime { docker: "x" cpu: t memory: "~{mem_gb} GiB" }
}
task sized {
  input { File f }
  Int disk = ceil(size(f, "GB"))
  Int mem_gb = 4
  Int sort_mem = mem_gb - 1
  command <<< echo ~{disk} ~{sort_mem} >>>
  output { String o = read_string(stdout()) }
  runtime { docker: "x" cpu: 1 memory: "~{mem_gb} GiB" disk: "~{disk} GB" }
}
task scaled {
  input { Int mem_gb }
  command <<< echo hi >>>
  output { String o = read_string(stdout()) }
  runtime { docker: "x" cpu: 2 memory: "~{mem_gb} GB" gpuCount: 1 }
}
workflow main {
  input { File f  Int wf_mem = 16  Array[String] ids  Int? mem_override }
  Int n_samples = length(ids)
  Int total = 64
  Int shards = 8
  Int per = total / shards
  Int big_mem = select_first([mem_override, 8 + ceil(n_samples * 0.1)])
  call lit
  call withinput { input: n = per }
  call withinput as w2 { input: n = wf_mem }
  scatter (i in ids) { call sized { input: f = f } }
  if (length(ids) > 1) { call scaled { input: mem_gb = big_mem } }
  output { String o = lit.o }
}
"""


def setUpModule():
    try:
        import WDL  # noqa: F401
    except ImportError:
        raise unittest.SkipTest("miniwdl not importable (activate .venv)")


def load_script():
    spec = importlib.util.spec_from_file_location("wdl_resources", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules["wdl_resources"] = mod
    spec.loader.exec_module(mod)
    return mod


class FixtureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import WDL
        cls.mod = load_script()
        cls.tmp = tempfile.TemporaryDirectory()
        path = pathlib.Path(cls.tmp.name) / "main.wdl"
        path.write_text(FIXTURE)
        doc = WDL.load(str(path))
        cls.tasks = {}
        cls.mod.walk(doc, "ugc_wgw_fixture", cls.tasks, cls.mod.graph.collect_docs(doc))
        for info in cls.tasks.values():
            cls.mod.evaluate(info)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_literals_and_arithmetic_fold(self):
        t = self.tasks["lit"]
        self.assertEqual((t.cpu, t.memory), (4, 8 * G))
        self.assertEqual(t.notes, ["command uses threads"])

    def test_call_site_and_input_default_with_varying_sites(self):
        t = self.tasks["withinput"]
        self.assertEqual(t.cpu, 2)
        self.assertIsNone(t.memory)  # 16G from the constant call site, 24G from the workflow input default
        self.assertIn("cpu from input default", t.notes)
        self.assertTrue(any(n.startswith("memory varies by call site: 16G, 24G") for n in t.notes), t.notes)
        self.assertIn("memory from call site", t.notes)
        self.assertIn("memory from workflow input default", t.notes)
        self.assertEqual(len(t.call_sites), 2)

    def test_size_is_dynamic_but_memory_folds(self):
        t = self.tasks["sized"]
        self.assertEqual((t.cpu, t.memory), (1, 4 * G))
        self.assertIn("command uses mem", t.notes)  # through sort_mem = mem_gb - 1
        self.assertNotIn("command uses threads", t.notes)

    def test_n_scaled_memory_shows_formula_and_gpu(self):
        t = self.tasks["scaled"]
        self.assertEqual(t.cpu, 2)
        self.assertIsNone(t.memory)
        joined = "; ".join(t.notes)
        self.assertIn('memory: "~{mem_gb} GB" where mem_gb = big_mem = select_first([mem_override, 8 + ceil(n_samples * 0.1)])',
                      joined)
        self.assertIn("memory scales with N", t.notes)
        self.assertIn("gpu", t.notes)

    def test_outputs_render(self):
        tsv = self.mod.tsv_text(self.tasks)
        lines = [l for l in tsv.splitlines() if not l.startswith("#")]
        self.assertEqual(lines[0].split("\t"), list(self.mod.TSV_COLUMNS))
        row = dict(zip(lines[0].split("\t"), next(l for l in lines if l.startswith("lit\t")).split("\t")))
        self.assertEqual((row["cpu"], row["memory"], row["time"], row["stages"]), ("4", "8G", "-", "fixture"))
        self.assertTrue(row["source"].endswith("main.wdl:2"))
        md = self.mod.md_text(self.tasks)
        self.assertIn("| `lit` | 4 | 8G | command uses threads |", md)
        self.assertIn("## Largest requests", md)

    def test_duplicate_task_name_is_an_error(self):
        import WDL
        dup = FIXTURE.replace("task sized {", "task lit2 {").replace("call sized", "call lit2")
        other = "version 1.0\ntask lit { command <<< echo >>> runtime { docker: \"y\" cpu: 1 memory: \"1 GiB\" } }\n"
        with tempfile.TemporaryDirectory() as tmp:
            (pathlib.Path(tmp) / "other.wdl").write_text(other)
            main = dup.replace("version 1.0\n", 'version 1.0\nimport "other.wdl" as O\n', 1)
            main = main.replace("  call lit\n", "  call lit\n  call O.lit as lit_other\n")
            p = pathlib.Path(tmp) / "main.wdl"
            p.write_text(main)
            doc = WDL.load(str(p))
            with self.assertRaises(SystemExit) as cm:
                self.mod.walk(doc, "ugc_wgw_fixture", {}, self.mod.graph.collect_docs(doc))
            self.assertIn("defined twice", str(cm.exception))


class RealClosureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mod = load_script()
        cls.tasks = cls.mod.inventory()

    def test_known_requests(self):
        dv = self.tasks["deepvariant_call_variants_cpu"]
        self.assertEqual((dv.cpu, dv.memory), (64, 28 * G))  # v4.0.0: 28 GiB, 64 threads capped by max_nproc
        self.assertIn("cpu from workflow input default", dv.notes)  # threads = min(max_nproc, 64) since v4.0.0
        asm = self.tasks["ugc_wgw_hifiasm_assemble"]
        self.assertEqual((asm.cpu, asm.memory), (48, 288 * G))
        self.assertIn("cpu from workflow input default", asm.notes)
        gl = self.tasks["glnexus"]
        self.assertEqual(gl.cpu, 32)
        self.assertIsNone(gl.memory)
        self.assertIn("memory scales with N", gl.notes)
        self.assertIn("command uses mem", gl.notes)
        self.assertEqual(gl.source, "workflows/overrides/glnexus.wdl:9")
        self.assertEqual(sorted(gl.stages), ["ugc_wgw_cohort_call", "ugc_wgw_cohort_merge"])
        self.assertTrue(all(len(t.stages) >= 1 for t in self.tasks.values()))
        self.assertGreaterEqual(len(self.tasks), 45)

    def test_generated_files_current_and_deterministic(self):
        self.assertEqual(self.mod.tsv_text(self.tasks), self.mod.TSV_OUTPUT.read_text())
        self.assertEqual(self.mod.md_text(self.tasks), self.mod.MD_OUTPUT.read_text())
        again = self.mod.inventory()
        self.assertEqual(self.mod.tsv_text(again), self.mod.tsv_text(self.tasks))

    def test_tsv_rows_parse_as_policy(self):
        from ugc_wgw_miniwdl import policy
        rules = policy.parse_policy(str(self.mod.TSV_OUTPUT))
        self.assertGreaterEqual(len(rules), 45)
        by = {r.pattern: r for r in rules}
        self.assertEqual(by["deepvariant_call_variants_cpu"].cpu, 64)
        self.assertEqual(by["deepvariant_call_variants_cpu"].memory, 28 * G)
        self.assertIsNone(by["glnexus"].memory)


if __name__ == "__main__":
    unittest.main()
