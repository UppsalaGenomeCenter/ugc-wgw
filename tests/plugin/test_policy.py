"""Tests for ugc_wgw_miniwdl.policy: the TSV format, units, matching and application."""
from __future__ import annotations

import os
import tempfile
import unittest

from ugc_wgw_miniwdl import policy

G = 1024 ** 3

HEADER = "task\tcpu\tmemory\ttime\tpartition\tconstraint\n"
SAMPLE = HEADER + (
    "# comment\n"
    "\n"
    "deepvariant_call_variants_cpu\t48\t-\t-\t-\t-\n"
    "ugc_wgw_hifiasm_assemble\t-\t-\t2-00:00:00\tfat\t-\n"
    "glnexus\t-\t200G\t\t\t\n"
    "deepvariant_*\t-\t-\t12:00:00\t-\t-\n"
)


class ParseTests(unittest.TestCase):
    def test_rows_and_unset_cells(self):
        rules = policy.parse_policy_text(SAMPLE, "p.tsv")
        self.assertEqual([r.pattern for r in rules],
                         ["deepvariant_call_variants_cpu", "ugc_wgw_hifiasm_assemble", "glnexus", "deepvariant_*"])
        self.assertEqual(rules[0].values(), {"cpu": 48})
        self.assertEqual(rules[1].values(), {"time": 2880, "partition": "fat"})
        self.assertEqual(rules[2].values(), {"memory": 200 * G})
        self.assertEqual(rules[3].values(), {"time": 720})
        self.assertEqual(rules[0].source, "p.tsv:4")   # line numbers count comments and blanks
        self.assertEqual(rules[3].source, "p.tsv:7")

    def test_short_rows_and_info_columns(self):
        text = "task\tcpu\tstages\tsource\tnotes\nmosdepth\t4\tsingleton\tx.wdl:5\tcommand uses threads\nbam_stats\n"
        rules = policy.parse_policy_text(text, "p")
        self.assertEqual(rules[0].values(), {"cpu": 4})
        self.assertEqual(rules[1].values(), {})

    def test_errors_name_file_and_line(self):
        cases = {
            "cpu\ttask\nx\t1\n": "p:1: the first column",
            "task\tcores\nx\t1\n": "p:1: unknown column(s) cores",
            "task\tcpu\tcpu\nx\t1\t1\n": "p:1: duplicate column",
            "task\tcpu\nx\t1\t2\n": "p:2: 3 cells for 2 columns",
            "task\tcpu\n\t1\n": "p:2: bad task name",
            "task\tcpu\nx\tfour\n": "p:2: column cpu: cpu must be a positive integer",
            "task\tcpu\nx\t0\n": "p:2: column cpu",
            "task\tmemory\nx\t256\n": "p:2: column memory: memory needs a number and a unit",
            "task\tmemory\nx\t256X\n": "p:2: column memory: unknown memory unit",
            "task\ttime\nx\t1:2:3:4\n": "p:2: column time: time must be",
            "task\ttime\nx\t0\n": "p:2: column time: time must be at least one minute",
            "x\t1\n": "p:1: the first column",
        }
        for text, needle in cases.items():
            with self.assertRaises(policy.PolicyError, msg=text) as cm:
                policy.parse_policy_text(text, "p")
            self.assertIn(needle, str(cm.exception), text)

    def test_comment_only_file_is_empty(self):
        self.assertEqual(policy.parse_policy_text("# nothing\n\n", "p"), [])
        self.assertEqual(policy.parse_policy_text("", "p"), [])
        self.assertEqual(policy.parse_policy_text(HEADER, "p"), [])

    def test_load_policy_caches_by_mtime_and_size(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "r.tsv")
            with open(path, "w") as fh:
                fh.write(HEADER + "a\t2\t-\t-\t-\t-\n")
            first = policy.load_policy(path)
            self.assertIs(policy.load_policy(path), first)
            with open(path, "w") as fh:
                fh.write(HEADER + "a\t3\t-\t-\t-\t-\n")
            os.utime(path, (1_700_000_000, 1_700_000_000))
            second = policy.load_policy(path)
            self.assertEqual(second[0].cpu, 3)
            self.assertEqual(policy.load_policy(os.path.join(tmp, "missing.tsv")), [])


class UnitTests(unittest.TestCase):
    def test_memory_units_are_binary(self):
        for text, expect in [("256G", 256 * G), ("256GB", 256 * G), ("256GiB", 256 * G), ("256 GiB", 256 * G),
                             ("1.5T", 1536 * G), ("512M", 512 * 1024 ** 2), ("1024K", 1024 ** 2), ("7B", 7)]:
            self.assertEqual(policy.parse_memory(text), expect, text)
        self.assertEqual(policy.format_memory(256 * G), "256G")
        self.assertEqual(policy.format_memory(1536 * G), "1536G")
        self.assertEqual(policy.format_memory(2048 * G), "2T")
        self.assertEqual(policy.format_memory(512 * 1024 ** 2), "512M")
        self.assertEqual(policy.format_memory(G + 1), "1.00G")

    def test_time_follows_sbatch(self):
        for text, expect in [("90", 90), ("30:30", 31), ("1:30:00", 90), ("2-00", 2880), ("2-12", 3600),
                             ("2-12:30", 3630), ("2-00:00:00", 2880), ("0-00:00:01", 1), ("00:00:30", 1)]:
            self.assertEqual(policy.parse_time(text), expect, text)
        self.assertEqual(policy.format_time(2880), "2-00:00:00")
        self.assertEqual(policy.format_time(90), "01:30:00")
        self.assertEqual(policy.format_time(3630), "2-12:30:00")


class ResolveTests(unittest.TestCase):
    def setUp(self):
        self.rules = policy.parse_policy_text(SAMPLE, "p")

    def test_all_matching_rows_apply_later_wins(self):
        r = policy.resolve(self.rules, "deepvariant_call_variants_cpu")
        self.assertEqual({k: v for k, (v, _) in r.items()}, {"cpu": 48, "time": 720})
        self.assertEqual(r["time"][1].source, "p:7")
        self.assertEqual(policy.resolve(self.rules, "mosdepth"), {})
        later = policy.parse_policy_text(HEADER + "*\t2\t-\t-\t-\t-\nmosdepth\t8\t-\t-\t-\t-\n", "q")
        self.assertEqual(policy.resolve(later, "mosdepth")["cpu"][0], 8)
        self.assertEqual(policy.resolve(later, "trgt")["cpu"][0], 2)

    def test_apply_maps_keys_and_reports_changes(self):
        rv = {"cpu": 64, "memory_reservation": 256 * G, "time_minutes": 4320}
        changes, sources = policy.apply(self.rules, "deepvariant_call_variants_cpu", rv)
        self.assertEqual(rv, {"cpu": 48, "memory_reservation": 256 * G, "time_minutes": 720})
        self.assertEqual(changes, {"cpu": (64, 48), "time_minutes": (4320, 720)})
        self.assertEqual(sources, ["p:4", "p:7"])
        rv = {"cpu": 48, "memory_reservation": 288 * G}
        changes, _ = policy.apply(self.rules, "ugc_wgw_hifiasm_assemble", rv)
        self.assertEqual(rv["slurm_partition"], "fat")
        self.assertEqual(rv["time_minutes"], 2880)
        self.assertEqual(changes, {"time_minutes": (None, 2880), "slurm_partition": (None, "fat")})
        self.assertEqual(policy.describe(changes), "time -→2-00:00:00, partition -→fat")

    def test_apply_no_change_when_equal(self):
        rv = {"cpu": 48}
        changes, sources = policy.apply(self.rules, "deepvariant_call_variants_cpu", rv)
        self.assertEqual(changes, {"time_minutes": (None, 720)})
        self.assertEqual(sources, ["p:7"])

    def test_apply_bounded_by_backend_limits_and_rescales_memory_limit(self):
        rv = {"cpu": 2, "memory_reservation": 4 * G, "memory_limit": 8 * G}
        rules = policy.parse_policy_text(HEADER + "x\t64\t200G\t-\t-\t-\n", "p")
        changes, _ = policy.apply(rules, "x", rv, {"cpu": 24, "mem_bytes": 100 * G})
        self.assertEqual(rv["cpu"], 24)
        self.assertEqual(rv["memory_reservation"], 100 * G)
        self.assertEqual(rv["memory_limit"], 200 * G)
        self.assertEqual(changes["cpu"], (2, 24))
        rv = {"cpu": 2}
        policy.apply(rules, "x", rv, {"cpu": 2 ** 62, "mem_bytes": 2 ** 62})  # miniwdl-slurm: no limit
        self.assertEqual(rv["cpu"], 64)
        self.assertEqual(rv["memory_reservation"], 200 * G)

    def test_describe_memory_and_cpu(self):
        self.assertEqual(policy.describe({"cpu": (64, 48), "memory_reservation": (256 * G, 200 * G)}),
                         "cpu 64→48, memory 256G→200G")


if __name__ == "__main__":
    unittest.main()
