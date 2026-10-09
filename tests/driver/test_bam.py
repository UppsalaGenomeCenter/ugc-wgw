import tempfile
import unittest
from pathlib import Path

from ugc_wgw import bam, preflight
from ugc_wgw.db import DB

from .helpers import make_project, write_bam, write_tsv


class BamInspectTest(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def test_small_file_is_read_whole_and_exact(self):
        info = bam.inspect(write_bam(self.tmp / "a.bam", reads=10, read_len=1000, movie="m1", sample="S1"))
        self.assertTrue(info.ok, info.problems)
        self.assertEqual((info.reads, info.bases, info.reads_exact, info.pbi), (10, 10_000, True, False))
        self.assertEqual((info.movies, info.samples, info.read_groups, info.aligned, info.n_ref), (["m1"], ["S1"], 1, False, 0))
        self.assertTrue(info.eof_ok)
        self.assertEqual(info.mean_read_length, 1000)
        d = info.to_dict()
        self.assertEqual((d["gb"], d["mean_read_length"], d["bytes"]), (0.0, 1000, info.bytes))
        self.assertEqual(bam.BamInfo.from_dict(d).reads, 10)

    def test_large_file_is_estimated_from_a_sample(self):
        info = bam.inspect(write_bam(self.tmp / "big.bam", reads=3000, read_len=15_000), sample_records=500)
        self.assertTrue(info.ok)
        self.assertEqual(info.records_sampled, 500)
        self.assertFalse(info.reads_exact)
        self.assertAlmostEqual(info.reads, 3000, delta=60)          # within 2 %
        self.assertAlmostEqual(info.bases / 1e6, 45.0, delta=1.0)
        self.assertAlmostEqual(info.coverage, 45e6 / 3.1e9, places=3)

    def test_pbi_gives_the_exact_count(self):
        info = bam.inspect(write_bam(self.tmp / "p.bam", reads=5, read_len=100, pbi_reads=12345))
        self.assertEqual((info.reads, info.reads_exact, info.pbi, info.bases), (12345, True, True, 1_234_500))
        self.assertIsNone(bam.pbi_reads(self.tmp / "none.bam"))

    def test_problems(self):
        info = bam.inspect(write_bam(self.tmp / "e.bam", empty=True))
        self.assertEqual(info.problems, ["empty file (0 bytes)"])
        info = bam.inspect(write_bam(self.tmp / "h.bam", header_only=True))
        self.assertEqual(info.problems, ["no reads (header only)"])
        self.assertEqual(info.movies, ["m84000_240101_000000_s1"])   # the header was still read
        info = bam.inspect(write_bam(self.tmp / "t.bam", reads=300, truncated=True))
        self.assertEqual(len(info.problems), 1)
        self.assertIn("truncated", info.problems[0])
        self.assertFalse(info.eof_ok)
        (self.tmp / "x.bam").write_text("not a bam\n")
        self.assertIn("truncated", bam.inspect(self.tmp / "x.bam").problems[0])
        self.assertIn("not readable", bam.inspect(self.tmp / "missing.bam").problems[0])
        # a complete BGZF stream that is not a BAM
        from .helpers import BGZF_EOF, bgzf_block
        (self.tmp / "g.bam").write_bytes(bgzf_block(b"BAX\x01" + bytes(40)) + BGZF_EOF)
        self.assertEqual(bam.inspect(self.tmp / "g.bam").problems, ["unreadable BAM: not a BAM file (bad magic)"])

    def test_warnings(self):
        info = bam.inspect(write_bam(self.tmp / "al.bam", aligned=True))
        self.assertTrue(info.ok)
        self.assertTrue(info.aligned)
        self.assertIn("aligned input (1 reference sequences)", info.warnings[0])
        info = bam.inspect(write_bam(self.tmp / "rg.bam", read_group=False))
        self.assertEqual(info.warnings, ["no @RG read group in the header"])
        self.assertEqual(info.movies, [])

    def test_quick_check(self):
        p = write_bam(self.tmp / "q.bam", reads=20)
        self.assertIsNone(bam.quick_check(p))
        self.assertIsNone(bam.quick_check(p, p.stat().st_size))
        self.assertIn("size changed since registration", bam.quick_check(p, p.stat().st_size + 1))
        self.assertIn("truncated", bam.quick_check(write_bam(self.tmp / "t.bam", truncated=True)))
        self.assertEqual(bam.quick_check(write_bam(self.tmp / "e.bam", empty=True)), "empty file (0 bytes)")
        self.assertIn("not readable", bam.quick_check(self.tmp / "nope.bam"))

    def test_formatters(self):
        self.assertEqual((bam.fmt_bases(3_300_000_000), bam.fmt_bases(2_500_000), bam.fmt_bases(900)), ("3.3 Gb", "2.5 Mb", "1 kb"))
        self.assertEqual((bam.fmt_bytes(0), bam.fmt_bytes(2048), bam.fmt_bytes(1.5 * 1024 ** 3)), ("0 B", "2.0 KB", "1.5 GB"))


class PreflightTest(unittest.TestCase):
    def test_read_paths_and_check(self):
        tmp = Path(tempfile.mkdtemp())
        cfg = make_project(tmp)
        from ugc_wgw import samples
        from ugc_wgw.log import Events
        db = DB(cfg.db_path)
        events = Events(db, cfg.logs_dir / "events.jsonl")
        tsv = write_tsv(tmp, [{"sample_id": "S1", "hifi_reads": "a.bam,b.bam", "fail_reads": "f.bam"}])
        samples.add_samples(db, events, tsv)
        rec = db.get_sample("S1")
        doc = {"ugc_wgw_singleton.sample_id": "S1", "ugc_wgw_singleton.hifi_reads": rec.hifi_reads,
               "ugc_wgw_singleton.fail_reads": rec.fail_reads, "ugc_wgw_singleton.threads": 4,
               "ugc_wgw_assembly.father_hifi_reads": ["/nonexistent/p.bam"]}
        self.assertEqual(preflight.read_paths(doc), rec.hifi_reads + rec.fail_reads + ["/nonexistent/p.bam"])
        n, total, problems = preflight.check(doc, db)
        self.assertEqual(n, 4)
        self.assertEqual(total, sum(Path(p).stat().st_size for p in rec.hifi_reads + rec.fail_reads))
        self.assertEqual(len(problems), 1)
        self.assertIn("/nonexistent/p.bam: not readable", problems[0])
        # the registered size is what is compared
        Path(rec.hifi_reads[0]).write_bytes(Path(rec.hifi_reads[0]).read_bytes() + b"\0")
        _, _, problems = preflight.check({"ugc_wgw_singleton.hifi_reads": rec.hifi_reads}, db)
        self.assertEqual(len(problems), 1)
        self.assertIn("size changed since registration", problems[0])
        self.assertEqual(preflight.check({"ugc_wgw_singleton.threads": 4}, db), (0, 0, []))
        db.close()


if __name__ == "__main__":
    unittest.main()
