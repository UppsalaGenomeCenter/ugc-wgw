"""The chrom_coverage task's Python (between its begin/end markers in the WDL) on synthetic PAFs: per-chromosome
coverage, the karyotype decision, and dipcall's regrouped haplotypes for a non-trio XY sample."""
import gzip
import json
import pathlib
import subprocess
import sys
import tempfile
import unittest

from .helpers import REPO

WDL = REPO / "workflows" / "ugc_wgw" / "assembly" / "chrom_coverage.wdl"
BEGIN, END = "# ---- chrom_coverage.py (begin) ----\n", "# ---- chrom_coverage.py (end) ----\n"
CHR = {"chr1": 248_956_422, "chrX": 156_040_895, "chrY": 57_227_415, "chr20": 64_444_167}
PAR = "chrX\t0\t2781479\nchrX\t155701383\t156030895\n"


def script() -> str:
    text = WDL.read_text()
    assert text.count(BEGIN) == 1 and text.count(END) == 1
    return text[text.index(BEGIN) + len(BEGIN):text.index(END)]


def paf_line(q: str, qlen: int, t: str, ts: int, te: int, mapq: int = 60, tp: str = "P") -> str:
    n = te - ts
    return "\t".join([q, str(qlen), "0", str(n), "+", t, str(CHR[t]), str(ts), str(te), str(n), str(n), str(mapq), f"tp:A:{tp}", "cs:Z::" + str(n)]) + "\n"


class ChromCoverageTest(unittest.TestCase):
    def setUp(self):
        self.tmp = pathlib.Path(tempfile.mkdtemp())
        (self.tmp / "chrom_coverage.py").write_text(script())
        (self.tmp / "ref.fai").write_text("".join(f"{c}\t{n}\t0\t60\t61\n" for c, n in CHR.items()))
        (self.tmp / "par.bed").write_text(PAR)

    def fasta(self, name: str, contigs: list[str]) -> pathlib.Path:
        p = self.tmp / f"{name}.fasta.gz"
        with gzip.open(p, "wt") as fh:
            for c in contigs:
                fh.write(f">{c} extra words\nACGT\nACGT\n")
        return p

    def run_script(self, haps: list[tuple[str, str, str | None]], trio: bool = False, x_min: str = "0.5", y_min: str = "0.1") -> dict:
        regroup = self.tmp / "regrouped"; regroup.mkdir(exist_ok=True)
        specs = [f"{label}={paf}" + (f"={fa}" if fa else "") for label, paf, fa in haps]
        res = subprocess.run([sys.executable, str(self.tmp / "chrom_coverage.py"), "S1", str(self.tmp / "ref.fai"), str(self.tmp / "par.bed"),
                              "true" if trio else "false", "5", x_min, y_min, str(self.tmp / "cov.tsv"), str(self.tmp / "kar.json"),
                              str(self.tmp / "kar.txt"), str(self.tmp / "regroup.tsv"), str(regroup), *specs],
                             capture_output=True, text=True)
        self.assertEqual(res.returncode, 0, res.stderr)
        doc = json.loads((self.tmp / "kar.json").read_text())
        self.assertEqual((self.tmp / "kar.txt").read_text().strip(), doc["karyotype"])
        return doc

    def write_paf(self, name: str, lines: list[str]) -> str:
        p = self.tmp / f"{name}.paf.gz"
        with gzip.open(p, "wt") as fh:
            fh.write("".join(lines))
        return str(p)

    def test_coverage_rows_and_unknown_for_a_chr20_only_assembly(self):
        h1 = self.write_paf("h1", [paf_line("c1", 30_000_000, "chr20", 0, 30_000_000), paf_line("c2", 20_000_000, "chr20", 30_000_000, 50_000_000),
                                   paf_line("low", 1_000_000, "chr20", 50_000_000, 51_000_000, mapq=1),
                                   paf_line("sec", 1_000_000, "chr20", 50_000_000, 51_000_000, tp="S")])
        h2 = self.write_paf("h2", [paf_line("d1", 60_000_000, "chr20", 0, 60_000_000)])
        doc = self.run_script([("hap1", h1, None), ("hap2", h2, None)])
        self.assertEqual(doc["karyotype"], "unknown")
        self.assertIn("below 0.5", doc["reason"])
        self.assertFalse(doc["regrouped"])
        rows = [l.split("\t") for l in (self.tmp / "cov.tsv").read_text().splitlines()[1:]]
        self.assertEqual([(r[0], r[1], int(r[2]), int(r[3])) for r in rows], [("hap1", "chr20", 50_000_000, 2), ("hap2", "chr20", 60_000_000, 1)])
        self.assertAlmostEqual(float(rows[1][5]), 60_000_000 / CHR["chr20"], places=3)

    def test_xx_and_xy(self):
        x_len, y_len = CHR["chrX"], CHR["chrY"]
        # XX: both haplotypes cover chrX (PAR bases excluded), nothing on chrY
        h1 = self.write_paf("xx1", [paf_line("x1", x_len, "chrX", 0, x_len)])
        h2 = self.write_paf("xx2", [paf_line("x2", x_len, "chrX", 0, x_len)])
        doc = self.run_script([("hap1", h1, None), ("hap2", h2, None)])
        self.assertEqual(doc["karyotype"], "XX")
        self.assertEqual(doc["chrx_bases_outside_par"], 2 * (x_len - 2781479 - (156030895 - 155701383)))
        self.assertEqual(doc["chry_bases_outside_par"], 0)
        # XY: one chrX, chrY above a tenth of its length
        h1 = self.write_paf("xy1", [paf_line("x1", x_len, "chrX", 0, x_len)])
        h2 = self.write_paf("xy2", [paf_line("y1", 20_000_000, "chrY", 10_000_000, 30_000_000)])
        doc = self.run_script([("hap1", h1, None), ("hap2", h2, None)], trio=True)
        self.assertEqual(doc["karyotype"], "XY")
        self.assertAlmostEqual(doc["chry_fraction"], 20_000_000 / y_len, places=3)
        self.assertFalse(doc["regrouped"])   # trio-binned: nothing moves

    def test_non_trio_xy_regroups_chrx_and_chry_contigs_for_dipcall(self):
        x_len = CHR["chrX"]
        # hap1 carries the chrX contig and an autosome; hap2 carries the chrY contig and a PAR-only chrX piece
        h1 = self.write_paf("r1", [paf_line("xcontig", x_len, "chrX", 2_800_000, x_len - 400_000), paf_line("a1", 50_000_000, "chr1", 0, 50_000_000)])
        h2 = self.write_paf("r2", [paf_line("ycontig", 20_000_000, "chrY", 10_000_000, 30_000_000), paf_line("parcontig", 2_000_000, "chrX", 0, 2_000_000),
                                   paf_line("a2", 50_000_000, "chr1", 50_000_000, 100_000_000)])
        f1 = self.fasta("h1", ["xcontig", "a1"]); f2 = self.fasta("h2", ["ycontig", "parcontig", "a2"])
        doc = self.run_script([("hap1", h1, str(f1)), ("hap2", h2, str(f2))])
        self.assertEqual(doc["karyotype"], "XY")
        self.assertTrue(doc["regrouped"])
        self.assertEqual(doc["moved_contigs"], 2)
        log = [l.split("\t") for l in (self.tmp / "regroup.tsv").read_text().splitlines()[1:]]
        self.assertEqual({(r[0], r[1], r[3]) for r in log}, {("hap1", "xcontig", "hap2"), ("hap2", "ycontig", "hap1")})
        def names(p):
            with gzip.open(p, "rt") as fh:
                return [l[1:].split()[0] for l in fh if l.startswith(">")]
        self.assertEqual(names(self.tmp / "regrouped" / "S1.dipcall.hap1.fasta.gz"), ["a1", "ycontig"])
        self.assertEqual(names(self.tmp / "regrouped" / "S1.dipcall.hap2.fasta.gz"), ["xcontig", "parcontig", "a2"])   # PAR-only chrX stays


if __name__ == "__main__":
    unittest.main()
