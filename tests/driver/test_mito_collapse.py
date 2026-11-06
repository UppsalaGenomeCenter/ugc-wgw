"""The chrM collapse script embedded in workflows/ugc_wgw/mito/collapse.wdl (ugc_wgw_mito_collapse), run as the task
runs it: nallo's genotype rule, the SV split, the per-haplotype columns, the summary, the messages."""
from __future__ import annotations

import gzip
import json
import os
import pathlib
import subprocess
import sys
import tempfile
import unittest

from .helpers import REPO

WDL = REPO / "workflows" / "ugc_wgw" / "mito" / "collapse.wdl"
BEGIN, END = "# ---- mito_collapse.py (begin) ----\n", "# ---- mito_collapse.py (end) ----\n"
HEADER = ('##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
          '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">\n##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Read depth for each allele">\n'
          '##FORMAT=<ID=VAF,Number=A,Type=Float,Description="Variant allele fractions.">\n##contig=<ID=chrM,length=16569>\n'
          '##mitorsaw_version="0.2.13-49bcf18"\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n')


def script() -> str:
    text = WDL.read_text()
    assert text.count(BEGIN) == 1 and text.count(END) == 1
    body = text[text.index(BEGIN) + len(BEGIN):text.index(END)]
    return "\n".join(line[4:] if line.startswith("    ") else line for line in body.splitlines()) + "\n"


def rec(pos: int, ref: str, alt: str, gt: str, vaf: str = "1", dp: int = 100) -> str:
    return f"chrM\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\tGT:DP:AD:VAF\t{gt}:{dp}:10,90:{vaf}\n"


class CollapseTest(unittest.TestCase):
    def setUp(self):
        self.tmp = pathlib.Path(tempfile.mkdtemp())
        (self.tmp / "collapse.py").write_text(script())

    def run_it(self, records: str, haplotypes: list[tuple[str, float]], sv_min: int = 50, sample: str = "S1") -> dict[str, str]:
        with gzip.open(self.tmp / "in.vcf.gz", "wt") as fh:
            fh.write(HEADER + records)
        (self.tmp / "in.json").write_text(json.dumps({"haplotypes": [{"label": l, "estimated_abundance": a, "num_ref_variants": 0} for l, a in haplotypes],
                                                      "fingerprint_stats": {}}))
        res = subprocess.run([sys.executable, "collapse.py", "in.vcf.gz", "in.json", sample, "out", str(sv_min)], cwd=self.tmp, capture_output=True, text=True)
        self.assertEqual(res.returncode, 0, res.stderr)
        out = {name: (self.tmp / f"out.mt.{name}").read_text() for name in ("small_variants.vcf", "svs.vcf", "haplotypes.vcf", "haplotypes.tsv", "summary.tsv")}
        out["messages"] = (self.tmp / "messages.txt").read_text()
        out["stderr"] = res.stderr
        return out

    @staticmethod
    def gts(vcf: str) -> list[str]:
        return [line.split("\t")[9].split(":")[0] for line in vcf.splitlines() if not line.startswith("#")]

    def test_rule_and_split(self):
        haps = [("hap_0", 0.6), ("hap_1", 0.3), ("hap_2", 0.1)]
        out = self.run_it(rec(73, "A", "G", "1|1|1") + rec(152, "T", "C", "0|1|1", vaf="0.4") + rec(300, "G", "A", "0|0|0")
                          + rec(400, "C", "T", ".|1|.") + rec(500, "A", "T", ".|.|.") + rec(600, "A", "T,C", "1|2|2", vaf="0.3,0.7")
                          + rec(700, "A", "T,C", "0|1|2", vaf="0.3,0.1") + rec(800, "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC", "A", "0|0|1", vaf="0.1")
                          + rec(900, "A", "<DEL>", "1|1|1"), haps)
        self.assertEqual(self.gts(out["small_variants.vcf"]), ["1/1", "0/1", "0/0", "1/1", "./.", "1/2", "0/1"])   # 600: both alternate, the more abundant first; 700: 0 first, then the most abundant ALT
        self.assertEqual(self.gts(out["svs.vcf"]), ["0/1", "1/1"])
        self.assertIn("1 record(s) with more than two alleles", out["messages"])
        cols = [line for line in out["haplotypes.vcf"].splitlines() if line.startswith("#CHROM")][0].split("\t")
        self.assertEqual(cols[9:], ["S1.hap_0", "S1.hap_1", "S1.hap_2"])
        rows = [line.split("\t") for line in out["haplotypes.vcf"].splitlines() if not line.startswith("#")]
        self.assertEqual([r[8:] for r in rows][1], ["GT", "0", "1", "1"])
        self.assertEqual([r[8:] for r in rows][3], ["GT", ".", "1", "."])
        self.assertEqual(len(rows), 9)
        self.assertIn("##ugc_wgw_mito_abundances=hap_0:0.6000,hap_1:0.3000,hap_2:0.1000", out["small_variants.vcf"])
        self.assertIn("##ugc_wgw_mito_source=mitorsaw 0.2.13-49bcf18", out["svs.vcf"])
        self.assertNotIn("##FORMAT=<ID=VAF", out["haplotypes.vcf"])
        summ = dict(zip(*[line.split("\t") for line in out["summary.tsv"].splitlines()]))
        self.assertEqual((summ["sample_id"], summ["vcf_sample"], summ["haplotypes"], summ["records"]), ("S1", "S1", "3", "9"))
        self.assertEqual((summ["homoplasmic"], summ["heteroplasmic"], summ["svs"]), ("2", "3", "2"))   # 73 and 400 homoplasmic; 152, 600 (two alternates) and 700 heteroplasmic
        self.assertEqual((summ["het_vaf_mean"], summ["het_vaf_max"], summ["major_haplotype"], summ["major_abundance"]), ("0.3333", "0.4000", "hap_0", "0.6000"))
        self.assertEqual(out["haplotypes.tsv"].splitlines()[1:], ["hap_0\t0.6000\t3", "hap_1\t0.3000\t6", "hap_2\t0.1000\t6"])

    def test_single_haplotype_and_empty(self):
        out = self.run_it(rec(263, "A", "G", "1") + rec(310, "T", "TC", "1", vaf="0.92"), [("hap_0", 1.0)])
        self.assertEqual(self.gts(out["small_variants.vcf"]), ["1/1", "1/1"])
        self.assertEqual(out["messages"], "")
        out = self.run_it("", [("hap_0", 1.0)])
        self.assertEqual(self.gts(out["small_variants.vcf"]), [])
        self.assertIn("no chrM variant called", out["messages"])
        self.assertTrue(out["haplotypes.vcf"].rstrip().endswith("FORMAT\tS1.hap_0"))
        summ = dict(zip(*[line.split("\t") for line in out["summary.tsv"].splitlines()]))
        self.assertEqual((summ["records"], summ["het_vaf_mean"], summ["major_haplotype"]), ("0", "", "hap_0"))

    def test_haplotype_count_mismatch(self):
        out = self.run_it(rec(73, "A", "G", "1|0"), [("hap_0", 1.0)])
        self.assertIn("1 haplotypes in the JSON but 2 GT entries", out["messages"])
        self.assertEqual(self.gts(out["small_variants.vcf"]), ["0/1"])
        self.assertEqual(out["haplotypes.tsv"].splitlines()[1:], ["hap_0\t1.0000\t1", "hap_1\t0.0000\t0"])

    def test_real_trio_when_present(self):
        """The three mitorsaw outputs of the design pass (chrM slices of the GIAB trio, not in the repository): the known
        counts. UGC_WGW_MITORSAW_TRIO_DIR names the directory holding <id>.mitorsaw.vcf.gz and <id>.mitorsaw.json."""
        where = os.environ.get("UGC_WGW_MITORSAW_TRIO_DIR")
        mito = pathlib.Path(where) if where else None
        if mito is None or not (mito / "HG003.mitorsaw.vcf.gz").exists():
            self.skipTest("UGC_WGW_MITORSAW_TRIO_DIR is unset or lacks the trio's mitorsaw outputs")
        for sid, expect in (("HG002", ("1", "15", "15", "0")), ("HG003", ("4", "41", "38", "3")), ("HG004", ("6", "19", "14", "5"))):
            res = subprocess.run([sys.executable, "collapse.py", str(mito / f"{sid}.mitorsaw.vcf.gz"), str(mito / f"{sid}.mitorsaw.json"), sid, sid, "50"],
                                 cwd=self.tmp, capture_output=True, text=True)
            self.assertEqual(res.returncode, 0, res.stderr)
            summ = dict(zip(*[line.split("\t") for line in (self.tmp / f"{sid}.mt.summary.tsv").read_text().splitlines()]))
            self.assertEqual((summ["haplotypes"], summ["records"], summ["homoplasmic"], summ["heteroplasmic"]), expect, sid)


if __name__ == "__main__":
    unittest.main()
