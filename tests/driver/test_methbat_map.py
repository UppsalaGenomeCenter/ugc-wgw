"""Tests for scripts/make-methbat-map.py: CpG islands joined with the genes near them and typed."""
from __future__ import annotations

import gzip
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "make-methbat-map.py"

GTF = """\
##description: toy
chr1\tHAVANA\tgene\t1001\t5000\t.\t+\t.\tgene_id "ENSG1"; gene_type "protein_coding"; gene_name "GENE1"; level 2; hgnc_id "HGNC:1";
chr1\tHAVANA\ttranscript\t1001\t5000\t.\t+\t.\tgene_id "ENSG1"; transcript_id "ENST1"; gene_name "GENE1";
chr1\tHAVANA\tgene\t1501\t4000\t.\t+\t.\tgene_id "ENSG3"; gene_type "protein_coding"; gene_name "GENE3"; level 2; hgnc_id "HGNC:3";
chr1\tHAVANA\tgene\t6001\t9000\t.\t-\t.\tgene_id "ENSG2"; gene_type "lncRNA"; gene_name "GENE2"; level 2; hgnc_id "HGNC:2";
chr1\tHAVANA\tgene\t20001\t20500\t.\t+\t.\tgene_id "ENSG4"; gene_type "processed_pseudogene"; gene_name "PSEUDO1"; level 2;
chr1\tHAVANA\tgene\t30001\t30800\t.\t+\t.\tgene_id "ENSG5"; gene_type "lncRNA"; gene_name "LINC-NAMED"; level 2;
chr2\tHAVANA\tgene\t101\t900\t.\t+\t.\tgene_id "ENSG6"; gene_type "protein_coding"; gene_name "GENE6"; level 2; hgnc_id "HGNC:6";
"""

REGIONS = """\
chrom\tstart\tend\tcpg_label
chr1\t900\t1100\tCpG:_10
chr1\t3000\t3200\tCpG:_11
chr1\t50000\t50200\tCpG:_12
chr1\t7000\t7100\tCpG:_13
chr1\t30900\t31000\tCpG:_14
chr3\t10\t20\tCpG:_15
"""

IMPRINTING = """\
chrom\tstart\tend\tcpg_label\texpected_category\tanomalous_categories
chr1\t6900\t7050\tTOY:IG-DMR\tAlleleSpecificMethylation\tMethylated;Unmethylated
"""


def run(tmp: Path, *extra: str, gz: bool = False) -> tuple[list[str], list[dict[str, str]], str]:
    regions = tmp / "islands.tsv"
    regions.write_text(REGIONS)
    imprinting = tmp / "imprinting.tsv"
    imprinting.write_text(IMPRINTING)
    if gz:
        gtf = tmp / "toy.gtf.gz"
        with gzip.open(gtf, "wt") as fh:
            fh.write(GTF)
    else:
        gtf = tmp / "toy.gtf"
        gtf.write_text(GTF)
    res = subprocess.run([sys.executable, str(SCRIPT), "--regions", str(regions), "--gtf", str(gtf), "--imprinting", str(imprinting), *extra],
                         capture_output=True, text=True)
    if res.returncode != 0:
        raise AssertionError(res.stderr)
    lines = res.stdout.splitlines()
    meta = [l for l in lines if l.startswith("##")]
    body = [l for l in lines if not l.startswith("##")]
    header = body[0].split("\t")
    rows = [dict(zip(header, l.split("\t"))) for l in body[1:]]
    return meta, rows, res.stderr


class MethbatMapTest(unittest.TestCase):
    def test_types_genes_and_header(self) -> None:
        with tempfile.TemporaryDirectory() as d:
            meta, rows, err = run(Path(d), gz=True)
        by_label = {r["region_label"]: r for r in rows}
        self.assertEqual(len(rows), 6, "every island keeps its row")
        self.assertEqual(list(rows[0]), ["chrom", "start", "end", "region_label", "hgnc_id", "gene_symbol", "region_type"])
        # a TSS (GENE1 at 1000) inside the island: promoter; GENE3's TSS 1500 and GENE2's 8999 (- strand) are within 10 kb
        self.assertEqual(by_label["CpG:_10"]["region_type"], "promoter")
        self.assertEqual(by_label["CpG:_10"]["gene_symbol"], "GENE1;GENE3;GENE2")
        self.assertEqual(by_label["CpG:_10"]["hgnc_id"], "HGNC:1;HGNC:3;HGNC:2")
        # inside GENE1 and GENE3, no TSS within 1 kb: gene_body; nearest TSS first (GENE3 1500 bp, GENE1 1900, GENE2 5800)
        self.assertEqual(by_label["CpG:_11"]["region_type"], "gene_body")
        self.assertEqual(by_label["CpG:_11"]["gene_symbol"], "GENE3;GENE1;GENE2")
        # nothing within 10 kb: intergenic with empty gene columns
        self.assertEqual(by_label["CpG:_12"]["region_type"], "intergenic")
        self.assertEqual((by_label["CpG:_12"]["hgnc_id"], by_label["CpG:_12"]["gene_symbol"]), ("", ""))
        # overlaps the imprinting region: imprinted wins over the gene body of GENE2 (- strand, TSS at 8999)
        self.assertEqual(by_label["CpG:_13"]["region_type"], "imprinted")
        self.assertEqual(by_label["CpG:_13"]["gene_symbol"], "GENE2;GENE3;GENE1")
        # the pseudogene and the unnamed lncRNA are not genes by the default rule
        self.assertEqual(by_label["CpG:_14"]["region_type"], "intergenic")
        self.assertEqual(by_label["CpG:_14"]["gene_symbol"], "")
        # a chromosome without genes
        self.assertEqual(by_label["CpG:_15"]["region_type"], "intergenic")
        self.assertTrue(any(m.startswith("##methbat_map_generator=") for m in meta))
        self.assertTrue(any(m.startswith("##methbat_map_input=gtf\ttoy.gtf.gz\tsha256:") for m in meta))
        self.assertIn("6 islands", err)
        self.assertIn("imprinted 1, promoter 1, gene_body 1, intergenic 3", err)

    def test_named_genes_and_windows(self) -> None:
        with tempfile.TemporaryDirectory() as d:
            _, rows, _ = run(Path(d), "--genes", "named", "--promoter", "200", "--tss-window", "300")
        by_label = {r["region_label"]: r for r in rows}
        # the named lncRNA counts now (TSS 30000, 900 bp from the island: outside the 300 bp window, and its body does not
        # overlap either) -> still intergenic; the pseudogene never counts
        self.assertEqual(by_label["CpG:_14"]["gene_symbol"], "")
        # a 300 bp window drops GENE3 (TSS 400 bp away) and GENE2 from the first island; GENE1 stays (TSS inside)
        self.assertEqual(by_label["CpG:_10"]["gene_symbol"], "GENE1")
        self.assertEqual(by_label["CpG:_10"]["region_type"], "promoter")
        # gene body overlap keeps genes whose TSS is outside the window
        self.assertEqual(by_label["CpG:_11"]["gene_symbol"], "GENE3;GENE1")
        self.assertEqual(by_label["CpG:_11"]["region_type"], "gene_body")


if __name__ == "__main__":
    unittest.main()
