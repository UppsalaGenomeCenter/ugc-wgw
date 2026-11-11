"""Tests for scripts/make-strchive-map.py: STRchive's thresholds joined to the TRGT catalog's named loci, the motifs
that count as pathogenic in the catalog's spelling, and the ordinary loci overlapping a named one."""
from __future__ import annotations

import gzip
import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "make-strchive-map.py"

CATALOG = """\
chr4\t100\t200\tID=HD_HTT;MOTIFS=CAG,CCG;STRUC=<TR>
chr4\t90\t260\tID=chr4_90_260_GCC;MOTIFS=GCC;STRUC=(GCC)n
chr4\t5000\t5060\tID=CANVAS_RFC1;MOTIFS=AAGGG,ACAGG,AAAAG,AAGAG;STRUC=<TR>;INCLUDE_FAIL_READS
chr4\t9000\t9040\tID=chr4_9000_9040_A;MOTIFS=A;STRUC=(A)n
chr20\t300\t420\tID=CJD_PRNP;MOTIFS=CCTCAGGG,CCTCATGG,GGCCTCAT;STRUC=<TR>
chrX\t700\t760\tID=PRTS_ARX;MOTIFS=NGC;STRUC=<TR>
chr1\t10\t30\tID=HMNR7_VWA1;MOTIFS=GGCGCGGAGC;STRUC=<TR>
"""


def locus(lid: str, gene: str, chrom: str, start: int, stop: int, pathogenic: list[str], **kw: object) -> dict[str, object]:
    d: dict[str, object] = {"id": lid, "disease_id": lid.split("_")[0], "gene": gene, "disease": f"{gene} disease", "inheritance": ["AD"],
                            "chrom": chrom, "start_hg38": start, "stop_hg38": stop, "motif_len": len(pathogenic[0]),
                            "reference_motif_reference_orientation": [pathogenic[0]], "pathogenic_motif_reference_orientation": pathogenic,
                            "benign_min": 5.0, "benign_max": 26.0, "intermediate_min": 27.0, "intermediate_max": 39.0, "pathogenic_min": 40.0,
                            "pathogenic_max": 250.0, "locus_structure": "(X)*"}
    d.update(kw)
    return d


LOCI = [
    locus("HD_HTT", "HTT", "chr4", 100, 160, ["CAG"]),                                          # the catalog's frame is wider
    locus("CANVAS_RFC1", "RFC1", "chr4", 5000, 5060, ["AAGGG", "ACAGG", "AGAGG"], inheritance=["AR"],   # AGAGG is not in the catalog
          benign_min=0.0, benign_max=11.0, intermediate_min=11.0, intermediate_max=200.0, pathogenic_min=400.0),
    locus("CJD_PRNP", "PRNP", "chr20", 300, 420, ["CCTCATGG"], benign_min=4.0, benign_max=4.0, intermediate_min=None,
          intermediate_max=None, pathogenic_min=5.0, pathogenic_max=16.0),
    locus("PRTS_ARX", "ARX", "chrX", 700, 760, ["GCG"], inheritance=["XR"], benign_max=None, intermediate_min=None, intermediate_max=None),
    locus("HMNR7_VWA1", "VWA1", "chr1", 10, 30, ["GGCGCGGAGC"], inheritance=["AR"], benign_min=2.0, benign_max=2.0, pathogenic_min=1.0,
          pathogenic_max=3.0, intermediate_min=None, intermediate_max=None),
]


def load_module():
    spec = importlib.util.spec_from_file_location("make_strchive_map", SCRIPT)
    mod = importlib.util.module_from_spec(spec)   # type: ignore[arg-type]
    spec.loader.exec_module(mod)                  # type: ignore[union-attr]
    return mod


class StrchiveMapTest(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.catalog = self.tmp / "catalog.bed.gz"
        with gzip.open(self.catalog, "wt") as fh:
            fh.write(CATALOG)
        self.json = self.tmp / "STRchive-loci.json"

    def run_script(self, loci: list[dict[str, object]]) -> subprocess.CompletedProcess:
        self.json.write_text(json.dumps(loci))
        return subprocess.run([sys.executable, str(SCRIPT), "--strchive", str(self.json), "--strchive-tag", "v0.0.0", "--catalog", str(self.catalog),
                               "--out-loci", str(self.tmp / "loci.tsv"), "--out-overlaps", str(self.tmp / "overlaps.tsv")],
                              capture_output=True, text=True)

    def table(self, name: str) -> tuple[list[str], list[dict[str, str]]]:
        lines = (self.tmp / name).read_text().splitlines()
        head = [line for line in lines if line.startswith("##")]
        body = [line for line in lines if not line.startswith("##")]
        cols = body[0].split("\t")
        return head, [dict(zip(cols, line.split("\t"))) for line in body[1:]]

    def test_tables(self):
        res = self.run_script(LOCI)
        self.assertEqual(res.returncode, 0, res.stderr)
        head, rows = self.table("loci.tsv")
        self.assertTrue(any("CC BY 4.0" in h for h in head))
        self.assertTrue(any(h.startswith("##tr_map_input=strchive\tSTRchive-loci.json@v0.0.0\tsha256:") for h in head))
        by = {r["trid"]: r for r in rows}
        self.assertEqual([r["trid"] for r in rows], ["HD_HTT", "CANVAS_RFC1", "CJD_PRNP", "PRTS_ARX", "HMNR7_VWA1"])   # the catalog's contig order
        htt = by["HD_HTT"]
        self.assertEqual((htt["pathogenic_catalog_motifs"], htt["motif_match"], htt["frame_differs"], htt["twin_trids"]),
                         ("CAG", "exact", "yes", "chr4_90_260_GCC"))
        self.assertEqual((htt["benign_max"], htt["intermediate_min"], htt["pathogenic_min"], htt["motif_len"], htt["start"], htt["end"]),
                         ("26", "27", "40", "3", "100", "200"))   # the catalog's interval, integral copies without a decimal
        rfc1 = by["CANVAS_RFC1"]
        self.assertEqual((rfc1["pathogenic_catalog_motifs"], rfc1["motif_match"], rfc1["inheritance"]), ("AAGGG,ACAGG", "partial:exact", "AR"))
        self.assertNotIn("AAGAG", rfc1["pathogenic_catalog_motifs"])   # AGAGG and the benign AAGAG are different rotation classes
        prnp = by["CJD_PRNP"]
        self.assertEqual((prnp["pathogenic_catalog_motifs"], prnp["motif_match"]), ("CCTCATGG,GGCCTCAT", "rotation"))   # exact and a rotation of it
        self.assertEqual((prnp["intermediate_min"], prnp["frame_differs"]), ("", "no"))
        arx = by["PRTS_ARX"]
        self.assertEqual((arx["pathogenic_catalog_motifs"], arx["motif_match"], arx["benign_max"], arx["inheritance"]), ("NGC", "iupac", "", "XR"))
        self.assertEqual((by["HMNR7_VWA1"]["benign_max"], by["HMNR7_VWA1"]["pathogenic_min"]), ("2", "1"))   # the accepted nested range
        self.assertEqual({r["assessable"] for r in rows}, {"yes"})                # none of the fixture's loci is on the list
        self.assertTrue(any("not assessable by copy count: ADTKD_MUC1" in h for h in head))
        _h, pairs = self.table("overlaps.tsv")
        self.assertEqual(pairs, [{"ordinary_trid": "chr4_90_260_GCC", "named_trid": "HD_HTT", "overlap_bp": "100"}])
        self.assertIn("5 named loci (1 with a twin), 1 overlap pairs", res.stderr)

    def test_shipped_table_marks_muc1(self):
        """The table the pipeline ships: 73 loci, one of them not assessable by copy count, with its reason."""
        text = (REPO / "references" / "tr_disease_loci.strchive-2.4.4.tsv").read_text().splitlines()
        cols = next(l for l in text if not l.startswith("##")).split("\t")
        rows = [dict(zip(cols, l.split("\t"))) for l in text if not l.startswith("##")][1:]
        self.assertEqual(len(rows), 73)
        no = [r for r in rows if r["assessable"] == "no"]
        self.assertEqual([r["trid"] for r in no], ["ADTKD_MUC1"])
        self.assertIn("7-cytosine tract", no[0]["assessable_note"])
        self.assertTrue(all(r["assessable_note"] == "" for r in rows if r["assessable"] == "yes"))

    def test_refusals(self):
        res = self.run_script(LOCI[:-1])                                   # a named catalog locus without a STRchive entry
        self.assertNotEqual(res.returncode, 0)
        self.assertIn("in the catalog only ['HMNR7_VWA1']", res.stderr)
        res = self.run_script(LOCI + [locus("NEW_GENE", "GENE", "chr1", 1, 2, ["A"])])
        self.assertIn("in STRchive only ['NEW_GENE']", res.stderr)
        nested = [dict(entry, benign_max=45.0) if entry["id"] == "HD_HTT" else entry for entry in LOCI]
        res = self.run_script(nested)                                      # a nested range nobody has looked at
        self.assertNotEqual(res.returncode, 0)
        self.assertIn("HD_HTT: benign_max 45.0 reaches pathogenic_min 40.0", res.stderr)

    def test_motif_matching(self):
        m = load_module()
        self.assertEqual(m.match_kind("GGCCTG", "GGCCTG"), "exact")
        self.assertEqual(m.match_kind("GCCTGG", "GGCCTG"), "rotation")
        self.assertEqual(m.match_kind("GCG", "GCN"), "iupac")
        self.assertEqual(m.match_kind("CGG", "GCN"), "iupac")              # a rotation that is compatible
        self.assertIsNone(m.match_kind("AGAGG", "AAGAG"))                  # RFC1's pathogenic and benign motifs
        self.assertIsNone(m.match_kind("CAG", "CTG"))                      # the other strand is not matched
        self.assertIsNone(m.match_kind("CAG", "CAGG"))
        self.assertEqual(m.match_motifs(["GAA"], ["A", "GAA"]), (["GAA"], "exact"))
        self.assertEqual(m.match_motifs(["TTT"], ["A"]), ([], "none"))
        self.assertEqual((m.num(None), m.num(26.0), m.num(7.5)), ("", "26", "7.5"))


if __name__ == "__main__":
    unittest.main()
