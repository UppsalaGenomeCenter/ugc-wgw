"""Tests for scripts/manifest-to-samples.py: a per-file manifest becomes a sheet that `ugc-wgw samples add` accepts."""
from __future__ import annotations

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from ugc_wgw.db import DB

from .helpers import REPO, make_project, run_cli

SCRIPT = REPO / "scripts" / "manifest-to-samples.py"
HEADER = "project,sample,file,family_id,paternal_id,maternal_id,sex,phenotype\n"


def convert(*args: str) -> tuple[int, str, str]:
    res = subprocess.run([sys.executable, str(SCRIPT), *args], capture_output=True, text=True)
    return res.returncode, res.stdout, res.stderr


class ManifestTest(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def write(self, text: str, name: str = "manifest.csv") -> Path:
        p = self.tmp / name
        p.write_text(text)
        return p

    def test_trio_with_fail_reads_and_metadata(self):
        m = self.write(HEADER
                       + "p1,K1,/d/K1/r1.hifi_reads.bc1.bam,F1,D1,M1,1,2\n"
                       + "p1,K1,/d/K1/r2.hifi_reads.bc1.bam,F1,D1,M1,1,2\n"
                       + "p1,K1,/d/K1/r1.fail_reads.bc1.bam,F1,D1,M1,1,2\n"
                       + "p1,K1,/d/K1/r1.hifi_reads.bc1.bam,F1,D1,M1,1,2\n"   # duplicate: dropped with a warning
                       + "p1,D1,/d/D1/r1.hifi_reads.bc2.bam,F1,0,0,1,1\n"
                       + "p1,M1,/d/M1/r1.hifi_reads.bc3.bam,F1,NA,NA,2,1\n"
                       + "p1,X1,/d/X1/r1.hifi_reads.bc4.bam,F2,0,0,0,0\n")
        code, out, err = convert(str(m))
        self.assertEqual(code, 0, err)
        rows = [line.split("\t") for line in out.splitlines()]
        self.assertEqual(rows[0], ["sample_id", "sex", "hifi_reads", "fail_reads", "father_id", "mother_id", "project", "family_id", "phenotype"])
        by = {r[0]: r for r in rows[1:]}
        self.assertEqual(by["K1"], ["K1", "MALE", "/d/K1/r1.hifi_reads.bc1.bam,/d/K1/r2.hifi_reads.bc1.bam", "/d/K1/r1.fail_reads.bc1.bam",
                                    "D1", "M1", "p1", "F1", "2"])
        self.assertEqual(by["M1"][1:6], ["FEMALE", "/d/M1/r1.hifi_reads.bc3.bam", "", "", ""])
        self.assertEqual(by["X1"][1], "")   # PLINK 0 = unknown sex
        self.assertIn("duplicate file dropped", err)
        self.assertIn("4 sample(s), 6 file(s), sex known for 3, parents for 1", err)
        # the sheet is accepted by the driver, parents and metadata included
        sheet = self.tmp / "samples.tsv"
        sheet.write_text(out)
        cfg = make_project(self.tmp)
        code, _, err = run_cli(["--project", str(cfg.project_dir), "samples", "add", str(sheet), "--no-check"])
        self.assertEqual(code, 0, err)
        db = DB(cfg.db_path)
        k1 = db.get_sample("K1")
        db.close()
        self.assertEqual((k1.sex, k1.father_id, k1.mother_id, len(k1.hifi_reads), len(k1.fail_reads)), ("MALE", "D1", "M1", 2, 1))
        self.assertEqual((k1.meta["family_id"], k1.meta["phenotype"], k1.meta["project"]), ("F1", "2", "p1"))

    def test_tsv_input_aliases_project_filter_and_check(self):
        m = self.write("sample_id\thifi_reads\tsex\tfather_id\tmother_id\tproject\n"
                       "A\t/d/A.bam\tF\t\t\tp1\nB\t/d/B.bam\tMALE\t\t\tp2\n", "manifest.tsv")
        code, out, err = convert(str(m), "--project", "p1")
        self.assertEqual(code, 0, err)
        self.assertEqual([line.split("\t")[0] for line in out.splitlines()], ["sample_id", "A"])
        self.assertIn("FEMALE", out)
        code, out, err = convert(str(m), "--check")
        self.assertEqual(code, 1)
        self.assertIn("file not found: /d/A.bam", err)
        real = self.write(HEADER + f"p1,S,{self.tmp / 'manifest.tsv'},F,0,0,1,0\n", "real.csv")
        self.assertEqual(convert(str(real), "--check")[0], 0)

    def test_errors_are_listed(self):
        m = self.write(HEADER
                       + "p1,S1,/d/a.bam,F1,0,0,1,0\n"
                       + "p1,S1,/d/b.bam,F1,0,0,2,0\n"        # sex conflict
                       + "p1,S2,/d/c.bam,F1,S2,0,3,0\n"        # bad sex code, own father
                       + "p1,S3,/d/d.fail_reads.bam,F1,0,0,1,0\n")   # fail reads only
        code, out, err = convert(str(m))
        self.assertEqual(code, 1)
        for text in ("sex 'FEMALE' differs from 'MALE'", "sex '3' is not", "is its own father_id", "S3: no hifi_reads file"):
            self.assertIn(text, err)
        self.assertEqual(out, "")
        m = self.write("project,name,path\np1,S1,/d/a.bam\n")
        code, out, err = convert(str(m))
        self.assertEqual(code, 1)
        self.assertIn("no column for the sample id", err)
        m = self.write(HEADER + "p1,S1,/d/a.bam,F1,0,0,1,0\n")
        code, out, err = convert(str(m), "--project", "nope")
        self.assertEqual(code, 1)
        self.assertIn("no samples of project 'nope'", err)


if __name__ == "__main__":
    unittest.main()
