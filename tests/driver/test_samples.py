import json
import tempfile
import unittest
from pathlib import Path

from ugc_wgw import samples
from ugc_wgw.db import DB
from ugc_wgw.log import Events
from ugc_wgw.util import UgcError

from .helpers import make_project, run_cli, write_tsv


class SamplesTest(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.cfg = make_project(self.tmp)
        self.db = DB(self.cfg.db_path)
        self.events = Events(self.db, self.cfg.logs_dir / "events.jsonl")

    def tearDown(self):
        self.db.close()

    def test_parse_and_add(self):
        tsv = write_tsv(self.tmp, [
            {"sample_id": "S1", "sex": "male", "hifi_reads": "s1a.bam,s1b.bam", "fail_reads": "s1f.bam", "batch": "b1"},
            {"sample_id": "S2", "sex": "", "hifi_reads": "s2.bam", "father_id": "S1"},
        ])
        n = samples.add_samples(self.db, self.events, tsv)
        self.assertEqual(n, 2)
        s1 = self.db.get_sample("S1")
        self.assertEqual(s1.sex, "MALE")
        self.assertEqual(len(s1.hifi_reads), 2)
        self.assertEqual(len(s1.fail_reads), 1)
        self.assertTrue(all(p.startswith("/") for p in s1.hifi_reads))
        self.assertEqual(s1.meta["batch"], "b1")
        s2 = self.db.get_sample("S2")
        self.assertIsNone(s2.sex)
        self.assertEqual(s2.father_id, "S1")
        self.assertEqual(s2.fail_reads, [])

    def test_guide_example_sheet(self):
        example = Path(__file__).resolve().parents[2] / "docs" / "guide" / "examples" / "samples.tsv"
        self.assertEqual(samples.add_samples(self.db, self.events, example, check_paths=False), 4)
        child = self.db.get_sample("HG002")
        self.assertEqual((child.father_id, child.mother_id, child.sex), ("HG003", "HG004", "MALE"))
        self.assertEqual(len(child.hifi_reads), 2)
        self.assertEqual(len(child.fail_reads), 1)
        self.assertEqual(child.meta["batch"], "2026-01")
        self.assertIsNone(self.db.get_sample("NA12878").sex)

    def test_remove_refusals_and_force(self):
        from ugc_wgw import cohorts
        from .helpers import run_cli, seed_success, write_ids
        tsv = write_tsv(self.tmp, [{"sample_id": "R1", "hifi_reads": "r1.bam"}, {"sample_id": "R2", "hifi_reads": "r2.bam"},
                                   {"sample_id": "R3", "hifi_reads": "r3.bam"}])
        samples.add_samples(self.db, self.events, tsv)
        cohorts.freeze(self.db, self.events, "CR", write_ids(self.tmp, ["R2"]))
        seed_success(self.cfg, self.db, "singleton", "sample", "R3", {"sample_id": "R3", "hifi_reads": ["/x"]})
        with self.assertRaisesRegex(UgcError, "unknown sample"):
            samples.remove_samples(self.db, self.events, ["nope"], force=False, results_dir=self.cfg.results_dir)
        with self.assertRaisesRegex(UgcError, "member of cohort"):
            samples.remove_samples(self.db, self.events, ["R2"], force=True, results_dir=self.cfg.results_dir)
        with self.assertRaisesRegex(UgcError, "use --force"):
            samples.remove_samples(self.db, self.events, ["R1", "R3"], force=False, results_dir=self.cfg.results_dir)
        self.assertTrue(self.db.sample_exists("R1"))  # nothing deleted when any id is refused
        out = samples.remove_samples(self.db, self.events, ["R1", "R3"], force=True, results_dir=self.cfg.results_dir)
        self.assertEqual([(o["sample_id"], o["runs_deleted"]) for o in out], [("R1", 0), ("R3", 1)])
        self.assertFalse(self.db.sample_exists("R3"))
        self.assertEqual(self.db.runs_for("sample", "R3"), [])
        self.assertTrue((self.cfg.results_dir / "samples" / "R3").exists())  # results on disk untouched
        code, _, err = run_cli(["--project", str(self.cfg.project_dir), "samples", "remove", "R2"])
        self.assertEqual(code, 1)
        self.assertIn("cohorts are immutable", err)

    def test_missing_path_and_no_check(self):
        tsv = self.tmp / "bad.tsv"
        tsv.write_text("sample_id\thifi_reads\nS9\t/nonexistent/x.bam\n")
        with self.assertRaises(UgcError):
            samples.add_samples(self.db, self.events, tsv)
        self.assertEqual(samples.add_samples(self.db, self.events, tsv, check_paths=False), 1)

    def test_duplicates_and_replace(self):
        tsv = write_tsv(self.tmp, [{"sample_id": "S1", "hifi_reads": "a.bam"}])
        samples.add_samples(self.db, self.events, tsv)
        with self.assertRaises(UgcError):
            samples.add_samples(self.db, self.events, tsv)
        tsv2 = write_tsv(self.tmp, [{"sample_id": "S1", "sex": "FEMALE", "hifi_reads": "b.bam"}], name="s2.tsv")
        samples.add_samples(self.db, self.events, tsv2, replace=True)
        self.assertEqual(self.db.get_sample("S1").sex, "FEMALE")
        dup = self.tmp / "dup.tsv"
        dup.write_text("sample_id\thifi_reads\nX\ta\nX\tb\n")
        with self.assertRaises(UgcError):
            samples.add_samples(self.db, self.events, dup, check_paths=False)

    def test_bad_values(self):
        for text in ("sample_id\thifi_reads\tsex\nS1\ta.bam\tother\n", "sample_id\thifi_reads\nbad id\ta.bam\n",
                     "sample_id\thifi_reads\nS1\t\n", "sample_id\nS1\n"):
            tsv = self.tmp / "v.tsv"
            tsv.write_text(text)
            with self.assertRaises(UgcError):
                samples.add_samples(self.db, self.events, tsv, check_paths=False)

    def test_list_via_cli(self):
        tsv = write_tsv(self.tmp, [{"sample_id": "S1", "hifi_reads": "a.bam"}])
        code, out, err = run_cli(["--project", str(self.cfg.project_dir), "samples", "add", str(tsv)])
        self.assertEqual(code, 0, err)
        code, out, _ = run_cli(["--project", str(self.cfg.project_dir), "samples", "list", "--json"])
        self.assertEqual(code, 0)
        rows = json.loads(out)
        self.assertEqual(rows, [{"sample_id": "S1", "sex": "", "singleton": "-"}])


if __name__ == "__main__":
    unittest.main()
