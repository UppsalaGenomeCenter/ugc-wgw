import json
import tempfile
import unittest
from pathlib import Path

from ugc_wgw import samples
from ugc_wgw.db import DB
from ugc_wgw.log import Events
from ugc_wgw.util import UgcError

from .helpers import make_project, run_cli, write_bam, write_tsv


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

    def warnings(self) -> list[str]:
        return [str(json.loads(l)["detail"]["message"]) for l in (self.cfg.logs_dir / "events.jsonl").read_text().splitlines()
                if json.loads(l)["event"] == "sample.warning"]

    def test_inspection_records_yield_and_warns(self):
        write_bam(self.tmp / "data" / "big.bam", reads=1200, read_len=15_000, movie="m84000_240101_000000_s2", pbi_reads=1200)
        tsv = write_tsv(self.tmp, [{"sample_id": "S1", "hifi_reads": "s1a.bam,big.bam", "fail_reads": "s1f.bam"}])
        self.assertEqual(samples.add_samples(self.db, self.events, tsv), 1)
        rec = self.db.get_sample("S1")
        small, big = rec.hifi_reads
        self.assertEqual((rec.input_info[small]["reads"], rec.input_info[small]["reads_exact"]), (10, True))
        self.assertEqual((rec.input_info[big]["reads"], rec.input_info[big]["pbi"], rec.input_info[big]["bytes"]),
                         (1200, True, Path(big).stat().st_size))
        self.assertEqual(rec.input_info[big]["movies"], ["m84000_240101_000000_s2"])
        self.assertIn(rec.fail_reads[0], rec.input_info)
        check = rec.meta["input_check"]
        self.assertEqual((check["hifi_files"], check["fail_files"], check["reads"], check["exact"], check["dropped"]),
                         (2, 1, 1210, True, []))
        self.assertAlmostEqual(check["gbases"], 0.018, places=3)
        warnings = self.warnings()
        self.assertTrue(any("s1a.bam: only 10 reads (10 kb), below file_reads_min 1000" in w for w in warnings), warnings)
        self.assertFalse(any("big.bam" in w for w in warnings), warnings)      # 1200 reads is above the floor
        self.assertTrue(any("S1: 18.0 Mb of HiFi bases in 2 file(s), about 0.0x of GRCh38, below sample_gbases_min 30" in w
                            for w in warnings), warnings)
        added = [json.loads(l) for l in (self.cfg.logs_dir / "events.jsonl").read_text().splitlines()
                 if json.loads(l)["event"] == "sample.added"][0]["detail"]
        self.assertEqual((added["reads"], added["gbases"]), (1210, 0.018))

    def test_thresholds(self):
        self.assertEqual(samples.input_thresholds(None), samples.INPUT_THRESHOLDS)
        self.assertEqual(samples.input_thresholds({"file_reads_min": 0})["file_reads_min"], 0.0)
        for bad in ({"nope": 1}, {"file_reads_min": "x"}, {"sample_gbases_min": -1}, {"file_reads_min": True}):
            with self.assertRaises(UgcError):
                samples.input_thresholds(bad)
        tsv = write_tsv(self.tmp, [{"sample_id": "S1", "hifi_reads": "a.bam"}])
        samples.add_samples(self.db, self.events, tsv, thresholds={"file_reads_min": 0, "sample_gbases_min": 0})
        self.assertEqual(self.warnings(), [])
        # the project's config carries them
        (self.tmp / "p2").mkdir()
        (self.tmp / "p3").mkdir()
        cfg2 = make_project(self.tmp / "p2", input_thresholds={"file_reads_min": 5})
        code, _, err = run_cli(["--project", str(cfg2.project_dir), "samples", "add", str(tsv)])
        self.assertEqual(code, 0, err)
        self.assertIn("below sample_gbases_min 30", err)
        self.assertNotIn("file_reads_min", err)
        cfg3 = make_project(self.tmp / "p3", input_thresholds={"typo": 5})
        code, _, err = run_cli(["--project", str(cfg3.project_dir), "samples", "add", str(tsv)])
        self.assertEqual(code, 1)
        self.assertIn("input_thresholds: unknown key 'typo'", err)

    def test_empty_truncated_and_readless_files_are_refused(self):
        write_bam(self.tmp / "data" / "e.bam", empty=True)
        write_bam(self.tmp / "data" / "t.bam", reads=300, truncated=True)
        write_bam(self.tmp / "data" / "h.bam", header_only=True)
        tsv = write_tsv(self.tmp, [{"sample_id": "S1", "hifi_reads": "ok.bam,e.bam"}, {"sample_id": "S2", "hifi_reads": "t.bam"},
                                   {"sample_id": "S3", "hifi_reads": "ok2.bam", "fail_reads": "h.bam"}])
        with self.assertRaises(UgcError) as ctx:
            samples.add_samples(self.db, self.events, tsv)
        msg = str(ctx.exception)
        self.assertIn("S1: hifi_reads", msg)
        self.assertIn("e.bam: empty file (0 bytes)", msg)
        self.assertIn("t.bam: truncated: the BGZF end-of-file marker is missing", msg)
        self.assertIn("S3: fail_reads", msg)
        self.assertIn("h.bam: no reads (header only)", msg)
        self.assertEqual(self.db.list_sample_ids(), [])
        # --no-inspect: existence only, nothing recorded about the files
        self.assertEqual(samples.add_samples(self.db, self.events, tsv, inspect=False), 3)
        rec = self.db.get_sample("S1")
        self.assertEqual((len(rec.hifi_reads), rec.input_info, rec.meta.get("input_check")), (2, {}, None))

    def test_drop_empty(self):
        write_bam(self.tmp / "data" / "e.bam", empty=True)
        write_bam(self.tmp / "data" / "h.bam", header_only=True)
        write_bam(self.tmp / "data" / "t.bam", reads=300, truncated=True)
        tsv = write_tsv(self.tmp, [{"sample_id": "S1", "hifi_reads": "a.bam,e.bam,h.bam", "fail_reads": "h.bam"},
                                   {"sample_id": "S2", "hifi_reads": "e.bam"}])
        with self.assertRaisesRegex(UgcError, "S2: no hifi_reads left after dropping empty files"):
            samples.add_samples(self.db, self.events, tsv, drop_empty=True)
        tsv = write_tsv(self.tmp, [{"sample_id": "S1", "hifi_reads": "a.bam,e.bam,h.bam", "fail_reads": "h.bam"}], "one.tsv")
        self.assertEqual(samples.add_samples(self.db, self.events, tsv, drop_empty=True), 1)
        rec = self.db.get_sample("S1")
        self.assertEqual(([Path(p).name for p in rec.hifi_reads], rec.fail_reads), (["a.bam"], []))
        self.assertEqual([d["reason"] for d in rec.meta["input_check"]["dropped"]],
                         ["empty file (0 bytes)", "no reads (header only)", "no reads (header only)"])
        dropped = [json.loads(l) for l in (self.cfg.logs_dir / "events.jsonl").read_text().splitlines()
                   if json.loads(l)["event"] == "sample.input_dropped"]
        self.assertEqual([(d["detail"]["kind"], Path(d["detail"]["path"]).name) for d in dropped],
                         [("hifi_reads", "e.bam"), ("hifi_reads", "h.bam"), ("fail_reads", "h.bam")])
        self.assertTrue(any("dropped hifi_reads e.bam" in w for w in self.warnings()))
        # damage is never dropped
        tsv = write_tsv(self.tmp, [{"sample_id": "S3", "hifi_reads": "a.bam,t.bam"}], "t.tsv")
        with self.assertRaisesRegex(UgcError, "t.bam: truncated"):
            samples.add_samples(self.db, self.events, tsv, drop_empty=True)

    def test_check_cli(self):
        proj = str(self.cfg.project_dir)
        write_bam(self.tmp / "data" / "big.bam", reads=1500, read_len=12_000, movie="mA", pbi_reads=1500)
        tsv = write_tsv(self.tmp, [{"sample_id": "S1", "hifi_reads": "big.bam,small.bam"}, {"sample_id": "S2", "hifi_reads": "s2.bam"}])
        code, _, err = run_cli(["--project", proj, "samples", "add", str(tsv)])
        self.assertEqual(code, 0, err)
        code, out, err = run_cli(["--project", proj, "samples", "check", "--json"])
        self.assertEqual(code, 0, err)
        rows = json.loads(out)
        self.assertEqual([(r["sample_id"], r["file"]) for r in rows], [("S1", "big.bam"), ("S1", "small.bam"), ("S2", "s2.bam")])
        self.assertEqual(rows[0]["status"], "ok")
        self.assertEqual((rows[0]["reads"], rows[0]["movie"], rows[0]["bases"]), ("1,500", "mA", "18.0 Mb"))
        self.assertTrue(rows[1]["status"].startswith("warning: only 10 reads"), rows[1])
        self.assertIn("3 file(s) of 2 sample(s): 0 problem(s)", err)
        self.assertIn("warning: S2: 10 kb of HiFi bases", err)
        code, out, err = run_cli(["--project", proj, "samples", "check", "S2", "--tsv"])
        self.assertEqual(code, 0, err)
        self.assertEqual(out.splitlines()[0], "sample_id\tkind\tfile\tsize\treads\tbases\tmovie\tstatus")
        self.assertEqual(len(out.splitlines()), 2)
        # the data changes under the project: check sees it, --stored shows what registration saw
        small = self.db.get_sample("S1").hifi_reads[1]
        write_bam(Path(small), reads=300, truncated=True)
        code, out, err = run_cli(["--project", proj, "samples", "check", "S1", "--json"])
        self.assertEqual(code, 1)
        self.assertTrue(json.loads(out)[1]["status"].startswith("problem: truncated"))
        self.assertIn("problem: S1: hifi_reads", err)
        self.assertEqual(self.db.get_sample("S1").input_info[small]["reads"], 10)   # a problem does not overwrite the record
        code, out, err = run_cli(["--project", proj, "samples", "check", "--stored", "--json"])
        self.assertEqual(code, 0, err)
        self.assertTrue(json.loads(out)[1]["status"].startswith("warning: only 10 reads"))
        code, _, err = run_cli(["--project", proj, "samples", "check", "nope"])
        self.assertEqual(code, 1)
        self.assertIn("unknown sample(s): nope", err)
        # a sample registered with --no-inspect has nothing stored
        tsv2 = write_tsv(self.tmp, [{"sample_id": "S3", "hifi_reads": "s3.bam"}], "s3.tsv")
        run_cli(["--project", proj, "samples", "add", str(tsv2), "--no-inspect"])
        code, out, _ = run_cli(["--project", proj, "samples", "check", "S3", "--stored", "--json"])
        self.assertEqual(json.loads(out), [{"sample_id": "S3", "kind": "hifi_reads", "file": "s3.bam", "status": "not inspected"}])
        code, out, _ = run_cli(["--project", proj, "samples", "check", "S3", "--json"])   # and check fills it in
        self.assertEqual(code, 0)
        self.assertEqual(self.db.get_sample("S3").input_info[self.db.get_sample("S3").hifi_reads[0]]["reads"], 10)

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
