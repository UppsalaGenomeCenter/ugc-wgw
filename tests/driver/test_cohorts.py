import hashlib
import tempfile
import unittest
from pathlib import Path

from ugc_wgw import cohorts, samples
from ugc_wgw.db import DB
from ugc_wgw.log import Events
from ugc_wgw.util import UgcError

from .helpers import make_project, write_ids, write_tsv


class CohortsTest(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.cfg = make_project(self.tmp)
        self.db = DB(self.cfg.db_path)
        self.events = Events(self.db, self.cfg.logs_dir / "events.jsonl")
        samples.add_samples(self.db, self.events, write_tsv(self.tmp, [
            {"sample_id": s, "hifi_reads": f"{s}.bam"} for s in ("S1", "S2", "S3")]))

    def tearDown(self):
        self.db.close()

    def test_freeze_order_sha_immutable(self):
        ids = write_ids(self.tmp, ["S3", "S1", "# comment", "", "S2"])
        cohort = cohorts.freeze(self.db, self.events, "C1", ids)
        self.assertEqual(cohort.members, ["S3", "S1", "S2"])
        self.assertEqual(cohort.sample_list_sha256, hashlib.sha256(b"S3\nS1\nS2\n").hexdigest())
        with self.assertRaises(UgcError):
            cohorts.freeze(self.db, self.events, "C1", ids)
        # adding a sample later does not change membership
        samples.add_samples(self.db, self.events, write_tsv(self.tmp, [{"sample_id": "S4", "hifi_reads": "s4.bam"}], "s4.tsv"))
        self.assertEqual(self.db.get_cohort("C1").members, ["S3", "S1", "S2"])

    def test_guide_example_list(self):
        example = Path(__file__).resolve().parents[2] / "docs" / "guide" / "examples" / "cohort.txt"
        self.assertEqual(cohorts.read_ids(example), ["HG002", "HG003", "HG004", "NA12878"])

    def test_unknown_and_duplicate_ids(self):
        with self.assertRaises(UgcError):
            cohorts.freeze(self.db, self.events, "C2", write_ids(self.tmp, ["S1", "NOPE"]))
        with self.assertRaises(UgcError):
            cohorts.freeze(self.db, self.events, "C3", write_ids(self.tmp, ["S1", "S1"]))
        self.assertIsNone(self.db.get_cohort("C2"))


if __name__ == "__main__":
    unittest.main()
