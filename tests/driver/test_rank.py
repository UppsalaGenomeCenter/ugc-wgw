import itertools
import tempfile
import unittest
from array import array
from pathlib import Path

from ugc_wgw import rank
from ugc_wgw.util import UgcError

from .helpers import HtmlChecker, make_project, run_cli, write_ids, write_tsv

HEADER = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{samples}\n"


def vcf_text(samples: list[str], records: list[tuple[str, dict[str, str]]], fmt: str = "GT") -> str:
    """records: (SVTYPE, {sample: genotype field}); a sample not named gets svx's `./.`."""
    lines = [HEADER.replace("{samples}", "\t".join(samples))]
    for i, (kind, gts) in enumerate(records):
        cols = "\t".join(gts.get(s, "./.") for s in samples)
        lines.append(f"chr1\t{1000 + i}\tsv{i}\tN\t<{kind}>\t.\tPASS\tEND={2000 + i};SVTYPE={kind};SVLEN=100\t{fmt}\t{cols}\n")
    return "".join(lines)


def matrix_of(members: list[str], carriers: list[list[int]]) -> rank.Matrix:
    return rank.Matrix(members=members, site_carriers=[array("H", c) for c in carriers])


class ReadMatrixTest(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def test_carriers_types_and_missing_columns(self):
        """A carrier has any non-reference allele; `./.`, `0/0` and `0|0` are not carriers whatever follows the GT;
        BND and CNV records are left out; a wanted member without a column is reported, an unwanted column ignored."""
        text = vcf_text(["A", "B", "C", "X"], [
            ("INS", {"A": "0/1:5", "B": "1|1:9", "C": "0/0:3", "X": "0/1:1"}),
            ("DEL", {"A": "./.:.", "B": "0|0:2", "C": "./1:4"}),
            ("BND", {"A": "0/1:5", "B": "0/1:5"}),
            ("CNV", {"C": "0/1:5"}),
            ("DUP", {"X": "0/1:5"}),          # carried only by a column nobody asked for: no site
            ("INV", {"B": "1:7"}),            # a haploid call
        ], fmt="GT:AD")
        p = self.tmp / "c.vcf"
        p.write_text(text)
        m = rank.read_sv_matrix(p, ["A", "B", "C", "D"])
        self.assertEqual(m.missing, ["D"])
        self.assertEqual(m.records, 4)                                    # INS, DEL, DUP, INV
        self.assertEqual([list(c) for c in m.site_carriers], [[0, 1], [2], [1]])
        self.assertEqual(m.by_type, {"INS": 1, "DEL": 1, "INV": 1})

    def test_gzip_and_a_format_that_does_not_start_with_gt(self):
        import gzip
        p = self.tmp / "c.vcf.gz"
        with gzip.open(p, "wt") as fh:
            fh.write(vcf_text(["A", "B"], [("INS", {"A": "0/1"})]))
        self.assertEqual([list(c) for c in rank.read_sv_matrix(p, ["A", "B"]).site_carriers], [[0]])
        bad = self.tmp / "bad.vcf"
        bad.write_text(vcf_text(["A"], [("INS", {"A": "5:0/1"})], fmt="AD:GT"))
        with self.assertRaisesRegex(UgcError, "does not start with GT"):
            rank.read_sv_matrix(bad, ["A"])
        none = self.tmp / "none.vcf"
        none.write_text(vcf_text(["X"], [("INS", {"X": "0/1"})]))
        with self.assertRaisesRegex(UgcError, "no #CHROM line with a column"):
            rank.read_sv_matrix(none, ["A"])


class GreedyTest(unittest.TestCase):
    def test_common_alleles_decide_and_ties_fall_to_all_alleles(self):
        """A carries three shared sites; B one shared site and five private ones. B has more alleles, A more common
        ones: A is first. Then the member adding the most uncovered shared sites, private ones breaking ties."""
        members = ["A", "B", "C", "D"]
        sites = [[0, 1], [0, 2], [0, 3],                     # shared: s0 A,B  s1 A,C  s2 A,D
                 [1], [1], [1], [1], [1],                   # B's five singletons
                 [2, 3], [3]]                               # s8 C,D shared; s9 D's singleton
        m = matrix_of(members, sites)
        klass = rank.site_classes(m, 0.01)
        self.assertEqual(list(klass), [0, 0, 0, 3, 3, 3, 3, 3, 0, 3])
        rk = rank.greedy_rank(m, klass, eligible=[0, 1, 2, 3], first=[], related={})
        # after A: s8 is the only uncovered shared site (C and D tie on it); D has one more allele of its own
        self.assertEqual([members[i] for i in rk.order], ["A", "D", "B", "C"])
        self.assertEqual(rk.new_common, [3, 1, 0, 0])
        self.assertEqual(rk.new_all, [3, 2, 5, 0])
        self.assertEqual(rk.totals, [4, 0, 0, 6])
        self.assertEqual(rk.cum[-1], [4, 0, 0, 6])
        self.assertEqual(list(rk.first_pos), [0, 0, 0, 2, 2, 2, 2, 2, 1, 1])

    def test_frequency_classes(self):
        """200 members: 10 carriers is 5%, 2 carriers 1%, 1 carrier a singleton; common_min moves the middle line."""
        members = [f"M{i}" for i in range(200)]
        m = matrix_of(members, [list(range(10)), list(range(9)), [0, 1], [0, 1, 2], [5]])
        self.assertEqual(list(rank.site_classes(m, 0.01)), [0, 1, 1, 1, 3])
        self.assertEqual(list(rank.site_classes(m, 0.012)), [0, 1, 2, 1, 3])

    def test_relatives_are_held_back_and_first_is_honoured(self):
        members = ["A", "B", "C", "D"]
        sites = [[0, 1], [0, 1], [0, 1], [1, 2], [2, 3], [3], [0, 2]]
        m = matrix_of(members, sites)
        klass = rank.site_classes(m, 0.01)
        related = {0: {1: 0.5}, 1: {0: 0.5}}
        rk = rank.greedy_rank(m, klass, eligible=[0, 1, 2, 3], first=[], related=related)
        order = [members[i] for i in rk.order]
        self.assertEqual(order[-1], "B")                                 # A and B tie; A is taken, B waits
        self.assertEqual(order[0], "A")
        self.assertEqual(rk.tier[1], "related")
        self.assertEqual(rk.related_to[1], (0, 0.5))
        self.assertEqual({rk.tier[i] for i in (0, 2, 3)}, {"unrelated"})
        # --first: D then B at the top, as given; A is then B's relative and waits
        rk = rank.greedy_rank(m, klass, eligible=[0, 1, 2, 3], first=[3, 1], related=related)
        self.assertEqual([members[i] for i in rk.order], ["D", "B", "C", "A"])
        self.assertEqual([rk.tier[i] for i in rk.order], ["first", "first", "unrelated", "related"])
        # two relatives both named first: both are placed, in the order given
        rk = rank.greedy_rank(m, klass, eligible=[0, 1, 2, 3], first=[0, 1], related=related)
        self.assertEqual([members[i] for i in rk.order[:2]], ["A", "B"])
        self.assertEqual([rk.tier[i] for i in rk.order[:2]], ["first", "first"])
        self.assertNotIn(1, rk.related_to)

    def test_expected_random_equals_the_average_over_every_subset(self):
        members = list("ABCDEF")
        sites = [[0], [0, 1], [1, 2, 3], [0, 1, 2, 3, 4, 5], [4, 5], [2]]
        n = len(members)
        by_count: dict[int, int] = {}
        for c in sites:
            by_count[len(c)] = by_count.get(len(c), 0) + 1
        exp = rank.expected_random([by_count], n)
        for k in range(1, n + 1):
            subsets = list(itertools.combinations(range(n), k))
            mean = sum(sum(1 for c in sites if set(c) & set(s)) for s in subsets) / len(subsets)
            self.assertAlmostEqual(exp[k - 1][0], mean, places=9, msg=f"k={k}")
        self.assertAlmostEqual(exp[-1][0], len(sites))

    def test_held_out_shares(self):
        """30 members: the first 20 carry sites 0-19 one each plus a shared one; member 20+i carries the shared site,
        one site of an earlier member and one of its own: two of its three alleles are already there."""
        n = 30
        members = [f"M{i}" for i in range(n)]
        sites = [[i] + ([20 + i] if i < 10 else []) for i in range(20)]      # site i: member i (and member 20+i)
        sites.append(list(range(n)))                                           # carried by everyone
        sites += [[20 + i] for i in range(10)]                                 # private to the late members
        m = matrix_of(members, sites)
        klass = rank.site_classes(m, 0.01)
        rk = rank.greedy_rank(m, klass, eligible=list(range(n)), first=list(range(20)), related={})
        rows = rank.held_out(rk, klass, checkpoints=(20, 25, 40))
        self.assertEqual([r["first"] for r in rows], [20, 25])                # 40 is beyond the ranking
        self.assertEqual(rows[0]["members_outside"], 10)
        self.assertAlmostEqual(rows[0]["all_median"], 2 / 3)
        self.assertAlmostEqual(rows[0]["all_min"], 2 / 3)
        self.assertAlmostEqual(rows[0]["common_median"], 1.0)                 # their shared alleles are all covered


class LinesChartTest(unittest.TestCase):
    def test_integer_axis_has_no_repeated_or_fractional_ticks(self):
        import re
        from ugc_wgw import summary_svg as sv
        svg = sv.lines([("a", "ln2", [(1, 0.5), (2, 0.8), (3, 0.9), (4, 1.0)]), ("b", "ln3", [(1, 0.4), (4, 1.0)])],
                       x_label="members", y_label="share", percent=True)
        ticks = re.findall(r'text-anchor="middle">([^<]+)</text>', svg)
        self.assertEqual(ticks, ["1", "2", "3", "4", "members"])
        self.assertEqual(svg.count("<polyline"), 2)
        self.assertIn("100%", svg)
        self.assertEqual(sv.lines([], percent=True), "<p class='muted'>No data.</p>")
        long = sv.lines([("a", "ln2", [(i, i / 2000) for i in range(1, 2001)])], percent=True)
        self.assertLess(long.count(","), 800)                             # thinned to a few hundred points


class EligibilityTest(unittest.TestCase):
    def test_read_based_fails_exclude_assembly_fails_do_not(self):
        def digest(flags, depth=30.0, n50=18000):
            return {"qc": {"status": "fail" if any(s == "fail" for _c, s in flags) else "pass",
                           "flags": [{"code": c, "severity": s} for c, s in flags]},
                    "coverage": {"mean": depth, "inferred_sex": "XX"}, "reads": {"read_length_n50": n50}, "identity": {"resolved": "XX"}}
        digests = {"A": digest([]), "B": digest([("low_depth", "fail"), ("asm_hap_length", "fail")]),
                   "C": digest([("asm_hap_length", "fail"), ("low_mapped", "warn")]), "D": digest([], depth=12.0), "E": digest([], n50=9000)}
        rows = {m.sample_id: m for m in rank.member_table(digests, ["A", "B", "C", "D", "E", "F"], min_depth=15, min_read_n50=10000)}
        self.assertEqual({k: v.eligible for k, v in rows.items()}, {"A": True, "B": False, "C": True, "D": False, "E": False, "F": False})
        self.assertEqual(rows["B"].reason, "QC fail: low_depth")
        self.assertEqual(rows["C"].qc_status, "fail")                     # shown as it is; the assembly flag does not decide here
        self.assertEqual(rows["D"].reason, "mean depth 12.0 below 15")
        self.assertEqual(rows["E"].reason, "read N50 9000 below 10000")
        self.assertEqual(rows["F"].reason, "no successful run to summarise")
        self.assertEqual((rows["A"].sex, rows["A"].depth_mean, rows["A"].read_n50), ("XX", 30.0, 18000.0))


class EmitTest(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        rows = [("1", "A", "first"), ("2", "B", "unrelated"), ("3", "C", "unrelated"), ("4", "D", "unrelated"),
                ("5", "E", "related"), ("", "F", "excluded")]
        text = "\t".join(rank.RANK_COLUMNS) + "\n"
        for r, sid, tier in rows:
            cells = {"rank": r, "sample_id": sid, "tier": tier}
            text += "\t".join(cells.get(c, "") for c in rank.RANK_COLUMNS) + "\n"
        self.path = self.tmp / "C.rank.tsv"
        self.path.write_text(text)

    def test_prefixes_are_nested_and_failures_are_replaced_from_the_ranking(self):
        two = rank.emit(self.path, 2, pass_list=None, out=self.tmp / "p2.txt")
        three = rank.emit(self.path, 3, pass_list=None, out=self.tmp / "p3.txt")
        self.assertEqual((two, three), (["A", "B"], ["A", "B", "C"]))
        self.assertEqual((self.tmp / "p3.txt").read_text(), "A\nB\nC\n")
        passing = write_ids(self.tmp, ["A", "C", "D", "E"])                  # B's assembly failed QC
        self.assertEqual(rank.emit(self.path, 3, pass_list=passing, out=self.tmp / "q3.txt"), ["A", "C", "D"])
        self.assertEqual(rank.emit(self.path, 2, pass_list=passing, out=self.tmp / "q2.txt"), ["A", "C"])

    def test_relatives_only_on_request_and_never_more_than_there_are(self):
        with self.assertRaisesRegex(UgcError, "only 4 ranked member.s. qualify, 5 asked for .*--include-related"):
            rank.emit(self.path, 5, pass_list=None, out=self.tmp / "x.txt")
        self.assertEqual(rank.emit(self.path, 5, pass_list=None, out=self.tmp / "x.txt", include_related=True), ["A", "B", "C", "D", "E"])
        self.assertEqual(rank.emit(self.path, 0, pass_list=None, out=self.tmp / "all.txt"), ["A", "B", "C", "D", "E"])   # the excluded member never
        with self.assertRaisesRegex(UgcError, "count of members"):
            rank.emit(self.path, -1, pass_list=None, out=self.tmp / "x.txt")
        other = self.tmp / "other.tsv"
        other.write_text("a\tb\n1\t2\n")
        with self.assertRaisesRegex(UgcError, "not a ranking"):
            rank.emit(other, 1, pass_list=None, out=self.tmp / "x.txt")


class RankEndToEnd(unittest.TestCase):
    """A four-member cohort through the fake engine (every member carries every fixture SV, S1 and S2 are related,
    S3 and S4 duplicates): the four files, the tiers, a nested emit and a freeze from it."""

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.cfg = make_project(self.tmp)
        self.proj = str(self.cfg.project_dir)
        tsv = write_tsv(self.tmp, [{"sample_id": f"S{i}", "sex": "MALE", "hifi_reads": f"s{i}.bam"} for i in range(1, 5)])
        self.assertEqual(run_cli(["--project", self.proj, "samples", "add", str(tsv)])[0], 0)
        self.assertEqual(run_cli(["--project", self.proj, "cohort", "freeze", "C1", "--samples", str(write_ids(self.tmp, ["S1", "S2", "S3", "S4"]))])[0], 0)
        code, _, err = run_cli(["--project", self.proj, "submit", "--mode", "standalone", "--cohort", "C1"],
                               env={"UGC_WGW_FAKE_RELATEDNESS": "S1:S2=0.52,S3:S4=0.97"})
        self.assertEqual(code, 0, err)
        self.out = self.cfg.results_dir / "reports" / "rank"

    def test_rank_emit_freeze(self):
        code, out, err = run_cli(["--project", self.proj, "cohort", "rank", "C1"])
        self.assertEqual(code, 0, err)
        self.assertEqual([Path(p).name for p in out.split()], ["C1.rank.tsv", "C1.rank.curves.tsv", "C1.rank.heldout.tsv", "C1.rank.html"])
        self.assertIn("ranked 4 of 4 members (0 placed first, 2 unrelated, 2 held back as relatives; 0 excluded)", err)
        rows = rank.read_rank(self.out / "C1.rank.tsv")
        self.assertEqual([(r["rank"], r["sample_id"], r["tier"], r["related_to"]) for r in rows],
                         [("1", "S1", "unrelated", ""), ("2", "S3", "unrelated", ""), ("3", "S2", "related", "S1"), ("4", "S4", "related", "S3")])
        self.assertEqual(rows[0]["cum_all"], "1.0000")                        # identical carriers: the first member holds everything
        self.assertEqual(rows[2]["relatedness"], "0.520")
        curves = (self.out / "C1.rank.curves.tsv").read_text().splitlines()
        self.assertEqual(curves[0].split("\t")[:2], ["n", "greedy_ge5"])
        self.assertEqual(len(curves), 5)
        page = (self.out / "C1.rank.html").read_text()
        checker = HtmlChecker()
        checker.feed(page)
        checker.close()
        self.assertEqual((checker.errors, checker.stack, checker.titles), ([], [], 1))
        self.assertIn("members ranked for a pangenome graph", page)
        self.assertIn("this ranking", page)
        # the first two as a cohort; three need the relatives
        code, out, err = run_cli(["--project", self.proj, "cohort", "rank", "C1", "--emit", "2"])
        self.assertEqual(code, 0, err)
        self.assertEqual(Path(out.strip()).read_text(), "S1\nS3\n")
        code, _, err = run_cli(["--project", self.proj, "cohort", "rank", "C1", "--emit", "3"])
        self.assertNotEqual(code, 0)
        self.assertIn("--include-related", err)
        code, out, err = run_cli(["--project", self.proj, "cohort", "rank", "C1", "--emit", "0"])
        self.assertEqual(Path(out.strip()).read_text(), "S1\nS3\nS2\nS4\n")
        self.assertEqual(run_cli(["--project", self.proj, "cohort", "freeze", "P2", "--samples", str(self.out / "C1.top2.txt")])[0], 0)

    def test_first_floors_and_refusals(self):
        first = self.tmp / "first.txt"
        first.write_text("S4\n")
        code, _, err = run_cli(["--project", self.proj, "cohort", "rank", "C1", "--first", str(first), "--out-dir", str(self.tmp / "r")])
        self.assertEqual(code, 0, err)
        rows = rank.read_rank(self.tmp / "r" / "C1.rank.tsv")
        self.assertEqual([(r["sample_id"], r["tier"]) for r in rows][:2], [("S4", "first"), ("S1", "unrelated")])
        self.assertEqual(rows[-1]["related_to"], "S4")                        # S3, the duplicate of the forced member
        code, _, err = run_cli(["--project", self.proj, "cohort", "rank", "C1", "--min-depth", "1000", "--out-dir", str(self.tmp / "r2")])
        self.assertNotEqual(code, 0)
        self.assertIn("no eligible member", err)
        code, _, err = run_cli(["--project", self.proj, "cohort", "rank", "NOPE"])
        self.assertNotEqual(code, 0)
        self.assertIn("unknown cohort NOPE", err)
        code, _, err = run_cli(["--project", self.proj, "cohort", "rank", "C1", "--mode", "joint"])
        self.assertNotEqual(code, 0)
        self.assertIn("no successful cohort_merge run in joint mode", err)
        code, _, err = run_cli(["--project", self.proj, "cohort", "rank", "C1", "--emit", "2", "--out-dir", str(self.tmp / "empty")])
        self.assertNotEqual(code, 0)
        self.assertIn("no saved ranking", err)


if __name__ == "__main__":
    unittest.main()
