"""cohort_graph_gt: the prerequisite in another cohort (the planner), the builder, the scripts embedded in
workflows/ugc_wgw/graph_gt/cohort.wdl (the site list, one shard's join, the gather), and the stage end to end with
the fake engine: the cohort page's section and the export."""
from __future__ import annotations

import gzip
import json
import os
import shutil
import stat
import subprocess
import sys
import tempfile
import textwrap
import unittest
from pathlib import Path

from ugc_wgw import cohorts, config, inputs, plan, samples
from ugc_wgw.db import DB
from ugc_wgw.log import Events

from .helpers import REPO, VERSION, HtmlChecker, make_project, run_cli, seed_success, write_ids, write_tsv

WDL = REPO / "workflows" / "ugc_wgw" / "graph_gt" / "cohort.wdl"
MEMBERS = ["S1", "S2", "S3"]

# Stand-ins for the three bcftools calls and the tabix call the join script makes, on plain gzip files: the test is
# of the script's arithmetic, not of htslib.
SHIM = r'''#!/usr/bin/env python3
import gzip, os, sys
tool, args = os.path.basename(sys.argv[0]), sys.argv[1:]
def regions(path):
    return [(f[0], int(f[1]), int(f[2])) for f in (l.split() for l in open(path)) if len(f) >= 3]
def hit(regs, chrom, pos):
    return any(c == chrom and s < pos <= e for c, s, e in regs)
def lines(path):
    with gzip.open(path, "rt") as fh:
        return [l.rstrip("\n") for l in fh]
if tool == "tabix":
    regs = regions(args[args.index("-R") + 1])
    for l in lines(args[-1]):
        f = l.split("\t")
        if not l.startswith("#") and hit(regs, f[0], int(f[1])):
            print(l)
elif args[0] == "--version":
    print("bcftools shim")
elif args[0] == "view" and "-h" in args:
    print("\n".join(l for l in lines(args[-1]) if l.startswith("#")))
elif args[0] == "view":
    regs = regions(args[args.index("-R") + 1])
    for l in lines(args[-1]):
        f = l.split("\t")
        if not l.startswith("#") and hit(regs, f[0], int(f[1])):
            print(l)
elif args[0] == "query":
    regs = regions(args[args.index("-R") + 1])
    for l in lines(args[-1]):
        f = l.split("\t")
        if l.startswith("#") or not hit(regs, f[0], int(f[1])):
            continue
        keys = f[8].split(":")
        vals = dict(zip(keys, f[9].split(":")))
        print("\t".join([f[2]] + [vals.get(k, ".") for k in ("GT", "GQ", "DP", "AD")]))
else:
    sys.exit("shim: unexpected call " + " ".join(sys.argv))
'''


def script(name: str) -> str:
    text = WDL.read_text()
    begin, end = f"# ---- {name} (begin) ----\n", f"# ---- {name} (end) ----\n"
    return textwrap.dedent(text.split(begin, 1)[1].split(end, 1)[0])


def read_tsv(path: Path) -> list[dict[str, str]]:
    lines = path.read_text().splitlines()
    cols = lines[0].split("\t")
    return [dict(zip(cols, line.split("\t"))) for line in lines[1:]]


def gz(path: Path, text: str) -> Path:
    with gzip.open(path, "wt") as fh:
        fh.write(text)
    return path


class ScriptsTest(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        shim = self.tmp / "bin"
        shim.mkdir()
        for tool in ("bcftools", "tabix"):
            (shim / tool).write_text(SHIM)
            (shim / tool).chmod((shim / tool).stat().st_mode | stat.S_IXUSR)
        self.env = dict(os.environ, PATH=f"{shim}{os.pathsep}{os.environ['PATH']}")

    def run_script(self, name: str, args: list[str]) -> subprocess.CompletedProcess:
        res = subprocess.run([sys.executable, "-c", script(name), *args], capture_output=True, text=True, cwd=self.tmp, env=self.env)
        self.assertEqual(res.returncode, 0, res.stderr)
        return res

    def test_site_list(self):
        ins, big = "A" + "C" * 60, "A" + "G" * 12000
        records = [
            ("chr1", 100, "A", ins, ["0|1", "0|0"]),
            ("chr1", 100, "A", ins, ["0|0", "1|."]),          # the same allele from a second record: one site, genotypes added
            ("chr1", 200, "T" * 1501, "T", ["1|1", "0|1"]),
            ("chr1", 300, "A", big, ["0|1", "0|0"]),          # longer than size_max: listed, marked, not genotyped
            ("chr1", 400, "A" * 80, "C" * 80, ["0|1", "0|0"]),  # long on both sides, no length change
            ("chr1", 500, "A", "AT", ["0|1", "0|0"]),         # a small variant
            ("chr1", 600, "A", "<INV>", ["0|1", "0|0"]),      # symbolic
        ]
        vcf = self.tmp / "norm.vcf"
        vcf.write_text("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000000>\n##INFO=<ID=LV,Number=1,Type=Integer,Description=\"x\">\n"
                       "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tZ\n"
                       + "".join(f"{c}\t{p}\t.\t{r}\t{a}\t60\t.\tLV=0\tGT\t" + "\t".join(g) + "\n" for c, p, r, a, g in records))
        # the two region lists: a tandem repeat over the first site (gzipped, as the TRGT catalogue is), an excluded region over the second
        gz(self.tmp / "tr.bed.gz", "chr1\t50\t100\tID=x;MOTIFS=A;STRUC=<TR>\n" "chr1\t90\t95\ty\n" "chr2\t0\t999\tz\n")
        (self.tmp / "excl.bed").write_text("chr1\t199\t200\tgap\n" "chr1\t100\t101\tgap\n")
        self.run_script("graph_gt_sites.py", [str(vcf), "50", "10000", "out", "tr.bed.gz", "excl.bed"])
        sites = [l.split("\t") for l in (self.tmp / "out.sites.vcf").read_text().splitlines() if not l.startswith("#")]
        self.assertEqual([(s[1], s[2].split("-")[2], s[2].split("-")[3]) for s in sites], [("100", "INS", "60"), ("200", "DEL", "1500"), ("300", "INS", "12000")])
        self.assertEqual(len({s[2] for s in sites}), 3)
        self.assertTrue(all(len(s) == 8 and s[6] == "PASS" for s in sites))
        info = {s[1]: dict(kv.split("=") if "=" in kv else (kv, True) for kv in s[7].split(";")) for s in sites}
        # position 100 is the last base of the repeat (0-based 50 to 100, a longer region before a shorter one); 101 is outside
        self.assertEqual([("IN_TR" in info[p], "IN_EXCL" in info[p]) for p in ("100", "200", "300")], [(True, False), (False, True), (False, False)])
        # the duplicate's carrier is added: A 0|1, Z 1|. -> two alternate alleles; Z's called alleles stay two
        self.assertEqual((info["100"]["GRAPH_AC"], info["100"]["GRAPH_AN"], info["100"]["GRAPH_NS"], info["100"]["SVLEN"]), ("2", "4", "2", "60"))
        self.assertEqual((info["200"]["GRAPH_AC"], info["200"]["SVTYPE"], info["200"]["SVLEN"]), ("3", "DEL", "-1500"))
        self.assertEqual(info["300"].get("SKIP"), "too_long")
        self.assertNotIn("SKIP", info["100"])
        header = [l for l in (self.tmp / "out.sites.vcf").read_text().splitlines() if l.startswith("#")]
        self.assertIn("##contig=<ID=chr1,length=1000000>", header)
        self.assertEqual(header[-1], "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
        self.assertFalse(any("ID=LV" in l for l in header))
        calls = (self.tmp / "out.graph_calls.tsv").read_text().splitlines()
        self.assertEqual(calls[0], "#chrom\tpos\tid\tA\tZ")
        self.assertEqual([l.split("\t")[3:] for l in calls[1:]], [["1:2", "1:2"], ["2:2", "1:2"], ["1:2", "0:2"]])
        skipped = read_tsv(self.tmp / "out.skipped.tsv")
        self.assertEqual([(r["pos"], r["reason"]) for r in skipped], [("300", "too_long"), ("400", "no_length_change"), ("600", "symbolic")])
        self.assertEqual(skipped[0]["id"], sites[2][2])
        counts = {(r["metric"], r["key"]): int(r["value"]) for r in read_tsv(self.tmp / "out.sites_counts.tsv")}
        self.assertEqual(counts[("records", "read")], 7)
        self.assertEqual((counts[("sites", "INS:50-999")], counts[("sites", "DEL:1k-9.9k")], counts[("sites", "INS:10k+")]), (1, 1, 1))
        self.assertEqual((counts[("left_out", "duplicate_merged")], counts[("left_out", "below_size_min")], counts[("left_out", "no_length_change")],
                          counts[("left_out", "symbolic")], counts[("not_genotyped", "too_long")]), (1, 1, 1, 1, 1))
        self.assertEqual((counts[("graph_sample", "A")], counts[("graph_sample", "Z")]), (1, 1))
        self.assertEqual((counts[("sites_in_tandem_repeats", "INS")], counts[("sites_in_excluded_regions", "DEL")]), (1, 1))
        # without the lists nothing is marked
        self.run_script("graph_gt_sites.py", [str(vcf), "50", "10000", "plain"])
        self.assertNotIn("IN_", "".join(l for l in (self.tmp / "plain.sites.vcf").read_text().splitlines() if not l.startswith("#")))

    def join_inputs(self) -> list[str]:
        head = "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000000>\n##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"x\">\n"
        sites = [("chr1", 100, "s1", "SVTYPE=INS;SVLEN=60;GRAPH_AC=1;GRAPH_AN=4;GRAPH_NS=2;IN_TR"),
                 ("chr1", 200, "s2", "SVTYPE=DEL;SVLEN=-1500;GRAPH_AC=1;GRAPH_AN=4;GRAPH_NS=2;IN_EXCL;IN_TR"),
                 ("chr1", 300, "s3", "SVTYPE=INS;SVLEN=12000;GRAPH_AC=1;GRAPH_AN=4;GRAPH_NS=2;SKIP=too_long"),
                 ("chr1", 400, "s4", "SVTYPE=INS;SVLEN=70;GRAPH_AC=2;GRAPH_AN=4;GRAPH_NS=2"),
                 ("chr1", 1500, "s5", "SVTYPE=INS;SVLEN=80;GRAPH_AC=1;GRAPH_AN=4;GRAPH_NS=2")]   # starts outside the shard
        gz(self.tmp / "sites.vcf.gz", head + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
           + "".join(f"{c}\t{p}\t{i}\tA\tAC\t.\tPASS\t{info}\n" for c, p, i, info in sites))
        gz(self.tmp / "graph_calls.tsv.gz", "#chrom\tpos\tid\tA\tZ\n" "chr1\t100\ts1\t1:2\t0:2\n" "chr1\t200\ts2\t1:2\t0:2\n"
           "chr1\t300\ts3\t1:2\t0:2\n" "chr1\t400\ts4\t1:1\t1:2\n" "chr1\t1500\ts5\t1:2\t0:2\n")
        (self.tmp / "shard.bed").write_text("chr1\t0\t1000\n")
        members = self.tmp / "members"
        shutil.rmtree(members, ignore_errors=True)   # the script renames the bundle's entries: a second run starts from a fresh bundle
        members.mkdir()
        fmt = "GT:FT:SQ:GQ:PS:DP:AD:KS"

        def member(sample: str, calls: dict[str, str]) -> str:
            return (head + f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n"
                    + "".join(f"{c}\t{p}\t{i}\tA\tAC\t.\tPASS\t{info}\t{fmt}\t{calls[i]}\n" for c, p, i, info in sites))

        gz(members / "000000__A.gt.vcf.gz", member("A", {"s1": "0|1:0:60:30:7:20:10,10:90", "s2": "1|1:0:60:40:7:18:0,18:95",
                                                          "s3": "./.:0:0:0:.:0:0,0:0", "s4": "0|1:0:60:30:7:20:10,10:90",
                                                          "s5": "0|1:0:60:30:7:20:10,10:90"}))
        gz(members / "000001__B.gt.vcf.gz", member("B", {"s1": "0|0:0:0:50:.:22:22,0:99", "s2": "./.:6:0:0:.:2:2,0:0",
                                                          "s3": "./.:0:0:0:.:0:0,0:0", "s4": "1:0:60:35:.:11:0,11:97",
                                                          "s5": "0|0:0:0:50:.:22:22,0:99"}))
        (members / "000002__A.gt.vcf.gz.tbi").write_text("")
        (members / "000003__B.gt.vcf.gz.tbi").write_text("")
        for name, values in (("ids.txt", ["A", "B", "C"]), ("sexes.txt", ["XX", "XY", ""]), ("read_ids.txt", ["A", "B"])):
            (self.tmp / name).write_text("\n".join(values) + "\n")
        return ["sites.vcf.gz", "graph_calls.tsv.gz", "shard.bed", "members", "ids.txt", "sexes.txt", "read_ids.txt", "shard"]

    def test_join_of_one_shard(self):
        res = self.run_script("graph_gt_join.py", self.join_inputs())
        lines = res.stdout.splitlines()
        header = [l for l in lines if l.startswith("##")]
        for tag in ("##INFO=<ID=SVTYPE", "##INFO=<ID=AC,", "##INFO=<ID=AF_XY,", "##FORMAT=<ID=GT,", "##FORMAT=<ID=AD,", "##contig=<ID=chr1"):
            self.assertTrue(any(l.startswith(tag) for l in header), tag)
        self.assertEqual([l for l in lines if l.startswith("#CHROM")][0].split("\t")[8:], ["FORMAT", "A", "B", "C"])
        rec = {l.split("\t")[2]: l.split("\t") for l in lines if not l.startswith("#")}
        self.assertEqual(sorted(rec), ["s1", "s2", "s3", "s4"])            # s5 starts outside the shard
        # a member's column: unphased genotype with GQ, DP, AD; uncalled without them; a member without reads uncalled
        self.assertEqual(rec["s1"][8:], ["GT:GQ:DP:AD", "0/1:30:20:10,10", "0/0:50:22:22,0", "./."])
        self.assertEqual(rec["s2"][9:], ["1/1:40:18:0,18", "./.", "./."])
        self.assertEqual(rec["s4"][9:], ["0/1:30:20:10,10", "1:35:11:0,11", "./."])   # a haploid call stays haploid
        info = {k: dict(kv.split("=") if "=" in kv else (kv, True) for kv in v[7].split(";")) for k, v in rec.items()}
        self.assertEqual([info["s1"][k] for k in ("AC", "AN", "AF", "NS", "AC_XX", "AN_XX", "AC_XY", "AN_XY", "AF_XY", "GRAPH_AC")],
                         ["1", "4", "0.25", "2", "1", "2", "0", "2", "0", "1"])
        self.assertEqual([info["s2"][k] for k in ("AC", "AN", "AF", "NS", "AN_XY", "AF_XY")], ["2", "2", "1", "1", "0", "."])
        self.assertEqual([info["s4"][k] for k in ("AC", "AN", "AC_XY", "AN_XY")], ["2", "3", "1", "1"])
        self.assertNotIn("AC", info["s3"])                                  # not genotyped: no cohort counts
        self.assertEqual(info["s3"]["SKIP"], "too_long")
        members = {(r["sample_id"], r["svtype"]): r for r in read_tsv(self.tmp / "shard.members.tsv")}
        self.assertEqual([members[("A", "INS")][k] for k in ("called", "missing", "hom_ref", "het", "hom_alt", "dp_sum", "dp_n")], ["2", "0", "0", "2", "0", "40", "2"])
        self.assertEqual([members[("A", "DEL")][k] for k in ("called", "het", "hom_alt")], ["1", "0", "1"])
        self.assertEqual([members[("B", "INS")][k] for k in ("called", "hom_ref", "hom_alt")], ["2", "1", "1"])   # the hemizygous call counts with hom_alt
        self.assertEqual([members[("B", "DEL")][k] for k in ("called", "missing")], ["0", "1"])
        self.assertEqual([members[("C", "INS")][k] for k in ("called", "missing")], ["0", "2"])                   # the site too long is not counted
        conc = {(r["sample_id"], r["svtype"], r["size_class"], r["context"]): r for r in read_tsv(self.tmp / "shard.concordance.tsv")}
        # Z is not in this cohort. A: s1 (in a repeat) het in both; s4 (elsewhere) the assembly has one haplotype there (1:1), the
        # reads two: not comparable; s2 (an excluded region wins over a repeat) het against hom
        self.assertEqual(sorted(conc), [("A", "DEL", "1k-9.9k", "excluded_region"), ("A", "INS", "50-999", "other"), ("A", "INS", "50-999", "tandem_repeat")])
        self.assertEqual([conc[("A", "INS", "50-999", "tandem_repeat")][k] for k in ("compared", "concordant", "not_comparable")], ["1", "1", "0"])
        self.assertEqual([conc[("A", "INS", "50-999", "other")][k] for k in ("compared", "concordant", "not_comparable")], ["0", "0", "1"])
        self.assertEqual([conc[("A", "DEL", "1k-9.9k", "excluded_region")][k] for k in ("compared", "concordant", "zygosity")], ["1", "0", "1"])
        counts = {(r["svtype"], r["size_class"]): r for r in read_tsv(self.tmp / "shard.site_counts.tsv")}
        self.assertEqual([counts[("INS", "50-999")][k] for k in ("sites", "genotyped", "variable", "graph_only", "no_call", "in_tandem_repeat",
                                                                 "in_excluded_region")], ["2", "2", "2", "0", "0", "1", "0"])
        self.assertEqual([counts[("DEL", "1k-9.9k")][k] for k in ("in_tandem_repeat", "in_excluded_region")], ["0", "1"])
        self.assertEqual([counts[("INS", "10k+")][k] for k in ("sites", "genotyped")], ["1", "0"])
        self.assertIn("4 sites, 2 members with reads of 3, 1 compared with their assembly", res.stderr)

    def test_join_of_an_empty_shard_and_gather(self):
        args = self.join_inputs()
        (self.tmp / "empty.bed").write_text("chr2\t0\t1000\n")
        first = self.run_script("graph_gt_join.py", args)
        self.join_inputs()
        empty = self.run_script("graph_gt_join.py", args[:2] + ["empty.bed"] + args[3:7] + ["empty"])
        self.assertEqual([l for l in empty.stdout.splitlines() if not l.startswith("#")], [])
        self.assertEqual([l for l in empty.stdout.splitlines() if l.startswith("#")], [l for l in first.stdout.splitlines() if l.startswith("#")])
        self.assertEqual(read_tsv(self.tmp / "empty.members.tsv"), [])
        for kind in ("members", "concordance", "site_counts"):
            (self.tmp / f"{kind}.list").write_text(f"{self.tmp}/shard.{kind}.tsv\n{self.tmp}/empty.{kind}.tsv\n")
        (self.tmp / "sites_counts.tsv").write_text("metric\tkey\tvalue\n" "graph_sample\tA\t1\n" "graph_sample\tZ\t1\n" "records\tread\t9\n"
                                                    "sites\tINS:50-999\t3\n" "sites\tDEL:1k-9.9k\t1\n" "sites\tINS:10k+\t1\n" "not_genotyped\ttoo_long\t1\n")
        self.run_script("graph_gt_gather.py", ["members.list", "concordance.list", "site_counts.list", "sites_counts.tsv", "ids.txt", "sexes.txt",
                                               "read_ids.txt", "G1", "size_min=50;size_max=10000;genotyper=kanpig", "out"])
        rows = {r["sample_id"]: r for r in read_tsv(self.tmp / "out.members.tsv")}
        self.assertEqual(list(rows), ["A", "B", "C"])
        self.assertEqual([rows["A"][k] for k in ("sex", "reads", "in_graph", "sites", "called", "missing", "call_rate", "het", "hom_alt", "ins_carried",
                                                 "del_carried", "mean_dp")], ["XX", "yes", "yes", "3", "3", "0", "1.0000", "2", "1", "2", "1", "19.3"])
        self.assertEqual([rows["B"][k] for k in ("reads", "in_graph", "sites", "called", "call_rate")], ["yes", "no", "3", "2", "0.6667"])
        self.assertEqual([rows["C"][k] for k in ("reads", "sites", "called", "call_rate", "mean_dp")], ["no", "3", "0", "", ""])
        conc = {(r["sample_id"], r["svtype"], r["size_class"], r["context"]): r for r in read_tsv(self.tmp / "out.concordance.tsv")}
        self.assertEqual([conc[("ALL", "ALL", "ALL", "ALL")][k] for k in ("compared", "concordant", "zygosity", "not_comparable", "concordance")],
                         ["2", "1", "1", "1", "0.5000"])
        self.assertEqual((conc[("A", "ALL", "ALL", "ALL")]["compared"], conc[("ALL", "ALL", "ALL", "tandem_repeat")]["concordance"],
                          conc[("ALL", "ALL", "ALL", "excluded_region")]["concordance"], conc[("ALL", "ALL", "ALL", "other")]["concordance"],
                          conc[("A", "ALL", "ALL", "tandem_repeat")]["compared"], conc[("ALL", "INS", "ALL", "ALL")]["compared"]),
                         ("2", "1.0000", "0.0000", "", "1", "1"))
        self.assertEqual(list(conc)[0], ("ALL", "ALL", "ALL", "ALL"))
        summ: dict[str, dict[str, str]] = {}
        for r in read_tsv(self.tmp / "out.summary.tsv"):
            self.assertEqual(r["resource"], "graph_gt")
            summ.setdefault(r["metric"], {})[r["key"]] = r["value"]
        self.assertEqual((summ["graph"]["cohort"], summ["graph"]["members"]), ("G1", "2"))
        self.assertEqual((summ["members"]["cohort"], summ["members"]["with_reads"], summ["members"]["without_reads"], summ["members"]["in_graph_with_reads"]),
                         ("3", "2", "1", "1"))
        # five alleles were listed, four start in a region shard: the fifth is reported, not lost silently
        self.assertEqual((summ["sites"]["listed"], summ["sites"]["in_region_shards"], summ["sites"]["outside_region_shards"]), ("5", "4", "1"))
        self.assertEqual((summ["sites_genotyped"]["ALL:ALL"], summ["sites_variable"]["INS:50-999"], summ["sites_sites"]["INS:ALL"]), ("3", "2", "3"))
        self.assertEqual((summ["genotypes"]["called"], summ["genotypes"]["possible"], summ["genotypes"]["call_rate"]), ("5", "6", "0.8333"))
        self.assertEqual((summ["concordance"]["share"], summ["parameter"]["genotyper"], summ["site_list_not_genotyped"]["too_long"]), ("0.5000", "kanpig", "1"))
        self.assertEqual((summ["concordance_tandem_repeat"]["share"], summ["concordance_excluded_region"]["compared"], summ["concordance_other"]["compared"],
                          summ["sites_in_tandem_repeat"]["ALL:ALL"], summ["sites_in_excluded_region"]["DEL:1k-9.9k"]), ("1.0000", "1", "0", "1", "1"))


class PlannerAndBuilderTest(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.cfg = make_project(self.tmp)
        self.db = DB(self.cfg.db_path)
        self.events = Events(self.db, self.cfg.logs_dir / "events.jsonl")
        samples.add_samples(self.db, self.events, write_tsv(self.tmp, [
            {"sample_id": sid, "sex": "MALE" if sid == "S1" else "FEMALE" if sid == "S2" else "", "hifi_reads": f"{sid}.bam"} for sid in MEMBERS]))
        self.cohort = cohorts.freeze(self.db, self.events, "C1", write_ids(self.tmp, MEMBERS))
        cohorts.freeze(self.db, self.events, "G1", write_ids(self.tmp, ["S1", "S2"], name="g.txt"))

    def tearDown(self):
        self.db.close()

    def blocked(self, **sel) -> dict[str, str]:
        p = plan.compute(self.db, self.cfg, plan.Selection(mode="standalone", cohort="C1", **sel), VERSION)
        self.runnable = [c.stage for c in p.runnable]
        return {b.candidate.stage: b.reason for b in p.blocked}

    def test_the_graph_is_a_prerequisite_in_another_cohort(self):
        # a project that names no graph cohort has no such stage; asked for by name, it says what is missing
        self.assertNotIn("cohort_graph_gt", self.blocked())
        self.assertEqual(self.blocked(stage="cohort_graph_gt")["cohort_graph_gt"], plan.NO_GRAPH_COHORT)
        self.cfg.graph_gt_cohort = "NOPE"
        self.assertIn("NOPE is not a frozen cohort", self.blocked()["cohort_graph_gt"])
        # the members' own per-sample stage comes first ...
        self.cfg.graph_gt_cohort = "G1"
        self.assertIn("waiting for singleton of S1, singleton of S2, singleton of S3", self.blocked()["cohort_graph_gt"])
        for sid in MEMBERS:
            seed_success(self.cfg, self.db, "singleton", "sample", sid, {"sample_id": sid})
        # ... then the graph, of the other cohort, in the other mode
        self.assertEqual(self.blocked()["cohort_graph_gt"], "waiting for cohort_pangenome of cohort G1 (config.json graph_gt_cohort)")
        seed_success(self.cfg, self.db, "cohort_pangenome", "cohort", "G1", {"cohort_id": "G1", "sample_ids": ["S1", "S2"]}, mode="assembly",
                     version="0.1.0")
        # runnable although no cohort stage of this cohort has run (it reads the members' alignments, nothing of
        # cohort_merge or cohort_methyl), with a graph built at another version
        blocked = self.blocked()
        self.assertNotIn("cohort_graph_gt", blocked)
        self.assertIn("cohort_graph_gt", self.runnable)
        self.assertIn("cohort_freq", blocked)

    def test_builder(self):
        for sid in MEMBERS:
            seed_success(self.cfg, self.db, "singleton", "sample", sid, {"sample_id": sid})
        graph = seed_success(self.cfg, self.db, "cohort_pangenome", "cohort", "G1", {"cohort_id": "G1", "sample_ids": ["S1", "S2"]}, mode="assembly",
                             version="0.1.0")
        ctx = inputs.BuildContext(self.cfg, self.db, "standalone", VERSION, "cohort", "C1", cohort=self.cohort)
        with self.assertRaisesRegex(Exception, "names no graph_gt_cohort"):
            inputs.generate(ctx, "cohort_graph_gt")
        self.cfg.graph_gt_cohort = "G1"
        doc, _ = inputs.generate(ctx, "cohort_graph_gt")
        d = {k.split(".", 1)[1]: v for k, v in doc.items()}
        self.assertEqual((d["cohort_id"], d["graph_cohort_id"], d["read_sample_ids"]), ("C1", "G1", d["sample_ids"]))
        self.assertEqual(d["sample_ids"], self.cohort.members)
        self.assertEqual(d["sample_sexes"], [{"S1": "XY", "S2": "XX", "S3": ""}[s] for s in d["sample_ids"]])
        self.assertEqual(d["read_sample_sexes"], d["sample_sexes"])
        self.assertTrue(d["graph_vcf"].endswith("/pangenome_wave_vcf/G1.pangenome.GRCh38_GIABv3.wave.vcf.gz"))
        self.assertIn("/cohorts/G1/0.1.0/cohort_pangenome/", d["graph_vcf"])
        self.assertTrue(d["graph_vcf_index"].endswith(".wave.vcf.gz.tbi"))
        self.assertEqual([Path(b).name for b in d["aligned_bams"]], [f"{s}.GRCh38_GIABv3.haplotagged.bam" for s in d["sample_ids"]])
        self.assertEqual([Path(b).name for b in d["aligned_bam_indices"]], [f"{s}.GRCh38_GIABv3.haplotagged.bam.bai" for s in d["sample_ids"]])
        self.assertTrue(d["ugc_wgw_container_registry"])
        # the graph's run is in the manifest's record of what the stage read
        self.assertIn({"cohort_id": "G1", "stage": "cohort_pangenome", "ugc_wgw_version": "0.1.0", "run_id": graph.run_id}, ctx.member_runs)
        # a member whose run has no haplotagged alignment is left out of the genotyping and stays a column
        run = self.db.latest_success("sample", "S3", "singleton", VERSION)
        out = Path(run.run_path) / "outputs.json"
        doc3 = json.loads(out.read_text())
        for key in ("ugc_wgw_singleton.merged_haplotagged_bam", "ugc_wgw_singleton.merged_haplotagged_bam_index"):
            doc3[key] = None
        out.write_text(json.dumps(doc3))
        doc, _ = inputs.generate(inputs.BuildContext(self.cfg, self.db, "standalone", VERSION, "cohort", "C1", cohort=self.cohort), "cohort_graph_gt")
        d = {k.split(".", 1)[1]: v for k, v in doc.items()}
        self.assertEqual((sorted(d["read_sample_ids"]), len(d["aligned_bams"]), len(d["sample_ids"]), len(d["read_sample_sexes"])), (["S1", "S2"], 2, 3, 2))


class EndToEndTest(unittest.TestCase):
    """The stage through the driver with the fake engine: submitted like any other once the project names a graph
    cohort, its section on the cohort page, its tables in the export."""

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.cfg = make_project(self.tmp)
        self.proj = str(self.cfg.project_dir)
        tsv = write_tsv(self.tmp, [{"sample_id": sid, "sex": "MALE", "hifi_reads": f"{sid}.bam"} for sid in MEMBERS])
        self.assertEqual(run_cli(["--project", self.proj, "samples", "add", str(tsv)])[0], 0)
        self.assertEqual(run_cli(["--project", self.proj, "cohort", "freeze", "C1", "--samples", str(write_ids(self.tmp, MEMBERS))])[0], 0)
        self.assertEqual(run_cli(["--project", self.proj, "cohort", "freeze", "G1", "--samples", str(write_ids(self.tmp, ["S1", "S2"], name="g.txt"))])[0], 0)

    def test_submit_page_and_export(self):
        code, out, err = run_cli(["--project", self.proj, "submit", "--mode", "standalone", "--cohort", "C1"])
        self.assertEqual(code, 0, err)
        self.assertNotIn("cohort_graph_gt", out + err)                      # no graph cohort named: the stage does not exist for the project
        code, out, err = run_cli(["--project", self.proj, "submit", "--mode", "standalone", "--cohort", "C1", "--stage", "cohort_graph_gt"])
        self.assertIn("no graph_gt_cohort in the project's config.json", out + err)
        cfg = config.load(self.cfg.project_dir)
        cfg.graph_gt_cohort = "G1"
        config.save(cfg)
        code, out, err = run_cli(["--project", self.proj, "submit", "--mode", "standalone", "--cohort", "C1"])
        self.assertIn("waiting for cohort_pangenome of cohort G1", out + err)
        code, out, err = run_cli(["--project", self.proj, "submit", "--mode", "assembly", "--cohort", "G1"])
        self.assertEqual(code, 0, err)
        code, out, err = run_cli(["--project", self.proj, "submit", "--mode", "standalone", "--cohort", "C1"])
        self.assertEqual(code, 0, err)
        db = DB(cfg.db_path)
        run = db.latest_success("cohort", "C1", "cohort_graph_gt", VERSION)
        db.close()
        self.assertIsNotNone(run)
        outs = json.loads((Path(run.run_path) / "outputs.json").read_text())
        self.assertTrue(outs["ugc_wgw_cohort_graph_gt.graph_gt_vcf"].endswith("C1.graph_gt.G1.GRCh38_GIABv3.vcf.gz"))
        self.assertEqual(len(outs["ugc_wgw_cohort_graph_gt.plups"]), 3)
        man = json.loads((Path(run.run_path) / "run_manifest.json").read_text())
        self.assertIn("cohort_pangenome", json.dumps(man))                    # which graph run the stage read
        # the cohort page
        code, out, err = run_cli(["--project", self.proj, "summary", "--cohort", "C1"])
        self.assertEqual(code, 0, err)
        out_dir = cfg.results_dir / "reports" / "summary"
        text = (out_dir / "C1.summary.html").read_text()
        checker = HtmlChecker()
        checker.feed(text)
        self.assertEqual((checker.errors, checker.stack), ([], []))
        for needle in ('id="sec-graph-gt"', "Structural variants of the graph in every member (cohort_graph_gt)", "graph of cohort <strong>G1</strong>",
                       "SV alleles in the list", "not genotyped: longer than 10 kb", "3 of 3", "Reads against assemblies",
                       "reads agree with the assembly, outside repeats and excluded regions", "97.0% of 600 genotypes",
                       "in a tandem repeat (TRGT catalogue)", "93.4%", "exclusion list", "80.0%", "by type and length",
                       "carrier missed by the reads", "Sites by type and length", "1k-9.9k", "graph_gt_missing_rate", "per member of the graph",
                       "Median call rate 97.0%"):
            self.assertIn(needle, text, needle)
        cd = json.loads((out_dir / "C1.summary.json").read_text())
        gg = cd["callsets"]["cohort_graph_gt"]
        self.assertEqual((gg["summary"]["graph"]["cohort"], gg["summary"]["sites"]["listed"], gg["summary"]["concordance"]["share"]), ("G1", 1000, 0.9368))
        self.assertEqual((gg["members"]["flagged"], round(gg["members"]["median_call_rate"], 4), len(gg["members"]["rows"])), (["S3"], 0.9697, 3))
        self.assertEqual({r["sample_id"] for r in gg["concordance"]}, {"ALL", "S2", "S1"})
        self.assertEqual(gg["summary"]["concordance_other"]["share"], 0.97)
        self.assertGreater(gg["graph_gt_vcf"]["bytes"], 0)
        # a cohort without the stage has no such section (the page of the graph's own cohort, in assembly mode)
        code, out, err = run_cli(["--project", self.proj, "summary", "--cohort", "G1", "--mode", "assembly"])
        self.assertEqual(code, 0, err)
        self.assertNotIn("sec-graph-gt", (out_dir / "G1.summary.html").read_text() if (out_dir / "G1.summary.html").exists() else "")
        # the export copies the stage's tables
        dest = self.tmp / "export"
        code, out, err = run_cli(["--project", self.proj, "export", "--cohort", "C1", "--out", str(dest), "--domains", "graph_gt"])
        self.assertEqual(code, 0, err)
        self.assertEqual(sorted(p.name for p in dest.glob("C1.graph_gt.*")),
                         ["C1.graph_gt.concordance.tsv", "C1.graph_gt.members.tsv", "C1.graph_gt.skipped.tsv", "C1.graph_gt.summary.tsv"])


if __name__ == "__main__":
    unittest.main()
