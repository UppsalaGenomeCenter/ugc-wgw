"""`ugc-wgw export`: the frequency tables, the per-member rows and the command end to end."""
import gzip
import json
import tempfile
import unittest
from pathlib import Path

from ugc_wgw import cli, export, summary
from ugc_wgw.db import DB

from .helpers import make_digest, make_project, run_cli, seed_files, write_ids, write_tsv


def _read_tsv(path: Path) -> tuple[str, list[str], list[dict[str, str]]]:
    lines = path.read_text().splitlines()
    comment = lines[0]
    cols = lines[1].split("\t")
    rows = [dict(zip(cols, line.split("\t"))) for line in lines[2:]]
    return comment, cols, rows


class FreqTableTest(unittest.TestCase):
    def test_counts_fractions_strata(self):
        rows = [{"sample_id": "A", "gene": "G", "dip": "*1/*1", "called": True},
                {"sample_id": "B", "gene": "G", "dip": "*1/*1", "called": True},
                {"sample_id": "C", "gene": "G", "dip": "*1/*2", "called": True},
                {"sample_id": "D", "gene": "G", "dip": "NO_READS", "called": False}]
        strata = export.strata_of({sid: {"sample": {"sex_sheet": s}, "coverage": {}} for sid, s in
                                   (("A", "MALE"), ("B", "FEMALE"), ("C", "MALE"), ("D", ""))}, ("sex",))
        self.assertEqual(strata["sex"], {"A": "XY", "B": "XX", "C": "XY"})
        out = export.freq_table(rows, "dip", strata, by=("gene",))
        every = [r for r in out if r["stratum_kind"] == "ALL"]
        self.assertEqual([(r["dip"], r["count"], r["fraction"], r["n_called"], r["n_members"]) for r in every],
                         [("*1/*1", 2, 0.6667, 3, 4), ("*1/*2", 1, 0.3333, 3, 4)])
        xy = [r for r in out if r["stratum"] == "XY"]
        self.assertEqual([(r["dip"], r["count"], r["n_called"], r["n_members"]) for r in xy], [("*1/*1", 1, 2, 2), ("*1/*2", 1, 2, 2)])
        self.assertEqual(sum(float(r["fraction"]) for r in xy), 1.0)
        self.assertEqual([r["stratum"] for r in out if r["stratum_kind"] == "sex" and r["stratum"] == "XX"], ["XX"])
        self.assertNotIn("D", strata["sex"])

    def test_empty_and_single(self):
        self.assertEqual(export.freq_table([], "x", {"ALL": {}}, by=("g",)), [])
        one = export.freq_table([{"sample_id": "A", "g": "G", "x": 3, "called": True}], "x", {"ALL": {"A": "ALL"}}, by=("g",))
        self.assertEqual((one[0]["x"], one[0]["fraction"], one[0]["n_called"]), (3, 1.0, 1))

    def test_parse_domains(self):
        self.assertEqual(export.parse_domains(None), list(export.DOMAINS))
        self.assertEqual(export.parse_domains(["cn,pgx", "cn", "members"]), ["cn", "pgx", "members"])
        with self.assertRaises(summary.UgcError):
            export.parse_domains(["nope"])


class RowsTest(unittest.TestCase):
    def test_cn_rows_member_without_paraphase(self):
        d1 = make_digest("S1")
        d2 = make_digest("S2", targeted={"paraphase": {}})
        rows = export.cn_rows({"S1": d1, "S2": d2})
        s1 = [r for r in rows if r["sample_id"] == "S1"]
        s2 = [r for r in rows if r["sample_id"] == "S2"]
        self.assertEqual(len(s1), 3 + 2)
        self.assertEqual([r["region"] for r in s2], ["KIV2", "D4Z4"])
        smn1 = next(r for r in s1 if r["region"] == "smn1")
        self.assertEqual((smn1["failed"], smn1["called"], smn1["total_cn"]), (True, False, None))
        cyp = next(r for r in s1 if r["region"] == "CYP2D6")
        self.assertEqual((cyp["called"], cyp["total_cn"], cyp["gene_cn"], cyp["haplotypes_by_gene"], cyp["fusions"]),
                         (True, 4, 2, "CYP2D6:2;CYP2D7:2", "CYP2D6-CYP2D7"))
        kiv2 = next(r for r in s1 if r["region"] == "KIV2")
        self.assertEqual((kiv2["called"], kiv2["failed"], kiv2["total_cn"]), (False, True, None))
        d4z4 = next(r for r in s1 if r["region"] == "D4Z4")
        self.assertEqual((d4z4["called"], d4z4["total_cn"], d4z4["genes"]), (True, 38.0, "DUX4"))
        tg = export.targeted_summary({"S1": d1, "S2": d2}, export.strata_of({"S1": d1, "S2": d2}, ("sex",)))
        freq = tg["paraphase"]["freq"]
        self.assertTrue(all(r["region"] != "KIV2" for r in freq))
        cyp_all = [r for r in freq if r["region"] == "CYP2D6" and r["stratum_kind"] == "ALL"]
        self.assertEqual([(r["metric"], r["value"], r["count"], r["n_called"], r["n_members"]) for r in cyp_all],
                         [("total_cn", 4, 1, 1, 1), ("gene_cn", 2, 1, 1, 1)])
        self.assertEqual([r["region"] for r in tg["paraphase"]["regions"]][:1], ["CYP2D6"])   # ties by n_called, then name
        self.assertEqual(tg["paraphase"]["n_members_with_failed"], 1)

    def test_pgx_and_hla_rows_uncalled_gene(self):
        d = make_digest("S1")
        pgx = export.pgx_rows({"S1": d})
        genes = {r["gene"]: r for r in pgx}
        self.assertNotIn("HLA-A", genes)
        self.assertFalse(genes["TPMT"]["called"])
        self.assertEqual((genes["CYP2D6"]["haplotype_1"], genes["CYP2D6"]["haplotype_2"], genes["CYP2D6"]["called"]), ("*1", "*4", True))
        freq = export.freq_table(pgx, "simple_diplotype", {"ALL": {"S1": "ALL"}}, by=("gene",))
        self.assertTrue(all(r["gene"] != "TPMT" for r in freq))
        hla = export.hla_rows({"S1": d})
        by = {r["gene"]: r for r in hla}
        self.assertEqual(sorted(by), ["HLA-A", "HLA-DRB3"])
        self.assertFalse(by["HLA-DRB3"]["called"])
        alleles = export.freq_table(export.hla_allele_rows(hla), "allele", {"ALL": {"S1": "ALL"}}, by=("gene",))
        self.assertEqual([(r["allele"], r["count"], r["n_called"]) for r in alleles if r["gene"] == "HLA-A"],
                         [("*02:01:01:01", 1, 2), ("*03:01:01:01", 1, 2)])
        # the TSV fallback when the JSON is absent
        d2 = make_digest("S2", targeted={"starphase_meta": {}})
        pgx2 = {r["gene"]: r for r in export.pgx_rows({"S2": d2})}
        self.assertEqual((pgx2["CYP2C19"]["haplotype_1"], pgx2["CYP2C19"]["called"], pgx2["TPMT"]["called"]), ("*1", True, False))
        self.assertEqual([r["gene"] for r in export.hla_rows({"S2": d2})], ["HLA-A", "HLA-DRB3"])

    def test_mito_rows(self):
        d = make_digest("S1", targeted={"mitorsaw": {"fingerprint": {"passing_explanation": False, "unexplained_fraction": 0.05}}})
        rows = export.mito_rows({"S1": d})
        self.assertEqual(len(rows), 2)
        self.assertEqual({r["unexplained_fraction"] for r in rows}, {0.05})
        self.assertTrue(all(r["mito_unexplained_high"] for r in rows))
        self.assertEqual((rows[0]["haplogroup"], rows[1]["n_haplotypes"]), ("", 2))
        tg = export.targeted_summary({"S1": d}, {"ALL": {"S1": "ALL"}})
        self.assertEqual((tg["mito"]["multi_haplotype"], tg["mito"]["minor_abundance_median"], tg["mito"]["flagged"]), (["S1"], 0.03, ["S1"]))
        # a run before 0.15 (no collapsed summary, no haplogroups of its own): the haplogroups come from cohort_merge's gathered
        # table and the most abundant haplotype is the major one, so the cohort counts it (the 864 cohort counted none)
        legacy = make_digest("S3")
        legacy["targeted"]["mitorsaw"]["collapsed"] = {}
        for h in legacy["targeted"]["mitorsaw"]["haplotypes"]:
            h["haplogroup"] = ""
        gathered = [{"sample_id": "S3", "haplotype": "hap_0", "haplogroup": "H1h1", "quality": 0.87, "made_in_stage": "yes"},
                    {"sample_id": "S3", "haplotype": "hap_1", "haplogroup": "H1", "quality": 0.8, "made_in_stage": "yes"}]
        rows = export.mito_rows({"S3": legacy}, gathered)
        self.assertEqual([(r["haplotype"], r["haplogroup"], r["haplogroup_source"], r["major_haplotype"]) for r in rows],
                         [("hap_0", "H1h1", "cohort_merge", True), ("hap_1", "H1", "cohort_merge", False)])
        self.assertEqual(export.haplogroup_rows(rows), [{"sample_id": "S3", "haplogroup": "H1h1", "macro_haplogroup": "H", "called": True}])
        tg = export.targeted_summary({"S3": legacy}, {"ALL": {"S3": "ALL"}}, cohort_haplogroups=gathered)
        self.assertEqual((tg["mito"]["n_haplogrouped"], tg["mito"]["haplogroups"]), (1, {"H1h1": 1}))
        self.assertEqual(export.mito_rows({"S3": legacy})[0]["major_haplotype"], True)   # abundance decides even without a gathered table

    def test_member_and_variant_rows(self):
        d = make_digest("S1")
        v = export.variants_rows({"S1": d})[0]
        self.assertEqual((v["snv"], v["n_hets"], v["roh_total_bp"], v["sv_total"], v["trgt_disease_called"], v["foo_extra"]),
                         (3800000.0, 2000000, 10000000, 25.0, 69, "7"))
        self.assertNotIn("SNV_count", v)
        ph = export.phasing_rows({"S1": d})[0]
        self.assertEqual((ph["phased_fraction"], ph["blocks"], ph["block_n50"]), (0.9, 3, 400000))
        cov = export.coverage_rows({"S1": d})[0]
        self.assertEqual((cov["depth_chr1"], cov["chrX_ratio"], cov["read_count"]), (30.0, 0.5, 1.0))


class ExportEndToEnd(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.cfg = make_project(self.tmp)
        self.proj = str(self.cfg.project_dir)
        tsv = write_tsv(self.tmp, [{"sample_id": "S1", "sex": "MALE", "hifi_reads": "s1.bam", "batch": "b1"},
                                   {"sample_id": "S2", "sex": "FEMALE", "hifi_reads": "s2.bam", "batch": "b2"}])
        self.assertEqual(run_cli(["--project", self.proj, "samples", "add", str(tsv)])[0], 0)
        self.assertEqual(run_cli(["--project", self.proj, "cohort", "freeze", "C1", "--samples", str(write_ids(self.tmp, ["S1", "S2"]))])[0], 0)
        code, out, err = run_cli(["--project", self.proj, "submit", "--mode", "standalone", "--cohort", "C1"],
                                 env={"UGC_WGW_FAKE_INFERRED_SEX": "S1=FEMALE,S2=FEMALE"})
        self.assertEqual(code, 0, err)
        self.export_dir = self.cfg.results_dir / "cohorts" / "C1" / cli.config.read_version(self.cfg.code_dir) / "export"

    def test_all_domains(self):
        code, out, err = run_cli(["--project", self.proj, "export", "--cohort", "C1"])
        self.assertEqual(code, 0, err)
        paths = [Path(p) for p in out.split()]
        self.assertTrue(all(p.parent == self.export_dir for p in paths), out)
        names = {p.name for p in paths}
        for n in ("C1.members.tsv", "C1.coverage.tsv", "C1.variants.tsv", "C1.phasing.tsv", "C1.cn.tsv", "C1.cn.freq.tsv",
                  "C1.pgx.tsv", "C1.pgx.freq.tsv", "C1.hla.tsv", "C1.hla.freq.tsv", "C1.mito.tsv", "C1.mito.freq.tsv"):
            self.assertIn(n, names)
        tr_names = ["C1.tr.loci.tsv.gz", "C1.tr.outliers.tsv.gz", "C1.tr.outliers_top.tsv", "C1.tr.samples.tsv", "C1.tr.screen.tsv", "C1.tr.screen_loci.tsv",
                    "C1.tr.spectrum.tsv.gz", "C1.tr.summary.tsv"]   # cohort_tr's tables, the three big ones still bgzipped
        self.assertEqual(sorted(n for n in names if n.startswith("C1.tr.")), tr_names)
        with gzip.open(self.export_dir / "C1.tr.spectrum.tsv.gz", "rt") as fh:
            self.assertEqual(fh.readline().rstrip("\n").split("\t"), ["trid", "al", "n_alleles", "n_members"])
        self.assertIn("assembly: no member has an assembly", err)
        comment, cols, rows = _read_tsv(self.export_dir / "C1.members.tsv")
        self.assertTrue(comment.startswith("# ugc-wgw "))
        self.assertEqual([r["sample_id"] for r in rows], ["S1", "S2"])
        self.assertEqual([(r["sex_sheet"], r["sex_inferred"], r["sex_group"], r["sheet_batch"], r["present"]) for r in rows],
                         [("MALE", "FEMALE", "XY", "b1", "true"), ("FEMALE", "FEMALE", "XX", "b2", "true")])
        self.assertIn("sex_mismatch", rows[0]["qc_flags"])
        self.assertTrue(rows[0]["run_singleton"])
        _c, _cols, cn = _read_tsv(self.export_dir / "C1.cn.tsv")
        self.assertEqual(len(cn), 2 * (3 + 2))
        _c, _cols, cnf = _read_tsv(self.export_dir / "C1.cn.freq.tsv")
        hit = [r for r in cnf if r["region"] == "CYP2D6" and r["metric"] == "total_cn" and r["stratum_kind"] == "ALL"]
        self.assertEqual([(r["value"], r["count"], r["fraction"], r["n_members"]) for r in hit], [("4", "2", "1", "2")])
        self.assertTrue(any(r["stratum"] == "XY" for r in cnf))
        _c, _cols, pf = _read_tsv(self.export_dir / "C1.pgx.freq.tsv")
        self.assertTrue(any(r["gene"] == "CYP2C19" and r["simple_diplotype"] == "*1/*17" and r["count"] == "2" for r in pf))
        _c, cols, mito = _read_tsv(self.export_dir / "C1.mito.tsv")
        self.assertIn("haplogroup", cols)
        man = json.loads((self.export_dir / "export_manifest.json").read_text())
        self.assertEqual((man["schema"], man["digest_schema"], man["cohort"]["n"], man["members"]["with_digest"]), (1, summary.DIGEST_SCHEMA, 2, 2))
        self.assertTrue(man["source_runs"]["members"]["S1"]["singleton"])
        self.assertEqual(sorted(man["domains"]["tr"]["files"]), tr_names)
        self.assertTrue(man["domains"]["tr"]["source_run"].startswith("C1-cohort_tr-"))
        self.assertEqual(man["strata"]["kinds"], ["ALL", "sex"])
        self.assertIsNone(man["strata"]["ancestry"])
        for f in man["files"]:
            p = self.export_dir / f["name"]
            self.assertTrue(p.exists())
            text = gzip.decompress(p.read_bytes()).decode() if p.name.endswith(".gz") else p.read_text()
            self.assertEqual(max(0, sum(1 for line in text.splitlines() if not line.startswith("#")) - 1), f["rows"])   # an empty table has 0 rows
        # the digests are shared with `summary`
        self.assertTrue((self.cfg.results_dir / "reports" / "summary" / "S1.summary.json").exists())

    def test_domains_selection_and_errors(self):
        code, out, err = run_cli(["--project", self.proj, "export", "--cohort", "C1", "--domains", "members,cn"])
        self.assertEqual(code, 0, err)
        self.assertEqual(sorted(Path(p).name for p in out.split()), ["C1.cn.freq.tsv", "C1.cn.tsv", "C1.members.tsv"])
        code, out, err = run_cli(["--project", self.proj, "export", "--cohort", "C1", "--domains", "nope"])
        self.assertEqual(code, 1)
        self.assertIn("unknown export domain", err)
        code, out, err = run_cli(["--project", self.proj, "export", "--cohort", "C1", "--domains", "tr"])
        self.assertEqual(code, 0, err)
        self.assertEqual(len(out.split()), 8)
        code, out, err = run_cli(["--project", self.proj, "export", "--cohort", "C1", "--domains", "assembly"])
        self.assertEqual(code, 1)
        self.assertIn("nothing written", err)
        code, out, err = run_cli(["--project", self.proj, "export", "--cohort", "NOPE"])
        self.assertEqual((code, "unknown cohort" in err), (1, True))
        code, out, err = run_cli(["--project", self.proj, "export", "--cohort", "C1", "--mode", "joint"])
        self.assertEqual(code, 1)
        self.assertIn("no member of C1 has a successful joint run", err)

    def test_strata_without_cohort_qc(self):
        code, out, err = run_cli(["--project", self.proj, "export", "--cohort", "C1", "--domains", "pgx", "--strata", "ancestry"])
        self.assertEqual(code, 0, err)
        self.assertIn("ancestry strata: no cohort_qc run", err)
        _c, _cols, pf = _read_tsv(self.export_dir / "C1.pgx.freq.tsv")
        self.assertEqual({r["stratum_kind"] for r in pf}, {"ALL"})
        code, out, err = run_cli(["--project", self.proj, "export", "--cohort", "C1", "--domains", "pgx", "--strata", "none", "--force"])
        self.assertEqual(code, 0, err)
        _c, _cols, pf = _read_tsv(self.export_dir / "C1.pgx.freq.tsv")
        self.assertEqual({r["stratum_kind"] for r in pf}, {"ALL"})

    def test_identity_copy_and_ancestry_strata(self):
        db = DB(self.cfg.db_path)
        try:
            pairs = "#sample_a\tsample_b\trelatedness\tibs0\tibs2\nS1\tS2\t0.01\t40\t200\n"
            anc = ("#sample_id\tpredicted_ancestry\tgiven_ancestry\tEAS_prob\tAFR_prob\tAMR_prob\tSAS_prob\tEUR_prob\tPC1\tPC2\n"
                   "HG00096\tEUR\tEUR\t0\t0\t0\t0\t1\t0.1\t0.2\nS1\tEUR\t\t0\t0\t0\t0\t1\t0.1\t0.2\nS2\tAFR\t\t0\t1\t0\t0\t0\t-3\t0.2\n")
            members = "#sample_id\tsomalier_sex\n S1\tMALE\nS2\tFEMALE\n".replace(" S1", "S1")
            qc_run = seed_files(self.cfg, db, "cohort_qc", "cohort", "C1",
                                {"qc_pairs": pairs, "qc_samples": "#sample_id\nS1\nS2\n", "qc_groups": "", "qc_summary": "section\tmetric\tkey\tvalue\n",
                                 "qc_members": members, "qc_related_pairs": pairs, "qc_ancestry": anc})
        finally:
            db.close()
        code, out, err = run_cli(["--project", self.proj, "export", "--cohort", "C1", "--domains", "identity,pgx,members", "--strata", "all"])
        self.assertEqual(code, 0, err)
        self.assertNotIn("ancestry strata", err)
        self.assertEqual((self.export_dir / "C1.identity.pairs.tsv").read_text(), pairs)
        self.assertTrue((self.export_dir / "C1.identity.ancestry.tsv").exists())
        _c, _cols, pf = _read_tsv(self.export_dir / "C1.pgx.freq.tsv")
        self.assertEqual({r["stratum"] for r in pf if r["stratum_kind"] == "ancestry"}, {"EUR", "AFR"})
        comment, _cols, rows = _read_tsv(self.export_dir / "C1.members.tsv")
        self.assertEqual([(r["ancestry"], r["sex_somalier"]) for r in rows], [("EUR", "MALE"), ("AFR", "FEMALE")])
        self.assertIn("1000 Genomes", comment)
        man = json.loads((self.export_dir / "export_manifest.json").read_text())
        self.assertEqual(man["source_runs"]["cohort"]["cohort_qc"], qc_run.run_id)
        self.assertEqual(man["strata"]["ancestry"]["unlabelled"], [])
        self.assertEqual(man["domains"]["identity"]["files"][0], "C1.identity.pairs.tsv")

    def test_out_and_force(self):
        other = self.tmp / "elsewhere"
        code, out, err = run_cli(["--project", self.proj, "export", "--cohort", "C1", "--domains", "members", "--out", str(other)])
        self.assertEqual(code, 0, err)
        self.assertTrue((other / "C1.members.tsv").exists())
        filled = self.tmp / "filled"
        filled.mkdir()
        (filled / "stray.txt").write_text("x")
        code, out, err = run_cli(["--project", self.proj, "export", "--cohort", "C1", "--domains", "members", "--out", str(filled)])
        self.assertEqual(code, 1)
        self.assertIn("no export_manifest.json", err)
        (filled / "C1.hla.tsv").write_text("stale")
        code, out, err = run_cli(["--project", self.proj, "export", "--cohort", "C1", "--domains", "members", "--out", str(filled), "--force"])
        self.assertEqual(code, 0, err)
        self.assertFalse((filled / "C1.hla.tsv").exists())
        self.assertIn("removed stale C1.hla.tsv", err)
        self.assertTrue((filled / "stray.txt").exists())
        # without --force a stale table survives, the manifest describes only this export
        (filled / "C1.hla.tsv").write_text("stale")
        code, out, err = run_cli(["--project", self.proj, "export", "--cohort", "C1", "--domains", "members", "--out", str(filled)])
        self.assertEqual(code, 0, err)
        self.assertTrue((filled / "C1.hla.tsv").exists())
        man = json.loads((filled / "export_manifest.json").read_text())
        self.assertEqual([f["name"] for f in man["files"]], ["C1.members.tsv"])

    def test_json_and_threshold(self):
        code, out, err = run_cli(["--project", self.proj, "export", "--cohort", "C1", "--domains", "mito", "--json",
                                  "--threshold", "mito_unexplained_max=0.5"])
        self.assertEqual(code, 0, err)
        man = json.loads(out)
        self.assertEqual(man["thresholds"]["mito_unexplained_max"], 0.5)
        self.assertEqual(man["domains"]["mito"]["files"], ["C1.mito.tsv", "C1.mito.freq.tsv"])


class DocsTest(unittest.TestCase):
    def test_export_flags_documented(self):
        parser = cli.build_parser()
        sub = next(a for a in parser._actions if isinstance(a, cli.argparse._SubParsersAction))
        ex = sub.choices["export"]
        guide = (Path(__file__).resolve().parents[2] / "docs" / "guide" / "11-reference.md").read_text()
        for action in ex._actions:
            for opt in action.option_strings:
                if opt.startswith("--"):
                    self.assertIn(opt, guide, f"{opt} is not documented in guide 11")
        self.assertIn("mito_unexplained_max", guide)
        self.assertIn("mito_unexplained_max", (Path(__file__).resolve().parents[2] / "docs" / "guide" / "13-qc-flags.md").read_text())

    def test_every_flag_and_threshold_is_in_the_flags_chapter(self):
        """Chapter 13 is the reference for flags: every flag code the driver can raise and every threshold key is named
        there, so a new flag cannot ship undocumented."""
        import re
        from ugc_wgw import summary
        root = Path(__file__).resolve().parents[2]
        chapter = (root / "docs" / "guide" / "13-qc-flags.md").read_text()
        code = (root / "bin" / "ugc_wgw" / "summary.py").read_text() + (root / "bin" / "ugc_wgw" / "summary_cohort.py").read_text()
        flags = set(re.findall(r'flag\((?:[a-z_]+, )?"([a-z_]+)", "(?:fail|warn|info)"', code)) | set(re.findall(r'"code": "([a-z_]+)", "severity"', code))
        flags |= {"low_depth", "low_mapped", "low_read_quality"}       # raised from a table of (code, key, threshold)
        self.assertGreaterEqual(len(flags), 25, sorted(flags))
        for flag in sorted(flags):
            self.assertIn(f"`{flag}`", chapter, f"flag {flag} is not in chapter 13")
        for key in summary.DEFAULT_THRESHOLDS:
            self.assertIn(f"`{key}`", chapter, f"threshold {key} is not in chapter 13")


if __name__ == "__main__":
    unittest.main()
