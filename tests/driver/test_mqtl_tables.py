"""Tests for the scripts embedded in workflows/ugc_wgw/mqtl/cohort.wdl (cohort_mqtl): the island-by-member matrix,
the cis tests of one shard (the statistics, a planted effect among nulls, a confounder), and the gather."""
from __future__ import annotations

import gzip
import math
import random
import statistics
import subprocess
import sys
import tempfile
import textwrap
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
WDL = REPO / "workflows" / "ugc_wgw" / "mqtl" / "cohort.wdl"


def script(name: str) -> str:
    text = WDL.read_text()
    begin, end = f"# ---- {name} (begin) ----\n", f"# ---- {name} (end) ----\n"
    return textwrap.dedent(text.split(begin, 1)[1].split(end, 1)[0])


def run_script(name: str, args: list[str], cwd: Path) -> subprocess.CompletedProcess:
    res = subprocess.run([sys.executable, "-c", script(name), *args], capture_output=True, text=True, cwd=cwd)
    if res.returncode != 0:
        raise AssertionError(res.stderr)
    return res


def read_tsv(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rt") as fh:  # type: ignore[operator]
        lines = [line.rstrip("\n") for line in fh if not line.startswith("##")]
    cols = lines[0].split("\t")
    return [dict(zip(cols, line.split("\t"))) for line in lines[1:]]


class StatsTest(unittest.TestCase):
    def setUp(self):
        self.ns: dict[str, object] = {}
        exec("import math\nfrom statistics import NormalDist\n" + script("mqtl_stats"), self.ns)

    def test_student_t_against_tabulated_values(self):
        t2 = self.ns["t_two_sided"]
        for t, df, p in ((2.0, 10, 0.073388), (2.228, 10, 0.05), (1.96, 1000, 0.050273), (3.0, 30, 0.005388), (12.706, 1, 0.05), (0.0, 20, 1.0)):
            got = t2(t, df)                                                       # type: ignore[operator]
            self.assertAlmostEqual(got / p, 1.0, delta=0.02, msg=f"t={t} df={df}: {got}")

        def tail(t, df, upper=80.0, steps=400_000):
            """Twice the density's integral from t on (Simpson), for the far tail no table prints."""
            c = math.exp(math.lgamma((df + 1) / 2) - math.lgamma(df / 2)) / math.sqrt(df * math.pi)
            h = (upper - t) / steps
            f = lambda x: c * (1 + x * x / df) ** (-(df + 1) / 2)       # noqa: E731
            s = f(t) + f(upper) + sum((4 if i % 2 else 2) * f(t + i * h) for i in range(1, steps))
            return 2 * s * h / 3

        for t, df in ((5.0, 100), (8.0, 150), (6.0, 40)):
            self.assertAlmostEqual(t2(t, df) / tail(t, df), 1.0, delta=1e-4, msg=f"t={t} df={df}")   # type: ignore[operator]
        self.assertAlmostEqual(t2(-2.0, 10), t2(2.0, 10))                         # type: ignore[operator]
        self.assertTrue(math.isnan(t2(1.0, 0)))                                   # type: ignore[operator, arg-type]

    def test_rank_normal_and_residuals(self):
        rn, basis, residual = self.ns["rank_normal"], self.ns["basis"], self.ns["residual"]
        self.assertAlmostEqual(sum(rn([5.0, 1.0, 3.0, 4.0, 9.0])), 0.0, places=9)   # type: ignore[operator, arg-type]
        z = rn([5.0, 1.0, 3.0, 3.0, 9.0])                                         # type: ignore[operator]
        self.assertEqual(z[2], z[3])                                              # ties share a rank
        self.assertLess(z[1], z[2])
        self.assertLess(z[0], z[4])
        ones, x = [1.0] * 6, [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
        q = basis([ones, x, [2 * v + 3 for v in x]])                              # type: ignore[operator]
        self.assertEqual(len(q), 2)                                               # the third column adds nothing
        r = residual([3 + 2 * v for v in x], q)                                   # type: ignore[operator]
        self.assertTrue(all(abs(v) < 1e-9 for v in r))
        r = residual([0.0, 1.0, 0.0, 1.0, 0.0, 1.0], q)                           # type: ignore[operator]
        self.assertAlmostEqual(sum(r), 0.0, places=9)
        self.assertAlmostEqual(sum(a * b for a, b in zip(r, x)), 0.0, places=9)

    def test_the_pair_test_recovers_a_slope_and_matches_ols(self):
        test, basis = self.ns["test"], self.ns["basis"]
        rnd = random.Random(3)
        n = 200
        g = [float(rnd.choice([0, 0, 0, 1, 1, 2])) for _ in range(n)]
        c = [rnd.gauss(0, 1) for _ in range(n)]
        y = [50 + 8 * gi + 5 * ci + rnd.gauss(0, 4) for gi, ci in zip(g, c)]
        q = basis([[1.0] * n, c])                                                 # type: ignore[operator]
        r, t, df, p, slope = test(y, g, q)                                        # type: ignore[operator, misc]
        self.assertEqual(df, n - 2 - 1)
        self.assertAlmostEqual(slope, 8.0, delta=1.0)
        self.assertGreater(r, 0.7)
        self.assertLess(p, 1e-30)
        self.assertIsNone(test(y, [1.0] * n, q))                                  # type: ignore[operator]  # no variation left in the genotype


class ShardTest(unittest.TestCase):
    N = 160

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        rnd = random.Random(11)
        n = self.N
        self.members = [f"M{i:03d}" for i in range(n)]
        self.pc1 = [rnd.gauss(0, 1) for _ in range(n)]
        self.sex = [rnd.choice(["XX", "XY"]) for _ in range(n)]
        cov = ["sample_id\tused\twhy_not\tsex\tPC1\tPC2\tPC3\tPC4\tPC5\tgranulocyte\tmonocyte\tT-cell\tB-cell\tNK-cell"]
        for i, m in enumerate(self.members):
            cov.append("\t".join([m, "yes", "", self.sex[i], f"{self.pc1[i]:.4f}"] + [f"{rnd.gauss(0, 1):.4f}" for _ in range(4)]
                                 + [f"{rnd.uniform(0.3, 0.6):.4f}", "0.1000", f"{rnd.uniform(0.1, 0.3):.4f}", "0.0500", "0.0300"]))
        cov.append("OUT\tno\texcluded (exclude_ids)\tXX\t\t\t\t\t\t\t\t\t\t")
        (self.tmp / "cov.tsv").write_text("\n".join(cov) + "\n")
        # genotypes
        self.sv_effect = [float(rnd.choice([0, 0, 0, 1, 1, 2])) for _ in range(n)]
        self.sv_confounded = [1.0 if self.pc1[i] + rnd.gauss(0, 0.3) > 0.4 else 0.0 for i in range(n)]
        self.nulls = [[float(rnd.choice([0, 0, 1])) for _ in range(n)] for _ in range(60)]
        self.tr_len = [rnd.choice([30.0, 30.0, 60.0, 90.0, 150.0]) for _ in range(n)]
        self.rare = [1.0 if i < 3 else 0.0 for i in range(n)]
        # islands
        y_effect = [50 + 10 * self.sv_effect[i] + rnd.gauss(0, 5) for i in range(n)]
        y_null = [40 + rnd.gauss(0, 6) for _ in range(n)]
        y_conf = [45 + 9 * self.pc1[i] + rnd.gauss(0, 3) for i in range(n)]
        y_tr = [20 + 0.05 * self.tr_len[i] + rnd.gauss(0, 2) for i in range(n)]
        y_flat = [1.0 + (0.1 if i % 2 else 0.0) for i in range(n)]
        y_gappy = [60 + rnd.gauss(0, 6) for _ in range(n)]
        rows = [("chr1", 100000, 100500, "CpG:_1", "GENEA", "promoter", y_effect), ("chr1", 200000, 200400, "CpG:_2", "GENEB", "gene_body", y_null),
                ("chr1", 300000, 300400, "CpG:_3", "GENEC", "promoter", y_conf), ("chr1", 400000, 400300, "CpG:_4", "", "intergenic", y_tr),
                ("chr1", 500000, 500300, "CpG:_5", "", "intergenic", y_flat), ("chr1", 600000, 600300, "CpG:_6", "", "intergenic", y_gappy),
                ("chr1", 2000000, 2000300, "CpG:_7", "", "intergenic", y_null), ("chrX", 100000, 100300, "CpG:_8", "", "promoter", y_effect)]
        with gzip.open(self.tmp / "matrix.tsv.gz", "wt") as fh:
            fh.write("\t".join(["chrom", "start", "end", "region_label", "gene_symbol", "region_type", "n_members", "mean", "sd"] + self.members + ["OUT"]) + "\n")
            for chrom, s, e, label, gene, rtype, y in rows:
                vals = [f"{v:.1f}" for v in y]
                if label == "CpG:_6":
                    vals = ["" if i < n // 2 else v for i, v in enumerate(vals)]      # half the members without coverage
                fh.write("\t".join([chrom, str(s), str(e), label, gene, rtype, str(n), "0", "0"] + vals + ["99.0"]) + "\n")
        (self.tmp / "shard.bed").write_text("chr1\t0\t1000000\n")
        # the VCF columns are in another order than the members: reversed, with the unused member first
        self.vcf_order = ["OUT"] + self.members[::-1]
        (self.tmp / "samples.txt").write_text("\n".join(self.vcf_order) + "\n")

        def gts(v):
            return ["0/0"] + [{0.0: "./.", 1.0: "0/1", 2.0: "1|1"}[x] for x in v[::-1]]

        sv = [["chr1", "100900", "101200", "svEffect", "DEL", "-300"] + gts(self.sv_effect),
              ["chr1", "320000", "320050", "svConfounded", "INS", "250"] + gts(self.sv_confounded),
              ["chr1", "150000", "150001", "svRare", "INS", "80"] + gts(self.rare),
              ["chr1", "480000", "480100", "svFarFromFlat", "DEL", "-100"] + gts(self.sv_effect),
              ["chrX", "100100", "100200", "svX", "DEL", "-100"] + gts(self.sv_effect)]
        for k, v in enumerate(self.nulls):
            sv.append(["chr1", str(180000 + 500 * k), str(180000 + 500 * k + 60), f"svNull{k}", "DEL", "-60"] + gts(v))
        (self.tmp / "sv.txt").write_text("".join("\t".join(r) + "\n" for r in sv))
        tr = [["chr1", "405000", "405090", "chr1_405000_405090_CAG", "CAG", "."] + ["."] + [f"{int(x)},{int(x)}" for x in self.tr_len[::-1]]]
        (self.tmp / "tr.txt").write_text("".join("\t".join(r) + "\n" for r in tr))

    def run_shard(self, prefix="s1", min_members=40, cov="cov.tsv"):
        run_script("mqtl_shard.py", ["matrix.tsv.gz", cov, "shard.bed", "sv.txt", "samples.txt", "tr.txt", "samples.txt", prefix,
                                     "100000", "10", str(min_members), "0.8", "2.0"], self.tmp)
        pairs = read_tsv(self.tmp / f"{prefix}.pairs.tsv")
        counts = {(r["metric"], r["key"]): int(r["value"]) for r in read_tsv(self.tmp / f"{prefix}.counts.tsv")}
        return pairs, counts

    def test_planted_effects_nulls_and_what_is_left_out(self):
        pairs, counts = self.run_shard()
        by = {(r["region_label"], r["variant_id"]): r for r in pairs}
        hit = by[("CpG:_1", "svEffect")]
        self.assertLess(float(hit["p"]), 1e-20)
        self.assertAlmostEqual(float(hit["effect"]), 10.0, delta=1.5)
        self.assertEqual((hit["effect_unit"], hit["gene_symbol"], hit["region_type"], hit["variant_type"]), ("points per allele", "GENEA", "promoter", "DEL"))
        self.assertEqual((hit["distance"], hit["n"], hit["variant_len"]), ("400", str(self.N), "300"))
        self.assertEqual(int(hit["n_minor"]), sum(1 for x in self.sv_effect if x > 0))
        # twelve columns named, nine in the model: the fixture's monocyte, B and NK fractions are the same for everybody
        # and add nothing beside the intercept
        self.assertEqual(int(hit["df"]), self.N - 9 - 1)
        tr = by[("CpG:_4", "chr1_405000_405090_CAG")]
        self.assertLess(float(tr["p"]), 1e-10)
        self.assertAlmostEqual(float(tr["effect"]), 5.0, delta=1.0)           # 0.05 points per bp
        self.assertEqual(tr["effect_unit"], "points per 100 bp")
        # the island that follows an ancestry component: its variant is no hit while the components are in the model
        self.assertGreater(float(by[("CpG:_3", "svConfounded")]["p"]), 1e-3)
        # nulls: 60 variants against the null island, p-values spread over the unit interval
        null_p = [float(r["p"]) for (label, vid), r in by.items() if label == "CpG:_2" and vid.startswith("svNull")]
        self.assertEqual(len(null_p), 60)
        self.assertGreater(statistics.mean(null_p), 0.35)
        self.assertLess(statistics.mean(null_p), 0.65)
        self.assertLess(sum(1 for p in null_p if p < 0.01), 4)
        # left out: the rare variant, the flat island, the island half the members do not cover, the island outside the shard, chrX
        self.assertNotIn(("CpG:_1", "svRare"), by)
        self.assertFalse(any(label in ("CpG:_5", "CpG:_6", "CpG:_7", "CpG:_8") for label, _ in by))
        self.assertEqual(counts[("islands", "not tested: little spread")], 1)
        self.assertEqual(counts[("islands", "not tested: too few members with coverage")], 1)
        self.assertEqual(counts[("islands", "in the shard")], 6)
        self.assertEqual(counts[("islands", "tested")], 4)
        self.assertGreaterEqual(counts[("pairs", "not tested: too few members off the commonest genotype")], 1)
        self.assertEqual((counts[("members", "used")], counts[("members", "cohort")]), (self.N, self.N + 1))
        self.assertEqual(sorted(k for (m, k) in counts if m == "covariate"),
                         sorted(["intercept", "sex", "PC1", "PC2", "PC3", "PC4", "PC5", "granulocyte", "monocyte", "T-cell", "B-cell", "NK-cell"]))
        self.assertEqual(counts[("variants_tested", "TR")], 1)
        self.assertEqual(counts[("pairs", "tested")], len(pairs))

    def test_without_the_components_the_confounded_variant_looks_like_a_hit(self):
        """One member without ancestry: the components leave the model for everybody, and the island that follows PC1
        is now associated with the variant that follows PC1 too. The page's inflation factor is there for this."""
        lines = (self.tmp / "cov.tsv").read_text().splitlines()
        c = lines[1].split("\t")
        c[4] = ""
        lines[1] = "\t".join(c)
        (self.tmp / "cov2.tsv").write_text("\n".join(lines) + "\n")
        pairs, counts = self.run_shard(prefix="s2", cov="cov2.tsv")
        by = {(r["region_label"], r["variant_id"]): r for r in pairs}
        self.assertLess(float(by[("CpG:_3", "svConfounded")]["p"]), 1e-6)
        self.assertNotIn(("covariate", "PC1"), counts)
        self.assertIn(("covariate", "granulocyte"), counts)
        self.assertEqual(int(by[("CpG:_1", "svEffect")]["df"]), self.N - 4 - 1)     # intercept, sex, granulocyte, T-cell

    def test_too_few_members_writes_headers_and_a_note(self):
        pairs, counts = self.run_shard(prefix="s3", min_members=500)
        self.assertEqual(pairs, [])
        self.assertEqual(counts[("note", "fewer than 500 members with a profile: nothing tested")], 1)
        self.assertEqual((self.tmp / "s3.pairs.tsv").read_text().split("\t")[0], "chrom")


class MatrixAndGatherTest(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def profile(self, sid: str, rows: list[tuple[str, str, str, str]]) -> str:
        """rows: (label, background, mean, coverage) for two islands' worth of lines."""
        head = ("chrom\tstart\tend\tregion_label\tsummary_label\tcompare_label\tbackground_category\tcategory_pop_count\tcategory_pop_freq\t"
                "asm_fishers_pvalue\tmean_hap1_methyl\tmean_hap2_methyl\tmean_meth_delta\tmean_abs_meth_delta_zscore\tmean_combined_methyl\t"
                "mean_combined_methyl_delta\tmean_combined_methyl_zscore\tnum_phased_sites\tnum_partial_sites\tnum_unphased_sites\tmedian_total_coverage\t"
                "median_hap1_coverage\tmedian_hap2_coverage\thgnc_id\tgene_symbol\tregion_type")
        lines = ["##methbat_version=1.1.0", "##datetime=x", head]
        coords = {"CpG:_1": ("chr1", "100", "200", "HGNC:1", "GENEA", "promoter"), "CpG:_2": ("chr2", "500", "900", "", "", "intergenic")}
        for label, background, mean, cov in rows:
            c, s, e, hg, gene, rt = coords[label]
            lines.append("\t".join([c, s, e, label, "Methylated", "", background, "", "", "", "", "", "", "", mean, "", "", "5", "0", "0", cov, "", "", hg, gene, rt]))
        p = self.tmp / f"{sid}.REF.methbat.scored.annotated.tsv"
        p.write_text("\n".join(lines) + "\n")
        return str(p)

    def test_matrix_and_covariates(self):
        profiles = [self.profile("A.1", [("CpG:_1", "ALL", "80.25", "30"), ("CpG:_1", "XX", "80.25", "30"), ("CpG:_2", "ALL", "10.0", "5")]),
                    self.profile("B", [("CpG:_1", "ALL", "60.0", "12"), ("CpG:_2", "ALL", "", "0")]),
                    self.profile("C", [("CpG:_1", "ALL", "55.0", "40"), ("CpG:_2", "ALL", "20.0", "40")])]
        (self.tmp / "profiles.txt").write_text("\n".join(profiles) + "\n")
        (self.tmp / "sexes.tsv").write_text("A.1\tXX\nB\tXY\nC\t\nD\tXX\n")
        (self.tmp / "deconv.tsv").write_text("sample_id\tatlas_md5\tgranulocyte\tmonocyte\tT-cell\tB-cell\tNK-cell\tneuron\nA.1\tx\t0.5\t0.1\t0.2\t0.05\t0.03\t0\nB\tx\t0.6\t0.1\t0.2\t0.0\t0.0\t0\n")
        (self.tmp / "anc.tsv").write_text("#sample_id\tpredicted_ancestry\tPC1\tPC2\tPC3\tPC4\tPC5\nA.1\tEUR\t1.5\t2\t3\t4\t5\nNA1\tEUR\t9\t9\t9\t9\t9\n")
        (self.tmp / "exclude.txt").write_text("C\n")
        run_script("mqtl_matrix.py", ["profiles.txt", "REF", "sexes.tsv", "deconv.tsv", "anc.tsv", "exclude.txt", "8", "out"], self.tmp)
        m = read_tsv(self.tmp / "out.matrix.tsv.gz")
        self.assertEqual(list(m[0])[:9], ["chrom", "start", "end", "region_label", "gene_symbol", "region_type", "n_members", "mean", "sd"])
        self.assertEqual(list(m[0])[9:], ["A.1", "B"])                         # C is excluded, D has no profile
        self.assertEqual((m[0]["region_label"], m[0]["gene_symbol"], m[0]["A.1"], m[0]["B"], m[0]["n_members"], m[0]["mean"]),
                         ("CpG:_1", "GENEA", "80.2", "60.0", "2", "70.10"))
        self.assertEqual((m[1]["A.1"], m[1]["B"], m[1]["n_members"]), ("", "", "0"))   # coverage 5 is below 8; B has no value
        cov = {r["sample_id"]: r for r in read_tsv(self.tmp / "out.covariates.tsv")}
        self.assertEqual({k: (v["used"], v["why_not"]) for k, v in cov.items()},
                         {"A.1": ("yes", ""), "B": ("yes", ""), "C": ("no", "excluded (exclude_ids)"), "D": ("no", "no scored profile in cohort_methyl")})
        self.assertEqual((cov["A.1"]["PC1"], cov["A.1"]["granulocyte"], cov["A.1"]["sex"], cov["B"]["PC1"], cov["B"]["T-cell"]),
                         ("1.5000", "0.5000", "XX", "", "0.2000"))
        # no ancestry and no deconvolution at all: the columns stay empty
        run_script("mqtl_matrix.py", ["profiles.txt", "REF", "sexes.tsv", "", "", "exclude.txt", "8", "out2"], self.tmp)
        self.assertEqual(read_tsv(self.tmp / "out2.covariates.tsv")[0]["PC1"], "")

    def test_gather_corrects_and_summarises(self):
        cols = ["chrom", "start", "end", "region_label", "gene_symbol", "region_type", "variant_type", "variant_id", "variant_start", "variant_end",
                "variant_len", "distance", "n", "n_minor", "r", "t", "df", "p", "effect", "effect_unit"]

        def row(chrom, start, vid, vtype, t, p, rtype="promoter"):
            return "\t".join([chrom, str(start), str(start + 100), "CpG", "G", rtype, vtype, vid, str(start + 500), str(start + 600), "100", "400", "100", "20",
                              "0.1", str(t), "90", str(p), "1.0", "points per allele"])

        (self.tmp / "a.pairs.tsv").write_text("\t".join(cols) + "\n" + row("chr2", 100, "v1", "DEL", 5.0, 0.001) + "\n" + row("chr10", 100, "v2", "TR", 1.0, 0.04) + "\n")
        (self.tmp / "b.pairs.tsv").write_text("\t".join(cols) + "\n" + row("chr1", 900, "v3", "INS", 2.0, 0.02, "gene_body") + "\n" + row("chr1", 100, "v4", "DEL", 0.5, 0.8) + "\n")
        (self.tmp / "a.counts.tsv").write_text("metric\tkey\tvalue\nislands\ttested\t2\nmembers\tused\t100\ncovariate\tsex\t1\npairs\ttested\t2\n")
        (self.tmp / "b.counts.tsv").write_text("metric\tkey\tvalue\nislands\ttested\t3\nmembers\tused\t100\ncovariate\tsex\t1\npairs\ttested\t2\n")
        (self.tmp / "pairs.txt").write_text("a.pairs.tsv\nb.pairs.tsv\n")
        (self.tmp / "counts.txt").write_text("a.counts.tsv\nb.counts.tsv\n")
        run_script("mqtl_gather.py", ["pairs.txt", "counts.txt", "0.05", "window_bp=100000;min_carriers=10", "out"], self.tmp)
        pairs = read_tsv(self.tmp / "out.pairs.tsv.gz")
        self.assertEqual([r["variant_id"] for r in pairs], ["v4", "v3", "v1", "v2"])          # genome order: chr1, chr2, chr10
        q = {r["variant_id"]: float(r["q"]) for r in pairs}
        # p = 0.001, 0.02, 0.04, 0.8 with m = 4: 0.004, 0.04, 0.0533, 0.8
        self.assertAlmostEqual(q["v1"], 0.004)
        self.assertAlmostEqual(q["v3"], 0.04)
        self.assertAlmostEqual(q["v2"], 0.04 * 4 / 3, places=4)
        self.assertAlmostEqual(q["v4"], 0.8)
        self.assertEqual([r["variant_id"] for r in read_tsv(self.tmp / "out.top.tsv")], ["v1", "v3"])
        s = {(r["metric"], r["key"]): r["value"] for r in read_tsv(self.tmp / "out.summary.tsv")}
        self.assertEqual((s[("pairs", "total")], s[("hits", "total")], s[("hits", "islands")], s[("islands", "tested")], s[("members", "used")]),
                         ("4", "2", "2", "5", "100"))
        self.assertEqual((s[("hits_by_variant_type", "DEL")], s[("hits_by_variant_type", "INS")], s[("hits_by_region_type", "gene_body")]), ("1", "1", "1"))
        self.assertAlmostEqual(float(s[("inflation", "lambda")]), ((1.0 + 4.0) / 2) / 0.4549364, places=2)
        self.assertEqual((s[("parameters", "window_bp")], s[("parameters", "fdr")], s[("covariate", "sex")]), ("100000", "0.05", "1"))
        # no pairs at all: headers, zero hits, no inflation row
        (self.tmp / "e.pairs.tsv").write_text("\t".join(cols) + "\n")
        (self.tmp / "e.counts.tsv").write_text("metric\tkey\tvalue\nnote\tfewer than 40 members with a profile: nothing tested\t1\nmembers\tused\t4\n")
        (self.tmp / "pairs2.txt").write_text("e.pairs.tsv\ne.pairs.tsv\n")
        (self.tmp / "counts2.txt").write_text("e.counts.tsv\ne.counts.tsv\n")
        run_script("mqtl_gather.py", ["pairs2.txt", "counts2.txt", "0.05", "", "out2"], self.tmp)
        s2 = {(r["metric"], r["key"]): r["value"] for r in read_tsv(self.tmp / "out2.summary.tsv")}
        self.assertEqual((s2[("pairs", "total")], s2[("hits", "total")], s2[("note", "fewer than 40 members with a profile: nothing tested")]), ("0", "0", "1"))
        self.assertNotIn(("inflation", "lambda"), s2)
        self.assertEqual(read_tsv(self.tmp / "out2.pairs.tsv.gz"), [])


if __name__ == "__main__":
    unittest.main()
