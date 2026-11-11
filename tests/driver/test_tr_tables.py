"""Tests for the scripts embedded in workflows/ugc_wgw/tr/cohort.wdl (cohort_tr): the per-locus statistics with the
expansion-outlier scan, the screen of the named loci against STRchive's thresholds, and the gather."""
from __future__ import annotations

import statistics
import subprocess
import sys
import tempfile
import textwrap
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
WDL = REPO / "workflows" / "ugc_wgw" / "tr" / "cohort.wdl"


def script(name: str) -> str:
    text = WDL.read_text()
    begin, end = f"# ---- {name} (begin) ----\n", f"# ---- {name} (end) ----\n"
    return textwrap.dedent(text.split(begin, 1)[1].split(end, 1)[0])


def run_script(name: str, args: list[str], cwd: Path, stdin: str = "", ok: bool = True) -> subprocess.CompletedProcess:
    res = subprocess.run([sys.executable, "-c", script(name), *args], capture_output=True, text=True, cwd=cwd, input=stdin)
    if ok and res.returncode != 0:
        raise AssertionError(res.stderr)
    return res


def read_tsv(path: Path) -> list[dict[str, str]]:
    lines = [line for line in path.read_text().splitlines() if not line.startswith("##")]
    cols = lines[0].split("\t")
    return [dict(zip(cols, line.split("\t"))) for line in lines[1:]]


LOCI_HEAD = ("trid\tgene\tdisease_id\tdisease\tinheritance\tchrom\tstart\tend\tmotif_len\tbenign_min\tbenign_max\tintermediate_min\tintermediate_max\t"
             "pathogenic_min\tpathogenic_max\treference_motifs\tpathogenic_motifs\tcatalog_motifs\tpathogenic_catalog_motifs\tmotif_match\tframe_differs\t"
             "twin_trids\tlocus_structure\n")


def locus_row(trid: str, inheritance: str, chrom: str, start: int, end: int, motif_len: int, thresholds: tuple, catalog: str, pathogenic: str,
              twin: str = "") -> str:
    bmin, bmax, imin, imax, pmin, pmax = ("" if t is None else str(t) for t in thresholds)
    gene = trid.split("_", 1)[1]
    return "\t".join([trid, gene, trid.split("_")[0], f"{gene} disease", inheritance, chrom, str(start), str(end), str(motif_len), bmin, bmax, imin, imax,
                      pmin, pmax, pathogenic.split(",")[0], pathogenic, catalog, pathogenic, "exact", "no", twin, "(X)*"]) + "\n"


class LocusStatsTest(unittest.TestCase):
    """tr_locus_stats.py: stdin is bcftools query's `CHROM POS END TRID MOTIFS [AL x N] [SD x N]`."""

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.samples = [f"S{i:02d}" for i in range(1, 25)]                 # 24 members: S01-S12 XY, S13-S24 XX
        (self.tmp / "samples.txt").write_text("\n".join(self.samples) + "\n")
        (self.tmp / "sexes.tsv").write_text("".join(f"{s}\t{'XY' if i < 12 else 'XX'}\n" for i, s in enumerate(self.samples)))
        (self.tmp / "overlaps.tsv").write_text("##header\nordinary_trid\tnamed_trid\toverlap_bp\nchr4_90_260_GCC\tHD_HTT\t100\n")

    def line(self, chrom: str, pos: int, end: int, trid: str, motifs: str, als: list[str], sds: list[str] | None = None) -> str:
        sds = sds or ["." if a == "." else ",".join("9" for _ in a.split(",")) for a in als]
        return "\t".join([chrom, str(pos), str(end), trid, motifs, *als, *sds]) + "\n"

    def run_stats(self, lines: list[str], *, min_bp=50, mad_min=5, min_called=20, fraction=0.01, carriers_min=2, ok=True):
        res = run_script("tr_locus_stats.py", ["samples.txt", "sexes.tsv", "overlaps.tsv", "out", str(min_bp), str(mad_min), str(min_called), str(fraction),
                                               str(carriers_min)], self.tmp, stdin="".join(lines), ok=ok)
        if not ok:
            return res, {}
        return res, {k: read_tsv(self.tmp / f"out.{k}.tsv") for k in ("loci", "spectrum", "outliers", "members")}

    def test_distribution_and_spectrum(self):
        als = ["30,30"] * 10 + ["30,36"] * 6 + ["33,39"] * 4 + ["."] * 4
        _res, t = self.run_stats([self.line("chr1", 100, 130, "chr1_100_130_CAG", "CAG", als)])
        row = t["loci"][0]
        flat = [30] * 20 + [30] * 6 + [36] * 6 + [33] * 4 + [39] * 4
        self.assertEqual((row["n_called"], row["n_alleles"], row["n_distinct"], row["n_expected"], row["ref_len"]), ("20", "40", "4", "24", "30"))
        self.assertEqual(float(row["al_median"]), statistics.median(flat))
        self.assertEqual(float(row["al_mad"]), statistics.median(abs(v - statistics.median(flat)) for v in flat))
        self.assertEqual((row["al_min"], row["al_max"], row["al_p95"], row["outlier_status"], row["named_locus"]), ("30", "39", "39", "none", ""))
        self.assertEqual({(r["al"], r["n_alleles"], r["n_members"]) for r in t["spectrum"]}, {("30", "26", "16"), ("33", "4", "4"), ("36", "6", "6"), ("39", "4", "4")})
        members = {r["sample_id"]: r for r in t["members"]}
        self.assertEqual((members["S01"]["n_loci"], members["S01"]["n_loci_uncalled"], members["S24"]["n_loci_uncalled"]), ("1", "0", "1"))
        # an odd number of alleles (haploid calls) and an even split: the median averages the two middle values
        _res, t = self.run_stats([self.line("chrX", 5, 25, "chrX_5_25_A", "A", ["20"] * 12 + ["20,22"] * 12)])
        self.assertEqual((t["loci"][0]["n_alleles"], t["loci"][0]["al_median"]), ("36", "20"))
        _res, t = self.run_stats([self.line("chr1", 5, 25, "chr1_5_25_A", "A", ["20,20"] * 12 + ["22,22"] * 12)])
        self.assertEqual(t["loci"][0]["al_median"], "21")

    def test_outliers_scale_and_rarity(self):
        base = ["300,300"] * 23
        # MAD 0 at a 60 bp motif: one extra copy (60 bp) is over the bp floor but is one repeat unit, not five
        _res, t = self.run_stats([self.line("chr1", 1000, 1300, "chr1_1000_1300_" + "A" * 60, "A" * 60, base + ["300,360"])])
        self.assertEqual((t["loci"][0]["outlier_status"], t["loci"][0]["n_carriers"], t["outliers"]), ("none", "0", []))
        # six extra copies: listed, with the other allele, the spanning reads of the long one, the excess in units and in scaled MADs
        _res, t = self.run_stats([self.line("chr1", 1000, 1300, "chr1_1000_1300_" + "A" * 60, "A" * 60, base + ["300,660"],
                                            sds=["9,9"] * 23 + ["11,4"])])
        (pair,) = t["outliers"]
        self.assertEqual((pair["sample_id"], pair["al"], pair["al_other"], pair["sd"], pair["excess_bp"], pair["excess_units"], pair["z"], pair["both_alleles"],
                          pair["n_carriers"], pair["cohort_median"], pair["scale"], pair["ref_len"]),
                         ("S24", "660", "300", "4", "360", "6", "6", "no", "1", "300", "60", "300"))
        self.assertEqual(t["loci"][0]["outlier_status"], "listed")
        self.assertEqual({r["sample_id"]: r["n_outliers"] for r in t["members"]}["S24"], "1")
        # a homopolymer with MAD 0: 49 bp over the median is under the floor, 50 is listed
        for extra, status in ((49, "none"), (50, "listed")):
            _res, t = self.run_stats([self.line("chr1", 10, 30, "chr1_10_30_A", "A", ["20,20"] * 23 + [f"20,{20 + extra}"])])
            self.assertEqual(t["loci"][0]["outlier_status"], status, extra)
        # both alleles long: one row, both_alleles yes, both spanning-read counts
        _res, t = self.run_stats([self.line("chr1", 10, 30, "chr1_10_30_A", "A", ["20,20"] * 23 + ["90,95"], sds=["9,9"] * 23 + ["5,6"])])
        self.assertEqual((t["outliers"][0]["al"], t["outliers"][0]["al_other"], t["outliers"][0]["both_alleles"], t["outliers"][0]["sd"]), ("95", "90", "yes", "5,6"))
        self.assertEqual({r["sample_id"]: r["n_outliers_both_alleles"] for r in t["members"]}["S24"], "1")
        # a common long allele is polymorphism: three carriers of 24 exceed max(2, 1% of the called), the locus lists nobody
        _res, t = self.run_stats([self.line("chr1", 10, 30, "chr1_10_30_A", "A", ["20,20"] * 21 + ["20,120"] * 3)])
        self.assertEqual((t["loci"][0]["outlier_status"], t["loci"][0]["n_carriers"], t["loci"][0]["al_max"], t["outliers"]), ("common", "3", "120", []))
        # ... two carriers are listed, each row saying so; a pedigree cohort raises the fraction
        _res, t = self.run_stats([self.line("chr1", 10, 30, "chr1_10_30_A", "A", ["20,20"] * 22 + ["20,120"] * 2)])
        self.assertEqual([(r["sample_id"], r["n_carriers"]) for r in t["outliers"]], [("S23", "2"), ("S24", "2")])
        _res, t = self.run_stats([self.line("chr1", 10, 30, "chr1_10_30_A", "A", ["20,20"] * 21 + ["20,120"] * 3)], fraction=0.2)
        self.assertEqual(t["loci"][0]["outlier_status"], "listed")
        # too few called members: not scanned, the carriers still counted in the table
        _res, t = self.run_stats([self.line("chr1", 10, 30, "chr1_10_30_A", "A", ["20,20"] * 18 + ["20,120"] + ["."] * 5)])
        self.assertEqual((t["loci"][0]["outlier_status"], t["loci"][0]["n_called"], t["outliers"]), ("below_min_called", "19", []))

    def test_named_loci_haploid_and_empty(self):
        hap = ["40"] * 11 + ["400"] + ["."] * 12                           # chrY: the XY members, haploid; the XX ones are not expected
        lines = [self.line("chrY", 500, 540, "chrY_500_540_AT", "AT", hap),
                 self.line("chr4", 100, 200, "HD_HTT", "CAG,CCG", ["60,60"] * 23 + ["60,400"]),
                 self.line("chr4", 90, 260, "chr4_90_260_GCC", "GCC", ["170,170"] * 23 + ["170,510"]),
                 self.line("chr2", 1, 9, "chr2_1_9_A", "A", ["."] * 24)]
        _res, t = self.run_stats(lines, min_called=10)
        by = {r["trid"]: r for r in t["loci"]}
        self.assertEqual((by["chrY_500_540_AT"]["n_expected"], by["chrY_500_540_AT"]["n_called"], by["chrY_500_540_AT"]["n_alleles"]), ("12", "12", "12"))
        self.assertEqual((by["HD_HTT"]["named_locus"], by["chr4_90_260_GCC"]["named_locus"]), ("HD_HTT", "HD_HTT"))   # the locus and its twin
        self.assertEqual((by["chr2_1_9_A"]["n_called"], by["chr2_1_9_A"]["outlier_status"], by["chr2_1_9_A"]["al_median"]), ("0", "below_min_called", ""))
        pairs = {(r["sample_id"], r["trid"]): r for r in t["outliers"]}
        self.assertEqual(pairs[("S12", "chrY_500_540_AT")]["al_other"], "")
        self.assertEqual(pairs[("S24", "chr4_90_260_GCC")]["named_locus"], "HD_HTT")
        members = {r["sample_id"]: r for r in t["members"]}
        # outliers at a named locus or its twin never feed the QC count
        self.assertEqual((members["S24"]["n_outliers"], members["S24"]["n_outliers_named"], members["S12"]["n_outliers"]), ("0", "2", "1"))
        self.assertEqual((members["S24"]["n_loci"], members["S24"]["n_loci_uncalled"], members["S01"]["n_loci_uncalled"]), ("4", "2", "1"))
        # an empty shard (the contigs without a locus): header-only tables, every sample with zero counts
        _res, t = self.run_stats([])
        self.assertEqual((t["loci"], t["spectrum"], t["outliers"], len(t["members"]), t["members"][0]["n_loci"]), ([], [], [], 24, "0"))
        # a line with the wrong number of fields stops the task
        res, _t = self.run_stats(["chr1\t1\t9\tchr1_1_9_A\tA\t20,20\n"], ok=False)
        self.assertNotEqual(res.returncode, 0)
        self.assertIn("fields for 24 samples", res.stderr)


VCF_HEAD = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
FMT = "GT:AL:ALLR:SD:MC:MS:AP:AM:PS"


def allele(n: int) -> str:
    """An allele sequence of AL n: the pad base and n bases."""
    return "N" + "A" * n


class ScreenTest(unittest.TestCase):
    """tr_screen.py: the named loci of a VCF text against the reference table, LPS re-paired through GT."""

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.members = ["M1", "M2", "M3", "F1", "F2"]
        (self.tmp / "members.txt").write_text("\n".join(self.members) + "\n")
        (self.tmp / "sexes.tsv").write_text("M1\tXY\nM2\tXY\nM3\tXY\nF1\tXX\nF2\tXX\n")
        (self.tmp / "overlaps.tsv").write_text("ordinary_trid\tnamed_trid\toverlap_bp\nchr4_90_260_GCC\tHD_HTT\t100\n")
        rows = [
            locus_row("HD_HTT", "AD", "chr4", 100, 200, 3, (6, 26, 27, 39, 40, 250), "CAG,CCG", "CAG", twin="chr4_90_260_GCC"),
            locus_row("CANVAS_RFC1", "AR", "chr4", 5000, 5060, 5, (0, 11, 11, 200, 400, 2750), "AAGGG,ACAGG,AAAAG", "AAGGG,ACAGG"),
            locus_row("FXS_FMR1", "XD", "chrX", 700, 760, 3, (5, 44, 45, 200, 201, 2000), "CGG", "CGG"),
            locus_row("PRTS_ARX", "XR", "chrX", 900, 960, 3, (None, None, None, None, 20, 30), "NGC", "NGC"),
            locus_row("HMNR7_VWA1", "AR", "chr1", 10, 30, 10, (2, 2, None, None, 1, 3), "GGCGCGGAGC", "GGCGCGGAGC"),
            locus_row("CPEO_POLG", "", "chr15", 300, 330, 3, (10, 10, None, None, None, None), "GCT", "GCT"),
            locus_row("SCA_MISSING", "AD", "chr9", 50, 80, 3, (5, 20, None, None, 40, 100), "CAG", "CAG"),
        ]
        (self.tmp / "loci.tsv").write_text("##tr_map_source=test\n" + LOCI_HEAD + "".join(rows))

    def record(self, chrom: str, pos: int, end: int, trid: str, motifs: str, alts: list[int], fields: dict[str, str], ref_len: int | None = None) -> str:
        ref = allele(end - pos if ref_len is None else ref_len)
        cols = [fields.get(m, ".:.:.:.:.:.:.:.:.") for m in self.vcf_samples]
        return "\t".join([chrom, str(pos), ".", ref, ",".join(allele(a) for a in alts) or ".", ".", ".",
                          f"TRID={trid};END={end};MOTIFS={motifs};STRUC=<TR>", FMT, *cols]) + "\n"

    def run_screen(self, records: list[str], lps: str = "", borderline_bp: int = 6, ok: bool = True):
        (self.tmp / "named.vcf").write_text(VCF_HEAD + "\t".join(self.vcf_samples) + "\n" + "".join(records))
        (self.tmp / "named.lps.tsv").write_text("trid\tmotif\t" + "\t".join(self.vcf_samples) + "\n" + lps)
        res = run_script("tr_screen.py", ["screen", "loci.tsv", "overlaps.tsv", "named.vcf", "named.lps.tsv", "members.txt", "sexes.tsv", "out", str(borderline_bp)],
                         self.tmp, ok=ok)
        if not ok:
            return res, {}, {}
        screen = {(r["sample_id"], r["trid"]): r for r in read_tsv(self.tmp / "out.screen.tsv")}
        return res, screen, {r["trid"]: r for r in read_tsv(self.tmp / "out.screen_loci.tsv")}

    vcf_samples = ["F2", "M1", "M2", "M3", "F1"]                            # not the cohort order: samples are matched by name

    def test_regions(self):
        res = run_script("tr_screen.py", ["regions", "loci.tsv", "overlaps.tsv"], self.tmp)
        rows = [line.split("\t") for line in res.stdout.splitlines()]
        self.assertIn(["chr4", "99", "200"], rows)                          # POS is the base before the repeat
        self.assertIn(["chr4", "89", "260"], rows)                          # the twin of HD_HTT
        self.assertEqual(len(rows), 8)

    def test_pairing_classes_and_inheritance(self):
        htt = self.record("chr4", 100, 200, "HD_HTT", "CAG,CCG", [60, 150, 66], {
            # 2|1: FORMAT values are in allele-index order (1 then 2), the LPS in GT order (2 then 1)
            "M1": "2|1:60,150:60-60,148-153:20,12:17_3,45_5:.:0.99,0.95:0.1,0.2:7",
            "M2": "1/1:60,60:60-60,60-60:15,14:17_3,17_3:.:0.99,0.99:.,.:.",
            "M3": "3|1:60,66:60-60,66-66:15,14:17_3,19_3:.:0.99,0.99:0.1,0.1:7",
            "F1": "0/2:100,150:99-101,150-150:9,8:29_4,45_5:.:0.98,0.97:0.3,0.4:.",
        })
        twin = self.record("chr4", 90, 260, "chr4_90_260_GCC", "GCC", [200], {"M1": "0/1:170,200:.:9,9:56,66:.:0.9,0.9:.,.:.", "F2": "0/0:170,170:.:9,9:56,56:.:0.9,0.9:.,.:."})
        lps = "HD_HTT\tCAG\t.\t40,17\t17,17\t19,17\t29,45\nHD_HTT\tCCG\t.\t5,3\t3,3\t3,3\t4,5\n"
        res, s, loci = self.run_screen([htt, twin], lps)
        m1 = s[("M1", "HD_HTT")]
        self.assertEqual((m1["gt"], m1["ploidy"], m1["class"], m1["n_alleles_pathogenic"], m1["mode_match"], m1["measure"]), ("2|1", "2", "pathogenic", "1", "consistent", "mc_pathogenic"))
        # a1 is the allele with more copies: 45 CAG on 150 bp, whose longest pure CAG run (GT slot 0) is 40
        self.assertEqual((m1["a1_copies"], m1["a1_class"], m1["a1_al"], m1["a1_lps"], m1["a1_sd"], m1["a1_am"], m1["a1_mc"]), ("45", "pathogenic", "150", "40", "12", "0.2", "45_5"))
        self.assertEqual((m1["a2_copies"], m1["a2_class"], m1["a2_al"], m1["a2_lps"], m1["a2_sd"]), ("17", "normal", "60", "17", "20"))
        self.assertEqual((m1["twin_trid"], m1["twin_al"], m1["notes"]), ("chr4_90_260_GCC", "170,200", ""))
        self.assertEqual((s[("M2", "HD_HTT")]["class"], s[("M2", "HD_HTT")]["a1_lps"], s[("M2", "HD_HTT")]["mode_match"]), ("normal", "17", "none"))
        m3 = s[("M3", "HD_HTT")]                                           # 3|1: position 0 is allele 1 (GT slot 1), position 1 allele 3 (slot 0)
        self.assertEqual((m3["a1_copies"], m3["a1_al"], m3["a1_lps"], m3["a2_copies"], m3["a2_lps"]), ("19", "66", "19", "17", "17"))
        f1 = s[("F1", "HD_HTT")]                                           # 0/2: the reference allele of 100 bp carries 29 CAG: intermediate
        self.assertEqual((f1["a1_class"], f1["a2_class"], f1["a2_copies"], f1["class"], f1["a2_lps"]), ("pathogenic", "intermediate", "29", "pathogenic", "29"))
        f2 = s[("F2", "HD_HTT")]
        self.assertEqual((f2["class"], f2["gt"], f2["ploidy"], f2["twin_al"], f2["a1_copies"]), ("uncalled", ".", "0", "170,170", ""))
        self.assertIn("0 allele-length mismatches", res.stderr)
        row = loci["HD_HTT"]
        self.assertEqual((row["n_expected"], row["n_called"], row["call_rate"], row["n_alleles"], row["n_pathogenic"], row["n_normal"], row["copies_max"],
                          row["call_note"], row["twin_trid"]), ("5", "4", "0.8", "8", "2", "2", "45", "threshold_or_frame_suspect", "chr4_90_260_GCC"))
        self.assertEqual((row["am_median_xy"], row["am_median_xx"]), ("0.1", "0.35"))
        self.assertEqual((loci["SCA_MISSING"]["n_called"], loci["SCA_MISSING"]["call_note"]), ("0", "mostly_uncalled"))   # no record at all
        self.assertEqual(s[("M1", "SCA_MISSING")]["class"], "uncalled")

    def test_motif_replacement_and_recessive(self):
        def rfc1(alts, fields):
            return self.record("chr4", 5000, 5060, "CANVAS_RFC1", "AAGGG,ACAGG,AAAAG", alts, fields)
        rec = rfc1([2500, 2100, 2050], {
            "M1": "1/1:2500,2500:.:9,9:0_0_500,0_0_500:.:0.9,0.9:.,.:.",      # long on the benign AAAAG only
            "M2": "2/3:2100,2050:.:9,9:420_0_0,0_410_0:.:0.9,0.9:.,.:.",      # AAGGG on one allele, ACAGG on the other
            "M3": "0/2:60,2100:.:9,9:0_0_12,420_0_0:.:0.9,0.9:.,.:.",         # one pathogenic allele
            "F1": "2/2:2100,2100:.:9,9:420_0_0,420_0_0:.:0.9,0.9:.,.:.",
        })
        _res, s, loci = self.run_screen([rec])
        m1 = s[("M1", "CANVAS_RFC1")]
        self.assertEqual((m1["class"], m1["a1_copies"], m1["n_alleles_pathogenic"], m1["mode_match"]), ("normal", "0", "0", "none"))
        self.assertEqual(m1["notes"], "lps_missing;expanded_other_motif;hom_call")
        m2 = s[("M2", "CANVAS_RFC1")]
        self.assertEqual((m2["class"], m2["a1_copies"], m2["a2_copies"], m2["n_alleles_pathogenic"], m2["mode_match"]), ("pathogenic", "420", "410", "2", "consistent"))
        m3 = s[("M3", "CANVAS_RFC1")]
        self.assertEqual((m3["class"], m3["n_alleles_pathogenic"], m3["mode_match"], m3["a2_class"]), ("pathogenic", "1", "carrier", "normal"))
        self.assertEqual((s[("F1", "CANVAS_RFC1")]["mode_match"], "hom_call" in s[("F1", "CANVAS_RFC1")]["notes"]), ("consistent", True))
        self.assertEqual(loci["CANVAS_RFC1"]["n_expanded_other_motif"], "2")
        # a reference-length allele is never an expansion of another motif, however many units its frame holds
        # (CJD_PRNP on real data: 123 bp = one 27-mer and four 24-mers, pathogenic from five)
        ref = rfc1([], {"M1": "0/0:60,60:.:9,9:0_0_12,0_0_12:.:0.9,0.9:.,.:."})
        table = self.tmp / "loci.tsv"
        table.write_text(table.read_text().replace("CANVAS_RFC1\tRFC1\tCANVAS\tRFC1 disease\tAR\tchr4\t5000\t5060\t5\t0\t11\t11\t200\t400",
                                                   "CANVAS_RFC1\tRFC1\tCANVAS\tRFC1 disease\tAR\tchr4\t5000\t5060\t5\t0\t11\t11\t200\t12"))
        _res, s, loci = self.run_screen([ref])
        self.assertEqual((s[("M1", "CANVAS_RFC1")]["class"], s[("M1", "CANVAS_RFC1")]["notes"], loci["CANVAS_RFC1"]["n_expanded_other_motif"]),
                         ("normal", "lps_missing;hom_call", "0"))

    def test_x_linked_nulls_and_nested_ranges(self):
        fmr1 = self.record("chrX", 700, 760, "FXS_FMR1", "CGG", [165, 630, 132], {
            "M1": "1:165:165-165:20:55:.:0.99:0.9:.",                         # haploid XY, 55 CGG: a premutation, never normal
            "M2": "2:630:600-650:9:210:.:0.9:0.95:.",                         # full mutation
            "F1": "0/3:60,132:60-60,132-132:9,9:20,44:.:0.99,0.99:0.5,0.5:.", # 44 CGG, one from the threshold
            "F2": "1:165:165-165:20:55:.:0.99:0.9:.",                         # an XX member called haploid
            "M3": "0/1:60,165:60-60,165-165:9,9:20,55:.:0.99,0.99:0.1,0.1:.", # an XY member called diploid
        })
        arx = self.record("chrX", 900, 960, "PRTS_ARX", "NGC", [75, 45], {"M1": "1:75:75-75:9:25:.:0.9:.:.", "F1": "0/1:60,75:60-60,75-75:9,9:15,25:.:0.9,0.9:.,.:.",
                                                                           "F2": "1/1:75,75:75-75,75-75:9,9:25,25:.:0.9,0.9:.,.:.", "M2": "2:45:45-45:9:15:.:0.9:.:."})
        vwa1 = self.record("chr1", 10, 30, "HMNR7_VWA1", "GGCGCGGAGC", [10, 30], {"M1": "0/0:20,20:20-20,20-20:9,9:2,2:.:0.9,0.9:.,.:.",
                                                                                 "M2": "0/1:20,10:20-20,10-10:9,9:2,1:.:0.9,0.9:.,.:.",
                                                                                 "F1": "2/2:30,30:30-30,30-30:9,9:3,3:.:0.9,0.9:.,.:."})
        polg = self.record("chr15", 300, 330, "CPEO_POLG", "GCT", [33, 27], {"M1": "0/1:30,33:30-30,33-33:9,9:10,11:.:0.9,0.9:.,.:.",
                                                                           "M2": "0/2:30,27:30-30,27-27:9,9:10,9:.:0.9,0.9:.,.:."})
        _res, s, loci = self.run_screen([fmr1, arx, vwa1, polg])   # no LPS table rows here: every called row says lps_missing

        def notes(member: str, trid: str) -> list[str]:
            return [n for n in s[(member, trid)]["notes"].split(";") if n and n != "lps_missing"]

        self.assertEqual((s[("M1", "FXS_FMR1")]["class"], s[("M1", "FXS_FMR1")]["ploidy"], s[("M1", "FXS_FMR1")]["notes"]), ("intermediate", "1", "lps_missing"))
        m2 = s[("M2", "FXS_FMR1")]
        self.assertEqual((m2["class"], m2["mode_match"], m2["borderline"], m2["a1_am"]), ("pathogenic", "consistent", "yes", "0.95"))   # ALLR 600-650 spans 201 CGG
        f1 = s[("F1", "FXS_FMR1")]
        self.assertEqual((f1["a1_copies"], f1["class"], f1["borderline"]), ("44", "normal", "yes"))   # within 6 bp of the intermediate range
        self.assertEqual((notes("F2", "FXS_FMR1"), notes("M3", "FXS_FMR1")), (["ploidy_sex_mismatch"], ["ploidy_sex_mismatch"]))
        self.assertEqual((loci["FXS_FMR1"]["am_median_xy"], loci["FXS_FMR1"]["am_median_xx"]), ("0.5", "0.5"))
        # XR: a haploid pathogenic allele is consistent, one of two a carrier, two consistent; no benign range: below_pathogenic
        self.assertEqual([s[(m, "PRTS_ARX")]["mode_match"] for m in ("M1", "F1", "F2")], ["consistent", "carrier", "consistent"])
        self.assertEqual((s[("M2", "PRTS_ARX")]["class"], s[("F1", "PRTS_ARX")]["a2_class"], s[("M1", "PRTS_ARX")]["measure"]), ("below_pathogenic", "below_pathogenic", "mc_pathogenic"))   # a1 is the longer allele
        self.assertIn("hom_call", s[("F2", "PRTS_ARX")]["notes"])
        # the benign value sits inside the pathogenic range: 2 copies are normal, 1 and 3 pathogenic
        self.assertEqual([s[(m, "HMNR7_VWA1")]["class"] for m in ("M1", "M2", "F1")], ["normal", "pathogenic", "pathogenic"])
        self.assertEqual((s[("M2", "HMNR7_VWA1")]["mode_match"], s[("F1", "HMNR7_VWA1")]["mode_match"], s[("M1", "HMNR7_VWA1")]["borderline"]), ("carrier", "consistent", "no"))
        # no pathogenic threshold at all: only the benign range is known
        self.assertEqual((s[("M1", "CPEO_POLG")]["class"], s[("M2", "CPEO_POLG")]["class"], s[("M1", "CPEO_POLG")]["mode_match"]), ("above_benign", "below_benign", "none"))
        self.assertEqual(loci["HMNR7_VWA1"]["call_note"], "threshold_or_frame_suspect")

    def test_fallbacks_and_refusals(self):
        # the record's motifs are not the table's: length over motif length, with the note
        changed = self.record("chr4", 100, 200, "HD_HTT", "CAG,CCG,CAA", [150], {"M1": "0/1:100,150:.:9,9:20_5_1,40_5_1:.:0.9,0.9:.,.:."})
        # an allele whose sequence is not AL long is counted
        bad = self.record("chr9", 50, 80, "SCA_MISSING", "CAG", [90], {"M1": "0/1:30,95:.:9,9:10,31:.:0.9,0.9:.,.:."})
        res, s, _loci = self.run_screen([changed, bad])
        m1 = s[("M1", "HD_HTT")]
        self.assertEqual((m1["measure"], m1["a1_copies"], m1["class"], m1["a1_lps"], m1["notes"]), ("al_div_motif_len", "50", "pathogenic", "", "catalog_motifs_changed"))
        self.assertEqual(s[("F1", "HD_HTT")]["notes"], "catalog_motifs_changed")
        self.assertIn("al_seq_mismatch", s[("M1", "SCA_MISSING")]["notes"])
        self.assertIn("1 allele-length mismatches", res.stderr)
        self.assertEqual(s[("M1", "SCA_MISSING")]["a1_class"], "above_benign")   # 31 copies: between benign_max 20 and pathogenic_min 40
        # a cohort member absent from the VCF stops the task
        (self.tmp / "members.txt").write_text("M1\nGHOST\n")
        res, _s, _l = self.run_screen([changed], ok=False)
        self.assertNotEqual(res.returncode, 0)
        self.assertIn("cohort members absent from the TRGT VCF: GHOST", res.stderr)


class GatherTest(unittest.TestCase):
    def test_tables(self):
        tmp = Path(tempfile.mkdtemp())
        loci_head = "chrom\tstart\tend\ttrid\tmotifs\tnamed_locus\tref_len\tn_expected\tn_called\tn_alleles\tn_distinct\tal_min\tal_median\tal_mad\tal_p95\tal_p99\tal_max\tn_carriers\toutlier_status\n"
        pair_head = ("sample_id\ttrid\tchrom\tstart\tend\tmotifs\tnamed_locus\tal\tal_other\tsd\texcess_bp\texcess_units\tz\tcohort_median\tcohort_mad\tscale\tref_len\t"
                     "n_called\tn_carriers\tboth_alleles\n")
        def pair(sample, trid, z, excess, named=""):
            return f"{sample}\t{trid}\tchr1\t1\t9\tA\t{named}\t99\t20\t9\t{excess}\t{excess}\t{z}\t20\t0\t1\t8\t30\t1\tno\n"
        (tmp / "a.loci.tsv").write_text(loci_head + "chr1\t1\t9\tL1\tA\t\t8\t3\t3\t6\t2\t8\t8\t0\t9\t9\t99\t1\tlisted\nchr1\t20\t29\tL2\tA\t\t9\t3\t3\t6\t1\t9\t9\t0\t9\t9\t9\t0\tnone\n")
        (tmp / "b.loci.tsv").write_text(loci_head)                              # an empty shard
        (tmp / "c.loci.tsv").write_text(loci_head + "chr2\t1\t9\tL3\tA\t\t8\t3\t3\t6\t2\t8\t8\t0\t9\t9\t99\t3\tcommon\n")
        for s in "abc":
            (tmp / f"{s}.spectrum.tsv").write_text("trid\tal\tn_alleles\tn_members\n" + ("L1\t8\t5\t3\nL1\t99\t1\t1\n" if s == "a" else ""))
        (tmp / "a.outliers.tsv").write_text(pair_head + pair("S2", "L1", 91, 91) + pair("S1", "L9", 7, 60) + pair("S1", "L8", 9, 70) + pair("S1", "L7", 9, 80, "HD_HTT"))
        (tmp / "b.outliers.tsv").write_text(pair_head)
        (tmp / "c.outliers.tsv").write_text(pair_head)
        m_head = "sample_id\tn_loci\tn_loci_uncalled\tn_outliers\tn_outliers_named\tn_outliers_both_alleles\n"
        (tmp / "a.members.tsv").write_text(m_head + "S1\t2\t0\t2\t1\t0\nS2\t2\t1\t1\t0\t0\nS3\t2\t0\t0\t0\t0\n")
        (tmp / "b.members.tsv").write_text(m_head + "S1\t0\t0\t0\t0\t0\nS2\t0\t0\t0\t0\t0\nS3\t0\t0\t0\t0\t0\n")
        (tmp / "c.members.tsv").write_text(m_head + "S1\t1\t1\t0\t0\t0\nS2\t1\t0\t0\t0\t0\nS3\t1\t0\t0\t0\t0\n")
        for kind in ("loci", "spectrum", "outliers", "members"):
            (tmp / f"{kind}.list").write_text("".join(f"{tmp}/{s}.{kind}.tsv\n" for s in "abc"))
        (tmp / "screen.tsv").write_text("sample_id\ttrid\tclass\tnotes\nS1\tHD_HTT\tpathogenic\t\nS1\tFXS_FMR1\tuncalled\t\nS2\tHD_HTT\tnormal\t\nS2\tFXS_FMR1\tintermediate\tploidy_sex_mismatch\n"
                                        "S3\tHD_HTT\tbelow_pathogenic\t\nS3\tFXS_FMR1\tnormal\t\n")
        (tmp / "screen_loci.tsv").write_text("trid\tcall_note\nHD_HTT\tok\nFXS_FMR1\tfew_uncalled\n")
        (tmp / "samples.txt").write_text("S1\nS2\nS3\n")
        (tmp / "sexes.tsv").write_text("S1\tXY\nS2\tXX\nS3\t\n")
        (tmp / "loci_ref.tsv").write_text("##tr_map_source=STRchive vX, licensed CC BY 4.0\n##tr_map_input=strchive\tfile\tsha256:abc\n" + LOCI_HEAD)
        res = run_script("tr_gather.py", ["loci.list", "spectrum.list", "outliers.list", "members.list", "screen.tsv", "screen_loci.tsv", "samples.txt", "sexes.tsv",
                                          "loci_ref.tsv", "out", "2", "outlier_min_bp=50;outlier_mad_min=5.0"], tmp)
        self.assertIn("3 members, 3 loci, 4 outlier pairs, 2 spectrum rows", res.stderr)
        self.assertEqual([r["trid"] for r in read_tsv(tmp / "out.loci.tsv")], ["L1", "L2", "L3"])            # shard order
        self.assertEqual(len(read_tsv(tmp / "out.outliers.tsv")), 4)
        top = read_tsv(tmp / "out.outliers_top.tsv")
        self.assertEqual([(r["sample_id"], r["trid"]) for r in top], [("S1", "L7"), ("S1", "L8"), ("S2", "L1")])   # two per member: by z, then excess
        samples = {r["sample_id"]: r for r in read_tsv(tmp / "out.samples.tsv")}
        self.assertEqual((samples["S1"]["sex"], samples["S1"]["n_loci_called"], samples["S1"]["n_loci_uncalled"], samples["S1"]["n_outliers"], samples["S1"]["n_outliers_named"]),
                         ("XY", "2", "1", "2", "1"))
        self.assertEqual((samples["S1"]["n_named_called"], samples["S1"]["n_named_uncalled"], samples["S1"]["n_loci_pathogenic"]), ("1", "1", "1"))
        self.assertEqual((samples["S2"]["n_loci_intermediate"], samples["S2"]["n_ploidy_sex_mismatch"], samples["S3"]["n_loci_above_benign"], samples["S3"]["sex"]), ("1", "1", "1", ""))
        summ = {(r["resource"], r["metric"], r["key"]): r["value"] for r in read_tsv(tmp / "out.summary.tsv")}
        self.assertEqual((summ[("cohort", "n_members", "")], summ[("cohort", "n_loci", "")], summ[("cohort", "n_outlier_pairs", "")]), ("3", "3", "4"))
        self.assertEqual((summ[("loci", "outlier_status", "listed")], summ[("loci", "outlier_status", "common")], summ[("loci", "outlier_status", "below_min_called")]), ("1", "1", "0"))
        self.assertEqual((summ[("parameters", "outlier_mad_min", "")], summ[("screen", "call_note", "few_uncalled")], summ[("screen_pathogenic", "n_members", "HD_HTT")]), ("5.0", "1", "1"))
        self.assertIn("CC BY 4.0", summ[("source", "tr_map_source", "")])
        self.assertEqual((summ[("outlier_loci", "n_carriers", "L1")], summ[("outlier_loci", "named_locus", "L7")]), ("1", "HD_HTT"))
        # a shard whose header differs stops the task
        (tmp / "c.spectrum.tsv").write_text("trid\tal\n")
        res = run_script("tr_gather.py", ["loci.list", "spectrum.list", "outliers.list", "members.list", "screen.tsv", "screen_loci.tsv", "samples.txt", "sexes.tsv",
                                          "loci_ref.tsv", "out", "2", ""], tmp, ok=False)
        self.assertIn("header differs from the first shard's", res.stderr)


if __name__ == "__main__":
    unittest.main()
