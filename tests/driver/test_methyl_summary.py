"""Tests of the Python embedded in workflows/ugc_wgw/methyl/cohort.wdl: the per-member annotation/events script of
ugc_wgw_methbat_profile_bg and the cohort summary of ugc_wgw_methyl_summary, run from the WDL text on synthetic
MethBat 1.1.0 outputs (the real column layouts of 2026-11-06)."""
from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import textwrap
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
WDL = REPO / "workflows" / "ugc_wgw" / "methyl" / "cohort.wdl"


def script(name: str) -> str:
    text = WDL.read_text()
    begin, end = f"# ---- {name} (begin) ----\n", f"# ---- {name} (end) ----\n"
    body = text.split(begin, 1)[1].split(end, 1)[0]
    return textwrap.dedent(body)


SCORED_COLS = ["chrom", "start", "end", "region_label", "summary_label", "compare_label", "background_category", "category_pop_count",
               "category_pop_freq", "asm_fishers_pvalue", "mean_hap1_methyl", "mean_hap2_methyl", "mean_meth_delta", "mean_abs_meth_delta_zscore",
               "mean_combined_methyl", "mean_combined_methyl_delta", "mean_combined_methyl_zscore", "num_phased_sites", "num_partial_sites",
               "num_unphased_sites", "median_total_coverage", "median_hap1_coverage", "median_hap2_coverage"]
BG_COLS = ["chrom", "start", "end", "region_label", "data_category", "num_phased", "num_unphased", "NoData", "Uncategorized", "Methylated",
           "Unmethylated", "AlleleSpecificMethylation", "avg_abs_meth_deltas", "stdev_abs_meth_deltas", "avg_combined_methyls", "stdev_combined_methyls"]
REPORT_COLS = ["chrom", "start", "end", "region_label", "report_summary", "qc_warnings", "expected_category", "summary_label", "mean_combined_methyl",
               "mean_meth_delta", "mean_hap1_methyl", "mean_hap2_methyl", "asm_fishers_pvalue", "num_phased_sites", "num_partial_sites",
               "num_unphased_sites", "median_total_coverage", "median_hap1_coverage", "median_hap2_coverage"]
COMPARE_COLS = ["chrom", "start", "end", "baseline_category", "compare_category", "summary_comparison", "zscore_avg_abs_meth_deltas",
                "delta_avg_abs_meth_deltas", "baseline_num_phased", "compare_num_phased", "zscore_avg_combined_methyls", "delta_avg_combined_methyls",
                "baseline_num_samples", "compare_num_samples"]
REGIONS = [("chr1", "1000", "1500", "CpG:_1"), ("chr1", "5000", "5600", "CpG:_2"), ("chr2", "100", "300", "CpG:_3")]
MAP = "##methbat_map_generator=test\nchrom\tstart\tend\tregion_label\thgnc_id\tgene_symbol\tregion_type\n" \
      "chr1\t1000\t1500\tCpG:_1\tHGNC:1;HGNC:2\tGENEA;GENEB\tpromoter\n" \
      "chr1\t5000\t5600\tCpG:_2\t\t\tintergenic\n" \
      "chr2\t100\t300\tCpG:_3\tHGNC:3\tGENEC\timprinted\n"


def scored_tsv(rows_by_region_label: dict[tuple[str, str], tuple[str, str]]) -> str:
    """rows_by_region_label: (region_label, background label) -> (summary_label, compare_label)."""
    lines = ["##methbat_version=1.1.0-test", "\t".join(SCORED_COLS)]
    for chrom, start, end, label in REGIONS:
        for bg_label in ("ALL", "XX"):
            summary, compare = rows_by_region_label.get((label, bg_label), ("Methylated", "Uncategorized"))
            lines.append("\t".join([chrom, start, end, label, summary, compare, bg_label, "9", "0.9", "0.5", "80", "82", "2", "0.1",
                                    "81", "-30" if compare.startswith("Hypo") else "25", "-3.5" if compare.startswith("Hypo") else "3.2",
                                    "20", "0", "0", "30", "15", "15"]))
    return "\n".join(lines) + "\n"


def run_script(name: str, args: list[str], cwd: Path) -> str:
    res = subprocess.run([sys.executable, "-c", script(name), *args], capture_output=True, text=True, cwd=cwd)
    if res.returncode != 0:
        raise AssertionError(res.stderr)
    return res.stderr


def read_tsv(path: Path) -> list[dict[str, str]]:
    lines = [l for l in path.read_text().splitlines() if not l.startswith("##")]
    header = lines[0].split("\t")
    return [dict(zip(header, l.split("\t"))) for l in lines[1:]]


class AnnotateTest(unittest.TestCase):
    def test_annotation_events_and_counts(self) -> None:
        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d)
            (tmp / "map.tsv").write_text(MAP)
            (tmp / "S1.scored.tsv").write_text(scored_tsv({("CpG:_1", "ALL"): ("Unmethylated", "HypoMethylated"),
                                                           ("CpG:_1", "XX"): ("Unmethylated", "HypoMethylated"),
                                                           ("CpG:_3", "XX"): ("AlleleSpecificMethylation", "HyperASM"),
                                                           ("CpG:_2", "ALL"): ("NoData", "InsufficientData")}))
            err = run_script("methbat_annotate.py", ["S1.scored.tsv", "map.tsv", "S1", "S1.GRCh38"], tmp)
            ann = read_tsv(tmp / "S1.GRCh38.methbat.scored.annotated.tsv")
            self.assertEqual(len(ann), 6)
            self.assertEqual(list(ann[0])[-3:], ["hgnc_id", "gene_symbol", "region_type"])
            self.assertEqual((ann[0]["gene_symbol"], ann[0]["region_type"]), ("GENEA;GENEB", "promoter"))
            self.assertTrue((tmp / "S1.GRCh38.methbat.scored.annotated.tsv").read_text().startswith("##methbat_version"), "MethBat's header lines stay")
            events = read_tsv(tmp / "S1.GRCh38.methbat.events.tsv")
            self.assertEqual([(e["region_label"], e["background_category"], e["compare_label"]) for e in events],
                             [("CpG:_1", "ALL", "HypoMethylated"), ("CpG:_1", "XX", "HypoMethylated"), ("CpG:_3", "XX", "HyperASM")])
            self.assertTrue(all(e["sample_id"] == "S1" for e in events))
            self.assertEqual(events[2]["gene_symbol"], "GENEC")
            counts = {(c["label"], c["metric"]): int(c["value"]) for c in read_tsv(tmp / "S1.GRCh38.methbat.counts.tsv")}
            self.assertEqual(counts[("ALL", "regions")], 3)
            self.assertEqual(counts[("ALL", "HypoMethylated")], 1)
            self.assertEqual(counts[("XX", "HyperASM")], 1)
            self.assertEqual(counts[("ALL", "InsufficientData")], 1)
            self.assertEqual(counts[("ALL", "summary_Unmethylated")], 1)
            self.assertEqual(counts[("ALL", "summary_NoData")], 1)
            self.assertNotIn(("XX", "summary_Unmethylated"), counts, "summary labels are counted once, on the ALL rows")
            self.assertIn("S1: ALL: 1 events; XX: 2 events", err)


class SummaryTest(unittest.TestCase):
    def test_cohort_tables(self) -> None:
        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d)
            (tmp / "map.tsv").write_text(MAP)
            ids = ["S1", "S2", "S3"]
            labels = ["XX", "XY", ""]
            # per-member events and counts through the annotate script, so the two scripts agree on the formats
            plans = {
                "S1": {("CpG:_1", "ALL"): ("Unmethylated", "HypoMethylated"), ("CpG:_3", "ALL"): ("AlleleSpecificMethylation", "HyperASM")},
                "S2": {("CpG:_1", "ALL"): ("Unmethylated", "HypoMethylated")},
                "S3": {},
            }
            for sid in ids:
                (tmp / f"{sid}.scored.tsv").write_text(scored_tsv(plans[sid]))
                run_script("methbat_annotate.py", [f"{sid}.scored.tsv", "map.tsv", sid, f"{sid}.GRCh38"], tmp)
            # the background: 3 regions x 3 labels
            bg = ["##methbat_version=1.1.0-test", "\t".join(BG_COLS)]
            for chrom, start, end, label in REGIONS:
                for cat, n in (("ALL", 3), ("XX", 1), ("XY", 1)):
                    bg.append("\t".join([chrom, start, end, label, cat, "0", str(n), "0", "0", str(n), "0", "0", "1.5", "0.5", "80.0", "4.0"]))
            (tmp / "bg.tsv").write_text("\n".join(bg) + "\n")
            # reports: S1 anomalous at one region, S2 all PASS, S3 made in the stage
            for sid, summaries in (("S1", ["PASS", "Anomalous"]), ("S2", ["PASS", "PASS"]), ("S3", ["Inconclusive", "AnomalousQcWarning"])):
                rows = ["##methbat_version=1.1.0-test", "\t".join(REPORT_COLS)]
                for (region, summ) in zip(("PLAGL1:alt-TSS-DMR", "GNAS-XL:Ex1-DMR"), summaries):
                    rows.append("\t".join(["chr6", "1", "2", region, summ, "PASS" if summ != "AnomalousQcWarning" else "WeakASM", "AlleleSpecificMethylation",
                                           "AlleleSpecificMethylation", "50", "70", "15", "85", "0.0", "100", "0", "0", "30", "15", "15"]))
                (tmp / f"{sid}.report.tsv").write_text("\n".join(rows) + "\n")
            # one comparison with a note, and a deconvolution for S1
            cmp_rows = ["##methbat_version=1.1.0-test", "\t".join(COMPARE_COLS)]
            for chrom, start, end, _ in REGIONS:
                cmp_rows.append("\t".join([chrom, start, end, "XX", "XY", "InsufficientData", "", "0.0", "0", "0", "", "0.0", "1", "1"]))
            (tmp / "cmp.tsv").write_text("\n".join(cmp_rows) + "\n")
            (tmp / "S1.deconv.json").write_text(json.dumps({"metrics": {"normalized_fraction": {"B-cell": 0.6, "T-cell": 0.4}}}))
            w = lambda name, items: ((tmp / name).write_text("".join(x + "\n" for x in items)), str(tmp / name))[1]
            args = ["C1.methyl.GRCh38", w("ids", ids), w("labels", labels),
                    w("events", [str(tmp / f"{s}.GRCh38.methbat.events.tsv") for s in ids]),
                    w("counts", [str(tmp / f"{s}.GRCh38.methbat.counts.tsv") for s in ids]),
                    w("reports", [str(tmp / f"{s}.report.tsv") for s in ids]), w("sources", ["run", "run", "stage"]),
                    w("compares", [str(tmp / "cmp.tsv")]), w("notes", ["compare XX:XY: 1 vs 1 members, below min_samples 10: every row is InsufficientData"]),
                    str(tmp / "bg.tsv"), str(tmp / "map.tsv"), w("deconv", [str(tmp / "S1.deconv.json")]), "10"]
            err = run_script("methyl_summary.py", args, tmp)
            self.assertIn("3 members, 3 regions, 3 events, 6 imprinting rows", err)

            events = read_tsv(tmp / "C1.methyl.GRCh38.events.tsv")
            self.assertEqual([(e["sample_id"], e["region_label"]) for e in events], [("S1", "CpG:_1"), ("S2", "CpG:_1"), ("S1", "CpG:_3")],
                             "sorted by gene (GENEA before GENEC), then member")
            regions = {(r["region_label"], r["label"]): r for r in read_tsv(tmp / "C1.methyl.GRCh38.regions.tsv")}
            self.assertEqual(len(regions), 9)
            self.assertEqual((regions[("CpG:_1", "ALL")]["n_hypo"], regions[("CpG:_1", "ALL")]["n_datasets"], regions[("CpG:_1", "ALL")]["gene_symbol"]), ("2", "3", "GENEA;GENEB"))
            self.assertEqual((regions[("CpG:_3", "ALL")]["n_hyper_asm"], regions[("CpG:_3", "XX")]["n_hyper_asm"]), ("1", "0"))
            self.assertEqual(regions[("CpG:_2", "XY")]["region_type"], "intergenic")
            imp = read_tsv(tmp / "C1.methyl.GRCh38.imprinting.tsv")
            self.assertEqual(len(imp), 6)
            self.assertEqual({r["sample_id"]: r["report_source"] for r in imp}, {"S1": "run", "S2": "run", "S3": "stage"})
            samples = {r["sample_id"]: r for r in read_tsv(tmp / "C1.methyl.GRCh38.samples.tsv")}
            self.assertEqual((samples["S1"]["n_hypo"], samples["S1"]["n_hyper_asm"], samples["S1"]["n_events"]), ("1", "1", "2"))
            self.assertEqual((samples["S1"]["imprinting_pass"], samples["S1"]["imprinting_anomalous"]), ("1", "1"))
            self.assertEqual((samples["S3"]["imprinting_anomalous_qc"], samples["S3"]["report_source"]), ("1", "stage"))
            self.assertEqual((samples["S1"]["deconvolution_top_cell_type"], samples["S1"]["deconvolution_top_fraction"]), ("B-cell", "0.600"))
            self.assertEqual(samples["S2"]["deconvolution_top_cell_type"], "")
            self.assertEqual(samples["S2"]["labels"], "XY")
            summary = {(r["resource"], r["metric"], r["key"]): r["value"] for r in read_tsv(tmp / "C1.methyl.GRCh38.summary.tsv")}
            self.assertEqual(summary[("cohort", "n_members", "")], "3")
            self.assertEqual(summary[("cohort", "n_events", "")], "3")
            self.assertEqual(summary[("cohort", "reports_made_in_stage", "")], "1")
            self.assertEqual((summary[("labels", "n_members", "ALL")], summary[("labels", "n_members", "XX")]), ("3", "1"))
            self.assertIn(("notes", "label_below_min", "XX"), summary)
            self.assertEqual(summary[("genes", "members_with_event", "GENEA")], "2")
            self.assertEqual(summary[("genes", "members_with_event", "GENEC")], "1")
            self.assertEqual(summary[("region_type", "n_events", "promoter")], "2")
            self.assertEqual(summary[("imprinting_regions", "n_anomalous", "GNAS-XL:Ex1-DMR")], "2")
            self.assertEqual(summary[("compare", "XX:XY", "InsufficientData")], "3")
            self.assertIn(("notes", "compare", "XX:XY"), summary)
            self.assertEqual(summary[("deconvolution", "top_cell_type", "S1")], "B-cell")
            self.assertEqual(summary[("members", "n_hypo", "S2")], "1")
            self.assertEqual(summary[("imprinting", "Anomalous", "S1")], "1")


if __name__ == "__main__":
    unittest.main()
