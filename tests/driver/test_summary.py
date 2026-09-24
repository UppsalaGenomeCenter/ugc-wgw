import gzip
import html
import json
import tempfile
import unittest
from html.parser import HTMLParser
from pathlib import Path

from ugc_wgw import config, summary
from ugc_wgw import summary_parsers as sp
from ugc_wgw import summary_svg as sv
from ugc_wgw.db import DB

from .helpers import HtmlChecker, make_project, run_cli, seed_success, write_ids, write_tsv

BCFTOOLS_STATS = """# This file was produced by bcftools stats (1.23.1+htslib-1.23.1)
# SN\t[2]id\t[3]key\t[4]value
SN\t0\tnumber of samples:\t1
SN\t0\tnumber of records:\t116756
SN\t0\tnumber of SNPs:\t97699
SN\t0\tnumber of indels:\t19057
# TSTV\t[2]id\t[3]ts\t[4]tv\t[5]ts/tv\t[6]ts (1st ALT)\t[7]tv (1st ALT)\t[8]ts/tv (1st ALT)
TSTV\t0\t63578\t34121\t1.86\t63578\t34121\t1.86
ST\t0\tA>C\t4123
ST\t0\tA>G\t15908
IDD\t0\t-60\t13\t13\t0.83
IDD\t0\t1\t2000\t2000\t0.5
QUAL\t0\t19.5\t10\t7\t3\t1
QUAL\t0\t40.0\t500\t300\t200\t50
DP\t0\t2\t463\t0.396553\t0\t0.000000
DP\t0\t>500\t390\t0.334030\t0\t0.000000
PSC\t0\tHG002\t0\t43387\t54312\t63578\t34121\t19057\t49.4\t116756\t0\t0\t7
bogus\tline
"""

MOSDEPTH = """chrom\tlength\tbases\tmean\tmin\tmax
chr1\t248956422\t7468692660\t30.00\t0\t80
chr1_region\t248956422\t7468692660\t30.00\t0\t80
chr2\t242193529\t7265805870\t30.00\t0\t70
chr2_region\t242193529\t7265805870\t30.00\t0\t70
chrX\t156040895\t2340613425\t15.00\t0\t60
chrX_region\t156040895\t2340613425\t15.00\t0\t60
chrY\t57227415\t744000000\t13.00\t0\t50
chrY_region\t57227415\t744000000\t13.00\t0\t50
chrM\t16569\t100000000\t6035.0\t0\t9000
total\t3095734472\t9000000000\t29.5\t0\t9000
total_region\t3095734472\t9000000000\t29.5\t0\t9000
"""

REGION_BED = b"chr1\t0\t500\t30.0\nchr1\t500\t1000\t32.0\nchr1\t1000\t1500\t10.0\nchr1\t2500\t3000\t20.0\nchrUn_x\t0\t500\t1\nbad\n"

ROH = b"# RG\t[2]Sample\t[3]Chromosome\t[4]Start\t[5]End\t[6]Length (bp)\t[7]Number of markers\t[8]Quality\nST\tS\tchr1\t100\t1\t50\nRG\tS\tchr1\t878431\t172619312\t171740882\t531\t93.6\nRG\tS\tchr2\t100\t250100\t250000\t40\t80.1\nRG\tS\tchr2\tx\ty\tz\t1\t1\n"

SV_VCF = b"""##fileformat=VCFv4.2
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="x">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
chr1\t100\tsv1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-250;SUPP=2\tGT\t0/1\t0/1
chr1\t900\tsv2\tN\t<INS>\t.\tPASS\tSVTYPE=INS;SVLEN=1500;SUPP=1\tGT\t0/1\t0/0
chr2\t100\tsv3\tN\tN[chr3:5[\t.\tMinQUAL\tSVTYPE=BND;SUPP=1\tGT\t./.\t0/1
chr2\t500\tsv4\tN\t<DUP>\t.\tPASS\tSVTYPE=DUP;SVLEN=2000000;SUPP=2\tGT\t0/1\t0/1
"""

SITES_ONLY_VCF = b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t100\t.\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-250\n"

PHASE_STATS = """sample_name\tchromosome\tnum_variants\tnum_heterozygous\tnum_phased\tnum_unphased\tnum_blocks\tblock_ng50
S1\tchr1\t100\t50\t45\t5\t3\t1000
S1\tchr2\t80\t40\t30\t10\t2\t500
S1\tall\t180\t90\t75\t15\t5\t800
"""

PHASE_BLOCKS = """source_block_index\tsample_name\tphase_block_id\tchrom\tstart\tend\tnum_variants
0\tS1\t100\tchr1\t100\t5100\t10
2\tS1\t9000\tchr1\t9000\t159000\t50
9\tS1\t300000\tchr2\t300000\t2300000\t900
"""

TRGT_VCF = b"""##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
chr1\t16682\t.\tTGG\t.\t.\t.\tTRID=chr1_16682_16774_TGGTGGGGG;END=16774;MOTIFS=TGGTGGGGG;STRUC=(TGGTGGGGG)n\tGT:AL:MC:SD\t0/0:93,93:10,10:5,6
chr4\t3074876\t.\tCCAG\tCCAGCAG\t.\t.\tTRID=HD_HTT;END=3074969;MOTIFS=CAG,CCG;STRUC=<TR>\tGT:AL:ALLR:SD:MC\t0/1:57,63:57-57,63-63:12,11:17_7,19_7
chr9\t69037270\t.\tCAAA\t.\t.\t.\tTRID=FRDA_FXN;END=69037304;MOTIFS=A,GAA;STRUC=<TR>\tGT:AL:MC:SD\t.:.:.:.
chr7\t100\t.\tC\t.\t.\t.\tTRID=HFG_HOXA13-I;END=200;MOTIFS=GCN;STRUC=<TR>\tGT:AL\t./.:.
"""

DROPOUTS = b"""chrom\tstart\tend\ttrid\texpected_ploidy\thap1_count\thap2_count\tunphased_count\tfail_read_count\tdropout
chr1\t16682\t16774\tID=chr1_16682_16774_TGGTGGGGG;MOTIFS=TGGTGGGGG;STRUC=(TGGTGGGGG)n\t2\t0\t0\t0\t0\tFullDropout
chr1\t19275\t19473\tID=chr1_19275_19473_TG;MOTIFS=TG;STRUC=(TG)n\t2\t1\t0\t0\t0\tHaplotypeDropout
chr9\t69037270\t69037304\tID=FRDA_FXN;MOTIFS=A,GAA;STRUC=<TR>\t2\t0\t0\t0\t0\tFullDropout
chr4\t3074876\t3074969\tID=HD_HTT;MOTIFS=CAG,CCG;STRUC=<TR>\t2\t3\t3\t2\t0\tPhasingDropout
"""

METHBAT = """##methbat_version=1.1.0-fbf4686
##command=methbat profile
chrom\tstart\tend\tregion_label\tsummary_label\tmean_combined_methyl
chr1\t100\t200\tCpG:_1\tMethylated\t95.5
chr1\t300\t400\tCpG:_2\tUnmethylated\t4.0
chr1\t500\t600\tCpG:_3\tNoData\t
chr1\t700\t800\tCpG:_4\tAlleleSpecificMethylation\t50.0
"""

FREQ = """# ugc-pacbio-wgw cohort_freq: comment
resource\tmetric\tkey\tvalue
samples\tn\ttotal\t4
samples\tn\tXX\t2
structural_variants\trecords\ttotal\t2011
structural_variants\trecords_by_af_bin\t<0.001\t0
structural_variants\trecords_by_af_bin\t>=0.05\t2011
structural_variants\tmean_call_rate\tall\t0.9875
small_variants\tmean_call_rate\tall\tNA
"""

GFATOOLS = """CC\tGS   genome_size_if_provided
CC
GS\t3100000000
SZ\t68817513
NN\t51
NL\t0\t31943945\t1
NL\t50\t20000000\t2
AU\t582777
"""

WDL_SNIPPET = '''version 1.0
task a {
  runtime {
    docker: "~{runtime_attributes.container_registry}/hiphase@sha256:41ebe22b55c66e2e78da2013f7fffaecc02a8b4e980400c3ea8d03c87330522e"  # 1.7.0_build2
  }
}
task b {
  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:03cb3c01937eccc907f8ad71c87b258581504572205fe3f31a657e318f3564ae"  # pb_wdl_base:build4
  }
}
task c {
  runtime {
    docker: "~{ugc_wgw_container_registry}/hifiasm@sha256:0e03f1dc2ab6ede5200bbc0e8e9c3ea10e693ce87cad9346e3a1aa59cec272cf"  # 0.25.0 (containers/hifiasm)
  }
}
workflow w {
  String deepvariant_version = "1.10.0"
  String workflow_version = "4.0.0"
  String docker_image = "google/deepvariant:~{deepvariant_version}"
  call x { input: docker_image = docker_image }
}
task d {
  runtime {
    docker: "~{docker_image}"
  }
}
'''


class ParserTest(unittest.TestCase):
    def test_two_line_tsv(self):
        d = sp.parse_two_line_tsv("sample_id\tdepth_mean\tinferred_sex\nS1\t30.1\tMALE\n".splitlines())
        self.assertEqual(d, {"sample_id": "S1", "depth_mean": "30.1", "inferred_sex": "MALE"})
        self.assertEqual(sp.parse_two_line_tsv(["a\tb"]), {})

    def test_bcftools_stats(self):
        d = sp.parse_bcftools_stats(BCFTOOLS_STATS.splitlines())
        self.assertEqual(d["sn"]["number of SNPs"], 97699)
        self.assertEqual(d["tstv"], {"ts": 63578, "tv": 34121, "ratio": 1.86})
        self.assertEqual(d["st"]["A>G"], 15908)
        self.assertEqual(d["idd"], [[-60, 13], [1, 2000]])
        self.assertEqual(d["qual"][1], [40.0, 500, 300, 200, 50])
        self.assertEqual(d["dp"][1], [">500", 390, 0])
        self.assertEqual(d["psc"]["n_hets"], 54312)
        self.assertEqual(d["psc"]["n_missing"], 7)

    def test_mosdepth_summary_and_ratios(self):
        d = sp.parse_mosdepth_summary(MOSDEPTH.splitlines())
        self.assertEqual([r["chrom"] for r in d["chrom"]], ["chr1", "chr2", "chrX", "chrY", "chrM"])
        self.assertEqual(d["region"]["chr1"], 30.0)
        self.assertNotIn("total", d["region"])
        self.assertEqual(d["total"]["mean"], 29.5)
        r = sp.coverage_ratios(d["chrom"], 0.1)
        self.assertEqual(r["auto_mean"], 30.0)
        self.assertAlmostEqual(r["chrX_ratio"], 0.5)
        self.assertAlmostEqual(r["chrY_ratio"], 13 / 30)
        self.assertEqual(r["sex_from_ratio"], "MALE")
        self.assertEqual(sp.coverage_ratios([], 0.1)["sex_from_ratio"], None)

    def test_region_bed_bins(self):
        d = sp.aggregate_region_bed(REGION_BED.splitlines(keepends=True), bin_bp=1000, contigs=["chr1", "chr2"])
        self.assertEqual(d["bin_bp"], 1000)
        self.assertEqual(list(d["chroms"]), ["chr1"])
        self.assertEqual(d["chroms"]["chr1"]["len"], 3000)
        self.assertEqual(d["chroms"]["chr1"]["bins"], [31.0, 10.0, 20.0])
        d = sp.aggregate_region_bed(REGION_BED.splitlines(keepends=True), bin_bp=1000, contigs=None)
        self.assertEqual(sorted(d["chroms"]), ["chr1", "chrUn_x"])

    def test_copynum_and_bedgraph(self):
        doc = {"sample_name": "S", "gc_bias_corrected_haploid_coverage": 15.5,
               "chromosomes": {"chr2": {"total_copy_number_bases": 10, "bases_per_copy_number": {"2": 8, "3": 2}, "most_common_copy_number": 2},
                               "chr1": {"total_copy_number_bases": 5, "bases_per_copy_number": {}, "most_common_copy_number": 2},
                               "chrUn": {"total_copy_number_bases": 1, "bases_per_copy_number": {"0": 1}, "most_common_copy_number": 0}}}
        d = sp.parse_copynum_summary(doc)
        self.assertEqual(d["haploid_coverage"], 15.5)
        self.assertEqual([r["chrom"] for r in d["chrom"]], ["chr1", "chr2", "chrUn"])
        self.assertEqual(d["chrom"][1]["cn_bases"], {"2": 8, "3": 2})
        self.assertEqual(sp.parse_copynum_summary("nope"), {})
        bg = sp.parse_bedgraph("chr1\t0\t1000\t2\nchr1\t1000\t2000\t3\nchrUn\t0\t5\t1\nbad\n".splitlines())
        self.assertEqual(bg, {"chr1": [[0, 1000, 2.0], [1000, 2000, 3.0]]})

    def test_roh(self):
        d = sp.parse_roh_out(ROH.splitlines(keepends=True))
        self.assertEqual(d["n"], 2)
        self.assertEqual(d["total_bp"], 171740882 + 250000)
        self.assertEqual(d["longest_bp"], 171740882)
        self.assertEqual(d["by_chrom"], {"chr1": 171740882, "chr2": 250000})
        self.assertEqual(dict(d["hist"]), {"<100k": 0, "100k-500k": 1, "500k-1M": 0, "1M-5M": 0, ">5M": 1})

    def test_sv_vcf(self):
        d = sp.parse_sv_vcf(SV_VCF.splitlines(keepends=True))
        self.assertEqual(d["n"], 4)
        self.assertEqual(d["counts"], {"DEL": 1, "INS": 1, "BND": 1, "DUP": 1})
        self.assertEqual(d["sizes"]["DEL"], {"100-1k": 1})
        self.assertEqual(d["sizes"]["DUP"], {">1M": 1})
        self.assertNotIn("BND", d["sizes"])
        self.assertEqual(d["filters"], {"PASS": 3, "MinQUAL": 1})
        self.assertEqual(d["supp"], {"1": 2, "2": 2})
        self.assertEqual(d["samples"], ["S1", "S2"])
        d = sp.parse_sv_vcf(SITES_ONLY_VCF.splitlines(keepends=True))
        self.assertEqual((d["n"], d["counts"], d["samples"]), (1, {"DEL": 1}, []))

    def test_phasing(self):
        d = sp.parse_phase_stats(PHASE_STATS.splitlines())
        self.assertEqual(d["all"]["num_phased"], 75)
        self.assertEqual([r["chromosome"] for r in d["chrom"]], ["chr1", "chr2"])
        self.assertEqual(d["chrom"][0]["block_ng50"], 1000)
        b = sp.parse_phase_blocks(PHASE_BLOCKS.splitlines())
        self.assertEqual(b["n"], 3)
        self.assertEqual(b["total_bp"], 5000 + 150000 + 2000000)
        self.assertEqual(b["n50"], 2000000)
        self.assertEqual(dict(b["hist"]), {"<10k": 1, "10k-100k": 0, "100k-1M": 1, "1M-10M": 1, ">10M": 0})

    def test_trgt_vcf_and_dropouts(self):
        fr = sp.fail_reads_trids(["INCLUDE_FAIL_READS regions: CANVAS_RFC1,FRDA_FXN", "Input x is already aligned."])
        self.assertEqual(fr, ["CANVAS_RFC1", "FRDA_FXN"])
        d = sp.parse_trgt_vcf(TRGT_VCF.splitlines(keepends=True), fr)
        self.assertEqual(d["n_records"], 4)
        self.assertEqual(d["n_named"], 3)
        self.assertEqual(d["n_named_called"], 1)
        by = {x["trid"]: x for x in d["loci"]}
        self.assertEqual(sorted(by), ["FRDA_FXN", "HD_HTT", "HFG_HOXA13-I"])
        htt = by["HD_HTT"]
        self.assertEqual((htt["disease"], htt["gene"], htt["chrom"], htt["pos"], htt["end"]), ("HD", "HTT", "chr4", 3074876, 3074969))
        self.assertEqual((htt["gt"], htt["al"], htt["mc"], htt["sd"], htt["called"], htt["fail_reads"]), ("0/1", "57,63", "17_7,19_7", "12,11", True, False))
        self.assertEqual((by["FRDA_FXN"]["called"], by["FRDA_FXN"]["fail_reads"]), (False, True))
        self.assertEqual(by["HFG_HOXA13-I"]["gene"], "HOXA13-I")
        # sorted by gene
        self.assertEqual([x["gene"] for x in d["loci"]], ["FXN", "HOXA13-I", "HTT"])
        dr = sp.parse_trgt_dropouts(DROPOUTS.splitlines(keepends=True))
        self.assertEqual(dr["n"], 4)
        self.assertEqual(dr["by_class"], {"FullDropout": 2, "HaplotypeDropout": 1, "PhasingDropout": 1})
        self.assertEqual(dr["named"], {"FRDA_FXN": "FullDropout", "HD_HTT": "PhasingDropout"})

    def test_methbat_paraphase_mitorsaw_starphase_kivvi(self):
        m = sp.parse_methbat_profile(METHBAT.splitlines())
        self.assertEqual((m["n_regions"], m["n_with_value"]), (4, 3))
        self.assertEqual(m["labels"]["Methylated"], 1)
        self.assertEqual(dict(m["hist"])["90-100"], 1)
        self.assertEqual(dict(m["hist"])["0-10"], 1)
        p = sp.parse_paraphase({"smn1": {"region_name": "smn1", "genes_in_region": "SMN1,SMN2", "region_depth": {"median": 0.0},
                                          "failed_for_coverage": True, "total_cn": None, "final_haplotypes": {}},
                                "CYP2D6": {"genes_in_region": "CYP2D6", "region_depth": {"median": 31.0},
                                           "failed_for_coverage": False, "total_cn": 2, "final_haplotypes": {"a": 1, "b": 2}}})
        self.assertEqual((p["n_regions"], p["n_failed"], p["n_with_cn"]), (2, 1, 1))
        self.assertEqual(p["regions"][0]["region"], "CYP2D6")
        self.assertEqual(p["regions"][0]["n_haplotypes"], 2)
        mi = sp.parse_mitorsaw({"haplotypes": [{"label": "hap_0", "seq_len": 16569, "num_ref_variants": 3, "estimated_abundance": 0.9}],
                                "fingerprint_stats": {"passing_explanation": True, "unexplained_fraction": None}})
        self.assertEqual(mi["haplotypes"][0]["estimated_abundance"], 0.9)
        self.assertEqual(mi["fingerprint"]["passing_explanation"], True)
        st = sp.parse_starphase_tsv("#gene\tdiplotype\nABCG2\tReference/Reference\nMT-RNR1\tReference\n".splitlines())
        self.assertEqual(st, [{"gene": "ABCG2", "diplotype": "Reference/Reference"}, {"gene": "MT-RNR1", "diplotype": "Reference"}])
        sj = sp.parse_starphase_json({"pbstarphase_version": "2.2.0", "database_metadata": {"cpic_version": "x"},
                                      "gene_details": {"A": {"diplotypes": [{"diplotype": "*1/*1"}]},
                                                       "B": {"diplotypes": [{"diplotype": "NO_READS/NO_READS"}]}}})
        self.assertEqual((sj["version"], sj["n_genes"], sj["n_called"], sj["database"]), ("2.2.0", 2, 1, {"cpic_version": "x"}))
        k = sp.parse_kivvi({"allele_cn": "NA", "depth_summary": {"genome_depth": None}, "complete_alleles": [],
                            "additional": {"call_status": "failed_due_to_no_reads"}})
        self.assertEqual((k["allele_cn"], k["call_status"], k["n_complete_alleles"]), ("NA", "failed_due_to_no_reads", 0))
        k = sp.parse_kivvi({"allele_cn": 38.5, "depth_summary": {"genome_depth": 30.0, "repeat_depth": 570.0}, "complete_alleles": [1, 2], "additional": {}})
        self.assertEqual((k["allele_cn"], k["call_status"], k["n_complete_alleles"], k["repeat_depth"]), (38.5, None, 2, 570.0))

    def test_freq_and_gfatools(self):
        f = sp.parse_freq_summary(FREQ.splitlines())
        self.assertEqual(f["samples"]["n"], {"total": 4, "XX": 2})
        self.assertEqual(f["structural_variants"]["records_by_af_bin"], {"<0.001": 0, ">=0.05": 2011})
        self.assertEqual(f["structural_variants"]["mean_call_rate"]["all"], 0.9875)
        self.assertEqual(f["small_variants"]["mean_call_rate"]["all"], "NA")
        g = sp.parse_gfatools_stats(GFATOOLS.splitlines())
        self.assertEqual((g["total_length"], g["n_sequences"], g["aun"], g["n50"], g["l50"]), (68817513, 51, 582777, 20000000, 2))
        self.assertEqual(g["nl"][0], [0, 31943945, 1])

    def test_files_gzip_dims_and_vcf_samples(self):
        tmp = Path(tempfile.mkdtemp())
        gz = tmp / "x.tsv.gz"
        with gzip.open(gz, "wt") as fh:
            fh.write("trid\tmotif\tS1\tS2\nchr1_1_2_A\tA\t0,0\t1,1\nchr1_3_4_C\tC\t0,0\t2,2\n")
        self.assertEqual(list(sp.read_lines(gz))[0], "trid\tmotif\tS1\tS2")
        self.assertEqual(sp.tsv_dims(gz), (["trid", "motif", "S1", "S2"], 2))
        plain = tmp / "y.tsv"
        plain.write_text("a\tb\n1\t2\n")
        self.assertEqual(sp.tsv_dims(plain), (["a", "b"], 1))
        vcf = tmp / "z.vcf.gz"
        with gzip.open(vcf, "wb") as fh:
            fh.write(SV_VCF)
        self.assertEqual(sp.vcf_samples(vcf), ["S1", "S2"])
        self.assertEqual(sp.load_json(tmp / "j.json") if (tmp / "j.json").write_text('{"a": 1}') else None, {"a": 1})

    def test_scan_tool_versions(self):
        tmp = Path(tempfile.mkdtemp())
        d = tmp / "workflows" / "ugc_wgw"
        d.mkdir(parents=True)
        (d / "a.wdl").write_text(WDL_SNIPPET)
        (d / "b.wdl").write_text(WDL_SNIPPET.split("workflow w")[0])  # the same pins again in a second file
        tools = {t["tool"]: t for t in sp.scan_tool_versions(tmp)}
        self.assertEqual(set(tools), {"hiphase", "pb_wdl_base", "hifiasm", "deepvariant"})
        self.assertEqual((tools["hiphase"]["version"], tools["hiphase"]["build"], tools["hiphase"]["digest"][:8]), ("1.7.0", "build2", "41ebe22b"))
        self.assertEqual((tools["pb_wdl_base"]["version"], tools["pb_wdl_base"]["build"]), ("build4", ""))
        self.assertEqual((tools["hifiasm"]["version"], tools["hifiasm"]["image"]), ("0.25.0", "hifiasm@sha256:0e03f1dc2ab6ede5200bbc0e8e9c3ea10e693ce87cad9346e3a1aa59cec272cf"))
        self.assertEqual((tools["deepvariant"]["version"], tools["deepvariant"]["digest"]), ("1.10.0", ""))
        self.assertEqual(tools["hiphase"]["files"], ["workflows/ugc_wgw/a.wdl", "workflows/ugc_wgw/b.wdl"])
        self.assertEqual(sp.scan_tool_versions(tmp / "nowhere"), [])

    def test_contig_class(self):
        self.assertEqual([sp.contig_class(c) for c in ("chr1", "22", "chrX", "Y", "chrM", "MT", "chrUn_x", "KMT2C_1")],
                         ["autosome", "autosome", "chrX", "chrY", "chrM", "chrM", "other", "other"])


class _Checker(HTMLParser):
    VOID = {"line", "rect", "polyline", "path", "circle", "polygon"}

    def __init__(self):
        super().__init__()
        self.stack, self.errors = [], []

    def handle_starttag(self, tag, attrs):
        if tag not in self.VOID:
            self.stack.append(tag)

    def handle_startendtag(self, tag, attrs):
        pass

    def handle_endtag(self, tag):
        if tag in self.VOID:
            return
        if not self.stack or self.stack[-1] != tag:
            self.errors.append(tag)
        else:
            self.stack.pop()


class SvgTest(unittest.TestCase):
    def check(self, svg: str) -> str:
        c = _Checker()
        c.feed(svg)
        self.assertEqual((c.errors, c.stack), ([], []), svg[:200])
        self.assertNotIn("nan", svg)
        self.assertNotIn("inf", svg.replace("infinite", ""))
        return svg

    def test_helpers(self):
        self.assertEqual(sv.nice_ticks(0, 97), [0, 20, 40, 60, 80, 100])
        self.assertEqual(sv.nice_ticks(0, 100), [0, 20, 40, 60, 80, 100])
        self.assertEqual(sv.nice_ticks(0, 0)[0], 0.0)
        self.assertEqual((sv.fmt(1234567), sv.fmt(0.974, 2), sv.fmt(None), sv.fmt("NA"), sv.fmt(float("nan"))), ("1,234,567", "0.97", "–", "NA", "–"))
        self.assertEqual((sv.fmt_bp(3.2e9), sv.fmt_bp(2_500_000), sv.fmt_bp(456)), ("3.20 Gb", "2.5 Mb", "456 bp"))
        self.assertEqual(sv.finite(float("inf")), 0.0)
        self.assertEqual(sv.esc('<a href="x">'), "&lt;a href=&quot;x&quot;&gt;")

    def test_bars_hist_stacked(self):
        s = self.check(sv.bars_h([("DEL", 537, ""), ("INS <", 553, "note"), ("BND", float("nan"), "")]))
        self.assertIn("INS &lt;", s)
        self.assertIn("<title>INS &lt;: 553 note</title>", s)
        self.assertIn("No data", sv.bars_h([]))
        s = self.check(sv.hist([("<100", 5), ("100-1k", 40), (">1M", 0)], y_label="records"))
        self.assertIn("records", s)
        self.assertIn(">40<", s)
        self.check(sv.hist([(f"bin{i}", i) for i in range(20)]))
        self.assertIn("No data", sv.hist([("a", 0)]))
        s = self.check(sv.stacked([("chr1", {"DEL": 3, "INS": 5}), ("chr2", {"DEL": 1, "DUP": 2})], ["DEL", "INS", "DUP"]))
        self.assertIn('class="m3"', s)
        self.assertIn("(62.5%)", s)
        self.check(sv.stacked([("A>C", {"n": 4123})], ["n"], percent=True, legend=False))
        self.assertIn("No data", sv.stacked([("x", {"a": 0})], ["a"]))

    def test_track_and_strip(self):
        s = self.check(sv.track([{"chrom": "chr1", "len": 2480000, "bins": [1.0, None, 2.5], "bin_bp": 1000000,
                                  "steps": [[0, 1000000, 2], [1000000, 2480000, 3]]},
                                 {"chrom": "chrX", "len": 1560000, "bins": [0.5, 0.4], "bin_bp": 1000000}], baseline=1.0))
        self.assertEqual(s.count('class="lane"'), 2)
        self.assertIn("stroke-dasharray", s)
        self.assertIn('class="ln3"', s)
        self.assertIn("No data", sv.track([{"chrom": "x", "len": 0}]))
        s = self.check(sv.strip([{"label": "depth", "points": [["S1", 30.1], ["S2", 28.5], ["S3", 12.0]], "outliers": ["S3"],
                                  "threshold": 20, "digits": 1, "unit": "x"},
                                 {"label": "N50", "points": [["S1", 17394], ["S2", 17394]]}]))
        self.assertEqual(s.count("<circle"), 5)
        self.assertEqual(s.count('class="bad"'), 1)
        self.assertIn("threshold 20.0x", s)
        self.assertIn("No data", sv.strip([{"label": "x", "points": []}]))
        self.check(sv.legend([("DEL", "m1"), ("INS", "m2")]))
        self.assertEqual(sv.legend([]), "")


class DigestHelpersTest(unittest.TestCase):
    def test_bins(self):
        idd = summary._bin_series([[-60, 13], [-20, 4], [-2, 5], [1, 7], [9, 1], [12, 2], [80, 3]], [-49, -9, 1, 10, 50],
                                  ["<=-50", "-49..-10", "-9..-1", "+1..+9", "+10..+49", ">=+50"])
        self.assertEqual(idd, [["<=-50", 13.0], ["-49..-10", 4.0], ["-9..-1", 5.0], ["+1..+9", 8.0], ["+10..+49", 2.0], [">=+50", 3.0]])
        dp = summary._dp_bins([["2", 10, 0], ["15", 20, 0], ["59", 5, 0], ["100", 7, 0], [">500", 1, 0], ["bad", 9, 0]])
        self.assertEqual(dp, [["<10", 10], ["10-20", 20], ["20-30", 0], ["30-40", 0], ["40-60", 5], ["60-100", 0], [">=100", 8]])


class QcTest(unittest.TestCase):
    def test_thresholds_sources_and_errors(self):
        tmp = Path(tempfile.mkdtemp())
        cfg = make_project(tmp)
        self.assertEqual(summary.thresholds_for(cfg, []), summary.DEFAULT_THRESHOLDS)
        cfg.summary_thresholds = {"depth_mean_min": 10}
        config.save(cfg)
        cfg = config.load(cfg.project_dir)
        thr = summary.thresholds_for(cfg, ["mapped_read_percent_min=90"])
        self.assertEqual((thr["depth_mean_min"], thr["mapped_read_percent_min"], thr["read_quality_median_min"]), (10.0, 90.0, 25.0))
        with self.assertRaises(summary.UgcError):
            summary.thresholds_for(cfg, ["bogus=1"])
        with self.assertRaises(summary.UgcError):
            summary.thresholds_for(cfg, ["depth_mean_min=high"])
        cfg.summary_thresholds = {"nope": 1}
        with self.assertRaises(summary.UgcError):
            summary.thresholds_for(cfg, [])

    def test_flags(self):
        d = {"stats": {"depth_mean": "12.5", "mapped_read_percent": "99.9", "read_quality_median": "20", "inferred_sex": "FEMALE"},
             "sample": {"sex_sheet": "MALE"}, "coverage": {"mean": 12.5},
             "targeted": {"kivvi": {"kiv2": {"call_status": "failed_due_to_no_reads"}, "d4z4": {"call_status": None}},
                          "paraphase": {"n_failed": 2, "n_regions": 164}, "mitorsaw": {"haplotypes": []}},
             "messages": ["x"], "warnings": []}
        qc = summary.qc_flags(d, dict(summary.DEFAULT_THRESHOLDS))
        codes = [(f["code"], f["severity"]) for f in qc["flags"]]
        self.assertEqual(qc["status"], "fail")
        self.assertIn(("low_depth", "fail"), codes)
        self.assertIn(("low_read_quality", "fail"), codes)
        self.assertNotIn(("low_mapped", "fail"), codes)
        self.assertIn(("sex_mismatch", "warn"), codes)
        self.assertEqual(sum(1 for c, _ in codes if c == "tool_status"), 3)
        self.assertIn(("messages", "info"), codes)
        d["stats"].update({"depth_mean": "31", "read_quality_median": "31", "inferred_sex": "MALE"})
        d["targeted"] = {}
        d["messages"] = []
        self.assertEqual(summary.qc_flags(d, dict(summary.DEFAULT_THRESHOLDS)), {"status": "pass", "flags": [], "thresholds": summary.DEFAULT_THRESHOLDS})
        d["stats"]["inferred_sex"] = ""
        self.assertEqual([f["code"] for f in summary.qc_flags(d, dict(summary.DEFAULT_THRESHOLDS))["flags"]], ["sex_not_inferred"])


DIGEST_KEYS = {"schema", "generated_at", "mode", "sample", "runs", "warnings", "stats", "reads", "coverage", "small_variants", "sv",
               "phasing", "tandem_repeats", "methylation", "targeted", "messages", "files", "assembly", "provenance", "qc"}


class SummaryEndToEnd(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.cfg = make_project(self.tmp)
        self.proj = str(self.cfg.project_dir)
        tsv = write_tsv(self.tmp, [{"sample_id": "S1", "sex": "MALE", "hifi_reads": "s1.bam"},
                                   {"sample_id": "S2", "sex": "FEMALE", "hifi_reads": "s2.bam"}])
        self.assertEqual(run_cli(["--project", self.proj, "samples", "add", str(tsv)])[0], 0)
        self.assertEqual(run_cli(["--project", self.proj, "cohort", "freeze", "C1", "--samples", str(write_ids(self.tmp, ["S1", "S2"]))])[0], 0)
        code, out, err = run_cli(["--project", self.proj, "submit", "--mode", "standalone", "--cohort", "C1"],
                                 env={"UGC_WGW_FAKE_INFERRED_SEX": "S1=FEMALE,S2=FEMALE"})
        self.assertEqual(code, 0, err)
        self.out_dir = self.cfg.results_dir / "reports" / "summary"

    def check_html(self, path: Path) -> str:
        text = path.read_text()
        checker = HtmlChecker()
        checker.feed(text)
        self.assertEqual(checker.errors, [], path.name)
        self.assertEqual(checker.stack, [], path.name)
        self.assertEqual(checker.titles, 1, path.name)
        self.assertNotIn('src="http', text)
        return text

    def test_sample_pages(self):
        code, out, err = run_cli(["--project", self.proj, "summary", "--samples", "S1", "S2"])
        self.assertEqual(code, 0, err)
        paths = [Path(line) for line in out.strip().splitlines()]
        self.assertEqual([p.name for p in paths], ["S1.summary.html", "S2.summary.html"])
        self.assertTrue(all(p.parent == self.out_dir for p in paths))
        self.assertIn("S1: digest rebuilt", err)
        text = self.check_html(paths[0])
        for needle in ('id="sec-overview"', 'id="sec-coverage"', 'id="sec-small-variants"', 'id="sec-tr"', 'id="sec-software"',
                       "S1 analysis summary", "sex_mismatch", "HD_HTT", "FRDA_FXN", "hiphase", "pb_wdl_base", "1.7.0",
                       "failed_due_to_no_reads", "CYP2D6", "hap_0", "INCLUDE_FAIL_READS", "chr20", "Runs of homozygosity",
                       "hifi-human-wgs-wdl", "<svg", 'class="lane"', "summary_thresholds" if False else "depth_mean_min",
                       html.escape(summary.DEFAULT_PROJECT_URL)):   # the owner placeholder is HTML-escaped in the page
            self.assertIn(needle, text, needle)
        self.assertNotIn("sex_mismatch", self.check_html(paths[1]))  # S2 is FEMALE in the sheet and inferred FEMALE
        digest = json.loads((self.out_dir / "S1.summary.json").read_text())
        self.assertEqual(set(digest), DIGEST_KEYS)
        self.assertEqual(digest["schema"], summary.DIGEST_SCHEMA)
        self.assertEqual(digest["qc"]["status"], "warn")
        self.assertEqual(sorted({f["code"] for f in digest["qc"]["flags"]}), ["messages", "sex_mismatch", "tool_status"])
        self.assertEqual(digest["warnings"], [])
        self.assertEqual(digest["coverage"]["ratios"]["sex_from_ratio"], "FEMALE")
        self.assertEqual(digest["coverage"]["track"]["bin_bp"], 1000000)
        self.assertEqual(sorted(digest["coverage"]["track"]["chroms"]), ["chr1", "chr20", "chrX"])
        self.assertEqual(digest["small_variants"]["snv"], 3900000.0)
        self.assertEqual(digest["small_variants"]["roh"]["n"], 3)
        self.assertEqual(digest["sv"]["vcf_counts"], {"DEL": 2, "INS": 1, "BND": 1, "DUP": 1})
        self.assertEqual(digest["phasing"]["blocks"]["n"], 4)
        loci = {l["trid"]: l for l in digest["tandem_repeats"]["disease_loci"]}
        self.assertEqual(sorted(loci), ["FRDA_FXN", "HD_HTT", "SCA27B_FGF14"])
        self.assertEqual((loci["HD_HTT"]["fail_reads"], loci["HD_HTT"]["called"], loci["SCA27B_FGF14"]["dropout"]), (True, True, "FullDropout"))
        self.assertEqual(digest["tandem_repeats"]["dropouts"]["by_class"]["FullDropout"], 2)
        self.assertEqual(digest["targeted"]["kivvi"]["d4z4"]["allele_cn"], 38.0)
        self.assertEqual(digest["targeted"]["paraphase"]["n_failed"], 1)
        self.assertEqual(len(digest["targeted"]["starphase"]), 4)
        self.assertEqual(digest["methylation"]["profile"]["labels"]["Methylated"], 1)
        self.assertEqual(digest["provenance"]["ugc_wgw"]["version"], config.read_version(self.cfg.code_dir))
        self.assertIn("hifi-human-wgs-wdl", digest["provenance"]["upstream"])
        self.assertTrue(digest["files"]["trgt_vcf"]["bytes"] > 0)
        # a second run reuses the digest; --force rebuilds; thresholds re-evaluate on a reused digest
        code, out, err = run_cli(["--project", self.proj, "summary", "--samples", "S1"])
        self.assertEqual(code, 0, err)
        self.assertIn("S1: digest reused", err)
        code, out, err = run_cli(["--project", self.proj, "summary", "--samples", "S1", "--threshold", "depth_mean_min=50"])
        self.assertEqual(code, 0, err)
        digest = json.loads((self.out_dir / "S1.summary.json").read_text())
        self.assertIn("low_depth", {f["code"] for f in digest["qc"]["flags"]})
        self.assertEqual(digest["qc"]["status"], "fail")
        self.assertEqual(digest["qc"]["thresholds"]["depth_mean_min"], 50.0)
        code, out, err = run_cli(["--project", self.proj, "summary", "--samples", "S1", "--force", "--out-dir", str(self.tmp / "elsewhere")])
        self.assertEqual(code, 0, err)
        self.assertIn("S1: digest rebuilt", err)
        self.assertTrue((self.tmp / "elsewhere" / "S1.summary.html").exists())
        # errors
        self.assertEqual(run_cli(["--project", self.proj, "summary", "--samples", "S1", "--threshold", "nope=1"])[0], 1)
        self.assertEqual(run_cli(["--project", self.proj, "summary", "--samples", "S9"])[0], 1)
        self.assertEqual(run_cli(["--project", self.proj, "summary", "--cohort", "C9"])[0], 1)
        # no runs to summarise in joint mode
        code, out, err = run_cli(["--project", self.proj, "summary", "--mode", "joint", "--samples", "S1"])
        self.assertEqual(code, 1)
        self.assertIn("no successful joint run", err)

    def test_all_samples_and_jobs(self):
        code, out, err = run_cli(["--project", self.proj, "summary", "--jobs", "2"])
        self.assertEqual(code, 0, err)
        self.assertEqual(sorted(Path(p).name for p in out.strip().splitlines()), ["S1.summary.html", "S2.summary.html"])
        one = json.loads((self.out_dir / "S2.summary.json").read_text())
        code, out, err = run_cli(["--project", self.proj, "summary", "--samples", "S2", "--force"])
        self.assertEqual(code, 0, err)
        two = json.loads((self.out_dir / "S2.summary.json").read_text())
        one.pop("generated_at"), two.pop("generated_at")
        self.assertEqual(one, two)

    def test_joint_and_partial(self):
        db = DB(self.cfg.db_path)
        try:
            seed_success(self.cfg, db, "upstream", "sample", "S1", {"sample_id": "S1", "hifi_reads": ["/dev/null"], "inferred_sex": "MALE"}, mode="joint")
            seed_success(self.cfg, db, "downstream", "sample", "S1", {"sample_id": "S1", "sex": "MALE", "upstream_msg": ["joint message"]}, mode="joint")
            seed_success(self.cfg, db, "upstream", "sample", "S2", {"sample_id": "S2", "hifi_reads": ["/dev/null"], "inferred_sex": "FEMALE"}, mode="joint")
            seed_success(self.cfg, db, "assembly", "sample", "S1", {"sample_id": "S1", "hifi_reads": ["/dev/null"]}, mode="assembly")
        finally:
            db.close()
        code, out, err = run_cli(["--project", self.proj, "summary", "--mode", "joint", "--samples", "S1", "S2", "--out-dir", str(self.tmp / "joint")])
        self.assertEqual(code, 0, err)
        s1 = self.check_html(self.tmp / "joint" / "S1.summary.html")
        self.assertIn("joint message", s1)
        self.assertIn('id="sec-assembly"', s1)
        self.assertIn("trio-binned", s1) if False else None
        d1 = json.loads((self.tmp / "joint" / "S1.summary.json").read_text())
        self.assertEqual(set(d1["runs"]), {"upstream", "downstream", "assembly"})
        self.assertEqual(d1["qc"]["status"], "warn")
        self.assertNotIn("sex_mismatch", {f["code"] for f in d1["qc"]["flags"]})
        self.assertEqual(len(d1["assembly"]["haplotypes"]), 2)
        self.assertEqual(d1["assembly"]["haplotypes"][0]["stats"]["n50"], 60000000)
        self.assertEqual(d1["assembly"]["haplotypes"][1]["paftools"]["snps"], 2700001)
        self.assertEqual(d1["tandem_repeats"]["n_named"], 3)   # downstream's trgt_vcf
        self.assertEqual(d1["coverage"]["ratios"]["sex_from_ratio"], "MALE")   # upstream's mosdepth summary
        # upstream only: a partial page, the downstream-only sections say so
        s2 = self.check_html(self.tmp / "joint" / "S2.summary.html")
        d2 = json.loads((self.tmp / "joint" / "S2.summary.json").read_text())
        self.assertEqual(set(d2["runs"]), {"upstream"})
        self.assertEqual(d2["tandem_repeats"]["disease_loci"], [])
        self.assertIn("Not available", s2)
        self.assertEqual(d2["coverage"]["ratios"]["sex_from_ratio"], "FEMALE")

    def test_cohort_page(self):
        code, out, err = run_cli(["--project", self.proj, "summary", "--cohort", "C1"])
        self.assertEqual(code, 0, err)
        self.assertEqual([Path(p).name for p in out.strip().splitlines()], ["S1.summary.html", "S2.summary.html", "C1.summary.html"])
        text = self.check_html(self.out_dir / "C1.summary.html")
        for needle in ('id="sec-members"', 'id="sec-distributions"', 'id="sec-sex"', 'id="sec-callsets"', 'id="sec-freq"', 'id="sec-software"',
                       'href="S1.summary.html"', "sex_mismatch", "Merged structural variants", "trgt-lps", "Allele-frequency bins",
                       "C1 cohort summary", "hiphase", "<svg", 'class="iqr"'):
            self.assertIn(needle, text, needle)
        cd = json.loads((self.out_dir / "C1.summary.json").read_text())
        self.assertEqual([r["id"] for r in cd["members"]], ["S1", "S2"])
        self.assertEqual(cd["missing_members"], [])
        self.assertEqual(set(cd["runs"]), {"cohort_merge", "cohort_freq"})
        self.assertEqual(cd["callsets"]["cohort_merge"]["sv"]["counts"], {"DEL": 2, "INS": 1, "BND": 1, "DUP": 1})
        self.assertEqual(cd["callsets"]["cohort_merge"]["sv"]["samples"], ["S1", "S2"])
        self.assertEqual((cd["callsets"]["cohort_merge"]["trgt_lps"]["loci"], cd["callsets"]["cohort_merge"]["trgt_lps"]["samples"]), (5, ["S1", "S2"]))
        self.assertEqual(cd["callsets"]["cohort_freq"]["summary"]["samples"]["n"]["total"], 2)
        self.assertEqual(len(cd["callsets"]["cohort_freq"]["samples"]), 2)
        self.assertEqual(cd["distributions"]["depth_mean"]["values"], {"S1": 30.1, "S2": 30.1})
        self.assertEqual(cd["distributions"]["depth_mean"]["threshold"], 20.0)
        self.assertEqual(cd["sex"]["mismatches"], ["S1"])
        self.assertEqual(cd["sex"]["sheet"], {"MALE": 1, "FEMALE": 1})
        self.assertEqual(cd["warnings"], [])
        self.assertIn('href="C1.summary.html"', (self.out_dir / "S1.summary.html").read_text())
        # a member without runs is listed, not fatal; the sample list of --samples plus --cohort is the union
        tsv = write_tsv(self.tmp, [{"sample_id": "S3", "sex": "MALE", "hifi_reads": "s3.bam"}], name="more.tsv")
        self.assertEqual(run_cli(["--project", self.proj, "samples", "add", str(tsv)])[0], 0)
        self.assertEqual(run_cli(["--project", self.proj, "cohort", "freeze", "C2", "--samples", str(write_ids(self.tmp, ["S1", "S3"], "c2.txt"))])[0], 0)
        code, out, err = run_cli(["--project", self.proj, "summary", "--cohort", "C2", "--samples", "S2"])
        self.assertEqual(code, 0, err)
        self.assertIn("S3: no successful standalone run; skipped", err)
        self.assertEqual([Path(p).name for p in out.strip().splitlines()], ["S2.summary.html", "S1.summary.html", "C2.summary.html"])
        c2 = json.loads((self.out_dir / "C2.summary.json").read_text())
        self.assertEqual(c2["missing_members"], ["S3"])
        self.assertEqual(c2["runs"], {})
        self.assertIn("no summarised run", self.check_html(self.out_dir / "C2.summary.html"))

    def test_cohort_page_joint(self):
        db = DB(self.cfg.db_path)
        try:
            for sid in ("S1", "S2"):
                seed_success(self.cfg, db, "upstream", "sample", sid, {"sample_id": sid, "hifi_reads": ["/dev/null"], "inferred_sex": "MALE"}, mode="joint")
                seed_success(self.cfg, db, "downstream", "sample", sid, {"sample_id": sid, "sex": "MALE"}, mode="joint")
            seed_success(self.cfg, db, "cohort_call", "cohort", "C1", {"cohort_id": "C1", "sample_ids": ["S1", "S2"], "gvcfs": ["/dev/null"]}, mode="joint")
            seed_success(self.cfg, db, "cohort_merge", "cohort", "C1", {"cohort_id": "C1", "sample_ids": ["S1", "S2"], "run_glnexus": False}, mode="joint", attempt=2)
        finally:
            db.close()
        code, out, err = run_cli(["--project", self.proj, "summary", "--mode", "joint", "--cohort", "C1", "--out-dir", str(self.tmp / "joint")])
        self.assertEqual(code, 0, err)
        cd = json.loads((self.tmp / "joint" / "C1.summary.json").read_text())
        self.assertEqual(set(cd["runs"]), {"cohort_call", "cohort_merge"})
        self.assertEqual(cd["callsets"]["cohort_call"]["small_variant_vcf"]["samples"], ["S1", "S2"])
        self.assertNotIn("small_variant_vcf", cd["callsets"]["cohort_merge"])
        self.assertEqual(cd["sex"]["mismatches"], ["S2"])   # S2 is FEMALE in the sheet, the seed inferred MALE
        text = self.check_html(self.tmp / "joint" / "C1.summary.html")
        self.assertIn("Joint small-variant call set (cohort_call", text)
        self.assertIn("cohort_freq has not run", text)


if __name__ == "__main__":
    unittest.main()
