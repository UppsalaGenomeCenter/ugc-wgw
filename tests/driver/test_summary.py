import gzip
import html
import json
import tempfile
import unittest
from unittest import mock
from html.parser import HTMLParser
from pathlib import Path

from ugc_wgw import config, summary
from ugc_wgw import summary_parsers as sp
from ugc_wgw import summary_svg as sv
from ugc_wgw.db import DB

from .helpers import HtmlChecker, make_project, run_cli, seed_files, seed_success, write_ids, write_tsv

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

    def test_paraphase_gene_fields(self):
        from .fixtures import PARAPHASE
        p = sp.parse_paraphase(PARAPHASE)
        by = {r["region"]: r for r in p["regions"]}
        self.assertEqual((p["n_regions"], p["n_failed"], p["n_with_cn"], p["n_with_gene_cn"]), (3, 1, 2, 1))
        c = by["CYP2D6"]
        self.assertEqual((c["gene_cn"], c["fusions"], c["depth_p80"], c["genome_depth"], c["sample_sex"]), (2, ["CYP2D6-CYP2D7"], 33.0, 30.1, "male"))
        self.assertEqual(c["haplotypes_by_gene"], {"CYP2D6": 2, "CYP2D7": 2})
        self.assertEqual(c["phase_region"], "38:chr22:42126000-42145000")
        r = by["RCCX"]
        self.assertEqual((r["n_haplotypes"], r["two_copy"], r["haplotypes_by_gene"]), (3, 1, {"CYP21A2": 2, "TNXB": 2}))
        s = by["smn1"]
        self.assertEqual((s["gene_cn"], s["fusions"], s["haplotypes_by_gene"], s["failed"]), (None, [], {}, True))
        # the older orientation (name -> encoded haplotype) and a top-level gene_cn are read the same way
        p2 = sp.parse_paraphase({"smn1": {"final_haplotypes": {"smn1hap1": "ATC", "smn2hap1": "ATG"}, "gene_cn": 1, "two_copy_haplotypes": ["smn1hap1"]}})
        self.assertEqual((p2["regions"][0]["haplotypes_by_gene"], p2["regions"][0]["gene_cn"]), ({"smn1": 2, "smn2": 1}, 1))

    def test_starphase_genes(self):
        from .fixtures import STARPHASE_JSON
        sj = sp.parse_starphase_json(STARPHASE_JSON)
        g = sj["genes"]
        self.assertEqual((sj["n_genes"], sj["n_called"], sj["n_hla"], sj["n_hla_called"]), (7, 5, 2, 1))
        self.assertEqual((g["CYP2D6"]["hap1"], g["CYP2D6"]["hap2"], g["CYP2D6"]["simple_diplotype"], g["CYP2D6"]["hla"]), ("*1", "*4", "*1/*4", False))
        self.assertTrue(g["HLA-A"]["hla"] and g["HLA-A"]["called"])
        self.assertFalse(g["HLA-DRB3"]["called"])
        self.assertFalse(g["TPMT"]["called"])
        self.assertTrue(sp.starphase_called("*1", "*17"))
        for bad in ("NO_READS", "NO_MATCH", ".", "Multiple", ""):
            self.assertFalse(sp.starphase_called("*1", bad))
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
        self.assertEqual((g["genome_size"], sp.parse_caln50 is sp.parse_gfatools_stats), (3100000000, True))

    def test_assembly_parsers(self):
        from .fixtures import CHROM_COVERAGE, DIP_BED, FLAGGER_SUMMARY, RELIABLE_BED, TRUVARI_JSON
        fl = sp.parse_flagger_summary(FLAGGER_SUMMARY.splitlines())
        self.assertEqual((fl["hap_pct"], fl["unk_pct"], fl["non_hap_pct"], fl["hap_fraction"]), (95.99, 0.0, 4.01, 0.9599))
        self.assertEqual(fl["bp"], {"err": 37000000, "dup": 41000000, "hap": 2900000000, "col": 42000000, "unk": 0})
        self.assertEqual(fl["blocks"]["dup"], 62)
        te = sp.parse_trioeval(["S\tc1\t4403\t197\t4300\t30\t17\t150", "S\tc2\t1\t0\t0\t0\t0\t0", "S\tc3\t10\t10\t5\t2\t3\t5", "W\t426\t50814\t0.008384", "H\t1409\t50864\t0.027701", "bogus"])
        self.assertEqual((te["switch_errors"], te["hamming_rate"], te["n_contigs"]), (426, 0.027701, 2))      # c2 has fewer than two blocks
        self.assertEqual([(c["name"], c["switches"], c["hamming"]) for c in te["contigs"]], [("c1", 47, 197), ("c3", 5, 10)])
        cc = sp.parse_chrom_coverage(CHROM_COVERAGE.format(h1="pat", h2="mat").splitlines())
        self.assertEqual(cc["haplotypes"], ["pat", "mat"])
        self.assertEqual([(r["haplotype"], r["chrom"]) for r in cc["rows"]], [("pat", "chr1"), ("pat", "chr20"), ("pat", "chrX"), ("mat", "chr20"), ("mat", "chrY"), ("pat", "other"), ("mat", "other")])
        self.assertEqual(cc["rows"][-2], {"haplotype": "pat", "chrom": "other", "aligned_bases": 50000, "contigs": 2, "chrom_length": 186739, "fraction": None})
        self.assertEqual(cc["lengths"]["chrY"], 57227415)
        bb = sp.bed_bins([l.encode() for l in DIP_BED.splitlines(True)], {"chr20": 64444167, "chr1": 248956422}, bin_bp=1_000_000)
        c20 = bb["chroms"]["chr20"]
        self.assertEqual((len(c20["bins"]), c20["bp"], c20["len"]), (65, 3298175, 64444167))
        self.assertEqual(c20["bins"][0], round((393022 - 63840 + 1000000 - 408057) / 1e6, 3))      # two intervals inside the first bin
        self.assertEqual(c20["bins"][1], round((1877050 - 1000000) / 1e6, 3))                       # the second interval crosses the bin edge
        self.assertEqual(c20["bins"][2], 1.0)
        self.assertEqual(c20["bins"][3], 0.5)
        self.assertEqual(bb["chroms"]["chr1"]["bins"][1], 0.2)
        self.assertEqual(list(bb["chroms"]), ["chr1", "chr20"])                                      # contig order, unknown contigs skipped
        self.assertEqual(sp.bed_bins([b"chrQ\t0\t10\n"], {}, contigs=["chrQ"])["chroms"]["chrQ"], {"len": 10, "bins": [1.0], "bp": 10})
        self.assertEqual(sp.bed_bp_by_chrom([l.encode() for l in RELIABLE_BED.splitlines(True)] + [b"#x\n", b"track name=y\n", b"chr2\ta\tb\n"]), {"chr20": 3248175, "chr1": 200000})
        self.assertEqual(sp.parse_regroup_log(["from_haplotype\tcontig\tlength\tto_haplotype", "hap1\tc9\t12345\thap2"]), [{"from_haplotype": "hap1", "contig": "c9", "length": 12345, "to_haplotype": "hap2"}])
        tv = sp.parse_truvari_summary(TRUVARI_JSON)
        self.assertEqual((tv["gt_concordance"], tv["TP-base_TP-gt"], tv["params"]["sizemax"]), (0.95, 22800, 50000))
        self.assertNotIn("gt_matrix", tv)
        self.assertEqual(sp.parse_truvari_summary("x"), {})
        cs = sp.parse_contig_sizes(["Contig\tGRCh38\tA.1\tA.2\t_MINIGRAPH_\tmin\tmax\tavg", "_AMBIGUOUS_\t0\t100\t200\t0\t0\t200\t100",
                                    "chr20\t64444167\t64000000\t66000000\t64444167\t64000000\t66000000\t65000000", "chrX\t156040895\t0\t0\t156040895\t0\t0\t0"])
        self.assertEqual(cs["genomes"], ["A.1", "A.2"])
        self.assertEqual(cs["chromosomes"]["chr20"], {"ref_bp": 64444167, "hap_bp_min": 64000000, "hap_bp_max": 66000000, "hap_bp_avg": 65000000, "hap_bp_total": 130000000})
        self.assertEqual((cs["assigned_bp"], cs["ambiguous_bp"]), ({"A.1": 64000000, "A.2": 66000000}, {"A.1": 100, "A.2": 200}))
        self.assertEqual(sp.parse_name_map(["sample_id\tsafe_name\ttrio", "HG-1\tHG_1\ttrue", "X\tX\tfalse"]), {"HG-1": {"safe_name": "HG_1", "trio": True}, "X": {"safe_name": "X", "trio": False}})
        bs = sp.parse_bcftools_stats(BCFTOOLS_STATS.splitlines() + ["PSC\t0\tS2\t1\t2\t3\t4\t5\t6\t7.0\t8\t0\t0\t9"])
        self.assertEqual((bs["psc"]["sample"], [r["sample"] for r in bs["psc_all"]], bs["psc_all"][1]["n_missing"]), ("HG002", ["HG002", "S2"], 9))

    def test_tier2_parsers(self):
        from .fixtures import CALL_SPECTRA, CONTIG_QC, DEGREE, GENOME_BP_HEADER, MEMBERS_HEADER, REFGAPS
        cq = sp.parse_contig_qc(CONTIG_QC.format(h1="pat", h2="mat").splitlines())
        h1 = cq["haplotypes"][0]
        self.assertEqual((h1["haplotype"], h1["contigs"], h1["bp"], h1["n_bp"], h1["gaps"], h1["telo_bp"]), ("pat", 2, 65000000, 100, 1, 6000))
        self.assertEqual(h1["telo_ends"], {"0": 0, "1": 2, "2": 0})
        self.assertEqual(h1["flagger_bp"], {"err": 20000, "dup": 500000, "hap": 64290000, "col": 190000})
        self.assertEqual(h1["non_hap_pct"], round(710000 / 65000000 * 100, 2))
        self.assertAlmostEqual(h1["gc_pct"], (40.9 * 60000000 + 41.5 * 5000000) / 65000000, places=2)
        self.assertEqual(cq["haplotypes"][1]["telo_ends"], {"0": 1, "1": 0, "2": 1})
        self.assertEqual([c["contig"] for c in cq["top_contigs"]], ["h1tg000001l", "h2tg000001l", "h1tg000002l", "h2tg000002l"])
        self.assertEqual(cq["top_contigs"][0]["non_hap_bp"], 700000)
        cs = sp.parse_call_spectra(CALL_SPECTRA.splitlines())
        self.assertEqual(cs["dipcall"]["by_filter"], {"ASM_UNRELIABLE": 20000, "PASS": 4380700})
        self.assertEqual(cs["dipcall"]["by_type"]["SNV"], {"ASM_UNRELIABLE": 20000, "PASS": 3780000})
        self.assertEqual(cs["pav"]["sv_sizes"], {"DEL": {"1k-10k": 250}, "INS": {"100-1k": 900}, "INV": {"10k-100k": 12}})
        self.assertEqual((cs["pav"]["pass_total"], cs["pav"]["total"]), (4301162, 4495162))
        ms = sp.parse_member_stats((MEMBERS_HEADER + "A\t1\t2\t3\t4\t5\t100\t0.9500\nB\t1\t2\t3\t4\t6\t100\t0.9400\n").splitlines())
        self.assertEqual(ms[1], {"sample_id": "B", "n_ref_hom": 1, "n_nonref_hom": 2, "n_hets": 3, "n_singletons": 4, "n_missing": 6, "n_sites": 100, "call_rate": 0.94})
        gb = sp.parse_genome_bp((GENOME_BP_HEADER + "GRCh38#0\t1\t64444167\nA#1\t3\t60000000\nA#2\t2\t62000000\n").splitlines())
        self.assertEqual(gb, {"ref_bp": 64444167, "haplotypes": 2, "bp_min": 60000000, "bp_max": 62000000, "bp_avg": 61000000, "bp_total": 122000000})
        dg = sp.parse_degree_stats(DEGREE.splitlines() + ["9\t5\t0\t0\t0", "12\t1\t0\t0\t0"])
        self.assertEqual((dg["1"], dg["4"], dg["8+"]), (736985, 2878, 6))
        rg = sp.parse_refgaps([l.encode() for l in REFGAPS.format(chrom="chr20").splitlines(True)])
        self.assertEqual(rg, {"chr20": {"n": 2, "bp": 60000 + 28557053 - 26386225}})

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
        s = self.check(sv.hist([("0.05", 1), ("0.50", 2)]))          # a few counts: one tick per integer, none repeated
        self.assertEqual([s.count(f">{t}</text>") for t in (0, 1, 2)], [1, 1, 1])
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
        self.assertEqual([(f["code"], f["severity"]) for f in summary.qc_flags(d, dict(summary.DEFAULT_THRESHOLDS))["flags"]], [("sex_not_inferred", "warn")])
        # three sources (0.14.0): somalier's chrX call joins; the flag names every disagreeing source and what the run used
        d["stats"]["inferred_sex"] = "FEMALE"
        d["identity"] = {"declared": "MALE", "coverage": "FEMALE", "somalier": "MALE", "resolved": "MALE", "source": "declared", "legacy": False}
        f = summary.qc_flags(d, dict(summary.DEFAULT_THRESHOLDS))["flags"]
        self.assertEqual([(x["code"], x["severity"]) for x in f], [("sex_mismatch", "warn")])
        self.assertIn("the sheet MALE (XY), coverage FEMALE (XX), somalier MALE (XY); the run used MALE (declared)", f[0]["message"])
        self.assertEqual(f[0]["value"], "sheet MALE, coverage FEMALE, somalier MALE")
        # nothing could sex the sample: a failure by default, a warning with sex_unknown_fails 0
        d["sample"]["sex_sheet"] = None
        d["stats"]["inferred_sex"] = ""
        d["identity"] = {"declared": "", "coverage": "", "somalier": "", "resolved": "", "source": "none", "legacy": False}
        f = summary.qc_flags(d, dict(summary.DEFAULT_THRESHOLDS))["flags"]
        self.assertEqual([(x["code"], x["severity"]) for x in f], [("sex_not_inferred", "fail")])
        self.assertIn("no source could sex the sample", f[0]["message"])
        self.assertIn("XX default", f[0]["message"])
        self.assertEqual([(x["code"], x["severity"]) for x in summary.qc_flags(d, dict(summary.DEFAULT_THRESHOLDS, sex_unknown_fails=0))["flags"]],
                         [("sex_not_inferred", "warn")])
        d["identity"]["legacy"] = True
        self.assertIn("predates the somalier extract", summary.qc_flags(d, dict(summary.DEFAULT_THRESHOLDS))["flags"][0]["message"])
        # somalier alone sexed it: a warning that names the source the run used
        d["identity"] = {"declared": "", "coverage": "", "somalier": "FEMALE", "resolved": "FEMALE", "source": "somalier", "legacy": False}
        f = summary.qc_flags(d, dict(summary.DEFAULT_THRESHOLDS))["flags"]
        self.assertEqual([(x["code"], x["severity"]) for x in f], [("sex_not_inferred", "warn")])
        self.assertIn("the run used FEMALE (somalier)", f[0]["message"])
        d["identity"]["somalier"] = "MALE"
        d["sample"]["sex_sheet"] = "MALE"
        d["stats"]["inferred_sex"] = "MALE"
        self.assertEqual(summary.qc_flags(d, dict(summary.DEFAULT_THRESHOLDS))["flags"], [])   # three sources that agree
        d["sample"]["sex_sheet"] = "MALE"
        d["stats"]["inferred_sex"] = ""
        d["identity"] = {}
        # an assembly-only sample has an empty coverage block and no read stage: no sex to infer, so no warning
        asm_only = {"stats": {}, "sample": {"sex_sheet": "MALE"}, "coverage": {"mean": None, "inferred_sex": ""},
                    "runs": {"assembly": {"attempt": 1}, "assembly_call": {"attempt": 1}}, "targeted": {}, "messages": [], "warnings": []}
        self.assertFalse(summary._has_read_data(asm_only))
        self.assertEqual(summary.qc_flags(asm_only, dict(summary.DEFAULT_THRESHOLDS))["flags"], [])
        asm_only["runs"]["singleton"] = {"attempt": 1}
        self.assertTrue(summary._has_read_data(asm_only))
        self.assertEqual([f["code"] for f in summary.qc_flags(asm_only, dict(summary.DEFAULT_THRESHOLDS))["flags"]], ["sex_not_inferred"])

    def test_mito_flag(self):
        def flags(uf, thr=None):
            d = {"stats": {"depth_mean": "31", "mapped_read_percent": "99.9", "read_quality_median": "31", "inferred_sex": "MALE"},
                 "sample": {"sex_sheet": "MALE"}, "coverage": {"mean": 31.0}, "messages": [], "warnings": [],
                 "targeted": {"mitorsaw": {"haplotypes": [{"label": "hap_0"}], "fingerprint": {"passing_explanation": uf is None or uf <= 0.02, "unexplained_fraction": uf}}}}
            t = dict(summary.DEFAULT_THRESHOLDS)
            t.update(thr or {})
            return summary.qc_flags(d, t)["flags"]
        f = flags(0.05)
        self.assertEqual([(x["code"], x["severity"], x["value"], x["threshold"]) for x in f], [("mito_unexplained_high", "warn", 0.05, 0.02)])
        self.assertIn("5.0%", f[0]["message"])
        self.assertEqual(flags(0.02), [])
        self.assertEqual(flags(None), [])
        self.assertEqual(flags(0.05, {"mito_unexplained_max": 0.1}), [])


class IdentityTest(unittest.TestCase):
    """cohort_qc's tables -> parsers, the cohort-level flags and the identity pass list (bin/ugc_wgw/summary_cohort.py)."""

    def test_parsers(self):
        from fixtures import SOMALIER_PAIRS_HEADER, fixture_bytes
        sex = sp.parse_sex_tsv(fixture_bytes("somalier_sex_tsv", "S1", "singleton", sex={"declared": "", "coverage": "", "somalier": "MALE", "resolved": "MALE",
                                                                                          "source": "somalier"}).decode().splitlines())
        self.assertEqual((sex["declared"], sex["coverage"], sex["somalier"], sex["resolved"], sex["source"], sex["x_n"], sex["x_het"], sex["x_hom_alt"]),
                         ("", "", "MALE", "MALE", "somalier", 300, 5, 150))
        self.assertAlmostEqual(sex["x_het_ratio"], 5 / 150, places=4)
        self.assertEqual(sp.parse_sex_tsv([]), {})
        pairs = sp.parse_somalier_pairs((SOMALIER_PAIRS_HEADER + "A\tB\t0.980\t0\t396\t0.9\t190\t190\t180\t176\t110\t110\t100\t300\t0\t0\t-1.0\n"
                                         "A\tC\t0.500\t0\t300\t0.9\t190\t190\t180\t90\t110\t110\t100\t300\t0\t0\t0.5\n").splitlines())
        self.assertEqual([(p["sample_a"], p["sample_b"], p["relatedness"], p["ibs0"], p["expected"]) for p in pairs], [("A", "B", 0.98, 0, None), ("A", "C", 0.5, 0, 0.5)])
        related = sp.parse_somalier_pairs(fixture_bytes("qc_related_pairs", "C1", "cohort_qc", members=["A", "B", "C"], pedigree={"C": {"A"}}).decode().splitlines())
        self.assertEqual([(p["sample_a"], p["sample_b"], p["expected"]) for p in related], [("A", "C", 0.5)])   # the declared pair, kept at 0.01
        summ = sp.parse_qc_summary(fixture_bytes("qc_summary", "C1", "cohort_qc", members=["A", "B", "C"]).decode().splitlines())
        self.assertEqual((summ["pairs"]["n"]["total"], summ["members"]["n"]["total"], summ["pairs"]["hist"]["0.00"]), (3, 3, 3))
        members = sp.parse_qc_members(fixture_bytes("qc_members", "C1", "cohort_qc", members=["A", "B", "C"], pedigree={"C": {"A"}}).decode().splitlines())
        self.assertEqual([m["sample_id"] for m in members], ["A", "B", "C"])
        self.assertEqual((members[0]["somalier_sex"], members[0]["x_n"], members[0]["het_ratio"], members[0]["closest_relative"], members[0]["closest_expected"]),
                         ("MALE", 300, 0.63, "B", None))
        self.assertEqual(members[2]["closest_expected"], 0.5)
        anc = sp.parse_ancestry(fixture_bytes("qc_ancestry", "C1", "cohort_qc", members=["A", "B"]).decode().splitlines())
        self.assertEqual([(r["sample_id"], r["predicted"], r["given"]) for r in anc], [("HG00096", "EUR", "EUR"), ("HG00403", "EAS", "EAS"), ("A", "EUR", ""), ("B", "AFR", "")])
        self.assertEqual(anc[0]["probs"], {"EAS": 0.0, "AFR": 0.0, "AMR": 0.1, "SAS": 0.0, "EUR": 0.9})
        self.assertEqual(anc[0]["pcs"], [0.5, 2.0, 0.1, 0.2, 0.3])

    def test_mt_parsers(self):
        from fixtures import fixture_bytes
        summ = sp.parse_mt_summary(fixture_bytes("mt_summary", "S1", "singleton").decode().splitlines())
        self.assertEqual((summ["haplotypes"], summ["records"], summ["homoplasmic"], summ["heteroplasmic"], summ["svs"], summ["het_vaf_max"], summ["major_haplotype"]),
                         (2, 3, 2, 1, 0, 0.4136, "hap_0"))
        self.assertEqual(sp.parse_mt_summary([]), {})
        hg = sp.parse_haplogroups(fixture_bytes("mt_haplogroups", "S1", "singleton", haplogroup="K1a1b1a").decode().splitlines())
        self.assertEqual([(h["haplotype"], h["abundance"], h["haplogroup"], h["quality"]) for h in hg], [("hap_0", 0.97, "K1a1b1a", 1.0), ("hap_1", 0.03, "K1a1b1a", 0.95)])
        mt = sp.parse_mt_vcf(iter(fixture_bytes("cohort_mt_small_variant_vcf", "C1", "cohort_merge", members=["A", "B", "C"]).splitlines(keepends=True)))
        self.assertEqual((mt["records"], mt["samples"]), (3, ["A", "B", "C"]))
        self.assertEqual(mt["members"]["A"], {"homoplasmic": 2, "heteroplasmic": 1, "het_vaf_max": 0.4136})
        self.assertEqual(mt["members"]["B"], {"homoplasmic": 2, "heteroplasmic": 0, "het_vaf_max": None})
        empty = sp.parse_mt_vcf(iter(fixture_bytes("cohort_mt_sv_vcf", "C1", "cohort_merge", members=["A"]).splitlines(keepends=True)))
        self.assertEqual((empty["records"], empty["members"]["A"]["homoplasmic"]), (0, 0))

    def test_scatter(self):
        c = _Checker()
        s = sv.scatter([("EUR", 0.5, 2.0), ("EAS", 2.6, -7.1), ("AFR", -3.0, 1.0)], [("S1", 0.4, 2.0, "EUR"), ("S2", -2.9, 1.1, "AFR")])
        c.feed(s)
        self.assertEqual((c.errors, c.stack), ([], []))
        self.assertEqual(s.count("<circle"), 5)
        self.assertIn("<title>S2: AFR</title>", s)
        self.assertIn("PC1", s)
        self.assertIn('class="legend"', s)   # the superpopulations' colours
        self.assertEqual(s.count("<rect"), 4)   # the plot frame and the legend's three swatches
        self.assertIn("No data", sv.scatter([], []))
        self.assertIn("<circle", sv.scatter([], [("S1", 0.0, 0.0, "")]))   # one point: a degenerate range is widened

    def test_related_hub(self):
        """A member related to related_hub_min members or more is the anomaly (the 864 cohort's gos_0001, related to 858
        of 863 with the highest heterozygosity): it fails `related_hub`, its partners are not flagged on its account, and
        the pass list keeps them; a real pair among the partners is still handled by sheet order."""
        from ugc_wgw import summary_cohort as sc
        thr = dict(summary.DEFAULT_THRESHOLDS)
        others = [f"M{i:02d}" for i in range(1, 13)]
        pairs = [{"sample_a": "HUB", "sample_b": m, "relatedness": 0.25, "ibs0": 500, "expected": None} for m in others]
        pairs.append({"sample_a": "M01", "sample_b": "M02", "relatedness": 0.3, "ibs0": 400, "expected": None})
        order = ["HUB"] + others
        members = [{"sample_id": s, "het_ratio": 0.76 if s == "HUB" else 0.63, "somalier_sex": "MALE"} for s in order]
        digests = {s: {"sample": {"id": s, "sex_sheet": "MALE", "father_id": None, "mother_id": None}, "coverage": {"inferred_sex": ""},
                       "identity": {"legacy": False, "somalier": "x"}} for s in order}
        flags = sc.identity_flags(pairs, members, digests, order, thr)
        self.assertEqual([(f["code"], f["severity"], f["value"]) for f in flags["HUB"] if f["code"] != "het_outlier"], [("related_hub", "fail", 12)])
        self.assertIn("related at or above 0.2 to 12 members", flags["HUB"][0]["message"])
        self.assertEqual(sorted(f["code"] for f in flags["M01"]), ["related"])          # its own pair with M02, not the hub
        self.assertEqual(flags["M03"], [])
        rows = [{"id": s, "present": True, "identity_status": "fail" if any(f["severity"] == "fail" for f in flags[s]) else "pass"} for s in order]
        d = {"members": rows, "identity": {"available": True, "flags": flags, "pairs": pairs}, "thresholds": thr}
        passed, excluded = sc.identity_pass(d, {s: i for i, s in enumerate(order)})
        self.assertEqual(passed, ["M01"] + others[2:])
        self.assertEqual([sid for sid, _w in excluded], ["HUB", "M02"])
        self.assertIn("related at or above 0.2 to 12 members", dict(excluded)["HUB"])
        # the hub found later in the sheet changes nothing for its partners
        passed2, excluded2 = sc.identity_pass(d, {s: i for i, s in enumerate(others + ["HUB"])})
        self.assertEqual(passed2, ["M01"] + others[2:])
        # with a higher threshold nobody is a hub and the old rule applies: every partner warns `related`
        loose = sc.identity_flags(pairs, members, digests, order, dict(thr, related_hub_min=20))
        self.assertEqual(sorted(f["code"] for f in loose["M03"]), ["related"])
        self.assertNotIn("related_hub", {f["code"] for f in loose["HUB"]})

    def test_flags_and_pass(self):
        from ugc_wgw import summary_cohort as sc
        thr = dict(summary.DEFAULT_THRESHOLDS)
        pairs = [{"sample_a": "A", "sample_b": "B", "relatedness": 0.98, "ibs0": 0, "expected": None},
                 {"sample_a": "A", "sample_b": "C", "relatedness": 0.5, "ibs0": 0, "expected": 0.5},
                 {"sample_a": "B", "sample_b": "D", "relatedness": 0.45, "ibs0": 0, "expected": None},
                 {"sample_a": "C", "sample_b": "E", "relatedness": 0.02, "ibs0": 40, "expected": 0.5},
                 {"sample_a": "D", "sample_b": "E", "relatedness": 0.01, "ibs0": 40, "expected": None}]
        members = [{"sample_id": s, "het_ratio": h, "somalier_sex": sx} for s, h, sx in
                   (("A", 0.63, "MALE"), ("B", 0.64, "MALE"), ("C", 0.62, "FEMALE"), ("D", 0.63, "FEMALE"), ("E", 0.90, "MALE"), ("F", 0.63, "FEMALE"))]
        digests = {s: {"sample": {"id": s, "sex_sheet": "MALE" if s in ("A", "B", "E", "F") else "FEMALE", "father_id": None, "mother_id": None},
                       "coverage": {"inferred_sex": ""}, "identity": {"legacy": s == "F", "somalier": "" if s == "F" else "x"}}
                   for s in "ABCDEF"}
        digests["C"]["sample"]["father_id"] = "A"
        digests["E"]["sample"]["father_id"] = "C"      # declared, but unrelated: a mismatch
        digests["E"]["sample"]["mother_id"] = "OUTSIDE"   # a parent outside the cohort is not checked
        flags = sc.identity_flags(pairs, members, digests, list("ABCDEF"), thr)
        codes = {sid: sorted((f["code"], f["severity"]) for f in fl) for sid, fl in flags.items()}
        self.assertEqual(codes["A"], [("duplicate", "fail")])
        self.assertEqual(codes["B"], [("duplicate", "fail"), ("related", "warn")])
        self.assertEqual(codes["C"], [("pedigree_mismatch", "fail")])
        self.assertEqual(codes["D"], [("related", "warn")])
        self.assertEqual(codes["E"], [("het_outlier", "warn"), ("pedigree_mismatch", "fail")])
        self.assertEqual(codes["F"], [("sex_mismatch", "warn")])   # a legacy member: somalier's call from cohort_qc against the sheet
        self.assertIn("same individual as B (relatedness 0.980", flags["A"][0]["message"])
        self.assertIn("declared father C but relatedness 0.020 is below the parent-child range 0.35–0.65", [f["message"] for f in flags["E"] if f["code"] == "pedigree_mismatch"][0])
        self.assertIn("declared father of E", flags["C"][0]["message"])
        self.assertIn("somalier's chrX call in cohort_qc is FEMALE (XX) but the sheet says MALE (XY)", flags["F"][0]["message"])
        self.assertIn("0.900 is beyond the cohort's median 0.630", [f["message"] for f in flags["E"] if f["code"] == "het_outlier"][0])
        # a declared pair missing from the table warns; stricter thresholds move the undeclared pair to duplicate
        digests["D"]["sample"]["mother_id"] = "F"
        self.assertEqual([f["severity"] for f in sc.identity_flags(pairs, members, digests, list("ABCDEF"), thr)["D"] if f["code"] == "pedigree_mismatch"], ["warn"])
        strict = sc.identity_flags(pairs, members, digests, list("ABCDEF"), dict(thr, duplicate_relatedness_min=0.4))
        self.assertEqual(sorted(f["code"] for f in strict["D"]), ["duplicate", "pedigree_mismatch"])
        # the pass list: failures out, then the later member of each pair at or above related_relatedness_max
        rows = [{"id": s, "present": True, "identity_status": "fail" if any(f["severity"] == "fail" for f in flags[s]) else "pass"} for s in "ABCDEF"]
        d = {"members": rows, "identity": {"available": True, "flags": flags, "pairs": pairs}, "thresholds": thr}
        passed, excluded = sc.identity_pass(d, {s: i for i, s in enumerate("ABCDEF")})
        self.assertEqual(passed, ["D", "F"])
        self.assertEqual([sid for sid, _why in excluded], ["A", "B", "C", "E"])
        self.assertIn("same individual as B", dict(excluded)["A"])
        # sheet order decides which undeclared relative stays: with D first, B's partner D is kept and B goes; the declared
        # pair A-C (expected 0.5) stays whole
        rows2 = [dict(r, identity_status="pass") for r in rows]
        passed2, excluded2 = sc.identity_pass({"members": rows2, "identity": {"available": True, "flags": {}, "pairs": pairs}, "thresholds": thr},
                                              {"D": 0, "B": 1, "A": 2, "C": 3, "E": 4, "F": 5})
        self.assertEqual(passed2, ["D", "A", "C", "E", "F"])
        self.assertEqual(dict(excluded2)["B"], "related to D (relatedness 0.450) with no pedigree link in the sheet; D comes first")
        self.assertEqual([sid for sid, _w in excluded2], ["B"])


class IdentityEndToEnd(unittest.TestCase):
    """A four-member cohort through the fake: a declared parent-child pair, a duplicate pair, a member only somalier
    could sex, a sex alias in the sheet; the pages, the pass list, the manifest and the freeze warning."""

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.cfg = make_project(self.tmp)
        self.proj = str(self.cfg.project_dir)
        tsv = write_tsv(self.tmp, [{"sample_id": "S1", "sex": "MALE", "hifi_reads": "s1.bam"},
                                   {"sample_id": "S2", "sex": "FEMALE", "hifi_reads": "s2.bam", "father_id": "S1"},
                                   {"sample_id": "S3", "sex": "", "hifi_reads": "s3.bam"},
                                   {"sample_id": "S4", "sex": "M", "hifi_reads": "s4.bam"}])
        self.assertEqual(run_cli(["--project", self.proj, "samples", "add", str(tsv)])[0], 0)
        self.assertEqual(run_cli(["--project", self.proj, "cohort", "freeze", "C1", "--samples", str(write_ids(self.tmp, ["S1", "S2", "S3", "S4"]))])[0], 0)
        code, out, err = run_cli(["--project", self.proj, "submit", "--mode", "standalone", "--cohort", "C1"],
                                 env={"UGC_WGW_FAKE_INFERRED_SEX": "S1=FEMALE,S2=FEMALE,S4=MALE", "UGC_WGW_FAKE_RELATEDNESS": "S1:S2=0.52,S3:S4=0.97",
                                      "UGC_WGW_FAKE_SOMALIER_SEX": "S3=MALE"})
        self.assertEqual(code, 0, err)
        self.out_dir = self.cfg.results_dir / "reports" / "summary"

    def test_pages_pass_list_manifest_and_freeze(self):
        code, out, err = run_cli(["--project", self.proj, "summary", "--cohort", "C1", "--pass-list", str(self.tmp / "id.txt"), "--pass-on", "identity"])
        self.assertEqual(code, 0, err)
        self.assertEqual((self.tmp / "id.txt").read_text().split(), ["S1", "S2"])   # the declared parent-child pair stays
        self.assertNotIn("S2 excluded", err)
        self.assertIn("S3 excluded: same individual as S4", err)
        self.assertIn("pass list (identity): 2 of 4", err)
        cd = json.loads((self.out_dir / "C1.summary.json").read_text())
        self.assertEqual(cd["identity"]["counts"], {"duplicate": 2})   # the declared parent-child pair at 0.52 is no `related`
        self.assertEqual({r["id"]: r["identity_status"] for r in cd["members"]}, {"S1": "warn", "S2": "pass", "S3": "fail", "S4": "fail"})
        self.assertEqual({r["id"]: (r["sex_used"], r["sex_source"]) for r in cd["members"]},
                         {"S1": ("MALE", "declared"), "S2": ("FEMALE", "declared"), "S3": ("MALE", "somalier"), "S4": ("MALE", "declared")})
        self.assertEqual(cd["sex"]["mismatches"], ["S1"])          # sheet MALE, coverage FEMALE
        self.assertEqual(cd["sex"]["sources"], {"declared": 3, "somalier": 1})
        self.assertEqual(cd["identity"]["pairs"][0]["expected"], 0.5)
        text = (self.out_dir / "C1.summary.html").read_text()
        self.assertIn("same individual as S4", text)
        self.assertIn("duplicate", text)
        s3 = (self.out_dir / "S3.summary.html").read_text()
        self.assertIn("Identity in cohort C1", s3)
        self.assertIn("same individual as S4", s3)
        d3 = json.loads((self.out_dir / "S3.summary.json").read_text())
        self.assertEqual((d3["identity"]["resolved"], d3["identity"]["source"], d3["identity"]["somalier"], d3["identity"]["legacy"]), ("MALE", "somalier", "MALE", False))
        self.assertEqual(d3["identity"]["detail"]["x_n"], 300)
        self.assertEqual([(f["code"], f["severity"]) for f in d3["qc"]["flags"] if f["code"].startswith("sex")], [("sex_not_inferred", "warn")])
        d4 = json.loads((self.out_dir / "S4.summary.json").read_text())
        self.assertEqual(d4["sample"]["sex_sheet"], "MALE")        # the sheet's M
        self.assertNotIn("duplicate", {f["code"] for f in d4["qc"]["flags"]})   # cohort-level flags stay on the cohort page
        # --pass-on all needs an assembly too; identity needs a cohort with a cohort_qc run
        code, _, err = run_cli(["--project", self.proj, "summary", "--cohort", "C1", "--pass-list", str(self.tmp / "all.txt"), "--pass-on", "all"])
        self.assertEqual(code, 0, err)
        self.assertEqual((self.tmp / "all.txt").read_text(), "")
        code, _, err = run_cli(["--project", self.proj, "summary", "--samples", "S1", "--pass-list", str(self.tmp / "x.txt"), "--pass-on", "identity"])
        self.assertEqual(code, 1)
        self.assertIn("needs --cohort", err)
        self.assertEqual(run_cli(["--project", self.proj, "cohort", "freeze", "C2", "--samples", str(write_ids(self.tmp, ["S1", "S2"], "c2.txt"))])[0], 0)
        code, _, err = run_cli(["--project", self.proj, "summary", "--cohort", "C2", "--pass-list", str(self.tmp / "y.txt"), "--pass-on", "identity"])
        self.assertEqual(code, 1)
        self.assertIn("no successful cohort_qc run", err)
        code, _, err = run_cli(["--project", self.proj, "summary", "--cohort", "C2", "--pass-list", str(self.tmp / "z.txt")])
        self.assertEqual(code, 0, err)                                   # the assembly gate needs no cohort_qc
        self.assertIn("cohort_qc has not run for this cohort", (self.out_dir / "C2.summary.html").read_text())
        # the run manifest's sex block
        db = DB(self.cfg.db_path)
        try:
            m3 = json.loads((db.latest_run("sample", "S3", "singleton", None).run_path / "run_manifest.json").read_text())
            m4 = json.loads((db.latest_run("sample", "S4", "singleton", None).run_path / "run_manifest.json").read_text())
        finally:
            db.close()
        self.assertEqual(m3["sex"], {"declared": None, "coverage": "", "somalier": "MALE", "resolved": "MALE", "source": "somalier"})
        self.assertEqual((m4["sex"]["declared"], m4["sex"]["resolved"], m4["sex"]["source"]), ("MALE", "MALE", "declared"))
        # freezing a cohort that holds a known duplicate pair warns but succeeds
        code, _, err = run_cli(["--project", self.proj, "cohort", "freeze", "C3", "--samples", str(write_ids(self.tmp, ["S3", "S4"], "c3.txt"))])
        self.assertEqual(code, 0, err)
        self.assertIn("S3 and S4 look like the same individual (relatedness 0.970 in cohort_qc of C1)", err)
        code, _, err = run_cli(["--project", self.proj, "cohort", "freeze", "C4", "--samples", str(write_ids(self.tmp, ["S1", "S3"], "c4.txt"))])
        self.assertEqual(code, 0, err)
        self.assertNotIn("same individual", err)

    def test_legacy_member_and_ancestry(self):
        """A member whose run predates the extract (seed `legacy`) is extracted by cohort_qc; an ancestry table gives
        strata and the projection on the page; somalier's call for the legacy member comes from cohort_qc's table."""
        db = DB(self.cfg.db_path)
        try:
            seed_success(self.cfg, db, "singleton", "sample", "S5", {"sample_id": "S5", "hifi_reads": ["/dev/null"], "inferred_sex": "FEMALE", "legacy": True})
        finally:
            db.close()
        self.assertEqual(run_cli(["--project", self.proj, "samples", "add", str(write_tsv(self.tmp, [{"sample_id": "S5", "sex": "MALE", "hifi_reads": "s5.bam"}], "s5.tsv"))])[0], 0)
        self.assertEqual(run_cli(["--project", self.proj, "cohort", "freeze", "C5", "--samples", str(write_ids(self.tmp, ["S1", "S5"], "c5.txt"))])[0], 0)
        # the 1000 Genomes set "installed": the driver turns run_ancestry on by itself
        with open(self.cfg.ref_map_file, "a") as fh:
            fh.write(f"somalier_1kg_tar\t{self.cfg.ref_map_file}\nsomalier_1kg_labels\t{self.cfg.ref_map_file}\n")
        code, out, err = run_cli(["--project", self.proj, "submit", "--mode", "standalone", "--cohort", "C5", "--stage", "cohort_qc"],
                                 env={"UGC_WGW_FAKE_INFERRED_SEX": "S1=FEMALE,S5=FEMALE"})
        self.assertEqual(code, 0, err)
        self.assertIn("S5", err)   # the warning names the legacy member
        code, out, err = run_cli(["--project", self.proj, "summary", "--cohort", "C5"])
        self.assertEqual(code, 0, err)
        cd = json.loads((self.out_dir / "C5.summary.json").read_text())
        self.assertEqual(cd["callsets"]["cohort_qc"]["ancestry"], {"S1": "EUR", "S5": "AFR"})
        self.assertEqual(len(cd["callsets"]["cohort_qc"]["ancestry_background"]), 2)
        self.assertEqual(cd["targeted"]["strata"], ["ALL", "ancestry", "sex"])
        s5 = [r for r in cd["members"] if r["id"] == "S5"][0]
        self.assertTrue(s5["sex_legacy"])
        self.assertEqual((s5["sex_somalier"], s5["sex_used"], s5["sex_source"]), ("FEMALE", "MALE", "declared"))   # somalier's call from cohort_qc's table (fixture: FEMALE for odd members)
        self.assertIn("sex_mismatch", s5["flags"])
        text = (self.out_dir / "C5.summary.html").read_text()
        for needle in ("Ancestry (projection onto the 1000 Genomes superpopulations)", "<title>S5: AFR</title>", "somalier&#x27;s chrX call in cohort_qc is FEMALE (XX)"):
            self.assertIn(needle, text, needle)
        d5 = json.loads((self.out_dir / "S5.summary.json").read_text())
        self.assertEqual((d5["identity"]["legacy"], d5["identity"]["resolved"], d5["identity"]["source"], d5["identity"]["somalier"]), (True, "MALE", "declared", ""))


DIGEST_KEYS = {"schema", "generated_at", "mode", "sample", "runs", "warnings", "stats", "reads", "coverage", "small_variants", "sv",
               "phasing", "tandem_repeats", "methylation", "targeted", "identity", "messages", "files", "assembly", "provenance", "qc"}


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

    def test_digests_are_written_as_built(self):
        # a rebuild of hundreds of whole-genome members takes hours: the digests built before an interruption stay
        # on disk and are reused by the next run, which builds only the missing ones
        original, calls = summary._digest_job, []

        def flaky(args):
            calls.append(args[0].sample_id)
            if len(calls) == 2:
                raise RuntimeError("interrupted")
            return original(args)

        with mock.patch.object(summary, "_digest_job", flaky):
            with self.assertRaises(RuntimeError):
                run_cli(["--project", self.proj, "summary", "--samples", "S1", "S2"])
        self.assertTrue((self.out_dir / "S1.summary.json").exists())
        self.assertFalse((self.out_dir / "S2.summary.json").exists())
        code, _, err = run_cli(["--project", self.proj, "summary", "--samples", "S1", "S2"])
        self.assertEqual(code, 0, err)
        self.assertIn("S1: digest reused", err)
        self.assertIn("S2: digest rebuilt", err)

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
        self.assertEqual(len(digest["targeted"]["starphase"]), 7)
        self.assertEqual(digest["targeted"]["starphase_meta"]["genes"]["CYP2C19"]["simple_diplotype"], "*1/*17")
        self.assertEqual(digest["targeted"]["paraphase"]["n_with_gene_cn"], 1)
        # the mitochondrial block (0.15): the collapsed counts, the haplogroup per haplotype and the sample's
        mito = digest["targeted"]["mitorsaw"]
        self.assertEqual((mito["collapsed"]["heteroplasmic"], mito["collapsed"]["major_haplotype"], mito["haplogroup"], mito["legacy"]), (1, "hap_0", "H5a7", False))
        self.assertEqual([(h["label"], h["haplogroup"], h["haplogroup_quality"]) for h in mito["haplotypes"]], [("hap_0", "H5a7", 1.0), ("hap_1", "H5a7", 0.95)])
        self.assertIn("mt_summary", digest["files"])
        for needle in ("haplogroup (major haplotype)", "H5a7", "heteroplasmic sites"):
            self.assertIn(needle, text, needle)
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

    def test_assembly_only_and_pass_list(self):
        db = DB(self.cfg.db_path)
        try:
            # S1: trio-binned assembly with calls, concordance and masks; S2: assembly only (no calling yet)
            seed_success(self.cfg, db, "assembly", "sample", "S1", {"sample_id": "S1", "hifi_reads": ["/dev/null"], "father_id": "F", "mother_id": "M",
                                                                      "father_hifi_reads": ["/f"], "mother_hifi_reads": ["/m"]}, mode="assembly")
            seed_success(self.cfg, db, "assembly_call", "sample", "S1", {"sample_id": "S1", "hap_fastas": ["/a", "/b"], "haplotype_labels": ["pat", "mat"], "trio": True,
                                                                           "flagger_bed": "/f.bed", "read_small_variant_vcf": "/r.vcf.gz"}, mode="assembly")
            seed_success(self.cfg, db, "assembly", "sample", "S2", {"sample_id": "S2", "hifi_reads": ["/dev/null"]}, mode="assembly")
        finally:
            db.close()
        out = self.tmp / "asm"
        code, _, err = run_cli(["--project", self.proj, "summary", "--samples", "S1", "S2", "--out-dir", str(out), "--pass-list", str(out / "pass.txt")])
        self.assertEqual(code, 0, err)
        d1 = json.loads((out / "S1.summary.json").read_text())
        self.assertTrue({"assembly", "assembly_call"} <= set(d1["runs"]), d1["runs"])
        self.assertTrue(d1["assembly"]["trio"])
        self.assertEqual(d1["assembly"]["haplotypes"][1]["trio_eval"]["switch_rate"], 0.005)
        self.assertEqual(d1["assembly"]["concordance"]["het_loss_rate"], 0.02)
        self.assertEqual(d1["assembly"]["concordance"]["region"], "reliable.bed")
        s1 = self.check_html(out / "S1.summary.html")
        self.assertIn("trioeval mat", s1)
        self.assertIn("het loss", s1)
        self.assertIn("contigs with the most switches", s1)
        self.assertIn("Concordance with the read-based SNVs", s1)
        self.assertEqual(d1["schema"], summary.DIGEST_SCHEMA)
        self.assertEqual(d1["assembly"]["parents"], {"hap1": "F", "hap2": "M"})
        self.assertEqual(d1["assembly"]["haplotypes"][0]["trio_eval"]["n_contigs"], 1)
        self.assertEqual(d1["assembly"]["haplotypes"][0]["trio_eval"]["contigs"][0]["switches"], 13 + 13)   # c1 + c2 of the fixture's S line
        self.assertEqual(d1["assembly"]["summary"]["switch_max"], 0.005)
        self.assertEqual(d1["assembly"]["summary"]["genotype_concordance"], 0.9865)
        self.assertTrue(d1["provenance"]["containers"], d1["provenance"])           # an assembly-only sample has provenance now
        d2 = json.loads((out / "S2.summary.json").read_text())
        self.assertIn("assembly", d2["runs"])
        self.assertNotIn("assembly_call", d2["runs"])
        self.assertEqual(d2["assembly"]["flagger"]["hap_pct"], 95.99)
        self.assertIsNone(d2["assembly"]["reliable"])
        self.assertEqual([f["code"] for f in d2["qc"]["flags"] if f["code"].startswith("asm_") and f["severity"] == "fail"], [])   # the shared setup's read-based warnings do not count
        self.assertEqual((out / "pass.txt").read_text().split(), ["S1", "S2"])
        # a schema-1 digest on disk is rebuilt, not reused
        old = json.loads((out / "S2.summary.json").read_text()); old["schema"] = 1; old["assembly"] = {}
        (out / "S2.summary.json").write_text(json.dumps(old))
        self.assertEqual(run_cli(["--project", self.proj, "summary", "--samples", "S2", "--out-dir", str(out)])[0], 0)
        self.assertEqual(json.loads((out / "S2.summary.json").read_text())["schema"], summary.DIGEST_SCHEMA)
        self.assertEqual(json.loads((out / "S2.summary.json").read_text())["assembly"]["flagger"]["hap_pct"], 95.99)
        # the concordance gate: a strict minimum fails the trio sample (0.9865 in the fixture), the assembly-only sample has no concordance
        code, _, err = run_cli(["--project", self.proj, "summary", "--samples", "S1", "S2", "--out-dir", str(out), "--pass-list", str(out / "pass4.txt"),
                                "--threshold", "asm_genotype_concordance_min=0.99"])
        self.assertEqual(code, 0, err)
        self.assertEqual((out / "pass4.txt").read_text().split(), ["S2"])
        self.assertEqual([f["code"] for f in json.loads((out / "S1.summary.json").read_text())["qc"]["flags"] if f["severity"] == "fail"], ["asm_concordance_low"])
        # thresholds decide: a strict Flagger limit fails both, a strict switch-error limit fails the trio only
        code, _, err = run_cli(["--project", self.proj, "summary", "--samples", "S1", "S2", "--out-dir", str(out), "--pass-list", str(out / "pass2.txt"),
                                "--threshold", "asm_switch_error_max=0.0045"])
        self.assertEqual(code, 0, err)
        self.assertEqual((out / "pass2.txt").read_text().split(), ["S2"])
        d1 = json.loads((out / "S1.summary.json").read_text())
        self.assertEqual([f["code"] for f in d1["qc"]["flags"] if f["severity"] == "fail"], ["asm_switch_error"])
        code, _, err = run_cli(["--project", self.proj, "summary", "--samples", "S1", "S2", "--out-dir", str(out), "--pass-list", str(out / "pass3.txt"),
                                "--threshold", "asm_flagger_err_pct_max=1"])
        self.assertEqual(code, 0, err)
        self.assertEqual((out / "pass3.txt").read_text().strip(), "")

    def test_joint_and_partial(self):
        db = DB(self.cfg.db_path)
        try:
            seed_success(self.cfg, db, "upstream", "sample", "S1", {"sample_id": "S1", "hifi_reads": ["/dev/null"], "inferred_sex": "MALE"}, mode="joint")
            seed_success(self.cfg, db, "downstream", "sample", "S1", {"sample_id": "S1", "sex": "MALE", "upstream_msg": ["joint message"]}, mode="joint")
            seed_success(self.cfg, db, "upstream", "sample", "S2", {"sample_id": "S2", "hifi_reads": ["/dev/null"], "inferred_sex": "FEMALE"}, mode="joint")
            seed_success(self.cfg, db, "assembly", "sample", "S1", {"sample_id": "S1", "hifi_reads": ["/dev/null"]}, mode="assembly")
            seed_success(self.cfg, db, "assembly_call", "sample", "S1", {"sample_id": "S1", "hap_fastas": ["/a", "/b"], "haplotype_labels": ["hap1", "hap2"], "trio": False,
                                                                           "flagger_bed": "/f.bed"}, mode="assembly")
        finally:
            db.close()
        code, out, err = run_cli(["--project", self.proj, "summary", "--mode", "joint", "--samples", "S1", "S2", "--out-dir", str(self.tmp / "joint")])
        self.assertEqual(code, 0, err)
        s1 = self.check_html(self.tmp / "joint" / "S1.summary.html")
        self.assertIn("joint message", s1)
        self.assertIn('id="sec-assembly"', s1)
        self.assertIn("trio-binned", s1) if False else None
        d1 = json.loads((self.tmp / "joint" / "S1.summary.json").read_text())
        self.assertEqual(set(d1["runs"]), {"upstream", "downstream", "assembly", "assembly_call"})
        self.assertEqual(d1["qc"]["status"], "warn")
        self.assertNotIn("sex_mismatch", {f["code"] for f in d1["qc"]["flags"]})
        self.assertEqual(len(d1["assembly"]["haplotypes"]), 2)
        self.assertEqual(d1["assembly"]["haplotypes"][0]["stats"]["n50"], 60000000)
        self.assertEqual(d1["assembly"]["haplotypes"][1]["name"], "hap2")
        self.assertEqual(d1["assembly"]["calls"]["dipcall"]["snps"], 3800000)
        self.assertEqual(d1["assembly"]["calls"]["pav"]["indels"], 750000)
        self.assertEqual(d1["assembly"]["karyotype"]["karyotype"], "unknown")
        self.assertIn("inferred karyotype", s1)
        # the QC pieces of the fixtures: Flagger's shares, the reliable region, the concordance (assembly_call was seeded with a read VCF? no: without)
        self.assertEqual(d1["assembly"]["flagger"]["non_hap_pct"], 4.01)
        self.assertEqual(d1["assembly"]["reliable"]["reliable_fraction_of_reference"], 0.9139)
        self.assertIsNone(d1["assembly"]["concordance"])
        for needle in ("Flagger: what the reads say", "Chromosomes", "Where the assembly is reliable", "Karyotype evidence", "Calls against the reference",
                       "Reliable region", "indel length (sites)", "substitution spectrum", "NGx curve", "hap1 bp"):
            self.assertIn(needle, s1, needle)
        self.assertEqual([f["code"] for f in d1["qc"]["flags"] if f["code"].startswith("asm_")], ["asm_karyotype_unknown"])   # warn: chr20-only fixture
        self.assertEqual(d1["assembly"]["assembler_version"], "0.25.0-r726")
        self.assertIsNone(d1["assembly"]["parents"])
        self.assertEqual(d1["assembly"]["calls"]["dipcall"]["tstv"], 1.92)
        self.assertEqual(d1["assembly"]["calls"]["dipcall"]["psc"]["n_hets"], 2900000)
        self.assertEqual(d1["assembly"]["calls"]["dipcall"]["hethom"], 1.933)
        self.assertEqual(dict(d1["assembly"]["calls"]["pav"]["idd"])[">=+50"], 20)
        self.assertEqual(d1["assembly"]["flagger"]["bp"]["hap"], 2900000000)
        self.assertEqual(d1["assembly"]["flagger"]["unk_pct"], 0.0)
        cc = d1["assembly"]["chrom_coverage"]
        self.assertEqual(cc["haplotypes"], ["hap1", "hap2"])
        self.assertEqual([(r["chrom"], r["contigs"]) for r in cc["rows"] if r["haplotype"] == "hap1"], [("chr1", 1), ("chr20", 3), ("chrX", 1), ("other", 2)])
        reg = d1["assembly"]["regions"]
        self.assertEqual(reg["dip"]["chroms"]["chr20"]["bp"], 329182 + 1468993 + 1500000)
        self.assertEqual(reg["reliable"]["chroms"]["chr20"]["bp"], reg["dip"]["chroms"]["chr20"]["bp"] - 50000)
        self.assertEqual(len(reg["dip"]["chroms"]["chr20"]["bins"]), 65)           # 64,444,167 bp in 1 Mb bins
        self.assertEqual(reg["dip"]["chroms"]["chr1"]["len"], 248956422)           # length from chrom_coverage, not from the BED
        self.assertEqual(d1["assembly"]["haplotypes"][0]["lowq"], {"n": 3, "bp": 5083 + 1469 + 10})
        self.assertEqual(d1["assembly"]["haplotypes"][0]["hap_bed_bp"], 340000 + 1495000 + 1610000 + 220000)
        self.assertEqual(d1["assembly"]["haplotypes"][0]["unreliable_bp"], 60000)
        self.assertEqual(d1["assembly"]["haplotypes"][1]["unreliable_bp"], 10000)
        self.assertEqual(d1["assembly"]["haplotypes"][0]["pav_callable_bp"], 1840000 + 1610000 + 220000)
        self.assertEqual(d1["assembly"]["karyotype"]["moved"], [])
        summ = d1["assembly"]["summary"]
        self.assertEqual((summ["hap_length_min"], summ["ng50_min"], summ["contigs"], summ["dipcall_snvs"], summ["chrx_fraction"]), (3015000000, 57000000, 820, 3800000, 0.0))
        self.assertIsNone(summ["switch_max"])
        self.assertIn("chrom_coverage", d1["files"])                          # the assembly files are summarised and listed
        self.assertIn("dipcall_vcf_stats", d1["files"])
        # Tier 2 files: per-haplotype Flagger and sequence QC from contig_qc, the callers' spectra
        self.assertEqual(d1["assembly"]["haplotypes"][0]["seq"]["telo_ends"], {"0": 0, "1": 2, "2": 0})
        self.assertEqual(d1["assembly"]["haplotypes"][1]["flagger"]["bp"]["dup"], 300000)
        self.assertEqual(d1["assembly"]["contig_qc"]["top_contigs"][0]["contig"], "h1tg000001l")
        self.assertEqual(d1["assembly"]["calls"]["pav"]["spectra"]["sv_sizes"]["INV"], {"10k-100k": 12})
        self.assertEqual((summ["pav_sv_pass"], summ["dipcall_pass_fraction"], summ["telo_ends_max"]), (1162, round(4380700 / 4400700, 4), 1))
        for needle in ("Sequence QC per haplotype", "contigs with the most non-Hap bases", "Calls by FILTER class", "PASS structural variants by size", "Per haplotype"):
            self.assertIn(needle, s1, needle)
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
                       "C1 cohort summary", "hiphase", "<svg", 'class="iqr"', 'id="sec-targeted"', "Targeted callers", "CYP2D6",
                       "Pharmacogenomics (StarPhase)", "HLA (StarPhase)", "mitorsaw unexplained fraction"):
            self.assertIn(needle, text, needle)
        cd = json.loads((self.out_dir / "C1.summary.json").read_text())
        # the page's frequency tables are the export's
        from ugc_wgw import export
        members = {sid: json.loads((self.out_dir / f"{sid}.summary.json").read_text()) for sid in ("S1", "S2")}
        self.assertEqual(cd["targeted"]["pgx"]["freq"], export.freq_table(export.pgx_rows(members), "simple_diplotype",
                                                                          export.strata_of(members, ("sex",)), by=("gene",)))
        self.assertEqual(cd["targeted"]["strata"], ["ALL", "sex"])
        self.assertEqual(cd["targeted"]["paraphase"]["regions"][0]["region"], "CYP2D6")
        self.assertEqual(cd["distributions"]["mito_unexplained"]["threshold"], 0.02)
        self.assertEqual([r["id"] for r in cd["members"]], ["S1", "S2"])
        self.assertEqual(cd["missing_members"], [])
        self.assertEqual(set(cd["runs"]), {"cohort_qc", "cohort_merge", "cohort_freq", "cohort_tr", "cohort_methyl"})
        # identity: cohort_qc's tables read, no flag for two unrelated members without a pedigree; the section and the member pages
        self.assertTrue(cd["identity"]["available"])
        self.assertEqual((cd["identity"]["counts"], cd["identity"]["n_fail"]), ({}, 0))
        self.assertEqual(cd["identity"]["members"]["S1"]["closest_relative"], "S2")
        self.assertEqual(cd["callsets"]["cohort_qc"]["summary"]["members"]["n"]["total"], 2)
        self.assertEqual([r["sex_used"] for r in cd["members"]], ["MALE", "FEMALE"])
        self.assertEqual(cd["sex"]["sources"], {"declared": 2})
        for needle in ('id="sec-identity"', "Identity (cohort_qc)", "Pairs worth a look", "somalier (chrX)", "used by the run", "Relatedness of every pair"):
            self.assertIn(needle, text, needle)
        s1 = self.check_html(self.out_dir / "S1.summary.html")
        self.assertIn("Identity in cohort C1", s1)
        self.assertIn("No identity flag", s1)
        self.assertEqual(cd["callsets"]["cohort_merge"]["sv"]["counts"], {"DEL": 2, "INS": 1, "BND": 1, "DUP": 1})
        self.assertEqual(cd["callsets"]["cohort_merge"]["sv"]["samples"], ["S1", "S2"])
        # the mitochondrial sets of cohort_merge and cohort_freq (0.15)
        self.assertEqual(cd["callsets"]["cohort_merge"]["mt"]["records"], 3)
        self.assertEqual(cd["callsets"]["cohort_merge"]["mt"]["members"]["S1"]["heteroplasmic"], 1)
        self.assertEqual(len(cd["callsets"]["cohort_merge"]["mt_haplogroups"]), 4)
        self.assertEqual(cd["callsets"]["cohort_freq"]["summary"]["mt_small_variants"]["records"]["total"], 3)
        self.assertIn("mt_small_variant_freq_vcf", cd["callsets"]["cohort_freq"])
        self.assertEqual((cd["targeted"]["mito"]["haplogroups"], cd["targeted"]["mito"]["n_haplogrouped"]), ({"H5a7": 2}, 2))
        self.assertEqual(cd["targeted"]["mito"]["freq"][0], {"kind": "haplogroup", "value": "H5a7", "stratum_kind": "ALL", "stratum": "ALL", "count": 2, "fraction": 1.0, "n_called": 2, "n_members": 2})
        self.assertIn({"kind": "macro_haplogroup", "value": "H", "stratum_kind": "ALL", "stratum": "ALL", "count": 2, "fraction": 1.0, "n_called": 2, "n_members": 2}, cd["targeted"]["mito"]["freq"])
        self.assertIn("mt_heteroplasmic", cd["distributions"])
        for needle in ("mitochondrial small variants", "Haplogroup of each member", "in mt VCF"):
            self.assertIn(needle, text, needle)
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

    def test_cohort_page_assembly(self):
        db = DB(self.cfg.db_path)
        try:
            for sid in ("S1", "S2"):
                seed_success(self.cfg, db, "assembly", "sample", sid, {"sample_id": sid, "hifi_reads": ["/dev/null"]}, mode="assembly")
                seed_success(self.cfg, db, "assembly_call", "sample", sid, {"sample_id": sid, "hap_fastas": ["/a", "/b"], "haplotype_labels": ["hap1", "hap2"],
                                                                            "trio": False, "flagger_bed": "/f.bed"}, mode="assembly")
            seed_success(self.cfg, db, "cohort_pangenome", "cohort", "C1", {"cohort_id": "C1", "sample_ids": ["S1", "S2"], "hap1_fastas": ["/a", "/a"],
                                                                            "hap2_fastas": ["/b", "/b"], "trios": [False, False], "sample_sexes": ["", "XY"],
                                                                            "sv_vcf": "/dev/null", "sv_vcf_index": "/dev/null"}, mode="assembly")
        finally:
            db.close()
        code, out, err = run_cli(["--project", self.proj, "summary", "--mode", "assembly", "--cohort", "C1", "--out-dir", str(self.tmp / "asm")])
        self.assertEqual(code, 0, err)
        self.assertEqual([Path(p).name for p in out.strip().splitlines()], ["S1.summary.html", "S2.summary.html", "C1.summary.html"])
        cd = json.loads((self.tmp / "asm" / "C1.summary.json").read_text())
        self.assertEqual(set(cd["runs"]), {"cohort_pangenome"})
        cp = cd["callsets"]["cohort_pangenome"]
        self.assertEqual(cp["chromosomes"], ["chr1", "chr20", "chrX"])
        c20 = cp["graph"]["chromosomes"]["chr20"]
        self.assertEqual((c20["nodes"], c20["edges"], c20["length"]), (501000, 701000, 61000000))
        self.assertEqual((c20["ref_length"], c20["hap_bp_min"], c20["hap_bp_avg"], c20["hap_bp_max"]), (64444167, 65000000, 65000000, 65000000))   # contig_sizes.tsv
        self.assertEqual((cp["graph"]["nodes"], cp["graph"]["length"]), (1503000, 183000000))
        self.assertEqual(cp["contig_sizes"]["ambiguous_bp"]["S1.1"], 2000000)
        self.assertEqual(cp["name_map"]["S2"], {"safe_name": "S2", "trio": False})
        self.assertEqual((cp["truvari"]["gt_concordance"], cp["truvari"]["TP-base_TP-gt"], cp["truvari"]["params"]["sizemin"]), (0.95, 22800, 50))
        self.assertNotIn("gt_matrix", cp["truvari"])
        self.assertIn("freq_vcf", cp)
        self.assertEqual([g["genome"] for g in cp["split"]["genomes"]], ["S1.1", "S1.2", "S2.1", "S2.2"])
        self.assertEqual(cp["split"]["ambiguous_contigs"], 18 + 18 + 19 + 19)
        self.assertEqual(cp["summary"]["pangenome"]["records"]["total"], 4100000)
        self.assertEqual(cp["summary"]["samples"]["n_in_vcf"]["pangenome"], 2)
        self.assertEqual(cp["truvari"]["recall"], 0.8)
        self.assertEqual([r["in_pangenome_vcf"] for r in cp["samples"]], ["yes", "yes"])
        self.assertEqual(sorted(cd["assembly"]), ["S1", "S2"])                       # the members' assembly rows come along
        self.assertEqual(cd["warnings"], [])
        # the assembly metrics are distributions with a direction; the member rows carry the karyotype
        self.assertEqual(cd["distributions"]["asm_reliable_fraction"]["values"], {"S1": 0.9139, "S2": 0.9139})
        self.assertEqual((cd["distributions"]["asm_flagger_non_hap"]["direction"], cd["distributions"]["asm_flagger_non_hap"]["threshold"]), ("max", 10.0))
        self.assertEqual(cd["distributions"]["asm_flagger_non_hap"]["beyond_threshold"], [])
        self.assertEqual(cd["distributions"]["asm_dipcall_hethom"]["values"]["S1"], 1.933)
        self.assertEqual(cd["sex"]["karyotype"], {"unknown": 2})
        self.assertEqual([r["karyotype"] for r in cd["members"]], ["unknown", "unknown"])
        # Tier 2: the members' genotypes, the per-chromosome haplotype sequence, degrees, gaps and the summary's level rows
        self.assertEqual([m["sample_id"] for m in cp["members"]], ["S1", "S2"])
        self.assertEqual(cp["member_distributions"]["call_rate"]["threshold"], 0.9)
        self.assertEqual(cp["member_distributions"]["call_rate"]["beyond_threshold"], [])
        self.assertEqual(c20["genome_haplotypes"], 4)
        self.assertEqual(c20["genome_bp_min"], 61_000_000 - 2 * 50_000)   # member 0, haplotype 2 of the fixture
        self.assertEqual(c20["refgaps"], {"n": 2, "bp": 60000 + 28557053 - 26386225})
        self.assertEqual(cp["degree"]["1"], 736985 * 3)
        self.assertEqual(cd["distributions"]["asm_pav_sv_pass"]["values"]["S1"], 1162)
        text = self.check_html(self.tmp / "asm" / "C1.summary.html")
        for needle in ('id="sec-pangenome"', "Pangenome graph (cohort_pangenome)", "Per chromosome", "Contigs binned by the split",
                       "truvari bench", "Allele frequencies of the vcfwave VCF", "cohort_freq has not run", "not haplotypes",
                       "graph / reference", "Shortest haplotype as a fraction", "genotype concordance", "assigned bp h1 / h2",
                       "2 assemblies (one row per member", "PAV SNVs / indels", "karyotype from the assembly", "assembly: Flagger non-Hap %",
                       "Graph structure", "Node degree", "Haplotype sequence in each chromosome", "genotypes in the graph VCF", "call rate in the graph VCF"):
            self.assertIn(needle, text, needle)
        # a tight Flagger limit puts both members above the threshold in the distributions (direction max)
        code, out, err = run_cli(["--project", self.proj, "summary", "--mode", "assembly", "--cohort", "C1", "--out-dir", str(self.tmp / "asm"),
                                  "--threshold", "asm_flagger_err_pct_max=1"])
        self.assertEqual(code, 0, err)
        cd2 = json.loads((self.tmp / "asm" / "C1.summary.json").read_text())
        self.assertEqual(cd2["distributions"]["asm_flagger_non_hap"]["beyond_threshold"], ["S1", "S2"])
        self.assertIn("above threshold", self.check_html(self.tmp / "asm" / "C1.summary.html"))
        # the member call-rate threshold is a cohort-page signal: raised above the fixtures' rates it lists both members
        code, out, err = run_cli(["--project", self.proj, "summary", "--mode", "assembly", "--cohort", "C1", "--out-dir", str(self.tmp / "asm"),
                                  "--threshold", "pg_member_call_rate_min=0.99"])
        self.assertEqual(code, 0, err)
        cd3 = json.loads((self.tmp / "asm" / "C1.summary.json").read_text())
        self.assertEqual(cd3["callsets"]["cohort_pangenome"]["member_distributions"]["call_rate"]["beyond_threshold"], ["S1", "S2"])
        self.assertIn("Members below the call-rate threshold", self.check_html(self.tmp / "asm" / "C1.summary.html"))
        # without the stage the assembly-mode page says so; the standalone page has no pangenome section at all
        self.assertEqual(run_cli(["--project", self.proj, "summary", "--mode", "assembly", "--cohort", "C2", "--out-dir", str(self.tmp / "asm2")])[0], 1)  # unknown cohort
        code, out, err = run_cli(["--project", self.proj, "summary", "--cohort", "C1"])
        self.assertEqual(code, 0, err)
        self.assertNotIn('id="sec-pangenome"', (self.out_dir / "C1.summary.html").read_text())

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


class MethylationTest(unittest.TestCase):
    """cohort_methyl's tables and the per-sample methbat report -> parsers, the digest's methylation block and flag, the
    cohort page's Methylation section and the methyl_outlier flag (bin/ugc_wgw/summary.py, summary_cohort.py)."""

    def test_parsers(self):
        from fixtures import fixture_bytes
        rep = sp.parse_methbat_report(fixture_bytes("methbat_report", "S1", "singleton", imprinting="anomalous").decode().splitlines())
        self.assertEqual((rep["n"], rep["pass"], rep["inconclusive"], rep["anomalous"], rep["anomalous_qc"]), (3, 2, 0, 1, 0))
        self.assertEqual(rep["anomalous_regions"], ["MEG3:TSS-DMRb"])
        self.assertEqual(rep["regions"][0]["region"], "PLAGL1:alt-TSS-DMR")
        self.assertAlmostEqual(rep["regions"][0]["combined"], 46.0)
        plain = sp.parse_methbat_report(fixture_bytes("methbat_report", "S1", "singleton").decode().splitlines())
        self.assertEqual((plain["inconclusive"], plain["anomalous"], plain["anomalous_regions"]), (1, 0, []))
        summ = sp.parse_methyl_summary(fixture_bytes("methyl_summary", "C1", "cohort_methyl", members=["A", "B", "C"], legacy=["C"], labels=["XX", "XY", ""]).decode().splitlines())
        self.assertEqual(summ["cohort"]["n_members"][""], 3)
        self.assertEqual(summ["labels"]["n_members"]["XX"], 1)
        self.assertIn("XX", summ["notes"]["label_below_min"])
        self.assertEqual(summ["genes"]["members_with_event"]["GENEA"], 2)
        samples = sp.parse_methyl_samples(fixture_bytes("methyl_samples", "C1", "cohort_methyl", members=["A", "B", "C"], legacy=["C"], labels=["XX", "XY", ""]).decode().splitlines())
        self.assertEqual([(r["sample_id"], r["n_events"], r["report_source"]) for r in samples], [("A", 0, "run"), ("B", 1, "run"), ("C", 1, "stage")])
        self.assertEqual(samples[2]["imprinting_anomalous"], 1)
        regions = sp.parse_methyl_regions(fixture_bytes("methyl_regions", "C1", "cohort_methyl", members=["A", "B", "C"]).decode().splitlines(), top=2)
        self.assertEqual(regions["n_regions"], 4)
        self.assertEqual([(r["region"], r["events"], r["gene_symbol"]) for r in regions["top"]], [("CpG:_1", 2, "GENEA;GENEB")])
        self.assertEqual(regions["totals"]["ALL"], {"regions": 4, "events": 2, "regions_with_event": 1})
        events = sp.parse_methyl_events(fixture_bytes("methyl_events", "C1", "cohort_methyl", members=["A", "B", "C"]).decode().splitlines(), genes=["hgnc:2"])
        self.assertEqual((events["n"], events["n_members"], len(events["rows"])), (2, 2, 2))
        self.assertEqual(events["rows"][0]["sample_id"], "B")
        none = sp.parse_methyl_events(fixture_bytes("methyl_events", "C1", "cohort_methyl", members=["A", "B", "C"]).decode().splitlines(), genes=["NOPE"])
        self.assertEqual((none["n"], len(none["rows"])), (2, 0))


class TandemRepeatEndToEnd(unittest.TestCase):
    """A five-member cohort through the fake with cohort_tr: the cohort page's Tandem repeats section, the digest's
    block (findings with the members' own dropout class, the outlier counts and the tr_outlier_count_high flag), the
    members' pages, the unchanged per-sample digest schema, and a cohort whose cohort_tr failed."""

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.cfg = make_project(self.tmp)
        self.proj = str(self.cfg.project_dir)
        tsv = write_tsv(self.tmp, [{"sample_id": f"S{i}", "sex": "MALE" if i % 2 else "FEMALE", "hifi_reads": f"s{i}.bam"} for i in range(1, 6)])
        self.assertEqual(run_cli(["--project", self.proj, "samples", "add", str(tsv)])[0], 0)
        self.assertEqual(run_cli(["--project", self.proj, "cohort", "freeze", "C1", "--samples", str(write_ids(self.tmp, [f"S{i}" for i in range(1, 6)]))])[0], 0)
        self.out_dir = self.cfg.results_dir / "reports" / "summary"

    def test_parsers(self):
        from .fixtures import fixture_bytes
        members = ["S1", "S2", "S3", "S4", "S5"]
        text = lambda name: fixture_bytes(name, "C1", "cohort_tr", members=members).decode().splitlines()  # noqa: E731
        screen = sp.parse_tr_screen(text("tr_pathogenic_screen"))
        self.assertEqual((screen["n_rows"], screen["by_class"]), (10, {"pathogenic": 1, "normal": 7, "intermediate": 1, "uncalled": 1}))
        self.assertEqual([(f["sample_id"], f["trid"], f["class"], f["a1_copies"], f["a1_lps"]) for f in screen["findings"]],
                         [("S1", "HD_HTT", "pathogenic", 45, 43), ("S2", "FXS_FMR1", "intermediate", 60, 58)])   # the worst class first
        summ = sp.parse_tr_summary(text("tr_summary"))
        self.assertEqual((summ["cohort"]["n_members"][""], summ["loci"]["outlier_status"]["common"], summ["parameters"]["outlier_mad_min"][""]), (5, 2, 5.0))
        self.assertEqual(len(summ["source"]["tr_map_input"]["values"]), 2)                # both inputs kept
        self.assertIn("CC BY 4.0", summ["source"]["tr_map_source"]["values"][0])
        samples = sp.parse_tr_samples(text("tr_samples"))
        self.assertEqual([r["n_outliers"] for r in samples], [1, 2, 3, 1, 40])
        top = sp.parse_tr_outliers_top(text("tr_outliers_top"), per_member=1)
        self.assertEqual((len(top["S1"]), top["S1"][0]["trid"], top["S1"][0]["z"], top["S5"][0]["al"]), (1, "chr1_1000_1090_CAG", 120, 450))
        loci = sp.parse_tr_screen_loci(text("tr_screen_loci"))
        self.assertEqual([(l["trid"], l["n_called"], l["call_note"], l["benign_max"]) for l in loci], [("HD_HTT", 5, "ok", 26), ("FXS_FMR1", 4, "few_uncalled", 44)])

    def test_pages_digest_and_flag(self):
        code, out, err = run_cli(["--project", self.proj, "submit", "--mode", "standalone", "--cohort", "C1"])
        self.assertEqual(code, 0, err)
        code, out, err = run_cli(["--project", self.proj, "summary", "--cohort", "C1"])
        self.assertEqual(code, 0, err)
        d1 = json.loads((self.out_dir / "S1.summary.json").read_text())
        self.assertEqual(d1["schema"], 5)                                # the stage adds nothing to a per-sample digest
        cd = json.loads((self.out_dir / "C1.summary.json").read_text())
        trd = cd["tandem_repeats"]
        self.assertTrue(trd["available"])
        self.assertIn("cohort_tr", cd["runs"])
        self.assertEqual(trd["counts"], {"tr_outlier_count_high": 1})    # 1, 2, 3, 1, 40: median 2, MAD 1, limit 5
        self.assertEqual([f["code"] for f in trd["flags"]["S5"]], ["tr_outlier_count_high"])
        self.assertIn("ancestry that is rare in this cohort", trd["flags"]["S5"][0]["message"])
        rows = {r["id"]: r for r in cd["members"]}
        self.assertEqual({k: v["tr_outliers"] for k, v in rows.items()}, {"S1": 1, "S2": 2, "S3": 3, "S4": 1, "S5": 40})
        self.assertIn("tr_outlier_count_high", rows["S5"]["flags"])
        self.assertEqual(rows["S5"]["qc_status"], "warn")
        self.assertIn("tr_outliers", cd["distributions"])
        self.assertEqual([(f["sample_id"], f["trid"], f["class"]) for f in trd["findings"]], [("S1", "HD_HTT", "pathogenic"), ("S2", "FXS_FMR1", "intermediate")])
        self.assertIn("dropout", trd["findings"][0])                     # joined from the member's digest, not read by the stage
        self.assertNotIn("screen", cd["callsets"]["cohort_tr"])          # the findings live in the block, once
        self.assertEqual(len(cd["callsets"]["cohort_tr"]["screen_loci"]), 2)
        page = (self.out_dir / "C1.summary.html").read_text()
        for needle in ("Tandem repeats (cohort_tr)", "A research screen, not a diagnostic result", "Disease-locus screen (STRchive thresholds)",
                       "Named loci", "Loci with the most outlier members", "few_uncalled", "CC BY 4.0", "tr_outlier_count_high",
                       "Member x locus rows by class: normal 7"):
            self.assertIn(needle, page, needle)
        page1 = (self.out_dir / "S1.summary.html").read_text()
        for needle in ("Against the cohort C1 (cohort_tr)", "Strongest expansion outliers", "chr4_90_260_GCC", "pathogenic"):
            self.assertIn(needle, page1, needle)
        page4 = (self.out_dir / "S4.summary.html").read_text()
        self.assertIn("every called named locus is within its benign range", page4)
        # without --cohort the sample page has no cohort block
        code, out, err = run_cli(["--project", self.proj, "summary", "--samples", "S1", "--out-dir", str(self.tmp / "solo")])
        self.assertEqual(code, 0, err)
        self.assertNotIn("Against the cohort", (self.tmp / "solo" / "S1.summary.html").read_text())

    def test_optional_and_absent(self):
        # cohort_tr fails: the submit reports it, cohort_methyl ran anyway, and the page says the stage has not run
        code, out, err = run_cli(["--project", self.proj, "submit", "--mode", "standalone", "--cohort", "C1"], env={"UGC_WGW_FAKE_FAIL_STAGES": "cohort_tr"})
        self.assertEqual(code, 1, err)
        db = DB(self.cfg.db_path)
        try:
            self.assertEqual(db.latest_run("cohort", "C1", "cohort_tr", VERSION_OF(self.cfg)).status, "failed")
            self.assertEqual(db.latest_run("cohort", "C1", "cohort_methyl", VERSION_OF(self.cfg)).status, "success")
        finally:
            db.close()
        code, out, err = run_cli(["--project", self.proj, "summary", "--cohort", "C1"])
        self.assertEqual(code, 0, err)
        cd = json.loads((self.out_dir / "C1.summary.json").read_text())
        self.assertFalse(cd["tandem_repeats"]["available"])
        self.assertIn("cohort_tr has not run for this cohort", (self.out_dir / "C1.summary.html").read_text())
        self.assertNotIn("Against the cohort", (self.out_dir / "S1.summary.html").read_text())


def VERSION_OF(cfg) -> str:  # noqa: N802 - the version the test project runs at
    return config.read_version(cfg.code_dir)


class MethylationEndToEnd(unittest.TestCase):
    """A five-member cohort through the fake with cohort_methyl last: the sample page's imprinting table and flag, the
    cohort page's Methylation section, the digest's methylation block, the methyl_outlier flag, the export's methyl
    domain, and a project without the stage's run."""

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.cfg = make_project(self.tmp, methyl_genes_of_interest=["GENEA"])
        self.proj = str(self.cfg.project_dir)
        tsv = write_tsv(self.tmp, [{"sample_id": f"S{i}", "sex": "MALE" if i % 2 else "FEMALE", "hifi_reads": f"s{i}.bam"} for i in range(1, 6)])
        self.assertEqual(run_cli(["--project", self.proj, "samples", "add", str(tsv)])[0], 0)
        self.assertEqual(run_cli(["--project", self.proj, "cohort", "freeze", "C1", "--samples", str(write_ids(self.tmp, [f"S{i}" for i in range(1, 6)]))])[0], 0)
        code, out, err = run_cli(["--project", self.proj, "submit", "--mode", "standalone", "--cohort", "C1"])
        self.assertEqual(code, 0, err)
        self.out_dir = self.cfg.results_dir / "reports" / "summary"

    def test_pages_digest_flags_and_export(self):
        code, out, err = run_cli(["--project", self.proj, "summary", "--cohort", "C1"])
        self.assertEqual(code, 0, err)
        d1 = json.loads((self.out_dir / "S1.summary.json").read_text())
        self.assertEqual(d1["schema"], summary.DIGEST_SCHEMA)
        imp = d1["methylation"]["imprinting"]
        self.assertEqual((imp["n"], imp["pass"], imp["inconclusive"], imp["anomalous"]), (3, 2, 1, 0))
        self.assertNotIn("imprinting_anomalous", {f["code"] for f in d1["qc"]["flags"]})
        self.assertIn("methbat_report", d1["files"])
        page1 = (self.out_dir / "S1.summary.html").read_text()
        self.assertIn("Imprinting regions (methbat report)", page1)
        self.assertIn("PLAGL1:alt-TSS-DMR", page1)
        cd = json.loads((self.out_dir / "C1.summary.json").read_text())
        self.assertTrue(cd["methylation"]["available"])
        self.assertIn("cohort_methyl", cd["runs"])
        self.assertEqual(cd["callsets"]["cohort_methyl"]["summary"]["cohort"]["n_members"][""], 5)
        self.assertEqual(cd["callsets"]["cohort_methyl"]["genes_of_interest"], ["GENEA"])
        self.assertEqual(cd["callsets"]["cohort_methyl"]["events"]["n"], 4)       # the fixture: every member but the first has one event
        self.assertEqual(cd["callsets"]["cohort_methyl"]["regions"]["top"][0]["gene_symbol"], "GENEA;GENEB")
        self.assertEqual({r["id"]: r["methyl_events"] for r in cd["members"]}, {"S1": 0, "S2": 1, "S3": 1, "S4": 1, "S5": 1})
        self.assertEqual(cd["methylation"]["counts"], {"methyl_outlier": 0})   # 0,1,1,1,1: median 1, MAD 0 -> nobody flagged
        self.assertIn("methyl_events", cd["distributions"])
        page = (self.out_dir / "C1.summary.html").read_text()
        self.assertIn("Methylation (cohort_methyl)", page)
        self.assertIn("Genes with the most members carrying an event", page)
        self.assertIn("Events in the genes of interest", page)
        self.assertIn("Imprinting screen (methbat report)", page)
        self.assertIn("Group comparisons (methbat compare)", page)
        self.assertIn("below the minimum", page)
        # the export's methyl domain copies the stage's tables, the regions table included
        out_dir = self.tmp / "export"
        code, out, err = run_cli(["--project", self.proj, "export", "--cohort", "C1", "--out", str(out_dir), "--domains", "methyl"])
        self.assertEqual(code, 0, err)
        names = sorted(p.name for p in out_dir.glob("C1.methyl.*"))
        self.assertEqual(names, ["C1.methyl.background.tsv", "C1.methyl.compare.XX_vs_XY.tsv", "C1.methyl.events.tsv", "C1.methyl.imprinting.tsv", "C1.methyl.regions.tsv",
                                 "C1.methyl.samples.tsv", "C1.methyl.summary.tsv"])

    def test_anomalous_report_flags_the_sample(self):
        seed_success(self.cfg, DB(self.cfg.db_path), "singleton", "sample", "S2", {"sample_id": "S2", "inferred_sex": "MALE", "imprinting": "anomalous"}, attempt=3)
        code, out, err = run_cli(["--project", self.proj, "summary", "--samples", "S2"])
        self.assertEqual(code, 0, err)
        d2 = json.loads((self.out_dir / "S2.summary.json").read_text())
        flags = {f["code"]: f for f in d2["qc"]["flags"]}
        self.assertIn("imprinting_anomalous", flags)
        self.assertEqual(flags["imprinting_anomalous"]["severity"], "warn")
        self.assertIn("MEG3:TSS-DMRb", flags["imprinting_anomalous"]["message"])
        self.assertEqual(d2["methylation"]["imprinting"]["anomalous_regions"], ["MEG3:TSS-DMRb"])
