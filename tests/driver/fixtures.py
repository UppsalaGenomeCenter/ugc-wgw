"""Small but realistic contents for the result files `ugc-wgw summary` parses; tests/driver/fake_miniwdl.py writes
them in place of its placeholder text (gzip when the file name ends in .gz). Shapes follow the real files
(checked against the smoke results of 2026-09-22)."""
from __future__ import annotations

import json

STATS_COLS = ["sample_id", "read_count", "read_length_mean", "read_length_median", "read_length_n50", "read_quality_mean",
              "read_quality_median", "mapped_read_count", "mapped_read_percent", "gap_compressed_identity_mean",
              "gap_compressed_identity_median", "depth_mean", "inferred_sex", "stat_phased_basepairs", "phase_block_ng50",
              "cpg_combined_count", "cpg_hap1_count", "cpg_hap2_count", "methbat_methylated_count",
              "methbat_unmethylated_count", "methbat_asm_count", "SNV_count", "TSTV_ratio", "HETHOM_ratio", "INDEL_count",
              "sv_DUP_count", "sv_DEL_count", "sv_INS_count", "sv_INV_count", "sv_SWAP_count", "sv_BND_count",
              "trgt_genotyped_count", "trgt_uncalled_count"]


def stats_values(sid: str, inferred_sex: str) -> dict[str, str]:
    """The stats table of a fake sample; the same numbers are exposed as `stat_*` scalar outputs."""
    return dict(zip(STATS_COLS, [sid, "3100000", "16687.5", "16052.0", "17394", "31.43", "31.0", "3099700", "99.99", "98.96", "99.7",
                                 "30.1", inferred_sex, "2900000000", "1800000", "28000000", "26000000", "26100000", "9000", "17000",
                                 "600", "3900000", "2.03", "1.55", "820000", "300", "9000", "9500", "40", "20", "700",
                                 "900000", "37000"]))


def _mosdepth(inferred_sex: str) -> str:
    x = {"MALE": 15.0, "FEMALE": 30.0}.get(inferred_sex, 0.0)
    y = {"MALE": 14.0, "FEMALE": 1.0}.get(inferred_sex, 0.0)
    auto = 30.0 if inferred_sex else 0.0
    rows = [("chr%d" % i, 200_000_000 - i * 5_000_000, auto) for i in range(1, 23)] + [("chrX", 156_040_895, x), ("chrY", 57_227_415, y), ("chrM", 16569, 3000.0)]
    lines = ["chrom\tlength\tbases\tmean\tmin\tmax"]
    for name, length, mean in rows:
        for suffix in ("", "_region"):
            lines.append(f"{name}{suffix}\t{length}\t{int(length * mean)}\t{mean:.2f}\t0\t{int(mean * 3)}")
    total = sum(r[1] for r in rows)
    lines += [f"total\t{total}\t{int(total * auto)}\t{auto:.2f}\t0\t9000", f"total_region\t{total}\t{int(total * auto)}\t{auto:.2f}\t0\t9000"]
    return "\n".join(lines) + "\n"


def _region_bed() -> str:
    lines = []
    for chrom, n, depth in (("chr1", 30, 30.0), ("chr20", 20, 31.0), ("chrX", 10, 15.0)):
        for i in range(n):
            lines.append(f"{chrom}\t{i * 500}\t{i * 500 + 500}\t{depth + (i % 3):.2f}")
    return "\n".join(lines) + "\n"


BCFTOOLS_STATS = """# This file was produced by bcftools stats (1.23.1+htslib-1.23.1) and can be plotted using plot-vcfstats.
# SN\t[2]id\t[3]key\t[4]value
SN\t0\tnumber of samples:\t1
SN\t0\tnumber of records:\t4720000
SN\t0\tnumber of no-ALTs:\t0
SN\t0\tnumber of SNPs:\t3900000
SN\t0\tnumber of MNPs:\t0
SN\t0\tnumber of indels:\t820000
SN\t0\tnumber of others:\t0
SN\t0\tnumber of multiallelic sites:\t0
SN\t0\tnumber of multiallelic SNP sites:\t0
# TSTV\t[2]id\t[3]ts\t[4]tv\t[5]ts/tv\t[6]ts (1st ALT)\t[7]tv (1st ALT)\t[8]ts/tv (1st ALT)
TSTV\t0\t2613000\t1287000\t2.03\t2613000\t1287000\t2.03
ST\t0\tA>C\t210000
ST\t0\tA>G\t650000
ST\t0\tA>T\t190000
ST\t0\tC>A\t220000
ST\t0\tC>G\t230000
ST\t0\tC>T\t660000
ST\t0\tG>A\t655000
ST\t0\tG>C\t228000
ST\t0\tG>T\t221000
ST\t0\tT>A\t189000
ST\t0\tT>C\t648000
ST\t0\tT>G\t199000
IDD\t0\t-60\t130\t130\t0.83
IDD\t0\t-2\t120000\t120000\t0.5
IDD\t0\t1\t300000\t300000\t0.5
IDD\t0\t12\t20000\t20000\t0.5
QUAL\t0\t19.5\t1000\t700\t300\t100
QUAL\t0\t45.0\t3800000\t2600000\t1200000\t800000
DP\t0\t8\t4000\t0.1\t0\t0.0
DP\t0\t30\t4000000\t84.7\t0\t0.0
DP\t0\t>500\t390\t0.01\t0\t0.0
PSC\t0\t{sid}\t0\t1850000\t2870000\t2613000\t1287000\t820000\t30.2\t150000\t0\t0\t12
"""

ROH_OUT = """# The command line was:\tbcftools roh --AF-dflt 0.4 x.vcf.gz
# RG\t[2]Sample\t[3]Chromosome\t[4]Start\t[5]End\t[6]Length (bp)\t[7]Number of markers\t[8]Quality (average fwd-bwd phred score)
# ST\t[2]Sample\t[3]Chromosome\t[4]Position\t[5]State (0:HW, 1:AZ)\t[6]Quality (fwd-bwd phred score)
ST\t{sid}\tchr1\t100\t0\t50
RG\t{sid}\tchr1\t878431\t2619312\t1740882\t531\t93.6
RG\t{sid}\tchr2\t100\t250100\t250000\t40\t80.1
RG\t{sid}\tchr7\t100\t60100\t60000\t10\t70.0
"""

SV_VCF = """##fileformat=VCFv4.2
##source="sawfish 2.2.1"
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=A,Type=Integer,Description="Length of structural variant">
##INFO=<ID=SUPP,Number=1,Type=Integer,Description="Number of samples supporting the variant">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{samples}
chr1\t100\tsv1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-250;SUPP=2\tGT\t{gts}
chr1\t900\tsv2\tN\t<INS>\t.\tPASS\tSVTYPE=INS;SVLEN=1500;SUPP=1\tGT\t{gts}
chr2\t100\tsv3\tN\tN[chr3:5[\t.\tMinQUAL\tSVTYPE=BND;SUPP=1\tGT\t{gts}
chr2\t500\tsv4\tN\t<DUP>\t.\tPASS\tSVTYPE=DUP;SVLEN=2000000;SUPP=3\tGT\t{gts}
chr20\t500\tsv5\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-25000;SUPP=3\tGT\t{gts}
"""

PHASE_STATS_HEADER = ("sample_name\tchromosome\tnum_variants\tnum_heterozygous\tnum_phased\tnum_unphased\tnum_het_snv\tnum_phased_snv\tnum_blocks\t"
                      "num_singletons\tvariants_per_block_median\tvariants_per_block_mean\tvariants_per_block_min\tvariants_per_block_max\t"
                      "variants_per_block_sum\tbasepairs_per_block_median\tbasepairs_per_block_mean\tbasepairs_per_block_min\tbasepairs_per_block_max\t"
                      "basepairs_per_block_sum\tblock_ng50")


def _phase_stats(sid: str) -> str:
    lines = [PHASE_STATS_HEADER]
    total = [0] * 6
    for i in range(1, 23):
        v, het, ph = 200000 - i * 5000, 120000 - i * 3000, 118000 - i * 3000
        blocks, bp = 300 - i * 5, 200_000_000 - i * 5_000_000
        lines.append(f"{sid}\tchr{i}\t{v}\t{het}\t{ph}\t{het - ph}\t{het - 100}\t{ph - 100}\t{blocks}\t3\t50\t400\t1\t3000\t{ph}\t100000\t300000\t1\t2000000\t{bp}\t1800000")
        for j, x in enumerate((v, het, ph, het - ph, blocks, bp)):
            total[j] += x
    lines.append(f"{sid}\tall\t{total[0]}\t{total[1]}\t{total[2]}\t{total[3]}\t{total[1] - 2200}\t{total[2] - 2200}\t{total[4]}\t66\t50\t400\t1\t3000\t{total[2]}\t100000\t300000\t1\t2000000\t{total[5]}\t1800000")
    return "\n".join(lines) + "\n"


def _phase_blocks(sid: str) -> str:
    lines = ["source_block_index\tsample_name\tphase_block_id\tchrom\tstart\tend\tnum_variants"]
    for i, (chrom, start, length) in enumerate((("chr1", 100, 5000), ("chr1", 9000, 150000), ("chr2", 300000, 2000000), ("chr20", 1000, 900000))):
        lines.append(f"{i}\t{sid}\t{start}\t{chrom}\t{start}\t{start + length}\t{10 * (i + 1)}")
    return "\n".join(lines) + "\n"


TRGT_VCF = """##fileformat=VCFv4.2
##trgtVersion=5.1.0-ec66463
##INFO=<ID=TRID,Number=1,Type=String,Description="Tandem repeat ID">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AL,Number=.,Type=Integer,Description="Length of each allele">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sid}
chr1\t16682\t.\tTGG\t.\t.\t.\tTRID=chr1_16682_16774_TGGTGGGGG;END=16774;MOTIFS=TGGTGGGGG;STRUC=(TGGTGGGGG)n\tGT:AL:ALLR:SD:MC:MS:AP:AM\t0/0:93,93:93-93,93-93:10,10:10,10:0(0-93),0(0-93):1,1:.,.
chr4\t3074876\t.\tCCAG\tCCAGCAG\t.\t.\tTRID=HD_HTT;END=3074969;MOTIFS=CAG,CCG;STRUC=<TR>\tGT:AL:ALLR:SD:MC:MS:AP:AM\t0/1:57,63:57-57,63-63:12,11:17_7,19_7:0(0-51)_1(51-72),0(0-57)_1(57-78):1,1:.,.
chr9\t69037270\t.\tCAAA\t.\t.\t.\tTRID=FRDA_FXN;END=69037304;MOTIFS=A,GAA;STRUC=<TR>\tGT:AL:ALLR:SD:MC:MS:AP:AM\t0/0:34,34:34-34,34-34:9,9:16_6,16_6:0(0-16)_1(16-34),0(0-16)_1(16-34):1,1:.,.
chr13\t102161574\t.\tGAA\t.\t.\t.\tTRID=SCA27B_FGF14;END=102161726;MOTIFS=GAA;STRUC=<TR>\tGT:AL:ALLR:SD:MC:MS:AP:AM\t.:.:.:.:.:.:.:.
chr20\t500\t.\tCAG\t.\t.\t.\tTRID=chr20_500_600_CAG;END=600;MOTIFS=CAG;STRUC=(CAG)n\tGT:AL:ALLR:SD:MC:MS:AP:AM\t0/0:100,100:100-100,100-100:8,8:33,33:0(0-100),0(0-100):1,1:.,.
"""

DROPOUTS = """chrom\tstart\tend\ttrid\texpected_ploidy\thap1_count\thap2_count\tunphased_count\tfail_read_count\tdropout
chr1\t19275\t19473\tID=chr1_19275_19473_TG;MOTIFS=TG;STRUC=(TG)n\t2\t1\t0\t0\t0\tHaplotypeDropout
chr1\t30000\t30100\tID=chr1_30000_30100_A;MOTIFS=A;STRUC=(A)n\t2\t0\t0\t0\t0\tFullDropout
chr13\t102161574\t102161726\tID=SCA27B_FGF14;MOTIFS=GAA;STRUC=<TR>\t2\t0\t0\t0\t0\tFullDropout
chr5\t1000\t1100\tID=chr5_1000_1100_AT;MOTIFS=AT;STRUC=(AT)n\t2\t3\t3\t2\t0\tPhasingDropout
"""

METHBAT = """##methbat_version=1.1.0-fbf4686
##datetime=2026-09-22 15:03:31
##command=methbat profile --input-pileup x.bed.gz --input-regions cpgIslandExt.sorted.hg38.tsv
##base_modification=5mC
##strand=combined
chrom\tstart\tend\tregion_label\tsummary_label\tcompare_label\tbackground_category\tcategory_pop_count\tcategory_pop_freq\tasm_fishers_pvalue\tmean_hap1_methyl\tmean_hap2_methyl\tmean_meth_delta\tmean_abs_meth_delta_zscore\tmean_combined_methyl\tmean_combined_methyl_delta\tmean_combined_methyl_zscore\tnum_phased_sites\tnum_partial_sites\tnum_unphased_sites\tmedian_total_coverage\tmedian_hap1_coverage\tmedian_hap2_coverage
chr1\t100\t200\tCpG:_1\tMethylated\tInsufficientData\t\t\t\t1.000000e0\t95.0\t96.0\t1.0\t\t95.5\t\t\t10\t0\t0\t30\t15\t15
chr1\t300\t400\tCpG:_2\tUnmethylated\tInsufficientData\t\t\t\t1.000000e0\t3.0\t5.0\t2.0\t\t4.0\t\t\t10\t0\t0\t30\t15\t15
chr1\t500\t600\tCpG:_3\tNoData\tInsufficientData\t\t\t\t\t\t\t\t\t\t\t\t0\t0\t0\t0\t0\t0
chr1\t700\t800\tCpG:_4\tAlleleSpecificMethylation\tInsufficientData\t\t\t\t1.000000e-9\t5.0\t95.0\t90.0\t\t50.0\t\t\t10\t0\t0\t30\t15\t15
"""

PARAPHASE = {
    "smn1": {"region_name": "smn1", "phase_region": "38:chr5:70917100-70961220", "genes_in_region": "SMN1,SMN2", "sample_sex": None,
             "genome_depth": None, "region_depth": {"median": 0.0, "percentile80": 0.0}, "failed_for_coverage": True, "total_cn": None,
             "final_haplotypes": {}, "two_copy_haplotypes": [], "region_specific_info": {}},
    "CYP2D6": {"region_name": "CYP2D6", "phase_region": "38:chr22:42126000-42145000", "genes_in_region": "CYP2D6,CYP2D7", "sample_sex": "male",
               "genome_depth": 30.1, "region_depth": {"median": 29.0, "percentile80": 33.0}, "failed_for_coverage": False, "total_cn": 4,
               "final_haplotypes": {"x": "CYP2D6_hap1", "y": "CYP2D6_hap2", "z": "CYP2D7_hap1", "w": "CYP2D7_hap2"}, "two_copy_haplotypes": [],
               "region_specific_info": {"gene_cn": 2, "fusions_called": ["CYP2D6-CYP2D7"]}},
    "RCCX": {"region_name": "rccx", "phase_region": "38:chr6:32000000-32100000", "genes_in_region": "CYP21A2,TNXB", "sample_sex": "male",
             "genome_depth": 30.1, "region_depth": {"median": 31.0, "percentile80": 35.0}, "failed_for_coverage": False, "total_cn": 4,
             "final_haplotypes": {"s1": "CYP21A2_hap1", "s2": "TNXB_hap1", "s3": "TNXB_hap2"}, "two_copy_haplotypes": ["s1"], "region_specific_info": {}},
}

MITORSAW = {"haplotypes": [{"label": "hap_0", "seq_len": 16569, "fingerprint_alleles": [], "num_ref_variants": 38, "estimated_abundance": 0.97},
                           {"label": "hap_1", "seq_len": 16569, "fingerprint_alleles": [], "num_ref_variants": 39, "estimated_abundance": 0.03}],
            "fingerprint_stats": {"fingerprint_counts": {}, "unexplained_fingerprints": [], "max_unexplained_cost": 0, "unexplained_cost": 0,
                                  "unexplained_fraction": 0.0, "passing_explanation": True}}

STARPHASE_TSV = ("#gene\tdiplotype\nABCG2\trs2231142 reference (G)/rs2231142 reference (G)\nCYP2C19\t*1/*17\nCYP2D6\t*1/*4\nMT-RNR1\tReference\n"
                 "TPMT\tNO_READS/NO_READS\nHLA-A\t*02:01:01:01/*03:01:01:01\nHLA-DRB3\t./.\n")
STARPHASE_JSON = {"pbstarphase_version": "2.2.0-5bf32fc",
                  "database_metadata": {"pbstarphase_version": "2.2.0-c56b9ba", "cpic_version": "API-2026-07-24", "hla_version": "v3.65.0-alpha",
                                        "pharmvar_version": "6.2.27", "build_time": "2026-07-24T12:44:47Z"},
                  "gene_details": {"ABCG2": {"diplotypes": [{"hap1": "ref", "hap2": "ref", "diplotype": "ref/ref"}]},
                                   "CYP2C19": {"diplotypes": [{"hap1": "*1", "hap2": "*17", "diplotype": "*1/*17"}]},
                                   "CYP2D6": {"diplotypes": [{"hap1": "*1", "hap2": "*4", "diplotype": "*1/*4"}]},
                                   "MT-RNR1": {"diplotypes": [{"hap1": "Reference", "hap2": "Reference", "diplotype": "Reference"}]},
                                   "TPMT": {"diplotypes": [{"hap1": "NO_READS", "hap2": "NO_READS", "diplotype": "NO_READS/NO_READS"}]},
                                   "HLA-A": {"diplotypes": [{"hap1": "*02:01:01:01", "hap2": "*03:01:01:01", "diplotype": "*02:01:01:01/*03:01:01:01"}],
                                             "simple_diplotypes": None, "mapping_details": {}},
                                   "HLA-DRB3": {"diplotypes": [{"hap1": ".", "hap2": ".", "diplotype": "./."}], "simple_diplotypes": None}}}
KIVVI_NA = {"allele_cn": "NA", "depth_summary": {"genome_depth": None, "repeat_depth": None}, "complete_alleles": [], "partial_alleles": [],
            "supporting_reads": {}, "allele_info": [], "complete_allele_variants": {}, "other_unit_variants": [], "methylation": {}, "read_info": {},
            "additional": {"call_status": "failed_due_to_no_reads"}}
KIVVI_OK = {"allele_cn": 38.0, "depth_summary": {"genome_depth": 30.1, "repeat_depth": 571.9}, "complete_alleles": [1, 2], "partial_alleles": [],
            "supporting_reads": {}, "allele_info": [], "complete_allele_variants": {}, "other_unit_variants": [], "methylation": {}, "read_info": {},
            "additional": {}}

COPYNUM = {"sample_name": "{sid}", "gc_bias_corrected_haploid_coverage": 15.2,
           "chromosomes": {"chr1": {"total_copy_number_bases": 230000000, "bases_per_copy_number": {"2": 229000000, "3": 1000000}, "most_common_copy_number": 2},
                           "chr20": {"total_copy_number_bases": 61639000, "bases_per_copy_number": {"2": 61000000, "1": 639000}, "most_common_copy_number": 2},
                           "chrX": {"total_copy_number_bases": 150000000, "bases_per_copy_number": {"1": 150000000}, "most_common_copy_number": 1}}}
BEDGRAPH = "chr1\t0\t100000000\t2\nchr1\t100000000\t101000000\t3\nchr1\t101000000\t230000000\t2\nchr20\t0\t61639000\t2\nchrX\t0\t150000000\t1\n"

FREQ_SUMMARY = """# ugc-pacbio-wgw cohort_freq: AF bins use the largest allele frequency of a record with AC > 0; call rate = NS / n_samples
resource\tmetric\tkey\tvalue
samples\tn\ttotal\t{n}
samples\tn\tXX\t{xx}
samples\tn\tXY\t{xy}
samples\tn\tunknown\t0
samples\tn_sex_source\tsheet\t{n}
samples\tn_sex_source\tinferred\t0
samples\tn_sex_source\tnone\t0
samples\tn_in_vcf\tsmall_variants\t{n}
samples\tn_in_vcf\tstructural_variants\t{n}
small_variants\trecords\ttotal\t5200000
small_variants\trecords_by_type\tINDEL\t900000
small_variants\trecords_by_type\tSNP\t4300000
small_variants\trecords_by_contig_class\tautosome\t5100000
small_variants\trecords_by_contig_class\tchrX\t100000
small_variants\trecords_by_filter\t.\t5200000
small_variants\trecords\tac_gt0\t5190000
small_variants\trecords\tsingletons\t1200000
small_variants\trecords_by_af_bin\t<0.001\t0
small_variants\trecords_by_af_bin\t<0.01\t0
small_variants\trecords_by_af_bin\t<0.05\t0
small_variants\trecords_by_af_bin\t>=0.05\t5190000
small_variants\trecords_ac_gt0_in\tXX\t4000000
small_variants\trecords_ac_gt0_in\tXY\t4100000
small_variants\tmean_call_rate\tall\t0.9950
small_variants\trecords\tcall_rate_below_0.9\t12000
structural_variants\trecords\ttotal\t30000
structural_variants\trecords_by_type\tDEL\t12000
structural_variants\trecords_by_type\tINS\t15000
structural_variants\trecords_by_type\tBND\t3000
structural_variants\trecords_by_contig_class\tautosome\t29000
structural_variants\trecords_by_contig_class\tchrX\t1000
structural_variants\trecords_by_filter\tPASS\t30000
structural_variants\trecords\tac_gt0\t30000
structural_variants\trecords\tsingletons\t9000
structural_variants\trecords_by_af_bin\t<0.001\t0
structural_variants\trecords_by_af_bin\t<0.01\t0
structural_variants\trecords_by_af_bin\t<0.05\t0
structural_variants\trecords_by_af_bin\t>=0.05\t30000
structural_variants\trecords_ac_gt0_in\tXX\t25000
structural_variants\trecords_ac_gt0_in\tXY\t26000
structural_variants\tmean_call_rate\tall\t1.0000
structural_variants\trecords\tcall_rate_below_0.9\t0
"""

# cohort_pangenome: the chromosomes the fake builds (chromosomes input absent), the graph sizes per chromosome, the split's
# contig binning per haplotype, the frequency summary of the vcfwave VCF and truvari's figures
PANGENOME_CHROMS = ["chr1", "chr20", "chrX"]
GRAPH_STATS = "#ref_chrom\tnodes\tedges\tlength\n{chrom}\t{nodes}\t{edges}\t{length}\n"
AMBIGUOUS_HEADER = "genome\tassigned\tassigned_bp\tambiguous\tambiguous_bp\n"
TRUVARI_JSON = {"TP-base": 24000, "TP-comp": 24500, "FP": 2500, "FN": 6000, "precision": 0.9074, "recall": 0.8, "f1": 0.8503,
                "base cnt": 30000, "comp cnt": 27000, "gt_concordance": 0.95, "TP-base_TP-gt": 22800, "TP-base_FP-gt": 1200,
                "TP-comp_TP-gt": 23275, "TP-comp_FP-gt": 1225, "gt_matrix": {"(0, 1)": {"(0, 1)": 20000, "(1, 1)": 800}},
                "params": {"sizemin": 50, "sizefilt": 50, "sizemax": 50000, "pctseq": 0.0, "pctsize": 0.7, "refdist": 500}}
PANGENOME_FREQ_SUMMARY = """# ugc-pacbio-wgw cohort_pangenome: AF bins use the largest allele frequency of a record with AC > 0; call rate = NS / n_samples
resource\tmetric\tkey\tvalue
samples\tn\ttotal\t{n}
samples\tn\tXX\t{xx}
samples\tn\tXY\t{xy}
samples\tn\tunknown\t0
samples\tn_sex_source\tkaryotype\t{n}
samples\tn_sex_source\tnone\t0
samples\tn_in_vcf\tpangenome\t{n}
pangenome\trecords\ttotal\t4100000
pangenome\trecords_by_type\tDEL\t12000
pangenome\trecords_by_type\tINDEL\t700000
pangenome\trecords_by_type\tINS\t14000
pangenome\trecords_by_type\tMIXED\t9000
pangenome\trecords_by_type\tMNV\t15000
pangenome\trecords_by_type\tSNV\t3350000
pangenome\trecords_by_contig_class\tautosome\t4000000
pangenome\trecords_by_contig_class\tchrX\t100000
pangenome\trecords_by_filter\t.\t4100000
pangenome\trecords\tac_gt0\t4100000
pangenome\trecords\tsingletons\t1500000
pangenome\trecords_by_af_bin\t<0.001\t0
pangenome\trecords_by_af_bin\t<0.01\t0
pangenome\trecords_by_af_bin\t<0.05\t0
pangenome\trecords_by_af_bin\t>=0.05\t4100000
pangenome\trecords_ac_gt0_in\tXX\t3000000
pangenome\trecords_ac_gt0_in\tXY\t3100000
pangenome\tmean_call_rate\tall\t0.9700
pangenome\trecords\tcall_rate_below_0.9\t150000
"""

# Tier 2 (0.12.0) files: per-contig QC, the callers' spectra, the graph's per-member statistics, per-chromosome path bp, degrees and gaps
CONTIG_QC = ("haplotype\tcontig\tlength\tgc_pct\tn_bp\tgaps\ttelo_5p_bp\ttelo_3p_bp\terr_bp\tdup_bp\thap_bp\tcol_bp\n"
             "{h1}\th1tg000001l\t60000000\t40.9\t0\t0\t3200\t0\t20000\t500000\t59300000\t180000\n"
             "{h1}\th1tg000002l\t5000000\t41.5\t100\t1\t0\t2800\t0\t0\t4990000\t10000\n"
             "{h2}\th2tg000001l\t61000000\t40.8\t0\t0\t0\t0\t15000\t300000\t60600000\t85000\n"
             "{h2}\th2tg000002l\t3000000\t42.0\t0\t0\t3000\t3100\t0\t0\t3000000\t0\n")
CALL_SPECTRA = ("caller\tfilter_class\ttype\tsize_bin\tcontig_class\tcount\n"
                "dipcall\tASM_UNRELIABLE\tSNV\t<50\tautosome\t20000\ndipcall\tPASS\tDEL\t100-1k\tautosome\t300\ndipcall\tPASS\tINDEL\t<50\tautosome\t600000\n"
                "dipcall\tPASS\tINS\t<100\tautosome\t400\ndipcall\tPASS\tSNV\t<50\tautosome\t3700000\ndipcall\tPASS\tSNV\t<50\tchrX\t80000\n"
                "pav\tASM_UNRELIABLE\tSNV\t<50\tautosome\t194000\npav\tPASS\tDEL\t1k-10k\tautosome\t250\npav\tPASS\tINDEL\t<50\tautosome\t700000\n"
                "pav\tPASS\tINS\t100-1k\tautosome\t900\npav\tPASS\tINV\t10k-100k\tautosome\t12\npav\tPASS\tSNV\t<50\tautosome\t3600000\n")
MEMBERS_HEADER = "sample_id\tn_ref_hom\tn_nonref_hom\tn_hets\tn_singletons\tn_missing\tn_sites\tcall_rate\n"
GENOME_BP_HEADER = "genome\tpaths\tbp\n"
DEGREE = "Degree\tSides\tNodes(min)\tNodes(max)\tNodes(total)\n0\t2\t2\t0\t0\n1\t736985\t392753\t344232\t2\n2\t376561\t169496\t207065\t344230\n3\t12480\t2094\t10386\t47132\n4\t2878\t91\t2787\t3900\n"
REFGAPS = "{chrom}\t0\t60000\trefgap\n{chrom}\t26386225\t28557053\trefgap\n"

GFATOOLS = "CC\tGS   genome_size_if_provided\nCC\tSZ   total_sequence_length\nCC\tNN   number_of_sequences\nCC\tNL   x   Nx   Lx\nCC\tAU   auN\nCC\nGS\t3100000000\nSZ\t{sz}\nNN\t{nn}\nNL\t0\t{n0}\t1\nNL\t50\t{n50}\t{l50}\nAU\t{aun}\n"
CALLER_SN = ("# bcftools stats\nSN\t0\tnumber of samples:\t1\nSN\t0\tnumber of records:\t{rec}\nSN\t0\tnumber of SNPs:\t{snps}\nSN\t0\tnumber of indels:\t{indels}\n"
             "SN\t0\tnumber of MNPs:\t0\nSN\t0\tnumber of others:\t2\nSN\t0\tnumber of multiallelic sites:\t3073\n"
             "TSTV\t0\t{ts}\t{tv}\t{tstv}\t{ts}\t{tv}\t{tstv}\n"
             "ST\t0\tA>C\t3450\nST\t0\tA>G\t14328\nST\t0\tC>T\t15000\n"
             "IDD\t0\t-60\t242\t232\t0.64\nIDD\t0\t-1\t3000\t2900\t0.5\nIDD\t0\t1\t3200\t3100\t0.5\nIDD\t0\t55\t20\t20\t0.9\n"
             "PSC\t0\t{sid}\t0\t{hom}\t{het}\t{ts}\t{tv}\t{indels}\t0.0\t{singletons}\t0\t0\t{missing}\n")
# assembly_call BEDs on chr20 (and a chr1 piece): dip.bed, reliable.bed (dip.bed minus two masked intervals), the per-haplotype dipcall
# BEDs, the lifted unreliable regions, PAV's callable regions, and hifiasm's low-quality regions in contig coordinates
DIP_BED = "chr20\t63840\t393022\nchr20\t408057\t1877050\nchr20\t2000000\t3500000\nchr1\t1000000\t1200000\n"
RELIABLE_BED = "chr20\t63840\t393022\nchr20\t408057\t1800000\nchr20\t1850000\t1877050\nchr20\t2000000\t3500000\nchr1\t1000000\t1200000\n"
HAP_BED = "chr20\t60000\t400000\nchr20\t405000\t1900000\nchr20\t1990000\t3600000\nchr1\t990000\t1210000\n"
UNRELIABLE_BED = "chr20\t1800000\t1850000\nchr20\t4000000\t4010000\n"
CALLABLE_BED = "#chrom\tpos\tend\nchr20\t60000\t1900000\nchr20\t1990000\t3600000\nchr1\t990000\t1210000\n"
LOWQ_BED = "h{hap}tg000001l\t0\t5083\t100\tread1,read2\nh{hap}tg000001l\t3914257\t3915726\t70\tread3\nh{hap}tg000002l\t10\t20\t60\tread4\n"
REGROUP_LOG = "from_haplotype\tcontig\tlength\tto_haplotype\thome_chrom\tchrY_bases\tchrX_nonpar_bases\n"
CHROM_COVERAGE = ("haplotype\tchrom\taligned_bases\tcontigs\tchrom_length\tfraction\n"
                  "{h1}\tchr1\t823\t1\t248956422\t0.0000\n{h1}\tchr20\t60000000\t3\t64444167\t0.9310\n{h1}\tchrX\t556\t1\t156040895\t0.0000\n"
                  "{h1}\tchrUn_KI270742v1\t50000\t2\t186739\t0.2677\n"
                  "{h2}\tchr20\t61000000\t2\t64444167\t0.9466\n{h2}\tchrY\t0\t0\t57227415\t0.0000\n{h2}\tchr1_KI270706v1_random\t20000\t1\t175055\t0.1142\n")
FLAGGER_SUMMARY = ("#metric\tlevel\tunit\tannotation_type\tannotation\tsize_bin\tlabel_filter\tErr\tDup\tHap\tCol\tUnk\n"
                   "PREDICTION\toverlap_based\tcount\tannotation\twhole_genome\tALL_SIZES\tALL_LABELS\t3.00\t62.00\t104.00\t8.00\t0.00\n"
                   "PREDICTION\tbase_level\tcount\tannotation\twhole_genome\tALL_SIZES\tALL_LABELS\t37000000\t41000000\t2900000000\t42000000\t0\n"
                   "PREDICTION\tbase_level\tpercentage\tannotation\twhole_genome\tALL_SIZES\tALL_LABELS\t1.24\t1.37\t95.99\t1.40\t0.00\n"
                   "PREDICTION\twindow_level\tpercentage\tannotation\twhole_genome\tALL_SIZES\tALL_LABELS\t1.30\t1.40\t95.90\t1.40\t0.00\n")
KARYOTYPE_JSON = {"karyotype": "unknown", "reason": "non-PAR chrX bases aligned are 0.000 of the chrX length, below 0.5", "trio": False,
                  "chrx_bases_outside_par": 556, "chry_bases_outside_par": 0, "chrx_fraction": 0.0, "chry_fraction": 0.0,
                  "chrx_length": 156040895, "chry_length": 57227415, "haplotypes": ["hap1", "hap2"],
                  "thresholds": {"min_mapq": 5, "chrx_min_fraction": 0.5, "chry_min_fraction": 0.1}, "regrouped": False, "moved_contigs": 0}


SOMALIER_SAMPLES_HEADER = ("#family_id\tsample_id\tpaternal_id\tmaternal_id\tsex\tphenotype\toriginal_pedigree_sex\tgt_depth_mean\tgt_depth_sd\t"
                           "depth_mean\tdepth_sd\tab_mean\tab_std\tn_hom_ref\tn_het\tn_hom_alt\tn_unknown\tp_middling_ab\tcontamination_charr\t"
                           "contamination_charr_n_sites\tX_depth_mean\tX_n\tX_hom_ref\tX_het\tX_hom_alt\tY_depth_mean\tY_n\n")
SOMALIER_PAIRS_HEADER = ("#sample_a\tsample_b\trelatedness\tibs0\tibs2\tconcordance\thets_a\thets_b\thets_ab\tshared_hets\thom_alts_a\thom_alts_b\t"
                         "shared_hom_alts\tn\tx_ibs0\tx_ibs2\texpected_relatedness\n")


def _relatedness_table() -> dict[tuple[str, str], float]:
    """UGC_WGW_FAKE_RELATEDNESS=S1:S2=0.98,S1:S3=0.5 -> {("S1","S2"): 0.98, ...}."""
    import os
    out: dict[tuple[str, str], float] = {}
    for item in os.environ.get("UGC_WGW_FAKE_RELATEDNESS", "").split(","):
        if "=" in item and ":" in item:
            pair, value = item.split("=", 1)
            a, b = pair.split(":", 1)
            try:
                out[(a, b)] = float(value)
            except ValueError:
                pass
    return out


def _declared(a: str, b: str, pedigree: dict[str, set[str]] | None) -> bool:
    return bool(pedigree) and (b in pedigree.get(a, set()) or a in pedigree.get(b, set()))  # type: ignore[union-attr]


def _qc_pairs(members: list[str], *, related_only: bool, pedigree: dict[str, set[str]] | None = None) -> str:
    """somalier's pairs table (or cohort_qc's related_pairs: the pairs at or above 0.1 or with an expectation, plus the
    `expected` column): relatedness from UGC_WGW_FAKE_RELATEDNESS, expected 0.5 for a declared parent-child pair."""
    rel = _relatedness_table()
    rows = []
    for i, a in enumerate(members):
        for b in members[i + 1:]:
            r = rel.get((a, b), rel.get((b, a), 0.01))
            declared = _declared(a, b, pedigree)
            if related_only and r < 0.1 and not declared:
                continue
            ibs0 = 0 if r >= 0.4 else 40
            expected = "0.5" if declared else "-1.0"
            rows.append(f"{a}\t{b}\t{r:.3f}\t{ibs0}\t{int(200 * r + 200)}\t0.9\t190\t190\t180\t{int(180 * r)}\t110\t110\t100\t300\t0\t0\t{expected}"
                        + (("\t" + ("0.5" if declared else "")) if related_only else "") + "\n")
    header = SOMALIER_PAIRS_HEADER.rstrip("\n") + ("\texpected\n" if related_only else "\n")
    return header + "".join(rows)


MT_HEADER = ('##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
             '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">\n##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Read depth for each allele">\n'
             '##FORMAT=<ID=VAF,Number=A,Type=Float,Description="Variant allele fractions.">\n##contig=<ID=chrM,length=16569>\n'
             '##mitorsaw_version="0.2.13-49bcf18"\n##ugc_wgw_mito_collapse_rule=every haplotype reference -> 0/0; every haplotype alternate -> a/a; mixed -> 0/a (heteroplasmic)\n'
             '##ugc_wgw_mito_haplotypes=2\n##ugc_wgw_mito_abundances=hap_0:0.9700,hap_1:0.0300\n##ugc_wgw_mito_sv_min_size=50\n')
MT_SITES = ((263, "A", "G"), (152, "T", "C"), (750, "A", "G"))   # 152 is the heteroplasmic one


def _mt_vcf(samples: list[str], het_for: set[str]) -> str:
    """A collapsed chrM VCF: 263 and 750 homoplasmic in every sample, 152 heteroplasmic (0/1, VAF 0.4) for `het_for`, 0/0 otherwise."""
    rows = []
    for pos, ref, alt in MT_SITES:
        cols = []
        for s in samples:
            if pos != 152:
                cols.append("1/1:6500:0,6490:1")
            elif s in het_for:
                cols.append("0/1:4057:2370,1672:0.4136")
            else:
                cols.append("0/0:4000:4000,0:0")
        rows.append(f"chrM\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\tGT:DP:AD:VAF\t" + "\t".join(cols))
    return MT_HEADER + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n" + "".join(r + "\n" for r in rows)


def _mt_haplogroups(haplogroup: str) -> list[str]:
    return [f"hap_0\t0.9700\t{haplogroup}\t1.0000\t1-16569\t263G 750G\t", f"hap_1\t0.0300\t{haplogroup}\t0.9500\t1-16569\t263G 750G\t152C (localPrivateMutation)"]


METHBAT_REPORT_COLS = ("chrom\tstart\tend\tregion_label\treport_summary\tqc_warnings\texpected_category\tsummary_label\tmean_combined_methyl\tmean_meth_delta\t"
                       "mean_hap1_methyl\tmean_hap2_methyl\tasm_fishers_pvalue\tnum_phased_sites\tnum_partial_sites\tnum_unphased_sites\tmedian_total_coverage\t"
                       "median_hap1_coverage\tmedian_hap2_coverage\n")


def _methbat_report(imprinting: str) -> str:
    """methbat report with three regions: two PASS, the third Inconclusive (low haplotype coverage) or, with
    imprinting = "anomalous", Anomalous without a QC warning."""
    third = ("Anomalous", "PASS", "Methylated", "92.0", "4.0") if imprinting == "anomalous" else ("Inconclusive", "LowHaplotypeCoverage", "NoData", "", "")
    return ("##methbat_version=1.1.0-fbf4686\n" + METHBAT_REPORT_COLS
            + "chr6\t144006940\t144008751\tPLAGL1:alt-TSS-DMR\tPASS\tPASS\tAlleleSpecificMethylation\tAlleleSpecificMethylation\t46.0\t83.5\t13.2\t96.7\t0.0\t143\t0\t0\t41\t25\t16\n"
            + "chr20\t58853849\t58856408\tGNAS-XL:Ex1-DMR\tPASS\tPASS\tAlleleSpecificMethylation\tAlleleSpecificMethylation\t54.0\t50.7\t28.5\t79.2\t0.0\t201\t0\t0\t34\t17\t17\n"
            + f"chr14\t100824186\t100827641\tMEG3:TSS-DMRb\t{third[0]}\t{third[1]}\tAlleleSpecificMethylation\t{third[2]}\t{third[3]}\t{third[4]}\t\t\t1.000000e0\t0\t0\t0\t0\t0\t0\n")


def _methyl_samples(members: list[str], legacy: list[str], labels: list[str]) -> str:
    rows = []
    for i, m in enumerate(members):
        lab = labels[i] if i < len(labels) else ""
        n_hypo = 0 if i == 0 else 1
        anomalous = 1 if m in legacy else 0
        rows.append(f"{m}\t{lab}\t{n_hypo}\t0\t0\t0\t{n_hypo}\t{3 - n_hypo}\t3\t1\t1\t1\t0\t0\t2\t{1 - anomalous}\t{anomalous}\t0\t{'stage' if m in legacy else 'run'}\t\t\n")
    return ("sample_id\tlabels\tn_hypo\tn_hyper\tn_hypo_asm\tn_hyper_asm\tn_events\tn_insufficient\tn_regions_scored\tn_methylated\tn_unmethylated\tn_asm\t"
            "n_uncategorized\tn_nodata\timprinting_pass\timprinting_inconclusive\timprinting_anomalous\timprinting_anomalous_qc\treport_source\t"
            "deconvolution_top_cell_type\tdeconvolution_top_fraction\n" + "".join(rows))


def _methyl_summary(members: list[str], legacy: list[str], labels: list[str]) -> str:
    rows = ["resource\tmetric\tkey\tvalue"]
    n_events = max(0, len(members) - 1)
    for i, m in enumerate(members):
        n = 0 if i == 0 else 1
        rows += [f"members\tn_events\t{m}\t{n}", f"members\tn_hypo\t{m}\t{n}", f"members\tn_regions_scored\t{m}\t3",
                 f"imprinting\tPASS\t{m}\t2", f"imprinting\tInconclusive\t{m}\t{0 if m in legacy else 1}", f"imprinting\tAnomalous\t{m}\t{1 if m in legacy else 0}",
                 f"imprinting\tAnomalousQcWarning\t{m}\t0"]
    rows += [f"cohort\tn_members\t\t{len(members)}", "cohort\tn_regions\t\t4", f"cohort\tn_events\t\t{n_events}", "cohort\tbackground_min_datasets\t\t10",
             f"cohort\treports_made_in_stage\t\t{sum(1 for m in members if m in legacy)}", f"labels\tn_members\tALL\t{len(members)}"]
    counts: dict[str, int] = {}
    for lab in labels:
        for part in lab.split(";"):
            if part:
                counts[part] = counts.get(part, 0) + 1
    for lab, n in sorted(counts.items()):
        rows.append(f"labels\tn_members\t{lab}\t{n}")
        if n < 10:
            rows.append(f"notes\tlabel_below_min\t{lab}\t{n} member(s): comparisons against {lab} are InsufficientData")
    if len(members) < 10:
        rows.append(f"notes\tlabel_below_min\tALL\t{len(members)} member(s): comparisons against ALL are InsufficientData")
    if n_events:
        rows += [f"genes\tmembers_with_event\tGENEA\t{n_events}", f"genes\tmembers_with_event\tGENEB\t{n_events}", f"region_type\tn_events\tpromoter\t{n_events}"]
    rows += ["compare\tXX:XY\tInsufficientData\t4", "notes\tcompare\tXX:XY\tcompare XX:XY: 1 vs 1 members, below min_samples 10: every row is InsufficientData"]
    return "\n".join(rows) + "\n"


def fixture_bytes(name: str, sid: str, stage: str, inferred_sex: str = "", index: int = 0,
                  members: list[str] | None = None, pedigree: dict[str, set[str]] | None = None,
                  sex: dict[str, str] | None = None, legacy: list[str] | None = None, mt: bool = False,
                  haplogroup: str = "H5a7", imprinting: str = "", labels: list[str] | None = None) -> bytes | None:
    """The content of output `name` for sample/cohort `sid`, or None for outputs that stay placeholders. `pedigree`
    (child -> parents in the cohort) shapes cohort_qc's expectations; `sex` (declared, coverage, somalier, resolved,
    source) the sex resolution table of a per-sample run."""
    text: str | None = None
    doc: object = None
    if name == "stats_file":
        vals = stats_values(sid, inferred_sex)
        text = "\t".join(STATS_COLS) + "\n" + "\t".join(vals[c] for c in STATS_COLS) + "\n"
    elif name == "msg_file":
        text = f"INCLUDE_FAIL_READS regions: HD_HTT,FRDA_FXN\nInput {sid}.hifi_reads.bam is already aligned.  Alignments and haplotype tags will be stripped, and chunking will be disabled.\n"
    elif name == "mosdepth_summary":
        text = _mosdepth(inferred_sex)
    elif name == "mosdepth_region_bed":
        text = _region_bed()
    elif name == "sv_copynum_summary":
        doc = json.loads(json.dumps(COPYNUM).replace("{sid}", sid))
    elif name == "sv_copynum_bedgraph":
        text = BEDGRAPH
    elif name == "small_variant_stats":
        text = BCFTOOLS_STATS.replace("{sid}", sid)
    elif name == "bcftools_roh_out":
        text = ROH_OUT.replace("{sid}", sid)
    elif name == "mt_small_variant_vcf":
        text = _mt_vcf([sid], {sid})
    elif name == "mt_sv_vcf":
        text = _mt_vcf([sid], set()).split("#CHROM")[0] + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sid + "\n"
    elif name == "mt_haplotypes_vcf":
        text = (MT_HEADER + f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sid}.hap_0\t{sid}.hap_1\n"
                + "".join(f"chrM\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\tGT\t{'0' if pos == 152 else '1'}\t1\n" for pos, ref, alt in MT_SITES))
    elif name == "mt_haplotypes_tsv":
        text = "haplotype\tabundance\tn_variants\nhap_0\t0.9700\t2\nhap_1\t0.0300\t3\n"
    elif name == "mt_summary":
        text = ("sample_id\tvcf_sample\thaplotypes\trecords\thomoplasmic\theteroplasmic\tsvs\thet_vaf_mean\thet_vaf_max\tmajor_haplotype\tmajor_abundance\n"
                f"{sid}\t{sid}\t2\t3\t2\t1\t0\t0.4136\t0.4136\thap_0\t0.9700\n")
    elif name == "mt_haplogroups":
        text = "haplotype\tabundance\thaplogroup\tquality\trange\tfound_polys\tremaining_polys\n" + "".join(r + "\n" for r in _mt_haplogroups(haplogroup))
    elif name == "cohort_mt_small_variant_vcf":
        ms = members or []
        text = _mt_vcf(ms, {m for i, m in enumerate(ms) if i % 2 == 0})
    elif name == "cohort_mt_sv_vcf":
        ms = members or []
        text = _mt_vcf(ms, set()).split("#CHROM")[0] + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(ms) + "\n"
    elif name == "cohort_mt_haplogroups":
        ms = members or []
        text = ("sample_id\thaplotype\tabundance\thaplogroup\tquality\trange\tfound_polys\tremaining_polys\tmade_in_stage\n"
                + "".join(f"{m}\t{r}\t{'yes' if m in (legacy or []) else 'no'}\n" for m in ms for r in _mt_haplogroups("K1a1b1a" if m in (legacy or []) else haplogroup)))
    elif name in ("phased_sv_vcf", "sv_vcf"):
        text = SV_VCF.replace("{samples}", sid).replace("{gts}", "0/1")
    elif name == "cohort_sv_vcf":
        ms = members or []
        text = SV_VCF.replace("{samples}", "\t".join(ms)).replace("{gts}", "\t".join(["0/1"] * len(ms)))
    elif name == "cohort_small_variant_vcf":
        ms = members or []
        text = ("##fileformat=VCFv4.2\n##source=GLnexus\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(ms) + "\n"
                + "chr1\t100\t.\tA\tG\t50\t.\tAF=0.5\tGT\t" + "\t".join(["0/1"] * len(ms)) + "\n")
    elif name == "phase_stats":
        text = _phase_stats(sid)
    elif name == "phase_blocks":
        text = _phase_blocks(sid)
    elif name in ("phased_trgt_vcf", "trgt_vcf"):
        text = TRGT_VCF.replace("{sid}", sid)
    elif name == "trgt_coverage_dropouts":
        text = DROPOUTS
    elif name == "methbat_profile":
        text = METHBAT
    elif name in ("paraphase_summary", "paraphase_output_json"):
        doc = PARAPHASE
    elif name == "mitorsaw_hap_stats":
        doc = MITORSAW
    elif name == "pbstarphase_tsv":
        text = STARPHASE_TSV
    elif name in ("pbstarphase_summary", "pbstarphase_json"):
        doc = STARPHASE_JSON
    elif name == "kivvi_kiv2_json":
        doc = KIVVI_NA
    elif name == "kivvi_d4z4_json":
        doc = KIVVI_OK
    elif name == "freq_summary":
        ms = members or []
        tpl = PANGENOME_FREQ_SUMMARY if stage == "cohort_pangenome" else FREQ_SUMMARY
        text = tpl.replace("{n}", str(len(ms))).replace("{xx}", str(len(ms) // 2)).replace("{xy}", str(len(ms) - len(ms) // 2))
        if mt and stage == "cohort_freq":
            text += f"samples\tn_in_vcf\tmt_small_variants\t{len(ms)}\n" + "".join(
                f"mt_small_variants\t{m}\t{k}\t{val}\n" for m, k, val in (("records", "total", 3), ("records_by_type", "SNP", 3), ("records_by_contig_class", "chrM", 3),
                                                               ("records_by_filter", "PASS", 3), ("records", "ac_gt0", 3), ("records", "singletons", 0),
                                                               ("records_by_af_bin", "<0.001", 0), ("records_by_af_bin", "<0.01", 0), ("records_by_af_bin", "<0.05", 0),
                                                               ("records_by_af_bin", ">=0.05", 3), ("records_ac_gt0_in", "XX", 0), ("records_ac_gt0_in", "XY", 0),
                                                               ("mean_call_rate", "all", "1.0000"), ("records", "call_rate_below_0.9", 0)))
            text += "mt_svs\trecords\ttotal\t0\nmt_svs\trecords\tac_gt0\t0\nmt_svs\trecords\tsingletons\t0\nmt_svs\tmean_call_rate\tall\tNA\nmt_svs\trecords\tcall_rate_below_0.9\t0\n"
    elif name == "freq_samples" and stage == "cohort_pangenome":
        ms = members or []
        text = "sample_id\tgroup\tsex_source\tin_pangenome_vcf\n" + "".join(
            f"{m}\t{'XX' if i % 2 else 'XY'}\tkaryotype\tyes\n" for i, m in enumerate(ms))
    elif name == "freq_samples":
        ms = members or []
        text = "sample_id\tgroup\tsex_source\tin_small_variant_vcf\tin_sv_vcf\tin_mt_vcf\n" + "".join(
            f"{m}\t{'XX' if i % 2 else 'XY'}\tdeclared\tyes\tyes\t{'yes' if mt else 'none'}\n" for i, m in enumerate(ms))
    elif name == "somalier_sex_tsv":
        sx = sex or {"declared": inferred_sex, "coverage": inferred_sex, "somalier": "", "resolved": inferred_sex, "source": "declared" if inferred_sex else "none"}
        som = sx.get("somalier") or ""
        x_n, x_het, x_hom = (300, 5, 150) if som.upper() == "MALE" else (300, 120, 80) if som.upper() == "FEMALE" else (0, 0, 0)
        ratio = x_het / x_hom if x_hom else 0.0
        reason = ("somalier: X het/hom-alt ratio {:.4f} at {} sites".format(ratio, x_n) if som else "somalier: 0 chrX sites genotyped, below 50")
        text = ("sample_id\tdeclared\tcoverage\tsomalier\tresolved\tsource\tx_n\tx_het\tx_hom_alt\tx_het_ratio\tx_depth_mean\ty_depth_norm\treason\n"
                f"{sid}\t{sx.get('declared', '')}\t{sx.get('coverage', '')}\t{som}\t{sx.get('resolved', '')}\t{sx.get('source', '')}\t{x_n}\t{x_het}\t{x_hom}\t{ratio:.4f}\t"
                f"{15.0 if x_n else 0.0:.2f}\t{0.5 if som.upper() == 'MALE' else 0.0:.3f}\t{reason}\n")
    elif name in ("qc_pairs", "qc_related_pairs"):
        text = _qc_pairs(members or [], related_only=name == "qc_related_pairs", pedigree=pedigree)
    elif name == "qc_samples":
        ms = members or []
        text = SOMALIER_SAMPLES_HEADER + "".join(
            f"C\t{m}\t-9\t-9\t{1 if i % 2 == 0 else 2}\t-9\t-9\t30.0\t5.0\t1.0\t5.0\t0.5\t0.4\t170\t{190 + i}\t{110 + i}\t16900\t0.000\t0.01\t300\t15.0\t300\t100\t{5 if i % 2 == 0 else 120}\t{150 if i % 2 == 0 else 80}\t{8.0 if i % 2 == 0 else 0.0}\t10\n"
            for i, m in enumerate(ms))
    elif name == "qc_groups":
        text = ""
    elif name == "qc_ped":
        ms = members or []
        text = "".join(f"C\t{m}\t-9\t-9\t{1 if i % 2 == 0 else 2}\t-9\n" for i, m in enumerate(ms))
    elif name == "qc_members":
        ms = members or []
        rel = _relatedness_table()
        rows = []
        for i, m in enumerate(ms):
            best = max(((o, rel.get((m, o), rel.get((o, m), 0.01))) for o in ms if o != m), key=lambda t: t[1], default=("", 0.0))
            expected = "0.50" if best[0] and _declared(m, best[0], pedigree) else ""
            rows.append(f"{m}\t{'XY' if i % 2 == 0 else 'XX'}\tdeclared\t{'MALE' if i % 2 == 0 else 'FEMALE'}\t300\t{0.03 if i % 2 == 0 else 1.5}\t"
                        f"{0.63 + 0.01 * i:.4f}\t30.0\t{best[0]}\t{best[1]:.3f}\t{40 if best[1] < 0.2 else 0}\t{expected}\t{'EUR' if i % 2 == 0 else 'AFR'}\t0.9\tno\n")
        text = ("sample_id\tsex_group\tsex_source\tsomalier_sex\tx_n\tx_het_ratio\thet_ratio\tdepth_mean\tclosest_relative\tclosest_relatedness\t"
                "closest_ibs0\tclosest_expected\tancestry\tancestry_prob\textracted_here\n" + "".join(rows))
    elif name == "qc_summary":
        ms = members or []
        n = len(ms) * (len(ms) - 1) // 2
        text = f"section\tmetric\tkey\tvalue\npairs\tn\ttotal\t{n}\npairs\tn\tkept\t0\npairs\thist\t0.00\t{n}\nmembers\tn\ttotal\t{len(ms)}\nmembers\tn\textracted_here\t0\nmembers\tn\twith_ancestry\t{len(ms)}\n"
    elif name == "methbat_report" or name == "member_reports":
        text = _methbat_report(imprinting if name == "methbat_report" else ("anomalous" if (members or []) and index < len(members or []) and (members or [])[index] in (legacy or []) else ""))
    elif name == "methyl_samples":
        text = _methyl_samples(members or [], legacy or [], labels or [])
    elif name == "methyl_summary":
        text = _methyl_summary(members or [], legacy or [], labels or [])
    elif name == "methyl_regions":
        text = ("chrom\tstart\tend\tregion_label\tlabel\tn_datasets\tn_nodata\tn_uncategorized\tn_methylated\tn_unmethylated\tn_asm\tavg_combined_methyl\t"
                "stdev_combined_methyl\tavg_abs_meth_delta\tstdev_abs_meth_delta\tn_hypo\tn_hyper\tn_hypo_asm\tn_hyper_asm\thgnc_id\tgene_symbol\tregion_type\n"
                + "".join(f"chr1\t{100 + 200 * i}\t{200 + 200 * i}\tCpG:_{i + 1}\t{lab}\t{len(members or [])}\t0\t0\t{len(members or [])}\t0\t0\t80.0\t4.0\t1.5\t0.5\t"
                          f"{(len(members or []) - 1 if i == 0 and lab == 'ALL' else 0)}\t0\t0\t0\t{'HGNC:1;HGNC:2' if i == 0 else ''}\t{'GENEA;GENEB' if i == 0 else ''}\t{'promoter' if i == 0 else 'intergenic'}\n"
                          for i in range(4) for lab in ("ALL",)))
    elif name == "methyl_events":
        ms = members or []
        text = ("sample_id\tchrom\tstart\tend\tregion_label\tsummary_label\tcompare_label\tbackground_category\tcategory_pop_count\tcategory_pop_freq\t"
                "asm_fishers_pvalue\tmean_hap1_methyl\tmean_hap2_methyl\tmean_meth_delta\tmean_abs_meth_delta_zscore\tmean_combined_methyl\tmean_combined_methyl_delta\t"
                "mean_combined_methyl_zscore\tnum_phased_sites\tnum_partial_sites\tnum_unphased_sites\tmedian_total_coverage\tmedian_hap1_coverage\tmedian_hap2_coverage\t"
                "hgnc_id\tgene_symbol\tregion_type\n"
                + "".join(f"{m}\tchr1\t100\t200\tCpG:_1\tUnmethylated\tHypoMethylated\tALL\t{len(ms)}\t0.9\t0.5\t10\t12\t2\t0.1\t11\t-60\t-4.2\t20\t0\t0\t30\t15\t15\t"
                          f"HGNC:1;HGNC:2\tGENEA;GENEB\tpromoter\n" for m in ms[1:]))
    elif name == "methyl_imprinting":
        ms = members or []
        rows = []
        for m in ms:
            for region, summ, warn in (("PLAGL1:alt-TSS-DMR", "PASS", "PASS"), ("GNAS-XL:Ex1-DMR", "PASS", "PASS"),
                                       ("MEG3:TSS-DMRb", "Anomalous" if m in (legacy or []) else "Inconclusive", "PASS" if m in (legacy or []) else "LowHaplotypeCoverage")):
                rows.append(f"{m}\tchr6\t1\t2\t{region}\t{summ}\t{warn}\tAlleleSpecificMethylation\tAlleleSpecificMethylation\t50.0\t70.0\t15.0\t85.0\t0.0\t100\t0\t0\t30\t15\t15\t"
                            f"{'stage' if m in (legacy or []) else 'run'}\n")
        text = ("sample_id\tchrom\tstart\tend\tregion_label\treport_summary\tqc_warnings\texpected_category\tsummary_label\tmean_combined_methyl\tmean_meth_delta\t"
                "mean_hap1_methyl\tmean_hap2_methyl\tasm_fishers_pvalue\tnum_phased_sites\tnum_partial_sites\tnum_unphased_sites\tmedian_total_coverage\t"
                "median_hap1_coverage\tmedian_hap2_coverage\treport_source\n" + "".join(rows))
    elif name in ("cohort_background", "methyl_collection", "scored_profiles", "scored_profiles_annotated", "comparisons"):
        text = None
    elif name == "input_checks":
        text = ("check\ttarget\tstatus\tdetail\nbam_index\tbam\tpass\tx.bam.bai\nbam_sorted\tbam\tpass\t@HD SO:coordinate\n"
                "bam_reference\tbam\tpass\t25 contigs match the reference\nbam_methylation_tags\tbam\tpass\t1990 of 2000 first reads carry MM/ML\n")
    elif name == "input_provenance":
        text = "source\tkey\tvalue\nbam\tPG.pbmm2.PN\tpbmm2\nbam\tPG.pbmm2.VN\t1.13.1\nbam\tRG.rg1.SM\t" + sid + "\n"
    elif name == "qc_ancestry":
        ms = members or []
        text = ("#sample_id\tpredicted_ancestry\tgiven_ancestry\tEAS_prob\tAFR_prob\tAMR_prob\tSAS_prob\tEUR_prob\tPC1\tPC2\tPC3\tPC4\tPC5\n"
                "HG00096\tEUR\tEUR\t0.0\t0.0\t0.1\t0.0\t0.9\t0.5\t2.0\t0.1\t0.2\t0.3\nHG00403\tEAS\tEAS\t1.0\t0.0\t0.0\t0.0\t0.0\t2.6\t-7.1\t0.1\t0.2\t0.3\n"
                + "".join(f"{m}\t{'EUR' if i % 2 == 0 else 'AFR'}\t\t0.0\t{0.0 if i % 2 == 0 else 0.9}\t0.1\t0.0\t{0.9 if i % 2 == 0 else 0.0}\t{0.4 + i * 0.1:.1f}\t{2.0 - i:.1f}\t0.1\t0.2\t0.3\n"
                          for i, m in enumerate(ms)))
    elif name == "chrom_graph_stats":
        chrom = PANGENOME_CHROMS[index] if index < len(PANGENOME_CHROMS) else f"chr{index}"
        text = GRAPH_STATS.format(chrom=chrom, nodes=500_000 + index * 1000, edges=700_000 + index * 1000, length=60_000_000 + index * 1_000_000)
    elif name == "ambiguous_contigs":
        ms = members or []
        text = AMBIGUOUS_HEADER + "".join(f"{m}.{h}\t{30 + i}\t{3_000_000_000 - i * 1_000_000}\t{18 + i}\t{2_000_000 + i * 500_000}\n"
                                          for i, m in enumerate(ms) for h in (1, 2))
    elif name == "contig_sizes":
        ms = members or []
        gens = [f"{m}.{h}" for m in ms for h in (1, 2)]
        text = ("Contig\tGRCh38\t" + "\t".join(gens) + "\t_MINIGRAPH_\tmin\tmax\tavg\n"
                + "_AMBIGUOUS_\t0\t" + "\t".join(["2000000"] * len(gens)) + "\t0\t0\t2000000\t1000000\n"
                + "chr20\t64444167\t" + "\t".join(["65000000"] * len(gens)) + "\t64444167\t64444167\t65000000\t64800000\n")
    elif name == "name_map":
        ms = members or []
        text = "sample_id\tsafe_name\ttrio\n" + "".join(f"{m}\t{m.replace('-', '_').replace('.', '_')}\tfalse\n" for m in ms)
    elif name == "truvari_summary":
        doc = TRUVARI_JSON
    elif name == "contig_qc":
        text = CONTIG_QC.format(h1="hap1", h2="hap2")
    elif name == "call_spectra":
        text = CALL_SPECTRA
    elif name == "member_stats":
        ms = members or []
        text = MEMBERS_HEADER + "".join(f"{m}\t{1000000 + i}\t{800000 - i * 1000}\t{1500000 + i * 2000}\t{30000 + i * 5000}\t{90000 + i * 20000}\t4100000\t{1 - (90000 + i * 20000) / 4100000:.4f}\n"
                                        for i, m in enumerate(ms))
    elif name == "member_stats_shards":
        text = "# bcftools stats\nSN\t0\tnumber of records:\t1366666\n"
    elif name == "chrom_genome_bp":
        ms = members or []
        chrom = PANGENOME_CHROMS[index] if index < len(PANGENOME_CHROMS) else f"chr{index}"
        text = GENOME_BP_HEADER + f"GRCh38#0\t1\t{60_000_000 + index * 1_000_000}\n" + "".join(f"{m}#{h}\t{3 + i}\t{61_000_000 + i * 100_000 - h * 50_000}\n" for i, m in enumerate(ms) for h in (1, 2))
    elif name == "chrom_degree_stats":
        text = DEGREE
    elif name == "chrom_refgaps":
        chrom = PANGENOME_CHROMS[index] if index < len(PANGENOME_CHROMS) else f"chr{index}"
        text = REFGAPS.format(chrom=chrom)
    elif name == "cohort_trgt_lps":
        ms = members or []
        text = "trid\tmotif\t" + "\t".join(ms) + "\n" + "".join(f"chr1_{i}_{i + 50}_CAG\tCAG\t" + "\t".join(["0,0"] * len(ms)) + "\n" for i in range(1, 6))
    elif name == "assembly_stats":
        text = GFATOOLS.format(sz=3_020_000_000 - index * 5_000_000, nn=400 + index * 20, n0=140_000_000, n50=60_000_000 - index * 3_000_000, l50=15 + index, aun=58_000_000)
    elif name == "dipcall_vcf_stats":
        text = CALLER_SN.format(rec=4_500_000, snps=3_800_000, indels=700_000, ts=2_500_000, tv=1_300_000, tstv=1.92, sid=sid, hom=1_500_000, het=2_900_000, singletons=4_300_000, missing=5000)
    elif name == "pav_vcf_stats":
        text = CALLER_SN.format(rec=4_600_000, snps=3_850_000, indels=750_000, ts=1_900_000, tv=1_950_000, tstv=0.97, sid=sid, hom=1_400_000, het=3_100_000, singletons=4_400_000, missing=120_000)
    elif name == "dipcall_bed":
        text = DIP_BED
    elif name == "reliable_bed":
        text = RELIABLE_BED
    elif name == "dipcall_hap_beds":
        text = HAP_BED
    elif name == "unreliable_beds":
        text = UNRELIABLE_BED if index == 0 else "chr20\t4000000\t4010000\n"
    elif name == "pav_callable_beds":
        text = CALLABLE_BED
    elif name == "assembly_lowQ_beds":
        text = LOWQ_BED.format(hap=index + 1)
    elif name == "regroup_log":
        text = REGROUP_LOG
    elif name == "assembler_version":
        text = "0.25.0-r726\n"
    elif name == "karyotype_json":
        doc = KARYOTYPE_JSON
    elif name == "flagger_summary":
        text = FLAGGER_SUMMARY
    elif name == "trio_eval":
        text = f"S\th{index + 1}tg000001l\t16624\t22\t16610\t13\t13\t9\nW\t{199 + index}\t49822\t{0.0040 + index * 0.001:.6f}\nH\t540\t49871\t0.010828\n"
    elif name == "reliable_summary":
        text = "reference_bp\t3031042417\ndip_bed_bp\t2800000000\nunreliable_bp\t40000000\nunreliable_in_dip_bed_bp\t30000000\nreliable_bp\t2770000000\nreliable_fraction_of_reference\t0.9139\nreliable_fraction_of_dip_bed\t0.9893\n"
    elif name == "concordance":
        text = "region\treliable.bed\nsites_merged\t3900000\nboth_called\t3700000\nconcordant\t3650000\ndiscordant\t50000\nread_het\t2200000\nhet_loss\t44000\nread_only\t120000\nasm_only\t80000\nasm_only_het\t30000\ngenotype_concordance\t0.9865\nhet_loss_rate\t0.0200\n"
    elif name == "chrom_coverage":
        text = CHROM_COVERAGE.format(h1="hap1", h2="hap2")
    if doc is not None:
        text = json.dumps(doc, indent=2) + "\n"
    return text.encode() if text is not None else None
