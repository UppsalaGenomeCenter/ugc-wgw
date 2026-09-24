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
               "final_haplotypes": {"CYP2D6_hap1": "x", "CYP2D6_hap2": "y", "CYP2D7_hap1": "z", "CYP2D7_hap2": "w"}, "two_copy_haplotypes": [], "region_specific_info": {}},
    "RCCX": {"region_name": "rccx", "phase_region": "38:chr6:32000000-32100000", "genes_in_region": "CYP21A2,TNXB", "sample_sex": "male",
             "genome_depth": 30.1, "region_depth": {"median": 31.0, "percentile80": 35.0}, "failed_for_coverage": False, "total_cn": 4,
             "final_haplotypes": {"a": 1, "b": 2, "c": 3, "d": 4}, "two_copy_haplotypes": [], "region_specific_info": {}},
}

MITORSAW = {"haplotypes": [{"label": "hap_0", "seq_len": 16569, "fingerprint_alleles": [], "num_ref_variants": 38, "estimated_abundance": 0.97},
                           {"label": "hap_1", "seq_len": 16569, "fingerprint_alleles": [], "num_ref_variants": 39, "estimated_abundance": 0.03}],
            "fingerprint_stats": {"fingerprint_counts": {}, "unexplained_fingerprints": [], "max_unexplained_cost": 0, "unexplained_cost": 0,
                                  "unexplained_fraction": 0.0, "passing_explanation": True}}

STARPHASE_TSV = "#gene\tdiplotype\nABCG2\trs2231142 reference (G)/rs2231142 reference (G)\nCYP2C19\t*1/*17\nCYP2D6\t*1/*4\nMT-RNR1\tReference\n"
STARPHASE_JSON = {"pbstarphase_version": "2.2.0-5bf32fc",
                  "database_metadata": {"pbstarphase_version": "2.2.0-c56b9ba", "cpic_version": "API-2026-07-24", "hla_version": "v3.65.0-alpha",
                                        "pharmvar_version": "6.2.27", "build_time": "2026-07-24T12:44:47Z"},
                  "gene_details": {"ABCG2": {"diplotypes": [{"hap1": "ref", "hap2": "ref", "diplotype": "ref/ref"}]},
                                   "CYP2C19": {"diplotypes": [{"hap1": "*1", "hap2": "*17", "diplotype": "*1/*17"}]},
                                   "CYP2D6": {"diplotypes": [{"hap1": "*1", "hap2": "*4", "diplotype": "*1/*4"}]},
                                   "MT-RNR1": {"diplotypes": [{"hap1": "Reference", "hap2": "Reference", "diplotype": "Reference"}]}}}
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

GFATOOLS = "CC\tGS   genome_size_if_provided\nCC\tSZ   total_sequence_length\nCC\tNN   number_of_sequences\nCC\tNL   x   Nx   Lx\nCC\tAU   auN\nCC\nGS\t3100000000\nSZ\t{sz}\nNN\t{nn}\nNL\t0\t{n0}\t1\nNL\t50\t{n50}\t{l50}\nAU\t{aun}\n"
PAFTOOLS_SN = "# bcftools stats\nSN\t0\tnumber of samples:\t1\nSN\t0\tnumber of records:\t{rec}\nSN\t0\tnumber of SNPs:\t{snps}\nSN\t0\tnumber of indels:\t{indels}\n"


def fixture_bytes(name: str, sid: str, stage: str, inferred_sex: str = "", index: int = 0,
                  members: list[str] | None = None) -> bytes | None:
    """The content of output `name` for sample/cohort `sid`, or None for outputs that stay placeholders."""
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
        text = FREQ_SUMMARY.replace("{n}", str(len(ms))).replace("{xx}", str(len(ms) // 2)).replace("{xy}", str(len(ms) - len(ms) // 2))
    elif name == "freq_samples":
        ms = members or []
        text = "sample_id\tgroup\tsex_source\tin_small_variant_vcf\tin_sv_vcf\n" + "".join(
            f"{m}\t{'XX' if i % 2 else 'XY'}\tsheet\tyes\tyes\n" for i, m in enumerate(ms))
    elif name == "cohort_trgt_lps":
        ms = members or []
        text = "trid\tmotif\t" + "\t".join(ms) + "\n" + "".join(f"chr1_{i}_{i + 50}_CAG\tCAG\t" + "\t".join(["0,0"] * len(ms)) + "\n" for i in range(1, 6))
    elif name == "assembly_stats":
        text = GFATOOLS.format(sz=3_020_000_000 - index * 5_000_000, nn=400 + index * 20, n0=140_000_000, n50=60_000_000 - index * 3_000_000, l50=15 + index, aun=58_000_000)
    elif name == "paftools_vcf_stats":
        text = PAFTOOLS_SN.format(rec=3_000_000 + index, snps=2_700_000 + index, indels=300_000)
    if doc is not None:
        text = json.dumps(doc, indent=2) + "\n"
    return text.encode() if text is not None else None
