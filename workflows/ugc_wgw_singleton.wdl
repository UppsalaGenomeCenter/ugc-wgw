version 1.0

# ugc-pacbio-wgw — stage `singleton`: humanwgs_singleton call for call, reference data from the installed map instead of the data container, plus ugc_wgw_manifest.json
# origin: none (mirrors vendor/hifi-human-wgs-wdl/workflows/singleton.wdl @ v4.0.0; re-derive on every sync, docs/ENTRYPOINTS.md §3)
# see docs/DESIGN.md §5, docs/ENTRYPOINTS.md §3

import "../vendor/hifi-human-wgs-wdl/workflows/downstream/downstream.wdl" as Downstream
import "../vendor/hifi-human-wgs-wdl/workflows/process_trgt_catalog/process_trgt_catalog.wdl" as ProcessTrgtCatalog
import "overrides/upstream.wdl" as Upstream
import "../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/tasks/utilities.wdl" as Utilities
import "../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/workflows/backend_configuration/backend_configuration.wdl" as BackendConfiguration
import "../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/workflows/filter_messages/filter_messages.wdl" as FilterMessages
import "ugc_wgw/provenance.wdl" as Provenance
import "ugc_wgw/methyl/report.wdl" as MethylReport

workflow ugc_wgw_singleton {
  meta {
    description: "ugc-pacbio-wgw stage singleton: the PacBio HiFi human WGS singleton pipeline (humanwgs_singleton 4.0.0) with unchanged semantics, its reference data read from the installed reference map rather than unpacked from the workflow data container in every run, plus the ugc provenance record"
  }

  parameter_meta {
    ugc_wgw_version: {
      description: "ugc-pacbio-wgw version the run was submitted with; echoed as ugc_wgw_workflow_version"
    }
    sample_id: {
      description: "Sample ID"
    }
    hifi_reads: {
      description: "Unaligned or aligned hifi_reads BAMs (aligned input is realigned as a single chunk, with a message)"
    }
    fail_reads: {
      description: "Unaligned fail_reads BAMs, used only for TRGT at INCLUDE_FAIL_READS loci of the catalog"
    }
    ref_map_file: {
      description: "ugc-wgw reference map (docs/DESIGN.md §10.4): the files and scalars of upstream's workflow data container for one reference build, plus scatter_regions and pbmm2_index (ugc-wgw-extras)"
    }
    use_alignment_chunking: {
      description: "Whether to chunk BAM files for alignment. If false, all reads will be aligned in a single chunk."
    }
    use_gpu: {
      description: "Use GPU when possible"
    }
    declared_sex: {
      description: "The sample sheet's sex (MALE or FEMALE), when it has one: the first source of the sex the sex-dependent callers use (declared, else mosdepth's coverage inference, else somalier's chrX heterozygosity)"
    }
    somalier_sex_min_x_sites: {
      description: "chrX sites somalier must have genotyped before its sex call is used"
    }
    somalier_sex_min_x_depth: {
      description: "Mean chrX depth somalier must have seen before its sex call is used"
    }
    ugc_wgw_container_registry: {
      description: "Registry holding the ugc-built images (somalier, haplogrep3)"
    }
    mito_sv_min_size: {
      description: "Allele-length difference (bp) from which a mitorsaw record is an SV in the collapsed chrM call sets (nallo's 50)"
    }
    haplogrep_tree: {
      description: "haplogrep3 tree for the mitochondrial haplogroups: phylotree-fu-rcrs@1.3 (PhyloTree 17 forensic update) or phylotree-rcrs@17.3"
    }
    run_haplogrep: {
      description: "Classify each mitorsaw haplotype with haplogrep3"
    }
    use_parabricks_deepvariant: {
      description: "Use Parabricks DeepVariant"
    }
    backend: {
      description: "Backend where the workflow will be executed",
      choices: [
        "GCP",
        "Azure",
        "AWS-HealthOmics",
        "HPC"
      ]
    }
    zones: {
      description: "Zones where compute will take place; required if backend is set to 'GCP'"
    }
    cpuPlatform: {
      description: "Optional minimum CPU platform to use for tasks on GCP"
    }
    gpuType: {
      description: "GPU type of the GPU tasks (sbatch --gres gpu:<type>:N under miniwdl-slurm). Leave unset: the ugc_wgw_resources plugin applies config.json gpu_type at submit, so the type stays out of the call cache; a value given here wins over it"
    }
    container_registry: {
      description: "Container registry for upstream images; default quay.io/pacbio"
    }
    preemptible: {
      description: "Where possible, run tasks preemptibly (no effect on HPC)"
    }
    run_methbat_report: {
      description: "methbat report: label the imprinting regions of the reference map's methbat_report_regions_tsv against their expected state (ugc_wgw_methbat_report, after the MethBat profile; off = no methbat_report output)"
    }
  }

  input {
    String ugc_wgw_version
    String sample_id
    Array[File] hifi_reads
    Array[File]? fail_reads
    File ref_map_file
    Boolean use_alignment_chunking = true
    Boolean use_gpu = false
    Boolean use_parabricks_deepvariant = false
    String? declared_sex
    Int somalier_sex_min_x_sites = 50
    Float somalier_sex_min_x_depth = 5.0
    String ugc_wgw_container_registry
    Int mito_sv_min_size = 50
    String haplogrep_tree = "phylotree-fu-rcrs@1.3"
    Boolean run_haplogrep = true
    Boolean run_methbat_report = true

    # Backend configuration
    String backend = "HPC"
    String? zones
    String? cpuPlatform
    String? gpuType
    String? container_registry
    Boolean preemptible = true
  }

  call BackendConfiguration.backend_configuration { input:
    backend = backend,
    zones = zones,
    cpuPlatform = cpuPlatform,
    gpuType = gpuType,
    container_registry = container_registry
  }

  RuntimeAttributes default_runtime_attributes = if preemptible
    then backend_configuration.spot_runtime_attributes
    else backend_configuration.on_demand_runtime_attributes

  # Reference data: singleton.wdl takes these from unpack_container_manifest over the
  # workflow data container; install-bundle.sh renders the same manifest into the map once.
  #@ except: DeclarationName
  Map[String, String] ref_map = read_map(ref_map_file)

  String ref_name = ref_map["name"]
  File ref_fasta = ref_map["fasta"]  # !FileCoercion
  File somalier_sites_vcf = ref_map["somalier_sites_vcf"]  # !FileCoercion
  File ref_index = ref_map["fasta_index"]  # !FileCoercion
  File pbmm2_index = ref_map["pbmm2_index"]  # !FileCoercion
  Float max_norm_female_chrY_depth = ref_map["max_norm_female_chrY_depth"]  # !StringCoercion
  File trgt_tandem_repeat_bed = ref_map["trgt_tandem_repeat_bed"]  # !FileCoercion
  File sawfish_exclude_bed = ref_map["sawfish_exclude_bed"]  # !FileCoercion
  File sawfish_exclude_bed_index = ref_map["sawfish_exclude_bed_index"]  # !FileCoercion
  File sawfish_expected_bed_male = ref_map["sawfish_expected_bed_male"]  # !FileCoercion
  File sawfish_expected_bed_female = ref_map["sawfish_expected_bed_female"]  # !FileCoercion
  File methbat_region_tsv = ref_map["methbat_region_tsv"]  # !FileCoercion
  File methbat_report_regions_tsv = ref_map["methbat_report_regions_tsv"]  # !FileCoercion
  String paraphase_genome_build = ref_map["paraphase_genome_build"]
  Boolean run_starphase = ref_map["run_starphase"] == "true"

  call ProcessTrgtCatalog.process_trgt_catalog { input:
    trgt_catalog = trgt_tandem_repeat_bed,
    ref_fasta = ref_fasta,
    ref_index = ref_index,
    default_runtime_attributes = default_runtime_attributes
  }

  call Upstream.upstream { input:
    sample_id = sample_id,
    hifi_reads = hifi_reads,
    fail_reads = fail_reads,
    fail_reads_bed = process_trgt_catalog.include_fail_reads_bed,
    fail_reads_bait_index = process_trgt_catalog.fail_reads_bait_index,
    ref_name = ref_name,
    ref_fasta = ref_fasta,
    pbmm2_index = pbmm2_index,
    ref_index = ref_index,
    max_norm_female_chrY_depth = max_norm_female_chrY_depth,
    paraphase_genome_build = paraphase_genome_build,
    sawfish_exclude_bed = sawfish_exclude_bed,
    sawfish_exclude_bed_index = sawfish_exclude_bed_index,
    sawfish_expected_bed_male = sawfish_expected_bed_male,
    sawfish_expected_bed_female = sawfish_expected_bed_female,
    use_alignment_chunking = use_alignment_chunking,
    single_sample = true,
    use_gpu = use_gpu,
    use_parabricks_deepvariant = use_parabricks_deepvariant,
    declared_sex = declared_sex,
    somalier_sites_vcf = somalier_sites_vcf,
    somalier_sex_min_x_sites = somalier_sex_min_x_sites,
    somalier_sex_min_x_depth = somalier_sex_min_x_depth,
    ugc_wgw_container_registry = ugc_wgw_container_registry,
    mito_sv_min_size = mito_sv_min_size,
    haplogrep_tree = haplogrep_tree,
    run_haplogrep = run_haplogrep,
    default_runtime_attributes = default_runtime_attributes
  }

  call Downstream.downstream { input:
    sample_id = sample_id,
    sex = upstream.resolved_sex,  # the resolved sex (declared, else coverage, else somalier); upstream passes mosdepth's inference
    aligned_hifi_reads = upstream.aligned_hifi_reads,
    aligned_hifi_reads_index = upstream.aligned_hifi_reads_index,
    aligned_fail_reads = upstream.aligned_fail_reads,
    aligned_fail_reads_index = upstream.aligned_fail_reads_index,
    trgt_catalog = process_trgt_catalog.full_catalog,
    small_variant_vcf = upstream.small_variant_vcf,
    small_variant_vcf_index = upstream.small_variant_vcf_index,
    sv_vcf = select_first([
      upstream.sv_vcf
    ]),
    sv_vcf_index = select_first([
      upstream.sv_vcf_index
    ]),
    ref_name = ref_name,
    ref_fasta = ref_fasta,
    ref_index = ref_index,
    sawfish_expected_bed_male = sawfish_expected_bed_male,
    sawfish_expected_bed_female = sawfish_expected_bed_female,
    methbat_region_tsv = methbat_region_tsv,
    run_starphase = run_starphase,
    default_runtime_attributes = default_runtime_attributes
  }

  # ugc-wgw: the imprinting report on the 5mC pileup (none when methbat_pileup wrote no 5mC sites)
  if (run_methbat_report && defined(downstream.cpg_pileup_bed)) {
    call MethylReport.ugc_wgw_methbat_report as methbat_report_task { input:
      cpg_pileup_bed = select_first([
        downstream.cpg_pileup_bed
      ]),
      report_regions_tsv = methbat_report_regions_tsv,
      out_prefix = "~{sample_id}.~{ref_name}",
      runtime_attributes = default_runtime_attributes
    }
  }

  Array[Array[String]] stats = [
    [
      "sample_id",
      sample_id
    ],
    [
      "read_count",
      downstream.stat_read_count
    ],
    [
      "read_length_mean",
      downstream.stat_read_length_mean
    ],
    [
      "read_length_median",
      downstream.stat_read_length_median
    ],
    [
      "read_length_n50",
      downstream.stat_read_length_n50
    ],
    [
      "read_quality_mean",
      downstream.stat_read_quality_mean
    ],
    [
      "read_quality_median",
      downstream.stat_read_quality_median
    ],
    [
      "mapped_read_count",
      downstream.stat_mapped_read_count
    ],
    [
      "mapped_read_percent",
      downstream.stat_mapped_read_percent
    ],
    [
      "gap_compressed_identity_mean",
      downstream.stat_gap_compressed_identity_mean
    ],
    [
      "gap_compressed_identity_median",
      downstream.stat_gap_compressed_identity_median
    ],
    [
      "depth_mean",
      upstream.stat_depth_mean
    ],
    [
      "inferred_sex",
      upstream.inferred_sex
    ],
    [
      "stat_phased_basepairs",
      downstream.stat_phased_basepairs
    ],
    [
      "phase_block_ng50",
      downstream.stat_phase_block_ng50
    ],
    [
      "cpg_combined_count",
      downstream.stat_cpg_combined_count
    ],
    [
      "cpg_hap1_count",
      downstream.stat_cpg_hap1_count
    ],
    [
      "cpg_hap2_count",
      downstream.stat_cpg_hap2_count
    ],
    [
      "methbat_methylated_count",
      downstream.stat_methbat_methylated_count
    ],
    [
      "methbat_unmethylated_count",
      downstream.stat_methbat_unmethylated_count
    ],
    [
      "methbat_asm_count",
      downstream.stat_methbat_asm_count
    ],
    [
      "SNV_count",
      downstream.stat_SNV_count
    ],
    [
      "TSTV_ratio",
      downstream.stat_TSTV_ratio
    ],
    [
      "HETHOM_ratio",
      downstream.stat_HETHOM_ratio
    ],
    [
      "INDEL_count",
      downstream.stat_INDEL_count
    ],
    [
      "sv_DUP_count",
      downstream.stat_sv_DUP_count
    ],
    [
      "sv_DEL_count",
      downstream.stat_sv_DEL_count
    ],
    [
      "sv_INS_count",
      downstream.stat_sv_INS_count
    ],
    [
      "sv_INV_count",
      downstream.stat_sv_INV_count
    ],
    [
      "sv_SWAP_count",
      downstream.stat_sv_SWAP_count
    ],
    [
      "sv_BND_count",
      downstream.stat_sv_BND_count
    ],
    [
      "trgt_genotyped_count",
      downstream.stat_trgt_genotyped_count
    ],
    [
      "trgt_uncalled_count",
      downstream.stat_trgt_uncalled_count
    ]
  ]

  call FilterMessages.filter_messages { input:
    message_arrays = [
      process_trgt_catalog.msg,
      upstream.msg,
      downstream.msg
    ]
  }

  call Utilities.consolidate_stats { input:
    out_prefix = sample_id,
    stats = stats,
    msg_array = filter_messages.messages,
    runtime_attributes = default_runtime_attributes
  }

  # The mirrored workflow's name and version (singleton.wdl's workflow_name/workflow_version outputs)
  call Provenance.ugc_wgw_manifest_write { input:
    ugc_wgw_version = ugc_wgw_version,
    stage = "singleton",
    subject_type = "sample",
    subject_id = sample_id,
    upstream_workflow_name = "humanwgs_singleton",
    upstream_workflow_version = "4.0.0",
    runtime_attributes = default_runtime_attributes
  }

  output {
    # consolidated stats
    File stats_file = consolidate_stats.stats_tsv
    File msg_file = consolidate_stats.messages

    # bam stats
    File read_length_plot = downstream.read_length_plot
    File read_quality_plot = downstream.read_quality_plot
    File mapq_distribution_plot = downstream.mapq_distribution_plot
    File mg_distribution_plot = downstream.mg_distribution_plot
    String stat_read_count = downstream.stat_read_count
    String stat_read_length_mean = downstream.stat_read_length_mean
    String stat_read_length_median = downstream.stat_read_length_median
    String stat_read_length_n50 = downstream.stat_read_length_n50
    String stat_read_quality_mean = downstream.stat_read_quality_mean
    String stat_read_quality_median = downstream.stat_read_quality_median
    String stat_mapped_read_count = downstream.stat_mapped_read_count
    String stat_mapped_read_percent = downstream.stat_mapped_read_percent
    String stat_gap_compressed_identity_mean = downstream.stat_gap_compressed_identity_mean
    String stat_gap_compressed_identity_median = downstream.stat_gap_compressed_identity_median

    # merged, haplotagged alignments
    File merged_haplotagged_bam = downstream.merged_haplotagged_bam
    File merged_haplotagged_bam_index = downstream.merged_haplotagged_bam_index

    # mosdepth outputs
    File mosdepth_summary = upstream.mosdepth_summary
    File mosdepth_region_bed = upstream.mosdepth_region_bed
    File mosdepth_region_bed_index = upstream.mosdepth_region_bed_index
    File mosdepth_depth_distribution_plot = upstream.mosdepth_depth_distribution_plot
    String stat_depth_mean = upstream.stat_depth_mean
    String inferred_sex = upstream.inferred_sex
    String resolved_sex = upstream.resolved_sex
    String resolved_sex_source = upstream.resolved_sex_source
    String somalier_sex = upstream.somalier_sex
    File somalier_extract = upstream.somalier_extract
    File somalier_sex_tsv = upstream.somalier_sex_tsv

    # mitorsaw's chrM calls collapsed to diploid, and the haplogroups
    File mt_small_variant_vcf = upstream.mt_small_variant_vcf
    File mt_small_variant_vcf_index = upstream.mt_small_variant_vcf_index
    File mt_sv_vcf = upstream.mt_sv_vcf
    File mt_sv_vcf_index = upstream.mt_sv_vcf_index
    File mt_haplotypes_vcf = upstream.mt_haplotypes_vcf
    File mt_haplotypes_vcf_index = upstream.mt_haplotypes_vcf_index
    File mt_haplotypes_tsv = upstream.mt_haplotypes_tsv
    File mt_summary = upstream.mt_summary
    File? mt_haplogroups = upstream.mt_haplogroups
    String? mt_haplogroup = upstream.mt_haplogroup

    # phasing stats
    File phase_stats = downstream.phase_stats
    File phase_blocks = downstream.phase_blocks
    File phase_haplotags = downstream.phase_haplotags
    String stat_phased_basepairs = downstream.stat_phased_basepairs
    String stat_phase_block_ng50 = downstream.stat_phase_block_ng50

    # methylation outputs and profile
    File? cpg_pileup_bed = downstream.cpg_pileup_bed
    File? cpg_pileup_bed_index = downstream.cpg_pileup_bed_index
    File? hmcpg_pileup_bed = downstream.hmcpg_pileup_bed
    File? hmcpg_pileup_bed_index = downstream.hmcpg_pileup_bed_index
    String stat_cpg_hap1_count = downstream.stat_cpg_hap1_count
    String stat_cpg_hap2_count = downstream.stat_cpg_hap2_count
    String stat_cpg_combined_count = downstream.stat_cpg_combined_count
    File? methbat_profile = downstream.methbat_profile
    File? methbat_report = methbat_report_task.report  # ugc-wgw: imprinting regions labelled (run_methbat_report)
    String stat_methbat_methylated_count = downstream.stat_methbat_methylated_count
    String stat_methbat_unmethylated_count = downstream.stat_methbat_unmethylated_count
    String stat_methbat_asm_count = downstream.stat_methbat_asm_count

    # sv outputs
    File phased_sv_vcf = downstream.phased_sv_vcf
    File phased_sv_vcf_index = downstream.phased_sv_vcf_index
    File sv_supporting_reads = select_first([
      upstream.sv_supporting_reads
    ])
    File sv_copynum_bedgraph = select_first([
      upstream.sv_copynum_bedgraph
    ])
    File sv_depth_bw = select_first([
      upstream.sv_depth_bw
    ])
    File sv_gc_bias_corrected_depth_bw = select_first([
      upstream.sv_gc_bias_corrected_depth_bw
    ])
    File sv_copynum_summary = select_first([
      upstream.sv_copynum_summary
    ])

    # sv stats
    String stat_sv_DUP_count = downstream.stat_sv_DUP_count
    String stat_sv_DEL_count = downstream.stat_sv_DEL_count
    String stat_sv_INS_count = downstream.stat_sv_INS_count
    String stat_sv_INV_count = downstream.stat_sv_INV_count
    String stat_sv_SWAP_count = downstream.stat_sv_SWAP_count
    String stat_sv_BND_count = downstream.stat_sv_BND_count
    File sv_stats_plot = downstream.sv_stats_plot

    # small variant outputs
    File phased_small_variant_vcf = downstream.phased_small_variant_vcf
    File phased_small_variant_vcf_index = downstream.phased_small_variant_vcf_index
    File? small_variant_gvcf = upstream.small_variant_gvcf
    File? small_variant_gvcf_index = upstream.small_variant_gvcf_index

    # small variant stats
    File small_variant_stats = downstream.small_variant_stats
    File bcftools_roh_out = downstream.bcftools_roh_out
    File bcftools_roh_bed = downstream.bcftools_roh_bed
    String stat_small_variant_SNV_count = downstream.stat_SNV_count
    String stat_small_variant_INDEL_count = downstream.stat_INDEL_count
    String stat_small_variant_TSTV_ratio = downstream.stat_TSTV_ratio
    String stat_small_variant_HETHOM_ratio = downstream.stat_HETHOM_ratio
    File snv_distribution_plot = downstream.snv_distribution_plot
    File indel_distribution_plot = downstream.indel_distribution_plot

    # trgt outputs
    File phased_trgt_vcf = downstream.trgt_vcf
    File phased_trgt_vcf_index = downstream.trgt_vcf_index
    File trgt_spanning_reads = downstream.trgt_spanning_reads
    File trgt_spanning_reads_index = downstream.trgt_spanning_reads_index
    File trgt_coverage_dropouts = downstream.trgt_coverage_dropouts
    String stat_trgt_genotyped_count = downstream.stat_trgt_genotyped_count
    String stat_trgt_uncalled_count = downstream.stat_trgt_uncalled_count

    # paraphase outputs
    File? paraphase_summary = upstream.paraphase_output_json
    File? paraphase_realigned_bam = upstream.paraphase_realigned_bam
    File? paraphase_realigned_bam_index = upstream.paraphase_realigned_bam_index
    File? paraphase_vcfs = upstream.paraphase_vcfs

    # per sample mitorsaw outputs
    File mitorsaw_vcf = upstream.mitorsaw_vcf
    File mitorsaw_vcf_index = upstream.mitorsaw_vcf_index
    File mitorsaw_hap_stats = upstream.mitorsaw_hap_stats

    # kivvi kiv2 outputs
    File? kivvi_kiv2_vcf = upstream.kivvi_kiv2_vcf
    File? kivvi_kiv2_vcf_index = upstream.kivvi_kiv2_vcf_index
    File? kivvi_kiv2_json = upstream.kivvi_kiv2_json
    File? kivvi_kiv2_realigned_bam = upstream.kivvi_kiv2_realigned_bam
    File? kivvi_kiv2_realigned_bam_index = upstream.kivvi_kiv2_realigned_bam_index
    File? kivvi_kiv2_allele_plot = upstream.kivvi_kiv2_allele_plot

    # kivvi d4z4 outputs
    File? kivvi_d4z4_vcf = upstream.kivvi_d4z4_vcf
    File? kivvi_d4z4_vcf_index = upstream.kivvi_d4z4_vcf_index
    File? kivvi_d4z4_json = upstream.kivvi_d4z4_json
    File? kivvi_d4z4_realigned_bam = upstream.kivvi_d4z4_realigned_bam
    File? kivvi_d4z4_realigned_bam_index = upstream.kivvi_d4z4_realigned_bam_index
    File? kivvi_d4z4_allele_plot = upstream.kivvi_d4z4_allele_plot

    # PGx outputs
    File? pbstarphase_summary = downstream.pbstarphase_json
    File? pbstarphase_tsv = downstream.pbstarphase_tsv

    # qc messages
    Array[String] msg = filter_messages.messages

    # ugc-wgw provenance
    String ugc_wgw_workflow_version = ugc_wgw_version
    File ugc_wgw_manifest = ugc_wgw_manifest_write.ugc_wgw_manifest
  }
}
