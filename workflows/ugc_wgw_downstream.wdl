version 1.0

# ugc-pacbio-wgw — stage `downstream`: phasing and everything after (downstream subworkflow, stats table, messages file), as singleton.wdl runs them
# origin: none
# see docs/DESIGN.md §5, docs/ENTRYPOINTS.md §5

import "../vendor/hifi-human-wgs-wdl/workflows/downstream/downstream.wdl" as Downstream
import "../vendor/hifi-human-wgs-wdl/workflows/process_trgt_catalog/process_trgt_catalog.wdl" as ProcessTrgtCatalog
import "../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/tasks/utilities.wdl" as Utilities
import "../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/workflows/backend_configuration/backend_configuration.wdl" as BackendConfiguration
import "../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/workflows/filter_messages/filter_messages.wdl" as FilterMessages
import "ugc_wgw/provenance.wdl" as Provenance

workflow ugc_wgw_downstream {
  meta {
    description: "ugc-pacbio-wgw stage downstream: HiPhase, pbjam BAM stats, TRGT, variant stats, MethBat methylation and StarPhase PGx for one sample from an aligned BAM plus a small-variant VCF and an SV VCF chosen by the driver, exactly as humanwgs_singleton runs them after the upstream subworkflow; reference data from the installed map"
  }

  parameter_meta {
    ugc_wgw_version: {
      description: "ugc-pacbio-wgw version the run was submitted with; echoed as ugc_wgw_workflow_version"
    }
    sample_id: {
      description: "Sample ID"
    }
    sex: {
      description: "Sex used by TRGT and the stats table: pass the inferred_sex output of the upstream stage (MALE, FEMALE, or empty for XX)",
      choices: [
        "MALE",
        "FEMALE",
        ""
      ]
    }
    aligned_hifi_reads: {
      description: "Aligned hifi_reads BAM (upstream stage)"
    }
    aligned_hifi_reads_index: {
      description: "Aligned hifi_reads BAM index"
    }
    aligned_fail_reads: {
      description: "Aligned fail_reads BAM (upstream stage, optional)"
    }
    aligned_fail_reads_index: {
      description: "Aligned fail_reads BAM index"
    }
    small_variant_vcf: {
      description: "Single-sample small variant VCF: the DeepVariant VCF, or the sample's split of the cohort_call GLnexus VCF"
    }
    small_variant_vcf_index: {
      description: "Small variant VCF index"
    }
    sv_vcf: {
      description: "Single-sample structural variant VCF: the per-sample sawfish VCF, or the sample's split of the cohort_call joint sawfish VCF"
    }
    sv_vcf_index: {
      description: "Structural variant VCF index"
    }
    ref_map_file: {
      description: "ugc-wgw reference map (docs/DESIGN.md §10.4): the files and scalars of upstream's workflow data container for one reference build, plus scatter_regions"
    }
    stat_depth_mean: {
      description: "stat_depth_mean output of the upstream stage, for the stats table (NA if absent)"
    }
    upstream_msg: {
      description: "msg output of the upstream stage, so that the messages file matches humanwgs_singleton"
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
  }

  input {
    String ugc_wgw_version
    String sample_id
    String sex
    File aligned_hifi_reads
    File aligned_hifi_reads_index
    File? aligned_fail_reads
    File? aligned_fail_reads_index
    File small_variant_vcf
    File small_variant_vcf_index
    File sv_vcf
    File sv_vcf_index
    File ref_map_file
    String? stat_depth_mean
    Array[String] upstream_msg = []

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

  # Reference data as singleton.wdl takes it from the workflow data container (ugc_wgw_singleton.wdl)
  #@ except: DeclarationName
  Map[String, String] ref_map = read_map(ref_map_file)

  String ref_name = ref_map["name"]
  File ref_fasta = ref_map["fasta"]  # !FileCoercion
  File ref_index = ref_map["fasta_index"]  # !FileCoercion
  File trgt_tandem_repeat_bed = ref_map["trgt_tandem_repeat_bed"]  # !FileCoercion
  File sawfish_expected_bed_male = ref_map["sawfish_expected_bed_male"]  # !FileCoercion
  File sawfish_expected_bed_female = ref_map["sawfish_expected_bed_female"]  # !FileCoercion
  File methbat_region_tsv = ref_map["methbat_region_tsv"]  # !FileCoercion
  Boolean run_starphase = ref_map["run_starphase"] == "true"

  # the same flag-stripped catalog process_trgt_catalog gives singleton.wdl (the bait index is not needed here)
  call ProcessTrgtCatalog.filter_trgt_catalog { input:
    trgt_catalog = trgt_tandem_repeat_bed,
    ref_fasta = ref_fasta,
    ref_index = ref_index,
    out_prefix = "fail_reads_subset",
    runtime_attributes = default_runtime_attributes
  }

  call Downstream.downstream { input:
    sample_id = sample_id,
    sex = sex,
    aligned_hifi_reads = aligned_hifi_reads,
    aligned_hifi_reads_index = aligned_hifi_reads_index,
    aligned_fail_reads = aligned_fail_reads,
    aligned_fail_reads_index = aligned_fail_reads_index,
    trgt_catalog = filter_trgt_catalog.full_catalog,
    small_variant_vcf = small_variant_vcf,
    small_variant_vcf_index = small_variant_vcf_index,
    sv_vcf = sv_vcf,
    sv_vcf_index = sv_vcf_index,
    ref_name = ref_name,
    ref_fasta = ref_fasta,
    ref_index = ref_index,
    sawfish_expected_bed_male = sawfish_expected_bed_male,
    sawfish_expected_bed_female = sawfish_expected_bed_female,
    methbat_region_tsv = methbat_region_tsv,
    run_starphase = run_starphase,
    default_runtime_attributes = default_runtime_attributes
  }

  # singleton.wdl's stats table; depth_mean and inferred_sex come from the upstream stage's outputs
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
      select_first([
        stat_depth_mean,
        "NA"
      ])
    ],
    [
      "inferred_sex",
      sex
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

  # upstream_msg is the upstream stage's already filtered messages; filtering again is a no-op for them
  call FilterMessages.filter_messages { input:
    message_arrays = [
      upstream_msg,
      downstream.msg
    ]
  }

  call Utilities.consolidate_stats { input:
    out_prefix = sample_id,
    stats = stats,
    msg_array = filter_messages.messages,
    runtime_attributes = default_runtime_attributes
  }

  call Provenance.ugc_wgw_manifest_write { input:
    ugc_wgw_version = ugc_wgw_version,
    stage = "downstream",
    subject_type = "sample",
    subject_id = sample_id,
    upstream_workflow_name = "downstream",
    runtime_attributes = default_runtime_attributes
  }

  output {
    # consolidated stats
    File stats_file = consolidate_stats.stats_tsv
    File msg_file = consolidate_stats.messages

    # hiphase outputs
    File merged_haplotagged_bam = downstream.merged_haplotagged_bam
    File merged_haplotagged_bam_index = downstream.merged_haplotagged_bam_index
    File phased_small_variant_vcf = downstream.phased_small_variant_vcf
    File phased_small_variant_vcf_index = downstream.phased_small_variant_vcf_index
    File phased_sv_vcf = downstream.phased_sv_vcf
    File phased_sv_vcf_index = downstream.phased_sv_vcf_index
    File phase_stats = downstream.phase_stats
    File phase_blocks = downstream.phase_blocks
    File phase_haplotags = downstream.phase_haplotags
    String stat_phased_basepairs = downstream.stat_phased_basepairs
    String stat_phase_block_ng50 = downstream.stat_phase_block_ng50

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
    File trgt_coverage_dropouts = downstream.trgt_coverage_dropouts

    # small variant stats
    File small_variant_stats = downstream.small_variant_stats
    File bcftools_roh_out = downstream.bcftools_roh_out
    File bcftools_roh_bed = downstream.bcftools_roh_bed
    String stat_SNV_count = downstream.stat_SNV_count
    String stat_INDEL_count = downstream.stat_INDEL_count
    String stat_TSTV_ratio = downstream.stat_TSTV_ratio
    String stat_HETHOM_ratio = downstream.stat_HETHOM_ratio
    File snv_distribution_plot = downstream.snv_distribution_plot
    File indel_distribution_plot = downstream.indel_distribution_plot

    # sv stats
    String stat_sv_DUP_count = downstream.stat_sv_DUP_count
    String stat_sv_DEL_count = downstream.stat_sv_DEL_count
    String stat_sv_INS_count = downstream.stat_sv_INS_count
    String stat_sv_INV_count = downstream.stat_sv_INV_count
    String stat_sv_BND_count = downstream.stat_sv_BND_count
    String stat_sv_SWAP_count = downstream.stat_sv_SWAP_count
    File sv_stats_plot = downstream.sv_stats_plot

    # trgt outputs
    File trgt_vcf = downstream.trgt_vcf
    File trgt_vcf_index = downstream.trgt_vcf_index
    File trgt_spanning_reads = downstream.trgt_spanning_reads
    File trgt_spanning_reads_index = downstream.trgt_spanning_reads_index
    String stat_trgt_genotyped_count = downstream.stat_trgt_genotyped_count
    String stat_trgt_uncalled_count = downstream.stat_trgt_uncalled_count

    # methylation outputs and profile
    File? cpg_pileup_bed = downstream.cpg_pileup_bed
    File? cpg_pileup_bed_index = downstream.cpg_pileup_bed_index
    File? hmcpg_pileup_bed = downstream.hmcpg_pileup_bed
    File? hmcpg_pileup_bed_index = downstream.hmcpg_pileup_bed_index
    String stat_cpg_hap1_count = downstream.stat_cpg_hap1_count
    String stat_cpg_hap2_count = downstream.stat_cpg_hap2_count
    String stat_cpg_combined_count = downstream.stat_cpg_combined_count
    File? methbat_profile = downstream.methbat_profile
    String stat_methbat_methylated_count = downstream.stat_methbat_methylated_count
    String stat_methbat_unmethylated_count = downstream.stat_methbat_unmethylated_count
    String stat_methbat_asm_count = downstream.stat_methbat_asm_count

    # pbstarphase outputs
    File? pbstarphase_json = downstream.pbstarphase_json
    File? pbstarphase_tsv = downstream.pbstarphase_tsv

    # qc messages
    Array[String] msg = filter_messages.messages

    # ugc-wgw provenance
    String ugc_wgw_workflow_version = ugc_wgw_version
    File ugc_wgw_manifest = ugc_wgw_manifest_write.ugc_wgw_manifest
  }
}
