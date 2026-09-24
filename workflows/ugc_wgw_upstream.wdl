version 1.0

# ugc-pacbio-wgw — stage `upstream`: per-sample work before phasing (process_trgt_catalog and the upstream subworkflow as singleton.wdl calls them; per-sample sawfish call by default)
# origin: none
# see docs/DESIGN.md §5, docs/ENTRYPOINTS.md §4

import "../vendor/hifi-human-wgs-wdl/workflows/process_trgt_catalog/process_trgt_catalog.wdl" as ProcessTrgtCatalog
import "../vendor/hifi-human-wgs-wdl/workflows/upstream/upstream.wdl" as Upstream
import "../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/workflows/backend_configuration/backend_configuration.wdl" as BackendConfiguration
import "../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/workflows/filter_messages/filter_messages.wdl" as FilterMessages
import "ugc_wgw/provenance.wdl" as Provenance

workflow ugc_wgw_upstream {
  meta {
    description: "ugc-pacbio-wgw stage upstream: alignment, DeepVariant, sawfish discover (and per-sample sawfish call unless run_sawfish_call is false), Paraphase, mitorsaw, kivvi, coverage and QC for one sample, exactly as humanwgs_singleton runs them before phasing; reference data from the installed map"
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
      description: "Unaligned fail_reads BAMs (used only for TRGT at INCLUDE_FAIL_READS loci, in the downstream stage)"
    }
    ref_map_file: {
      description: "ugc-wgw reference map (docs/DESIGN.md §10.4): the files and scalars of upstream's workflow data container for one reference build, plus scatter_regions"
    }
    use_alignment_chunking: {
      description: "Whether to chunk BAM files for alignment. If false, all reads will be aligned in a single chunk."
    }
    run_sawfish_call: {
      description: "Run the single-sample sawfish call here (upstream single_sample). false: emit only the discover tarball, for cohorts that run the multi-sample call in ugc_wgw_cohort_call"
    }
    use_gpu: {
      description: "Use GPU when possible"
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
  }

  input {
    String ugc_wgw_version
    String sample_id
    Array[File] hifi_reads
    Array[File]? fail_reads
    File ref_map_file
    Boolean use_alignment_chunking = true
    Boolean run_sawfish_call = true
    Boolean use_gpu = false
    Boolean use_parabricks_deepvariant = false

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
  Float max_norm_female_chrY_depth = ref_map["max_norm_female_chrY_depth"]  # !StringCoercion
  File trgt_tandem_repeat_bed = ref_map["trgt_tandem_repeat_bed"]  # !FileCoercion
  File sawfish_exclude_bed = ref_map["sawfish_exclude_bed"]  # !FileCoercion
  File sawfish_exclude_bed_index = ref_map["sawfish_exclude_bed_index"]  # !FileCoercion
  File sawfish_expected_bed_male = ref_map["sawfish_expected_bed_male"]  # !FileCoercion
  File sawfish_expected_bed_female = ref_map["sawfish_expected_bed_female"]  # !FileCoercion
  String paraphase_genome_build = ref_map["paraphase_genome_build"]

  # as singleton.wdl: strip INCLUDE_FAIL_READS flags and build the fail-reads bait index
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
    ref_index = ref_index,
    max_norm_female_chrY_depth = max_norm_female_chrY_depth,
    paraphase_genome_build = paraphase_genome_build,
    sawfish_exclude_bed = sawfish_exclude_bed,
    sawfish_exclude_bed_index = sawfish_exclude_bed_index,
    sawfish_expected_bed_male = sawfish_expected_bed_male,
    sawfish_expected_bed_female = sawfish_expected_bed_female,
    use_alignment_chunking = use_alignment_chunking,
    single_sample = run_sawfish_call,
    use_gpu = use_gpu,
    use_parabricks_deepvariant = use_parabricks_deepvariant,
    default_runtime_attributes = default_runtime_attributes
  }

  # singleton.wdl filters the three message arrays together; the downstream stage filters
  # this stage's messages with its own, so the concatenation is the same
  call FilterMessages.filter_messages { input:
    message_arrays = [
      process_trgt_catalog.msg,
      upstream.msg
    ]
  }

  call Provenance.ugc_wgw_manifest_write { input:
    ugc_wgw_version = ugc_wgw_version,
    stage = "upstream",
    subject_type = "sample",
    subject_id = sample_id,
    upstream_workflow_name = "upstream",
    runtime_attributes = default_runtime_attributes
  }

  output {
    # alignments
    File aligned_hifi_reads = upstream.aligned_hifi_reads
    File aligned_hifi_reads_index = upstream.aligned_hifi_reads_index
    File? aligned_fail_reads = upstream.aligned_fail_reads
    File? aligned_fail_reads_index = upstream.aligned_fail_reads_index

    # mosdepth outputs
    File mosdepth_summary = upstream.mosdepth_summary
    File mosdepth_region_bed = upstream.mosdepth_region_bed
    File mosdepth_region_bed_index = upstream.mosdepth_region_bed_index
    File mosdepth_depth_distribution_plot = upstream.mosdepth_depth_distribution_plot
    String inferred_sex = upstream.inferred_sex
    String stat_depth_mean = upstream.stat_depth_mean

    # per sample sv signatures
    File discover_tar = upstream.discover_tar

    # sawfish outputs for single sample
    File? sv_vcf = upstream.sv_vcf
    File? sv_vcf_index = upstream.sv_vcf_index
    File? sv_supporting_reads = upstream.sv_supporting_reads
    File? sv_copynum_bedgraph = upstream.sv_copynum_bedgraph
    File? sv_depth_bw = upstream.sv_depth_bw
    File? sv_gc_bias_corrected_depth_bw = upstream.sv_gc_bias_corrected_depth_bw
    File? sv_copynum_summary = upstream.sv_copynum_summary

    # small variant outputs
    File small_variant_vcf = upstream.small_variant_vcf
    File small_variant_vcf_index = upstream.small_variant_vcf_index
    File? small_variant_gvcf = upstream.small_variant_gvcf
    File? small_variant_gvcf_index = upstream.small_variant_gvcf_index

    # paraphase outputs
    File? paraphase_output_json = upstream.paraphase_output_json
    File? paraphase_realigned_bam = upstream.paraphase_realigned_bam
    File? paraphase_realigned_bam_index = upstream.paraphase_realigned_bam_index
    File? paraphase_vcfs = upstream.paraphase_vcfs

    # per sample mitorsaw outputs
    File mitorsaw_vcf = upstream.mitorsaw_vcf
    File mitorsaw_vcf_index = upstream.mitorsaw_vcf_index
    File mitorsaw_hap_stats = upstream.mitorsaw_hap_stats

    # per sample kivvi kiv2 outputs
    File? kivvi_kiv2_vcf = upstream.kivvi_kiv2_vcf
    File? kivvi_kiv2_vcf_index = upstream.kivvi_kiv2_vcf_index
    File? kivvi_kiv2_json = upstream.kivvi_kiv2_json
    File? kivvi_kiv2_realigned_bam = upstream.kivvi_kiv2_realigned_bam
    File? kivvi_kiv2_realigned_bam_index = upstream.kivvi_kiv2_realigned_bam_index
    File? kivvi_kiv2_allele_plot = upstream.kivvi_kiv2_allele_plot

    # per sample kivvi d4z4 outputs
    File? kivvi_d4z4_vcf = upstream.kivvi_d4z4_vcf
    File? kivvi_d4z4_vcf_index = upstream.kivvi_d4z4_vcf_index
    File? kivvi_d4z4_json = upstream.kivvi_d4z4_json
    File? kivvi_d4z4_realigned_bam = upstream.kivvi_d4z4_realigned_bam
    File? kivvi_d4z4_realigned_bam_index = upstream.kivvi_d4z4_realigned_bam_index
    File? kivvi_d4z4_allele_plot = upstream.kivvi_d4z4_allele_plot

    # qc messages
    Array[String] msg = filter_messages.messages

    # ugc-wgw provenance
    String ugc_wgw_workflow_version = ugc_wgw_version
    File ugc_wgw_manifest = ugc_wgw_manifest_write.ugc_wgw_manifest
  }
}
