version 1.0

# ugc-pacbio-wgw — stage `cohort_call`: GLnexus joint calling scattered by region, one-pass split per sample, optional multi-sample sawfish call
# origin: none
# see docs/DESIGN.md §5, §7.2; docs/ENTRYPOINTS.md §6

import "../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/tasks/bcftools.wdl" as Bcftools
import "../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/tasks/sawfish.wdl" as Sawfish
import "../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/workflows/backend_configuration/backend_configuration.wdl" as BackendConfiguration
import "ugc_wgw/cohort/glnexus_scatter.wdl" as GlnexusScatter
import "ugc_wgw/cohort/regions.wdl" as Regions
import "ugc_wgw/cohort/split_by_sample.wdl" as Split
import "ugc_wgw/provenance.wdl" as Provenance

workflow ugc_wgw_cohort_call {
  meta {
    description: "ugc-pacbio-wgw stage cohort_call: joint-call small variants across a cohort with GLnexus (gVCF slicing per sample, GLnexus per shard, concat) and split the result per sample for phasing in ugc_wgw_downstream; optionally joint-call structural variants with sawfish for small cohorts, as upstream joint.wdl does"
    outputs: {
      cohort_small_variant_vcf: {
        description: "Joint-called multi-sample small variant VCF"
      },
      cohort_small_variant_vcf_index: {
        description: "Index for joint-called small variant VCF"
      },
      split_small_variant_vcfs: {
        description: "Joint-called small variant VCFs split by sample, in sample_ids order"
      },
      split_small_variant_vcf_indices: {
        description: "Indices for split small variant VCFs"
      },
      cohort_sv_vcf: {
        description: "Joint-called structural variant VCF (if run_sawfish_joint_call)"
      },
      cohort_sv_vcf_index: {
        description: "Index for joint-called structural variant VCF"
      },
      sv_supporting_reads: {
        description: "Supporting reads for structural variants"
      },
      split_sv_vcfs: {
        description: "Joint-called structural variant VCFs split by sample, in sample_ids order"
      },
      split_sv_vcf_indices: {
        description: "Indices for split structural variant VCFs"
      },
      sv_copynum_bedgraph: {
        description: "CNV copy number BEDGraph per sample"
      },
      sv_depth_bw: {
        description: "CNV depth BigWig per sample"
      },
      sv_gc_bias_corrected_depth_bw: {
        description: "CNV GC-bias corrected depth BigWig per sample"
      },
      sv_copynum_summary: {
        description: "CNV copy number summary JSON per sample"
      },
      ugc_wgw_workflow_version: {
        description: "ugc_wgw_version echoed"
      },
      ugc_wgw_manifest: {
        description: "Workflow-side provenance record"
      }
    }
  }

  parameter_meta {
    ugc_wgw_version: {
      description: "ugc-pacbio-wgw version the run was submitted with"
    }
    cohort_id: {
      description: "Cohort ID; output prefix is <cohort_id>.joint.<ref_name>"
    }
    sample_ids: {
      description: "Sample IDs, in the order of every per-sample array"
    }
    gvcfs: {
      description: "Per-sample DeepVariant gVCFs (small_variant_gvcf from ugc_wgw_upstream)"
    }
    gvcf_indices: {
      description: "Per-sample gVCF indices"
    }
    ref_map_file: {
      description: "ugc-wgw reference map (docs/DESIGN.md §10.4); keys used here: name, fasta, fasta_index, scatter_regions"
    }
    scatter_regions_file: {
      description: "Override for the scatter_regions TSV (default: the ugc-wgw reference map entry); finer shards are allowed here"
    }
    glnexus_mem_gb: {
      description: "Memory (GB) per GLnexus shard; default 32 + 0.1 per sample"
    }
    glnexus_threads: {
      description: "Threads per GLnexus shard"
    }
    split_keep_homref: {
      description: "true: per-sample VCFs keep hom-ref sites and drop only uncalled genotypes (upstream parity); false: keep only sites where the sample carries an alt allele"
    }
    split_mem_gb: {
      description: "Memory (GB) for the split; default 8 + 0.02 per sample"
    }
    run_sawfish_joint_call: {
      description: "Joint-call structural variants across the cohort with sawfish (one job over all BAMs and discover tarballs; small cohorts only); needs discover_tars, aligned_bams, aligned_bam_indices"
    }
    discover_tars: {
      description: "Per-sample sawfish discover tarballs (discover_tar from ugc_wgw_upstream)"
    }
    aligned_bams: {
      description: "Per-sample aligned BAMs (aligned_hifi_reads from ugc_wgw_upstream)"
    }
    aligned_bam_indices: {
      description: "Per-sample aligned BAM indices"
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
    String cohort_id
    Array[String] sample_ids
    Array[File] gvcfs
    Array[File] gvcf_indices
    File ref_map_file
    File? scatter_regions_file
    Int? glnexus_mem_gb
    Int glnexus_threads = 32
    Boolean split_keep_homref = true
    Int? split_mem_gb

    Boolean run_sawfish_joint_call = false
    Array[File]? discover_tars
    Array[File]? aligned_bams
    Array[File]? aligned_bam_indices

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

  #@ except: DeclarationName
  Map[String, String] ref_map = read_map(ref_map_file)

  Int n_samples = length(sample_ids)
  String joint_prefix = "~{cohort_id}.joint.~{ref_map["name"]}"

  # Memory defaults scale with cohort size; first estimates until the cohort smoke test (DESIGN §7.2)
  Int glnexus_mem = select_first([glnexus_mem_gb, 32 + ceil(n_samples * 0.1)])
  Int split_mem = select_first([split_mem_gb, 8 + ceil(n_samples * 0.02)])

  File scatter_regions = select_first([scatter_regions_file, ref_map["scatter_regions"]])  # !FileCoercion

  # GLnexus scopes by BED, so shards need not be whole contigs here
  call Regions.ugc_wgw_regions_prepare { input:
    scatter_regions = scatter_regions,
    runtime_attributes = default_runtime_attributes
  }

  # ---- small variants: GLnexus per shard, concat, split per sample ---------

  call GlnexusScatter.ugc_wgw_glnexus_scatter { input:
    cohort_id = "~{cohort_id}.joint",
    sample_ids = sample_ids,
    gvcfs = gvcfs,
    gvcf_indices = gvcf_indices,
    ref_name = ref_map["name"],
    shard_names = ugc_wgw_regions_prepare.shard_names,
    shard_beds = ugc_wgw_regions_prepare.shard_beds,
    out_prefix = "~{joint_prefix}.small_variants",
    glnexus_mem_gb = glnexus_mem,
    glnexus_threads = glnexus_threads,
    default_runtime_attributes = default_runtime_attributes
  }

  # Output names as upstream joint.wdl: <sample_id>.<cohort_id>.joint.<ref>.small_variants.vcf.gz
  String glnexus_vcf_basename = basename(ugc_wgw_glnexus_scatter.vcf, ".vcf.gz")
  scatter (sample_id in sample_ids) {
    String split_glnexus_vcf_name = "~{sample_id}.~{glnexus_vcf_basename}.vcf.gz"
    String split_glnexus_vcf_index_name = "~{sample_id}.~{glnexus_vcf_basename}.vcf.gz.tbi"
  }

  call Split.ugc_wgw_bcftools_split as split_glnexus { input:
    vcf = ugc_wgw_glnexus_scatter.vcf,
    vcf_index = ugc_wgw_glnexus_scatter.vcf_index,
    sample_ids = sample_ids,
    split_vcf_names = split_glnexus_vcf_name,
    split_vcf_index_names = split_glnexus_vcf_index_name,
    keep_homref = split_keep_homref,
    mem_gb = split_mem,
    runtime_attributes = default_runtime_attributes
  }

  # ---- optional: multi-sample sawfish call, as upstream joint.wdl ----------

  if (run_sawfish_joint_call) {
    # In order to properly delocalize the outputs of sawfish_call in cloud engines
    # we need to generate the names of the outputs and pass these to the call.
    scatter (sample_id in sample_ids) {
      String copynum_bedgraph_name = "~{sample_id}.~{joint_prefix}.structural_variants.copynum.bedgraph"
      String depth_bw_name = "~{sample_id}.~{joint_prefix}.structural_variants.depth.bw"
      String gc_bias_corrected_depth_bw_name = "~{sample_id}.~{joint_prefix}.structural_variants.gc_bias_corrected_depth.bw"
      String copynum_summary_name = "~{sample_id}.~{joint_prefix}.structural_variants.copynum.summary.json"
    }

    call Sawfish.sawfish_call { input:
      sample_ids = sample_ids,
      discover_tars = select_first([discover_tars]),
      aligned_bams = select_first([aligned_bams]),
      aligned_bam_indices = select_first([aligned_bam_indices]),
      ref_fasta = ref_map["fasta"],  # !FileCoercion
      ref_index = ref_map["fasta_index"],  # !FileCoercion
      out_prefix = "~{joint_prefix}.structural_variants",
      copynum_bedgraph_names = copynum_bedgraph_name,
      depth_bw_names = depth_bw_name,
      gc_bias_corrected_depth_bw_names = gc_bias_corrected_depth_bw_name,
      copynum_summary_names = copynum_summary_name,
      runtime_attributes = default_runtime_attributes
    }

    String sv_vcf_basename = basename(sawfish_call.vcf, ".vcf.gz")
    scatter (sample_id in sample_ids) {
      String split_sv_vcf_name = "~{sample_id}.~{sv_vcf_basename}.vcf.gz"
      String split_sv_vcf_index_name = "~{sample_id}.~{sv_vcf_basename}.vcf.gz.tbi"
    }

    # upstream's task (N small passes over a small VCF), keeping every site as joint.wdl does
    call Bcftools.split_vcf_by_sample as split_sawfish { input:
      sample_ids = sample_ids,
      vcf = sawfish_call.vcf,
      vcf_index = sawfish_call.vcf_index,
      split_vcf_names = split_sv_vcf_name,
      split_vcf_index_names = split_sv_vcf_index_name,
      exclude_uncalled = false,
      runtime_attributes = default_runtime_attributes
    }
  }

  # ---- provenance ----------------------------------------------------------

  call Provenance.ugc_wgw_manifest_write { input:
    ugc_wgw_version = ugc_wgw_version,
    stage = "cohort_call",
    subject_type = "cohort",
    subject_id = cohort_id,
    member_ids = sample_ids,
    upstream_workflow_name = "joint",
    runtime_attributes = default_runtime_attributes
  }

  output {
    File cohort_small_variant_vcf = ugc_wgw_glnexus_scatter.vcf
    File cohort_small_variant_vcf_index = ugc_wgw_glnexus_scatter.vcf_index
    Array[File] split_small_variant_vcfs = split_glnexus.split_vcfs
    Array[File] split_small_variant_vcf_indices = split_glnexus.split_vcf_indices

    File? cohort_sv_vcf = sawfish_call.vcf
    File? cohort_sv_vcf_index = sawfish_call.vcf_index
    File? sv_supporting_reads = sawfish_call.supporting_reads
    Array[File]? split_sv_vcfs = split_sawfish.split_vcfs
    Array[File]? split_sv_vcf_indices = split_sawfish.split_vcf_indices
    Array[File]? sv_copynum_bedgraph = sawfish_call.copynum_bedgraph
    Array[File]? sv_depth_bw = sawfish_call.depth_bw
    Array[File]? sv_gc_bias_corrected_depth_bw = sawfish_call.gc_bias_corrected_depth_bw
    Array[File]? sv_copynum_summary = sawfish_call.copynum_summary

    String ugc_wgw_workflow_version = ugc_wgw_version
    File ugc_wgw_manifest = ugc_wgw_manifest_write.ugc_wgw_manifest
  }
}
