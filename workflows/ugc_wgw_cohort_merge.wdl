version 1.0

# ugc-pacbio-wgw — stage `cohort_merge`: SV merge (svx or bcftools), trgt merge + trgt-lps, optional phased small-variant merge, optional GLnexus
# origin: none
# see docs/DESIGN.md §5, §7; docs/ENTRYPOINTS.md §7

import "../vendor/hifi-human-wgs-wdl/workflows/process_trgt_catalog/process_trgt_catalog.wdl" as ProcessTrgtCatalog
import "../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/workflows/backend_configuration/backend_configuration.wdl" as BackendConfiguration
import "ugc_wgw/cohort/concat.wdl" as Concat
import "ugc_wgw/cohort/glnexus_scatter.wdl" as GlnexusScatter
import "ugc_wgw/cohort/regions.wdl" as Regions
import "ugc_wgw/cohort/sv_merge_bcftools.wdl" as BcftoolsMerge
import "ugc_wgw/cohort/svx.wdl" as Svx
import "ugc_wgw/cohort/trgt_merge_lps.wdl" as Trgt
import "ugc_wgw/provenance.wdl" as Provenance

workflow ugc_wgw_cohort_merge {
  meta {
    description: "ugc-pacbio-wgw stage cohort_merge: merge per-sample structural variants (svx or bcftools), tandem repeats (trgt merge, then trgt-lps) and optionally phased small variants across a cohort, and optionally joint-call small variants with GLnexus. Every cohort-wide step is scattered over scatter_regions."
    outputs: {
      cohort_sv_vcf: {
        description: "Cohort structural variant VCF"
      },
      cohort_sv_vcf_index: {
        description: "Index for cohort structural variant VCF"
      },
      cohort_trgt_vcf: {
        description: "Cohort TRGT VCF (trgt merge)"
      },
      cohort_trgt_vcf_index: {
        description: "Index for cohort TRGT VCF"
      },
      cohort_trgt_lps: {
        description: "Longest pure segment table (trgt-lps)"
      },
      cohort_phased_small_variant_vcf: {
        description: "bcftools merge of the per-sample phased small-variant VCFs (if merge_phased_small_variants)"
      },
      cohort_phased_small_variant_vcf_index: {
        description: "Index for merged phased small-variant VCF"
      },
      cohort_small_variant_vcf: {
        description: "GLnexus joint-called small-variant VCF, unphased (if run_glnexus)"
      },
      cohort_small_variant_vcf_index: {
        description: "Index for GLnexus VCF"
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
      description: "Cohort ID; output prefix is <cohort_id>.merged.<ref_name>"
    }
    sample_ids: {
      description: "Sample IDs, in the order of every per-sample array"
    }
    sv_vcfs: {
      description: "Per-sample sawfish VCFs (phased_sv_vcf from downstream/singleton, or the unphased sv_vcf)"
    }
    sv_vcf_indices: {
      description: "Per-sample sawfish VCF indices"
    }
    sv_merge_method: {
      description: "SV merge method",
      choices: [
        "svx",
        "bcftools"
      ]
    }
    svx_threads: {
      description: "svx --threads per shard"
    }
    svx_min_supp: {
      description: "svx --min-supp"
    }
    svx_mem_gb: {
      description: "Memory (GB) per svx shard; default 8 + 0.05 per sample"
    }
    trgt_vcfs: {
      description: "Per-sample TRGT VCFs (trgt_vcf from downstream, phased_trgt_vcf from singleton)"
    }
    trgt_vcf_indices: {
      description: "Per-sample TRGT VCF indices"
    }
    trgt_merge_threads: {
      description: "trgt merge --threads per shard"
    }
    trgt_merge_mem_gb: {
      description: "Memory (GB) per trgt merge shard; default 8 + 0.02 per sample"
    }
    trgt_lps_threads: {
      description: "trgt-lps --threads"
    }
    run_glnexus: {
      description: "Joint-call small variants with GLnexus (standalone mode); needs gvcfs and gvcf_indices"
    }
    gvcfs: {
      description: "Per-sample DeepVariant gVCFs (small_variant_gvcf)"
    }
    gvcf_indices: {
      description: "Per-sample gVCF indices"
    }
    glnexus_mem_gb: {
      description: "Memory (GB) per GLnexus shard; default 32 + 0.1 per sample"
    }
    glnexus_threads: {
      description: "Threads per GLnexus shard"
    }
    merge_phased_small_variants: {
      description: "bcftools merge of the per-sample phased small-variant VCFs, as upstream family.wdl does"
    }
    phased_small_variant_vcfs: {
      description: "Per-sample phased small-variant VCFs (phased_small_variant_vcf)"
    }
    phased_small_variant_vcf_indices: {
      description: "Per-sample phased small-variant VCF indices"
    }
    bcftools_merge_mem_gb: {
      description: "Memory (GB) per bcftools merge shard; default 8 + 0.02 per sample"
    }
    ref_map_file: {
      description: "ugc-wgw reference map (docs/DESIGN.md §10.4); keys used here: name, fasta, fasta_index, trgt_tandem_repeat_bed, scatter_regions"
    }
    ugc_wgw_container_registry: {
      description: "Registry holding the ugc-built images (svx, trgt-lps)"
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
      description: "Optional type of GPU/Accelerator to use"
    }
    container_registry: {
      description: "Container registry for upstream images; default quay.io/pacbio"
    }
    preemptible: {
      description: "Where possible, run tasks preemptibly (no effect on HPC)"
    }
  }

  input {
    String ugc_wgw_version
    String cohort_id
    Array[String] sample_ids

    Array[File] sv_vcfs
    Array[File] sv_vcf_indices
    String sv_merge_method = "svx"
    Int svx_threads = 8
    Int svx_min_supp = 1
    Int? svx_mem_gb

    Array[File] trgt_vcfs
    Array[File] trgt_vcf_indices
    Int trgt_merge_threads = 2
    Int? trgt_merge_mem_gb
    Int trgt_lps_threads = 8

    Boolean run_glnexus = true
    Array[File]? gvcfs
    Array[File]? gvcf_indices
    Int? glnexus_mem_gb
    Int glnexus_threads = 32

    Boolean merge_phased_small_variants = false
    Array[File]? phased_small_variant_vcfs
    Array[File]? phased_small_variant_vcf_indices
    Int? bcftools_merge_mem_gb

    File ref_map_file
    String ugc_wgw_container_registry

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
  String out_prefix = "~{cohort_id}.merged.~{ref_map["name"]}"

  # Memory defaults scale with cohort size; every one is overridable. Formulas are
  # first estimates to be replaced after the cohort smoke test (DESIGN §7).
  Int svx_mem = select_first([svx_mem_gb, 8 + ceil(n_samples * 0.05)])
  Int trgt_merge_mem = select_first([trgt_merge_mem_gb, 8 + ceil(n_samples * 0.02)])
  Int glnexus_mem = select_first([glnexus_mem_gb, 32 + ceil(n_samples * 0.1)])
  Int bcftools_merge_mem = select_first([bcftools_merge_mem_gb, 8 + ceil(n_samples * 0.02)])

  # Shards: whole contigs required, because trgt merge scopes by --contig only
  call Regions.ugc_wgw_regions_prepare { input:
    scatter_regions = ref_map["scatter_regions"],  # !FileCoercion
    require_whole_contigs = true,
    runtime_attributes = default_runtime_attributes
  }

  Int n_shards = length(ugc_wgw_regions_prepare.shard_names)

  # ---- structural variants -------------------------------------------------

  if (sv_merge_method == "svx") {
    # svx --trs wants the catalog without upstream's INCLUDE_FAIL_READS flags; reuse upstream's task
    call ProcessTrgtCatalog.filter_trgt_catalog { input:
      trgt_catalog = ref_map["trgt_tandem_repeat_bed"],  # !FileCoercion
      ref_fasta = ref_map["fasta"],  # !FileCoercion
      ref_index = ref_map["fasta_index"],  # !FileCoercion
      out_prefix = "fail_reads_subset",
      runtime_attributes = default_runtime_attributes
    }

    # non-BND types per shard; svx drops orphaned breakends under scoping (DESIGN §7.3)
    scatter (shard_index in range(n_shards)) {
      call Svx.ugc_wgw_svx_merge as svx_shard { input:
        vcfs = sv_vcfs,
        vcf_indices = sv_vcf_indices,
        out_prefix = "~{out_prefix}.structural_variants.~{ugc_wgw_regions_prepare.shard_names[shard_index]}",
        target_positions = ugc_wgw_regions_prepare.shard_targets[shard_index],
        svtypes = "INS,DEL,INV,DUP,CNV",
        trs_bed = filter_trgt_catalog.full_catalog,
        min_supp = svx_min_supp,
        threads = svx_threads,
        mem_gb = svx_mem,
        ugc_wgw_container_registry = ugc_wgw_container_registry,
        runtime_attributes = default_runtime_attributes
      }
    }

    # BND genome-wide, once, so that both mates of every breakend are in scope
    call Svx.ugc_wgw_svx_merge as svx_bnd { input:
      vcfs = sv_vcfs,
      vcf_indices = sv_vcf_indices,
      out_prefix = "~{out_prefix}.structural_variants.BND",
      svtypes = "BND",
      min_supp = svx_min_supp,
      threads = svx_threads,
      mem_gb = svx_mem,
      ugc_wgw_container_registry = ugc_wgw_container_registry,
      runtime_attributes = default_runtime_attributes
    }

    call Concat.ugc_wgw_bcftools_concat as concat_svx { input:
      vcfs = flatten([svx_shard.vcf, [svx_bnd.vcf]]),
      vcf_indices = flatten([svx_shard.vcf_index, [svx_bnd.vcf_index]]),
      out_prefix = "~{out_prefix}.structural_variants",
      allow_overlaps = true,
      runtime_attributes = default_runtime_attributes
    }
  }

  if (sv_merge_method == "bcftools") {
    scatter (shard_index in range(n_shards)) {
      call BcftoolsMerge.ugc_wgw_bcftools_merge as bcftools_sv_shard { input:
        vcfs = sv_vcfs,
        vcf_indices = sv_vcf_indices,
        regions_bed = ugc_wgw_regions_prepare.shard_beds[shard_index],
        out_prefix = "~{out_prefix}.structural_variants.~{ugc_wgw_regions_prepare.shard_names[shard_index]}",
        mem_gb = bcftools_merge_mem,
        runtime_attributes = default_runtime_attributes
      }
    }

    call Concat.ugc_wgw_bcftools_concat as concat_bcftools_sv { input:
      vcfs = bcftools_sv_shard.vcf,
      vcf_indices = bcftools_sv_shard.vcf_index,
      out_prefix = "~{out_prefix}.structural_variants",
      runtime_attributes = default_runtime_attributes
    }
  }

  # ---- tandem repeats ------------------------------------------------------

  scatter (shard_index in range(n_shards)) {
    call Trgt.ugc_wgw_trgt_merge { input:
      vcfs = trgt_vcfs,
      vcf_indices = trgt_vcf_indices,
      ref_fasta = ref_map["fasta"],  # !FileCoercion
      ref_index = ref_map["fasta_index"],  # !FileCoercion
      contigs = ugc_wgw_regions_prepare.shard_contigs[shard_index],
      out_prefix = "~{out_prefix}.trgt.~{ugc_wgw_regions_prepare.shard_names[shard_index]}",
      threads = trgt_merge_threads,
      mem_gb = trgt_merge_mem,
      runtime_attributes = default_runtime_attributes
    }
  }

  call Concat.ugc_wgw_bcftools_concat as concat_trgt { input:
    vcfs = ugc_wgw_trgt_merge.vcf,
    vcf_indices = ugc_wgw_trgt_merge.vcf_index,
    out_prefix = "~{out_prefix}.trgt",
    runtime_attributes = default_runtime_attributes
  }

  call Trgt.ugc_wgw_trgt_lps { input:
    vcf = concat_trgt.vcf,
    vcf_index = concat_trgt.vcf_index,
    out_prefix = "~{out_prefix}.trgt",
    threads = trgt_lps_threads,
    ugc_wgw_container_registry = ugc_wgw_container_registry,
    runtime_attributes = default_runtime_attributes
  }

  # ---- optional: phased small variants (bcftools merge, as family.wdl) ----

  if (merge_phased_small_variants) {
    scatter (shard_index in range(n_shards)) {
      call BcftoolsMerge.ugc_wgw_bcftools_merge as bcftools_small_variant_shard { input:
        vcfs = select_first([phased_small_variant_vcfs]),
        vcf_indices = select_first([phased_small_variant_vcf_indices]),
        regions_bed = ugc_wgw_regions_prepare.shard_beds[shard_index],
        out_prefix = "~{out_prefix}.small_variants.phased.~{ugc_wgw_regions_prepare.shard_names[shard_index]}",
        mem_gb = bcftools_merge_mem,
        runtime_attributes = default_runtime_attributes
      }
    }

    call Concat.ugc_wgw_bcftools_concat as concat_phased_small_variants { input:
      vcfs = bcftools_small_variant_shard.vcf,
      vcf_indices = bcftools_small_variant_shard.vcf_index,
      out_prefix = "~{out_prefix}.small_variants.phased",
      runtime_attributes = default_runtime_attributes
    }
  }

  # ---- optional: GLnexus (standalone mode) ---------------------------------

  if (run_glnexus) {
    call GlnexusScatter.ugc_wgw_glnexus_scatter { input:
      cohort_id = cohort_id,
      sample_ids = sample_ids,
      gvcfs = select_first([gvcfs]),
      gvcf_indices = select_first([gvcf_indices]),
      ref_name = ref_map["name"],
      shard_names = ugc_wgw_regions_prepare.shard_names,
      shard_beds = ugc_wgw_regions_prepare.shard_beds,
      out_prefix = "~{out_prefix}.small_variants",
      glnexus_mem_gb = glnexus_mem,
      glnexus_threads = glnexus_threads,
      default_runtime_attributes = default_runtime_attributes
    }
  }

  # ---- provenance ----------------------------------------------------------

  call Provenance.ugc_wgw_manifest_write { input:
    ugc_wgw_version = ugc_wgw_version,
    stage = "cohort_merge",
    subject_type = "cohort",
    subject_id = cohort_id,
    member_ids = sample_ids,
    runtime_attributes = default_runtime_attributes
  }

  output {
    File cohort_sv_vcf = select_first([
      concat_svx.vcf,
      concat_bcftools_sv.vcf
    ])
    File cohort_sv_vcf_index = select_first([
      concat_svx.vcf_index,
      concat_bcftools_sv.vcf_index
    ])

    File cohort_trgt_vcf = concat_trgt.vcf
    File cohort_trgt_vcf_index = concat_trgt.vcf_index
    File cohort_trgt_lps = ugc_wgw_trgt_lps.lps

    File? cohort_phased_small_variant_vcf = concat_phased_small_variants.vcf
    File? cohort_phased_small_variant_vcf_index = concat_phased_small_variants.vcf_index

    File? cohort_small_variant_vcf = ugc_wgw_glnexus_scatter.vcf
    File? cohort_small_variant_vcf_index = ugc_wgw_glnexus_scatter.vcf_index

    String ugc_wgw_workflow_version = ugc_wgw_version
    File ugc_wgw_manifest = ugc_wgw_manifest_write.ugc_wgw_manifest
  }
}
