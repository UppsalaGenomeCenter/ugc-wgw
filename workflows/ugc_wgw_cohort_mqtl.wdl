version 1.0

# ugc-pacbio-wgw — stage `cohort_mqtl`: structural variants and tandem repeats of cohort_merge tested in cis against the CpG islands' methylation of cohort_methyl; a research screen
# origin: none
# see docs/DESIGN.md §5, §18 (2026-11-13); docs/ENTRYPOINTS.md §9.5; docs/expansion/cohort-mqtl.design.md

import "../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/workflows/backend_configuration/backend_configuration.wdl" as BackendConfiguration
import "ugc_wgw/cohort/regions.wdl" as Regions
import "ugc_wgw/mqtl/cohort.wdl" as Mqtl
import "ugc_wgw/provenance.wdl" as Provenance

workflow ugc_wgw_cohort_mqtl {
  meta {
    description: "ugc-pacbio-wgw stage cohort_mqtl: does carrying a structural variant or a longer tandem repeat go with a different methylation of the CpG island next to it? The members' scored profiles of cohort_methyl become one island-by-member matrix; per region shard every autosomal island with enough covered members and some spread is tested against every SV (carrier dosage) and repeat (mean allele length) of cohort_merge within window_bp, methylation rank-normalised, both sides residualised on sex, the ancestry components of cohort_qc and the blood cell fractions of the deconvolution where every member has them; Benjamini-Hochberg over the run. Optional: no stage reads it. A research screen in cis, in one tissue, at CpG islands; an association is not a mechanism and nothing here is a diagnostic result. Below min_members the tables hold their headers and a note."
    outputs: {
      mqtl_matrix: {
        description: "<cohort_id>.mqtl.<ref>.matrix.tsv.gz: a row per CpG island, a column per member with the island's mean methylation (0 to 100)"
      },
      mqtl_covariates: {
        description: "<cohort_id>.mqtl.<ref>.covariates.tsv: per cohort member whether it is used, sex, ancestry components, blood cell fractions"
      },
      mqtl_pairs: {
        description: "<cohort_id>.mqtl.<ref>.pairs.tsv.gz: every tested island-variant pair with r, t, p, q and the effect in methylation points"
      },
      mqtl_top: {
        description: "<cohort_id>.mqtl.<ref>.top.tsv: the pairs at or below the false-discovery rate, by p-value"
      },
      mqtl_summary: {
        description: "<cohort_id>.mqtl.<ref>.summary.tsv: members, islands, variants, pairs, hits, the covariates in the model, the inflation factor, the parameters"
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
      description: "ugc-pacbio-wgw version"
    }
    cohort_id: {
      description: "Cohort ID; output prefix"
    }
    sample_ids: {
      description: "The cohort's members, in order"
    }
    sample_sexes: {
      description: "XX, XY or empty per member (the sex the read-based track resolved)"
    }
    scored_profiles: {
      description: "cohort_methyl's scored_profiles_annotated: one per member that stage scored (the sample ID is read from the file name)"
    }
    methyl_deconvolution: {
      description: "cohort_methyl's cell-type table, when its run has one"
    }
    sv_vcf: {
      description: "cohort_merge's SV VCF"
    }
    sv_vcf_index: {
      description: "Its index"
    }
    trgt_vcf: {
      description: "cohort_merge's merged TRGT VCF"
    }
    trgt_vcf_index: {
      description: "Its index"
    }
    ancestry: {
      description: "cohort_qc's somalier ancestry table, when the cohort has a cohort_qc run with one"
    }
    window_bp: {
      description: "A variant and an island are tested as a pair within this distance"
    }
    min_carriers: {
      description: "A variant is tested when at least this many members differ from the commonest genotype and at least as many have it"
    }
    min_members: {
      description: "Below this many members with a profile nothing is tested"
    }
    min_island_members_fraction: {
      description: "An island is tested when this share of the used members has a value"
    }
    min_island_sd: {
      description: "An island is tested when the standard deviation of its methylation over the members is at least this many points"
    }
    min_island_coverage: {
      description: "A member's value for an island counts when its median coverage there is at least this"
    }
    fdr: {
      description: "The false-discovery rate of the top table"
    }
    exclude_ids: {
      description: "Members left out (relatives, a failed sample)"
    }
    test_svs: {
      description: "Test structural variants"
    }
    test_trs: {
      description: "Test tandem repeats"
    }
    matrix_mem_gb: {
      description: "Memory (GB) of the matrix task"
    }
    shard_mem_gb: {
      description: "Memory (GB) of each shard's tests"
    }
    gather_mem_gb: {
      description: "Memory (GB) of the gather"
    }
    ref_map_file: {
      description: "Reference map: name and scatter_regions are used"
    }
    backend: {
      description: "Backend where the workflow will be executed ['GCP', 'Azure', 'AWS-HealthOmics', 'HPC']"
    }
    zones: {
      description: "Zones where compute will take place; required if backend is set to 'GCP'"
    }
    cpuPlatform: {
      description: "Optional minimum CPU platform to use for tasks on GCP"
    }
    gpuType: {
      description: "GPU type to use; unused in this stage"
    }
    container_registry: {
      description: "Container registry where workflow images are hosted"
    }
    preemptible: {
      description: "Where possible, run tasks preemptibly"
    }
  }

  input {
    String ugc_wgw_version
    String cohort_id
    Array[String] sample_ids
    Array[String] sample_sexes
    Array[File] scored_profiles
    File? methyl_deconvolution
    File sv_vcf
    File sv_vcf_index
    File trgt_vcf
    File trgt_vcf_index
    File? ancestry

    Int window_bp = 100000
    Int min_carriers = 10
    Int min_members = 40
    Float min_island_members_fraction = 0.8
    Float min_island_sd = 2.0
    Int min_island_coverage = 8
    Float fdr = 0.05
    Array[String] exclude_ids = []
    Boolean test_svs = true
    Boolean test_trs = true
    Int matrix_mem_gb = 8
    Int shard_mem_gb = 8
    Int gather_mem_gb = 4

    File ref_map_file

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

  Map[String, String] ref_map = read_map(ref_map_file)

  String ref_name = ref_map["name"]
  String out_prefix = "~{cohort_id}.mqtl.~{ref_name}"

  call Mqtl.ugc_wgw_mqtl_matrix { input:
    profiles = scored_profiles,
    ref_name = ref_name,
    sample_ids = sample_ids,
    sample_sexes = sample_sexes,
    deconvolution = methyl_deconvolution,
    ancestry = ancestry,
    exclude_ids = exclude_ids,
    min_island_coverage = min_island_coverage,
    out_prefix = out_prefix,
    mem_gb = matrix_mem_gb,
    runtime_attributes = default_runtime_attributes
  }

  # per shard: the islands that start in it, against the variants of the shard widened by the window
  call Regions.ugc_wgw_regions_prepare { input:
    scatter_regions = ref_map["scatter_regions"],  # !FileCoercion
    require_whole_contigs = false,
    runtime_attributes = default_runtime_attributes
  }

  scatter (shard_index in range(length(ugc_wgw_regions_prepare.shard_names))) {
    call Mqtl.ugc_wgw_mqtl_shard { input:
      matrix = ugc_wgw_mqtl_matrix.matrix,
      covariates = ugc_wgw_mqtl_matrix.covariates,
      sv_vcf = sv_vcf,
      sv_vcf_index = sv_vcf_index,
      trgt_vcf = trgt_vcf,
      trgt_vcf_index = trgt_vcf_index,
      regions_bed = ugc_wgw_regions_prepare.shard_beds[shard_index],
      window_bp = window_bp,
      min_carriers = min_carriers,
      min_members = min_members,
      min_island_members_fraction = min_island_members_fraction,
      min_island_sd = min_island_sd,
      test_svs = test_svs,
      test_trs = test_trs,
      out_prefix = "~{out_prefix}.~{ugc_wgw_regions_prepare.shard_names[shard_index]}",
      mem_gb = shard_mem_gb,
      runtime_attributes = default_runtime_attributes
    }
  }

  call Mqtl.ugc_wgw_mqtl_gather { input:
    pairs = ugc_wgw_mqtl_shard.pairs,
    counts = ugc_wgw_mqtl_shard.counts,
    fdr = fdr,
    parameters = "window_bp=~{window_bp};min_carriers=~{min_carriers};min_members=~{min_members};min_island_members_fraction=~{min_island_members_fraction};min_island_sd=~{min_island_sd};min_island_coverage=~{min_island_coverage};test_svs=~{test_svs};test_trs=~{test_trs}",
    out_prefix = out_prefix,
    mem_gb = gather_mem_gb,
    runtime_attributes = default_runtime_attributes
  }

  call Provenance.ugc_wgw_manifest_write { input:
    ugc_wgw_version = ugc_wgw_version,
    stage = "cohort_mqtl",
    subject_type = "cohort",
    subject_id = cohort_id,
    member_ids = sample_ids,
    runtime_attributes = default_runtime_attributes
  }

  output {
    File mqtl_matrix = ugc_wgw_mqtl_matrix.matrix
    File mqtl_covariates = ugc_wgw_mqtl_matrix.covariates
    File mqtl_pairs = ugc_wgw_mqtl_gather.pairs_gz
    File mqtl_top = ugc_wgw_mqtl_gather.top
    File mqtl_summary = ugc_wgw_mqtl_gather.summary
    String ugc_wgw_workflow_version = ugc_wgw_version
    File ugc_wgw_manifest = ugc_wgw_manifest_write.ugc_wgw_manifest
  }
}
