version 1.0

# ugc-pacbio-wgw — stage `cohort_freq`: in-house allele-frequency resource from the cohort VCFs (sites-only, sex strata, summary)
# origin: none
# see docs/DESIGN.md §5, §7.7; docs/ENTRYPOINTS.md §9

import "../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/workflows/backend_configuration/backend_configuration.wdl" as BackendConfiguration
import "ugc_wgw/cohort/concat.wdl" as Concat
import "ugc_wgw/cohort/freq.wdl" as Freq
import "ugc_wgw/cohort/regions.wdl" as Regions
import "ugc_wgw/provenance.wdl" as Provenance

workflow ugc_wgw_cohort_freq {
  meta {
    description: "ugc-pacbio-wgw stage cohort_freq: allele counts and frequencies (AC, AN, AF, NS, AC_Hom, AC_Het, AC_Hemi; overall and per sex group) of the cohort's joint small-variant VCF and merged SV VCF, written as sites-only VCFs with a summary. Scattered over scatter_regions; an SV absent from a sample's call set counts as reference."
    outputs: {
      small_variant_freq_vcf: {
        description: "Sites-only small-variant frequency VCF (absent when the cohort has no joint small-variant VCF)"
      },
      small_variant_freq_vcf_index: {
        description: "Index for small_variant_freq_vcf"
      },
      sv_freq_vcf: {
        description: "Sites-only structural-variant frequency VCF (SUPP and CF kept, IDLIST dropped)"
      },
      sv_freq_vcf_index: {
        description: "Index for sv_freq_vcf"
      },
      freq_summary: {
        description: "Long-format summary TSV (records by type, contig class, filter, AF bin; singletons; call rate; strata)"
      },
      freq_samples: {
        description: "Members with the sex group used and their presence in each VCF"
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
      description: "Cohort ID"
    }
    sample_ids: {
      description: "Cohort members; every one must be a sample of both VCFs"
    }
    sample_sexes: {
      description: "XX, XY or empty per member, aligned with sample_ids (the driver: sheet sex, else inferred sex, else empty)"
    }
    sample_sex_sources: {
      description: "sheet, inferred or empty per member, for the samples table"
    }
    small_variant_vcf: {
      description: "Cohort joint small-variant VCF (GLnexus): cohort_call in joint mode, cohort_merge with run_glnexus in standalone mode; omit when there is none"
    }
    small_variant_vcf_index: {
      description: "Index for small_variant_vcf"
    }
    sv_vcf: {
      description: "Cohort merged SV VCF (cohort_merge)"
    }
    sv_vcf_index: {
      description: "Index for sv_vcf"
    }
    freq_threads: {
      description: "Threads per shard task"
    }
    freq_mem_gb: {
      description: "Memory (GB) per shard task; the pipeline streams"
    }
    ref_map_file: {
      description: "ugc-wgw reference map (docs/DESIGN.md §10.4); keys used here: name, scatter_regions"
    }
    backend: {
      description: "Backend where the workflow will be executed (HPC on the cluster)"
    }
    zones: {
      description: "Zones where compute will take place; required if backend is set to 'AWS' or 'GCP'"
    }
    cpuPlatform: {
      description: "CPU platform to use; optional if backend is set to 'GCP'"
    }
    gpuType: {
      description: "GPU type to use; optional if backend is set to 'GCP'"
    }
    container_registry: {
      description: "Container registry where workflow images are hosted"
    }
    preemptible: {
      description: "Where possible, run tasks preemptibly (no effect on HPC)"
    }
  }

  input {
    String ugc_wgw_version
    String cohort_id
    Array[String] sample_ids
    Array[String] sample_sexes
    Array[String] sample_sex_sources = []

    File? small_variant_vcf
    File? small_variant_vcf_index
    File sv_vcf
    File sv_vcf_index

    Int freq_threads = 2
    Int freq_mem_gb = 4

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

  #@ except: DeclarationName
  Map[String, String] ref_map = read_map(ref_map_file)

  String out_prefix = "~{cohort_id}.freq.~{ref_map["name"]}"

  # Provenance of the frequency files, written into their headers
  Array[String] header_common = [
    "##ugc_wgw_cohort_freq_version=~{ugc_wgw_version}",
    "##ugc_wgw_cohort_freq_cohort=~{cohort_id}",
    "##ugc_wgw_cohort_freq_strata=XX,XY (sheet sex, else inferred sex; members of unknown sex are counted in the unstratified tags only)",
    "##ugc_wgw_cohort_freq_sex_chromosomes=genotypes are diploid on every contig, so hemizygous ALT alleles of XY samples count under AC_Hom_XY",
    "##ugc_wgw_cohort_freq_half_calls=a half-called genotype (./1) contributes its called allele only, counted under AC_Hemi with AN += 1 (GLnexus emits ./1 when one allele lacks depth; sawfish encodes depth-based CNV carriers as ./1)"
  ]

  call Regions.ugc_wgw_regions_prepare { input:
    scatter_regions = ref_map["scatter_regions"],  # !FileCoercion
    require_whole_contigs = false,
    runtime_attributes = default_runtime_attributes
  }

  Int n_shards = length(ugc_wgw_regions_prepare.shard_names)

  # ---- small variants (GLnexus joint VCF) -----------------------------------

  if (defined(small_variant_vcf)) {
    Array[String] header_small_variants = flatten([header_common, [
      "##ugc_wgw_cohort_freq_source=~{basename(select_first([small_variant_vcf]))}",
      "##ugc_wgw_cohort_freq_counting=GLnexus AF, AC and AN removed and recounted by bcftools +fill-tags; missing genotypes (./.) are excluded from AN"
    ]])

    scatter (shard_index in range(n_shards)) {
      call Freq.ugc_wgw_bcftools_freq as small_variant_shard { input:
        vcf = select_first([small_variant_vcf]),
        vcf_index = select_first([small_variant_vcf_index]),
        regions_bed = ugc_wgw_regions_prepare.shard_beds[shard_index],
        sample_ids = sample_ids,
        sample_sexes = sample_sexes,
        header_lines = header_small_variants,
        remove_tags = "INFO/AF,INFO/AC,INFO/AN",
        out_prefix = "~{out_prefix}.small_variants.~{ugc_wgw_regions_prepare.shard_names[shard_index]}",
        threads = freq_threads,
        mem_gb = freq_mem_gb,
        runtime_attributes = default_runtime_attributes
      }
    }

    call Concat.ugc_wgw_bcftools_concat as concat_small_variants { input:
      vcfs = small_variant_shard.shard_vcf,
      vcf_indices = small_variant_shard.shard_vcf_index,
      out_prefix = "~{out_prefix}.small_variants",
      runtime_attributes = default_runtime_attributes
    }
  }

  # ---- structural variants (svx merge) ---------------------------------------

  Array[String] header_svs = flatten([header_common, [
    "##ugc_wgw_cohort_freq_source=~{basename(sv_vcf)}",
    "##ugc_wgw_cohort_freq_missing_gt=./. counted as hom-ref before tagging (the sample's own sawfish call set has no such variant); SUPP = carrier samples, CF = SUPP / N_SAMPLES"
  ]])

  scatter (shard_index in range(n_shards)) {
    call Freq.ugc_wgw_bcftools_freq as sv_shard { input:
      vcf = sv_vcf,
      vcf_index = sv_vcf_index,
      regions_bed = ugc_wgw_regions_prepare.shard_beds[shard_index],
      sample_ids = sample_ids,
      sample_sexes = sample_sexes,
      header_lines = header_svs,
      missing_to_homref = true,
      remove_tags = "INFO/IDLIST",
      extra_tags = "CF:1=INFO/SUPP/N_SAMPLES",
      out_prefix = "~{out_prefix}.structural_variants.~{ugc_wgw_regions_prepare.shard_names[shard_index]}",
      threads = freq_threads,
      mem_gb = freq_mem_gb,
      runtime_attributes = default_runtime_attributes
    }
  }

  call Concat.ugc_wgw_bcftools_concat as concat_svs { input:
    vcfs = sv_shard.shard_vcf,
    vcf_indices = sv_shard.shard_vcf_index,
    out_prefix = "~{out_prefix}.structural_variants",
    runtime_attributes = default_runtime_attributes
  }

  # ---- summary and provenance ------------------------------------------------

  call Freq.ugc_wgw_freq_summary { input:
    small_variant_vcf = concat_small_variants.vcf,
    small_variant_vcf_index = concat_small_variants.vcf_index,
    sv_vcf = concat_svs.vcf,
    sv_vcf_index = concat_svs.vcf_index,
    small_variant_samples = if defined(small_variant_vcf) then select_first([small_variant_shard.vcf_samples])[0] else small_variant_vcf,
    sv_samples = sv_shard.vcf_samples[0],
    sample_ids = sample_ids,
    sample_sexes = sample_sexes,
    sample_sex_sources = sample_sex_sources,
    out_prefix = out_prefix,
    runtime_attributes = default_runtime_attributes
  }

  call Provenance.ugc_wgw_manifest_write { input:
    ugc_wgw_version = ugc_wgw_version,
    stage = "cohort_freq",
    subject_type = "cohort",
    subject_id = cohort_id,
    member_ids = sample_ids,
    runtime_attributes = default_runtime_attributes
  }

  output {
    File? small_variant_freq_vcf = concat_small_variants.vcf
    File? small_variant_freq_vcf_index = concat_small_variants.vcf_index
    File sv_freq_vcf = concat_svs.vcf
    File sv_freq_vcf_index = concat_svs.vcf_index
    File freq_summary = ugc_wgw_freq_summary.summary
    File freq_samples = ugc_wgw_freq_summary.samples
    String ugc_wgw_workflow_version = ugc_wgw_version
    File ugc_wgw_manifest = ugc_wgw_manifest_write.ugc_wgw_manifest
  }
}
