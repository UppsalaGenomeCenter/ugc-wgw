version 1.0

# ugc-pacbio-wgw — stage `intake`: a sample registered from an aligned BAM and pre-called VCFs (the vcf entry kind) gets what `downstream` needs from `upstream`, without alignment or calling
# origin: none
# see docs/DESIGN.md §5, §8.1, §18 (2026-11-06); docs/ENTRYPOINTS.md §11

import "../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/tasks/mosdepth.wdl" as Mosdepth
import "../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/workflows/backend_configuration/backend_configuration.wdl" as BackendConfiguration
import "ugc_wgw/intake/checks.wdl" as Checks
import "ugc_wgw/qc/somalier.wdl" as Somalier
import "ugc_wgw/provenance.wdl" as Provenance

workflow ugc_wgw_intake {
  meta {
    description: "ugc-pacbio-wgw stage intake (per sample, the vcf entry kind): the input checks against the installed reference (contigs, index, one sample per VCF, methylation tags), mosdepth's depth and sex inference, the somalier extract and the sex resolution of upstream (declared, else coverage, else somalier), and the given VCFs with their sample column named after the sample id. Its outputs carry upstream's names for the fields downstream reads, so the driver runs downstream next as after upstream. Minutes."
    outputs: {
      aligned_hifi_reads: {
        description: "The given aligned BAM, passed through"
      },
      aligned_hifi_reads_index: {
        description: "Its index"
      },
      mosdepth_summary: {
        description: "mosdepth summary"
      },
      mosdepth_region_bed: {
        description: "mosdepth regions BED"
      },
      mosdepth_region_bed_index: {
        description: "Index"
      },
      mosdepth_depth_distribution_plot: {
        description: "mosdepth depth distribution plot"
      },
      inferred_sex: {
        description: "Sex inferred from coverage (mosdepth)"
      },
      stat_depth_mean: {
        description: "Mean depth"
      },
      resolved_sex: {
        description: "The sex the sex-dependent steps use: declared, else coverage, else somalier (empty when none could sex the sample)"
      },
      resolved_sex_source: {
        description: "declared, coverage, somalier or none"
      },
      somalier_sex: {
        description: "somalier's chrX call (MALE, FEMALE or empty)"
      },
      sex_tsv: {
        description: "The sex resolution table (every source and the decision)"
      },
      somalier_extract: {
        description: "The somalier extract (cohort_qc reads it)"
      },
      small_variant_vcf: {
        description: "The given small-variant VCF with its sample column named after the sample id"
      },
      small_variant_vcf_index: {
        description: "Index"
      },
      sv_vcf: {
        description: "The given SV VCF, sample column renamed"
      },
      sv_vcf_index: {
        description: "Index"
      },
      small_variant_gvcf: {
        description: "The given gVCF, sample column renamed (when given: joint mode needs it)"
      },
      small_variant_gvcf_index: {
        description: "Index"
      },
      input_checks: {
        description: "The checks table: check, target, status, detail"
      },
      input_provenance: {
        description: "The BAM header's @PG/@RG fields and the VCF headers' source lines"
      },
      msg: {
        description: "Messages (the checks' warnings, the sex resolution's)"
      },
      ugc_wgw_workflow_version: {
        description: "The ugc-wgw version this run was started with (the engine records nothing else about the caller)"
      },
      ugc_wgw_manifest: {
        description: "Run manifest (stage, subject, version)"
      }
    }
  }

  parameter_meta {
    ugc_wgw_version: {
      description: "ugc-wgw version (driver-supplied; recorded in the manifest)"
    }
    sample_id: {
      description: "Sample ID"
    }
    declared_sex: {
      description: "The sample sheet's sex (MALE or FEMALE), when it has one: the first source of the resolved sex"
    }
    aligned_bam: {
      description: "The aligned BAM from the sample sheet (aligned_bam), to the installed reference"
    }
    aligned_bam_index: {
      description: "Its index"
    }
    small_variant_vcf_in: {
      description: "The small-variant VCF from the sheet (one sample)"
    }
    small_variant_vcf_in_index: {
      description: "Its index"
    }
    sv_vcf_in: {
      description: "The SV VCF from the sheet (one sample)"
    }
    sv_vcf_in_index: {
      description: "Its index"
    }
    gvcf_in: {
      description: "The gVCF from the sheet, when given (joint mode needs it)"
    }
    gvcf_in_index: {
      description: "Its index"
    }
    require_methylation_tags: {
      description: "Fail the intake when the BAM's first reads carry no MM/ML tags (true); false notes it instead"
    }
    somalier_sex_min_x_sites: {
      description: "Minimum chrX sites for somalier's sex call to count"
    }
    somalier_sex_min_x_depth: {
      description: "Minimum chrX depth for somalier's sex call to count"
    }
    ref_map_file: {
      description: "Reference map (ugc-wgw install): name, fasta, fasta_index, max_norm_female_chrY_depth, somalier_sites_vcf"
    }
    ugc_wgw_container_registry: {
      description: "Registry of the images ugc-wgw builds (somalier)"
    }
    backend: {
      description: "Backend where the workflow will be executed; HPC on the cluster"
    }
    zones: {
      description: "Zones where compute will take place; required if backend is set to 'AWS' or 'GCP'"
    }
    cpuPlatform: {
      description: "CPU platform; optional for GCP"
    }
    gpuType: {
      description: "GPU type to use; required if gpu is set to true for cloud backends; must match the hardware available"
    }
    container_registry: {
      description: "Container registry where workflow images are hosted (the installed SIF cache on the HPC)"
    }
    preemptible: {
      description: "Where possible, run tasks preemptibly"
    }
  }

  input {
    String ugc_wgw_version
    String sample_id
    String? declared_sex
    File aligned_bam
    File aligned_bam_index
    File small_variant_vcf_in
    File small_variant_vcf_in_index
    File sv_vcf_in
    File sv_vcf_in_index
    File? gvcf_in
    File? gvcf_in_index
    Boolean require_methylation_tags = true
    Int somalier_sex_min_x_sites = 50
    Float somalier_sex_min_x_depth = 5.0

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

  Map[String, String] ref_map = read_map(ref_map_file)

  String ref_name = ref_map["name"]
  File ref_fasta = ref_map["fasta"]  # !FileCoercion
  File ref_index = ref_map["fasta_index"]  # !FileCoercion
  Float max_norm_female_chrY_depth = ref_map["max_norm_female_chrY_depth"]  # !StringCoercion
  File somalier_sites_vcf = ref_map["somalier_sites_vcf"]  # !FileCoercion

  # ---- the checks first: a failure stops the stage before any other job ----

  call Checks.ugc_wgw_input_checks { input:
    sample_id = sample_id,
    aligned_bam = aligned_bam,
    aligned_bam_index = aligned_bam_index,
    small_variant_vcf = small_variant_vcf_in,
    small_variant_vcf_index = small_variant_vcf_in_index,
    sv_vcf = sv_vcf_in,
    sv_vcf_index = sv_vcf_in_index,
    gvcf = gvcf_in,
    gvcf_index = gvcf_in_index,
    ref_index = ref_index,
    require_methylation_tags = require_methylation_tags,
    runtime_attributes = default_runtime_attributes
  }

  # ---- what upstream would have made: depth and sex ----

  call Mosdepth.mosdepth { input:
    sample_id = sample_id,
    ref_name = ref_name,
    aligned_bam = aligned_bam,
    aligned_bam_index = aligned_bam_index,
    infer_sex = true,
    max_norm_female_chrY_depth = max_norm_female_chrY_depth,
    runtime_attributes = default_runtime_attributes
  }

  call Somalier.ugc_wgw_somalier_extract { input:
    sample_id = sample_id,
    aligned_bam = aligned_bam,
    aligned_bam_index = aligned_bam_index,
    sites_vcf = somalier_sites_vcf,
    ref_fasta = ref_fasta,
    ref_index = ref_index,
    ugc_wgw_container_registry = ugc_wgw_container_registry,
    runtime_attributes = default_runtime_attributes
  }

  call Somalier.ugc_wgw_somalier_relate_one { input:
    sample_id = sample_id,
    somalier = ugc_wgw_somalier_extract.somalier,
    sites_vcf = somalier_sites_vcf,
    ugc_wgw_container_registry = ugc_wgw_container_registry,
    runtime_attributes = default_runtime_attributes
  }

  call Somalier.ugc_wgw_resolve_sex { input:
    sample_id = sample_id,
    declared_sex = declared_sex,
    coverage_sex = mosdepth.inferred_sex,
    somalier_samples_tsv = ugc_wgw_somalier_relate_one.samples_tsv,
    min_x_sites = somalier_sex_min_x_sites,
    min_x_depth = somalier_sex_min_x_depth,
    runtime_attributes = default_runtime_attributes
  }

  # ---- the VCFs under the sample id and upstream's names ----

  call Checks.ugc_wgw_vcf_sample_rename as rename_small { input:
    vcf = small_variant_vcf_in,
    vcf_index = small_variant_vcf_in_index,
    sample_id = sample_id,
    out_name = "~{sample_id}.~{ref_name}.small_variants.vcf.gz",
    runtime_attributes = default_runtime_attributes
  }

  call Checks.ugc_wgw_vcf_sample_rename as rename_sv { input:
    vcf = sv_vcf_in,
    vcf_index = sv_vcf_in_index,
    sample_id = sample_id,
    out_name = "~{sample_id}.~{ref_name}.structural_variants.vcf.gz",
    runtime_attributes = default_runtime_attributes
  }

  if (defined(gvcf_in)) {
    call Checks.ugc_wgw_vcf_sample_rename as rename_gvcf { input:
      vcf = select_first([gvcf_in]),
      vcf_index = select_first([gvcf_in_index]),
      sample_id = sample_id,
      out_name = "~{sample_id}.~{ref_name}.small_variants.g.vcf.gz",
      runtime_attributes = default_runtime_attributes
    }
  }

  call Provenance.ugc_wgw_manifest_write { input:
    ugc_wgw_version = ugc_wgw_version,
    stage = "intake",
    subject_type = "sample",
    subject_id = sample_id,
    runtime_attributes = default_runtime_attributes
  }

  output {
    File aligned_hifi_reads = aligned_bam
    File aligned_hifi_reads_index = aligned_bam_index
    File mosdepth_summary = mosdepth.summary
    File mosdepth_region_bed = mosdepth.region_bed
    File mosdepth_region_bed_index = mosdepth.region_bed_index
    File mosdepth_depth_distribution_plot = mosdepth.depth_distribution_plot
    String inferred_sex = mosdepth.inferred_sex
    String stat_depth_mean = mosdepth.stat_depth_mean
    String resolved_sex = ugc_wgw_resolve_sex.resolved_sex
    String resolved_sex_source = ugc_wgw_resolve_sex.resolved_sex_source
    String somalier_sex = ugc_wgw_resolve_sex.somalier_sex
    File sex_tsv = ugc_wgw_resolve_sex.sex_tsv
    File somalier_extract = ugc_wgw_somalier_extract.somalier
    File small_variant_vcf = rename_small.vcf_out
    File small_variant_vcf_index = rename_small.vcf_out_index
    File sv_vcf = rename_sv.vcf_out
    File sv_vcf_index = rename_sv.vcf_out_index
    File? small_variant_gvcf = rename_gvcf.vcf_out
    File? small_variant_gvcf_index = rename_gvcf.vcf_out_index
    File input_checks = ugc_wgw_input_checks.checks
    File input_provenance = ugc_wgw_input_checks.provenance
    Array[String] msg = flatten([ugc_wgw_input_checks.msg, ugc_wgw_resolve_sex.msg])
    String ugc_wgw_workflow_version = ugc_wgw_version
    File ugc_wgw_manifest = ugc_wgw_manifest_write.ugc_wgw_manifest
  }
}
