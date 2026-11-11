version 1.0

# ugc-pacbio-wgw — stage `cohort_tr`: the cohort's tandem repeats from cohort_merge's TRGT VCF (the disease-locus screen against STRchive's thresholds, the genome-wide expansion-outlier scan, the allele-length resource)
# origin: none
# see docs/DESIGN.md §5, §7.6, §18 (2026-11-11); docs/ENTRYPOINTS.md §10

import "../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/workflows/backend_configuration/backend_configuration.wdl" as BackendConfiguration
import "ugc_wgw/cohort/regions.wdl" as Regions
import "ugc_wgw/tr/cohort.wdl" as Tr
import "ugc_wgw/provenance.wdl" as Provenance

workflow ugc_wgw_cohort_tr {
  meta {
    description: "ugc-pacbio-wgw stage cohort_tr: what a cohort makes answerable about tandem repeats, from the merged TRGT VCF and the trgt-lps table of cohort_merge alone. The 73 named disease loci of the catalog are screened per member against STRchive's thresholds (copies of the pathogenic motif per allele, classes, inheritance fit, methylation and purity carried along); every locus of the catalog gets its allele-length distribution and spectrum, and the members whose longer allele is far above the cohort's and rare in it are listed as expansion outliers. Reads no per-sample file; optional in the planner (cohort_methyl never waits for it); no de novo calling (trgt-denovo is not part of this stage)."
    outputs: {
      tr_pathogenic_screen: {
        description: "One row per member and named locus: copies and class per allele, inheritance fit, notes (a research screen, not a diagnostic result)"
      },
      tr_screen_loci: {
        description: "One row per named locus: thresholds, call rate, the copies distribution, members per class"
      },
      tr_locus_summary: {
        description: "Every locus of the catalog: called members, allele-length median, MAD and quantiles, carriers, outlier status (bgzipped)"
      },
      tr_allele_spectrum: {
        description: "trid, al, n_alleles, n_members: the cohort's allele-length resource (bgzipped)"
      },
      tr_outlier_pairs: {
        description: "Every member and locus pair listed as an expansion outlier (bgzipped)"
      },
      tr_outliers_top: {
        description: "The strongest outliers per member"
      },
      tr_samples: {
        description: "One row per member: loci called, outliers, the screen's classes"
      },
      tr_summary: {
        description: "Long-format summary: resource, metric, key, value"
      },
      ugc_wgw_workflow_version: {
        description: "The ugc-wgw version this run was started with (the engine records nothing else about the caller)"
      },
      ugc_wgw_manifest: {
        description: "Run manifest (stage, subject, version, member list)"
      }
    }
  }

  parameter_meta {
    ugc_wgw_version: {
      description: "ugc-wgw version (driver-supplied; recorded in the manifest)"
    }
    cohort_id: {
      description: "Cohort identifier (file prefix)"
    }
    sample_ids: {
      description: "Members, cohort order; every one must be a sample of the TRGT VCF"
    }
    sample_sexes: {
      description: "Per member, aligned: XX, XY or empty (the sex group as cohort_freq derives it); decides how many members a chrY locus expects and the ploidy check on chrX"
    }
    trgt_vcf: {
      description: "The cohort's merged TRGT VCF (cohort_merge's cohort_trgt_vcf)"
    }
    trgt_vcf_index: {
      description: "Index for trgt_vcf"
    }
    trgt_lps: {
      description: "The cohort's trgt-lps table (cohort_merge's cohort_trgt_lps)"
    }
    outlier_min_bp: {
      description: "An outlier's longer allele is at least this many bp above the cohort median"
    }
    outlier_mad_min: {
      description: "... and at least this many scaled MADs above it (the scale never below one repeat unit)"
    }
    outlier_min_called: {
      description: "Loci with fewer called members are not scanned"
    }
    outlier_max_carrier_fraction: {
      description: "A locus where more than this fraction of the called members qualify is common and lists nobody (raise it for a cohort of relatives)"
    }
    outlier_max_carriers_min: {
      description: "... but this many carriers are always allowed"
    }
    screen_borderline_bp: {
      description: "A screened allele is borderline when its class changes within this many bp of sequence"
    }
    top_outliers_per_member: {
      description: "Rows per member in tr_outliers_top (what the pages read)"
    }
    tr_disease_loci_override: {
      description: "A named-loci table to use instead of the reference map's (a newer STRchive release made with scripts/make-strchive-map.py)"
    }
    stats_mem_gb: {
      description: "Memory of the per-shard statistics task (GB)"
    }
    screen_mem_gb: {
      description: "Memory of the screen task (GB)"
    }
    gather_mem_gb: {
      description: "Memory of the gather task (GB)"
    }
    ref_map_file: {
      description: "Reference map (ugc-wgw install): name, scatter_regions, tr_disease_loci_tsv, tr_named_overlaps_tsv"
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
    String cohort_id
    Array[String] sample_ids
    Array[String] sample_sexes
    File trgt_vcf
    File trgt_vcf_index
    File trgt_lps

    Int outlier_min_bp = 50
    Float outlier_mad_min = 5.0
    Int outlier_min_called = 20
    Float outlier_max_carrier_fraction = 0.01
    Int outlier_max_carriers_min = 2
    Int screen_borderline_bp = 6
    Int top_outliers_per_member = 20
    File? tr_disease_loci_override
    Int stats_mem_gb = 2
    Int screen_mem_gb = 2
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
  String out_prefix = "~{cohort_id}.tr.~{ref_name}"
  File disease_loci_reference = ref_map["tr_disease_loci_tsv"]  # !FileCoercion
  File disease_loci = select_first([tr_disease_loci_override, disease_loci_reference])
  File named_overlaps = ref_map["tr_named_overlaps_tsv"]  # !FileCoercion

  # ---- the named disease loci: one task, 73 loci and their twins ----

  call Tr.ugc_wgw_tr_screen { input:
    trgt_vcf = trgt_vcf,
    trgt_vcf_index = trgt_vcf_index,
    trgt_lps = trgt_lps,
    disease_loci_tsv = disease_loci,
    named_overlaps = named_overlaps,
    sample_ids = sample_ids,
    sample_sexes = sample_sexes,
    borderline_bp = screen_borderline_bp,
    out_prefix = out_prefix,
    mem_gb = screen_mem_gb,
    runtime_attributes = default_runtime_attributes
  }

  # ---- every locus: per shard, a record belonging to the shard that holds its POS ----

  call Regions.ugc_wgw_regions_prepare { input:
    scatter_regions = ref_map["scatter_regions"],  # !FileCoercion
    require_whole_contigs = false,
    runtime_attributes = default_runtime_attributes
  }

  scatter (shard_index in range(length(ugc_wgw_regions_prepare.shard_names))) {
    call Tr.ugc_wgw_tr_locus_stats { input:
      trgt_vcf = trgt_vcf,
      trgt_vcf_index = trgt_vcf_index,
      regions_bed = ugc_wgw_regions_prepare.shard_beds[shard_index],
      sample_ids = sample_ids,
      sample_sexes = sample_sexes,
      named_overlaps = named_overlaps,
      outlier_min_bp = outlier_min_bp,
      outlier_mad_min = outlier_mad_min,
      outlier_min_called = outlier_min_called,
      outlier_max_carrier_fraction = outlier_max_carrier_fraction,
      outlier_max_carriers_min = outlier_max_carriers_min,
      out_prefix = "~{out_prefix}.~{ugc_wgw_regions_prepare.shard_names[shard_index]}",
      mem_gb = stats_mem_gb,
      runtime_attributes = default_runtime_attributes
    }
  }

  # ---- the tables the driver reads, and provenance ----

  call Tr.ugc_wgw_tr_gather { input:
    loci = ugc_wgw_tr_locus_stats.loci,
    spectra = ugc_wgw_tr_locus_stats.spectrum,
    outliers = ugc_wgw_tr_locus_stats.outliers,
    member_counts = ugc_wgw_tr_locus_stats.member_counts,
    screen = ugc_wgw_tr_screen.screen,
    screen_loci = ugc_wgw_tr_screen.screen_loci,
    sample_ids = sample_ids,
    sample_sexes = sample_sexes,
    disease_loci_tsv = disease_loci,
    top_outliers_per_member = top_outliers_per_member,
    parameters = "outlier_min_bp=~{outlier_min_bp};outlier_mad_min=~{outlier_mad_min};outlier_min_called=~{outlier_min_called};outlier_max_carrier_fraction=~{outlier_max_carrier_fraction};outlier_max_carriers_min=~{outlier_max_carriers_min};screen_borderline_bp=~{screen_borderline_bp}",
    out_prefix = out_prefix,
    mem_gb = gather_mem_gb,
    runtime_attributes = default_runtime_attributes
  }

  call Provenance.ugc_wgw_manifest_write { input:
    ugc_wgw_version = ugc_wgw_version,
    stage = "cohort_tr",
    subject_type = "cohort",
    subject_id = cohort_id,
    member_ids = sample_ids,
    runtime_attributes = default_runtime_attributes
  }

  output {
    File tr_pathogenic_screen = ugc_wgw_tr_screen.screen
    File tr_screen_loci = ugc_wgw_tr_screen.screen_loci
    File tr_locus_summary = ugc_wgw_tr_gather.loci_gz
    File tr_allele_spectrum = ugc_wgw_tr_gather.spectrum
    File tr_outlier_pairs = ugc_wgw_tr_gather.outliers_gz
    File tr_outliers_top = ugc_wgw_tr_gather.outliers_top
    File tr_samples = ugc_wgw_tr_gather.samples
    File tr_summary = ugc_wgw_tr_gather.summary
    String ugc_wgw_workflow_version = ugc_wgw_version
    File ugc_wgw_manifest = ugc_wgw_manifest_write.ugc_wgw_manifest
  }
}
