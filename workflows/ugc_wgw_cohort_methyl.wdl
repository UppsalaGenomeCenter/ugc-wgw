version 1.0

# ugc-pacbio-wgw — stage `cohort_methyl`: the cohort's methylation layer from the members' MethBat outputs (background, scores, comparisons, imprinting, optional joint segmentation and deconvolution)
# origin: none
# see docs/DESIGN.md §5, §8.1, §18 (2026-11-06); docs/ENTRYPOINTS.md §10

import "../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/workflows/backend_configuration/backend_configuration.wdl" as BackendConfiguration
import "ugc_wgw/methyl/report.wdl" as MethylReport
import "ugc_wgw/methyl/cohort.wdl" as Methyl
import "ugc_wgw/provenance.wdl" as Provenance

workflow ugc_wgw_cohort_methyl {
  meta {
    description: "ugc-pacbio-wgw stage cohort_methyl: a MethBat background profile of the cohort (every member, or N leave-one-out backgrounds), each member's 5mC pileup scored against it under the labels it carries (ALL, its sex group, the sheet columns the project names), group comparisons (XX against XY by default), the imprinting reports tabulated (made here for a member whose run predates the per-sample task), optionally the cohort's joint segmentation and per-member cell-type deconvolution, and the tables the driver reads: the long summary, the per-region table and the events file with every region's genes and type from the CpG-island map. Reads only per-sample outputs; last in the standalone and joint modes."
    outputs: {
      cohort_background: {
        description: "The cohort background profile (every member)"
      },
      cohort_backgrounds_loo: {
        description: "With leave_one_out: one background per member without that member, cohort order"
      },
      scored_profiles: {
        description: "Per member, cohort order: the scored profile as MethBat writes it"
      },
      scored_profiles_annotated: {
        description: "Per member: the scored profile with hgnc_id, gene_symbol and region_type appended"
      },
      comparisons: {
        description: "One comparison file per compare_pairs entry"
      },
      joint_segments: {
        description: "With run_joint_segment: the cohort's segmented methylome (BED)"
      },
      joint_combined_methyl: {
        description: "With run_joint_segment: cohort-averaged combined methylation per segment (bedGraph)"
      },
      joint_asm: {
        description: "With run_joint_segment: the haplotype spread per segment (bedGraph)"
      },
      deconvolution: {
        description: "With run_deconvolve: one estimate JSON per member, cohort order"
      },
      methyl_deconvolution: {
        description: "<prefix>.deconvolution.tsv: the estimates as one table, a row per member and a column per cell type of the atlas, with the fit (residual, raw fraction sum); header only without run_deconvolve"
      },
      methyl_summary: {
        description: "Long-format summary: resource, metric, key, value"
      },
      methyl_regions: {
        description: "Per region and label: dataset counts, event counts across members, the gene columns"
      },
      methyl_events: {
        description: "Every member's Hypo/Hyper rows with the gene columns, sorted by gene then member"
      },
      methyl_imprinting: {
        description: "Member x imprinting region from methbat report, with report_source"
      },
      methyl_samples: {
        description: "One row per member: labels, event counts, summary-label counts, imprinting counts"
      },
      methyl_collection: {
        description: "The collection file the background was built from (identifier, filename, labels)"
      },
      member_reports: {
        description: "The methbat report used per member, cohort order (the run's or the one made here)"
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
      description: "Members, cohort order"
    }
    sample_labels: {
      description: "Per member, aligned: the labels for the background, ';'-joined (the sex group as cohort_freq derives it, then the values of the sheet columns named in methyl_label_columns; may be empty). ALL is implicit."
    }
    cpg_pileup_beds: {
      description: "Per member, aligned: the 5mCpG pileup BED (singleton or downstream)"
    }
    methbat_profiles: {
      description: "Per member, aligned: the MethBat profile (singleton or downstream)"
    }
    methbat_reports: {
      description: "Per member, aligned: the methbat report of its run, null for a run from before the task (the report is then made here from the pileup)"
    }
    compare_pairs: {
      description: "Label pairs to compare, baseline:compare (default XX:XY); a side absent from the background gives a header-only file and a note"
    }
    leave_one_out: {
      description: "Score each member against a background of the other members (N builds) instead of the whole cohort"
    }
    background_min_datasets: {
      description: "Minimum members in a label group for a Hypo/Hyper call and for a comparison (MethBat's --min-datasets / --min-samples)"
    }
    zscore_min: {
      description: "Minimum absolute z-score for a Hypo/Hyper call (methbat profile --min-abs-zscore; compare --min-zscore)"
    }
    methyl_delta_min: {
      description: "Minimum absolute methylation delta in percentage points (methbat profile --min-abs-delta; compare --min-delta)"
    }
    run_joint_segment: {
      description: "Segment the cohort's methylome jointly from every pileup (one whole-genome task; its memory grows with the cohort)"
    }
    run_deconvolve: {
      description: "Estimate cell-type proportions per member against the shipped atlas (methbat deconvolve, nanomix's 39 cell types from bisulfite data; a minute and 2.5 GB per member). On by default since 0.20.0: blood samples resolve into granulocytes, T cells and monocytes as expected, minor fractions are noise. false skips it"
    }
    build_mem_gb: {
      description: "Memory of methbat build (GB)"
    }
    profile_mem_gb: {
      description: "Memory of each scoring task (GB)"
    }
    joint_segment_threads: {
      description: "Threads of methbat joint-segment"
    }
    joint_segment_mem_gb: {
      description: "Memory of methbat joint-segment (GB)"
    }
    summary_mem_gb: {
      description: "Memory of the summary task (GB)"
    }
    ref_map_file: {
      description: "Reference map (ugc-wgw install): methbat_map, methbat_report_regions_tsv, methbat_atlas_tsv, name"
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
    Array[String] sample_labels
    Array[File] cpg_pileup_beds
    Array[File] methbat_profiles
    Array[File?] methbat_reports

    Array[String] compare_pairs = ["XX:XY"]
    Boolean leave_one_out = false
    Int background_min_datasets = 10
    Float zscore_min = 3.0
    Float methyl_delta_min = 20
    Boolean run_joint_segment = false
    Boolean run_deconvolve = true
    Int build_mem_gb = 8
    Int profile_mem_gb = 4
    Int joint_segment_threads = 8
    Int joint_segment_mem_gb = 32
    Int summary_mem_gb = 4

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
  String out_prefix = "~{cohort_id}.methyl.~{ref_name}"
  File methbat_map = ref_map["methbat_map"]  # !FileCoercion
  File report_regions_tsv = ref_map["methbat_report_regions_tsv"]  # !FileCoercion

  # ---- the background of every member ----

  call Methyl.ugc_wgw_methbat_build { input:
    sample_ids = sample_ids,
    profiles = methbat_profiles,
    labels = sample_labels,
    out_prefix = out_prefix,
    mem_gb = build_mem_gb,
    runtime_attributes = default_runtime_attributes
  }

  if (leave_one_out) {
    scatter (j in range(length(sample_ids))) {
      call Methyl.ugc_wgw_methbat_build as build_loo { input:
        sample_ids = sample_ids,
        profiles = methbat_profiles,
        labels = sample_labels,
        exclude_index = j,
        out_prefix = "~{out_prefix}.loo.~{sample_ids[j]}",
        mem_gb = build_mem_gb,
        runtime_attributes = default_runtime_attributes
      }
    }
  }

  # ---- per member: the imprinting report (the run's, or made here), the scoring, the deconvolution ----

  scatter (i in range(length(sample_ids))) {
    if (!defined(methbat_reports[i])) {
      call MethylReport.ugc_wgw_methbat_report { input:
        cpg_pileup_bed = cpg_pileup_beds[i],
        report_regions_tsv = report_regions_tsv,
        out_prefix = "~{sample_ids[i]}.~{ref_name}",
        runtime_attributes = default_runtime_attributes
      }
    }
    File member_report = select_first([methbat_reports[i], ugc_wgw_methbat_report.report])
    String report_source = if defined(methbat_reports[i]) then "run" else "stage"
    File background_full = ugc_wgw_methbat_build.background
  }

  Array[File] member_backgrounds = if leave_one_out then select_first([build_loo.background]) else background_full

  scatter (i in range(length(sample_ids))) {
    call Methyl.ugc_wgw_methbat_profile_bg { input:
      sample_id = sample_ids[i],
      cpg_pileup_bed = cpg_pileup_beds[i],
      background = member_backgrounds[i],
      labels = sample_labels[i],
      methbat_map = methbat_map,
      zscore_min = zscore_min,
      methyl_delta_min = methyl_delta_min,
      background_min_datasets = background_min_datasets,
      out_prefix = "~{sample_ids[i]}.~{ref_name}",
      mem_gb = profile_mem_gb,
      runtime_attributes = default_runtime_attributes
    }

    if (run_deconvolve) {
      call Methyl.ugc_wgw_methbat_deconvolve { input:
        sample_id = sample_ids[i],
        cpg_pileup_bed = cpg_pileup_beds[i],
        atlas_tsv = ref_map["methbat_atlas_tsv"],  # !FileCoercion
        out_prefix = "~{sample_ids[i]}.~{ref_name}",
        runtime_attributes = default_runtime_attributes
      }
    }
  }

  # ---- comparisons and the joint segmentation ----

  scatter (pair in compare_pairs) {
    call Methyl.ugc_wgw_methbat_compare { input:
      background = ugc_wgw_methbat_build.background,
      pair = pair,
      zscore_min = zscore_min,
      methyl_delta_min = methyl_delta_min,
      min_samples = background_min_datasets,
      out_prefix = out_prefix,
      runtime_attributes = default_runtime_attributes
    }
  }

  if (run_joint_segment) {
    call Methyl.ugc_wgw_methbat_joint_segment { input:
      sample_ids = sample_ids,
      cpg_pileup_beds = cpg_pileup_beds,
      out_prefix = out_prefix,
      threads = joint_segment_threads,
      mem_gb = joint_segment_mem_gb,
      runtime_attributes = default_runtime_attributes
    }
  }

  # ---- the tables the driver reads, and provenance ----

  call Methyl.ugc_wgw_methyl_summary { input:
    sample_ids = sample_ids,
    labels = sample_labels,
    events = ugc_wgw_methbat_profile_bg.events,
    counts = ugc_wgw_methbat_profile_bg.counts,
    reports = member_report,
    report_sources = report_source,
    comparisons = ugc_wgw_methbat_compare.comparison,
    compare_notes = ugc_wgw_methbat_compare.note,
    background = ugc_wgw_methbat_build.background,
    methbat_map = methbat_map,
    deconvolution = select_all(ugc_wgw_methbat_deconvolve.estimate),
    background_min_datasets = background_min_datasets,
    out_prefix = out_prefix,
    mem_gb = summary_mem_gb,
    runtime_attributes = default_runtime_attributes
  }

  call Provenance.ugc_wgw_manifest_write { input:
    ugc_wgw_version = ugc_wgw_version,
    stage = "cohort_methyl",
    subject_type = "cohort",
    subject_id = cohort_id,
    member_ids = sample_ids,
    runtime_attributes = default_runtime_attributes
  }

  output {
    File cohort_background = ugc_wgw_methbat_build.background
    Array[File]? cohort_backgrounds_loo = build_loo.background
    Array[File] scored_profiles = ugc_wgw_methbat_profile_bg.scored
    Array[File] scored_profiles_annotated = ugc_wgw_methbat_profile_bg.scored_annotated
    Array[File] comparisons = ugc_wgw_methbat_compare.comparison
    File? joint_segments = ugc_wgw_methbat_joint_segment.regions_bed
    File? joint_combined_methyl = ugc_wgw_methbat_joint_segment.combined_methyl_bedgraph
    File? joint_asm = ugc_wgw_methbat_joint_segment.asm_bedgraph
    Array[File?] deconvolution = ugc_wgw_methbat_deconvolve.estimate
    File methyl_summary = ugc_wgw_methyl_summary.summary
    File methyl_regions = ugc_wgw_methyl_summary.regions
    File methyl_events = ugc_wgw_methyl_summary.cohort_events
    File methyl_imprinting = ugc_wgw_methyl_summary.imprinting
    File methyl_samples = ugc_wgw_methyl_summary.samples
    File methyl_deconvolution = ugc_wgw_methyl_summary.deconvolution_table
    File methyl_collection = ugc_wgw_methbat_build.collection
    Array[File] member_reports = member_report
    String ugc_wgw_workflow_version = ugc_wgw_version
    File ugc_wgw_manifest = ugc_wgw_manifest_write.ugc_wgw_manifest
  }
}
