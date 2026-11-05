version 1.0

# ugc-pacbio-wgw — stage `cohort_qc`: identity and relatedness of the cohort's members (somalier), with an ancestry projection
# origin: none
# see docs/DESIGN.md §5, §8.1; docs/ENTRYPOINTS.md §9.1

import "../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/workflows/backend_configuration/backend_configuration.wdl" as BackendConfiguration
import "ugc_wgw/qc/somalier.wdl" as Somalier
import "ugc_wgw/provenance.wdl" as Provenance

workflow ugc_wgw_cohort_qc {
  meta {
    description: "ugc-pacbio-wgw stage cohort_qc: pairwise relatedness, each member's sex from chrX heterozygosity and a check of the sample sheet's pedigree, from one small somalier extract per member (made by upstream or singleton; made here from the aligned BAM for members whose run predates the extract), plus an ancestry projection onto the 1000 Genomes superpopulations when the reference set is installed. Minutes at any cohort size; the stage never blocks the stages after it, its flags feed `ugc-wgw summary --pass-list --pass-on identity`."
    outputs: {
      qc_pairs: {
        description: "somalier relate pairs table: every pair with relatedness, IBS0, IBS2, concordance and the pedigree's expected relatedness"
      },
      qc_samples: {
        description: "somalier relate samples table: depth, heterozygosity, chrX and chrY statistics per member"
      },
      qc_groups: {
        description: "somalier relate groups: members that look identical"
      },
      qc_html: {
        description: "somalier's interactive relatedness report"
      },
      qc_ped: {
        description: "The PED file made from the sheet"
      },
      qc_summary: {
        description: "Long-format summary (relatedness histogram, counts)"
      },
      qc_members: {
        description: "One row per member: closest relative, heterozygosity ratio, somalier's chrX call, ancestry, whether the extract was made here"
      },
      qc_related_pairs: {
        description: "The pairs with relatedness at or above 0.1 or an expectation from the pedigree"
      },
      qc_ancestry: {
        description: "somalier ancestry table (absent when run_ancestry is false)"
      },
      qc_ancestry_html: {
        description: "somalier ancestry report (absent when run_ancestry is false)"
      },
      somalier_extracts: {
        description: "The extract of every member, cohort order (the per-sample run's, or the one made here)"
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
      description: "Cohort members, cohort order"
    }
    somalier_files: {
      description: "The members' somalier extracts, aligned with sample_ids; null for a member whose per-sample run predates the extract (then its aligned BAM is given and the extract is made here)"
    }
    aligned_bams: {
      description: "Aligned BAM per member without an extract, null otherwise (aligned with sample_ids; empty when every member has an extract)"
    }
    aligned_bam_indices: {
      description: "Index per entry of aligned_bams"
    }
    sample_sexes: {
      description: "XX, XY or empty per member, aligned with sample_ids (the resolved sex of the member's run)"
    }
    sample_sex_sources: {
      description: "declared, coverage, somalier or empty per member"
    }
    father_ids: {
      description: "Father per member from the sheet, empty where none"
    }
    mother_ids: {
      description: "Mother per member from the sheet, empty where none"
    }
    run_ancestry: {
      description: "Project the members onto the 1000 Genomes PCA (needs the reference map's somalier_1kg_tar and somalier_1kg_labels; the driver sets it when the tar is installed)"
    }
    min_depth: {
      description: "somalier relate --min-depth"
    }
    min_ab: {
      description: "somalier relate --min-ab"
    }
    report_all_pairs: {
      description: "Report every pair in the pairs table (somalier subsamples unrelated pairs to 200 000 otherwise)"
    }
    relate_threads: {
      description: "Threads of the relate task"
    }
    relate_mem_gb: {
      description: "Memory (GB) of the relate task"
    }
    ancestry_mem_gb: {
      description: "Memory (GB) of the ancestry task (a PCA of 2 504 + N extracts: about 1 GB at N = 4)"
    }
    extract_mem_gb: {
      description: "Memory (GB) of an extract made here"
    }
    ref_map_file: {
      description: "ugc-wgw reference map (docs/DESIGN.md §10.4); keys used here: name, fasta, fasta_index, somalier_sites_vcf, and with run_ancestry somalier_1kg_tar, somalier_1kg_labels"
    }
    ugc_wgw_container_registry: {
      description: "Registry holding the ugc-built images (somalier)"
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
    Array[File?] somalier_files
    Array[File?] aligned_bams = []
    Array[File?] aligned_bam_indices = []
    Array[String] sample_sexes
    Array[String] sample_sex_sources = []
    Array[String] father_ids
    Array[String] mother_ids

    Boolean run_ancestry = false
    Int min_depth = 7
    Float min_ab = 0.3
    Boolean report_all_pairs = true
    Int relate_threads = 2
    Int relate_mem_gb = 4
    Int ancestry_mem_gb = 8
    Int extract_mem_gb = 2

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

  String out_prefix = "~{cohort_id}.qc.~{ref_map["name"]}"
  File sites_vcf = ref_map["somalier_sites_vcf"]  # !FileCoercion

  # ---- one extract per member: the per-sample run's, or made here from the aligned BAM ----

  scatter (i in range(length(sample_ids))) {
    if (!defined(somalier_files[i])) {
      call Somalier.ugc_wgw_somalier_extract { input:
        sample_id = sample_ids[i],
        aligned_bam = select_first([aligned_bams[i]]),
        aligned_bam_index = select_first([aligned_bam_indices[i]]),
        sites_vcf = sites_vcf,
        ref_fasta = ref_map["fasta"],  # !FileCoercion
        ref_index = ref_map["fasta_index"],  # !FileCoercion
        mem_gb = extract_mem_gb,
        ugc_wgw_container_registry = ugc_wgw_container_registry,
        runtime_attributes = default_runtime_attributes
      }
    }
    File member_somalier = select_first([somalier_files[i], ugc_wgw_somalier_extract.somalier])
    String legacy_flag = if defined(somalier_files[i]) then "" else sample_ids[i]
  }

  # ---- pedigree and relatedness ----

  call Somalier.ugc_wgw_somalier_ped { input:
    cohort_id = cohort_id,
    sample_ids = sample_ids,
    sample_sexes = sample_sexes,
    father_ids = father_ids,
    mother_ids = mother_ids,
    runtime_attributes = default_runtime_attributes
  }

  call Somalier.ugc_wgw_somalier_relate { input:
    out_prefix = out_prefix,
    somalier_files = member_somalier,
    ped = ugc_wgw_somalier_ped.ped,
    sites_vcf = sites_vcf,
    min_depth = min_depth,
    min_ab = min_ab,
    report_all_pairs = report_all_pairs,
    threads = relate_threads,
    mem_gb = relate_mem_gb,
    ugc_wgw_container_registry = ugc_wgw_container_registry,
    runtime_attributes = default_runtime_attributes
  }

  # ---- ancestry (the map keys are read only here, so a map without them still runs the stage) ----

  if (run_ancestry) {
    call Somalier.ugc_wgw_somalier_ancestry { input:
      out_prefix = out_prefix,
      somalier_files = member_somalier,
      labelled_tar = ref_map["somalier_1kg_tar"],  # !FileCoercion
      labels = ref_map["somalier_1kg_labels"],  # !FileCoercion
      mem_gb = ancestry_mem_gb,
      ugc_wgw_container_registry = ugc_wgw_container_registry,
      runtime_attributes = default_runtime_attributes
    }
  }

  # ---- the tables the driver reads, and provenance ----

  call Somalier.ugc_wgw_qc_summary { input:
    out_prefix = out_prefix,
    pairs = ugc_wgw_somalier_relate.pairs,
    samples = ugc_wgw_somalier_relate.samples,
    ancestry_tsv = ugc_wgw_somalier_ancestry.ancestry_tsv,
    sample_ids = sample_ids,
    sample_sexes = sample_sexes,
    sample_sex_sources = sample_sex_sources,
    legacy_sample_ids = legacy_flag,
    runtime_attributes = default_runtime_attributes
  }

  call Provenance.ugc_wgw_manifest_write { input:
    ugc_wgw_version = ugc_wgw_version,
    stage = "cohort_qc",
    subject_type = "cohort",
    subject_id = cohort_id,
    member_ids = sample_ids,
    runtime_attributes = default_runtime_attributes
  }

  output {
    File qc_pairs = ugc_wgw_somalier_relate.pairs
    File qc_samples = ugc_wgw_somalier_relate.samples
    File qc_groups = ugc_wgw_somalier_relate.groups
    File qc_html = ugc_wgw_somalier_relate.html
    File qc_ped = ugc_wgw_somalier_ped.ped
    File qc_summary = ugc_wgw_qc_summary.summary
    File qc_members = ugc_wgw_qc_summary.members
    File qc_related_pairs = ugc_wgw_qc_summary.related_pairs
    File? qc_ancestry = ugc_wgw_somalier_ancestry.ancestry_tsv
    File? qc_ancestry_html = ugc_wgw_somalier_ancestry.ancestry_html
    Array[File] somalier_extracts = member_somalier
    String ugc_wgw_workflow_version = ugc_wgw_version
    File ugc_wgw_manifest = ugc_wgw_manifest_write.ugc_wgw_manifest
  }
}
