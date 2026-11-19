version 1.0

# ugc-pacbio-wgw — stage `cohort_graph_gt`: the structural variants of another cohort's pangenome graph genotyped in every member of this cohort from its reads' alignment (kanpig)
# origin: none
# see docs/DESIGN.md §5, §18 (2026-11-19); docs/ENTRYPOINTS.md §9.6; docs/expansion/pangenome-genotyping.design.md

import "../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/workflows/backend_configuration/backend_configuration.wdl" as BackendConfiguration
import "ugc_wgw/cohort/bundle.wdl" as Bundle
import "ugc_wgw/cohort/regions.wdl" as Regions
import "ugc_wgw/graph_gt/cohort.wdl" as GraphGt
import "ugc_wgw/provenance.wdl" as Provenance

workflow ugc_wgw_cohort_graph_gt {
  meta {
    description: "ugc-pacbio-wgw stage cohort_graph_gt: a graph is built from a subset of a cohort, so its structural variants have genotypes for the assembled members only. This stage takes the graph's sequence-resolved SV alleles (cohort_pangenome's vcfwave VCF, biallelic, 50 bp and longer) as a fixed list of sites and genotypes every member of a cohort at them from the alignment to the reference the member already has: a pileup index per member (kanpig plup; it depends on the member's BAM alone and is replayed by the call cache for the next graph), the genotypes from it (kanpig gt, ploidy by the member's sex), the members' columns joined per region shard into one cohort VCF with allele frequencies over the whole cohort beside the graph's own, and for the members that are in the graph the genotype from the reads against the genotype from the assembly. Genotyping, not discovery: no new variant appears. Alleles longer than size_max are listed and left uncalled. A research-grade genotype table with a measured error, not a truth set."
    outputs: {
      graph_gt_vcf: {
        description: "<cohort_id>.graph_gt.<graph_cohort_id>.<ref>.vcf.gz: every site with every member's genotype (GT unphased, GQ, DP, AD), AC/AN/AF/NS and the XX/XY strata over the cohort, GRAPH_AC/GRAPH_AN of the graph's members"
      },
      graph_gt_vcf_index: {
        description: "Its index"
      },
      graph_gt_sites_vcf: {
        description: "The same without the member columns: the frequency resource"
      },
      graph_gt_sites_vcf_index: {
        description: "Its index"
      },
      graph_gt_members: {
        description: "<prefix>.members.tsv: per member its sex, whether it has reads and is in the graph, genotyped sites, call rate, genotype counts, mean depth"
      },
      graph_gt_concordance: {
        description: "<prefix>.concordance.tsv: per member of the graph, SV type, size class and context (inside a region of the SV caller's exclusion list, inside a tandem repeat of the TRGT catalogue, or neither), the read genotype against the assembly genotype"
      },
      graph_gt_skipped: {
        description: "<prefix>.skipped.tsv: the alleles of the graph that are listed and not genotyped (too long) or left out of the list (symbolic, no length change), with the reason"
      },
      graph_gt_summary: {
        description: "<prefix>.summary.tsv: the graph and its members, the cohort's members, sites by type and size class, call rate, the comparison's totals, the parameters"
      },
      member_vcfs: {
        description: "<sample_id>.graph_gt.<graph_cohort_id>.<ref>.vcf.gz per genotyped member: every field kanpig writes (phase set included)"
      },
      member_vcf_indices: {
        description: "Their indexes"
      },
      plups: {
        description: "<sample_id>.<ref>.plup.gz per genotyped member: the pileup index, kept for the next graph"
      },
      plup_indices: {
        description: "Their indexes"
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
      description: "The cohort to genotype; output prefix"
    }
    sample_ids: {
      description: "The cohort's members, in order: the columns of the cohort VCF"
    }
    sample_sexes: {
      description: "XX, XY or empty per member (the sex the read-based track resolved)"
    }
    read_sample_ids: {
      description: "The members that have an alignment, in cohort order; the others are uncalled at every site"
    }
    read_sample_sexes: {
      description: "Their sexes (XX, XY or empty), in the same order: the ploidy of chrX and chrY"
    }
    aligned_bams: {
      description: "Their alignments to the reference (merged_haplotagged_bam of singleton or downstream), in the same order"
    }
    aligned_bam_indices: {
      description: "Their indexes"
    }
    graph_cohort_id: {
      description: "The cohort whose cohort_pangenome graph gives the sites (it may be this cohort)"
    }
    graph_vcf: {
      description: "That run's pangenome_wave_vcf"
    }
    graph_vcf_index: {
      description: "Its index"
    }
    size_min: {
      description: "Smallest length change (bp) of an allele in the site list"
    }
    size_max: {
      description: "Alleles longer than this are listed and not genotyped (kanpig reads pileups, not split alignments)"
    }
    plup_threads: {
      description: "Threads of the pileup task"
    }
    plup_mem_gb: {
      description: "Memory (GB) of the pileup task"
    }
    gt_threads: {
      description: "Threads of the genotyping task"
    }
    gt_mem_gb: {
      description: "Memory (GB) of the genotyping task"
    }
    join_mem_gb: {
      description: "Memory (GB) of each shard's join"
    }
    ref_map_file: {
      description: "Reference map: name, fasta, fasta_index, scatter_regions, trgt_tandem_repeat_bed, sawfish_exclude_bed and the two sawfish expected-copy-number BEDs are used"
    }
    ugc_wgw_container_registry: {
      description: "Registry holding the ugc-pacbio-wgw images (kanpig)"
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
    Array[String] read_sample_ids
    Array[String] read_sample_sexes
    Array[File] aligned_bams
    Array[File] aligned_bam_indices
    String graph_cohort_id
    File graph_vcf
    File graph_vcf_index

    Int size_min = 50
    Int size_max = 10000
    Int plup_threads = 8
    Int plup_mem_gb = 16
    Int gt_threads = 4
    Int gt_mem_gb = 8
    Int join_mem_gb = 8

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
  String out_prefix = "~{cohort_id}.graph_gt.~{graph_cohort_id}.~{ref_name}"

  call GraphGt.ugc_wgw_graph_gt_sites { input:
    graph_vcf = graph_vcf,
    graph_vcf_index = graph_vcf_index,
    reference = ref_map["fasta"],  # !FileCoercion
    reference_index = ref_map["fasta_index"],  # !FileCoercion
    tandem_repeat_bed = ref_map["trgt_tandem_repeat_bed"],  # !FileCoercion
    excluded_regions_bed = ref_map["sawfish_exclude_bed"],  # !FileCoercion
    size_min = size_min,
    size_max = size_max,
    out_prefix = out_prefix,
    runtime_attributes = default_runtime_attributes
  }

  # one member, two tasks: the pileup index is keyed by the member's BAM alone and so outlives the graph
  scatter (member_index in range(length(read_sample_ids))) {
    call GraphGt.ugc_wgw_kanpig_plup { input:
      sample_id = read_sample_ids[member_index],
      bam = aligned_bams[member_index],
      bam_index = aligned_bam_indices[member_index],
      ref_name = ref_name,
      size_min = size_min,
      size_max = size_max,
      threads = plup_threads,
      mem_gb = plup_mem_gb,
      ugc_wgw_container_registry = ugc_wgw_container_registry,
      runtime_attributes = default_runtime_attributes
    }

    call GraphGt.ugc_wgw_kanpig_gt { input:
      sample_id = read_sample_ids[member_index],
      sex = read_sample_sexes[member_index],
      plup = ugc_wgw_kanpig_plup.plup,
      plup_index = ugc_wgw_kanpig_plup.plup_index,
      sites = ugc_wgw_graph_gt_sites.sites,
      sites_index = ugc_wgw_graph_gt_sites.sites_index,
      reference = ref_map["fasta"],  # !FileCoercion
      reference_index = ref_map["fasta_index"],  # !FileCoercion
      expected_cn_xx = ref_map["sawfish_expected_bed_female"],  # !FileCoercion
      expected_cn_xy = ref_map["sawfish_expected_bed_male"],  # !FileCoercion
      size_min = size_min,
      size_max = size_max,
      out_prefix = "~{read_sample_ids[member_index]}.graph_gt.~{graph_cohort_id}.~{ref_name}",
      threads = gt_threads,
      mem_gb = gt_mem_gb,
      ugc_wgw_container_registry = ugc_wgw_container_registry,
      runtime_attributes = default_runtime_attributes
    }
  }

  # every shard's join reads every member's VCF: one tar instead of two files per member on its command line
  call Bundle.ugc_wgw_files_bundle as bundle_members { input:
    files = flatten([ugc_wgw_kanpig_gt.vcf, ugc_wgw_kanpig_gt.vcf_index]),
    name = "~{out_prefix}.members",
    runtime_attributes = default_runtime_attributes
  }

  call Regions.ugc_wgw_regions_prepare { input:
    scatter_regions = ref_map["scatter_regions"],  # !FileCoercion
    require_whole_contigs = false,
    runtime_attributes = default_runtime_attributes
  }

  scatter (shard_index in range(length(ugc_wgw_regions_prepare.shard_names))) {
    call GraphGt.ugc_wgw_graph_gt_join { input:
      members_bundle = bundle_members.bundle,
      read_sample_ids = read_sample_ids,
      sample_ids = sample_ids,
      sample_sexes = sample_sexes,
      sites = ugc_wgw_graph_gt_sites.sites,
      sites_index = ugc_wgw_graph_gt_sites.sites_index,
      graph_calls = ugc_wgw_graph_gt_sites.graph_calls,
      graph_calls_index = ugc_wgw_graph_gt_sites.graph_calls_index,
      regions_bed = ugc_wgw_regions_prepare.shard_beds[shard_index],
      out_prefix = "~{out_prefix}.~{ugc_wgw_regions_prepare.shard_names[shard_index]}",
      mem_gb = join_mem_gb,
      runtime_attributes = default_runtime_attributes
    }
  }

  call GraphGt.ugc_wgw_graph_gt_gather { input:
    vcfs = ugc_wgw_graph_gt_join.vcf,
    member_counts = ugc_wgw_graph_gt_join.members,
    concordance_counts = ugc_wgw_graph_gt_join.concordance,
    site_counts = ugc_wgw_graph_gt_join.site_counts,
    sites_counts = ugc_wgw_graph_gt_sites.counts,
    sample_ids = sample_ids,
    sample_sexes = sample_sexes,
    read_sample_ids = read_sample_ids,
    graph_cohort_id = graph_cohort_id,
    parameters = "size_min=~{size_min};size_max=~{size_max};genotyper=kanpig",
    out_prefix = out_prefix,
    runtime_attributes = default_runtime_attributes
  }

  call Provenance.ugc_wgw_manifest_write { input:
    ugc_wgw_version = ugc_wgw_version,
    stage = "cohort_graph_gt",
    subject_type = "cohort",
    subject_id = cohort_id,
    member_ids = sample_ids,
    runtime_attributes = default_runtime_attributes
  }

  output {
    File graph_gt_vcf = ugc_wgw_graph_gt_gather.vcf
    File graph_gt_vcf_index = ugc_wgw_graph_gt_gather.vcf_index
    File graph_gt_sites_vcf = ugc_wgw_graph_gt_gather.sites_vcf
    File graph_gt_sites_vcf_index = ugc_wgw_graph_gt_gather.sites_vcf_index
    File graph_gt_members = ugc_wgw_graph_gt_gather.members
    File graph_gt_concordance = ugc_wgw_graph_gt_gather.concordance
    File graph_gt_skipped = ugc_wgw_graph_gt_sites.skipped
    File graph_gt_summary = ugc_wgw_graph_gt_gather.summary
    Array[File] member_vcfs = ugc_wgw_kanpig_gt.vcf
    Array[File] member_vcf_indices = ugc_wgw_kanpig_gt.vcf_index
    Array[File] plups = ugc_wgw_kanpig_plup.plup
    Array[File] plup_indices = ugc_wgw_kanpig_plup.plup_index
    String ugc_wgw_workflow_version = ugc_wgw_version
    File ugc_wgw_manifest = ugc_wgw_manifest_write.ugc_wgw_manifest
  }
}
