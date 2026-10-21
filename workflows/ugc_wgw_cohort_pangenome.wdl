version 1.0
# ugc-pacbio-wgw — stage `cohort_pangenome`: a Minigraph-Cactus graph of the cohort's haplotype assemblies, built per chromosome, its cohort VCFs and their allele frequencies
# origin: none
# see docs/DESIGN.md §5, §13; docs/ENTRYPOINTS.md §8.2
import "../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/workflows/backend_configuration/backend_configuration.wdl" as BackendConfiguration
import "ugc_wgw/cohort/concat.wdl" as Concat
import "ugc_wgw/cohort/freq.wdl" as Freq
import "ugc_wgw/pangenome/cactus_chrom.wdl" as Chrom
import "ugc_wgw/pangenome/cactus_split.wdl" as Split
import "ugc_wgw/pangenome/names.wdl" as Names
import "ugc_wgw/pangenome/truvari.wdl" as Truvari
import "ugc_wgw/pangenome/vcf_post.wdl" as Post
import "ugc_wgw/provenance.wdl" as Provenance

workflow ugc_wgw_cohort_pangenome {
  meta {
    description: "ugc-pacbio-wgw stage cohort_pangenome: the members' haplotype assemblies (the assembly stage's contract) as a Minigraph-Cactus graph on the decoy-free GRCh38. One genome-wide split step maps every haplotype to a reference-only graph and bins contigs by chromosome (Cactus's --mgSplit, with its ambiguity rule); one task per requested chromosome constructs, maps, aligns and joins (GBZ, GFA, VCF, vcfwave VCF); the VCFs are gathered, renamed to the sample IDs and unphased for members that are not trio-binned; allele frequencies come from the cohort_freq tasks on the vcfwave VCF with missing genotypes kept missing (no assembled sequence there) and the XX/XY strata from the members' inferred karyotypes; optionally truvari bench against the cohort's read-based SV merge. The graph's paths of non-trio members are not haplotypes (hifiasm's local phasing): no LD, imputation or PanGenie-style use without statistical phasing."
    outputs: {
      pangenome_vcf: {
        description: "<cohort_id>.pangenome.<ref_name>.vcf.gz: the cohort's graph VCF (vg deconstruct) over the requested chromosomes, sample IDs restored, non-trio members unphased"
      },
      pangenome_vcf_index: {
        description: "Tabix index"
      },
      pangenome_wave_vcf: {
        description: "The same after vcfwave: alleles realigned to the reference, complex records decomposed"
      },
      pangenome_wave_vcf_index: {
        description: "Tabix index"
      },
      chromosomes_built: {
        description: "The chromosomes built, the order of every per-chromosome array"
      },
      chrom_gbzs: {
        description: "The per-chromosome GBZ graphs"
      },
      chrom_gfas: {
        description: "The per-chromosome GFAs (gzipped)"
      },
      chrom_graph_stats: {
        description: "graph-stats.tsv per chromosome (nodes, edges, length)"
      },
      chrom_logs: {
        description: "The Cactus logs per chromosome, tarred"
      },
      name_map: {
        description: "sample_id, Cactus-safe name, trio per member"
      },
      split_log: {
        description: "cactus-graphmap-split's minigraph.split.log"
      },
      contig_sizes: {
        description: "The split step's contig_sizes.tsv"
      },
      ambiguous_contigs: {
        description: "genome, assigned, ambiguous contig counts from the split"
      },
      freq_vcf: {
        description: "<cohort_id>.pangenome.<ref_name>.freq.vcf.gz: sites-only frequencies of the vcfwave VCF (AC, AN, AF, NS, per-stratum tags)"
      },
      freq_vcf_index: {
        description: "Tabix index"
      },
      freq_summary: {
        description: "The frequency summary (resource, metric, key, value)"
      },
      freq_samples: {
        description: "Members with their stratum and whether the VCF has them"
      },
      truvari_summary: {
        description: "truvari bench summary.json of the graph's SVs against the cohort_merge SV VCF (with sv_vcf)"
      },
      truvari_details: {
        description: "truvari's output directory, tarred (with sv_vcf)"
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
      description: "Frozen cohort ID; output prefix"
    }
    sample_ids: {
      description: "Members in the frozen cohort's order"
    }
    hap1_fastas: {
      description: "Each member's first haplotype (the assembly stage's hap_fastas[0]), same order"
    }
    hap2_fastas: {
      description: "Each member's second haplotype, same order"
    }
    trios: {
      description: "Whether each member's assembly is trio-binned (its genotypes keep their phase)"
    }
    sample_sexes: {
      description: "XX, XY or empty per member, from assembly_call's inferred karyotype: the frequency strata"
    }
    sample_sex_sources: {
      description: "Where each member's sex came from (karyotype, or empty)"
    }
    chromosomes: {
      description: "Chromosomes to build; default: the map's asm_graph_chromosomes (chr1-22, X, Y). The split step is genome-wide regardless"
    }
    build_other_contigs: {
      description: "Also build the graph of the other contig (chrM and the unplaced contigs lumped together)"
    }
    sv_vcf: {
      description: "The cohort_merge SV VCF to compare the graph's SVs with (truvari bench), when given"
    }
    sv_vcf_index: {
      description: "Its index"
    }
    split_cpu: {
      description: "Cores of the genome-wide split step (every haplotype mapped)"
    }
    split_mem_gb: {
      description: "Memory (GB) of the split step"
    }
    split_disk_gb: {
      description: "Disk (GB) of the split step; default from the input sizes"
    }
    chrom_cpu: {
      description: "Cores of each chromosome task"
    }
    chrom_mem_gb_floor: {
      description: "Memory (GB) of a chromosome task at least"
    }
    chrom_mem_gb_per_100mb: {
      description: "Memory (GB) per 100 Mb of chromosome; the larger of the two rules applies"
    }
    freq_threads: {
      description: "Threads of each frequency shard"
    }
    freq_mem_gb: {
      description: "Memory (GB) of each frequency shard"
    }
    ref_map_file: {
      description: "ugc-wgw reference map (docs/DESIGN.md §10.4); keys used here: name, asm_fasta, asm_fasta_index, asm_graph_chromosomes"
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
      description: "Optional type of GPU/Accelerator to use (no GPU task in this stage)"
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
    Array[File] hap1_fastas
    Array[File] hap2_fastas
    Array[Boolean] trios
    Array[String] sample_sexes
    Array[String] sample_sex_sources = []
    Array[String]? chromosomes
    Boolean build_other_contigs = false
    File? sv_vcf
    File? sv_vcf_index
    Int split_cpu = 32
    Int split_mem_gb = 128
    Int? split_disk_gb
    Int chrom_cpu = 16
    Int chrom_mem_gb_floor = 64
    Int chrom_mem_gb_per_100mb = 32
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
  String ref_name = ref_map["name"]
  String out_prefix = "~{cohort_id}.pangenome.~{ref_name}"
  Array[String] graph_chromosomes = read_lines(ref_map["asm_graph_chromosomes"])  # !FileCoercion
  Array[String] requested = select_first([chromosomes, graph_chromosomes])
  Array[String] other = if build_other_contigs then ["chrOther"] else []
  Array[String] to_build = flatten([requested, other])

  call Names.ugc_wgw_pangenome_names { input:
    cohort_id = cohort_id,
    sample_ids = sample_ids,
    trios = trios,
    ref_index = ref_map["asm_fasta_index"],  # !FileCoercion
    chromosomes = requested,
    chrom_mem_gb_floor = chrom_mem_gb_floor,
    chrom_mem_gb_per_100mb = chrom_mem_gb_per_100mb,
    runtime_attributes = default_runtime_attributes
  }

  # hap1, hap2 of each member in cohort order, matching the genome names <safe>.1, <safe>.2
  Array[File] hap_fastas = flatten(transpose([hap1_fastas, hap2_fastas]))

  call Split.ugc_wgw_cactus_split { input:
    cohort_id = cohort_id,
    ref_name = ref_name,
    ref_fasta = ref_map["asm_fasta"],  # !FileCoercion
    genome_names = ugc_wgw_pangenome_names.genome_names,
    hap_fastas = hap_fastas,
    ref_contigs = ref_map["asm_graph_chromosomes"],  # !FileCoercion
    cpu = split_cpu,
    mem_gb = split_mem_gb,
    disk_gb = split_disk_gb,
    runtime_attributes = default_runtime_attributes
  }

  scatter (chrom in to_build) {
    Int chrom_mem = if chrom == "chrOther" then chrom_mem_gb_floor else ugc_wgw_pangenome_names.chrom_mem_gb[chrom]
    call Chrom.ugc_wgw_cactus_chrom { input:
      cohort_id = cohort_id,
      ref_name = ref_name,
      chromosome = chrom,
      chrom_tars = ugc_wgw_cactus_split.chrom_tars,
      chromfile = ugc_wgw_cactus_split.chromfile,
      split_root = ugc_wgw_cactus_split.split_root,
      cpu = chrom_cpu,
      mem_gb = chrom_mem,
      runtime_attributes = default_runtime_attributes
    }
  }

  Array[String] header_common = [
    "##ugc_wgw_cohort_pangenome_version=~{ugc_wgw_version}",
    "##ugc_wgw_cohort_pangenome_cohort=~{cohort_id}",
    "##ugc_wgw_cohort_pangenome_graph=Minigraph-Cactus v3.3.0, per-chromosome construction (--mgSplit) on the decoy-free GRCh38",
    "##ugc_wgw_cohort_pangenome_ambiguity=contigs enter a chromosome's graph only when at least 75/50/25 % of their length (below 100 kb, below 1 Mb, from 1 Mb) aligns to one reference contig with a best/second ratio of 2 or more at MAPQ 5 or more (cactus-graphmap-split); other contigs are in no graph",
    "##ugc_wgw_cohort_pangenome_phase=genotypes of trio-binned members are phased (paternal|maternal); the others are unphased: hifiasm's hap1/hap2 are local phase blocks, so the graph's paths of those members are not haplotypes (no LD, imputation or PanGenie-style use without statistical phasing)"
  ]

  call Post.ugc_wgw_pangenome_vcf_post as post_raw { input:
    chromosomes = to_build,
    vcfs = ugc_wgw_cactus_chrom.vcf,
    vcf_indices = ugc_wgw_cactus_chrom.vcf_index,
    rename_samples = ugc_wgw_pangenome_names.rename_samples,
    unphase_samples = ugc_wgw_pangenome_names.unphase_samples,
    header_lines = header_common,
    out_prefix = out_prefix,
    runtime_attributes = default_runtime_attributes
  }

  call Post.ugc_wgw_pangenome_vcf_post as post_wave { input:
    chromosomes = to_build,
    vcfs = ugc_wgw_cactus_chrom.wave_vcf,
    vcf_indices = ugc_wgw_cactus_chrom.wave_vcf_index,
    rename_samples = ugc_wgw_pangenome_names.rename_samples,
    unphase_samples = ugc_wgw_pangenome_names.unphase_samples,
    header_lines = flatten([header_common, ["##ugc_wgw_cohort_pangenome_vcfwave=alleles realigned to the reference and complex records decomposed by vcfwave"]]),
    out_prefix = "~{out_prefix}.wave",
    runtime_attributes = default_runtime_attributes
  }

  # frequencies of the vcfwave VCF, one shard per requested chromosome; a missing genotype stays missing (no assembled sequence)
  Array[String] header_freq = flatten([header_common, [
    "##ugc_wgw_cohort_freq_strata=XX,XY (the karyotype assembly_call inferred; members of unknown karyotype are counted in the unstratified tags only)",
    "##ugc_wgw_cohort_freq_missing_gt=./. is kept: the haplotype has no assembled sequence there, so the call rate (NS over members) reports assembly coverage"
  ]])
  scatter (i in range(length(requested))) {
    call Freq.ugc_wgw_bcftools_freq as freq_shard { input:
      vcf = post_wave.vcf,
      vcf_index = post_wave.vcf_index,
      regions_bed = ugc_wgw_pangenome_names.chrom_beds[i],
      sample_ids = sample_ids,
      sample_sexes = sample_sexes,
      header_lines = header_freq,
      missing_to_homref = false,
      out_prefix = "~{out_prefix}.freq.~{requested[i]}",
      threads = freq_threads,
      mem_gb = freq_mem_gb,
      runtime_attributes = default_runtime_attributes
    }
  }
  call Concat.ugc_wgw_bcftools_concat as concat_freq { input:
    vcfs = freq_shard.shard_vcf,
    vcf_indices = freq_shard.shard_vcf_index,
    out_prefix = "~{out_prefix}.freq",
    runtime_attributes = default_runtime_attributes
  }
  call Post.ugc_wgw_pangenome_freq_summary { input:
    freq_vcf = concat_freq.vcf,
    freq_vcf_index = concat_freq.vcf_index,
    vcf_samples = freq_shard.vcf_samples[0],
    sample_ids = sample_ids,
    sample_sexes = sample_sexes,
    sample_sex_sources = sample_sex_sources,
    out_prefix = "~{out_prefix}.freq",
    runtime_attributes = default_runtime_attributes
  }

  if (defined(sv_vcf)) {
    # truvari refuses multiallelic records: the comparison copy is split, the primary VCFs are not
    call Truvari.ugc_wgw_vcf_split_multiallelic { input:
      vcf = post_wave.vcf,
      vcf_index = post_wave.vcf_index,
      out_prefix = "~{out_prefix}.wave",
      runtime_attributes = default_runtime_attributes
    }
    call Truvari.ugc_wgw_truvari_bench { input:
      base_vcf = select_first([sv_vcf]),
      base_vcf_index = select_first([sv_vcf_index]),
      comp_vcf = ugc_wgw_vcf_split_multiallelic.biallelic_vcf,
      comp_vcf_index = ugc_wgw_vcf_split_multiallelic.biallelic_vcf_index,
      out_prefix = out_prefix,
      runtime_attributes = default_runtime_attributes
    }
  }

  call Provenance.ugc_wgw_manifest_write { input:
    ugc_wgw_version = ugc_wgw_version,
    stage = "cohort_pangenome",
    subject_type = "cohort",
    subject_id = cohort_id,
    member_ids = sample_ids,
    upstream_workflow_name = "Minigraph-Cactus v3.3.0 (cactus-minigraph, cactus-graphmap, cactus-graphmap-split, cactus-align, cactus-graphmap-join)",
    runtime_attributes = default_runtime_attributes
  }

  output {
    File pangenome_vcf = post_raw.vcf
    File pangenome_vcf_index = post_raw.vcf_index
    File pangenome_wave_vcf = post_wave.vcf
    File pangenome_wave_vcf_index = post_wave.vcf_index
    Array[String] chromosomes_built = to_build
    Array[File] chrom_gbzs = ugc_wgw_cactus_chrom.gbz
    Array[File] chrom_gfas = ugc_wgw_cactus_chrom.gfa
    Array[File] chrom_graph_stats = ugc_wgw_cactus_chrom.graph_stats
    Array[File] chrom_logs = ugc_wgw_cactus_chrom.logs
    File name_map = ugc_wgw_pangenome_names.name_map
    File split_log = ugc_wgw_cactus_split.split_log
    File contig_sizes = ugc_wgw_cactus_split.contig_sizes
    File ambiguous_contigs = ugc_wgw_cactus_split.ambiguous_contigs
    File freq_vcf = concat_freq.vcf
    File freq_vcf_index = concat_freq.vcf_index
    File freq_summary = ugc_wgw_pangenome_freq_summary.summary
    File freq_samples = ugc_wgw_pangenome_freq_summary.samples
    File? truvari_summary = ugc_wgw_truvari_bench.summary
    File? truvari_details = ugc_wgw_truvari_bench.details
    String ugc_wgw_workflow_version = ugc_wgw_version
    File ugc_wgw_manifest = ugc_wgw_manifest_write.ugc_wgw_manifest
  }
}
