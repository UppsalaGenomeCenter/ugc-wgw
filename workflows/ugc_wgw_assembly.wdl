version 1.0

# ugc-pacbio-wgw — stage `assembly`: hifiasm de novo assembly of one sample, trio-binned with yak when both parents' reads are given
# origin: none
# see docs/DESIGN.md §5, §13; docs/ENTRYPOINTS.md §8

import "../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/tasks/samtools.wdl" as Samtools
import "../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/workflows/backend_configuration/backend_configuration.wdl" as BackendConfiguration
import "ugc_wgw/assembly/align_hifiasm.wdl" as Align
import "ugc_wgw/assembly/bcftools_zip_stats.wdl" as ZipStats
import "ugc_wgw/assembly/gfatools.wdl" as Gfatools
import "ugc_wgw/assembly/hifiasm.wdl" as Hifiasm
import "ugc_wgw/assembly/yak.wdl" as Yak
import "ugc_wgw/provenance.wdl" as Provenance

workflow ugc_wgw_assembly {
  meta {
    description: "ugc-pacbio-wgw stage assembly: samtools fasta, hifiasm (trio-binned with parental yak tables when father_hifi_reads and mother_hifi_reads are given), gfatools + calN50 per haplotype, minimap2 asm5 alignment to the reference, samtools merge of the haplotypes, paftools variant calls with bcftools stats. Ported from HiFi-human-assembly-WDL v1.0.2 into the WGS tree's conventions; aligns to the ref_map reference only."
    outputs: {
      assembly_noseq_gfas: {
        description: "Sequence-less GFAs (haplotype contigs and unitig graphs)"
      },
      assembly_lowQ_beds: {
        description: "Low-quality regions of the contigs and unitig graphs"
      },
      haplotypes: {
        description: "hap1, hap2: the order of every per-haplotype array below"
      },
      zipped_assembly_fastas: {
        description: "Haplotype contigs, bgzipped FASTA"
      },
      assembly_stats: {
        description: "calN50 statistics per haplotype"
      },
      asm_bams: {
        description: "Contigs aligned to the reference, per haplotype"
      },
      asm_bam_indices: {
        description: "Indices for asm_bams"
      },
      merged_asm_bam: {
        description: "Both haplotypes' alignments merged"
      },
      merged_asm_bam_index: {
        description: "Index for merged_asm_bam"
      },
      paftools_vcfs: {
        description: "paftools variant calls per haplotype (bgzipped)"
      },
      paftools_vcf_indices: {
        description: "Indices for paftools_vcfs"
      },
      paftools_vcf_stats: {
        description: "bcftools stats per haplotype VCF"
      },
      trio: {
        description: "Whether the assembly was trio-binned"
      },
      hap1_parent: {
        description: "Sample whose k-mers define hap1 (the father) in trio mode"
      },
      hap2_parent: {
        description: "Sample whose k-mers define hap2 (the mother) in trio mode"
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
    sample_id: {
      description: "Sample ID; output prefix is <sample_id>.asm"
    }
    hifi_reads: {
      description: "Unaligned HiFi BAMs of the sample"
    }
    father_id: {
      description: "Father's sample ID (trio mode; required with father_hifi_reads)"
    }
    mother_id: {
      description: "Mother's sample ID (trio mode; required with mother_hifi_reads)"
    }
    father_hifi_reads: {
      description: "Father's unaligned HiFi BAMs; with mother_hifi_reads this enables trio binning"
    }
    mother_hifi_reads: {
      description: "Mother's unaligned HiFi BAMs"
    }
    hifiasm_extra_params: {
      description: "Additional hifiasm parameters, appended in both modes"
    }
    hifiasm_threads: {
      description: "hifiasm threads"
    }
    hifiasm_mem_gb: {
      description: "hifiasm memory (GB)"
    }
    yak_threads: {
      description: "yak count threads (trio mode)"
    }
    yak_mem_gb: {
      description: "yak count memory (GB); default 70 at low parental depth, else 50"
    }
    yak_params: {
      description: "yak count parameters; default -b0 at low parental depth, else -b37"
    }
    hifiasm_trio_params: {
      description: "hifiasm trio binning thresholds; default -c1 -d1 at low parental depth, else -c2 -d5"
    }
    trio_low_depth_gb: {
      description: "A parent counts as low depth when its uncompressed FASTA is smaller than this many GB (60 GB is about 10x)"
    }
    ref_map_file: {
      description: "ugc-wgw reference map (docs/DESIGN.md §10.4); keys used here: name, fasta, fasta_index"
    }
    ugc_wgw_container_registry: {
      description: "Registry holding the ugc-built images (hifiasm)"
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
    String sample_id
    Array[File] hifi_reads
    String? father_id
    String? mother_id
    Array[File]? father_hifi_reads
    Array[File]? mother_hifi_reads
    String? hifiasm_extra_params
    Int hifiasm_threads = 48
    Int hifiasm_mem_gb = 288
    Int yak_threads = 24
    Int? yak_mem_gb
    String? yak_params
    String? hifiasm_trio_params
    Int trio_low_depth_gb = 90
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
  String ref_name = ref_map["name"]

  Boolean is_trio = defined(father_hifi_reads) && defined(mother_hifi_reads)

  scatter (bam in hifi_reads) {
    call Samtools.samtools_fasta as samtools_fasta_child { input:
      bam = bam,
      runtime_attributes = default_runtime_attributes
    }
  }

  if (is_trio) {
    String father = select_first([father_id])
    String mother = select_first([mother_id])

    scatter (bam in select_first([father_hifi_reads])) {
      call Samtools.samtools_fasta as samtools_fasta_father { input:
        bam = bam,
        runtime_attributes = default_runtime_attributes
      }
    }
    scatter (bam in select_first([mother_hifi_reads])) {
      call Samtools.samtools_fasta as samtools_fasta_mother { input:
        bam = bam,
        runtime_attributes = default_runtime_attributes
      }
    }

    # Upstream heuristic: below ~15x parental coverage keep singleton k-mers (yak -b0) and bin on
    # k-mers seen >= 1x in one parent and 0x in the other (-c1 -d1); otherwise bloom filter (-b37)
    # and the hifiasm defaults (-c2 -d5). 60 GB of uncompressed FASTA is about 10x.
    Boolean low_depth = size(samtools_fasta_father.fasta, "GB") < trio_low_depth_gb
      && size(samtools_fasta_mother.fasta, "GB") < trio_low_depth_gb
    String yak_params_used = select_first([yak_params, if low_depth then "-b0" else "-b37"])
    Int yak_mem_gb_used = select_first([yak_mem_gb, if low_depth then 70 else 50])
    String trio_params = select_first([hifiasm_trio_params, if low_depth then "-c1 -d1" else "-c2 -d5"])
    String hifiasm_params_trio = trio_params + " " + select_first([hifiasm_extra_params, ""])

    call Yak.ugc_wgw_yak_count as yak_count_father { input:
      sample_id = father,
      reads_fastas = samtools_fasta_father.fasta,
      yak_params = yak_params_used,
      threads = yak_threads,
      mem_gb = yak_mem_gb_used,
      runtime_attributes = default_runtime_attributes
    }
    call Yak.ugc_wgw_yak_count as yak_count_mother { input:
      sample_id = mother,
      reads_fastas = samtools_fasta_mother.fasta,
      yak_params = yak_params_used,
      threads = yak_threads,
      mem_gb = yak_mem_gb_used,
      runtime_attributes = default_runtime_attributes
    }
  }

  String? hifiasm_params = if is_trio then hifiasm_params_trio else hifiasm_extra_params

  call Hifiasm.ugc_wgw_hifiasm_assemble { input:
    sample_id = sample_id,
    reads_fastas = samtools_fasta_child.fasta,
    extra_params = hifiasm_params,
    father_yak = yak_count_father.yak,
    mother_yak = yak_count_mother.yak,
    threads = hifiasm_threads,
    mem_gb = hifiasm_mem_gb,
    ugc_wgw_container_registry = ugc_wgw_container_registry,
    runtime_attributes = default_runtime_attributes
  }

  # <sample_id>.asm.<bp|dip>.hap1.p_ctg.gfa -> hap1; the glob is lexicographic, so hap1 comes first
  scatter (gfa in ugc_wgw_hifiasm_assemble.assembly_hap_gfas) {
    String haplotype = sub(basename(gfa, ".p_ctg.gfa"), "^.*\\.", "")

    call Gfatools.ugc_wgw_gfatools_gfa2fa { input:
      gfa = gfa,
      runtime_attributes = default_runtime_attributes
    }

    call Align.ugc_wgw_minimap2_align_asm { input:
      sample_id = sample_id,
      haplotype = haplotype,
      query_fasta = ugc_wgw_gfatools_gfa2fa.zipped_fasta,
      ref_fasta = ref_map["fasta"],  # !FileCoercion
      ref_name = ref_name,
      runtime_attributes = default_runtime_attributes
    }

    call Align.ugc_wgw_paftools_call { input:
      bam = ugc_wgw_minimap2_align_asm.bam,
      bam_index = ugc_wgw_minimap2_align_asm.bam_index,
      ref_fasta = ref_map["fasta"],  # !FileCoercion
      sample_id = sample_id,
      runtime_attributes = default_runtime_attributes
    }

    call ZipStats.ugc_wgw_bcftools_zip_stats { input:
      vcf = ugc_wgw_paftools_call.vcf,
      ref_fasta = ref_map["fasta"],  # !FileCoercion
      ref_index = ref_map["fasta_index"],  # !FileCoercion
      sample_id = sample_id,
      runtime_attributes = default_runtime_attributes
    }
  }

  call Samtools.samtools_merge as merge_haps { input:
    bams = ugc_wgw_minimap2_align_asm.bam,
    out_prefix = "~{sample_id}.asm.~{ref_name}",
    runtime_attributes = default_runtime_attributes
  }

  Array[String] no_members = []

  call Provenance.ugc_wgw_manifest_write { input:
    ugc_wgw_version = ugc_wgw_version,
    stage = "assembly",
    subject_type = "sample",
    subject_id = sample_id,
    member_ids = if is_trio then select_all([father_id, mother_id]) else no_members,
    upstream_workflow_name = if is_trio then "de_novo_assembly_trio" else "de_novo_assembly_sample",
    runtime_attributes = default_runtime_attributes
  }

  output {
    Array[File] assembly_noseq_gfas = ugc_wgw_hifiasm_assemble.assembly_noseq_gfas
    Array[File] assembly_lowQ_beds = ugc_wgw_hifiasm_assemble.assembly_lowQ_beds
    Array[String] haplotypes = haplotype
    Array[File] zipped_assembly_fastas = ugc_wgw_gfatools_gfa2fa.zipped_fasta
    Array[File] assembly_stats = ugc_wgw_gfatools_gfa2fa.stats
    Array[File] asm_bams = ugc_wgw_minimap2_align_asm.bam
    Array[File] asm_bam_indices = ugc_wgw_minimap2_align_asm.bam_index
    File merged_asm_bam = merge_haps.merged_bam
    File merged_asm_bam_index = merge_haps.merged_bam_index
    Array[File] paftools_vcfs = ugc_wgw_bcftools_zip_stats.zipped_vcf
    Array[File] paftools_vcf_indices = ugc_wgw_bcftools_zip_stats.zipped_vcf_index
    Array[File] paftools_vcf_stats = ugc_wgw_bcftools_zip_stats.stats

    Boolean trio = is_trio
    String? hap1_parent = father
    String? hap2_parent = mother

    String ugc_wgw_workflow_version = ugc_wgw_version
    File ugc_wgw_manifest = ugc_wgw_manifest_write.ugc_wgw_manifest
  }
}
