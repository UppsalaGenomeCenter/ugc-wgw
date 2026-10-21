version 1.0

# ugc-pacbio-wgw — stage `assembly`: hifiasm de novo assembly of one sample, trio-binned with yak when both parents' reads are given; the output contract every later stage reads; HMM-Flagger and trio QC
# origin: none
# see docs/DESIGN.md §5, §13; docs/ENTRYPOINTS.md §8

import "../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/tasks/samtools.wdl" as Samtools
import "../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/workflows/backend_configuration/backend_configuration.wdl" as BackendConfiguration
import "ugc_wgw/assembly/flagger.wdl" as Flagger
import "ugc_wgw/assembly/gfatools.wdl" as Gfatools
import "ugc_wgw/assembly/hifiasm.wdl" as Hifiasm
import "ugc_wgw/assembly/yak.wdl" as Yak
import "ugc_wgw/fail.wdl" as Fail
import "ugc_wgw/provenance.wdl" as Provenance

workflow ugc_wgw_assembly {
  meta {
    description: "ugc-pacbio-wgw stage assembly: samtools fasta, hifiasm (trio-binned with parental yak tables when father_hifi_reads and mother_hifi_reads are given; --telo-m CCCTAA by default), gfatools gfa2fa and calN50 per haplotype. The assembler sits behind an output contract (assembler_used, hap_fastas, haplotype_labels, assembler_version); only hifiasm is implemented and any other value fails at the start. QC on the contract: HMM-Flagger (the reads aligned back to hap1+hap2, correct_bam, bam2cov, hmm_flagger) flags Err/Dup/Hap/Col regions; trio-binned samples get yak trioeval per haplotype. Calling against the reference is the assembly_call stage."
    outputs: {
      assembler_used: {
        description: "The assembler that ran (hifiasm)"
      },
      assembler_version: {
        description: "Its version line, one file"
      },
      haplotype_labels: {
        description: "Labels of the two haplotypes in the order of every per-haplotype array: hap1, hap2; pat, mat when trio-binned"
      },
      hap_fastas: {
        description: "Haplotype contigs, bgzipped FASTA, in haplotype order"
      },
      assembly_stats: {
        description: "calN50 statistics per haplotype"
      },
      assembly_noseq_gfas: {
        description: "Sequence-less GFAs (haplotype contigs and unitig graphs)"
      },
      assembly_lowQ_beds: {
        description: "Low-quality regions of the contigs and unitig graphs"
      },
      flagger_bed: {
        description: "HMM-Flagger's prediction for the diploid assembly: Err, Dup, Hap or Col per interval, in assembly coordinates"
      },
      flagger_summary: {
        description: "Flagger's summary table (prediction_summary_final.tsv)"
      },
      trio_eval: {
        description: "yak trioeval per haplotype (switch and Hamming error), trio-binned samples only"
      },
      read_alignment_bam: {
        description: "The reads aligned to the diploid assembly, only with keep_read_alignment"
      },
      read_alignment_bam_index: {
        description: "Its index, only with keep_read_alignment"
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
    assembler: {
      description: "The assembler behind the output contract; only hifiasm is implemented, anything else fails at the start"
    }
    hifiasm_telomere_motif: {
      description: "hifiasm --telo-m motif (CCCTAA for human); an empty string disables telomere identification. Changing it changes the hifiasm command and so the call-cache key of every sample"
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
    keep_read_alignment: {
      description: "Export the reads-to-assembly BAM Flagger worked on (tens of GB per sample; off by default)"
    }
    flagger_align_threads: {
      description: "minimap2 threads for the reads-to-assembly alignment"
    }
    flagger_align_mem_gb: {
      description: "Memory (GB) of that alignment (the -I8g index of the diploid assembly)"
    }
    flagger_threads: {
      description: "bam2cov and hmm_flagger threads"
    }
    flagger_mem_gb: {
      description: "Memory (GB) of bam2cov and hmm_flagger"
    }
    trioeval_c: {
      description: "yak trioeval -c; default: the -c the binning used"
    }
    trioeval_d: {
      description: "yak trioeval -d; default: the -d the binning used"
    }
    ref_map_file: {
      description: "ugc-wgw reference map (docs/DESIGN.md §10.4); keys used here: name"
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
    String sample_id
    String assembler = "hifiasm"
    String hifiasm_telomere_motif = "CCCTAA"
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
    Boolean keep_read_alignment = false
    Int flagger_align_threads = 32
    Int flagger_align_mem_gb = 64
    Int flagger_threads = 16
    Int flagger_mem_gb = 32
    Int? trioeval_c
    Int? trioeval_d
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

  # read for the input contract every entrypoint shares; this stage needs no key of it
  #@ except: DeclarationName
  Map[String, String] ref_map = read_map(ref_map_file)  # !UnusedDeclaration

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

  # the assembler seam: the contract below is all the rest of the stage and the later stages read
  if (assembler != "hifiasm") {
    call Fail.ugc_wgw_fail { input:  # !UnusedCall
      message = "assembler '" + assembler + "' is not implemented (only hifiasm)",
      runtime_attributes = default_runtime_attributes
    }
  }
  if (assembler == "hifiasm") {
    call Hifiasm.ugc_wgw_hifiasm_assemble { input:
      sample_id = sample_id,
      reads_fastas = samtools_fasta_child.fasta,
      extra_params = hifiasm_params,
      telomere_motif = hifiasm_telomere_motif,
      father_yak = yak_count_father.yak,
      mother_yak = yak_count_mother.yak,
      threads = hifiasm_threads,
      mem_gb = hifiasm_mem_gb,
      ugc_wgw_container_registry = ugc_wgw_container_registry,
      runtime_attributes = default_runtime_attributes
    }
  }
  Array[File] hap_gfas = select_first([ugc_wgw_hifiasm_assemble.assembly_hap_gfas])

  # <sample_id>.asm.<bp|dip>.hap1.p_ctg.gfa -> hap1; the glob is lexicographic, so hap1 comes first
  scatter (gfa in hap_gfas) {
    call Gfatools.ugc_wgw_gfatools_gfa2fa { input:
      gfa = gfa,
      runtime_attributes = default_runtime_attributes
    }
  }
  Array[String] labels = if is_trio then ["pat", "mat"] else ["hap1", "hap2"]

  # HMM-Flagger: the sample's reads back on its own diploid assembly
  call Flagger.ugc_wgw_flagger_align { input:
    sample_id = sample_id,
    reads_fastas = samtools_fasta_child.fasta,
    hap_fastas = ugc_wgw_gfatools_gfa2fa.zipped_fasta,
    threads = flagger_align_threads,
    mem_gb = flagger_align_mem_gb,
    ugc_wgw_container_registry = ugc_wgw_container_registry,
    runtime_attributes = default_runtime_attributes
  }
  call Flagger.ugc_wgw_flagger_correct_bam { input:
    sample_id = sample_id,
    bam = ugc_wgw_flagger_align.bam,
    bam_index = ugc_wgw_flagger_align.bam_index,
    runtime_attributes = default_runtime_attributes
  }
  call Flagger.ugc_wgw_flagger_hmm { input:
    sample_id = sample_id,
    bam = ugc_wgw_flagger_correct_bam.corrected_bam,
    bam_index = ugc_wgw_flagger_correct_bam.corrected_bam_index,
    dip_fasta = ugc_wgw_flagger_align.dip_fasta,
    threads = flagger_threads,
    mem_gb = flagger_mem_gb,
    runtime_attributes = default_runtime_attributes
  }
  if (keep_read_alignment) {
    File kept_bam = ugc_wgw_flagger_align.bam
    File kept_bam_index = ugc_wgw_flagger_align.bam_index
  }

  # trio QC: each haplotype against the parental k-mer tables, with the binning's thresholds
  if (is_trio) {
    Boolean trio_low = select_first([low_depth])
    Int eval_c = select_first([trioeval_c, if trio_low then 1 else 2])
    Int eval_d = select_first([trioeval_d, if trio_low then 1 else 5])
    scatter (i in range(length(labels))) {
      call Yak.ugc_wgw_yak_trioeval { input:
        sample_id = sample_id,
        haplotype = labels[i],
        hap_fasta = ugc_wgw_gfatools_gfa2fa.zipped_fasta[i],
        father_yak = select_first([yak_count_father.yak]),
        mother_yak = select_first([yak_count_mother.yak]),
        min_occurrence = eval_c,
        mid_occurrence = eval_d,
        runtime_attributes = default_runtime_attributes
      }
    }
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
    String assembler_used = assembler
    File assembler_version = select_first([ugc_wgw_hifiasm_assemble.version])
    Array[String] haplotype_labels = labels
    Array[File] hap_fastas = ugc_wgw_gfatools_gfa2fa.zipped_fasta
    Array[File] assembly_stats = ugc_wgw_gfatools_gfa2fa.stats
    Array[File] assembly_noseq_gfas = select_first([ugc_wgw_hifiasm_assemble.assembly_noseq_gfas])
    Array[File] assembly_lowQ_beds = select_first([ugc_wgw_hifiasm_assemble.assembly_lowQ_beds])
    File flagger_bed = ugc_wgw_flagger_hmm.prediction_bed
    File flagger_summary = ugc_wgw_flagger_hmm.summary_tsv
    Array[File]? trio_eval = ugc_wgw_yak_trioeval.trioeval
    File? read_alignment_bam = kept_bam
    File? read_alignment_bam_index = kept_bam_index
    Boolean trio = is_trio
    String? hap1_parent = father
    String? hap2_parent = mother
    String ugc_wgw_workflow_version = ugc_wgw_version
    File ugc_wgw_manifest = ugc_wgw_manifest_write.ugc_wgw_manifest
  }
}
