version 1.0
# ugc-pacbio-wgw — HMM-Flagger on one sample's diploid assembly: the HiFi reads aligned back to hap1+hap2, the alignment cleaned, coverage modelled per base into Err/Dup/Hap/Col regions
# origin: Flagger v1.2.0 (mobinasri/flagger) wdls/workflows/hmm_flagger_end_to_end_with_mapping.wdl, tasks/alignment/long_read_aligner.wdl, tasks/coverage/bam_coverage.wdl, tasks/hmm_flagger/hmm_flagger.wdl and ext/secphase correct_bam: commands ported into three tasks around the upstream images, not vendored (Cromwell runtime keys, relative submodule imports)
# see docs/DESIGN.md §5, §13; docs/ENTRYPOINTS.md §8; DERIVED_FILES.md
import "../../../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/structs.wdl"

task ugc_wgw_flagger_align {
  meta {
    description: "The sample's HiFi reads (the FASTAs samtools fasta made) aligned to its own diploid assembly, hap1 and hap2 concatenated, with minimap2 -k 19 -x map-hifi --eqx --cs -Y -L -y -I8g (Flagger's mapping options; -I8g indexes the whole diploid assembly in one batch), sorted and indexed. An intermediate of the stage unless keep_read_alignment."
    outputs: {
      bam: {
        description: "<sample_id>.asm.reads.bam: reads aligned to hap1+hap2"
      },
      bam_index: {
        description: "BAM index"
      },
      dip_fasta: {
        description: "<sample_id>.asm.dip.fasta.gz: hap1 and hap2 concatenated (the gzip concatenation of the two bgzipped FASTAs)"
      }
    }
  }

  parameter_meta {
    sample_id: {
      description: "Sample ID; output prefix"
    }
    reads_fastas: {
      description: "The sample's HiFi reads as FASTA, one per movie"
    }
    hap_fastas: {
      description: "The two haplotype FASTAs (bgzipped)"
    }
    threads: {
      description: "Threads (minimap2 gets threads - 4, samtools sort 4)"
    }
    mem_gb: {
      description: "Memory (GB): the -I8g index of a 6 Gb diploid assembly plus the sort"
    }
    ugc_wgw_container_registry: {
      description: "Registry holding the ugc-built images (minimap2)"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    String sample_id
    Array[File] reads_fastas
    Array[File] hap_fastas
    Int threads = 32
    Int mem_gb = 64
    String ugc_wgw_container_registry
    RuntimeAttributes runtime_attributes
  }

  String out_prefix = "~{sample_id}.asm.reads"
  Int minimap2_threads = if threads > 4 then threads - 4 else 1
  Int disk_size = ceil(size(reads_fastas, "GB") * 2 + size(hap_fastas, "GB") * 4 + 50)

  command <<<
    set -euo pipefail
    minimap2 --version >&2
    samtools --version | head -1 >&2

    cat ~{sep=" " hap_fastas} > "~{sample_id}.asm.dip.fasta.gz"
    mkdir -p TMP
    minimap2 \
      -t ~{minimap2_threads} \
      -k 19 \
      -a \
      -x map-hifi \
      --eqx \
      --cs \
      -Y \
      -L \
      -y \
      -I8g \
      "~{sample_id}.asm.dip.fasta.gz" \
      ~{sep=" " reads_fastas} \
    | samtools sort \
      -@ 4 \
      -T ./TMP/sort \
      -m 2G \
      -O BAM \
      -o "~{out_prefix}.bam"
    samtools index "~{out_prefix}.bam"
  >>>

  output {
    File bam = "~{out_prefix}.bam"
    File bam_index = "~{out_prefix}.bam.bai"
    File dip_fasta = "~{sample_id}.asm.dip.fasta.gz"
  }

  runtime {
    docker: "~{ugc_wgw_container_registry}/minimap2@sha256:00dd2dfbb7f395672eb0e71317899f913cd64af80c02864655dc9d848202604a"  # 2.31 + samtools 1.24 (containers/minimap2)
    cpu: threads
    memory: "~{mem_gb} GiB"
    disk: "~{disk_size} GB"
    disks: "local-disk ~{disk_size} HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
    cpuPlatform: runtime_attributes.cpuPlatform
  }
}

task ugc_wgw_flagger_correct_bam {
  meta {
    description: "secphase's correct_bam on the read alignment: primary alignments only, reads of 5 kb or more, alignments of 5 kb or more, gap-compressed divergence at most 0.1 (Flagger's correctBamOptions); indexed."
    outputs: {
      corrected_bam: {
        description: "<sample_id>.asm.reads.corrected.bam"
      },
      corrected_bam_index: {
        description: "BAM index"
      }
    }
  }

  parameter_meta {
    sample_id: {
      description: "Sample ID; output prefix"
    }
    bam: {
      description: "Reads aligned to the diploid assembly"
    }
    bam_index: {
      description: "Its index"
    }
    min_read_length: {
      description: "correct_bam --minReadLen"
    }
    min_alignment_length: {
      description: "correct_bam --minAlignmentLen"
    }
    max_divergence: {
      description: "correct_bam --maxDiv (gap-compressed divergence, the de tag)"
    }
    threads: {
      description: "Threads (BAM I/O)"
    }
    mem_gb: {
      description: "Memory (GB)"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    String sample_id
    File bam
    File bam_index
    Int min_read_length = 5000
    Int min_alignment_length = 5000
    Float max_divergence = 0.1
    Int threads = 8
    Int mem_gb = 16
    RuntimeAttributes runtime_attributes
  }

  String out_prefix = "~{sample_id}.asm.reads.corrected"
  Int disk_size = ceil(size(bam, "GB") * 2.5 + 20)

  command <<<
    set -euo pipefail
    ln -s "~{bam}" in.bam
    ln -s "~{bam_index}" in.bam.bai
    correct_bam \
      --primaryOnly \
      --minReadLen ~{min_read_length} \
      --minAlignmentLen ~{min_alignment_length} \
      --maxDiv ~{max_divergence} \
      --threads ~{threads} \
      -i in.bam \
      -o "~{out_prefix}.bam"
    samtools index "~{out_prefix}.bam"
  >>>

  output {
    File corrected_bam = "~{out_prefix}.bam"
    File corrected_bam_index = "~{out_prefix}.bam.bai"
  }

  runtime {
    docker: "mobinasri/secphase@sha256:0c5b3f0dd0bf6351fbdf5638ff67b421df1ee2a610d89cdc2bdd6b0653075470"  # v0.4.4
    cpu: threads
    memory: "~{mem_gb} GiB"
    disk: "~{disk_size} GB"
    disks: "local-disk ~{disk_size} HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
    cpuPlatform: runtime_attributes.cpuPlatform
  }
}

task ugc_wgw_flagger_hmm {
  meta {
    description: "bam2cov then hmm_flagger --preset hifi on the corrected alignment: per-base coverage against a whole-assembly annotation (every contig minus its N runs, as Flagger's own task builds it), then the HMM's prediction of Err, Dup, Hap and Col regions and its summary table."
    outputs: {
      prediction_bed: {
        description: "<sample_id>.asm.flagger.bed: Flagger's final prediction, one label per interval (column 4: Err, Dup, Hap or Col), in assembly coordinates"
      },
      summary_tsv: {
        description: "<sample_id>.asm.flagger.summary.tsv: Flagger's prediction_summary_final.tsv (the base_level percentage row of whole_genome gives the Err, Dup, Hap and Col shares)"
      }
    }
  }

  parameter_meta {
    sample_id: {
      description: "Sample ID; output prefix"
    }
    bam: {
      description: "The corrected read alignment"
    }
    bam_index: {
      description: "Its index"
    }
    dip_fasta: {
      description: "The diploid assembly the reads were aligned to (hap1+hap2)"
    }
    mapq_threshold: {
      description: "bam2cov --mapqThreshold"
    }
    min_alignment_length: {
      description: "bam2cov --minAlignmentLength"
    }
    threads: {
      description: "Threads"
    }
    mem_gb: {
      description: "Memory (GB)"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    String sample_id
    File bam
    File bam_index
    File dip_fasta
    Int mapq_threshold = 10
    Int min_alignment_length = 5000
    Int threads = 16
    Int mem_gb = 32
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(bam, "GB") + size(dip_fasta, "GB") * 3 + 64)

  command <<<
    set -euo pipefail
    hmm_flagger --version >&2 || true
    ln -s "~{bam}" reads.bam
    ln -s "~{bam_index}" reads.bam.bai

    # the whole-assembly annotation: every contig minus its runs of N (Flagger's WHOLE_GENOME_DEFAULT)
    gunzip -c "~{dip_fasta}" > dip.fa
    samtools faidx dip.fa
    python3 /home/scripts/get_N_coords.py --inputFasta dip.fa > N_coords.bed
    awk '{print $1"\t0\t"$2}' dip.fa.fai | bedtools sort -i - | bedtools subtract -a - -b N_coords.bed > whole_genome.bed
    printf '{"whole_genome": "%s/whole_genome.bed"}\n' "$PWD" > annotations.json

    bam2cov \
      --bam reads.bam \
      --mapqThreshold ~{mapq_threshold} \
      --clipRatioThreshold 0.1 \
      --output reads.cov.gz \
      --annotationJson annotations.json \
      --threads ~{threads} \
      --format all \
      --baselineAnnotation whole_genome \
      --downsampleRate 1.0 \
      --minAlignmentLength ~{min_alignment_length}

    mkdir -p flagger
    hmm_flagger \
      --input reads.cov.gz \
      --preset hifi \
      --outputDir flagger \
      --maxHighMapqRatio 0.25 \
      --minHighMapqRatio 0.75 \
      --iterations 100 \
      --trackName hmm_flagger \
      --convergenceTol 0.001 \
      --labelNames Err,Dup,Hap,Col \
      --threads ~{threads}
    cp flagger/final_flagger_prediction.bed "~{sample_id}.asm.flagger.bed"
    cp flagger/prediction_summary_final.tsv "~{sample_id}.asm.flagger.summary.tsv"
    grep base_level "~{sample_id}.asm.flagger.summary.tsv" | grep whole_genome | grep percentage >&2 || true
    rm -f dip.fa reads.cov.gz
  >>>

  output {
    File prediction_bed = "~{sample_id}.asm.flagger.bed"
    File summary_tsv = "~{sample_id}.asm.flagger.summary.tsv"
  }

  runtime {
    docker: "mobinasri/flagger@sha256:3571bbe1c5415126edd0e8e646c853cc17f390bf1e2eb7c03b283ae54bc8262c"  # v1.2.0
    cpu: threads
    memory: "~{mem_gb} GiB"
    disk: "~{disk_size} GB"
    disks: "local-disk ~{disk_size} HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
    cpuPlatform: runtime_attributes.cpuPlatform
  }
}
