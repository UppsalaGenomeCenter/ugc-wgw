version 1.0

# ugc-pacbio-wgw — align haplotype contigs to the reference with minimap2 and call variants from the alignment with paftools
# origin: HiFi-human-assembly-WDL v1.0.2 (025b3910) workflows/assemble_genome/assemble_genome.wdl; ported, not vendored
# see docs/DESIGN.md §5, §13; docs/ENTRYPOINTS.md §8

import "../../../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/structs.wdl"

task ugc_wgw_minimap2_align_asm {
  meta {
    description: "minimap2 -x asm5 alignment of one haplotype's contigs to the reference, sorted and indexed. Upstream's task lacked set -euo pipefail and located its outputs by glob; here the haplotype is an input and the output names are explicit"
    outputs: {
      bam: {
        description: "Aligned contigs, <sample_id>.<haplotype>.asm.<ref_name>.bam"
      },
      bam_index: {
        description: "BAM index"
      }
    }
  }

  parameter_meta {
    sample_id: {
      description: "Sample ID (read group)"
    }
    haplotype: {
      description: "hap1 or hap2"
    }
    query_fasta: {
      description: "Contigs (bgzipped FASTA)"
    }
    ref_fasta: {
      description: "Reference FASTA"
    }
    ref_name: {
      description: "Reference name"
    }
    threads: {
      description: "Threads (minimap2 gets threads - 4, samtools sort 3)"
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
    String haplotype
    File query_fasta
    File ref_fasta
    String ref_name
    Int threads = 16
    Int mem_gb = 128
    RuntimeAttributes runtime_attributes
  }

  String out_prefix = "~{sample_id}.~{haplotype}.asm.~{ref_name}"
  Int minimap2_threads = if threads > 4 then threads - 4 else 1
  Int disk_size = ceil((size(query_fasta, "GB") + size(ref_fasta, "GB")) * 2 + 20)

  command <<<
    set -euo pipefail

    minimap2 --version >&2
    samtools --version >&2

    ln --symbolic --verbose "~{ref_fasta}" .
    ln --symbolic --verbose "~{query_fasta}" .
    mkdir -p TMP

    minimap2 \
      -t ~{minimap2_threads} \
      -L \
      --secondary=no \
      --eqx \
      --cs \
      -a \
      -x asm5 \
      -R "@RG\\tID:~{sample_id}_hifiasm\\tSM:~{sample_id}" \
      "~{basename(ref_fasta)}" \
      "~{basename(query_fasta)}" \
    | samtools sort \
      -@ 3 \
      -T ./TMP/sort \
      -m 8G \
      -O BAM \
      -o "~{out_prefix}.bam"

    samtools index "~{out_prefix}.bam"
  >>>

  output {
    File bam = "~{out_prefix}.bam"
    File bam_index = "~{out_prefix}.bam.bai"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/align_hifiasm@sha256:0e8ad680b0e89376eb94fa8daa1a0269a4abe695ba39523a5c56a59d5c0e3953"  # 2.17_1.14_build2
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

task ugc_wgw_paftools_call {
  meta {
    description: "paftools.js call on one haplotype's contig alignment: variants from the cs tag and callable regions covered by exactly one contig; VCF against the reference"
    outputs: {
      vcf: {
        description: "paftools VCF (uncompressed)"
      }
    }
  }

  parameter_meta {
    bam: {
      description: "Contig alignment BAM"
    }
    bam_index: {
      description: "BAM index"
    }
    ref_fasta: {
      description: "Reference FASTA"
    }
    sample_id: {
      description: "Sample name for the VCF header"
    }
    threads: {
      description: "Threads (sort --parallel)"
    }
    mem_gb: {
      description: "Memory (GB); the whole-genome PAF is sorted in memory (sort -S mem_gb - 8)"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    File bam
    File bam_index
    File ref_fasta
    String sample_id
    Int threads = 4
    Int mem_gb = 32
    RuntimeAttributes runtime_attributes
  }

  String out_prefix = basename(bam, ".bam")
  Int sort_mem_gb = if mem_gb > 12 then mem_gb - 8 else 4
  Int disk_size = ceil((size(bam, "GB") + size(ref_fasta, "GB")) * 3 + 20)

  command <<<
    set -euo pipefail

    samtools --version >&2
    k8 /opt/minimap2-2.17/misc/paftools.js version >&2

    ln --symbolic --verbose "~{bam}" .
    ln --symbolic --verbose "~{bam_index}" .
    ln --symbolic --verbose "~{ref_fasta}" .
    mkdir -p TMP

    samtools view -h "~{basename(bam)}" \
    | k8 /opt/minimap2-2.17/misc/paftools.js sam2paf - \
    | sort -k6,6 -k8,8n -S ~{sort_mem_gb}G --parallel=~{threads} -T ./TMP \
    | k8 /opt/minimap2-2.17/misc/paftools.js call \
      -L5000 \
      -f "~{basename(ref_fasta)}" \
      -s "~{sample_id}" \
      - \
    > "~{out_prefix}.paftools.vcf"
  >>>

  output {
    File vcf = "~{out_prefix}.paftools.vcf"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/align_hifiasm@sha256:0e8ad680b0e89376eb94fa8daa1a0269a4abe695ba39523a5c56a59d5c0e3953"  # 2.17_1.14_build2
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
