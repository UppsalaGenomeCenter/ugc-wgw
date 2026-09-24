version 1.0

# ugc-pacbio-wgw — yak k-mer table of a parent's reads, for hifiasm trio binning
# origin: HiFi-human-assembly-WDL v1.0.2 (025b3910) workflows/de_novo_assembly_trio/de_novo_assembly_trio.wdl; ported, not vendored
# see docs/DESIGN.md §5, §13; docs/ENTRYPOINTS.md §8

import "../../../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/structs.wdl"

task ugc_wgw_yak_count {
  meta {
    description: "yak count over one parent's HiFi reads. -b0 keeps singleton k-mers (low parental depth), -b37 uses a bloom filter (recommended for human at normal depth)"
    outputs: {
      yak: {
        description: "k-mer table"
      }
    }
  }

  parameter_meta {
    sample_id: {
      description: "Parent sample ID (names the table)"
    }
    reads_fastas: {
      description: "Parent reads in FASTA"
    }
    yak_params: {
      description: "yak count parameters, e.g. -b0 or -b37"
    }
    threads: {
      description: "Threads"
    }
    mem_gb: {
      description: "Memory (GB): 70 without bloom filter (<= 30x), 50 with (<= 50x), for 24 threads"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    String sample_id
    Array[File] reads_fastas
    String yak_params
    Int threads = 24
    Int mem_gb = 70
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(reads_fastas, "GB") * 2 + 20)

  command <<<
    set -euo pipefail

    yak version >&2

    # shellcheck disable=SC2086
    yak count \
      -t ~{threads} \
      -o "~{sample_id}.yak" \
      ~{yak_params} \
      ~{sep=" " reads_fastas}
  >>>

  output {
    File yak = "~{sample_id}.yak"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/yak@sha256:7809845f03cb04e3a686178b006c42b0f61768797543d4cb46f62940ba44ffb7"  # 0.1_build3 (build4 ships a dangling yak symlink)
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
