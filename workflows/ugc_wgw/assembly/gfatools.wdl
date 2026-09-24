version 1.0

# ugc-pacbio-wgw — GFA to bgzipped FASTA plus calN50 assembly statistics
# origin: HiFi-human-assembly-WDL v1.0.2 (025b3910) workflows/assemble_genome/assemble_genome.wdl; ported, not vendored
# see docs/DESIGN.md §5, §13; docs/ENTRYPOINTS.md §8

import "../../../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/structs.wdl"

task ugc_wgw_gfatools_gfa2fa {
  meta {
    description: "Convert a hifiasm contig GFA to a bgzipped FASTA and compute size/NG50 statistics with calN50 (NG relative to 3.1 Gb)"
    outputs: {
      zipped_fasta: {
        description: "Contigs in bgzipped FASTA"
      },
      stats: {
        description: "calN50 statistics"
      }
    }
  }

  parameter_meta {
    gfa: {
      description: "Contig graph (GFA)"
    }
    threads: {
      description: "Threads (bgzip)"
    }
    mem_gb: {
      description: "Memory (GB)"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    File gfa
    Int threads = 2
    Int mem_gb = 4
    RuntimeAttributes runtime_attributes
  }

  String out_prefix = basename(gfa, ".gfa")
  Int disk_size = ceil(size(gfa, "GB") * 3 + 20)

  command <<<
    set -euo pipefail

    gfatools version >&2
    bgzip --version >&2
    k8 /opt/calN50/calN50.js -v >&2

    ln --symbolic --verbose "~{gfa}" .

    gfatools gfa2fa "~{basename(gfa)}" > "~{out_prefix}.fasta"

    bgzip \
      --threads ~{threads} \
      --stdout \
      "~{out_prefix}.fasta" \
    > "~{out_prefix}.fasta.gz"

    k8 /opt/calN50/calN50.js \
      -L3.1g \
      "~{out_prefix}.fasta.gz" \
    > "~{out_prefix}.fasta.stats.txt"
  >>>

  output {
    File zipped_fasta = "~{out_prefix}.fasta.gz"
    File stats = "~{out_prefix}.fasta.stats.txt"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/gfatools@sha256:0ec85987c043e84fd492e9e0908261e365b70a6d4eaf797c39a9c1da695a10ce"  # 0.5_34e0fcf_build2 (build3 ships a dangling gfatools symlink)
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
