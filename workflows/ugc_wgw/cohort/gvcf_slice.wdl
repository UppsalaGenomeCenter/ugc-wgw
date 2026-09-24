version 1.0

# ugc-pacbio-wgw — slice one sample's gVCF into the scatter regions (index-based), so GLnexus shards read only their slice
# origin: none
# see docs/DESIGN.md §7.2

import "../../../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/structs.wdl"

task ugc_wgw_bcftools_slice {
  meta {
    description: "bcftools view --regions-file per shard on one gVCF; GLnexus does not use gVCF indices, so this keeps total gVCF I/O at about one pass"
    outputs: {
      slices: {
        description: "gVCF slice per shard, in shard order"
      },
      slice_indices: {
        description: "Indices for the slices, same order"
      }
    }
  }

  parameter_meta {
    gvcf: {
      description: "Single-sample gVCF (bgzipped, indexed)"
    }
    gvcf_index: {
      description: "gVCF index"
    }
    shard_names: {
      description: "Shard names"
    }
    shard_beds: {
      description: "Shard BEDs, same order"
    }
    out_prefix: {
      description: "Output prefix (sample ID)"
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
    File gvcf
    File gvcf_index
    Array[String] shard_names
    Array[File] shard_beds
    String out_prefix
    Int threads = 2
    Int mem_gb = 4
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(gvcf, "GB") * 2 + 20)

  command <<<
    set -euo pipefail

    bcftools --version >&2

    ln --symbolic --verbose "~{gvcf}" .
    ln --symbolic --verbose "~{gvcf_index}" .

    : > slices.txt
    : > slice_indices.txt
    # --regions-overlap 1 keeps reference blocks that start before a shard but overlap it
    paste "~{write_lines(shard_names)}" "~{write_lines(shard_beds)}" \
    | while IFS=$'\t' read -r shard bed; do
      bcftools view \
        --threads ~{threads} \
        --regions-file "$bed" \
        --regions-overlap 1 \
        --output-type z \
        --output "~{out_prefix}.${shard}.g.vcf.gz" \
        "~{basename(gvcf)}"
      bcftools index \
        --threads ~{threads} \
        --tbi \
        "~{out_prefix}.${shard}.g.vcf.gz"
      echo "~{out_prefix}.${shard}.g.vcf.gz" >> slices.txt
      echo "~{out_prefix}.${shard}.g.vcf.gz.tbi" >> slice_indices.txt
    done
  >>>

  output {
    Array[File] slices = read_lines("slices.txt")
    Array[File] slice_indices = read_lines("slice_indices.txt")
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:03cb3c01937eccc907f8ad71c87b258581504572205fe3f31a657e318f3564ae"  # pb_wdl_base:build4
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
