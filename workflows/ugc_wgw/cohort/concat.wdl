version 1.0

# ugc-pacbio-wgw — bcftools concat of region shards back into one indexed VCF
# origin: none
# see docs/DESIGN.md §7.1

import "../../../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/structs.wdl"

task ugc_wgw_bcftools_concat {
  meta {
    description: "Concatenate per-shard VCFs (given in genome order) into one bgzipped, tabix-indexed VCF"
    outputs: {
      vcf: {
        description: "Concatenated VCF"
      },
      vcf_index: {
        description: "Index for concatenated VCF"
      }
    }
  }

  parameter_meta {
    vcfs: {
      description: "Shard VCFs in genome order"
    }
    vcf_indices: {
      description: "Shard VCF indices"
    }
    out_prefix: {
      description: "Output prefix"
    }
    allow_overlaps: {
      description: "Use bcftools concat --allow-overlaps (index-driven merge of the sorted inputs) when shards can interleave, e.g. a genome-wide BND shard next to per-contig shards"
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
    Array[File] vcfs
    Array[File] vcf_indices
    String out_prefix
    Boolean allow_overlaps = false
    Int threads = 4
    Int mem_gb = 8
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(vcfs, "GB") * 2 + 20)

  command <<<
    set -euo pipefail

    bcftools --version >&2

    # colocate every VCF with its index, keeping the given order
    mkdir -p in
    : > list.txt
    while read -r f; do
      ln --symbolic "$f" "in/$(basename "$f")"
      echo "in/$(basename "$f")" >> list.txt
    done < "~{write_lines(vcfs)}"
    while read -r f; do
      ln --symbolic "$f" "in/$(basename "$f")"
    done < "~{write_lines(vcf_indices)}"

    bcftools concat \
      --threads ~{threads} \
      ~{true="--allow-overlaps" false="" allow_overlaps} \
      --file-list list.txt \
      --output-type z \
      --output "~{out_prefix}.vcf.gz"
    bcftools index \
      --threads ~{threads} \
      --tbi \
      "~{out_prefix}.vcf.gz"
  >>>

  output {
    File vcf = "~{out_prefix}.vcf.gz"
    File vcf_index = "~{out_prefix}.vcf.gz.tbi"
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
