version 1.0

# ugc-pacbio-wgw — bcftools merge of per-sample VCFs restricted to one scatter region (site union, as upstream family.wdl)
# origin: none
# see docs/DESIGN.md §7.4, §7.5

import "../../../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/structs.wdl"

task ugc_wgw_bcftools_merge {
  meta {
    description: "bcftools merge of many single-sample VCFs within the regions of one shard; mirrors vendor/.../tasks/bcftools.wdl bcftools_merge plus --regions-file"
    outputs: {
      vcf: {
        description: "Merged VCF for the shard"
      },
      vcf_index: {
        description: "Index for merged VCF"
      }
    }
  }

  parameter_meta {
    vcfs: {
      description: "Per-sample VCFs (bgzipped, indexed, unique sample names)"
    }
    vcf_indices: {
      description: "Per-sample VCF indices"
    }
    regions_bed: {
      description: "BED of the shard (bcftools --regions-file)"
    }
    out_prefix: {
      description: "Output prefix"
    }
    threads: {
      description: "Threads"
    }
    mem_gb: {
      description: "Memory (GB); one reader per input, scale with sample count"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    Array[File] vcfs
    Array[File] vcf_indices
    File regions_bed
    String out_prefix
    Int threads = 4
    Int mem_gb = 8
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(vcfs, "GB") * 2 + 20)

  command <<<
    set -euo pipefail

    bcftools --version >&2

    # one open file per input VCF plus its index
    ulimit -Sn 65536 || true

    mkdir -p in
    : > list.txt
    while read -r f; do
      ln --symbolic "$f" "in/$(basename "$f")"
      echo "in/$(basename "$f")" >> list.txt
    done < "~{write_lines(vcfs)}"
    while read -r f; do
      ln --symbolic "$f" "in/$(basename "$f")"
    done < "~{write_lines(vcf_indices)}"

    if [ "~{length(vcfs)}" -eq 1 ]; then
      # bcftools merge needs two or more inputs; a single-sample cohort is passed through
      bcftools view \
        --threads ~{threads} \
        --regions-file "~{regions_bed}" \
        --output-type z \
        --output "~{out_prefix}.vcf.gz" \
        "$(head -n 1 list.txt)"
    else
      bcftools merge \
        --threads ~{threads} \
        --regions-file "~{regions_bed}" \
        --file-list list.txt \
        --output-type z \
        --output "~{out_prefix}.vcf.gz"
    fi
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
