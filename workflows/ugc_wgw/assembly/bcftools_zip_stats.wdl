version 1.0

# ugc-pacbio-wgw — bgzip + tabix a VCF and compute bcftools stats against the reference (one task, pb_wdl_base)
# origin: HiFi-human-assembly-WDL v1.0.2 (025b3910), wdl-common e37b3274, wdl/tasks/zip_index_vcf.wdl and bcftools_stats.wdl; ported, not vendored
# see docs/DESIGN.md §5; docs/ENTRYPOINTS.md §8

import "../../../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/structs.wdl"

task ugc_wgw_bcftools_zip_stats {
  meta {
    description: "Compress and index a plain VCF, then bcftools stats for the sample against the reference; folds upstream's zip_index_vcf and bcftools_stats into one job"
    outputs: {
      zipped_vcf: {
        description: "bgzipped VCF"
      },
      zipped_vcf_index: {
        description: "Tabix index"
      },
      stats: {
        description: "bcftools stats"
      }
    }
  }

  parameter_meta {
    vcf: {
      description: "Uncompressed VCF"
    }
    ref_fasta: {
      description: "Reference FASTA"
    }
    ref_index: {
      description: "Reference FASTA index"
    }
    sample_id: {
      description: "Sample to report stats for"
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
    File vcf
    File ref_fasta
    File ref_index
    String sample_id
    Int threads = 4
    Int mem_gb = 4
    RuntimeAttributes runtime_attributes
  }

  String out_name = basename(vcf)
  Int disk_size = ceil((size(vcf, "GB") + size(ref_fasta, "GB")) * 2 + 20)

  command <<<
    set -euo pipefail

    bgzip --version >&2
    tabix --version >&2
    bcftools --version >&2

    ln --symbolic --verbose "~{vcf}" .
    ln --symbolic --verbose "~{ref_fasta}" .
    ln --symbolic --verbose "~{ref_index}" .

    bgzip \
      --threads ~{threads} \
      --stdout \
      "~{out_name}" \
    > "~{out_name}.gz"
    tabix --preset vcf "~{out_name}.gz"

    bcftools stats \
      ~{if threads > 1
        then "--threads '" + (threads - 1) + "'"
        else ""
      } \
      --samples "~{sample_id}" \
      --fasta-ref "~{basename(ref_fasta)}" \
      "~{out_name}.gz" \
    > "~{out_name}.stats.txt"
  >>>

  output {
    File zipped_vcf = "~{out_name}.gz"
    File zipped_vcf_index = "~{out_name}.gz.tbi"
    File stats = "~{out_name}.stats.txt"
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
