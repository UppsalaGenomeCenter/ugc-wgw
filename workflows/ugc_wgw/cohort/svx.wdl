version 1.0

# ugc-pacbio-wgw — svx merge of per-sample sawfish VCFs, scoped to one shard (or genome-wide for BND)
# origin: none
# see docs/DESIGN.md §7.3

import "../../../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/structs.wdl"

task ugc_wgw_svx_merge {
  meta {
    description: "Merge sawfish SV/CNV VCFs across samples with svx. Scoped with --target-positions (POS-based) and --svtype; orphaned breakends are dropped silently under scoping (svx 0.6.0), so BND must be merged in one genome-wide call"
    outputs: {
      vcf: {
        description: "Merged, sorted VCF (svx writes the index with --sort-output)"
      },
      vcf_index: {
        description: "Index for merged VCF"
      }
    }
  }

  parameter_meta {
    vcfs: {
      description: "sawfish VCFs (single- or multi-sample, bgzipped, indexed, unique sample names)"
    }
    vcf_indices: {
      description: "VCF indices"
    }
    out_prefix: {
      description: "Output prefix"
    }
    target_positions: {
      description: "svx --target-positions value, e.g. (chr1:1-248956422),(chr1_KI270706v1_random:1-175055); omit for genome-wide"
    }
    svtypes: {
      description: "svx --svtype list: ALL, or a subset of INS,DEL,INV,DUP,BND,CNV"
    }
    trs_bed: {
      description: "TRGT-style repeat catalog BED for TR containment (svx --trs); flags must already be stripped"
    }
    min_supp: {
      description: "svx --min-supp: minimum number of supporting samples for a merged record"
    }
    threads: {
      description: "svx --threads"
    }
    mem_gb: {
      description: "Memory (GB); the entrypoint scales it with sample count"
    }
    ugc_wgw_container_registry: {
      description: "Registry holding the ugc-built images"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    Array[File] vcfs
    Array[File] vcf_indices
    String out_prefix
    String? target_positions
    String svtypes = "ALL"
    File? trs_bed
    Int min_supp = 1
    Int threads = 8
    Int mem_gb = 16
    String ugc_wgw_container_registry
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(vcfs, "GB") * 3 + 20)

  command <<<
    set -euo pipefail

    svx --version >&2

    mkdir -p in sort_tmp
    : > list.txt
    while read -r f; do
      ln --symbolic "$f" "in/$(basename "$f")"
      echo "in/$(basename "$f")" >> list.txt
    done < "~{write_lines(vcfs)}"
    while read -r f; do
      ln --symbolic "$f" "in/$(basename "$f")"
    done < "~{write_lines(vcf_indices)}"

    # shellcheck disable=SC2086
    svx merge \
      --vcf-list list.txt \
      ~{if length(vcfs) == 1 then "--force-single" else ""} \
      ~{if defined(target_positions) then "--target-positions '" + select_first([target_positions]) + "'" else ""} \
      --svtype "~{svtypes}" \
      ~{"--trs " + trs_bed} \
      --min-supp ~{min_supp} \
      --threads ~{threads} \
      --sort-output \
      --sort-tmp-dir sort_tmp \
      --output-type z \
      --output "~{out_prefix}.vcf.gz"

    rm -rf sort_tmp
  >>>

  output {
    File vcf = "~{out_prefix}.vcf.gz"
    File vcf_index = "~{out_prefix}.vcf.gz.tbi"
  }

  runtime {
    docker: "~{ugc_wgw_container_registry}/svx@sha256:5f5577e979e4a1e629c9630ebb647e3efe8ebe3978cd20dfb81a712807007a94"  # 0.6.0
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
