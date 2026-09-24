version 1.0

# ugc-pacbio-wgw — trgt merge scoped to one shard's contigs, and trgt-lps on the concatenated cohort TRGT VCF
# origin: none
# see docs/DESIGN.md §7.6

import "../../../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/structs.wdl"

task ugc_wgw_trgt_merge {
  meta {
    description: "trgt merge of per-sample TRGT VCFs restricted to the contigs of one shard; mirrors vendor/.../tasks/trgt.wdl trgt_merge with --vcf-list and --contig"
    outputs: {
      vcf: {
        description: "Merged TRGT VCF for the shard"
      },
      vcf_index: {
        description: "Index for merged TRGT VCF"
      }
    }
  }

  parameter_meta {
    vcfs: {
      description: "Per-sample TRGT VCFs (sorted, bgzipped, indexed, same catalog, unique sample names)"
    }
    vcf_indices: {
      description: "TRGT VCF indices"
    }
    ref_fasta: {
      description: "Reference FASTA"
    }
    ref_index: {
      description: "Reference FASTA index"
    }
    contigs: {
      description: "Comma-separated contigs of the shard (trgt merge --contig); shards must be whole contigs"
    }
    out_prefix: {
      description: "Output prefix"
    }
    no_index: {
      description: "Stream inputs without loading their indices (lower memory, reads every input fully)"
    }
    threads: {
      description: "trgt merge --threads ((de)compression only)"
    }
    mem_gb: {
      description: "Memory (GB); the entrypoint scales it with sample count"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    Array[File] vcfs
    Array[File] vcf_indices
    File ref_fasta
    File ref_index
    String contigs
    String out_prefix
    Boolean no_index = false
    Int threads = 2
    Int mem_gb = 8
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil((size(vcfs, "GB") + size(ref_fasta, "GB")) * 2 + 20)

  command <<<
    set -euo pipefail

    trgt --version >&2

    ln --symbolic --verbose "~{ref_fasta}" .
    ln --symbolic --verbose "~{ref_index}" .

    mkdir -p in
    : > list.txt
    while read -r f; do
      ln --symbolic "$f" "in/$(basename "$f")"
      echo "in/$(basename "$f")" >> list.txt
    done < "~{write_lines(vcfs)}"
    while read -r f; do
      ln --symbolic "$f" "in/$(basename "$f")"
    done < "~{write_lines(vcf_indices)}"

    # --quit-on-errors: a reference-allele mismatch means the inputs were genotyped with
    # different catalogs, which must fail the cohort rather than skip records
    # shellcheck disable=SC2086
    trgt merge \
      --vcf-list list.txt \
      ~{if length(vcfs) == 1 then "--force-single" else ""} \
      --contig "~{contigs}" \
      --genome "~{basename(ref_fasta)}" \
      --threads ~{threads} \
      --quit-on-errors \
      ~{true="--no-index" false="" no_index} \
      --output "~{out_prefix}.vcf.gz" \
      --write-index
  >>>

  output {
    File vcf = "~{out_prefix}.vcf.gz"
    File vcf_index = "~{out_prefix}.vcf.gz.tbi"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/trgt@sha256:648aee4a2c9d7371a48e454a7143861a242b853d81ff5453924cc0095d207824"  # 5.1.0_build2
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

task ugc_wgw_trgt_lps {
  meta {
    description: "Longest pure segment table (trgt-lps) from a multi-sample TRGT VCF"
    outputs: {
      lps: {
        description: "TSV: trid, motif, one column per sample with per-allele LPS values ('.' when missing)"
      }
    }
  }

  parameter_meta {
    vcf: {
      description: "Merged TRGT VCF"
    }
    vcf_index: {
      description: "Merged TRGT VCF index"
    }
    out_prefix: {
      description: "Output prefix"
    }
    threads: {
      description: "trgt-lps --threads"
    }
    mem_gb: {
      description: "Memory (GB)"
    }
    ugc_wgw_container_registry: {
      description: "Registry holding the ugc-built images"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    File vcf
    File vcf_index
    String out_prefix
    Int threads = 8
    Int mem_gb = 16
    String ugc_wgw_container_registry
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(vcf, "GB") * 2 + 20)

  command <<<
    set -euo pipefail

    trgt-lps --version >&2

    ln --symbolic --verbose "~{vcf}" .
    ln --symbolic --verbose "~{vcf_index}" .

    trgt-lps \
      --vcf "~{basename(vcf)}" \
      --threads ~{threads} \
    > "~{out_prefix}.lps.tsv"
  >>>

  output {
    File lps = "~{out_prefix}.lps.tsv"
  }

  runtime {
    docker: "~{ugc_wgw_container_registry}/trgt-lps@sha256:2a7b3c6c78a8d09983698d73f4fc140b89c15ae89d912203eed80b24a113cae7"  # 0.10.0
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
