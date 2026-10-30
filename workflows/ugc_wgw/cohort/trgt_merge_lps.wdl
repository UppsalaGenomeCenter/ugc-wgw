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
    description: "Longest pure segment table (trgt-lps) from a multi-sample TRGT VCF (one merge shard), rows in the VCF's record order"
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
    Int threads = 4
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
    > unordered.tsv

    # trgt-lps emits its rows in a thread-dependent order. Put them into the
    # VCF's record order (by TRID) so the table is reproducible and the
    # per-shard tables gather into the whole-cohort table byte for byte.
    zcat "~{basename(vcf)}" | { grep -v '^#' || true; } \
      | awk -F'\t' '{ if (match($8, /(^|;)TRID=[^;]*/)) { k = substr($8, RSTART, RLENGTH); sub(/^;?TRID=/, "", k); print k } else print $3 }' \
      > order.txt
    head -1 unordered.tsv > "~{out_prefix}.lps.tsv"
    awk -F'\t' -v OFS='\t' '
      FILENAME == ARGV[1] { ord[$1] = FNR; next }
      FNR == 1 { next }
      { if (!($1 in ord)) { print "trgt-lps row for a TRID not in the VCF: " $1 > "/dev/stderr"; exit 1 } print ord[$1], $0 }
    ' order.txt unordered.tsv \
      | sort -t "$(printf '\t')" -k1,1n | cut -f2- >> "~{out_prefix}.lps.tsv"
    echo "$(($(wc -l < "~{out_prefix}.lps.tsv") - 1)) loci with an LPS row of $(wc -l < order.txt) in the VCF" >&2
    rm -f unordered.tsv order.txt
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


task ugc_wgw_trgt_lps_gather {
  meta {
    description: "The per-shard trgt-lps tables as one: the first shard's header, then every shard's rows in shard order (one locus per row; trgt-lps computes per locus, so sharding by contig is exact). Replaces one trgt-lps over the whole cohort VCF, which took 13h52m for 864 samples (2026-10-30)"
    outputs: {
      lps: {
        description: "<out_prefix>.lps.tsv"
      }
    }
  }

  parameter_meta {
    shards: {
      description: "Per-shard trgt-lps TSVs, in shard order (the cohort VCF's contig order)"
    }
    out_prefix: {
      description: "Output prefix"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    Array[File] shards
    String out_prefix
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(shards, "GB") * 2 + 10)

  command <<<
    set -euo pipefail
    first=1
    while read -r f; do
      if [ "$first" = 1 ]; then
        head -1 "$f"
        first=0
      else
        # the headers must agree (same samples, same order): the merge used one sample list for every shard
        if ! cmp -s <(head -1 "$f") <(head -1 "$(head -1 '~{write_lines(shards)}')"); then
          echo "trgt-lps header of $f differs from the first shard's" >&2
          exit 1
        fi
      fi
      tail -n +2 "$f"
    done < '~{write_lines(shards)}' > "~{out_prefix}.lps.tsv"
    echo "$(($(wc -l < "~{out_prefix}.lps.tsv") - 1)) loci from $(wc -l < '~{write_lines(shards)}') shards" >&2
  >>>

  output {
    File lps = "~{out_prefix}.lps.tsv"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:03cb3c01937eccc907f8ad71c87b258581504572205fe3f31a657e318f3564ae"  # pb_wdl_base:build4
    cpu: 1
    memory: "2 GiB"
    disk: "~{disk_size} GB"
    disks: "local-disk ~{disk_size} HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
    cpuPlatform: runtime_attributes.cpuPlatform
  }
}

