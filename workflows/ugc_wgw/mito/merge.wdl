version 1.0

# ugc-pacbio-wgw — the cohort's mitochondrial call sets: bcftools merge of the members' collapsed chrM VCFs, and the haplogroups gathered
# origin: none
# see docs/DESIGN.md §5, §7, §18 (2026-11-05); docs/ENTRYPOINTS.md §7

import "../../../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/structs.wdl"

task ugc_wgw_mito_merge {
  meta {
    description: "bcftools merge of the members' collapsed chrM VCFs (one contig, no scatter), a member without a record at a site set to 0/0 (mitorsaw reports every site where it found an alternate allele, so absence is reference), header lines added, the members' per-sample haplotype-count and abundance header lines dropped, indexed. One member passes through."
    outputs: {
      vcf: {
        description: "<out_prefix>.vcf.gz"
      }
      vcf_index: {
        description: "Index"
      }
    }
  }

  parameter_meta {
    vcfs: {
      description: "The members' collapsed chrM VCFs (ugc_wgw_mito_collapse), cohort order; distinct sample names"
    }
    vcf_indices: {
      description: "Their indices"
    }
    header_lines: {
      description: "##key=value lines added to the header (provenance, the collapse rule, the 0/0 rule)"
    }
    out_prefix: {
      description: "Output prefix"
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
    Array[String] header_lines
    String out_prefix
    Int threads = 2
    Int mem_gb = 4
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(vcfs, "GB") * 2 + 5)

  command <<<
    set -euo pipefail
    export LC_ALL=C

    bcftools --version | head -1 >&2

    mkdir -p in
    : > list.txt
    while read -r f; do
      ln --symbolic "$f" "in/$(basename "$f")"
      echo "in/$(basename "$f")" >> list.txt
    done < "~{write_lines(vcfs)}"
    while read -r f; do
      ln --symbolic "$f" "in/$(basename "$f")"
    done < "~{write_lines(vcf_indices)}"
    ulimit -Sn 65536 || true

    if [ "$(wc -l < list.txt)" -eq 1 ]; then
      bcftools view --output-type u "$(cat list.txt)"
    else
      bcftools merge --threads ~{threads} --file-list list.txt --missing-to-ref --output-type u
    fi \
    | bcftools annotate --header-lines "~{write_lines(header_lines)}" --output-type b --output merged.bcf -
    # the members' own ##ugc_wgw_mito_haplotypes and ##ugc_wgw_mito_abundances lines describe one member each
    # and bcftools merge keeps the distinct ones without saying whose: dropped (they stay in each member's VCF
    # and in the haplogroups table); the collapse rule, the SV size and the mitorsaw version lines stay
    {
      bcftools view --header-only merged.bcf | grep -v -E '^##ugc_wgw_mito_(haplotypes|abundances)='
      bcftools view --no-header merged.bcf
    } | bcftools view --output-type z --output "~{out_prefix}.vcf.gz" -
    bcftools index --threads ~{threads} --tbi "~{out_prefix}.vcf.gz"
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

task ugc_wgw_mito_haplogroups_gather {
  meta {
    description: "One table of every member's mitochondrial haplogroups (ugc_wgw_haplogrep per member: a row per mitorsaw haplotype), with whether the classification was made in this stage (a member whose run predates it) or taken from the member's own run."
    outputs: {
      table: {
        description: "<out_prefix>.mt_haplogroups.tsv: sample_id, haplotype, abundance, haplogroup, quality, range, found_polys, remaining_polys, made_in_stage"
      }
    }
  }

  parameter_meta {
    sample_ids: {
      description: "The members with a haplogroup table, cohort order"
    }
    haplogroup_tsvs: {
      description: "Their tables (ugc_wgw_haplogrep), aligned with sample_ids"
    }
    made_in_stage: {
      description: "yes or no per member, aligned with sample_ids"
    }
    out_prefix: {
      description: "Output prefix"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    Array[String] sample_ids
    Array[File] haplogroup_tsvs
    Array[String] made_in_stage
    String out_prefix
    RuntimeAttributes runtime_attributes
  }

  command <<<
    set -euo pipefail

    python3 - "~{out_prefix}" "~{write_lines(sample_ids)}" "~{write_lines(haplogroup_tsvs)}" "~{write_lines(made_in_stage)}" <<'EOF'
    import sys
    prefix, ids_f, tsvs_f, made_f = sys.argv[1:5]
    read = lambda p: [l.rstrip("\n") for l in open(p) if l.strip()]
    ids, tsvs, made = read(ids_f), read(tsvs_f), read(made_f)
    with open(prefix + ".mt_haplogroups.tsv", "w") as out:
        out.write("sample_id\thaplotype\tabundance\thaplogroup\tquality\trange\tfound_polys\tremaining_polys\tmade_in_stage\n")
        for sid, path, flag in zip(ids, tsvs, made):
            with open(path) as fh:
                header = fh.readline().rstrip("\n").split("\t")
                for line in fh:
                    row = dict(zip(header, line.rstrip("\n").split("\t")))
                    out.write("\t".join([sid] + [row.get(c, "") for c in ("haplotype", "abundance", "haplogroup", "quality", "range", "found_polys", "remaining_polys")] + [flag]) + "\n")
    print(f"{len(ids)} members", file=sys.stderr)
    EOF
  >>>

  output {
    File table = "~{out_prefix}.mt_haplogroups.tsv"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:03cb3c01937eccc907f8ad71c87b258581504572205fe3f31a657e318f3564ae"  # pb_wdl_base:build4
    cpu: 1
    memory: "1 GiB"
    disk: "5 GB"
    disks: "local-disk 5 HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
    cpuPlatform: runtime_attributes.cpuPlatform
  }
}
