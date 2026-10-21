version 1.0
# ugc-pacbio-wgw — per-member genotype statistics of the cohort's graph VCF: bcftools stats per chromosome, summed per sample
# origin: none
# see docs/DESIGN.md §5, §13; docs/ENTRYPOINTS.md §8.2
import "../../../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/structs.wdl"

task ugc_wgw_pangenome_member_stats {
  meta {
    description: "bcftools stats --samples - of the vcfwave VCF inside one chromosome's BED: the PSC rows give every member's hom-ref, hom-alt and het genotypes, singletons and missing genotypes, the SN row the sites. One shard per requested chromosome keeps the per-member pass at genome scale within minutes."
    outputs: {
      stats: {
        description: "<out_prefix>.stats.txt: bcftools stats of the shard"
      }
    }
  }

  parameter_meta {
    vcf: {
      description: "The vcfwave VCF with genotypes"
    }
    vcf_index: {
      description: "Its index"
    }
    regions_bed: {
      description: "The chromosome's BED"
    }
    out_prefix: {
      description: "Output prefix"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    File vcf
    File vcf_index
    File regions_bed
    String out_prefix
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(vcf, "GB") + 10)

  command <<<
    set -euo pipefail
    bcftools --version | head -1 >&2
    bcftools stats --samples - --regions-file "~{regions_bed}" "~{vcf}" > "~{out_prefix}.stats.txt"
    grep -c "^PSC" "~{out_prefix}.stats.txt" | sed 's/$/ members/' >&2
  >>>

  output {
    File stats = "~{out_prefix}.stats.txt"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:03cb3c01937eccc907f8ad71c87b258581504572205fe3f31a657e318f3564ae"  # pb_wdl_base:build4
    cpu: 1
    memory: "4 GiB"
    disk: "~{disk_size} GB"
    disks: "local-disk ~{disk_size} HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
    cpuPlatform: runtime_attributes.cpuPlatform
  }
}

task ugc_wgw_pangenome_member_stats_gather {
  meta {
    description: "The per-chromosome bcftools stats of the graph VCF summed per member: hom-ref, hom-alt and het genotypes, singletons (alleles seen in one member only), missing genotypes, the sites counted, and the call rate (1 - missing / sites: a missing genotype means no assembled sequence at the site)."
    outputs: {
      members: {
        description: "<out_prefix>.members.tsv: sample_id, n_ref_hom, n_nonref_hom, n_hets, n_singletons, n_missing, n_sites, call_rate"
      }
    }
  }

  parameter_meta {
    stats: {
      description: "The shards' bcftools stats files"
    }
    out_prefix: {
      description: "Output prefix"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    Array[File] stats
    String out_prefix
    RuntimeAttributes runtime_attributes
  }

  command <<<
    set -euo pipefail
    python3 - "~{out_prefix}" ~{sep=" " stats} <<'PY'
import collections, sys
out_prefix, paths = sys.argv[1], sys.argv[2:]
tot = collections.defaultdict(lambda: collections.Counter())
order = []
sites = 0
for p in paths:
    for line in open(p):
        f = line.rstrip("\n").split("\t")
        if f[0] == "SN" and len(f) >= 4 and f[2].startswith("number of records"):
            sites += int(f[3])
        elif f[0] == "PSC" and len(f) >= 14:
            s = f[2]
            if s not in order:
                order.append(s)
            c = tot[s]
            c["n_ref_hom"] += int(f[3]); c["n_nonref_hom"] += int(f[4]); c["n_hets"] += int(f[5])
            c["n_singletons"] += int(f[10]); c["n_missing"] += int(f[13])
with open(out_prefix + ".members.tsv", "w") as out:
    out.write("sample_id\tn_ref_hom\tn_nonref_hom\tn_hets\tn_singletons\tn_missing\tn_sites\tcall_rate\n")
    for s in order:
        c = tot[s]
        rate = 1 - c["n_missing"] / sites if sites else 0.0
        out.write("\t".join(str(x) for x in [s, c["n_ref_hom"], c["n_nonref_hom"], c["n_hets"], c["n_singletons"], c["n_missing"], sites, f"{rate:.4f}"]) + "\n")
print(f"{len(order)} members over {sites} sites", file=sys.stderr)
PY
  >>>

  output {
    File members = "~{out_prefix}.members.tsv"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:03cb3c01937eccc907f8ad71c87b258581504572205fe3f31a657e318f3564ae"  # pb_wdl_base:build4
    cpu: 1
    memory: "2 GiB"
    disk: "10 GB"
    disks: "local-disk 10 HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
    cpuPlatform: runtime_attributes.cpuPlatform
  }
}
