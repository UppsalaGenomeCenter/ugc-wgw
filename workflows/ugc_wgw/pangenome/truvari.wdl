version 1.0
# ugc-pacbio-wgw — the graph VCF split to biallelic records, then truvari bench of its structural variants against the cohort's read-based SV merge
# origin: none
# see docs/DESIGN.md §5, §13; docs/ENTRYPOINTS.md §8.2
import "../../../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/structs.wdl"

task ugc_wgw_vcf_split_multiallelic {
  meta {
    description: "bcftools norm -m -any: every multiallelic record of the graph VCF split into biallelic ones (truvari bench refuses multiallelic records; the primary VCFs keep them as Cactus wrote them), and the INFO/SVTYPE and INFO/SVLEN definitions added to the header (truvari sets them on the copies it decomposes when matching a resolved variant against a breakend of the base)."
    outputs: {
      biallelic_vcf: {
        description: "<out_prefix>.biallelic.vcf.gz"
      },
      biallelic_vcf_index: {
        description: "Tabix index"
      }
    }
  }

  parameter_meta {
    vcf: {
      description: "The VCF to split"
    }
    vcf_index: {
      description: "Its index"
    }
    out_prefix: {
      description: "Output prefix"
    }
    threads: {
      description: "Threads"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    File vcf
    File vcf_index
    String out_prefix
    Int threads = 2
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(vcf, "GB") * 3 + 10)

  command <<<
    set -euo pipefail
    bcftools --version | head -1 >&2
    printf '%s\n' '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant (declared for truvari; set by its decomposition only)">' \
      '##INFO=<ID=SVLEN,Number=A,Type=Integer,Description="Length of structural variant (declared for truvari; set by its decomposition only)">' > header.txt
    bcftools norm --threads ~{threads} --multiallelics -any --output-type u "~{vcf}" \
      | bcftools annotate --header-lines header.txt --output-type z --output "~{out_prefix}.biallelic.vcf.gz" -
    tabix --preset vcf "~{out_prefix}.biallelic.vcf.gz"
  >>>

  output {
    File biallelic_vcf = "~{out_prefix}.biallelic.vcf.gz"
    File biallelic_vcf_index = "~{out_prefix}.biallelic.vcf.gz.tbi"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:03cb3c01937eccc907f8ad71c87b258581504572205fe3f31a657e318f3564ae"  # pb_wdl_base:build4
    cpu: threads
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

task ugc_wgw_truvari_bench {
  meta {
    description: "truvari bench with the cohort_merge SV VCF as the base and the pangenome's vcfwave VCF (split to biallelic records) as the comparison, variants of 50 bp or more, sequence comparison off (the merged SV VCF carries symbolic alleles): recall, precision and F1 of the graph's SVs against the read-based ones, as truvari's summary.json, plus the tp/fp/fn VCFs tarred."
    outputs: {
      summary: {
        description: "<out_prefix>.truvari.json: truvari bench's summary.json"
      },
      details: {
        description: "<out_prefix>.truvari.tar.gz: truvari's output directory"
      }
    }
  }

  parameter_meta {
    base_vcf: {
      description: "The cohort_merge SV VCF (truth-side)"
    }
    base_vcf_index: {
      description: "Its index"
    }
    comp_vcf: {
      description: "The pangenome's vcfwave VCF, biallelic"
    }
    comp_vcf_index: {
      description: "Its index"
    }
    sizemin: {
      description: "truvari --sizemin (and --sizefilt): smallest variant compared"
    }
    out_prefix: {
      description: "Output prefix"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    File base_vcf
    File base_vcf_index
    File comp_vcf
    File comp_vcf_index
    Int sizemin = 50
    String out_prefix
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil((size(base_vcf, "GB") + size(comp_vcf, "GB")) * 4 + 10)

  command <<<
    set -euo pipefail
    truvari version >&2
    ln -s "~{base_vcf}" base.vcf.gz; ln -s "~{base_vcf_index}" "base.vcf.gz.$(basename "~{base_vcf_index}" | sed 's/.*\.//')"
    ln -s "~{comp_vcf}" comp.vcf.gz; ln -s "~{comp_vcf_index}" "comp.vcf.gz.$(basename "~{comp_vcf_index}" | sed 's/.*\.//')"
    truvari bench -b base.vcf.gz -c comp.vcf.gz -o truvari --sizemin ~{sizemin} --sizefilt ~{sizemin} --pctseq 0
    cp truvari/summary.json "~{out_prefix}.truvari.json"
    tar -czf "~{out_prefix}.truvari.tar.gz" truvari
    python3 -c 'import json,sys; d=json.load(open(sys.argv[1])); print({k: d.get(k) for k in ("TP-base","TP-comp","FP","FN","precision","recall","f1")}, file=sys.stderr)' "~{out_prefix}.truvari.json"
  >>>

  output {
    File summary = "~{out_prefix}.truvari.json"
    File details = "~{out_prefix}.truvari.tar.gz"
  }

  runtime {
    docker: "quay.io/biocontainers/truvari@sha256:28a193205760fad402674f6011425cee09f45085ab82d57f5ff3df4bc59c1264"  # 5.5.0--pyhdfd78af_0
    cpu: 2
    memory: "8 GiB"
    disk: "~{disk_size} GB"
    disks: "local-disk ~{disk_size} HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
    cpuPlatform: runtime_attributes.cpuPlatform
  }
}
