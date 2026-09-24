version 1.0

# ugc-pacbio-wgw — hifiasm de novo assembly of one sample, optionally trio-binned with parental yak k-mer tables
# origin: HiFi-human-assembly-WDL v1.0.2 (025b3910) workflows/assemble_genome/assemble_genome.wdl; ported, not vendored
# see docs/DESIGN.md §5, §13; docs/ENTRYPOINTS.md §8

import "../../../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/structs.wdl"

task ugc_wgw_hifiasm_assemble {
  meta {
    description: "hifiasm assembly of HiFi reads; with father_yak and mother_yak the child's reads are trio-binned (hap1 = father, hap2 = mother). HiFi-only runs write <prefix>.bp.*, trio runs <prefix>.dip.*"
    outputs: {
      assembly_hap_gfas: {
        description: "Haplotype-resolved primary contig graphs, hap1 then hap2"
      },
      assembly_noseq_gfas: {
        description: "Sequence-less GFAs of the haplotype contigs and the unitig graphs"
      },
      assembly_lowQ_beds: {
        description: "Low-quality regions of the haplotype contigs and the unitig graphs"
      }
    }
  }

  parameter_meta {
    sample_id: {
      description: "Sample ID; output prefix is <sample_id>.asm"
    }
    reads_fastas: {
      description: "HiFi reads in FASTA (one per movie)"
    }
    extra_params: {
      description: "Additional hifiasm parameters (trio: -c/-d binning thresholds)"
    }
    father_yak: {
      description: "yak k-mer table of the father (trio mode)"
    }
    mother_yak: {
      description: "yak k-mer table of the mother (trio mode)"
    }
    threads: {
      description: "Threads"
    }
    mem_gb: {
      description: "Memory (GB); ~135 GB observed for a 30x human genome, 288 GB leaves headroom"
    }
    ugc_wgw_container_registry: {
      description: "Registry holding the ugc-built images"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    String sample_id
    Array[File] reads_fastas
    String? extra_params
    File? father_yak
    File? mother_yak
    Int threads = 48
    Int mem_gb = 288
    String ugc_wgw_container_registry
    RuntimeAttributes runtime_attributes
  }

  String prefix = "~{sample_id}.asm"
  Int disk_size = ceil(size(reads_fastas, "GB") * 4 + 20)

  command <<<
    set -euo pipefail

    hifiasm --version >&2

    # shellcheck disable=SC2086
    hifiasm \
      -o "~{prefix}" \
      -t ~{threads} \
      ~{default="" extra_params} \
      ~{"-1 " + father_yak} \
      ~{"-2 " + mother_yak} \
      ~{sep=" " reads_fastas}

    # exactly one primary contig graph per haplotype is expected (upstream collected them by glob only)
    for hap in hap1 hap2; do
      n=$(find . -maxdepth 1 -name "~{prefix}.*.${hap}.p_ctg.gfa" | wc -l)
      if [ "$n" -ne 1 ]; then
        echo "expected exactly one ${hap} p_ctg.gfa for ~{prefix}, found ${n}" >&2
        exit 1
      fi
    done
  >>>

  output {
    Array[File] assembly_hap_gfas = glob("~{prefix}.*.hap[12].p_ctg.gfa")
    Array[File] assembly_noseq_gfas = flatten([
      glob("~{prefix}.*.hap[12].p_ctg.noseq.gfa"),
      glob("~{prefix}.*.[pr]_utg.noseq.gfa")
    ])
    Array[File] assembly_lowQ_beds = flatten([
      glob("~{prefix}.*.hap[12].p_ctg.lowQ.bed"),
      glob("~{prefix}.*.[pr]_utg.lowQ.bed")
    ])
  }

  runtime {
    docker: "~{ugc_wgw_container_registry}/hifiasm@sha256:0e03f1dc2ab6ede5200bbc0e8e9c3ea10e693ce87cad9346e3a1aa59cec272cf"  # 0.25.0 (containers/hifiasm)
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
