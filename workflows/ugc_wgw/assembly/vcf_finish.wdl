version 1.0
# ugc-pacbio-wgw — finish an assembly caller's VCF: rename the sample, unphase the genotypes unless trio-binned, bgzip, tabix, bcftools stats
# origin: none
# see docs/DESIGN.md §5, §13; docs/ENTRYPOINTS.md §8.1
import "../../../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/structs.wdl"

task ugc_wgw_vcf_finish {
  meta {
    description: "The stage's primary VCF from a caller's raw one: optional sample rename (dipcall writes `syndip`), genotypes unphased with bcftools +setGT when the haplotypes are not trio-binned (hifiasm's hap1/hap2 phase is local, so column order means nothing beyond it), then bgzip, tabix and bcftools stats against the reference. With reliable_bed, every record outside it gets FILTER ASM_UNRELIABLE added (bcftools filter --mask-file ^bed --soft-filter, mode +): calls are annotated, never dropped. The callers' FILTER values are kept."
    outputs: {
      vcf: {
        description: "<sample_id>.<caller>.<ref_name>.vcf.gz"
      },
      vcf_index: {
        description: "Tabix index"
      },
      stats: {
        description: "bcftools stats of the finished VCF"
      }
    }
  }

  parameter_meta {
    vcf: {
      description: "The caller's VCF (bgzipped or plain)"
    }
    sample_id: {
      description: "Sample ID: the sample column and the output prefix"
    }
    caller: {
      description: "dipcall or pav: part of the output name"
    }
    ref_name: {
      description: "Reference name: part of the output name"
    }
    ref_fasta: {
      description: "Reference FASTA for bcftools stats"
    }
    ref_index: {
      description: "Its .fai"
    }
    unphase: {
      description: "Turn every genotype's | into / (true unless the haplotypes are trio-binned)"
    }
    rename_sample_from: {
      description: "When set, the single sample column with this name becomes sample_id"
    }
    reliable_bed: {
      description: "When given, records outside these regions get FILTER ASM_UNRELIABLE (added to any existing filter)"
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
    String sample_id
    String caller
    String ref_name
    File ref_fasta
    File ref_index
    Boolean unphase
    String? rename_sample_from
    File? reliable_bed
    Int threads = 4
    Int mem_gb = 8
    RuntimeAttributes runtime_attributes
  }

  String out_name = "~{sample_id}.~{caller}.~{ref_name}.vcf.gz"
  Int disk_size = ceil(size(vcf, "GB") * 6 + size(ref_fasta, "GB") + 10)

  command <<<
    set -euo pipefail
    bcftools --version | head -1 >&2
    ln -s "~{ref_fasta}" ref.fa
    ln -s "~{ref_index}" ref.fa.fai

    bcftools view -Oz -o step0.vcf.gz "~{vcf}"
    if [ -n "~{default="" rename_sample_from}" ]; then
      n=$(bcftools query --list-samples step0.vcf.gz | wc -l)
      [ "$n" -eq 1 ] || { echo "expected one sample column to rename, found $n" >&2; exit 1; }
      printf '%s\n' "~{sample_id}" > samples.txt
      bcftools reheader --samples samples.txt --output step1.vcf.gz step0.vcf.gz
    else
      mv step0.vcf.gz step1.vcf.gz
    fi
    if ~{true="true" false="false" unphase}; then
      bcftools +setGT --threads ~{threads} -Oz -o step2.vcf.gz step1.vcf.gz -- --target-gt a --new-gt u
    else
      mv step1.vcf.gz step2.vcf.gz
    fi
    if [ -n "~{default="" reliable_bed}" ]; then
      sort -k1,1 -k2,2n "~{reliable_bed}" > reliable.bed
      printf '##FILTER=<ID=ASM_UNRELIABLE,Description="Outside the reliable region: not covered by exactly one contig of each haplotype (dipcall dip.bed) or flagged by HMM-Flagger (Err, Dup or Col), lifted to the reference">\n' > filter.hdr
      bcftools annotate --header-lines filter.hdr -Ou step2.vcf.gz \
        | bcftools filter --mask-file ^reliable.bed --mask-overlap 1 --soft-filter ASM_UNRELIABLE --mode + --threads ~{threads} -Oz -o "~{out_name}"
    else
      bcftools view --threads ~{threads} -Oz -o "~{out_name}" step2.vcf.gz
    fi
    tabix --preset vcf "~{out_name}"
    bcftools stats --threads ~{threads} --fasta-ref ref.fa --samples - "~{out_name}" > "~{out_name}.stats.txt"
    rm -f step0.vcf.gz step1.vcf.gz step2.vcf.gz
  >>>

  output {
    File vcf_out = out_name
    File vcf_index = "~{out_name}.tbi"
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
