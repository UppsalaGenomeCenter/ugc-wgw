version 1.0
# ugc-pacbio-wgw — dipcall v0.3 (the official kit, as HPRC packages it): diploid small-variant calls and the confidently covered region from two haplotype assemblies
# origin: none (the kit's own run-dipcall + make, unchanged)
# see docs/DESIGN.md §5, §13; docs/ENTRYPOINTS.md §8.1
import "../../../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/structs.wdl"

task ugc_wgw_dipcall {
  meta {
    description: "run-dipcall [-x PAR] <sample_id> <ref> <hap1> <hap2> | make -j2: minimap2 2.20 (the kit's) aligns each haplotype, dipcall-aux pairs the calls. Male mode (-x, when xy) assumes chrY contigs in hap1 and non-PAR chrX contigs in hap2: the caller passes trio-binned haplotypes as they are and regrouped copies otherwise. The sample column is the kit's `syndip`; vcf_finish renames it. dipcall calls chr1-22, chrX and, with -x, chrY."
    outputs: {
      vcf: {
        description: "<sample_id>.dip.vcf.gz: the diploid calls, sample `syndip`, dipcall's FILTER values (GAP1/2, HET1/2, DIPX/Y)"
      },
      bed: {
        description: "<sample_id>.dip.bed: regions covered by exactly one contig of each haplotype (the confident region)"
      },
      hap_beds: {
        description: "<sample_id>.hap1.bed and .hap2.bed: per-haplotype covered regions"
      }
    }
  }

  parameter_meta {
    sample_id: {
      description: "Sample ID; dipcall's output prefix"
    }
    ref_fasta: {
      description: "Assembly reference FASTA (the map's asm_fasta)"
    }
    ref_index: {
      description: "Its .fai"
    }
    hap1_fasta: {
      description: "Haplotype 1 contigs (gzipped FASTA); paternal, or the chrY-bearing regrouped copy, when xy"
    }
    hap2_fasta: {
      description: "Haplotype 2 contigs (gzipped FASTA)"
    }
    xy: {
      description: "Run dipcall's male mode (-x par_bed)"
    }
    par_bed: {
      description: "Pseudoautosomal regions of chrX (dipcall's hs38.PAR.bed)"
    }
    threads: {
      description: "Threads per minimap2 (-t); make -j2 aligns both haplotypes at once, so the task takes twice this"
    }
    mem_gb: {
      description: "Memory (GB) for two concurrent minimap2 asm5 runs against GRCh38"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    String sample_id
    File ref_fasta
    File ref_index
    File hap1_fasta
    File hap2_fasta
    Boolean xy
    File par_bed
    Int threads = 8
    Int mem_gb = 48
    RuntimeAttributes runtime_attributes
  }

  Int cpu = threads * 2
  Int disk_size = ceil((size(hap1_fasta, "GB") + size(hap2_fasta, "GB")) * 8 + size(ref_fasta, "GB") * 2 + 20)

  command <<<
    set -euo pipefail
    export LC_ALL=C
    KIT=/opt/dipcall/dipcall.kit
    "$KIT/minimap2" --version >&2
    "$KIT/samtools" --version | head -1 >&2
    "$KIT/k8" "$KIT/paftools.js" version >&2

    # the kit indexes ref.fa next to itself and wants plain FASTA haplotypes
    ln -s "~{ref_fasta}" ref.fa
    ln -s "~{ref_index}" ref.fa.fai
    gzip -dc "~{hap1_fasta}" > hap1.fa
    gzip -dc "~{hap2_fasta}" > hap2.fa

    "$KIT/run-dipcall" \
      -t ~{threads} \
      ~{if xy then "-x '" + par_bed + "'" else ""} \
      "~{sample_id}" ref.fa hap1.fa hap2.fa \
    > "~{sample_id}.mak"
    make -j2 -f "~{sample_id}.mak"
    rm -f hap1.fa hap2.fa
    ls -la "~{sample_id}".* >&2
  >>>

  output {
    File vcf = "~{sample_id}.dip.vcf.gz"
    File bed = "~{sample_id}.dip.bed"
    Array[File] hap_beds = ["~{sample_id}.hap1.bed", "~{sample_id}.hap2.bed"]
  }

  runtime {
    docker: "humanpangenomics/hpp_dipcall_v0.3@sha256:2c9293158b5dffcc66a39a687ef3b3b3492886a027e195d56b2ce6434bc0f2b6"  # dipcall v0.3 kit, HPRC image (latest, 2021-08)
    cpu: cpu
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
