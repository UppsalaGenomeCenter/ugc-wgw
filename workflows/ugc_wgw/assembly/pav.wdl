version 1.0
# ugc-pacbio-wgw — PAV 3 (Phased Assembly Variant caller): small variants, SVs and inversions from two haplotype assemblies against the assembly reference
# origin: none (pav3 batch, the upstream image; the Snakemake profile is the package's default with align_map's threads set to align_threads and its memory lowered to 16 GiB)
# see docs/DESIGN.md §5, §13; docs/ENTRYPOINTS.md §8.1
import "../../../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/structs.wdl"

task ugc_wgw_pav {
  meta {
    description: "pav3 batch --cores threads -p profile in a directory holding pav.json (reference only: PAV's defaults apply, minimap2 -x asm20) and assemblies.tsv (name, hap_h1/hap_h2, or hap_pat/hap_mat when trio). The profile is PAV's default (data/workflow/profiles/pav_default) with align_map's threads set to align_threads. Writes <sample_id>.vcf.gz (.csi; sample column = sample_id, haplotypes phased in column order, PAV's own FILTER values) and the callable reference regions per haplotype, converted from Parquet to BED with pav3 maketable. The image has no tabix or bcftools: vcf_finish indexes and unphases."
    outputs: {
      vcf: {
        description: "<sample_id>.pav.raw.vcf.gz: PAV's VCF as written (csi-indexed)"
      },
      vcf_index: {
        description: "Its .csi"
      },
      callable_beds: {
        description: "callable_ref_<hap>.bed per haplotype: reference regions PAV could call, from its Parquet tables"
      }
    }
  }

  parameter_meta {
    sample_id: {
      description: "Sample ID: PAV's assembly name and the VCF's sample column"
    }
    ref_fasta: {
      description: "Assembly reference FASTA (the map's asm_fasta; PAV indexes it itself)"
    }
    hap1_fasta: {
      description: "Haplotype 1 contigs (gzipped FASTA)"
    }
    hap2_fasta: {
      description: "Haplotype 2 contigs (gzipped FASTA)"
    }
    trio: {
      description: "Trio-binned haplotypes: columns hap_pat/hap_mat instead of hap_h1/hap_h2"
    }
    threads: {
      description: "Cores Snakemake may use at once (pav3 batch --cores)"
    }
    align_threads: {
      description: "Threads of each align_map job (the profile's set-threads; PAV's default is 4). Both haplotypes align at once when threads allows"
    }
    mem_gb: {
      description: "Memory (GB): about 16 GB per concurrent align_map against GRCh38 (the asm20 index), plus headroom"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    String sample_id
    File ref_fasta
    File hap1_fasta
    File hap2_fasta
    Boolean trio
    Int threads = 8
    Int align_threads = 4
    Int mem_gb = 40
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil((size(hap1_fasta, "GB") + size(hap2_fasta, "GB")) * 10 + size(ref_fasta, "GB") * 3 + 20)
  String h1 = if trio then "hap_pat" else "hap_h1"
  String h2 = if trio then "hap_mat" else "hap_h2"

  command <<<
    set -euo pipefail
    pav3 --version >&2
    export HOME="$PWD"

    mkdir -p profile
    cat > profile/config.yaml <<'YAML'
---
latency-wait: 60
keep-going: False
rerun-incomplete: True
printshellcmds: False
rerun-triggers: "mtime"
default-resources:
  - mem_mb=2048
  - runtime=120
set-threads:
  - data_ref_info_table=4
  - align_map=~{align_threads}
  - align_tables=4
  - align_depth=4
  - call_intra_snv_insdel=4
  - call_intra_inv_flag=4
  - call_intra_inv=4
  - call_inter=4
  - call_integrate_sources=4
set-resources:
  - data_ref_info_table:runtime=240
  - align_map:mem_mb=16384
  - align_map:runtime=720
  - align_depth:mem_mb=8192
  - align_tables:mem_mb=8192
  - call_intra_inv:mem_mb=8192
  - call_intra_inv:runtime=240
  - call_intra_snv_insdel:mem_mb=8192
  - call_intra_snv_insdel:runtime=240
  - call_inter:mem_mb=8192
  - call_inter:runtime=240
  - call_integrate_sources:mem_mb=16384
YAML
    printf '{"reference": "%s"}\n' "~{ref_fasta}" > pav.json
    printf 'name\t~{h1}\t~{h2}\n~{sample_id}\t%s\t%s\n' "~{hap1_fasta}" "~{hap2_fasta}" > assemblies.tsv
    cat pav.json assemblies.tsv >&2

    pav3 batch --cores ~{threads} -p profile

    mv "~{sample_id}.vcf.gz" "~{sample_id}.pav.raw.vcf.gz"
    mv "~{sample_id}.vcf.gz.csi" "~{sample_id}.pav.raw.vcf.gz.csi"
    for t in results/"~{sample_id}"/call_hap/callable_ref_*.parquet; do
      hap=$(basename "$t" .parquet); hap=${hap#callable_ref_}
      pav3 maketable -i "$t" -o "~{sample_id}.pav.callable_ref_${hap}.bed"
    done
    ls -la "~{sample_id}".pav.* >&2
  >>>

  output {
    File vcf = "~{sample_id}.pav.raw.vcf.gz"
    File vcf_index = "~{sample_id}.pav.raw.vcf.gz.csi"
    Array[File] callable_beds = glob("~{sample_id}.pav.callable_ref_*.bed")
  }

  runtime {
    docker: "ghcr.io/becklaboratory/pav3@sha256:729964f51e6c2d6f84c42f56065d218bfd1f275026ad50b821898dc4874925c9"  # 3.0.0.dev26
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
