version 1.0

# ugc-pacbio-wgw — split a multi-sample VCF into per-sample VCFs in one pass (bcftools +split)
# origin: none
# see docs/DESIGN.md §7.2

import "../../../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/structs.wdl"

task ugc_wgw_bcftools_split {
  meta {
    description: "One-pass bcftools +split of a cohort VCF into single-sample VCFs. keep_homref=true reproduces upstream split_vcf_by_sample (--exclude-uncalled: only sites where the sample's genotype is missing are dropped); keep_homref=false keeps only sites where the sample carries an alt allele"
    outputs: {
      split_vcfs: {
        description: "Single-sample VCFs, in sample_ids order"
      },
      split_vcf_indices: {
        description: "Indices for the single-sample VCFs, same order"
      }
    }
  }

  parameter_meta {
    vcf: {
      description: "Multi-sample VCF (bgzipped, indexed)"
    }
    vcf_index: {
      description: "VCF index"
    }
    sample_ids: {
      description: "Samples to extract; every ID must exist in the VCF header"
    }
    split_vcf_names: {
      description: "Output file names, same order as sample_ids (must end in .vcf.gz)"
    }
    split_vcf_index_names: {
      description: "Output index names, same order"
    }
    keep_homref: {
      description: "true: drop only uncalled genotypes (upstream parity); false: keep only GT=\"alt\" sites"
    }
    threads: {
      description: "Threads"
    }
    mem_gb: {
      description: "Memory (GB); one bgzip writer per output, scale with sample count"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    File vcf
    File vcf_index
    Array[String] sample_ids
    Array[String] split_vcf_names
    Array[String] split_vcf_index_names
    Boolean keep_homref = true
    Int threads = 4
    Int mem_gb = 8
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(vcf, "GB") * 3 + 20)

  command <<<
    set -euo pipefail

    bcftools --version >&2

    # one output stream per sample
    ulimit -Sn 65536 || true

    ln --symbolic --verbose "~{vcf}" .
    ln --symbolic --verbose "~{vcf_index}" .

    # samples file: sample_id <TAB> - (keep name) <TAB> output base name (without .vcf.gz)
    paste "~{write_lines(sample_ids)}" "~{write_lines(split_vcf_names)}" \
    | awk -F'\t' -v OFS='\t' '{ sub(/\.vcf\.gz$/, "", $2); print $1, "-", $2 }' \
    > samples.tsv

    if [ "~{keep_homref}" = "true" ]; then
      filter=(--exclude 'GT="mis"')
    else
      filter=(--include 'GT="alt"')
    fi

    # bcftools +split (1.23) has no --threads option; `threads` only sizes the cpu reservation
    bcftools +split \
      --samples-file samples.tsv \
      "${filter[@]}" \
      --output-type z \
      --write-index=tbi \
      --output . \
      "~{basename(vcf)}"
  >>>

  output {
    Array[File] split_vcfs = split_vcf_names
    Array[File] split_vcf_indices = split_vcf_index_names
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
