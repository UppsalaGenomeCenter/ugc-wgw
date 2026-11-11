version 1.0

# ugc-pacbio-wgw — pblab's preview tools as tasks: pbweave alignment and pbdv small-variant calling (GPU and CPU), with the shims that make their outputs what pbmm2's and upstream's DeepVariant tasks deliver
# origin: none
# see docs/DESIGN.md §17 (pblab 20261021), §18 (2026-11-11); docs/ENTRYPOINTS.md §3, §4

import "../../../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/structs.wdl"

task ugc_wgw_pbweave_align {
  meta {
    description: "Align all of a sample's HiFi read files with pbweave (pblab, preview) into one coordinate-sorted BAM, as pbmm2 would with --preset HIFI --sample --sort --unmapped --strip: the read group's SM is the sample id, unmapped reads are kept, MM/ML and the other read tags are kept, kinetics are stripped. pbweave has no --strip-tags, so when an input is a haplotagged alignment the HP, PS and PC tags are removed in a second pass. The reference index is built in the task (seconds)."
    outputs: {
      aligned_bam: {
        description: "<sample_id>.<ref_name>.hifi_reads.bam"
      },
      aligned_bam_index: {
        description: "Its BAI"
      },
      msg: {
        description: "What the input files hold, in the words of upstream's index_input_bam: already aligned, no base modification tags, kinetics"
      }
    }
  }

  parameter_meta {
    sample_id: {
      description: "Sample ID (the SM of every read group)"
    }
    bams: {
      description: "HiFi reads, unaligned or aligned BAMs"
    }
    ref_fasta: {
      description: "Reference genome FASTA"
    }
    ref_index: {
      description: "Reference genome FASTA index"
    }
    ref_name: {
      description: "Reference name (part of the output name)"
    }
    strip_kinetics: {
      description: "Strip kinetics tags (pbweave --strip)"
    }
    threads: {
      description: "Threads for the index, the alignment and the sort"
    }
    mem_gb: {
      description: "Memory (GB); pbweave sorts in memory up to this less a margin, then spills to TMPDIR"
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
    Array[File] bams
    File ref_fasta
    File ref_index
    String ref_name
    Boolean strip_kinetics = true
    Int threads = 24
    Int mem_gb = 96
    String ugc_wgw_container_registry
    RuntimeAttributes runtime_attributes
  }

  String out_prefix = "~{sample_id}.~{ref_name}.hifi_reads"
  Int disk_size = ceil(size(bams, "GB") * 3 + size(ref_fasta, "GB") * 3 + 20)

  command <<<
    set -euo pipefail

    pbweave --version
    samtools --version | head -1

    # what each input holds, from its first reads (upstream's index_input_bam, without pysam)
    : > messages.txt
    hap_tags=false
    while read -r bam; do
      read -r aligned basemods kinetics haplotags <<< "$(samtools head --headers 0 --records 2000 "${bam}" | awk -F'\t' '
        { if ($3 != "*") a = 1
          for (i = 12; i <= NF; i++) { t = substr($i, 1, 2)
            if (t == "MM" || t == "ML") m = 1
            if (t == "fi" || t == "ri" || t == "fp" || t == "rp" || t == "ip" || t == "pw") k = 1
            if (t == "HP" || t == "PS" || t == "PC") h = 1 } }
        END { print a + 0, m + 0, k + 0, h + 0 }')"
      name=$(basename "${bam}")
      if [ "${aligned}" = 1 ]; then
        echo "Input ${name} is already aligned.  Alignments and haplotype tags will be stripped." >> messages.txt
      fi
      if [ "${basemods}" = 0 ]; then
        echo "Input ${name} does not contain base modification tags.  5mCpG pileups will not be generated." >> messages.txt
      fi
      if [ "${kinetics}" = 1 ] && [ "~{strip_kinetics}" = true ]; then
        echo "Input ${name} contains consensus kinetics tags. Kinetics will be stripped from the output." >> messages.txt
      fi
      if [ "${haplotags}" = 1 ]; then hap_tags=true; fi
    done < "~{write_lines(bams)}"
    cat messages.txt >&2

    # the fofn is recognised by its extension
    cp "~{write_lines(bams)}" reads.fofn

    pbweave index --num-threads ~{threads} --log-level INFO "~{ref_fasta}" reference.wvi

    pbweave align \
      --num-threads ~{threads} \
      --max-memory ~{if mem_gb > 12 then mem_gb - 4 else mem_gb}G \
      --preset HIFI \
      --sample "~{sample_id}" \
      --sort \
      --unmapped \
      ~{true="--strip" false="" strip_kinetics} \
      --min-length 50 \
      --log-level INFO \
      reference.wvi reads.fofn aligned.bam
    rm reference.wvi

    if [ "${hap_tags}" = true ]; then
      # pbmm2 --strip-tags HP,PS,PC: a haplotagged input must not hand its old phase to HiPhase's consumers
      samtools view --threads ~{threads} --remove-tag HP,PS,PC --bam --write-index \
        --output "~{out_prefix}.bam##idx##~{out_prefix}.bam.bai" aligned.bam
      rm aligned.bam aligned.bam.bai
    else
      mv aligned.bam "~{out_prefix}.bam"
      mv aligned.bam.bai "~{out_prefix}.bam.bai"
    fi
    samtools quickcheck "~{out_prefix}.bam"
  >>>

  output {
    File aligned_bam = "~{out_prefix}.bam"
    File aligned_bam_index = "~{out_prefix}.bam.bai"
    Array[String] msg = read_lines("messages.txt")
  }

  runtime {
    docker: "~{ugc_wgw_container_registry}/pblab@sha256:01eb8964524fae161f10568ddfc0551843883e50ba1fd38d19ca20b817c8ba5b"
    cpu: threads
    memory: "~{mem_gb} GiB"
    disk: "~{disk_size} GB"
    disks: "local-disk ~{disk_size} HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    zones: runtime_attributes.zones
    cpuPlatform: runtime_attributes.cpuPlatform
  }
}

task ugc_wgw_pbdv_call_gpu {
  meta {
    description: "Small variants with pbdv (pblab, preview; DeepVariant 1.10.0's PACBIO model in C++ and CUDA) on one NVIDIA GPU, delivered as upstream's DeepVariant tasks deliver theirs: genotypes unphased without PS (pbdv writes its read phasing into GT), the VCF without uncalled records, both files bgzipped and tabix-indexed under upstream's names. The headers are DeepVariant's own: no line names pbdv, the run manifest does."
    outputs: {
      vcf: {
        description: "<sample_id>.<ref_name>.small_variants.vcf.gz"
      },
      vcf_index: {
        description: "Its tabix index"
      },
      gvcf: {
        description: "<sample_id>.<ref_name>.small_variants.g.vcf.gz"
      },
      gvcf_index: {
        description: "Its tabix index"
      }
    }
  }

  parameter_meta {
    sample_id: {
      description: "Sample ID (the VCF's sample column)"
    }
    aligned_bam: {
      description: "Aligned, coordinate-sorted BAM"
    }
    aligned_bam_index: {
      description: "Its BAI"
    }
    ref_fasta: {
      description: "Reference genome FASTA"
    }
    ref_index: {
      description: "Reference genome FASTA index"
    }
    ref_name: {
      description: "Reference name (part of the output names)"
    }
    threads: {
      description: "Threads (make_examples is the CPU-bound part on a GPU too)"
    }
    mem_gb: {
      description: "Memory (GB)"
    }
    gpuCount: {
      description: "GPUs (pbdv uses the first)"
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
    File aligned_bam
    File aligned_bam_index
    File ref_fasta
    File ref_index
    String ref_name
    Int threads = 16
    Int mem_gb = 64
    Int gpuCount = 1
    String ugc_wgw_container_registry
    RuntimeAttributes runtime_attributes
  }

  String out_prefix = "~{sample_id}.~{ref_name}.small_variants"
  Int disk_size = ceil(size(aligned_bam, "GB") + size(ref_fasta, "GB") + 20)

  command <<<
    set -euo pipefail

    pbdv --version
    bcftools --version | head -1
    nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader >&2 || true

    # pbdv opens the index next to the BAM and the .fai next to the FASTA
    ln --symbolic "~{aligned_bam}" reads.bam
    ln --symbolic "~{aligned_bam_index}" reads.bam.bai
    ln --symbolic "~{ref_fasta}" reference.fasta
    ln --symbolic "~{ref_index}" reference.fasta.fai

    pbdv run \
      --device gpu \
      --num-threads ~{threads} \
      --ref reference.fasta \
      --sample-name "~{sample_id}" \
      --output-vcf raw.vcf.gz \
      --output-gvcf raw.g.vcf.gz \
      --log-level INFO \
      reads.bam

    # as upstream's DeepVariant tasks deliver them: unphased, no PS; the VCF without uncalled records
    bcftools annotate --no-version --remove FORMAT/PS --output-type u raw.vcf.gz \
      | bcftools +setGT --no-version --output-type u -- --target-gt a --new-gt u \
      | bcftools view --threads ~{threads} --exclude-uncalled --output-type z --output "~{out_prefix}.vcf.gz"
    bcftools annotate --no-version --remove FORMAT/PS --output-type u raw.g.vcf.gz \
      | bcftools +setGT --no-version --threads ~{threads} --output-type z --output "~{out_prefix}.g.vcf.gz" -- --target-gt a --new-gt u
    bcftools index --threads ~{threads} --tbi "~{out_prefix}.vcf.gz"
    bcftools index --threads ~{threads} --tbi "~{out_prefix}.g.vcf.gz"
    rm raw.vcf.gz raw.vcf.gz.tbi raw.g.vcf.gz raw.g.vcf.gz.tbi
  >>>

  output {
    File vcf = "~{out_prefix}.vcf.gz"
    File vcf_index = "~{out_prefix}.vcf.gz.tbi"
    File gvcf = "~{out_prefix}.g.vcf.gz"
    File gvcf_index = "~{out_prefix}.g.vcf.gz.tbi"
  }

  runtime {
    docker: "~{ugc_wgw_container_registry}/pblab@sha256:01eb8964524fae161f10568ddfc0551843883e50ba1fd38d19ca20b817c8ba5b"
    cpu: threads
    memory: "~{mem_gb} GiB"
    disk: "~{disk_size} GB"
    disks: "local-disk ~{disk_size} HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    gpuCount: gpuCount
    gpuType: runtime_attributes.gpuType
    zones: runtime_attributes.zones
    cpuPlatform: runtime_attributes.cpuPlatform
  }
}

task ugc_wgw_pbdv_call_cpu {
  meta {
    description: "Small variants with pbdv (pblab, preview; DeepVariant 1.10.0's PACBIO model in C++ and CUDA) on CPUs alone (the same networks on the CPU workers; several times slower than the GPU task and not bit-identical to it), delivered as upstream's DeepVariant tasks deliver theirs: genotypes unphased without PS (pbdv writes its read phasing into GT), the VCF without uncalled records, both files bgzipped and tabix-indexed under upstream's names. The headers are DeepVariant's own: no line names pbdv, the run manifest does."
    outputs: {
      vcf: {
        description: "<sample_id>.<ref_name>.small_variants.vcf.gz"
      },
      vcf_index: {
        description: "Its tabix index"
      },
      gvcf: {
        description: "<sample_id>.<ref_name>.small_variants.g.vcf.gz"
      },
      gvcf_index: {
        description: "Its tabix index"
      }
    }
  }

  parameter_meta {
    sample_id: {
      description: "Sample ID (the VCF's sample column)"
    }
    aligned_bam: {
      description: "Aligned, coordinate-sorted BAM"
    }
    aligned_bam_index: {
      description: "Its BAI"
    }
    ref_fasta: {
      description: "Reference genome FASTA"
    }
    ref_index: {
      description: "Reference genome FASTA index"
    }
    ref_name: {
      description: "Reference name (part of the output names)"
    }
    threads: {
      description: "Threads"
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
    String sample_id
    File aligned_bam
    File aligned_bam_index
    File ref_fasta
    File ref_index
    String ref_name
    Int threads = 32
    Int mem_gb = 64
    String ugc_wgw_container_registry
    RuntimeAttributes runtime_attributes
  }

  String out_prefix = "~{sample_id}.~{ref_name}.small_variants"
  Int disk_size = ceil(size(aligned_bam, "GB") + size(ref_fasta, "GB") + 20)

  command <<<
    set -euo pipefail

    pbdv --version
    bcftools --version | head -1

    # pbdv opens the index next to the BAM and the .fai next to the FASTA
    ln --symbolic "~{aligned_bam}" reads.bam
    ln --symbolic "~{aligned_bam_index}" reads.bam.bai
    ln --symbolic "~{ref_fasta}" reference.fasta
    ln --symbolic "~{ref_index}" reference.fasta.fai

    pbdv run \
      --device cpu \
      --num-threads ~{threads} \
      --ref reference.fasta \
      --sample-name "~{sample_id}" \
      --output-vcf raw.vcf.gz \
      --output-gvcf raw.g.vcf.gz \
      --log-level INFO \
      reads.bam

    # as upstream's DeepVariant tasks deliver them: unphased, no PS; the VCF without uncalled records
    bcftools annotate --no-version --remove FORMAT/PS --output-type u raw.vcf.gz \
      | bcftools +setGT --no-version --output-type u -- --target-gt a --new-gt u \
      | bcftools view --threads ~{threads} --exclude-uncalled --output-type z --output "~{out_prefix}.vcf.gz"
    bcftools annotate --no-version --remove FORMAT/PS --output-type u raw.g.vcf.gz \
      | bcftools +setGT --no-version --threads ~{threads} --output-type z --output "~{out_prefix}.g.vcf.gz" -- --target-gt a --new-gt u
    bcftools index --threads ~{threads} --tbi "~{out_prefix}.vcf.gz"
    bcftools index --threads ~{threads} --tbi "~{out_prefix}.g.vcf.gz"
    rm raw.vcf.gz raw.vcf.gz.tbi raw.g.vcf.gz raw.g.vcf.gz.tbi
  >>>

  output {
    File vcf = "~{out_prefix}.vcf.gz"
    File vcf_index = "~{out_prefix}.vcf.gz.tbi"
    File gvcf = "~{out_prefix}.g.vcf.gz"
    File gvcf_index = "~{out_prefix}.g.vcf.gz.tbi"
  }

  runtime {
    docker: "~{ugc_wgw_container_registry}/pblab@sha256:01eb8964524fae161f10568ddfc0551843883e50ba1fd38d19ca20b817c8ba5b"
    cpu: threads
    memory: "~{mem_gb} GiB"
    disk: "~{disk_size} GB"
    disks: "local-disk ~{disk_size} HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    zones: runtime_attributes.zones
    cpuPlatform: runtime_attributes.cpuPlatform
  }
}
