version 1.0

# ugc-pacbio-wgw — GLnexus over one scatter region (override of upstream task: --bed required, --list input, ulimit, scratch dir, sized per shard)
# origin: vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/tasks/glnexus.wdl @ v4.0.0
# see docs/DESIGN.md §7.2, DERIVED_FILES.md

import "../../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/structs.wdl"

task glnexus {
  meta {
    description: "Merge gVCFs from multiple samples and joint call genotypes with GLnexus, restricted to the regions of one shard (ugc-wgw override of the upstream task; config unchanged)"
    outputs: {
      vcf: {
        description: "Joint-called multi-sample VCF"
      },
      vcf_index: {
        description: "Index for joint-called multi-sample VCF"
      }
    }
  }

  parameter_meta {
    cohort_id: {
      description: "Cohort ID"
    }
    gvcfs: {
      description: "GVCFs"
    }
    gvcf_indices: {
      description: "GVCF indices"
    }
    ref_name: {
      description: "Reference name"
    }
    regions_bed: {
      description: "Regions BED of the shard (glnexus_cli --bed)"
    }
    shard: {
      description: "Shard name, used in the output file name"
    }
    mem_gb: {
      description: "Memory (GB); passed as --mem-gbytes and requested from the scheduler"
    }
    threads: {
      description: "Threads (glnexus_cli --threads)"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    String cohort_id
    Array[File] gvcfs
    Array[File] gvcf_indices
    String ref_name
    File regions_bed
    String shard
    Int threads = 32
    Int mem_gb = 60
    RuntimeAttributes runtime_attributes
  }

  String out_prefix = "~{cohort_id}.~{ref_name}.~{shard}"
  Int disk_size = ceil(size(gvcfs, "GB") * 2 + 100)

  command <<<
    set -euo pipefail

    # we use a custom config file based on DeepVariant_unfiltered, but with required_dp: 1
    cat << EOF > config.yml
    unifier_config:
      min_AQ1: 0
      min_AQ2: 0
      min_GQ: 0
      monoallelic_sites_for_lost_alleles: true
      max_alleles_per_site: 32
    genotyper_config:
      required_dp: 1
      revise_genotypes: false
      allow_partial_data: true
      more_PL: true
      trim_uncalled_alleles: true
      liftover_fields:
        - orig_names: [MIN_DP, DP]
          name: DP
          description: '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">'
          type: int
          combi_method: min
          number: basic
          count: 1
          ignore_non_variants: true
        - orig_names: [AD]
          name: AD
          description: '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">'
          type: int
          number: alleles
          combi_method: min
          default_type: zero
          count: 0
        - orig_names: [GQ]
          name: GQ
          description: '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">'
          type: int
          number: basic
          combi_method: min
          count: 1
          ignore_non_variants: true
        - orig_names: [PL]
          name: PL
          description: '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype Likelihoods">'
          type: int
          number: genotype
          combi_method: missing
          count: 0
          ignore_non_variants: true
    EOF

    # glneux_cli has no version option
    glnexus_cli --help 2>&1 | grep -Eo 'glnexus_cli release v[0-9a-f.-]+'
    bcftools --version

    # RocksDB opens many files in large cohorts (GLnexus wiki, Performance)
    ulimit -Sn 65536 || true
    # scratch on node-local disk when the scheduler provides one, else the task dir
    scratch="${TMPDIR:-$PWD}/~{out_prefix}.GLnexus.DB"

    glnexus_cli \
      --threads ~{threads} \
      --mem-gbytes ~{mem_gb} \
      --dir "${scratch}" \
      --config ./config.yml \
      --bed "~{regions_bed}" \
      --list "~{write_lines(gvcfs)}" \
    > "~{out_prefix}.small_variants.bcf"
    bcftools view \
      ~{if threads > 1
        then "--threads '" + (threads - 1) + "'"
        else ""
      } \
      --output-type z \
      --output-file "~{out_prefix}.small_variants.vcf.gz" \
      "~{out_prefix}.small_variants.bcf"
    bcftools index \
      ~{if threads > 1
        then "--threads '" + (threads - 1) + "'"
        else ""
      } \
      --tbi \
      "~{out_prefix}.small_variants.vcf.gz"
    # cleanup intermediate files
    rm --recursive --force "${scratch}" "~{out_prefix}.small_variants.bcf"
  >>>

  output {
    File vcf = "~{out_prefix}.small_variants.vcf.gz"
    File vcf_index = "~{out_prefix}.small_variants.vcf.gz.tbi"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/glnexus@sha256:ce6fecf59dddc6089a8100b31c29c1e6ed50a0cf123da9f2bc589ee4b0c69c8e"  # 1.4.3
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
