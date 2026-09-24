version 1.0

# ugc-pacbio-wgw — GLnexus over a cohort: slice gVCFs per sample, GLnexus per shard, concat
# origin: none
# see docs/DESIGN.md §7.2

import "../../../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/structs.wdl"
import "../../overrides/glnexus.wdl" as Glnexus
import "concat.wdl" as Concat
import "gvcf_slice.wdl" as GvcfSlice

workflow ugc_wgw_glnexus_scatter {
  meta {
    description: "Joint-call small variants across a cohort with GLnexus, scattered by shard: one slicing task per sample, one GLnexus task per shard, one concat"
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
    sample_ids: {
      description: "Sample IDs, same order as gvcfs"
    }
    gvcfs: {
      description: "Per-sample gVCFs"
    }
    gvcf_indices: {
      description: "Per-sample gVCF indices"
    }
    ref_name: {
      description: "Reference name"
    }
    shard_names: {
      description: "Shard names (ugc_wgw_regions_prepare)"
    }
    shard_beds: {
      description: "Shard BEDs, same order"
    }
    out_prefix: {
      description: "Output prefix of the concatenated VCF"
    }
    glnexus_mem_gb: {
      description: "Memory (GB) per GLnexus shard"
    }
    glnexus_threads: {
      description: "Threads per GLnexus shard"
    }
    default_runtime_attributes: {
      description: "Default runtime attribute structure"
    }
  }

  input {
    String cohort_id
    Array[String] sample_ids
    Array[File] gvcfs
    Array[File] gvcf_indices
    String ref_name
    Array[String] shard_names
    Array[File] shard_beds
    String out_prefix
    Int glnexus_mem_gb
    Int glnexus_threads = 32
    RuntimeAttributes default_runtime_attributes
  }

  scatter (sample_index in range(length(gvcfs))) {
    call GvcfSlice.ugc_wgw_bcftools_slice { input:
      gvcf = gvcfs[sample_index],
      gvcf_index = gvcf_indices[sample_index],
      shard_names = shard_names,
      shard_beds = shard_beds,
      out_prefix = sample_ids[sample_index],
      runtime_attributes = default_runtime_attributes
    }
  }

  # [sample][shard] -> [shard][sample]
  Array[Array[File]] slices_by_shard = transpose(ugc_wgw_bcftools_slice.slices)
  Array[Array[File]] slice_indices_by_shard = transpose(ugc_wgw_bcftools_slice.slice_indices)

  scatter (shard_index in range(length(shard_names))) {
    call Glnexus.glnexus { input:
      cohort_id = cohort_id,
      gvcfs = slices_by_shard[shard_index],
      gvcf_indices = slice_indices_by_shard[shard_index],
      ref_name = ref_name,
      regions_bed = shard_beds[shard_index],
      shard = shard_names[shard_index],
      mem_gb = glnexus_mem_gb,
      threads = glnexus_threads,
      runtime_attributes = default_runtime_attributes
    }
  }

  call Concat.ugc_wgw_bcftools_concat { input:
    vcfs = glnexus.vcf,
    vcf_indices = glnexus.vcf_index,
    out_prefix = out_prefix,
    runtime_attributes = default_runtime_attributes
  }

  output {
    File vcf = ugc_wgw_bcftools_concat.vcf
    File vcf_index = ugc_wgw_bcftools_concat.vcf_index
  }
}
