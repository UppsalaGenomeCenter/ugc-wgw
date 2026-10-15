version 1.0

# ugc-pacbio-wgw — pbmm2 alignment subworkflow with an optional prebuilt reference index (override of upstream subworkflow: create_pbmm2_index only when no index is given; the tasks are upstream's, imported, not copied)
# origin: vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/workflows/pbmm2/pbmm2.wdl @ v4.0.0 (workflow block only)
# see docs/DESIGN.md §18 (2026-10-15), DERIVED_FILES.md

import "../../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/structs.wdl"
import "../../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/workflows/pbmm2/pbmm2.wdl" as Tasks

workflow pbmm2 {
  meta {
    description: "Align HiFi reads to a reference genome with pbmm2 (ugc-wgw override: a prebuilt reference index skips create_pbmm2_index)"
    category: "Utility - Alignment"
    outputs: {
      aligned_bams: {
        description: "Array of aligned BAMs"
      },
      aligned_bam_indices: {
        description: "Array of aligned BAM indices"
      },
      msg: {
        description: "Array of messages"
      }
    }
  }

  parameter_meta {
    sample_id: {
      description: "Sample ID"
    }
    bam: {
      description: "Unaligned BAM"
    }
    use_alignment_chunking: {
      description: "Whether to split the input BAM into chunks for parallel alignment"
    }
    ref_fasta: {
      description: "Reference FASTA"
    }
    pbmm2_index: {
      description: "Prebuilt pbmm2 index (.mmi) of ref_fasta; when absent the index is built by create_pbmm2_index as upstream does (ugc-wgw)"
    }
    ref_name: {
      description: "Reference name"
    }
    default_runtime_attributes: {
      description: "Default runtime attribute structure"
    }
  }

  input {
    String sample_id
    File bam
    Boolean use_alignment_chunking = true
    File ref_fasta
    File? pbmm2_index
    String ref_name
    RuntimeAttributes default_runtime_attributes
  }

  # ugc-wgw: upstream builds the 5.8 GB index once per input BAM; a prebuilt one from the reference map skips it
  if (!defined(pbmm2_index)) {
    call Tasks.create_pbmm2_index { input:
      ref_fasta = ref_fasta,
      runtime_attributes = default_runtime_attributes
    }
  }

  call Tasks.index_input_bam { input:
    bam = bam,
    use_alignment_chunking = use_alignment_chunking,
    runtime_attributes = default_runtime_attributes
  }

  scatter (chunk_index in if (index_input_bam.num_chunks > 0)
    then range(index_input_bam.num_chunks)
    else [
      -1
    ]
  ) {
    call Tasks.pbmm2_align_wgs { input:
      sample_id = sample_id,
      bam = bam,
      pbindex = index_input_bam.pbindex,
      pbmm2_index = select_first([
        pbmm2_index,
        create_pbmm2_index.index
      ]),
      ref_name = ref_name,
      chunk_index = chunk_index + 1,  # wdl is 0-indexed, pbmm2 chunking is 1-indexed
      num_chunks = index_input_bam.num_chunks,
      runtime_attributes = default_runtime_attributes
    }
  }

  output {
    Array[File] aligned_bams = pbmm2_align_wgs.aligned_bam
    Array[File] aligned_bam_indices = pbmm2_align_wgs.aligned_bam_index
    Array[String] msg = index_input_bam.msg
  }
}
