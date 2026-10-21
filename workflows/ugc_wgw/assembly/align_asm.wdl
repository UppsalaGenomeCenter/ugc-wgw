version 1.0
# ugc-pacbio-wgw — align one haplotype's contigs to the assembly reference with minimap2 2.31: a sorted BAM for viewing and a PAF (M/I/D CIGAR, cs tag) for coverage, karyotype and liftover
# origin: the ported align_hifiasm task (HiFi-human-assembly-WDL v1.0.2, 025b3910), rewritten for our minimap2 image; paftools call retired
# see docs/DESIGN.md §5, §13; docs/ENTRYPOINTS.md §8.1
import "../../../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/structs.wdl"

task ugc_wgw_minimap2_align_asm {
  meta {
    description: "minimap2 -x asm5 --cs of one haplotype's contigs to the reference. One alignment run writes a SAM; samtools sort makes the BAM people view, paftools.js sam2paf the PAF the stage reads. No --eqx and no -L: sam2paf copies the SAM CIGAR into the cg tag as it is, and paftools liftover only understands M/I/D operations (=/X make it fail with an inconsistent CIGAR); -L would put long CIGARs into a CG tag that sam2paf does not read. The cs tag carries the differences."
    outputs: {
      bam: {
        description: "Aligned contigs, <sample_id>.<haplotype>.asm.<ref_name>.bam"
      },
      bam_index: {
        description: "BAM index"
      },
      paf: {
        description: "The same alignments as PAF with cs tags, gzipped"
      }
    }
  }

  parameter_meta {
    sample_id: {
      description: "Sample ID (read group SM)"
    }
    haplotype: {
      description: "Haplotype label (hap1/hap2, or pat/mat when trio-binned); part of the output names"
    }
    query_fasta: {
      description: "Contigs (bgzipped FASTA)"
    }
    ref_fasta: {
      description: "Assembly reference FASTA (the map's asm_fasta)"
    }
    ref_name: {
      description: "Reference name"
    }
    threads: {
      description: "Threads (minimap2 gets threads - 3, samtools sort 3)"
    }
    mem_gb: {
      description: "Memory (GB); the asm5 index of GRCh38 takes about 12 GB, the sort 6 GB"
    }
    ugc_wgw_container_registry: {
      description: "Registry holding the ugc-built images (minimap2)"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    String sample_id
    String haplotype
    File query_fasta
    File ref_fasta
    String ref_name
    Int threads = 16
    Int mem_gb = 32
    String ugc_wgw_container_registry
    RuntimeAttributes runtime_attributes
  }

  String out_prefix = "~{sample_id}.~{haplotype}.asm.~{ref_name}"
  Int minimap2_threads = if threads > 3 then threads - 3 else 1
  Int disk_size = ceil((size(query_fasta, "GB") + size(ref_fasta, "GB")) * 6 + 20)

  command <<<
    set -euo pipefail

    minimap2 --version >&2
    samtools --version | head -1 >&2
    paftools.js version >&2

    mkdir -p TMP
    minimap2 \
      -t ~{minimap2_threads} \
      -a \
      -x asm5 \
      --cs \
      --secondary=no \
      -R "@RG\tID:~{sample_id}_~{haplotype}\tSM:~{sample_id}" \
      "~{ref_fasta}" \
      "~{query_fasta}" \
    > aln.sam
    samtools sort \
      -@ 3 \
      -T ./TMP/sort \
      -m 2G \
      -O BAM \
      -o "~{out_prefix}.bam" \
      aln.sam
    samtools index "~{out_prefix}.bam"
    paftools.js sam2paf aln.sam | gzip -1 > "~{out_prefix}.paf.gz"
    rm -f aln.sam
  >>>

  output {
    File bam = "~{out_prefix}.bam"
    File bam_index = "~{out_prefix}.bam.bai"
    File paf = "~{out_prefix}.paf.gz"
  }

  runtime {
    docker: "~{ugc_wgw_container_registry}/minimap2@sha256:00dd2dfbb7f395672eb0e71317899f913cd64af80c02864655dc9d848202604a"  # 2.31 + samtools 1.24 + k8 (containers/minimap2)
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
