version 1.0
# ugc-pacbio-wgw — Flagger's unreliable regions lifted from assembly to reference coordinates, and dipcall's confident region without them
# origin: none
# see docs/DESIGN.md §5, §13; docs/ENTRYPOINTS.md §8.1
import "../../../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/structs.wdl"

task ugc_wgw_liftover_unreliable {
  meta {
    description: "The non-Hap intervals of Flagger's prediction BED (Err, Dup, Col; column 4) lifted through one haplotype's contig-to-reference PAF with paftools.js liftover (primary alignments of 50 kb or more at MAPQ 5 or more, dipcall's own thresholds; =/X CIGAR operations are rewritten to M first, since liftover reads M/I/D only). An interval spanning several alignments is lifted once per alignment, clipped; intervals of contigs without such an alignment are dropped, as they are outside dip.bed anyway."
    outputs: {
      unreliable_bed: {
        description: "<sample_id>.<haplotype>.asm.unreliable.bed: lifted non-Hap intervals in reference coordinates, sorted and merged"
      }
    }
  }

  parameter_meta {
    sample_id: {
      description: "Sample ID; output prefix"
    }
    haplotype: {
      description: "Haplotype label; part of the output name"
    }
    flagger_bed: {
      description: "Flagger's prediction BED for the diploid assembly (both haplotypes' contigs)"
    }
    paf: {
      description: "This haplotype's contig-to-reference PAF (gzipped), with cg tags of M/I/D operations"
    }
    min_alignment_length: {
      description: "paftools.js liftover -l"
    }
    min_mapq: {
      description: "paftools.js liftover -q"
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
    File flagger_bed
    File paf
    Int min_alignment_length = 50000
    Int min_mapq = 5
    String ugc_wgw_container_registry
    RuntimeAttributes runtime_attributes
  }

  String out_name = "~{sample_id}.~{haplotype}.asm.unreliable.bed"
  Int disk_size = ceil(size(paf, "GB") * 4 + 10)

  command <<<
    set -euo pipefail
    paftools.js version >&2
    # Flagger's BED starts with a `track` line; column 4 is the label
    awk 'BEGIN {OFS="\t"} /^(track|#)/ {next} $4 != "Hap" {print $1, $2, $3}' "~{flagger_bed}" > nonhap.bed
    echo "non-Hap intervals: $(wc -l < nonhap.bed)" >&2
    # liftover reads M/I/D CIGARs only: =/X (minimap2 --eqx) are rewritten to M in the cg tag
    gzip -dc "~{paf}" | awk 'BEGIN {FS = OFS = "\t"} {for (i = 13; i <= NF; i++) if ($i ~ /^cg:Z:/) gsub(/[=X]/, "M", $i); print}' > aln.paf
    paftools.js liftover -l ~{min_alignment_length} -q ~{min_mapq} aln.paf nonhap.bed \
      | cut -f1-3 \
      | sort -k1,1 -k2,2n \
      | awk 'BEGIN {OFS="\t"} NR == 1 {c=$1; s=$2; e=$3; next} $1 == c && $2 <= e {if ($3 > e) e = $3; next} {print c, s, e; c=$1; s=$2; e=$3} END {if (NR > 0) print c, s, e}' \
      > "~{out_name}"
    echo "lifted intervals: $(wc -l < "~{out_name}")" >&2
    rm -f aln.paf
  >>>

  output {
    File unreliable_bed = out_name
  }

  runtime {
    docker: "~{ugc_wgw_container_registry}/minimap2@sha256:00dd2dfbb7f395672eb0e71317899f913cd64af80c02864655dc9d848202604a"  # 2.31 + k8/paftools.js (containers/minimap2)
    cpu: 2
    memory: "8 GiB"
    disk: "~{disk_size} GB"
    disks: "local-disk ~{disk_size} HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
    cpuPlatform: runtime_attributes.cpuPlatform
  }
}

task ugc_wgw_reliable_regions {
  meta {
    description: "reliable.bed = dipcall's dip.bed minus the lifted unreliable regions of both haplotypes (bedtools subtract), with a summary of the bases: dip.bed, unreliable within it, reliable, and the reliable fraction of the reference's primary contigs."
    outputs: {
      reliable_bed: {
        description: "<sample_id>.asm.reliable.bed"
      },
      unreliable_merged_bed: {
        description: "<sample_id>.asm.unreliable.bed: both haplotypes' unreliable regions, merged"
      },
      summary_tsv: {
        description: "<sample_id>.asm.reliable.summary.tsv: metric, value (reference_bp, dip_bed_bp, unreliable_bp, unreliable_in_dip_bed_bp, reliable_bp, reliable_fraction_of_reference, reliable_fraction_of_dip_bed)"
      }
    }
  }

  parameter_meta {
    sample_id: {
      description: "Sample ID; output prefix"
    }
    dip_bed: {
      description: "dipcall's dip.bed"
    }
    unreliable_beds: {
      description: "The per-haplotype lifted unreliable regions"
    }
    ref_index: {
      description: "The assembly reference's .fai (chromosome lengths for the fraction)"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    String sample_id
    File dip_bed
    Array[File] unreliable_beds
    File ref_index
    RuntimeAttributes runtime_attributes
  }

  command <<<
    set -euo pipefail
    bedtools --version >&2
    sort -k1,1 -k2,2n ~{sep=" " unreliable_beds} | bedtools merge -i - > "~{sample_id}.asm.unreliable.bed"
    sort -k1,1 -k2,2n "~{dip_bed}" > dip.sorted.bed
    bedtools subtract -a dip.sorted.bed -b "~{sample_id}.asm.unreliable.bed" > "~{sample_id}.asm.reliable.bed"
    bedtools intersect -a dip.sorted.bed -b "~{sample_id}.asm.unreliable.bed" > unreliable_in_dip.bed
    ref=$(awk '$1 ~ /^chr([0-9]+|X|Y)$/ {s+=$2} END {print s+0}' "~{ref_index}")
    dip=$(awk '{s+=$3-$2} END {print s+0}' dip.sorted.bed)
    unr=$(awk '{s+=$3-$2} END {print s+0}' "~{sample_id}.asm.unreliable.bed")
    unrd=$(awk '{s+=$3-$2} END {print s+0}' unreliable_in_dip.bed)
    rel=$(awk '{s+=$3-$2} END {print s+0}' "~{sample_id}.asm.reliable.bed")
    {
      printf 'reference_bp\t%s\n' "$ref"
      printf 'dip_bed_bp\t%s\n' "$dip"
      printf 'unreliable_bp\t%s\n' "$unr"
      printf 'unreliable_in_dip_bed_bp\t%s\n' "$unrd"
      printf 'reliable_bp\t%s\n' "$rel"
      awk -v r="$ref" -v x="$rel" 'BEGIN {printf "reliable_fraction_of_reference\t%.4f\n", (r > 0 ? x / r : 0)}'
      awk -v d="$dip" -v x="$rel" 'BEGIN {printf "reliable_fraction_of_dip_bed\t%.4f\n", (d > 0 ? x / d : 0)}'
    } > "~{sample_id}.asm.reliable.summary.tsv"
    cat "~{sample_id}.asm.reliable.summary.tsv" >&2
  >>>

  output {
    File reliable_bed = "~{sample_id}.asm.reliable.bed"
    File unreliable_merged_bed = "~{sample_id}.asm.unreliable.bed"
    File summary_tsv = "~{sample_id}.asm.reliable.summary.tsv"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:03cb3c01937eccc907f8ad71c87b258581504572205fe3f31a657e318f3564ae"  # pb_wdl_base:build4
    cpu: 1
    memory: "4 GiB"
    disk: "10 GB"
    disks: "local-disk 10 HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
    cpuPlatform: runtime_attributes.cpuPlatform
  }
}
