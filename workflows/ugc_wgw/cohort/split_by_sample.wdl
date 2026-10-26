version 1.0

# ugc-pacbio-wgw — split a multi-sample VCF into per-sample VCFs in one pass (bcftools +split)
# origin: none
# see docs/DESIGN.md §7.2

import "../../../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/structs.wdl"

task ugc_wgw_bcftools_split {
  meta {
    description: "bcftools +split of a cohort VCF into single-sample VCFs, one process per contig that holds records (threads at a time, pieces on the node-local scratch), then one concatenation per sample in contig order. keep_homref=true reproduces upstream split_vcf_by_sample (--exclude-uncalled: only sites where the sample's genotype is missing are dropped); keep_homref=false keeps only sites where the sample carries an alt allele"
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
      description: "+split processes, then per-sample concatenations, running at once (bcftools +split itself has no threads)"
    }
    mem_gb: {
      description: "Memory (GB); every +split process holds one output stream per sample, so scale with samples × threads (cohort_call: 8 + 0.04 per sample at 8 threads)"
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

    # one output stream per sample in every +split process
    ulimit -Sn 65536 || true

    ln --symbolic --verbose "~{vcf}" .
    ln --symbolic --verbose "~{vcf_index}" .

    # samples file: sample_id <TAB> - (keep name) <TAB> output base name (without .vcf.gz)
    paste "~{write_lines(sample_ids)}" "~{write_lines(split_vcf_names)}" \
    | awk -F'\t' -v OFS='\t' '{ sub(/\.vcf\.gz$/, "", $2); print $1, "-", $2 }' \
    > samples.tsv

    # bcftools +split (1.23) parses every record and has no --threads: one pass over 800 samples took four
    # hours (about a microsecond per genotype, then one single-threaded bgzip stream per sample; 2026-10-26).
    # It runs here once per contig that holds records, `threads` at a time, each writing its per-sample
    # pieces at compression level 1 to the node-local scratch when the backend bound one (TMPDIR); then every
    # sample's pieces are concatenated in the header's contig order and recompressed, `threads` samples at a
    # time. A contig's region is a BED line, not a -r string: GRCh38's HLA contigs carry colons in their names.
    PIECES=${TMPDIR:-$PWD}/split_pieces
    VCF=~{basename(vcf)}
    KEEP_HOMREF=~{keep_homref}
    export KEEP_HOMREF PIECES VCF
    mkdir -p "$PIECES" regions
    bcftools index --stats "$VCF" \
    | awk -F'\t' -v OFS='\t' '$3 > 0 { print $1, 0, $2 > ("regions/" NR ".bed"); print NR, $1 }' > contigs.txt
    echo "$(wc -l < contigs.txt) contigs with records, $(wc -l < samples.tsv) samples, ~{threads} at a time" >&2
    cat > split_region.sh <<'SH'
#!/usr/bin/env bash
# <index> <contig>: one contig's per-sample pieces under $PIECES/<index>/
set -euo pipefail
mkdir -p "$PIECES/$1"
# keep_homref: drop only uncalled genotypes (upstream parity); else keep the sites where the sample carries an alt
if [ "$KEEP_HOMREF" = "true" ]; then set -- "$1" --exclude 'GT="mis"'; else set -- "$1" --include 'GT="alt"'; fi
bcftools +split --regions-file "regions/$1.bed" --regions-overlap 0 --samples-file samples.tsv \
  "$2" "$3" --output-type z1 --output "$PIECES/$1" "$VCF"
SH
    cat > concat_sample.sh <<'SH'
#!/usr/bin/env bash
# <output base name>: the sample's pieces in contig order, recompressed and indexed in the task directory
set -euo pipefail
: > "$PIECES/$1.list"
while read -r i _; do
  if [ -e "$PIECES/$i/$1.vcf.gz" ]; then echo "$PIECES/$i/$1.vcf.gz" >> "$PIECES/$1.list"; fi
done < contigs.txt
bcftools concat --file-list "$PIECES/$1.list" --output-type z --output "$1.vcf.gz" 2> "$PIECES/$1.concat.log"
bcftools index --tbi "$1.vcf.gz"
SH
    chmod +x split_region.sh concat_sample.sh
    xargs -P ~{threads} -L 1 ./split_region.sh < contigs.txt
    cut -f3 samples.tsv | xargs -P ~{threads} -I{} ./concat_sample.sh {}
    rm -rf "$PIECES"
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
