version 1.0
# ugc-pacbio-wgw — per-contig QC of one sample's haplotype assemblies: composition, gaps and telomeres (seqtk) joined with HMM-Flagger's labels
# origin: none
# see docs/DESIGN.md §5, §13; docs/ENTRYPOINTS.md §8
import "../../../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/structs.wdl"

task ugc_wgw_assembly_contig_qc {
  meta {
    description: "One row per contig of both haplotypes: length, GC %, N bases and N runs of 10 bp or more (seqtk comp, seqtk gap), the telomeric bases at the 5' and 3' ends (seqtk telo with the assembly's telomere motif; none when the motif is empty), and the bases HMM-Flagger labelled Err, Dup, Hap and Col (summed from flagger_bed). Runs in the Flagger image, which carries seqtk 1.4; minutes per haplotype."
    outputs: {
      contig_qc: {
        description: "<sample_id>.asm.contig_qc.tsv: haplotype, contig, length, gc_pct, n_bp, gaps, telo_5p_bp, telo_3p_bp, err_bp, dup_bp, hap_bp, col_bp"
      }
    }
  }

  parameter_meta {
    sample_id: {
      description: "Sample ID (output prefix)"
    }
    hap_fastas: {
      description: "The haplotype FASTAs of the assembly contract, in haplotype_labels order"
    }
    haplotype_labels: {
      description: "The haplotype labels (hap1/hap2, or pat/mat)"
    }
    flagger_bed: {
      description: "HMM-Flagger's prediction BED (assembly coordinates, label in column 4)"
    }
    telomere_motif: {
      description: "The telomere motif hifiasm was given (CCCTAA); empty disables the telomere columns"
    }
    cpu: {
      description: "Cores"
    }
    mem_gb: {
      description: "Memory (GB)"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    String sample_id
    Array[File] hap_fastas
    Array[String] haplotype_labels
    File flagger_bed
    String telomere_motif = "CCCTAA"
    Int cpu = 2
    Int mem_gb = 8
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(hap_fastas, "GB") * 2 + 10)

  command <<<
    set -euo pipefail
    SEQTK=/home/apps/seqtk/seqtk
    mkdir -p parts
    i=0
    for fa in ~{sep=" " hap_fastas}; do
      $SEQTK comp "$fa" > "parts/$i.comp.tsv"
      $SEQTK gap -l 10 "$fa" > "parts/$i.gap.bed"
      if [ -n "~{telomere_motif}" ]; then $SEQTK telo -m "~{telomere_motif}" "$fa" > "parts/$i.telo.bed"; else : > "parts/$i.telo.bed"; fi
      i=$((i + 1))
    done
    # Flagger's bases per contig and label (the BED has a track line; column 4 is Err, Dup, Hap or Col)
    awk -F'\t' '$0 !~ /^(track|#)/ && NF >= 4 { bp[$1 "\t" $4] += $3 - $2 } END { for (k in bp) print k "\t" bp[k] }' "~{flagger_bed}" > flagger_bp.tsv
    python3 - "~{sample_id}" "~{write_lines(haplotype_labels)}" <<'PY'
import collections, glob, sys
sample_id, labels_f = sys.argv[1:3]
labels = [l.strip() for l in open(labels_f) if l.strip()]
flag = collections.defaultdict(dict)
for line in open("flagger_bp.tsv"):
    f = line.rstrip("\n").split("\t")
    if len(f) == 3:
        flag[f[0]][f[1].lower()] = int(f[2])
with open(sample_id + ".asm.contig_qc.tsv", "w") as out:
    out.write("haplotype\tcontig\tlength\tgc_pct\tn_bp\tgaps\ttelo_5p_bp\ttelo_3p_bp\terr_bp\tdup_bp\thap_bp\tcol_bp\n")
    for i, label in enumerate(labels):
        gaps = collections.Counter()
        for line in open(f"parts/{i}.gap.bed"):
            f = line.split("\t")
            if len(f) >= 3:
                gaps[f[0]] += 1
        telo5, telo3 = collections.Counter(), collections.Counter()
        lengths = {}
        rows = []
        for line in open(f"parts/{i}.comp.tsv"):
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            name, length = f[0], int(f[1])
            a, c, g, t, n4 = int(f[2]), int(f[3]), int(f[4]), int(f[5]), int(f[8])
            acgt = a + c + g + t
            lengths[name] = length
            rows.append((name, length, round((c + g) / acgt * 100, 2) if acgt else 0.0, n4))
        for line in open(f"parts/{i}.telo.bed"):
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            try:
                start, end = int(f[1]), int(f[2])
            except ValueError:
                continue
            if start == 0:
                telo5[f[0]] += end - start
            if end >= lengths.get(f[0], -1):
                telo3[f[0]] += end - start
        for name, length, gc, n_bp in rows:
            fl = flag.get(name, {})
            out.write("\t".join(str(x) for x in [label, name, length, gc, n_bp, gaps[name], telo5[name], telo3[name],
                                                   fl.get("err", 0), fl.get("dup", 0), fl.get("hap", 0), fl.get("col", 0)]) + "\n")
PY
    awk -F'\t' 'NR > 1 { n[$1]++; telo[$1] += ($7 > 0) + ($8 > 0); nonhap[$1] += $9 + $10 + $12 } END { for (h in n) print h ": " n[h] " contigs, " telo[h] " telomeric ends, " nonhap[h] " non-Hap bp" }' "~{sample_id}.asm.contig_qc.tsv" >&2
  >>>

  output {
    File contig_qc = "~{sample_id}.asm.contig_qc.tsv"
  }

  runtime {
    docker: "mobinasri/flagger@sha256:3571bbe1c5415126edd0e8e646c853cc17f390bf1e2eb7c03b283ae54bc8262c"  # v1.2.0 (seqtk 1.4 at /home/apps/seqtk)
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
