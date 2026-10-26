version 1.0
# ugc-pacbio-wgw — what each haplotype assembled per chromosome, the karyotype inferred from chrX/chrY, and dipcall's regrouped haplotypes for non-trio XY samples
# origin: none
# see docs/DESIGN.md §5, §13; docs/ENTRYPOINTS.md §8.1
import "../../../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/structs.wdl"

task ugc_wgw_chrom_coverage {
  meta {
    description: "From the haplotypes' PAFs (primary alignments at or above min_mapq): aligned bases and contig counts per haplotype and chromosome; the karyotype (XX, XY or unknown) from chrX and chrY bases outside the PARs against two fractions of the chromosome lengths; and, for a non-trio XY sample, copies of the haplotype FASTAs with chrY contigs in hap1 and non-PAR chrX contigs in hap2, which is what dipcall's male mode assumes. The Python between the begin and end markers is tested on a synthetic PAF (tests/driver/test_chrom_coverage.py)."
    outputs: {
      coverage_tsv: {
        description: "haplotype, chrom, aligned_bases, contigs, chrom_length, fraction (one row per haplotype and chromosome)"
      },
      karyotype_json: {
        description: "The inferred karyotype with the chrX/chrY bases, fractions, thresholds and the reason; regrouped and moved_contigs"
      },
      inferred_karyotype: {
        description: "XX, XY or unknown"
      },
      regrouped_fastas: {
        description: "hap1 and hap2 copies for dipcall with chrX/chrY contigs regrouped; empty unless the sample is XY and not trio-binned"
      },
      regroup_log: {
        description: "One line per contig moved between the dipcall copies (empty when nothing moved)"
      }
    }
  }

  parameter_meta {
    sample_id: {
      description: "Sample ID; output prefix"
    }
    pafs: {
      description: "Contig-to-reference PAFs, one per haplotype, in haplotype order (gzipped or plain)"
    }
    hap_fastas: {
      description: "The haplotype FASTAs in the same order (regrouped for dipcall when needed)"
    }
    haplotype_labels: {
      description: "Haplotype labels in the same order"
    }
    ref_index: {
      description: "FASTA index of the assembly reference (chromosome lengths)"
    }
    par_bed: {
      description: "Pseudoautosomal regions (dipcall's hs38.PAR.bed): bases inside them count for neither chrX nor chrY"
    }
    trio: {
      description: "Whether the haplotypes are trio-binned (hap1 paternal): then nothing is regrouped"
    }
    min_mapq: {
      description: "Alignments below this MAPQ are ignored"
    }
    chrx_min_fraction: {
      description: "Below this fraction of chrX (non-PAR bases aligned, both haplotypes summed, over the chrX length) the karyotype is unknown"
    }
    chry_min_fraction: {
      description: "At or above this fraction of chrY (non-PAR bases aligned over the chrY length) the karyotype is XY, else XX"
    }
    threads: {
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
    Array[File] pafs
    Array[File] hap_fastas
    Array[String] haplotype_labels
    File ref_index
    File par_bed
    Boolean trio
    Int min_mapq = 5
    Float chrx_min_fraction = 0.5
    Float chry_min_fraction = 0.1
    Int threads = 2
    Int mem_gb = 8
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(hap_fastas, "GB") * 2 + size(pafs, "GB") + 10)

  command <<<
    set -euo pipefail
    python3 --version >&2
    mkdir -p regrouped
    cat > chrom_coverage.py <<'PY'
# ---- chrom_coverage.py (begin) ----
"""Per-chromosome assembly coverage, karyotype and dipcall regrouping from contig-to-reference PAFs (stdlib only).

usage: chrom_coverage.py <sample_id> <ref.fai> <par.bed> <trio:true|false> <min_mapq> <chrx_min_fraction>
                         <chry_min_fraction> <coverage.tsv> <karyotype.json> <karyotype.txt> <regroup.log>
                         <regrouped dir> <label>=<paf>[=<fasta>] ...
Primary alignments (tp:A:P or tp:A:I) at or above min_mapq count. Bases inside the PAR intervals count for
neither chrX nor chrY. A contig's home is the chromosome with most of its aligned bases; for a non-trio XY
sample, contigs mostly on chrY go to the hap1 copy and contigs mostly on non-PAR chrX to the hap2 copy.
"""
import collections
import gzip
import json
import sys


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_par(path):
    par = collections.defaultdict(list)
    with open(path) as fh:
        for line in fh:
            f = line.split()
            if len(f) >= 3 and not line.startswith("#"):
                par[f[0]].append((int(f[1]), int(f[2])))
    return par


def inside_par(par, chrom, start, end):
    """Bases of [start, end) that fall inside the PAR intervals of chrom."""
    return sum(max(0, min(end, e) - max(start, s)) for s, e in par.get(chrom, ()))


def read_paf(path, min_mapq):
    """Yield (qname, qlen, tname, tlen, tstart, tend) of the primary alignments."""
    with opener(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            tp = next((t[5:] for t in f[12:] if t.startswith("tp:A:")), "P")
            if tp not in ("P", "I") or int(f[11]) < min_mapq:
                continue
            yield f[0], int(f[1]), f[5], int(f[6]), int(f[7]), int(f[8])


def main(argv):
    (sample, fai, par_path, trio, min_mapq, x_min, y_min, cov_out, json_out, kar_out, log_out, regroup_dir) = argv[:12]
    trio, min_mapq, x_min, y_min = trio == "true", int(min_mapq), float(x_min), float(y_min)
    haps = []
    for spec in argv[12:]:
        parts = spec.split("=")
        haps.append((parts[0], parts[1], parts[2] if len(parts) > 2 else None))
    lengths = {}
    with open(fai) as fh:
        for line in fh:
            f = line.split("\t")
            lengths[f[0]] = int(f[1])
    par = read_par(par_path)
    cov = []                                   # (label, chrom, bases, contigs)
    contig_home = []                           # per haplotype: contig -> (home chrom, share, qlen, y_bases, x_nonpar)
    x_bases = y_bases = 0
    for label, paf, _ in haps:
        bases = collections.Counter()
        contigs = collections.defaultdict(set)
        per_contig = collections.defaultdict(collections.Counter)
        qlen = {}
        nonpar = collections.defaultdict(collections.Counter)   # contig -> chrom -> bases outside PAR
        for qname, ql, tname, tlen, ts, te in read_paf(paf, min_mapq):
            lengths.setdefault(tname, tlen)
            n = te - ts
            bases[tname] += n
            contigs[tname].add(qname)
            per_contig[qname][tname] += n
            nonpar[qname][tname] += n - inside_par(par, tname, ts, te)
            qlen[qname] = ql
        for chrom in sorted(bases, key=lambda c: (len(c), c)):
            cov.append((label, chrom, bases[chrom], len(contigs[chrom])))
        homes = {}
        for q, counts in per_contig.items():
            home, best = counts.most_common(1)[0]
            homes[q] = (home, best / sum(counts.values()), qlen[q], nonpar[q].get("chrY", 0), nonpar[q].get("chrX", 0))
            x_bases += nonpar[q].get("chrX", 0)
            y_bases += nonpar[q].get("chrY", 0)
        contig_home.append(homes)
    x_len, y_len = lengths.get("chrX"), lengths.get("chrY")
    x_frac = x_bases / x_len if x_len else 0.0
    y_frac = y_bases / y_len if y_len else 0.0
    if not x_len or not y_len:
        karyotype, reason = "unknown", "the reference index names no chrX or chrY"
    elif x_frac < x_min:
        karyotype, reason = "unknown", f"non-PAR chrX bases aligned are {x_frac:.3f} of the chrX length, below {x_min}"
    elif y_frac >= y_min:
        karyotype, reason = "XY", f"non-PAR chrY bases aligned are {y_frac:.3f} of the chrY length, at or above {y_min}"
    else:
        karyotype, reason = "XX", f"non-PAR chrY bases aligned are {y_frac:.3f} of the chrY length, below {y_min}"
    moved = []
    regrouped = karyotype == "XY" and not trio and len(haps) == 2 and all(h[2] for h in haps)
    if regrouped:
        dest = [{}, {}]                         # per source haplotype: contig -> destination index
        for i, homes in enumerate(contig_home):
            for q, (home, share, ql, yb, xb) in homes.items():
                if share <= 0.5:
                    continue
                if home == "chrY" and i != 0:
                    dest[i][q] = 0
                elif home == "chrX" and xb > 0 and i != 1:
                    dest[i][q] = 1
                if q in dest[i]:
                    moved.append((haps[i][0], q, ql, haps[dest[i][q]][0], home, yb, xb))
        # gzip level 4, bytes unparsed: the default level 9 took an hour and a half for two 3 Gb haplotypes
        # on the cluster (2026-10-26); level 4 is 25 times faster for a fifth more bytes (chr20 pair: 5.7 s,
        # 43 MB against 135 s, 39 MB; bgzip 10 s, 36 MB, but it is not in every environment the tests run in)
        outs = [gzip.open(f"{regroup_dir}/{sample}.dipcall.hap{k + 1}.fasta.gz", "wb", compresslevel=4) for k in range(2)]
        for i, (_, _, fasta) in enumerate(haps):
            with (gzip.open(fasta, "rb") if fasta.endswith(".gz") else open(fasta, "rb")) as fh:
                target = outs[i]
                for line in fh:
                    if line.startswith(b">"):
                        name = line[1:].split()[0].decode()
                        target = outs[dest[i].get(name, i)]
                    target.write(line)
        for o in outs:
            o.close()
    with open(cov_out, "w") as fh:
        fh.write("haplotype\tchrom\taligned_bases\tcontigs\tchrom_length\tfraction\n")
        for label, chrom, b, n in cov:
            ln = lengths.get(chrom, 0)
            fh.write(f"{label}\t{chrom}\t{b}\t{n}\t{ln}\t{(b / ln if ln else 0):.4f}\n")
    with open(log_out, "w") as fh:
        fh.write("from_haplotype\tcontig\tlength\tto_haplotype\thome_chrom\tchrY_bases\tchrX_nonpar_bases\n")
        for row in moved:
            fh.write("\t".join(str(x) for x in row) + "\n")
    doc = {"sample_id": sample, "karyotype": karyotype, "reason": reason, "trio": trio,
           "chrx_bases_outside_par": x_bases, "chry_bases_outside_par": y_bases, "chrx_length": x_len, "chry_length": y_len,
           "chrx_fraction": round(x_frac, 4), "chry_fraction": round(y_frac, 4),
           "thresholds": {"min_mapq": min_mapq, "chrx_min_fraction": x_min, "chry_min_fraction": y_min},
           "regrouped": regrouped, "moved_contigs": len(moved), "haplotypes": [h[0] for h in haps]}
    with open(json_out, "w") as fh:
        json.dump(doc, fh, indent=2, sort_keys=True)
        fh.write("\n")
    with open(kar_out, "w") as fh:
        fh.write(karyotype + "\n")
    print(f"karyotype {karyotype}: {reason}; regrouped={regrouped} moved={len(moved)}", file=sys.stderr)


if __name__ == "__main__":
    main(sys.argv[1:])
# ---- chrom_coverage.py (end) ----
PY
    paste -d '=' '~{write_lines(haplotype_labels)}' '~{write_lines(pafs)}' '~{write_lines(hap_fastas)}' > specs.txt
    # shellcheck disable=SC2046
    python3 chrom_coverage.py "~{sample_id}" "~{ref_index}" "~{par_bed}" ~{true="true" false="false" trio} ~{min_mapq} \
      ~{chrx_min_fraction} ~{chry_min_fraction} "~{sample_id}.asm.chrom_coverage.tsv" "~{sample_id}.asm.karyotype.json" \
      karyotype.txt "~{sample_id}.asm.regroup.tsv" regrouped $(cat specs.txt)
  >>>

  output {
    File coverage_tsv = "~{sample_id}.asm.chrom_coverage.tsv"
    File karyotype_json = "~{sample_id}.asm.karyotype.json"
    String inferred_karyotype = read_string("karyotype.txt")
    Array[File] regrouped_fastas = glob("regrouped/*.fasta.gz")
    File regroup_log = "~{sample_id}.asm.regroup.tsv"
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
