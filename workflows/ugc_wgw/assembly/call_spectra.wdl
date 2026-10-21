version 1.0
# ugc-pacbio-wgw — the dipcall and PAV calls of one sample counted by type, size, FILTER class and contig class
# origin: none
# see docs/DESIGN.md §5, §13; docs/ENTRYPOINTS.md §8.1
import "../../../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/structs.wdl"

task ugc_wgw_asm_call_spectra {
  meta {
    description: "bcftools query over both primary VCFs, one count per (caller, filter class, type, size bin, contig class): type from REF/ALT and INFO/SVTYPE (SNV, MNV, INDEL below 50 bp, INS and DEL from 50 bp, INV, OTHER), size bins as the SV sections use them, filter class PASS, ASM_UNRELIABLE (Flagger's mask alone), OTHER (another filter) or BOTH; the summaries draw the spectra and the PASS-versus-masked split from it."
    outputs: {
      call_spectra: {
        description: "<sample_id>.asm.call_spectra.tsv: caller, filter_class, type, size_bin, contig_class, count"
      }
    }
  }

  parameter_meta {
    sample_id: {
      description: "Sample ID (output prefix)"
    }
    dipcall_vcf: {
      description: "The primary dipcall VCF"
    }
    dipcall_vcf_index: {
      description: "Its index"
    }
    pav_vcf: {
      description: "The primary PAV VCF"
    }
    pav_vcf_index: {
      description: "Its index"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    String sample_id
    File dipcall_vcf
    File dipcall_vcf_index
    File pav_vcf
    File pav_vcf_index
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil((size(dipcall_vcf, "GB") + size(pav_vcf, "GB")) * 2 + 10)

  command <<<
    set -euo pipefail
    bcftools --version | head -1 >&2
    query() {  # query <vcf> <caller>: CHROM FILTER REF ALT SVTYPE SVLEN (dots when the tag is not in the header)
      fmt='%CHROM\t%FILTER\t%REF\t%ALT'
      for t in SVTYPE SVLEN; do
        if bcftools view --header-only "$1" | grep -q "^##INFO=<ID=$t,"; then fmt="$fmt\t%INFO/$t"; else fmt="$fmt\t."; fi
      done
      bcftools query --format "$fmt\n" "$1" | sed "s/^/$2\t/"
    }
    { query "~{dipcall_vcf}" dipcall; query "~{pav_vcf}" pav; } > calls.tsv
    python3 - "~{sample_id}" <<'PY'
import collections, sys
sample_id = sys.argv[1]
edges = [(100, "<100"), (1000, "100-1k"), (10000, "1k-10k"), (100000, "10k-100k"), (1000000, "100k-1M")]
def contig_class(c):
    c = c[3:] if c.startswith("chr") else c
    return "autosome" if c.isdigit() else "chrX" if c == "X" else "chrY" if c == "Y" else "chrM" if c in ("M", "MT") else "other"
def size_bin(n):
    for edge, label in edges:
        if n < edge:
            return label
    return ">1M"
counts = collections.Counter()
for line in open("calls.tsv"):
    f = line.rstrip("\n").split("\t")
    if len(f) < 7:
        continue
    caller, chrom, filt, ref, alt, svtype, svlen = f[:7]
    filters = set(x for x in filt.split(";") if x and x != ".")
    masked = "ASM_UNRELIABLE" in filters
    others = filters - {"PASS", "ASM_UNRELIABLE"}
    fclass = "BOTH" if masked and others else "ASM_UNRELIABLE" if masked else "OTHER" if others else "PASS"
    alts = [a for a in alt.split(",") if a]
    types = svtype.split(",") if svtype not in ("", ".") else []
    lens = []
    for x in svlen.split(","):
        try:
            lens.append(abs(int(x)))
        except ValueError:
            pass
    kinds = set()
    for i, a in enumerate(alts):
        t = types[i] if i < len(types) else (types[0] if types else "")
        n = lens[i] if i < len(lens) else (lens[0] if lens else None)
        if a.startswith("<") or t in ("INV", "DUP", "BND", "CPX"):
            sym = a.strip("<>").split(":")[0] if a.startswith("<") else t
            kind = "INV" if sym == "INV" else "DEL" if sym == "DEL" else "INS" if sym in ("INS", "DUP") else "OTHER"
            size = n if n is not None else 0
        else:
            d = len(a) - len(ref)
            size = abs(d)
            if d == 0:
                kind = "SNV" if len(ref) == 1 else "MNV"
            elif size < 50:
                kind = "INDEL"
            else:
                kind = "INS" if d > 0 else "DEL"
        kinds.add((kind, size))
    if not kinds:
        continue
    kind, size = max(kinds, key=lambda ks: ks[1])
    bin_label = "<50" if kind in ("SNV", "MNV", "INDEL") else size_bin(size)
    counts[(caller, fclass, kind, bin_label, contig_class(chrom))] += 1
with open(sample_id + ".asm.call_spectra.tsv", "w") as out:
    out.write("caller\tfilter_class\ttype\tsize_bin\tcontig_class\tcount\n")
    for key in sorted(counts):
        out.write("\t".join(key) + "\t" + str(counts[key]) + "\n")
by_caller = collections.Counter()
for (caller, fclass, kind, _b, _c), n in counts.items():
    by_caller[(caller, fclass)] += n
print({f"{c} {fc}": n for (c, fc), n in sorted(by_caller.items())}, file=sys.stderr)
PY
  >>>

  output {
    File call_spectra = "~{sample_id}.asm.call_spectra.tsv"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:03cb3c01937eccc907f8ad71c87b258581504572205fe3f31a657e318f3564ae"  # pb_wdl_base:build4
    cpu: 1
    memory: "4 GiB"
    disk: "~{disk_size} GB"
    disks: "local-disk ~{disk_size} HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
    cpuPlatform: runtime_attributes.cpuPlatform
  }
}
