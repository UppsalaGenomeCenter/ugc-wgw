version 1.0
# ugc-pacbio-wgw — the cohort pangenome's VCFs gathered over chromosomes, renamed to the sample IDs, unphased for non-trio members, and their frequency summary
# origin: none
# see docs/DESIGN.md §5, §13; docs/ENTRYPOINTS.md §8.2
import "../../../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/structs.wdl"

task ugc_wgw_pangenome_vcf_post {
  meta {
    description: "The per-chromosome VCFs of one kind (raw or vcfwave) brought to the cohort's sample columns (a member without sequence in a chromosome's graph, an XX member on chrY, has no column in that VCF: it is added with missing genotypes, every column in cohort order), bcftools concat, the Cactus-safe sample names changed back to the sample IDs (bcftools reheader), genotypes of the members that are not trio-binned unphased (| to / in those columns only: hifiasm's hap1/hap2 are local phase blocks, so the .1|.2 order means nothing beyond them), header lines about the cohort, the ambiguity rule and the phase policy added, bgzip, tabix."
    outputs: {
      vcf: {
        description: "<out_prefix>.vcf.gz"
      },
      vcf_index: {
        description: "Tabix index"
      }
    }
  }

  parameter_meta {
    vcfs: {
      description: "Per-chromosome VCFs of one kind, any order (concat sorts by the header's contig order)"
    }
    vcf_indices: {
      description: "Their indices"
    }
    rename_samples: {
      description: "safe_name<TAB>sample_id per member (bcftools reheader --samples)"
    }
    unphase_samples: {
      description: "Sample IDs whose genotypes are unphased"
    }
    header_lines: {
      description: "Header lines to add"
    }
    chromosomes: {
      description: "The chromosomes built (a header line)"
    }
    out_prefix: {
      description: "Output prefix"
    }
    threads: {
      description: "Threads"
    }
    mem_gb: {
      description: "Memory (GB)"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    Array[File] vcfs
    Array[File] vcf_indices
    File rename_samples
    Array[String] unphase_samples
    Array[String] header_lines
    Array[String] chromosomes
    String out_prefix
    Int threads = 4
    Int mem_gb = 8
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(vcfs, "GB") * 6 + size(vcf_indices, "GB") + 10)

  command <<<
    set -euo pipefail
    bcftools --version | head -1 >&2
    mkdir -p in
    for f in ~{sep=" " vcfs} ~{sep=" " vcf_indices}; do ln -s "$f" "in/$(basename "$f")"; done
    # a chromosome's graph holds only the members with contigs assigned to it, and Cactus writes no sample column
    # for the others (an XX member on chrY; a member whose assembly left a chromosome out); bcftools concat wants
    # identical sample columns, so such a VCF gets the absent members with missing genotypes, every column in
    # cohort order (the Cactus-safe names of rename_samples, which is that order)
    cut -f1 "~{rename_samples}" > samples.txt
    cat > pad.py <<'PY'
import sys
want = [l.strip() for l in open(sys.argv[1]) if l.strip()]
order = []
out = sys.stdout
for line in sys.stdin:
    if line.startswith("##"):
        out.write(line)
        continue
    if line.startswith("#CHROM"):
        names = line.rstrip("\n").split("\t")
        have = {n: i for i, n in enumerate(names[9:], start=9)}
        unknown = sorted(set(have) - set(want))
        if unknown:
            sys.exit("samples in the VCF but not in the cohort: " + ",".join(unknown))
        order = [have.get(n) for n in want]
        print(f"{len(have)} of {len(want)} members in this VCF; {len(want) - len(have)} added with missing genotypes", file=sys.stderr)
        out.write("\t".join(names[:9] + want) + "\n")
        continue
    f = line.rstrip("\n").split("\t")
    empty = ":".join("." for _ in f[8].split(":"))
    out.write("\t".join(f[:9] + [f[i] if i is not None else empty for i in order]) + "\n")
PY
    for f in in/*.vcf.gz; do
      if ! cmp -s <(bcftools query --list-samples "$f") samples.txt; then
        echo "$(basename "$f"): sample columns differ from the cohort's, padding" >&2
        bcftools view --output-type v "$f" | python3 pad.py samples.txt | bgzip --threads ~{threads} > "$f.padded"
        rm -f "$f" "$f.tbi" "$f.csi"
        mv "$f.padded" "$f"
        tabix --preset vcf "$f"
      fi
    done
    ls in/*.vcf.gz > list.txt
    cp "~{write_lines(header_lines)}" header.txt
    echo "##ugc_wgw_cohort_pangenome_chromosomes=~{sep=',' chromosomes}" >> header.txt
    # concat, the sample IDs back, the header lines; then the genotypes of the members that are not trio-binned
    # unphased (| to / in GT, those columns only): bcftools +setGT has no per-sample selection (its -s/-S subset
    # the output), so a stream filter does it
    cat > unphase.py <<'PY'
import sys
want = {l.strip() for l in open(sys.argv[1]) if l.strip()}
cols = []
out = sys.stdout
for line in sys.stdin:
    if line.startswith("#"):
        if line.startswith("#CHROM"):
            names = line.rstrip("\n").split("\t")
            cols = [i for i, n in enumerate(names) if i >= 9 and n in want]
            missing = want - set(names[9:])
            if missing:
                sys.exit("samples to unphase not in the VCF: " + ",".join(sorted(missing)))
            print(f"unphasing {len(cols)} of {len(names) - 9} samples", file=sys.stderr)
        out.write(line)
        continue
    if not cols:
        out.write(line)
        continue
    f = line.rstrip("\n").split("\t")
    for i in cols:
        gt, sep, rest = f[i].partition(":")
        f[i] = gt.replace("|", "/") + sep + rest
    out.write("\t".join(f) + "\n")
PY
    bcftools concat --threads ~{threads} --file-list list.txt --output-type u \
      | bcftools reheader --samples "~{rename_samples}" \
      | bcftools annotate --header-lines header.txt --output-type v - \
      | python3 unphase.py "~{write_lines(unphase_samples)}" \
      | bgzip --threads ~{threads} > "~{out_prefix}.vcf.gz"
    tabix --preset vcf "~{out_prefix}.vcf.gz"
    bcftools query --list-samples "~{out_prefix}.vcf.gz" | tr '\n' ' ' >&2; echo >&2
  >>>

  output {
    File vcf = "~{out_prefix}.vcf.gz"
    File vcf_index = "~{out_prefix}.vcf.gz.tbi"
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

task ugc_wgw_pangenome_freq_summary {
  meta {
    description: "The frequency VCF of the graph (sites only, AC/AN/AF/NS and the XX/XY strata from bcftools +fill-tags) read once into the cohort_freq summary shape (resource pangenome; metric, key, value): members by stratum and sex source, records by type (SNV, MNV, INDEL below 50 bp, INS/DEL from 50 bp, MIXED; by the allele length difference), by contig class, by FILTER, by nesting level (INFO/LV of vg deconstruct: 0 top level, 1, 2+), by vcfwave's TYPE and by its allele LEN bins (the records vcfwave realigned carry them; the rest count as '.'), call rate (NS over members; a missing genotype means no assembled sequence at the site), the cohort_freq allele-frequency bins, and records with carriers per stratum; plus the samples table."
    outputs: {
      summary: {
        description: "<out_prefix>.summary.tsv"
      },
      samples: {
        description: "<out_prefix>.samples.tsv: sample_id, group, sex_source, in_pangenome_vcf"
      }
    }
  }

  parameter_meta {
    freq_vcf: {
      description: "The concatenated frequency VCF (sites only)"
    }
    freq_vcf_index: {
      description: "Its index"
    }
    vcf_samples: {
      description: "The sample names of the graph VCF the frequencies were computed from"
    }
    sample_ids: {
      description: "Members in cohort order"
    }
    sample_sexes: {
      description: "XX, XY or empty per member"
    }
    sample_sex_sources: {
      description: "Where each member's sex came from (karyotype, or empty)"
    }
    out_prefix: {
      description: "Output prefix"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    File freq_vcf
    File freq_vcf_index
    File vcf_samples
    Array[String] sample_ids
    Array[String] sample_sexes
    Array[String] sample_sex_sources
    String out_prefix
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(freq_vcf, "GB") * 2 + 10)

  command <<<
    set -euo pipefail
    bcftools --version | head -1 >&2
    has_tag() { bcftools view --header-only "~{freq_vcf}" | grep -q "^##INFO=<ID=$1,"; }
    fmt='%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO/AC\t%INFO/AN\t%INFO/AF\t%INFO/NS'
    for t in AC_XX AC_XY LV TYPE LEN; do if has_tag $t; then fmt="$fmt\t%INFO/$t"; else fmt="$fmt\t."; fi; done
    bcftools query --format "$fmt\n" "~{freq_vcf}" > sites.tsv
    python3 - "~{out_prefix}" sites.tsv "~{vcf_samples}" "~{write_lines(sample_ids)}" "~{write_lines(sample_sexes)}" "~{write_lines(sample_sex_sources)}" <<'PY'
import collections, sys
out_prefix, sites, vcf_samples, ids_f, sexes_f, sources_f = sys.argv[1:7]
def lines(p): return [l.rstrip("\n") for l in open(p)]
ids, sexes, sources = lines(ids_f), lines(sexes_f), lines(sources_f)
sources += [""] * (len(ids) - len(sources))
in_vcf = set(lines(vcf_samples))
rows = []
def add(resource, metric, key, value): rows.append((resource, metric, key, value))
n = len(ids)
add("samples", "n", "total", n)
for g in ("XX", "XY"): add("samples", "n", g, sum(1 for s in sexes if s == g))
add("samples", "n", "unknown", sum(1 for s in sexes if s not in ("XX", "XY")))
for src in ("karyotype", "none"):
    add("samples", "n_sex_source", src, sum(1 for s in sources if (s or "none") == src))
add("samples", "n_in_vcf", "pangenome", sum(1 for s in ids if s in in_vcf))
def contig_class(c):
    if c in ("chrX", "X"): return "chrX"
    if c in ("chrY", "Y"): return "chrY"
    if c in ("chrM", "MT", "M"): return "chrM"
    return "autosome" if c.replace("chr", "").isdigit() else "other"
def nums(field):
    out = []
    for x in field.split(","):
        try: out.append(float(x))
        except ValueError: pass
    return out
by_type = collections.Counter(); by_class = collections.Counter(); by_filter = collections.Counter(); bins = collections.Counter()
by_level = collections.Counter(); by_vg_type = collections.Counter(); by_len = collections.Counter()
len_edges = [(2, "1"), (10, "2-9"), (50, "10-49"), (1000, "50-999"), (10000, "1k-9.9k")]
strata = {"XX": 0, "XY": 0}
total = ac_gt0 = singletons = low_ns = ns_n = 0; ns_sum = 0.0
for line in open(sites):
    f = line.rstrip("\n").split("\t")
    if len(f) < 11: continue
    chrom, pos, ref, alt, filt, ac, an, af, ns, acxx, acxy = f[:11]
    lv, vg_type, vg_len = (f[11:14] + [".", ".", "."])[:3]
    by_level["." if lv in ("", ".") else ("2+" if lv.isdigit() and int(lv) >= 2 else lv)] += 1
    for t in (vg_type.split(",") if vg_type not in ("", ".") else ["."]): by_vg_type[t] += 1
    lens = nums(vg_len)
    if lens:
        m = max(abs(x) for x in lens)
        by_len[next((label for edge, label in len_edges if m < edge), ">=10k")] += 1
    total += 1
    alts = alt.split(",")
    kinds = set()
    for a in alts:
        d = len(a) - len(ref)
        kinds.add("SNV" if d == 0 and len(ref) == 1 else ("MNV" if d == 0 else ("INDEL" if abs(d) < 50 else ("INS" if d > 0 else "DEL"))))
    by_type["MIXED" if len(kinds) > 1 else kinds.pop()] += 1
    by_class[contig_class(chrom)] += 1
    by_filter[filt] += 1
    total_ac = sum(nums(ac)); afs = nums(af)
    if total_ac > 0: ac_gt0 += 1
    if total_ac == 1: singletons += 1
    if afs and total_ac > 0:
        m = max(afs)
        bins["<0.001" if m < 0.001 else "<0.01" if m < 0.01 else "<0.05" if m < 0.05 else ">=0.05"] += 1
    nsv = nums(ns)
    if nsv:
        ns_sum += nsv[0]; ns_n += 1
        if n and nsv[0] < 0.9 * n: low_ns += 1
    if sum(nums(acxx)) > 0: strata["XX"] += 1
    if sum(nums(acxy)) > 0: strata["XY"] += 1
r = "pangenome"
add(r, "records", "total", total)
for k, v in sorted(by_type.items()): add(r, "records_by_type", k, v)
for k, v in sorted(by_class.items()): add(r, "records_by_contig_class", k, v)
for k, v in sorted(by_filter.items()): add(r, "records_by_filter", k, v)
for k, v in sorted(by_level.items()): add(r, "records_by_level", k, v)
for k, v in sorted(by_vg_type.items()): add(r, "records_by_vg_type", k, v)
for k in [label for _e, label in len_edges] + [">=10k"]:
    if by_len.get(k): add(r, "records_by_len_bin", k, by_len[k])
add(r, "records", "ac_gt0", ac_gt0)
add(r, "records", "singletons", singletons)
for k in ("<0.001", "<0.01", "<0.05", ">=0.05"): add(r, "records_by_af_bin", k, bins[k])
for g in ("XX", "XY"): add(r, "records_ac_gt0_in", g, strata[g])
add(r, "mean_call_rate", "all", "%.4f" % (ns_sum / ns_n / n) if ns_n and n else "NA")
add(r, "records", "call_rate_below_0.9", low_ns)
with open(out_prefix + ".summary.tsv", "w") as fh:
    fh.write("# ugc-pacbio-wgw cohort_pangenome: AF bins use the largest allele frequency of a record with AC > 0; "
             "call rate = NS / n_samples, a missing genotype means no assembled sequence at the site; "
             "types by the allele length difference (INDEL below 50 bp, INS/DEL from 50 bp)\n")
    fh.write("resource\tmetric\tkey\tvalue\n")
    for row in rows: fh.write("\t".join(str(x) for x in row) + "\n")
with open(out_prefix + ".samples.tsv", "w") as fh:
    fh.write("sample_id\tgroup\tsex_source\tin_pangenome_vcf\n")
    for sid, sex, src in zip(ids, sexes, sources):
        fh.write("\t".join([sid, sex if sex in ("XX", "XY") else "unknown", src or "none", "yes" if sid in in_vcf else "no"]) + "\n")
print(f"{total} records, {ac_gt0} with carriers, {dict(by_type)}", file=sys.stderr)
PY
  >>>

  output {
    File summary = "~{out_prefix}.summary.tsv"
    File samples = "~{out_prefix}.samples.tsv"
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
