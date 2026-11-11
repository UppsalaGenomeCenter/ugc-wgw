version 1.0

# ugc-pacbio-wgw — the cohort's tandem repeats from the merged TRGT VCF: per-locus allele-length statistics with the expansion-outlier scan, the screen of the named disease loci against STRchive's thresholds, and the gather of the shard tables
# origin: none (thresholds from STRchive, CC BY 4.0, through references/tr_disease_loci.*.tsv)
# see docs/DESIGN.md §5, §7.6, §18 (2026-11-11); docs/ENTRYPOINTS.md §10

import "../../../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/structs.wdl"

task ugc_wgw_tr_locus_stats {
  meta {
    description: "One shard of the cohort's merged TRGT VCF: per locus the allele-length distribution over the called alleles (median, MAD, quantiles), the allele spectrum, and the expansion outliers: a member's longer allele at least outlier_min_bp above the cohort median and at least outlier_mad_min scaled MADs (the scale never below one repeat unit), at a locus with enough called members, and only while the carriers are few (a common long allele is polymorphism: the locus is marked common and lists nobody). Contractions are not scanned."
    outputs: {
      loci: {
        description: "Per locus: chrom, start, end, trid, motifs, named_locus, ref_len, n_expected, n_called, n_alleles, n_distinct, al_min, al_median, al_mad, al_p95, al_p99, al_max, n_carriers, outlier_status (listed, common, none, below_min_called)"
      },
      spectrum: {
        description: "trid, al, n_alleles, n_members: every allele length seen at a locus"
      },
      outliers: {
        description: "One row per member and locus listed: the longer allele, the other one, spanning reads, excess in bp, repeat units and scaled MADs, the cohort's median and MAD, carriers, both_alleles"
      },
      member_counts: {
        description: "Per sample of the VCF: loci in the shard, loci uncalled, outliers at ordinary loci, outliers at named loci and their twins, outliers on both alleles"
      }
    }
  }

  parameter_meta {
    trgt_vcf: {
      description: "The cohort's merged TRGT VCF (cohort_merge)"
    }
    trgt_vcf_index: {
      description: "Index for trgt_vcf"
    }
    regions_bed: {
      description: "BED of the shard; a record belongs to the shard holding its POS (--regions-overlap 0)"
    }
    sample_ids: {
      description: "Cohort members; every one must be a sample of the VCF"
    }
    sample_sexes: {
      description: "XX, XY or empty per member, aligned with sample_ids; decides how many members a chrY locus expects"
    }
    named_overlaps: {
      description: "ordinary_trid, named_trid, overlap_bp: the ordinary loci that overlap a named disease locus (ugc-wgw-tr)"
    }
    outlier_min_bp: {
      description: "An outlier's longer allele is at least this many bp above the cohort median"
    }
    outlier_mad_min: {
      description: "... and at least this many scaled MADs above it (scale = max(1.4826 x MAD, the locus's shortest motif))"
    }
    outlier_min_called: {
      description: "Loci with fewer called members are not scanned"
    }
    outlier_max_carrier_fraction: {
      description: "A locus where more than this fraction of the called members qualify is common, not an outlier locus"
    }
    outlier_max_carriers_min: {
      description: "... but this many carriers are always allowed"
    }
    out_prefix: {
      description: "Output prefix"
    }
    threads: {
      description: "Threads (bcftools and the script, one each)"
    }
    mem_gb: {
      description: "Memory (GB); one record at a time"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    File trgt_vcf
    File trgt_vcf_index
    File regions_bed
    Array[String] sample_ids
    Array[String] sample_sexes
    File named_overlaps
    Int outlier_min_bp = 50
    Float outlier_mad_min = 5.0
    Int outlier_min_called = 20
    Float outlier_max_carrier_fraction = 0.01
    Int outlier_max_carriers_min = 2
    String out_prefix
    Int threads = 2
    Int mem_gb = 2
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(trgt_vcf, "GB") + 20)

  command <<<
    set -euo pipefail
    export LC_ALL=C

    bcftools --version >&2

    mkdir -p in
    ln --symbolic "~{trgt_vcf}" "in/~{basename(trgt_vcf)}"
    ln --symbolic "~{trgt_vcf_index}" "in/~{basename(trgt_vcf_index)}"

    # every cohort member must be a sample of the VCF, by name; the columns are read in the VCF's order
    bcftools query --list-samples "in/~{basename(trgt_vcf)}" > vcf_samples.txt
    sort vcf_samples.txt > vcf_samples.sorted.txt
    sort "~{write_lines(sample_ids)}" | comm -23 - vcf_samples.sorted.txt > missing.txt
    if [ -s missing.txt ]; then
      echo "cohort members absent from ~{basename(trgt_vcf)}: $(tr '\n' ' ' < missing.txt)" >&2
      exit 1
    fi
    paste "~{write_lines(sample_ids)}" "~{write_lines(sample_sexes)}" > sexes.tsv

    # the script goes to a file: its stdin is the query's output
    cat > tr_locus_stats.py <<'EOF'
    # ---- tr_locus_stats.py (begin) ----
    import math, sys
    from collections import Counter
    (samples_path, sexes_path, overlaps_path, prefix, min_bp, mad_min, min_called, max_fraction, max_carriers_min) = sys.argv[1:10]
    min_bp, mad_min, min_called, max_fraction, max_carriers_min = int(min_bp), float(mad_min), int(min_called), float(max_fraction), int(max_carriers_min)
    samples = [l.strip() for l in open(samples_path) if l.strip()]
    sex = dict(l.rstrip("\n").split("\t")[:2] for l in open(sexes_path) if "\t" in l)
    n = len(samples)
    n_not_xx = sum(1 for s in samples if sex.get(s, "") != "XX")
    twin_of = {}
    for l in open(overlaps_path):
        if l.startswith("#") or l.startswith("ordinary_trid"):
            continue
        c = l.rstrip("\n").split("\t")
        if len(c) >= 2:
            twin_of[c[0]] = c[1]

    def is_named(trid):
        p = trid.split("_")
        return not (len(p) >= 4 and p[0].startswith("chr") and p[1].isdigit() and p[2].isdigit())

    def wquantile(pairs, total, q):
        """Nearest-rank quantile of sorted (value, weight) pairs."""
        rank = max(1, math.ceil(q * total))
        acc = 0
        for v, w in pairs:
            acc += w
            if acc >= rank:
                return v
        return pairs[-1][0]

    def wmedian(pairs, total):
        """The median of sorted (value, weight) pairs, averaging the two middle values of an even total."""
        lo = wquantile(pairs, total, 0.5)
        if total % 2:
            return lo
        acc = 0
        for v, w in pairs:
            acc += w
            if acc >= total // 2 + 1:
                return (lo + v) / 2
        return lo

    def fmt(x):
        if isinstance(x, float):
            return str(int(x)) if x == int(x) else f"{x:.2f}".rstrip("0").rstrip(".")
        return str(x)

    LOCI = ["chrom", "start", "end", "trid", "motifs", "named_locus", "ref_len", "n_expected", "n_called", "n_alleles", "n_distinct", "al_min", "al_median",
            "al_mad", "al_p95", "al_p99", "al_max", "n_carriers", "outlier_status"]
    PAIRS = ["sample_id", "trid", "chrom", "start", "end", "motifs", "named_locus", "al", "al_other", "sd", "excess_bp", "excess_units", "z", "cohort_median",
             "cohort_mad", "scale", "ref_len", "n_called", "n_carriers", "both_alleles"]
    uncalled = [0] * n
    n_out = [0] * n
    n_out_named = [0] * n
    n_out_both = [0] * n
    n_loci = 0
    status_counts = Counter()
    with open(prefix + ".loci.tsv", "w") as loci, open(prefix + ".spectrum.tsv", "w") as spec, open(prefix + ".outliers.tsv", "w") as pairs_out:
        loci.write("\t".join(LOCI) + "\n")
        spec.write("trid\tal\tn_alleles\tn_members\n")
        pairs_out.write("\t".join(PAIRS) + "\n")
        for line in sys.stdin:
            f = line.rstrip("\n").split("\t")
            if len(f) != 5 + 2 * n:
                sys.exit(f"{f[3] if len(f) > 3 else line[:60]!r}: {len(f)} fields for {n} samples (expected {5 + 2 * n})")
            chrom, pos, end, trid, motifs = f[0], int(f[1]), int(f[2]), f[3], f[4]
            als = f[5:5 + n]
            n_loci += 1
            named = trid if is_named(trid) else twin_of.get(trid, "")
            expected = n_not_xx if chrom == "chrY" else n
            counts = Counter(als)
            if "." in counts:
                for i in [i for i, s in enumerate(als) if s == "."]:
                    uncalled[i] += 1
            n_called = n - counts.get(".", 0)
            row = {"chrom": chrom, "start": pos, "end": end, "trid": trid, "motifs": motifs, "named_locus": named, "ref_len": end - pos,
                   "n_expected": expected, "n_called": n_called}
            if n_called == 0:
                status_counts["below_min_called"] += 1
                row["outlier_status"] = "below_min_called"
                loci.write("\t".join(fmt(row.get(k, "")) for k in LOCI) + "\n")
                continue
            parsed = {s: [int(x) for x in s.split(",")] for s in counts if s != "."}
            alleles = Counter()
            members = Counter()
            for s, vals in parsed.items():
                w = counts[s]
                for v in vals:
                    alleles[v] += w
                for v in set(vals):
                    members[v] += w
            dist = sorted(alleles.items())
            total = sum(alleles.values())
            med = wmedian(dist, total)
            devw = Counter()
            for v, w in alleles.items():
                devw[abs(v - med)] += w
            mad = wmedian(sorted(devw.items()), total)
            motif_len = min(len(m) for m in motifs.split(","))
            scale = max(1.4826 * mad, motif_len)
            for v, w in dist:
                spec.write(f"{trid}\t{v}\t{w}\t{members[v]}\n")
            row.update({"n_alleles": total, "n_distinct": len(dist), "al_min": dist[0][0], "al_median": med, "al_mad": mad,
                        "al_p95": wquantile(dist, total, 0.95), "al_p99": wquantile(dist, total, 0.99), "al_max": dist[-1][0]})

            def qualifies(v):
                return v - med >= min_bp and (v - med) / scale >= mad_min

            cand = {s for s, vals in parsed.items() if qualifies(max(vals))}
            carriers = sum(counts[s] for s in cand)
            row["n_carriers"] = carriers
            if n_called < min_called:
                status = "below_min_called"
            elif not carriers:
                status = "none"
            elif carriers > max(max_carriers_min, int(max_fraction * n_called)):
                status = "common"
            else:
                status = "listed"
                sds = f[5 + n:]
                for i, s in enumerate(als):
                    if s not in cand:
                        continue
                    vals = parsed[s]
                    k = max(range(len(vals)), key=lambda j: vals[j])
                    longer = vals[k]
                    other = vals[1 - k] if len(vals) == 2 else ""
                    both = len(vals) == 2 and qualifies(other)
                    sd = sds[i].split(",")
                    excess = longer - med
                    pairs_out.write("\t".join(fmt(x) for x in (
                        samples[i], trid, chrom, pos, end, motifs, named, longer, other, sds[i] if both or len(sd) != len(vals) else sd[k], excess,
                        round(excess / motif_len, 1), round(excess / scale, 2), med, mad, round(scale, 2), end - pos, n_called, carriers,
                        "yes" if both else "no")) + "\n")
                    if is_named(trid) or named:
                        n_out_named[i] += 1
                    else:
                        n_out[i] += 1
                        if both:
                            n_out_both[i] += 1
            status_counts[status] += 1
            row["outlier_status"] = status
            loci.write("\t".join(fmt(row.get(k, "")) for k in LOCI) + "\n")
    with open(prefix + ".members.tsv", "w") as out:
        out.write("sample_id\tn_loci\tn_loci_uncalled\tn_outliers\tn_outliers_named\tn_outliers_both_alleles\n")
        for i, s in enumerate(samples):
            out.write(f"{s}\t{n_loci}\t{uncalled[i]}\t{n_out[i]}\t{n_out_named[i]}\t{n_out_both[i]}\n")
    print(f"{n_loci} loci, {n} samples: " + ", ".join(f"{k} {v}" for k, v in sorted(status_counts.items())), file=sys.stderr)
    # ---- tr_locus_stats.py (end) ----
    EOF

    bcftools query \
      --regions-file "~{regions_bed}" \
      --regions-overlap 0 \
      --format '%CHROM\t%POS\t%INFO/END\t%INFO/TRID\t%INFO/MOTIFS[\t%AL][\t%SD]\n' \
      "in/~{basename(trgt_vcf)}" \
    | python3 tr_locus_stats.py vcf_samples.txt sexes.tsv "~{named_overlaps}" "~{out_prefix}" \
        ~{outlier_min_bp} ~{outlier_mad_min} ~{outlier_min_called} ~{outlier_max_carrier_fraction} ~{outlier_max_carriers_min}
  >>>

  output {
    File loci = "~{out_prefix}.loci.tsv"
    File spectrum = "~{out_prefix}.spectrum.tsv"
    File outliers = "~{out_prefix}.outliers.tsv"
    File member_counts = "~{out_prefix}.members.tsv"
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

task ugc_wgw_tr_screen {
  meta {
    description: "The named disease loci of the catalog (STRchive's set) for every member: per allele the copies of the pathogenic motif(s) (TRGT's motif counts summed over the catalog motifs the reference table marks pathogenic; allele length over motif length when the table matched none), the class against STRchive's thresholds (normal, intermediate, pathogenic, the gaps by name, uncalled), the pathogenic motif's longest pure segment re-paired with the allele through the genotype, methylation, spanning reads and purity per allele, coded inheritance fit, and the overlapping ordinary locus's lengths. A research screen: a class is a reason to look, not a result."
    outputs: {
      screen: {
        description: "One row per member and named locus (uncalled ones included): sample_id, sex, trid, gene, disease, inheritance, chrom, start, end, gt, ploidy, class, n_alleles_pathogenic, mode_match, measure, borderline, a1_* and a2_* (copies, class, al, allr, lps, sd, ap, am, mc; a1 is the allele with more copies), twin_trid, twin_al, notes"
      },
      screen_loci: {
        description: "One row per named locus: the thresholds, members expected and called, call rate, the copies distribution, members per class (below_pathogenic: short of the pathogenic minimum at a locus without a benign maximum; not_assessable: the reference table says the thresholds do not count this repeat unit), the fractions of alleles above and below the benign range, median methylation per sex group, the twin, call_note (ok, few_uncalled, mostly_uncalled, not_assessable, threshold_or_frame_suspect: more than a fifth of the alleles above or below the benign range)"
      }
    }
  }

  parameter_meta {
    trgt_vcf: {
      description: "The cohort's merged TRGT VCF (cohort_merge)"
    }
    trgt_vcf_index: {
      description: "Index for trgt_vcf"
    }
    trgt_lps: {
      description: "The cohort's trgt-lps table (cohort_merge): trid, motif, one column per sample, longest pure segment per allele in genotype order"
    }
    disease_loci_tsv: {
      description: "The named loci with STRchive's thresholds and the catalog motifs that count as pathogenic (ugc-wgw-tr, scripts/make-strchive-map.py)"
    }
    named_overlaps: {
      description: "ordinary_trid, named_trid, overlap_bp (ugc-wgw-tr)"
    }
    sample_ids: {
      description: "Cohort members; every one must be a sample of the VCF"
    }
    sample_sexes: {
      description: "XX, XY or empty per member, aligned with sample_ids"
    }
    borderline_bp: {
      description: "An allele is borderline when its class changes within this many bp of sequence (or within TRGT's own allele-length range, when wider)"
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
    File trgt_vcf
    File trgt_vcf_index
    File trgt_lps
    File disease_loci_tsv
    File named_overlaps
    Array[String] sample_ids
    Array[String] sample_sexes
    Int borderline_bp = 6
    String out_prefix
    Int threads = 1
    Int mem_gb = 2
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(trgt_vcf, "GB") + size(trgt_lps, "GB") + 20)

  command <<<
    set -euo pipefail
    export LC_ALL=C

    bcftools --version >&2

    mkdir -p in
    ln --symbolic "~{trgt_vcf}" "in/~{basename(trgt_vcf)}"
    ln --symbolic "~{trgt_vcf_index}" "in/~{basename(trgt_vcf_index)}"
    paste "~{write_lines(sample_ids)}" "~{write_lines(sample_sexes)}" > sexes.tsv

    cat > tr_screen.py <<'EOF'
    # ---- tr_screen.py (begin) ----
    import math, re, statistics, sys

    def read_table(path):
        """A TSV with `##` comment lines and a header row: the list of row dicts."""
        rows, cols = [], None
        for line in open(path):
            if line.startswith("##") or not line.strip():
                continue
            c = line.rstrip("\n").split("\t")
            if cols is None:
                cols = c
                continue
            rows.append(dict(zip(cols, c)))
        return rows

    mode = sys.argv[1]
    loci = {r["trid"]: r for r in read_table(sys.argv[2])}
    twin_of = {r["ordinary_trid"]: r["named_trid"] for r in read_table(sys.argv[3])}

    if mode == "regions":
        # the intervals bcftools view -R needs: POS is the base before the repeat, so the BED start moves one to the left
        rows = [(r["chrom"], int(r["start"]), int(r["end"])) for r in loci.values()]
        for t in twin_of:
            chrom, start, end = t.split("_")[:3]
            rows.append((chrom, int(start), int(end)))
        for chrom, start, end in sorted(set(rows)):
            print(f"{chrom}\t{max(0, start - 1)}\t{end}")
        sys.exit(0)

    vcf_path, lps_path, samples_path, sexes_path, prefix, borderline_bp = sys.argv[4:10]
    borderline_bp = int(borderline_bp)
    members = [l.strip() for l in open(samples_path) if l.strip()]
    sex = dict(l.rstrip("\n").split("\t")[:2] for l in open(sexes_path) if "\t" in l)
    SEVERITY = ["uncalled", "not_assessable", "normal", "below_pathogenic", "below_benign", "above_benign", "intermediate", "pathogenic"]

    def num(s):
        return float(s) if s not in ("", None) else None

    def classify(c, r):
        """The class of an allele with c copies at a locus (STRchive's thresholds, in copies): the benign range first, since at
        two loci it sits inside the pathogenic one; then pathogenic, intermediate, and the names of the gaps."""
        if r.get("assessable", "yes") == "no":
            return "not_assessable"   # the thresholds do not count this repeat unit (the table's assessable_note says why)
        bmin, bmax, imin, imax, pmin = (num(r[k]) for k in ("benign_min", "benign_max", "intermediate_min", "intermediate_max", "pathogenic_min"))
        if bmax is not None and (bmin or 0) <= c <= bmax:
            return "normal"
        if pmin is not None and c >= pmin:
            return "pathogenic"
        if imin is not None and imax is not None and imin <= c <= imax:
            return "intermediate"
        if bmin is not None and c < bmin:
            return "below_benign"
        return "above_benign" if bmax is not None else "below_pathogenic"

    def mode_match(inheritance, n_path, ploidy):
        if n_path == 0:
            return "none"
        if inheritance == "AR" or (inheritance == "XR" and ploidy == 2):
            return "consistent" if n_path >= 2 else "carrier"
        return "consistent"

    def fmt(x):
        if x is None:
            return ""
        if isinstance(x, float):
            return str(int(x)) if x == int(x) else f"{x:.4f}".rstrip("0").rstrip(".")
        return str(x)

    # ---- the LPS rows of the named loci: {(trid, motif): {sample: [per GT slot]}}
    lps = {}
    lps_samples = None
    for line in open(lps_path):
        c = line.rstrip("\n").split("\t")
        if lps_samples is None:
            lps_samples = c[2:]
            continue
        if c[0] in loci:
            lps[(c[0], c[1])] = dict(zip(lps_samples, c[2:]))

    ALLELE = ["copies", "class", "al", "allr", "lps", "sd", "ap", "am", "mc"]
    COLS = (["sample_id", "sex", "trid", "gene", "disease", "inheritance", "chrom", "start", "end", "gt", "ploidy", "class", "n_alleles_pathogenic",
             "mode_match", "measure", "borderline"] + [f"a1_{k}" for k in ALLELE] + [f"a2_{k}" for k in ALLELE] + ["twin_trid", "twin_al", "notes"])
    LOCUS_COLS = ["trid", "gene", "disease", "inheritance", "chrom", "start", "end", "motifs", "pathogenic_motifs", "benign_min", "benign_max",
                  "intermediate_min", "intermediate_max", "pathogenic_min", "pathogenic_max", "measure", "n_expected", "n_called", "call_rate", "n_alleles",
                  "copies_median", "copies_p95", "copies_max", "n_normal", "n_intermediate", "n_above_benign", "n_below_benign", "n_below_pathogenic",
                  "n_not_assessable", "n_pathogenic", "n_expanded_other_motif", "frac_alleles_above_benign", "frac_alleles_below_benign", "am_median_xx",
                  "am_median_xy", "twin_trid", "call_note"]
    vcf_samples = None
    records = {}     # trid -> (chrom, pos, end, ref, alts, motifs, fmt keys, sample fields)
    for line in open(vcf_path):
        if line.startswith("##"):
            continue
        c = line.rstrip("\n").split("\t")
        if line.startswith("#CHROM"):
            vcf_samples = c[9:]
            continue
        info = dict(x.split("=", 1) for x in c[7].split(";") if "=" in x)
        trid = info.get("TRID", "")
        if trid in loci or trid in twin_of:
            records[trid] = (c[0], int(c[1]), int(info.get("END", 0)), c[3], c[4].split(","), info.get("MOTIFS", "").split(","), c[8].split(":"), c[9:])
    if vcf_samples is None:
        sys.exit("no #CHROM line in the VCF text")
    missing = [m for m in members if m not in vcf_samples]
    if missing:
        sys.exit("cohort members absent from the TRGT VCF: " + " ".join(missing))
    col = {s: i for i, s in enumerate(vcf_samples)}
    twin_by_named = {}
    for o, nm in twin_of.items():
        twin_by_named.setdefault(nm, o)

    def sample_field(rec, member):
        keys, fields = rec[6], rec[7]
        return dict(zip(keys, fields[col[member]].split(":")))

    n_mismatch = 0
    with open(prefix + ".screen.tsv", "w") as out, open(prefix + ".screen_loci.tsv", "w") as lout:
        out.write("\t".join(COLS) + "\n")
        lout.write("\t".join(LOCUS_COLS) + "\n")
        for trid, r in loci.items():
            rec = records.get(trid)
            twin = twin_by_named.get(trid, "")
            twin_rec = records.get(twin) if twin else None
            chrom = r["chrom"]
            motif_len = int(float(r["motif_len"] or 1)) or 1
            pcm = [m for m in r["pathogenic_catalog_motifs"].split(",") if m]
            rec_motifs = rec[5] if rec else r["catalog_motifs"].split(",")
            changed = rec is not None and rec_motifs != r["catalog_motifs"].split(",")
            use_mc = bool(pcm) and not changed and all(m in rec_motifs for m in pcm)
            measure = "mc_pathogenic" if use_mc else "al_div_motif_len"
            idx = [rec_motifs.index(m) for m in pcm] if use_mc else []
            expected = [m for m in members if not (chrom == "chrY" and sex.get(m, "") == "XX")]
            stats = {"copies": [], "classes": [], "row_classes": [], "am_xx": [], "am_xy": [], "other_motif": 0, "called": 0}
            for member in members:
                sx = sex.get(member, "")
                row = {"sample_id": member, "sex": sx, "trid": trid, "gene": r["gene"], "disease": r["disease"], "inheritance": r["inheritance"],
                       "chrom": chrom, "start": r["start"], "end": r["end"], "measure": measure, "twin_trid": twin}
                notes = []
                if twin_rec is not None:
                    row["twin_al"] = sample_field(twin_rec, member).get("AL", ".")
                d = sample_field(rec, member) if rec else {}
                gt_s = d.get("GT", ".")
                if not rec or gt_s in (".", "./.", "") or d.get("AL", ".") == ".":
                    row.update({"gt": ".", "ploidy": 0, "class": "uncalled", "n_alleles_pathogenic": 0, "mode_match": "none", "borderline": ""})
                    if changed:
                        notes.append("catalog_motifs_changed")
                    row["notes"] = ";".join(notes)
                    out.write("\t".join(fmt(row.get(k, "")) for k in COLS) + "\n")
                    continue
                gt = [int(x) for x in re.split(r"[/|]", gt_s)]
                ploidy = len(gt)
                order = sorted(range(ploidy), key=lambda i: (gt[i], i))   # FORMAT position p holds the allele of GT slot order[p]
                def vals(key):
                    v = d.get(key, ".").split(",")
                    return v if len(v) == ploidy else ["."] * ploidy
                al, allr, sd, mc, ap, am = (vals(k) for k in ("AL", "ALLR", "SD", "MC", "AP", "AM"))
                alleles = []
                for p in range(ploidy):
                    a = {"al": int(al[p]), "allr": allr[p], "sd": sd[p], "ap": ap[p], "am": am[p], "mc": mc[p]}
                    seq = rec[3] if gt[order[p]] == 0 else rec[4][gt[order[p]] - 1]
                    if len(seq) - 1 != a["al"]:
                        n_mismatch += 1
                        notes.append("al_seq_mismatch")
                    counts = [int(x) for x in mc[p].split("_")] if mc[p] != "." and all(x.isdigit() for x in mc[p].split("_")) else []
                    if use_mc and len(counts) == len(rec_motifs):
                        a["copies"] = sum(counts[i] for i in idx)
                    else:
                        a["copies"] = a["al"] // motif_len
                    best = None
                    for m in (pcm if use_mc else []):
                        v = (lps.get((trid, m)) or {}).get(member)
                        if v is None or v == ".":
                            continue
                        slots = v.split(",")
                        if len(slots) == ploidy and slots[order[p]].isdigit():
                            best = max(best or 0, int(slots[order[p]]))
                    if use_mc and best is None and "lps_missing" not in notes:
                        notes.append("lps_missing")
                    a["lps"] = best
                    a["class"] = classify(a["copies"], r)
                    pmin, pmax = num(r["pathogenic_min"]), num(r["pathogenic_max"])
                    # long by length, short by the pathogenic motif: only an allele longer than the reference one can be an
                    # expansion of another motif (TRGT's frame holds other motifs' copies in the reference allele too)
                    if use_mc and a["class"] != "not_assessable" and pmin is not None and a["al"] > rec[2] - rec[1] and a["al"] // motif_len >= pmin and a["copies"] < pmin:
                        notes.append("expanded_other_motif")
                        stats["other_motif"] += 1
                    if a["class"] == "pathogenic" and pmax is not None and a["copies"] > pmax:
                        notes.append("beyond_pathogenic_max")
                    # borderline: the class changes within the measurement's reach, borderline_bp of sequence or TRGT's own
                    # allele-length range, whichever is wider (a 24 bp unit is not mistaken by one copy; a 3 bp unit can be);
                    # the lower edge of the benign range is not a boundary anybody acts on
                    tol = borderline_bp // motif_len
                    lo, hi = a["copies"] - tol, a["copies"] + tol
                    rng = re.fullmatch(r"(\d+)-(\d+)", allr[p])
                    if rng:
                        lo = min(lo, a["copies"] - (a["al"] - int(rng.group(1))) // motif_len)
                        hi = max(hi, a["copies"] + (int(rng.group(2)) - a["al"]) // motif_len)
                    near = {classify(max(0, lo), r), a["class"], classify(hi, r)}
                    a["borderline"] = len({"normal" if c == "below_benign" else c for c in near}) > 1
                    alleles.append(a)
                    stats["copies"].append(a["copies"])
                    stats["classes"].append(a["class"])
                    if am[p] not in (".", "") and sx in ("XX", "XY"):
                        stats["am_xx" if sx == "XX" else "am_xy"].append(float(am[p]))
                alleles.sort(key=lambda a: (-a["copies"], -a["al"]))
                n_path = sum(1 for a in alleles if a["class"] == "pathogenic")
                # a homozygous call of an expanded allele at a recessive locus is two alleles or one with the other dropped out:
                # said only where it matters (it was on every homozygous row, 5,043 of them on an 864-member cohort)
                if r["inheritance"] in ("AR", "XR") and ploidy == 2 and gt[0] == gt[1] and any(a["class"] in ("pathogenic", "intermediate") for a in alleles):
                    notes.append("hom_call")
                if chrom == "chrX" and ((sx == "XY" and ploidy == 2) or (sx == "XX" and ploidy == 1)):
                    notes.append("ploidy_sex_mismatch")
                if changed:
                    notes.append("catalog_motifs_changed")
                row_class = max((a["class"] for a in alleles), key=SEVERITY.index)
                row.update({"gt": gt_s, "ploidy": ploidy, "class": row_class, "n_alleles_pathogenic": n_path,
                            "mode_match": mode_match(r["inheritance"], n_path, ploidy), "borderline": "yes" if any(a["borderline"] for a in alleles) else "no",
                            "notes": ";".join(dict.fromkeys(notes))})
                for name, a in zip(("a1", "a2"), alleles):
                    for k in ALLELE:
                        row[f"{name}_{k}"] = a[k]
                stats["called"] += 1
                stats["row_classes"].append(row_class)
                out.write("\t".join(fmt(row.get(k, "")) for k in COLS) + "\n")
            cp = sorted(stats["copies"])
            n_all = len(cp)
            # outside the benign range, where the locus has one: below_pathogenic (no benign maximum: every allele short of
            # the pathogenic minimum) and not_assessable are neither
            above = sum(1 for c in stats["classes"] if c in ("above_benign", "intermediate", "pathogenic"))
            below = sum(1 for c in stats["classes"] if c == "below_benign")
            frac = above / n_all if n_all else None
            frac_below = below / n_all if n_all else None
            rate = stats["called"] / len(expected) if expected else None
            if r.get("assessable", "yes") == "no":
                note = "not_assessable"
            elif not stats["called"] or (rate is not None and rate < 0.5):
                note = "mostly_uncalled"
            elif (frac is not None and frac > 0.2) or (frac_below is not None and frac_below > 0.2):
                note = "threshold_or_frame_suspect"
            elif stats["called"] < len(expected):
                note = "few_uncalled"
            else:
                note = "ok"
            rc = stats["row_classes"]
            lrow = {"trid": trid, "gene": r["gene"], "disease": r["disease"], "inheritance": r["inheritance"], "chrom": chrom, "start": r["start"],
                    "end": r["end"], "motifs": ",".join(rec_motifs), "pathogenic_motifs": r["pathogenic_catalog_motifs"], "measure": measure,
                    "n_expected": len(expected), "n_called": stats["called"], "call_rate": round(rate, 4) if rate is not None else None, "n_alleles": n_all,
                    "copies_median": statistics.median(cp) if cp else None, "copies_p95": cp[max(0, math.ceil(0.95 * n_all) - 1)] if cp else None,
                    "copies_max": cp[-1] if cp else None, "n_normal": rc.count("normal"), "n_intermediate": rc.count("intermediate"),
                    "n_above_benign": rc.count("above_benign"), "n_below_benign": rc.count("below_benign"),
                    "n_below_pathogenic": rc.count("below_pathogenic"), "n_not_assessable": rc.count("not_assessable"),
                    "n_pathogenic": rc.count("pathogenic"), "n_expanded_other_motif": stats["other_motif"],
                    "frac_alleles_above_benign": round(frac, 4) if frac is not None else None,
                    "frac_alleles_below_benign": round(frac_below, 4) if frac_below is not None else None,
                    "am_median_xx": round(statistics.median(stats["am_xx"]), 3) if stats["am_xx"] else None,
                    "am_median_xy": round(statistics.median(stats["am_xy"]), 3) if stats["am_xy"] else None, "twin_trid": twin, "call_note": note}
            for k in ("benign_min", "benign_max", "intermediate_min", "intermediate_max", "pathogenic_min", "pathogenic_max"):
                lrow[k] = r[k]
            lout.write("\t".join(fmt(lrow.get(k, "")) for k in LOCUS_COLS) + "\n")
    print(f"{len(loci)} named loci, {len(members)} members, {sum(1 for t in loci if t in records)} loci with a record, {n_mismatch} allele-length mismatches",
          file=sys.stderr)
    # ---- tr_screen.py (end) ----
    EOF

    # the named loci and their twins, with the base before each repeat (a record's POS)
    python3 tr_screen.py regions "~{disease_loci_tsv}" "~{named_overlaps}" > named.bed
    bcftools view --no-version --regions-file named.bed "in/~{basename(trgt_vcf)}" > named.vcf
    # the LPS rows of the named loci: every other TRID is chr<N>_<start>_<end>_<motif>
    awk -F'\t' 'NR == 1 || $1 !~ /^chr[^_]+_[0-9]+_[0-9]+_/' "~{trgt_lps}" > named.lps.tsv

    python3 tr_screen.py screen "~{disease_loci_tsv}" "~{named_overlaps}" named.vcf named.lps.tsv \
      "~{write_lines(sample_ids)}" sexes.tsv "~{out_prefix}" ~{borderline_bp}
  >>>

  output {
    File screen = "~{out_prefix}.screen.tsv"
    File screen_loci = "~{out_prefix}.screen_loci.tsv"
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

task ugc_wgw_tr_gather {
  meta {
    description: "The shard tables of ugc_wgw_tr_locus_stats in shard order (the big ones bgzipped), the strongest outliers per member, one row per member (loci called, outliers, the screen's classes), and the long-format summary the driver reads"
    outputs: {
      loci_gz: {
        description: "<out_prefix>.loci.tsv.gz: every locus of the catalog, the columns of the shard table"
      },
      spectrum: {
        description: "<out_prefix>.spectrum.tsv.gz: trid, al, n_alleles, n_members"
      },
      outliers_gz: {
        description: "<out_prefix>.outliers.tsv.gz: every outlier pair"
      },
      outliers_top: {
        description: "<out_prefix>.outliers_top.tsv: per member the top_outliers_per_member pairs by scaled MADs, then excess"
      },
      samples: {
        description: "<out_prefix>.samples.tsv: sample_id, sex, n_loci_called, n_loci_uncalled, n_outliers (ordinary loci), n_outliers_named, n_outliers_both_alleles, n_named_called, n_named_uncalled, n_loci_pathogenic, n_loci_intermediate, n_loci_above_benign, n_loci_below_benign, n_ploidy_sex_mismatch"
      },
      summary: {
        description: "<out_prefix>.summary.tsv: resource, metric, key, value (cohort counts, loci per outlier status, parameters, the reference table's header, the 50 loci with the most carriers, the screen's call notes and non-normal classes per locus)"
      }
    }
  }

  parameter_meta {
    loci: {
      description: "Shard loci tables, shard order"
    }
    spectra: {
      description: "Shard spectrum tables, shard order"
    }
    outliers: {
      description: "Shard outlier tables, shard order"
    }
    member_counts: {
      description: "Shard per-sample counts"
    }
    screen: {
      description: "The screen (ugc_wgw_tr_screen)"
    }
    screen_loci: {
      description: "The screen's per-locus table"
    }
    sample_ids: {
      description: "Cohort members, cohort order"
    }
    sample_sexes: {
      description: "XX, XY or empty per member, aligned with sample_ids"
    }
    disease_loci_tsv: {
      description: "The named-loci table, for its header lines (source and attribution)"
    }
    top_outliers_per_member: {
      description: "Rows per member in outliers_top"
    }
    parameters: {
      description: "key=value pairs joined by ';', written into the summary"
    }
    out_prefix: {
      description: "Output prefix"
    }
    threads: {
      description: "Threads (bgzip)"
    }
    mem_gb: {
      description: "Memory (GB); the outlier pairs are held in memory"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    Array[File] loci
    Array[File] spectra
    Array[File] outliers
    Array[File] member_counts
    File screen
    File screen_loci
    Array[String] sample_ids
    Array[String] sample_sexes
    File disease_loci_tsv
    Int top_outliers_per_member = 20
    String parameters = ""
    String out_prefix
    Int threads = 2
    Int mem_gb = 4
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(3 * (size(loci, "GB") + size(spectra, "GB") + size(outliers, "GB")) + 10)

  command <<<
    set -euo pipefail
    export LC_ALL=C

    paste "~{write_lines(sample_ids)}" "~{write_lines(sample_sexes)}" > sexes.tsv

    python3 - "~{write_lines(loci)}" "~{write_lines(spectra)}" "~{write_lines(outliers)}" "~{write_lines(member_counts)}" \
      "~{screen}" "~{screen_loci}" "~{write_lines(sample_ids)}" sexes.tsv "~{disease_loci_tsv}" "~{out_prefix}" \
      ~{top_outliers_per_member} "~{parameters}" <<'EOF'
    # ---- tr_gather.py (begin) ----
    import sys
    from collections import Counter
    (loci_list, spectrum_list, outliers_list, members_list, screen_path, screen_loci_path, samples_path, sexes_path, disease_path, prefix, top_n,
     params) = sys.argv[1:13]
    top_n = int(top_n)
    members = [l.strip() for l in open(samples_path) if l.strip()]
    sex = dict(l.rstrip("\n").split("\t")[:2] for l in open(sexes_path) if "\t" in l)

    def paths(list_file):
        return [l.strip() for l in open(list_file) if l.strip()]

    def concat(list_file, dest):
        """Shard tables one after the other under the first one's header; returns (columns, data rows written)."""
        cols, n = None, 0
        with open(dest, "w") as out:
            for p in paths(list_file):
                with open(p) as fh:
                    head = fh.readline()
                    if cols is None:
                        cols = head.rstrip("\n").split("\t")
                        out.write(head)
                    elif head.rstrip("\n").split("\t") != cols:
                        sys.exit(f"{p}: header differs from the first shard's")
                    for line in fh:
                        out.write(line)
                        n += 1
        return cols or [], n

    rows_out = []
    def add(resource, metric, key, value):
        rows_out.append((resource, metric, str(key), str(value)))

    # ---- loci: concatenated; the statuses counted on the way
    status = Counter()
    named_listed = []
    n_loci = 0
    loci_cols = None
    with open(prefix + ".loci.tsv", "w") as out:
        for p in paths(loci_list):
            with open(p) as fh:
                head = fh.readline()
                cols = head.rstrip("\n").split("\t")
                if loci_cols is None:
                    loci_cols = cols
                    out.write(head)
                    si, ti = cols.index("outlier_status"), cols.index("trid")
                elif cols != loci_cols:
                    sys.exit(f"{p}: header differs from the first shard's")
                for line in fh:
                    out.write(line)
                    n_loci += 1
                    status[line.rstrip("\n").split("\t")[si]] += 1
    _cols, n_spectrum = concat(spectrum_list, prefix + ".spectrum.tsv")

    # ---- outliers: concatenated; the top ones per member (by z, then excess) and the loci with the most carriers
    pair_cols, n_pairs = concat(outliers_list, prefix + ".outliers.tsv")
    per_member = {m: [] for m in members}
    locus_carriers = {}
    if n_pairs:
        ix = {k: i for i, k in enumerate(pair_cols)}
        with open(prefix + ".outliers.tsv") as fh:
            fh.readline()
            for line in fh:
                c = line.rstrip("\n").split("\t")
                per_member.setdefault(c[ix["sample_id"]], []).append((float(c[ix["z"]]), float(c[ix["excess_bp"]]), c))
                locus_carriers[c[ix["trid"]]] = (int(c[ix["n_carriers"]]), c[ix["chrom"]], c[ix["start"]], c[ix["motifs"]], c[ix["named_locus"]])
    with open(prefix + ".outliers_top.tsv", "w") as out:
        out.write("\t".join(pair_cols) + "\n")
        for m in members:
            for _z, _e, c in sorted(per_member.get(m, []), key=lambda t: (-t[0], -t[1], t[2][1]))[:top_n]:
                out.write("\t".join(c) + "\n")

    # ---- members: the shards' counts summed, the screen's classes per member
    counts = {m: Counter() for m in members}
    for p in paths(members_list):
        with open(p) as fh:
            cols = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                c = dict(zip(cols, line.rstrip("\n").split("\t")))
                if c["sample_id"] in counts:
                    for k in cols[1:]:
                        counts[c["sample_id"]][k] += int(c[k])
    screen_class = Counter()
    with open(screen_path) as fh:
        cols = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            c = dict(zip(cols, line.rstrip("\n").split("\t")))
            m = c["sample_id"]
            if m not in counts:
                continue
            cls = c["class"]
            counts[m]["n_named_uncalled" if cls == "uncalled" else "n_named_called"] += 1
            key = {"pathogenic": "n_loci_pathogenic", "intermediate": "n_loci_intermediate", "above_benign": "n_loci_above_benign",
                   "below_benign": "n_loci_below_benign"}.get(cls)   # below_pathogenic and not_assessable are not findings
            if key:
                counts[m][key] += 1
            if "ploidy_sex_mismatch" in c.get("notes", "").split(";"):
                counts[m]["n_ploidy_sex_mismatch"] += 1
            if cls not in ("uncalled", "normal"):
                screen_class[(cls, c["trid"])] += 1
    SAMPLE_COLS = ["sample_id", "sex", "n_loci_called", "n_loci_uncalled", "n_outliers", "n_outliers_named", "n_outliers_both_alleles", "n_named_called",
                   "n_named_uncalled", "n_loci_pathogenic", "n_loci_intermediate", "n_loci_above_benign", "n_loci_below_benign", "n_ploidy_sex_mismatch"]
    with open(prefix + ".samples.tsv", "w") as out:
        out.write("\t".join(SAMPLE_COLS) + "\n")
        for m in members:
            c = counts[m]
            c["n_loci_called"] = c["n_loci"] - c["n_loci_uncalled"]
            out.write("\t".join([m, sex.get(m, "")] + [str(c[k]) for k in SAMPLE_COLS[2:]]) + "\n")

    # ---- the long summary
    add("cohort", "n_members", "", len(members))
    add("cohort", "n_loci", "", n_loci)
    add("cohort", "n_spectrum_rows", "", n_spectrum)
    add("cohort", "n_outlier_pairs", "", n_pairs)
    for k in ("listed", "common", "none", "below_min_called"):
        add("loci", "outlier_status", k, status.get(k, 0))
    for kv in params.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            add("parameters", k, "", v)
    for line in open(disease_path):
        if not line.startswith("##"):
            break
        k, _, v = line[2:].rstrip("\n").partition("=")
        add("source", k, "", v.replace("\t", " "))
    for trid, (n, chrom, start, motifs, named) in sorted(locus_carriers.items(), key=lambda kv: (-kv[1][0], kv[0]))[:50]:
        add("outlier_loci", "n_carriers", trid, n)
        if named:
            add("outlier_loci", "named_locus", trid, named)
    with open(screen_loci_path) as fh:
        cols = fh.readline().rstrip("\n").split("\t")
        notes = Counter()
        for line in fh:
            c = dict(zip(cols, line.rstrip("\n").split("\t")))
            notes[c["call_note"]] += 1
        add("screen", "n_loci", "", sum(notes.values()))
        for k, v in sorted(notes.items()):
            add("screen", "call_note", k, v)
    for (cls, trid), n in sorted(screen_class.items(), key=lambda kv: (kv[0][0], -kv[1], kv[0][1])):
        add("screen_" + cls, "n_members", trid, n)
    with open(prefix + ".summary.tsv", "w") as out:
        out.write("resource\tmetric\tkey\tvalue\n")
        for r in rows_out:
            out.write("\t".join(r) + "\n")
    print(f"{len(members)} members, {n_loci} loci, {n_pairs} outlier pairs, {n_spectrum} spectrum rows", file=sys.stderr)
    # ---- tr_gather.py (end) ----
    EOF

    for table in loci spectrum outliers; do
      bgzip --force --threads ~{threads} "~{out_prefix}.$table.tsv"
    done
  >>>

  output {
    File loci_gz = "~{out_prefix}.loci.tsv.gz"
    File spectrum = "~{out_prefix}.spectrum.tsv.gz"
    File outliers_gz = "~{out_prefix}.outliers.tsv.gz"
    File outliers_top = "~{out_prefix}.outliers_top.tsv"
    File samples = "~{out_prefix}.samples.tsv"
    File summary = "~{out_prefix}.summary.tsv"
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
