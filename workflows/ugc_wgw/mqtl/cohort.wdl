version 1.0
# ugc-pacbio-wgw — cohort_mqtl's tasks: the island-by-member methylation matrix, the cis tests of structural variants and tandem repeats per region shard, the gather
# origin: none
# see docs/DESIGN.md §5, §18 (2026-11-13); docs/ENTRYPOINTS.md §9.5; docs/expansion/cohort-mqtl.design.md
import "../../../vendor/hifi-human-wgs-wdl/workflows/wdl-common/wdl/structs.wdl"

task ugc_wgw_mqtl_matrix {
  meta {
    description: "The members' annotated scored profiles of cohort_methyl as one table: a row per CpG island (the rows of the cohort-wide label ALL), a column per member with the island's mean methylation (mean_combined_methyl, 0 to 100), empty where the member's median coverage over the island is below min_island_coverage. Beside it the covariates of the model, one row per cohort member: whether it is used (a profile, not excluded), sex, the ancestry components of cohort_qc and the five blood cell fractions of the deconvolution, where those exist."
    outputs: {
      matrix: {
        description: "<out_prefix>.matrix.tsv.gz: chrom, start, end, region_label, gene_symbol, region_type, n_members, mean, sd, then the used members in cohort order"
      },
      covariates: {
        description: "<out_prefix>.covariates.tsv: sample_id, used, why_not, sex, PC1 to PC5, granulocyte, monocyte, T-cell, B-cell, NK-cell"
      }
    }
  }

  parameter_meta {
    profiles: {
      description: "cohort_methyl's scored_profiles_annotated, <sample_id>.<ref_name>.methbat.scored.annotated.tsv per member with a profile"
    }
    ref_name: {
      description: "The reference map's name: part of the profiles' file names, from which the sample ID is read"
    }
    sample_ids: {
      description: "The cohort's members, in order"
    }
    sample_sexes: {
      description: "XX, XY or empty per member"
    }
    deconvolution: {
      description: "cohort_methyl's methyl_deconvolution table, when the stage ran it"
    }
    ancestry: {
      description: "cohort_qc's somalier ancestry table, when that run exists"
    }
    exclude_ids: {
      description: "Members left out of the matrix and the tests"
    }
    min_island_coverage: {
      description: "A member's value for an island is used when its median coverage there is at least this"
    }
    out_prefix: {
      description: "Output prefix"
    }
    mem_gb: {
      description: "Memory (GB): the matrix is held in memory, 32,038 islands by the members"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    Array[File] profiles
    String ref_name
    Array[String] sample_ids
    Array[String] sample_sexes
    File? deconvolution
    File? ancestry
    Array[String] exclude_ids = []
    Int min_island_coverage = 8
    String out_prefix
    Int mem_gb = 8
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(profiles, "GB") * 2 + 20)

  command <<<
    set -euo pipefail
    export LC_ALL=C
    paste "~{write_lines(sample_ids)}" "~{write_lines(sample_sexes)}" > sexes.tsv

    cat > mqtl_matrix.py <<'EOF'
    # ---- mqtl_matrix.py (begin) ----
    import gzip, math, os, sys
    (profiles_path, ref_name, sexes_path, deconv_path, ancestry_path, exclude_path, min_cov, prefix) = sys.argv[1:9]
    min_cov = float(min_cov)
    CELLS = ["granulocyte", "monocyte", "T-cell", "B-cell", "NK-cell"]
    PCS = ["PC1", "PC2", "PC3", "PC4", "PC5"]
    cohort, sex = [], {}
    for l in open(sexes_path):
        c = l.rstrip("\n").split("\t")
        if c and c[0]:
            cohort.append(c[0])
            sex[c[0]] = c[1] if len(c) > 1 else ""
    excluded = {l.strip() for l in open(exclude_path) if l.strip()}
    suffix = "." + ref_name + ".methbat.scored.annotated.tsv"
    profile_of = {}
    for p in (l.strip() for l in open(profiles_path) if l.strip()):
        name = os.path.basename(p)
        if not name.endswith(suffix):
            sys.exit(f"{name}: not a <sample_id>{suffix} file")
        profile_of[name[:-len(suffix)]] = p

    def table(path):
        cols, rows = None, {}
        if not path or not os.path.exists(path):
            return rows
        for l in open(path):
            c = l.rstrip("\n").split("\t")
            if cols is None:
                cols = [x.lstrip("#") for x in c]
                continue
            r = dict(zip(cols, c))
            if r.get("sample_id"):
                rows[r["sample_id"]] = r
        return rows

    deconv, anc = table(deconv_path), table(ancestry_path)
    used, why = [], {}
    for s in cohort:
        if s in excluded:
            why[s] = "excluded (exclude_ids)"
        elif s not in profile_of:
            why[s] = "no scored profile in cohort_methyl"
        else:
            used.append(s)

    def num(r, k):
        try:
            v = float(r.get(k, ""))
        except (TypeError, ValueError):
            return ""
        return f"{v:.4f}" if math.isfinite(v) else ""

    with open(prefix + ".covariates.tsv", "w") as out:
        out.write("\t".join(["sample_id", "used", "why_not", "sex"] + PCS + CELLS) + "\n")
        for s in cohort:
            a, d = anc.get(s, {}), deconv.get(s, {})
            out.write("\t".join([s, "yes" if s in used else "no", why.get(s, ""), sex.get(s, "")]
                                + [num(a, k) for k in PCS] + [num(d, k) for k in CELLS]) + "\n")

    keys, meta, values = [], {}, {}
    for i, s in enumerate(used):
        cols = None
        col = [""] * len(keys)
        with open(profile_of[s]) as fh:
            for l in fh:
                if l.startswith("##"):
                    continue
                c = l.rstrip("\n").split("\t")
                if cols is None:
                    cols = {k: j for j, k in enumerate(c)}
                    need = ["chrom", "start", "end", "region_label", "background_category", "mean_combined_methyl", "median_total_coverage"]
                    missing = [k for k in need if k not in cols]
                    if missing:
                        sys.exit(f"{profile_of[s]}: no column {', '.join(missing)}")
                    continue
                if c[cols["background_category"]] != "ALL":
                    continue
                key = (c[cols["chrom"]], c[cols["start"]], c[cols["end"]])
                if key not in meta:
                    if i > 0:
                        continue          # the islands are the first profile's: every member is scored on the same background
                    meta[key] = (c[cols["region_label"]], c[cols["gene_symbol"]] if "gene_symbol" in cols else "",
                                 c[cols["region_type"]] if "region_type" in cols else "")
                    keys.append(key)
                    values[key] = [""] * len(used)
                v, cov = c[cols["mean_combined_methyl"]], c[cols["median_total_coverage"]]
                try:
                    if v != "" and float(cov or 0) >= min_cov:
                        values[key][i] = f"{float(v):.1f}"
                except ValueError:
                    pass
    with gzip.open(prefix + ".matrix.tsv.gz", "wt", compresslevel=4) as out:
        out.write("\t".join(["chrom", "start", "end", "region_label", "gene_symbol", "region_type", "n_members", "mean", "sd"] + used) + "\n")
        for key in keys:
            vals = [float(v) for v in values[key] if v != ""]
            n = len(vals)
            mean = sum(vals) / n if n else float("nan")
            sd = math.sqrt(sum((v - mean) ** 2 for v in vals) / (n - 1)) if n > 1 else float("nan")
            out.write("\t".join(list(key) + list(meta[key]) + [str(n), f"{mean:.2f}" if n else "", f"{sd:.2f}" if n > 1 else ""] + values[key]) + "\n")
    print(f"{len(keys)} islands, {len(used)} of {len(cohort)} members", file=sys.stderr)
    # ---- mqtl_matrix.py (end) ----
    EOF

    python3 mqtl_matrix.py "~{write_lines(profiles)}" "~{ref_name}" sexes.tsv "~{default="" deconvolution}" "~{default="" ancestry}" \
      "~{write_lines(exclude_ids)}" ~{min_island_coverage} "~{out_prefix}"
  >>>

  output {
    File matrix = "~{out_prefix}.matrix.tsv.gz"
    File covariates = "~{out_prefix}.covariates.tsv"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:03cb3c01937eccc907f8ad71c87b258581504572205fe3f31a657e318f3564ae"  # pb_wdl_base:build4
    cpu: 1
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

task ugc_wgw_mqtl_shard {
  meta {
    description: "One region shard of the screen: every autosomal island of the shard with enough covered members and some spread is tested against every structural variant (carrier dosage 0, 1, 2; svx's ./. is a non-carrier) and every tandem repeat (the mean of the member's two allele lengths) within window_bp of it. Methylation is rank-normalised across the members; both sides are residualised on the covariates (an intercept, sex, the ancestry components and the blood cell fractions when every used member has them); the test is the correlation of the residuals with a two-sided p-value from Student's t. The effect is also given on the methylation scale: points per allele, or per 100 bp of repeat length. A variant is tested when at least min_carriers members differ from the commonest genotype and at least as many have it. Standard-library Python."
    outputs: {
      pairs: {
        description: "<out_prefix>.pairs.tsv: one row per tested pair (no q-value yet: the gather corrects over the whole run)"
      },
      counts: {
        description: "<out_prefix>.counts.tsv: metric, key, value: islands seen, tested and left out by reason, variants tested by type, the model's covariates, members"
      }
    }
  }

  parameter_meta {
    matrix: {
      description: "The island-by-member matrix (ugc_wgw_mqtl_matrix)"
    }
    covariates: {
      description: "The covariates table (ugc_wgw_mqtl_matrix)"
    }
    sv_vcf: {
      description: "cohort_merge's SV VCF"
    }
    sv_vcf_index: {
      description: "Its index"
    }
    trgt_vcf: {
      description: "cohort_merge's merged TRGT VCF"
    }
    trgt_vcf_index: {
      description: "Its index"
    }
    regions_bed: {
      description: "The shard: islands starting inside it are this shard's"
    }
    window_bp: {
      description: "A variant and an island are a pair within this distance"
    }
    min_carriers: {
      description: "Members that must differ from the commonest genotype, and members that must have it"
    }
    min_members: {
      description: "Below this many used members nothing is tested"
    }
    min_island_members_fraction: {
      description: "An island is tested when this share of the used members has a value"
    }
    min_island_sd: {
      description: "An island is tested when its methylation's standard deviation over the members is at least this (points)"
    }
    test_svs: {
      description: "Test structural variants"
    }
    test_trs: {
      description: "Test tandem repeats"
    }
    out_prefix: {
      description: "Output prefix"
    }
    mem_gb: {
      description: "Memory (GB)"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    File matrix
    File covariates
    File sv_vcf
    File sv_vcf_index
    File trgt_vcf
    File trgt_vcf_index
    File regions_bed
    Int window_bp = 100000
    Int min_carriers = 10
    Int min_members = 40
    Float min_island_members_fraction = 0.8
    Float min_island_sd = 2.0
    Boolean test_svs = true
    Boolean test_trs = true
    String out_prefix
    Int mem_gb = 4
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil((size(sv_vcf, "GB") + size(trgt_vcf, "GB")) * 2 + 20)

  command <<<
    set -euo pipefail
    export LC_ALL=C
    bcftools --version | head -1 >&2
    mkdir -p in
    ln --symbolic "~{sv_vcf}" "in/~{basename(sv_vcf)}"
    ln --symbolic "~{sv_vcf_index}" "in/~{basename(sv_vcf_index)}"
    ln --symbolic "~{trgt_vcf}" "in/~{basename(trgt_vcf)}"
    ln --symbolic "~{trgt_vcf_index}" "in/~{basename(trgt_vcf_index)}"

    # the shard widened by the window: the variants an island at the shard's edge is tested against
    awk -F'\t' -v OFS='\t' -v w=~{window_bp} '!/^#/ && NF >= 3 { s = $2 - w; if (s < 0) s = 0; print $1, s, $3 + w }' "~{regions_bed}" > window.bed

    : > sv.txt; : > tr.txt
    bcftools query --list-samples "in/~{basename(sv_vcf)}" > sv_samples.txt
    bcftools query --list-samples "in/~{basename(trgt_vcf)}" > tr_samples.txt
    if [ "~{test_svs}" = true ] && [ -s window.bed ]; then
      bcftools query --regions-file window.bed --regions-overlap 1 \
        --include 'INFO/SVTYPE="INS" || INFO/SVTYPE="DEL" || INFO/SVTYPE="DUP" || INFO/SVTYPE="INV"' \
        --format '%CHROM\t%POS\t%INFO/END\t%ID\t%INFO/SVTYPE\t%INFO/SVLEN[\t%GT]\n' "in/~{basename(sv_vcf)}" | sort -u > sv.txt
    fi
    if [ "~{test_trs}" = true ] && [ -s window.bed ]; then
      bcftools query --regions-file window.bed --regions-overlap 1 \
        --format '%CHROM\t%POS\t%INFO/END\t%INFO/TRID\t%INFO/MOTIFS\t.[\t%AL]\n' "in/~{basename(trgt_vcf)}" | sort -u > tr.txt
    fi

    cat > mqtl_shard.py <<'EOF'
    # ---- mqtl_shard.py (begin) ----
    import gzip, math, re, sys
    from array import array
    from statistics import NormalDist
    (matrix_path, cov_path, bed_path, sv_path, sv_samples_path, tr_path, tr_samples_path, prefix,
     window, min_carriers, min_members, min_frac, min_sd) = sys.argv[1:14]
    window, min_carriers, min_members, min_frac, min_sd = int(window), int(min_carriers), int(min_members), float(min_frac), float(min_sd)
    PAIR_COLS = ["chrom", "start", "end", "region_label", "gene_symbol", "region_type", "variant_type", "variant_id", "variant_start",
                 "variant_end", "variant_len", "distance", "n", "n_minor", "r", "t", "df", "p", "effect", "effect_unit"]
    AUTOSOME = re.compile(r"^(chr)?[0-9]+$")
    NON_CARRIER = frozenset({"./.", ".|.", ".", "0/0", "0|0", "0"})
    counts = {}

    def bump(metric, key, n=1):
        counts[(metric, key)] = counts.get((metric, key), 0) + n

    # ---- mqtl_stats (begin) ----
    # statistics without third-party code
    def betacf(a, b, x):
        """The continued fraction of the incomplete beta function (modified Lentz)."""
        tiny = 1e-300
        qab, qap, qam = a + b, a + 1.0, a - 1.0
        c, d = 1.0, 1.0 - qab * x / qap
        d = 1.0 / (d if abs(d) > tiny else tiny)
        h = d
        for m in range(1, 400):
            m2 = 2 * m
            aa = m * (b - m) * x / ((qam + m2) * (a + m2))
            d = 1.0 + aa * d
            d = 1.0 / (d if abs(d) > tiny else tiny)
            c = 1.0 + aa / c
            c = c if abs(c) > tiny else tiny
            h *= d * c
            aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
            d = 1.0 + aa * d
            d = 1.0 / (d if abs(d) > tiny else tiny)
            c = 1.0 + aa / c
            c = c if abs(c) > tiny else tiny
            delta = d * c
            h *= delta
            if abs(delta - 1.0) < 3e-16:
                break
        return h

    def betainc(a, b, x):
        """The regularised incomplete beta function I_x(a, b)."""
        if x <= 0.0:
            return 0.0
        if x >= 1.0:
            return 1.0
        ln = math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b) + a * math.log(x) + b * math.log(1.0 - x)
        if x < (a + 1.0) / (a + b + 2.0):
            return math.exp(ln) * betacf(a, b, x) / a
        return 1.0 - math.exp(ln) * betacf(b, a, 1.0 - x) / b

    def t_two_sided(t, df):
        """P(|T| >= |t|) for Student's t with df degrees of freedom."""
        if df <= 0 or not math.isfinite(t):
            return float("nan")
        return betainc(df / 2.0, 0.5, df / (df + t * t))

    def rank_normal(values):
        """Rank-based inverse normal transform (ties get their mean rank)."""
        n = len(values)
        order = sorted(range(n), key=lambda i: values[i])
        ranks = [0.0] * n
        i = 0
        while i < n:
            j = i
            while j + 1 < n and values[order[j + 1]] == values[order[i]]:
                j += 1
            r = (i + j) / 2.0 + 1.0
            for k in range(i, j + 1):
                ranks[order[k]] = r
            i = j + 1
        nd = NormalDist()
        return [nd.inv_cdf((r - 0.5) / n) for r in ranks]

    def basis(columns):
        """An orthonormal basis of the columns' span (modified Gram-Schmidt); a column that adds nothing is dropped."""
        q = []
        for col in columns:
            v = list(col)
            scale = math.sqrt(sum(x * x for x in v)) or 1.0
            for u in q:
                d = sum(a * b for a, b in zip(u, v))
                v = [a - d * b for a, b in zip(v, u)]
            norm = math.sqrt(sum(x * x for x in v))
            if norm > 1e-8 * scale:
                q.append([x / norm for x in v])
        return q

    def residual(v, q):
        for u in q:
            d = sum(a * b for a, b in zip(u, v))
            v = [a - d * b for a, b in zip(v, u)]
        return v

    def prepare(y_raw, q):
        """An island's side of every pair, once: the residuals of the rank-normalised and of the raw methylation."""
        ry = residual(rank_normal(y_raw), q)
        return ry, residual(list(y_raw), q), sum(x * x for x in ry)

    def test_residuals(ry, rraw, syy, rg, sgg, n, k):
        """r, t, df, p and the slope of raw methylation on the variant, from both sides' residuals on k covariate columns."""
        df = n - k - 1
        if sgg <= 1e-12 or syy <= 1e-12 or df < 1:
            return None
        r = sum(a * b for a, b in zip(ry, rg)) / math.sqrt(sgg * syy)
        r = max(-0.999999999, min(0.999999999, r))
        t = r * math.sqrt(df / (1.0 - r * r))
        return r, t, df, t_two_sided(t, df), sum(a * b for a, b in zip(rraw, rg)) / sgg

    def test(y_raw, g, q):
        """One pair from the members' values, given as parallel lists."""
        ry, rraw, syy = prepare(y_raw, q)
        rg = residual(list(g), q)
        return test_residuals(ry, rraw, syy, rg, sum(x * x for x in rg), len(y_raw), len(q))

    # ---- mqtl_stats (end) ----

    # ---- members and covariates ----
    cov_rows, cov_cols = [], None
    for l in open(cov_path):
        c = l.rstrip("\n").split("\t")
        if cov_cols is None:
            cov_cols = c
            continue
        cov_rows.append(dict(zip(cov_cols, c)))
    used = [r for r in cov_rows if r.get("used") == "yes"]
    members = [r["sample_id"] for r in used]
    n_used = len(members)
    bump("members", "used", n_used)
    bump("members", "cohort", len(cov_rows))

    def complete(keys):
        return n_used > 0 and all(r.get(k, "") != "" for r in used for k in keys)

    columns, names = [[1.0] * n_used], ["intercept"]
    if n_used and all(r.get("sex") in ("XX", "XY") for r in used) and len({r["sex"] for r in used}) == 2:
        columns.append([1.0 if r["sex"] == "XY" else 0.0 for r in used])
        names.append("sex")
    pcs = ["PC1", "PC2", "PC3", "PC4", "PC5"]
    if complete(pcs):
        columns += [[float(r[k]) for r in used] for k in pcs]
        names += pcs
    cells = ["granulocyte", "monocyte", "T-cell", "B-cell", "NK-cell"]
    if complete(cells):
        columns += [[float(r[k]) for r in used] for k in cells]
        names += cells
    for nm in names:
        bump("covariate", nm)

    def write_out(rows):
        with open(prefix + ".pairs.tsv", "w") as out:
            out.write("\t".join(PAIR_COLS) + "\n")
            for r in rows:
                out.write("\t".join(r) + "\n")
        with open(prefix + ".counts.tsv", "w") as out:
            out.write("metric\tkey\tvalue\n")
            for (m, k), v in sorted(counts.items()):
                out.write(f"{m}\t{k}\t{v}\n")

    if n_used < min_members:
        bump("note", f"fewer than {min_members} members with a profile: nothing tested")
        write_out([])
        sys.exit(0)

    # ---- the shard's islands ----
    shard = {}
    for l in open(bed_path):
        c = l.rstrip("\n").split("\t")
        if len(c) >= 3 and not l.startswith("#"):
            shard.setdefault(c[0], []).append((int(c[1]), int(c[2])))
    islands = []
    with gzip.open(matrix_path, "rt") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        col_of = {s: i for i, s in enumerate(head)}
        idx = [col_of[s] for s in members]
        for l in fh:
            c = l.rstrip("\n").split("\t")
            chrom, start = c[0], int(c[1])
            if not any(a <= start < b for a, b in shard.get(chrom, ())):
                continue
            bump("islands", "in the shard")
            if not AUTOSOME.match(chrom):
                bump("islands", "not tested: not an autosome")
                continue
            vals = [c[i] if i < len(c) else "" for i in idx]
            present = [i for i, v in enumerate(vals) if v != ""]
            if len(present) < min_frac * n_used or len(present) < min_members:
                bump("islands", "not tested: too few members with coverage")
                continue
            y = [float(vals[i]) for i in present]
            mean = sum(y) / len(y)
            if math.sqrt(sum((v - mean) ** 2 for v in y) / (len(y) - 1)) < min_sd:
                bump("islands", "not tested: little spread")
                continue
            islands.append((chrom, start, int(c[2]), c[3], c[4], c[5], present, y))
    islands.sort(key=lambda r: (r[0], r[1]))

    # ---- the variants, in the members' order ----
    def reorder(samples_path):
        names_in = [l.strip() for l in open(samples_path) if l.strip()]
        pos = {s: i for i, s in enumerate(names_in)}
        missing = [s for s in members if s not in pos]
        if missing:
            sys.exit(f"members absent from {samples_path}: {' '.join(missing[:10])}")
        return [pos[s] for s in members]

    variants = {}   # chrom -> [(start, end, type, id, length, values or None per member)]
    if sum(1 for _ in open(sv_path)):
        order = reorder(sv_samples_path)
        for l in open(sv_path):
            c = l.rstrip("\n").split("\t")
            if not AUTOSOME.match(c[0]):
                continue
            gts = c[6:]
            g = bytearray(n_used)          # dosage 0, 1, 2 per member: a byte each
            for k, i in enumerate(order):
                a = gts[i]
                if a not in NON_CARRIER:
                    g[k] = min(2, sum(1 for x in re.split(r"[/|]", a) if x not in ("0", ".", "")))
            try:
                end = int(c[2])
            except ValueError:
                end = int(c[1])
            try:
                length = abs(int(c[5].split(",")[0]))
            except ValueError:
                length = max(0, end - int(c[1]))
            variants.setdefault(c[0], []).append((int(c[1]), max(end, int(c[1])), c[4], c[3], length, g))
    if sum(1 for _ in open(tr_path)):
        order = reorder(tr_samples_path)
        for l in open(tr_path):
            c = l.rstrip("\n").split("\t")
            if not AUTOSOME.match(c[0]):
                continue
            als = c[6:]
            g = array("f", [math.nan]) * n_used   # mean allele length per member; NaN without a call
            for k, i in enumerate(order):
                parts = [x for x in als[i].split(",") if x not in (".", "")]
                try:
                    if parts:
                        g[k] = sum(float(x) for x in parts) / len(parts)
                except ValueError:
                    pass
            try:
                end = int(c[2])
            except ValueError:
                end = int(c[1])
            variants.setdefault(c[0], []).append((int(c[1]), max(end, int(c[1])), "TR", c[3], max(0, end - int(c[1])), g))
    for chrom in variants:
        variants[chrom].sort(key=lambda v: v[0])
        bump("variants", "read", len(variants[chrom]))

    def usable(vtype, g):
        """The genotype vector for the model, or None with the reason: too few members off (or on) the commonest
        genotype, or too many without a call. A repeat's missing calls take the members' median length."""
        if vtype == "TR":
            called = sorted(x for x in g if x == x)
            if len(called) < 0.9 * len(g):
                return None, "too many members without a call"
            med = called[len(called) // 2]
            g = [x if x == x else med for x in g]
            keyed = [round(x) for x in g]
        else:
            g = [float(x) for x in g]
            keyed = g
        tally = {}
        for x in keyed:
            tally[x] = tally.get(x, 0) + 1
        top = max(tally.values())
        if len(g) - top < min_carriers or top < min_carriers:
            return None, "too few members off the commonest genotype"
        return g, ""

    # ---- the pairs ----
    # islands and variants are in genome order: a window of variants moves along with the islands, and what was computed
    # for a variant (its usable vector, its residual on the covariates of all members) is dropped when the window leaves it
    q_full = basis(columns)
    q_cache = {}
    rows = []
    tested_variants = set()
    full = tuple(range(n_used))
    cur_chrom, vlist, nxt, active, g_cache = None, [], 0, [], {}
    for chrom, start, end, label, gene, rtype, present, y in islands:
        if chrom != cur_chrom:
            cur_chrom, vlist, nxt, active, g_cache = chrom, variants.get(chrom, []), 0, [], {}
        while nxt < len(vlist) and vlist[nxt][0] <= end + window:
            active.append(vlist[nxt])
            nxt += 1
        if any(v[1] < start - window for v in active):
            for v in active:
                if v[1] < start - window:
                    g_cache.pop((v[2], v[3]), None)
            active = [v for v in active if v[1] >= start - window]
        key = tuple(present)
        is_full = key == full
        if is_full:
            q = q_full
        else:
            q = q_cache.get(key)
            if q is None:
                q = basis([[col[i] for i in present] for col in columns])
                if len(q_cache) < 2000:
                    q_cache[key] = q
        ry, rraw, syy = prepare(y, q)
        island_tested = False
        for vstart, vend, vtype, vid, vlen, g in active:
            if vstart > end + window:
                continue
            hit = g_cache.get((vtype, vid)) if is_full else None
            if hit is None:
                gv, why = usable(vtype, g if is_full else [g[i] for i in present])
                if gv is None:
                    hit = (None, why, None, 0.0, 0)
                else:
                    rg = array("d", residual(gv, q))
                    if vtype == "TR":
                        tally = {}
                        for x in gv:
                            tally[round(x)] = tally.get(round(x), 0) + 1
                        minor = len(gv) - max(tally.values())
                    else:
                        minor = sum(1 for x in gv if x > 0)
                    hit = (True, "", rg, sum(x * x for x in rg), minor)
                if is_full:
                    g_cache[(vtype, vid)] = hit
            ok, why, rg, sgg, minor = hit
            if ok is None:
                bump("pairs", "not tested: " + why)
                continue
            res = test_residuals(ry, rraw, syy, rg, sgg, len(y), len(q))
            if res is None:
                bump("pairs", "not tested: no variation left after the covariates")
                continue
            r, t, df, p, slope = res
            dist = 0 if (vstart <= end and vend >= start) else (start - vend if vend < start else vstart - end)
            if vtype == "TR":
                effect, unit = slope * 100.0, "points per 100 bp"
            else:
                effect, unit = slope, "points per allele"
            rows.append([chrom, str(start), str(end), label, gene, rtype, vtype, vid, str(vstart), str(vend), str(vlen), str(dist), str(len(y)),
                         str(minor), f"{r:.4f}", f"{t:.3f}", str(df), f"{p:.3e}", f"{effect:.3f}", unit])
            tested_variants.add((vtype, vid))
            bump("pairs", "tested")
            bump("pairs_by_type", vtype)
            island_tested = True
        bump("islands", "tested" if island_tested else "not tested: no usable variant in the window")
    for vtype, _vid in tested_variants:
        bump("variants_tested", vtype)
    write_out(rows)
    print(f"{len(islands)} islands, {len(rows)} pairs, {n_used} members, covariates: {', '.join(names)}", file=sys.stderr)
    # ---- mqtl_shard.py (end) ----
    EOF

    python3 mqtl_shard.py "~{matrix}" "~{covariates}" "~{regions_bed}" sv.txt sv_samples.txt tr.txt tr_samples.txt "~{out_prefix}" \
      ~{window_bp} ~{min_carriers} ~{min_members} ~{min_island_members_fraction} ~{min_island_sd}
  >>>

  output {
    File pairs = "~{out_prefix}.pairs.tsv"
    File counts = "~{out_prefix}.counts.tsv"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:03cb3c01937eccc907f8ad71c87b258581504572205fe3f31a657e318f3564ae"  # pb_wdl_base:build4
    cpu: 1
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

task ugc_wgw_mqtl_gather {
  meta {
    description: "The shards' pairs as one table with Benjamini-Hochberg q-values over every test of the run, the pairs at or below the false-discovery rate as the top table, and the summary: members, islands and variants tested and why the others were not, pairs and hits by variant type and region type, the covariates in the model, and the inflation factor (the median of the squared t statistics over 0.455, the median of a chi-square with one degree of freedom)."
    outputs: {
      pairs_gz: {
        description: "<out_prefix>.pairs.tsv.gz: every tested pair with its q-value, in genome order"
      },
      top: {
        description: "<out_prefix>.top.tsv: the pairs with q at or below fdr, by p-value"
      },
      summary: {
        description: "<out_prefix>.summary.tsv: resource, metric, key, value"
      }
    }
  }

  parameter_meta {
    pairs: {
      description: "The shards' pairs tables"
    }
    counts: {
      description: "The shards' counts tables"
    }
    fdr: {
      description: "The false-discovery rate of the top table"
    }
    parameters: {
      description: "The run's parameters, key=value separated by semicolons, written into the summary"
    }
    out_prefix: {
      description: "Output prefix"
    }
    mem_gb: {
      description: "Memory (GB): every pair's p-value is held for the correction"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    Array[File] pairs
    Array[File] counts
    Float fdr = 0.05
    String parameters
    String out_prefix
    Int mem_gb = 4
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(pairs, "GB") * 3 + 20)

  command <<<
    set -euo pipefail
    export LC_ALL=C

    cat > mqtl_gather.py <<'EOF'
    # ---- mqtl_gather.py (begin) ----
    import gzip, sys
    pairs_list, counts_list, fdr, parameters, prefix = sys.argv[1:6]
    fdr = float(fdr)
    head, rows = None, []
    for path in (l.strip() for l in open(pairs_list) if l.strip()):
        with open(path) as fh:
            h = fh.readline().rstrip("\n").split("\t")
            if head is None:
                head = h
            elif h != head:
                sys.exit(f"{path}: columns differ from the first shard's")
            for l in fh:
                rows.append(l.rstrip("\n").split("\t"))
    head = head or []
    col = {k: i for i, k in enumerate(head)}

    def chrom_key(c):
        d = c[3:] if c.startswith("chr") else c
        return (0, int(d)) if d.isdigit() else (1, d)

    if rows:
        rows.sort(key=lambda r: (chrom_key(r[col["chrom"]]), int(r[col["start"]]), int(r[col["variant_start"]]), r[col["variant_id"]]))
    # Benjamini-Hochberg: q(i) = min over j >= i of p(j) * m / j, in order of p
    m = len(rows)
    order = sorted(range(m), key=lambda i: float(rows[i][col["p"]])) if rows else []
    q = [1.0] * m
    best = 1.0
    for rank in range(m, 0, -1):
        i = order[rank - 1]
        best = min(best, float(rows[i][col["p"]]) * m / rank)
        q[i] = best
    with gzip.open(prefix + ".pairs.tsv.gz", "wt", compresslevel=4) as out:
        out.write("\t".join(head + ["q"]) + "\n")
        for i, r in enumerate(rows):
            out.write("\t".join(r + [f"{q[i]:.3e}"]) + "\n")
    hits = [i for i in order if q[i] <= fdr]
    with open(prefix + ".top.tsv", "w") as out:
        out.write("\t".join(head + ["q"]) + "\n")
        for i in hits:
            out.write("\t".join(rows[i] + [f"{q[i]:.3e}"]) + "\n")

    counts = {}
    for path in (l.strip() for l in open(counts_list) if l.strip()):
        with open(path) as fh:
            fh.readline()
            for l in fh:
                c = l.rstrip("\n").split("\t")
                if len(c) == 3:
                    if c[0] in ("members", "covariate"):          # the same in every shard
                        counts[(c[0], c[1])] = int(c[2])
                    elif c[0] == "note":
                        counts[(c[0], c[1])] = 1
                    else:
                        counts[(c[0], c[1])] = counts.get((c[0], c[1]), 0) + int(c[2])
    out_rows = [("mqtl", m_, k, v) for (m_, k), v in sorted(counts.items())]
    out_rows.append(("mqtl", "pairs", "total", m))
    out_rows.append(("mqtl", "hits", "total", len(hits)))
    by_type, by_region, islands_hit = {}, {}, set()
    for i in hits:
        r = rows[i]
        by_type[r[col["variant_type"]]] = by_type.get(r[col["variant_type"]], 0) + 1
        by_region[r[col["region_type"]] or "unannotated"] = by_region.get(r[col["region_type"]] or "unannotated", 0) + 1
        islands_hit.add((r[col["chrom"]], r[col["start"]]))
    out_rows += [("mqtl", "hits_by_variant_type", k, v) for k, v in sorted(by_type.items())]
    out_rows += [("mqtl", "hits_by_region_type", k, v) for k, v in sorted(by_region.items())]
    out_rows.append(("mqtl", "hits", "islands", len(islands_hit)))
    if m:
        chi = sorted(float(r[col["t"]]) ** 2 for r in rows)
        median = chi[m // 2] if m % 2 else (chi[m // 2 - 1] + chi[m // 2]) / 2.0
        out_rows.append(("mqtl", "inflation", "lambda", f"{median / 0.4549364:.3f}"))
    out_rows.append(("mqtl", "parameters", "fdr", f"{fdr:g}"))
    for kv in parameters.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out_rows.append(("mqtl", "parameters", k, v))
    with open(prefix + ".summary.tsv", "w") as out:
        out.write("resource\tmetric\tkey\tvalue\n")
        for r in out_rows:
            out.write("\t".join(str(x) for x in r) + "\n")
    print(f"{m} pairs, {len(hits)} at q <= {fdr:g}", file=sys.stderr)
    # ---- mqtl_gather.py (end) ----
    EOF

    python3 mqtl_gather.py "~{write_lines(pairs)}" "~{write_lines(counts)}" ~{fdr} "~{parameters}" "~{out_prefix}"
  >>>

  output {
    File pairs_gz = "~{out_prefix}.pairs.tsv.gz"
    File top = "~{out_prefix}.top.tsv"
    File summary = "~{out_prefix}.summary.tsv"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:03cb3c01937eccc907f8ad71c87b258581504572205fe3f31a657e318f3564ae"  # pb_wdl_base:build4
    cpu: 1
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
