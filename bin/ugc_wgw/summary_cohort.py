"""The cohort page of `ugc-wgw summary`: rendered from the members' per-sample digests (never from their files) plus
the cohort stages' own small outputs (merged SV VCF header and INFO, trgt-lps dimensions, cohort_freq summary).
See docs/guide/07-results.md ("Analysis summary")."""
from __future__ import annotations

import statistics
from pathlib import Path

from . import export
from . import summary_parsers as sp
from . import summary_svg as sv
from .db import Cohort, RunRecord
from .inputs import SEX_TO_GROUP
from .layout import run_files
from .outputs import Outputs
from .report import Raw, _bytes, _esc, _table, big_table
from .summary import DIGEST_SCHEMA, _cards, _has_read_data, _kv, _nav, _page, _sec, _software_sections, _status_badge
from .util import UgcError, read_json, utc_now

# (key, label, digits, unit, path into the sample digest, threshold key or None)
# (key, label, digits, unit, digest path, threshold key, direction): direction "min" means values below the threshold fail,
# "max" values above it; the assembly metrics read the flat `assembly.summary` block of the sample digest
MEMBER_METRICS: list[tuple[str, str, int, str, tuple[str, ...], str | None, str]] = [
    ("depth_mean", "mean depth", 1, "x", ("coverage", "mean"), "depth_mean_min", "min"),
    ("mapped_pct", "mapped reads %", 2, "%", ("reads", "mapped_read_percent"), "mapped_read_percent_min", "min"),
    ("read_q_median", "read quality (median)", 1, "", ("reads", "read_quality_median"), "read_quality_median_min", "min"),
    ("read_n50", "read N50 (bp)", 0, "", ("reads", "read_length_n50"), None, "min"),
    ("snv", "SNVs", 0, "", ("small_variants", "snv"), None, "min"),
    ("indel", "indels", 0, "", ("small_variants", "indel"), None, "min"),
    ("tstv", "Ts/Tv", 2, "", ("small_variants", "tstv"), None, "min"),
    ("hethom", "het/hom", 2, "", ("small_variants", "hethom"), None, "min"),
    ("sv_total", "SVs", 0, "", ("sv", "_total"), None, "min"),
    ("phased_bp", "phased bases", 0, "", ("phasing", "phased_bp"), None, "min"),
    ("block_ng50", "phase block NG50", 0, "", ("phasing", "block_ng50"), None, "min"),
    ("trgt_genotyped", "TRGT loci genotyped", 0, "", ("tandem_repeats", "genotyped"), None, "min"),
    ("cpg", "CpG sites", 0, "", ("methylation", "cpg", "combined"), None, "min"),
    ("roh_total", "ROH total (bp)", 0, "", ("small_variants", "roh", "total_bp"), None, "min"),
    ("chrx_ratio", "chrX / autosomes", 2, "", ("coverage", "ratios", "chrX_ratio"), None, "min"),
    ("chry_ratio", "chrY / autosomes", 3, "", ("coverage", "ratios", "chrY_ratio"), None, "min"),
    ("mito_unexplained", "mitorsaw unexplained fraction", 3, "", ("targeted", "mitorsaw", "fingerprint", "unexplained_fraction"), "mito_unexplained_max", "max"),
    ("mt_heteroplasmic", "chrM heteroplasmic sites", 0, "", ("targeted", "mitorsaw", "collapsed", "heteroplasmic"), None, "max"),
    ("asm_hap_length_min", "assembly: shorter haplotype (bp)", 0, "", ("assembly", "summary", "hap_length_min"), "asm_hap_length_min", "min"),
    ("asm_hap_length_max", "assembly: longer haplotype (bp)", 0, "", ("assembly", "summary", "hap_length_max"), "asm_hap_length_max", "max"),
    ("asm_ng50_min", "assembly: NG50, smaller haplotype (bp)", 0, "", ("assembly", "summary", "ng50_min"), None, "min"),
    ("asm_contigs", "assembly: contigs (both haplotypes)", 0, "", ("assembly", "summary", "contigs"), None, "max"),
    ("asm_flagger_non_hap", "assembly: Flagger non-Hap %", 2, "%", ("assembly", "summary", "flagger_non_hap_pct"), "asm_flagger_err_pct_max", "max"),
    ("asm_reliable_fraction", "assembly: reliable fraction of the reference", 3, "", ("assembly", "summary", "reliable_fraction"), "asm_reliable_fraction_min", "min"),
    ("asm_concordance", "assembly: genotype concordance with reads", 4, "", ("assembly", "summary", "genotype_concordance"), "asm_genotype_concordance_min", "min"),
    ("asm_het_loss", "assembly: het loss rate", 4, "", ("assembly", "summary", "het_loss"), "asm_het_loss_max", "max"),
    ("asm_switch_max", "assembly: switch error (trio, worse haplotype)", 4, "", ("assembly", "summary", "switch_max"), "asm_switch_error_max", "max"),
    ("asm_dipcall_snvs", "assembly: dipcall SNVs", 0, "", ("assembly", "summary", "dipcall_snvs"), None, "min"),
    ("asm_dipcall_hethom", "assembly: dipcall het/hom", 2, "", ("assembly", "summary", "dipcall_hethom"), None, "min"),
    ("asm_pav_records", "assembly: PAV records", 0, "", ("assembly", "summary", "pav_records"), None, "min"),
    ("asm_pav_sv_pass", "assembly: PAV PASS SVs (50 bp+)", 0, "", ("assembly", "summary", "pav_sv_pass"), None, "min"),
    ("asm_dipcall_pass_fraction", "assembly: dipcall PASS fraction", 3, "", ("assembly", "summary", "dipcall_pass_fraction"), None, "min"),
    ("asm_telo_contigs", "assembly: contigs telomeric at both ends (max over haplotypes)", 0, "", ("assembly", "summary", "telo_ends_max"), None, "min"),
    ("asm_chrx_fraction", "assembly: non-PAR chrX fraction", 3, "", ("assembly", "summary", "chrx_fraction"), None, "min"),
    ("asm_chry_fraction", "assembly: non-PAR chrY fraction", 3, "", ("assembly", "summary", "chry_fraction"), None, "min"),
]


def _dig(doc: object, path: tuple[str, ...]) -> object:
    cur = doc
    for key in path:
        if key == "_total" and isinstance(cur, dict):
            vals = [sv.finite(v) for v in cur.get("counts", {}).values()]
            return sum(vals) if vals else None
        if not isinstance(cur, dict):
            return None
        cur = cur.get(key)
    return cur


def _num(v: object) -> float | None:
    if v is None or isinstance(v, bool):
        return None
    try:
        x = float(v)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None
    return x if x == x else None


def supp_bins(supp: dict[str, object], n_samples: int) -> list[tuple[str, float]]:
    """The svx SUPP counts (samples carrying a merged SV -> records) in bins that read at any cohort size:
    1, 2, 3-5, 6-10, then shares of the cohort up to 1 %, 5 %, 10 %, 50 %, below all, and all members.
    Degenerate bins (equal edges at small cohorts) are dropped."""
    counts = {int(k): float(v) for k, v in supp.items() if str(k).isdigit()}
    if not counts:
        return []
    n = max(n_samples, max(counts))
    edges = [1, 2, 5, 10]
    for frac in (0.01, 0.05, 0.10, 0.50):
        edges.append(int(n * frac))
    edges.append(n - 1)
    uniq: list[int] = []
    for e in edges:
        if e >= 1 and (not uniq or e > uniq[-1]) and e < n:
            uniq.append(e)
    out: list[tuple[str, float]] = []
    lo = 1
    for hi in uniq:
        label = str(lo) if lo == hi else f"{lo}-{hi}"
        if hi > 10 and lo > 10:
            label += f" (≤{100 * hi / n:.0f}%)" if hi < n - 1 else " (<all)"
        out.append((label, sum(v for k, v in counts.items() if lo <= k <= hi)))
        lo = hi + 1
    out.append((f"all {n}", sum(v for k, v in counts.items() if k >= max(lo, n))))
    return [b for b in out if b[1] > 0 or b[0].startswith("1") or b[0].startswith("all")]


def _outliers(values: dict[str, float]) -> tuple[float | None, float | None, list[str]]:
    """Median, MAD and the ids beyond median ± 3 × 1.4826 × MAD (needs at least 4 values, a non-zero MAD)."""
    if len(values) < 4:
        return (statistics.median(values.values()) if values else None), None, []
    med = statistics.median(values.values())
    mad = statistics.median(abs(v - med) for v in values.values())
    if mad <= 0:
        return med, mad, []
    lim = 3 * 1.4826 * mad
    return med, mad, sorted(k for k, v in values.items() if abs(v - med) > lim)


def _manifest_of(run: RunRecord | None) -> dict[str, object]:
    if run is None:
        return {}
    p = run_files(run.run_path).manifest
    if p.exists():
        doc = read_json(p)
        if isinstance(doc, dict):
            return doc
    return {}


def _provenance(man: dict[str, object]) -> dict[str, object]:
    upstream = man.get("upstream") if isinstance(man.get("upstream"), dict) else {}
    references = man.get("references") if isinstance(man.get("references"), dict) else {}
    ugc_wgw_man = man.get("ugc_wgw_manifest") if isinstance(man.get("ugc_wgw_manifest"), dict) else {}
    return {"ugc_wgw": man.get("ugc_pacbio_wgw", {}), "engine": man.get("engine", {}), "host": man.get("host"),
            "upstream": {k: {"tag": v.get("tag"), "commit": v.get("commit")} for k, v in upstream.items() if isinstance(v, dict)},
            "workflow": ugc_wgw_man.get("upstream", {}) if isinstance(ugc_wgw_man, dict) else {},
            "references": {k: {"version": v.get("version"), "build": v.get("build"), "source": v.get("source"), "files": v.get("files", [])}
                           for k, v in references.items() if isinstance(v, dict)},
            "containers": sorted({str(c.get("image")) for c in man.get("containers", []) if isinstance(c, dict)})}


def _pangenome_entry(outs: Outputs, thr: dict[str, float] | None = None) -> dict[str, object]:
    """cohort_pangenome's small outputs: the chromosomes built with their graph sizes, the split's contig binning per
    genome, the frequency summary, the members table, truvari's figures, and the sizes of the VCFs."""
    entry: dict[str, object] = {}
    chroms = outs.opt("chromosomes_built")
    entry["chromosomes"] = [str(c) for c in chroms] if isinstance(chroms, list) else []
    per_chrom: dict[str, dict[str, object]] = {}
    stats = outs.opt("chrom_graph_stats")
    for p in (stats if isinstance(stats, list) else []):
        if not (isinstance(p, str) and Path(p).exists()):
            continue
        _cols, rows = sp.table(sp.read_lines(Path(p)))
        for r in rows:
            per_chrom[r.get("ref_chrom", "?")] = {k: sp._num(r.get(k, "")) for k in ("nodes", "edges", "length")}
    totals = {k: sum(float(c.get(k) or 0) for c in per_chrom.values()) for k in ("nodes", "edges", "length")}
    entry["graph"] = {"chromosomes": per_chrom, "nodes": int(totals["nodes"]), "edges": int(totals["edges"]), "length": int(totals["length"])}
    p = outs.opt("ambiguous_contigs")
    if isinstance(p, str) and Path(p).exists():
        _cols, rows = sp.table(sp.read_lines(Path(p)))
        contigs = [{k: (sp._num(v) if k != "genome" else v) for k, v in r.items()} for r in rows]
        entry["split"] = {"genomes": contigs,
                          "assigned_bp": sum(int(float(r.get("assigned_bp") or 0)) for r in contigs),  # type: ignore[arg-type]
                          "ambiguous_bp": sum(int(float(r.get("ambiguous_bp") or 0)) for r in contigs),  # type: ignore[arg-type]
                          "ambiguous_contigs": sum(int(float(r.get("ambiguous") or 0)) for r in contigs)}  # type: ignore[arg-type]
    p = outs.opt("freq_summary")
    if isinstance(p, str) and Path(p).exists():
        entry["summary"] = sp.parse_freq_summary(sp.read_lines(Path(p)))
    p = outs.opt("freq_samples")
    if isinstance(p, str) and Path(p).exists():
        entry["samples"] = sp.parse_freq_samples(sp.read_lines(Path(p)))
    p = outs.opt("truvari_summary")
    if isinstance(p, str) and Path(p).exists():
        tv = sp.parse_truvari_summary(read_json(Path(p)))
        if tv:
            entry["truvari"] = tv
    p = outs.opt("contig_sizes")
    if isinstance(p, str) and Path(p).exists():
        cs = sp.parse_contig_sizes(sp.read_lines(Path(p)))
        if cs:
            for chrom, g in per_chrom.items():
                c = cs["chromosomes"].get(chrom) if isinstance(cs.get("chromosomes"), dict) else None  # type: ignore[union-attr]
                if c:
                    g.update({"ref_length": c["ref_bp"], "hap_bp_min": c["hap_bp_min"], "hap_bp_avg": c["hap_bp_avg"], "hap_bp_max": c["hap_bp_max"]})
            entry["contig_sizes"] = {"genomes": cs.get("genomes", []), "assigned_bp": cs.get("assigned_bp", {}), "ambiguous_bp": cs.get("ambiguous_bp", {})}
    p = outs.opt("name_map")
    if isinstance(p, str) and Path(p).exists():
        entry["name_map"] = sp.parse_name_map(sp.read_lines(Path(p)))
    for name in ("pangenome_vcf", "pangenome_wave_vcf", "freq_vcf"):
        p = outs.opt(name)
        if isinstance(p, str) and Path(p).exists():
            entry[name] = {"bytes": Path(p).stat().st_size}
    # Tier 2 (0.12.0): per-member genotypes, the graph's per-chromosome genome sequence, degrees and reference gaps
    p = outs.opt("member_stats")
    if isinstance(p, str) and Path(p).exists():
        members = sp.parse_member_stats(sp.read_lines(Path(p)))
        entry["members"] = members
        limit = (thr or {}).get("pg_member_call_rate_min")
        dists: dict[str, object] = {}
        for key, label, digits in (("call_rate", "call rate in the graph VCF", 4), ("n_singletons", "singleton alleles", 0), ("n_hets", "het genotypes", 0),
                                   ("n_nonref_hom", "hom-alt genotypes", 0), ("n_missing", "missing genotypes", 0)):
            values = {str(m["sample_id"]): float(m[key]) for m in members if isinstance(m.get(key), (int, float))}
            if not values:
                continue
            med, mad, outliers = _outliers(values)
            beyond = sorted(k for k, v in values.items() if key == "call_rate" and limit is not None and v < limit)
            dists[key] = {"label": label, "digits": digits, "values": values, "median": med, "mad": mad, "outliers": outliers,
                          "threshold": limit if key == "call_rate" else None, "direction": "min", "beyond_threshold": beyond}
        entry["member_distributions"] = dists
    chroms_built: list[str] = entry.get("chromosomes") or []  # type: ignore[assignment]
    for name, parser, key in (("chrom_genome_bp", lambda q: sp.parse_genome_bp(sp.read_lines(q)), "genome_bp"),
                              ("chrom_refgaps", lambda q: sp.parse_refgaps(sp.read_binary_lines(q)), "refgaps")):
        paths = outs.opt(name)
        if isinstance(paths, list):
            for i, q in enumerate(paths):
                if not (isinstance(q, str) and Path(q).exists()):
                    continue
                chrom = chroms_built[i] if i < len(chroms_built) else Path(q).name
                parsed = parser(Path(q))
                if key == "refgaps":
                    g = per_chrom.setdefault(chrom, {})
                    tot = {"n": sum(v["n"] for v in parsed.values()), "bp": sum(v["bp"] for v in parsed.values())}  # type: ignore[index]
                    g["refgaps"] = tot
                else:
                    per_chrom.setdefault(chrom, {}).update({"genome_" + k: v for k, v in parsed.items()})  # type: ignore[union-attr]
    paths = outs.opt("chrom_degree_stats")
    if isinstance(paths, list):
        degree: dict[str, int] = {}
        for q in paths:
            if isinstance(q, str) and Path(q).exists():
                for k, v in sp.parse_degree_stats(sp.read_lines(Path(q))).items():
                    degree[k] = degree.get(k, 0) + v
        if degree:
            entry["degree"] = degree
    return entry


def _freq_resource_html(title: str, r: dict[str, object]) -> str:
    """The cards and charts of one resource of a frequency summary (cohort_freq's small_variants or structural_variants,
    cohort_pangenome's pangenome)."""
    rec: dict[str, object] = r.get("records", {})  # type: ignore[assignment]
    return (f"<h3>{_esc(title)}</h3>"
            + _cards([("records", sv.fmt(rec.get("total"))), ("with AC > 0", sv.fmt(rec.get("ac_gt0"))), ("singletons", sv.fmt(rec.get("singletons"))),
                      ("mean call rate", sv.fmt((r.get("mean_call_rate") or {}).get("all"), 4)), ("call rate < 0.9", sv.fmt(rec.get("call_rate_below_0.9"))),  # type: ignore[union-attr]
                      ("with AC > 0 in XX", sv.fmt((r.get("records_ac_gt0_in") or {}).get("XX"))), ("in XY", sv.fmt((r.get("records_ac_gt0_in") or {}).get("XY")))])  # type: ignore[union-attr]
            + '<div class="grid"><div><h4>Allele-frequency bins</h4>' + sv.hist([(k, sv.finite(v)) for k, v in (r.get("records_by_af_bin") or {}).items()], width=420, y_label="records")  # type: ignore[union-attr]
            + "</div><div><h4>By type</h4>" + sv.bars_h([(k, sv.finite(v), "") for k, v in sorted((r.get("records_by_type") or {}).items())], width=300, label_w=70)  # type: ignore[union-attr]
            + "</div><div><h4>By contig class</h4>" + sv.bars_h([(k, sv.finite(v), "") for k, v in sorted((r.get("records_by_contig_class") or {}).items())], width=300, label_w=70)  # type: ignore[union-attr]
            + "</div><div><h4>By filter</h4>" + sv.bars_h([(k, sv.finite(v), "") for k, v in sorted((r.get("records_by_filter") or {}).items())], width=300, label_w=110) + "</div></div>")  # type: ignore[union-attr]


def _assemblies_html(asm: dict[str, dict[str, object]], link, open_below: int = 50) -> str:  # type: ignore[no-untyped-def]
    """One row per member: both haplotypes' contiguity side by side, the trio QC, the karyotype, Flagger, the reliable
    fraction, the concordance and both callers' counts; collapsed above `open_below` members."""
    def pair(vals: list[object], fmt) -> str:  # type: ignore[no-untyped-def]
        return " / ".join(fmt(v) for v in vals) if vals else "–"

    rows = []
    for sid, a in asm.items():
        haps: list[dict[str, object]] = a.get("haplotypes") or []  # type: ignore[assignment]
        st = [h.get("stats") or {} for h in haps]
        te = [h.get("trio_eval") or {} for h in haps]
        calls: dict[str, dict[str, object]] = a.get("calls") or {}  # type: ignore[assignment]
        dip, pav = calls.get("dipcall", {}), calls.get("pav", {})
        kar: dict[str, object] = a.get("karyotype") or {}  # type: ignore[assignment]
        fl: dict[str, object] = a.get("flagger") or {}  # type: ignore[assignment]
        rel: dict[str, object] = a.get("reliable") or {}  # type: ignore[assignment]
        con: dict[str, object] = a.get("concordance") or {}  # type: ignore[assignment]
        rows.append([link(sid), "trio" if a.get("trio") else "single", kar.get("karyotype") or "–",
                     pair([s.get("total_length") for s in st], sv.fmt_bp), pair([s.get("n_sequences") for s in st], sv.fmt),
                     pair([s.get("n50") for s in st], lambda v: sv.fmt_bp(v) if v else "–"), pair([s.get("aun") for s in st], sv.fmt_bp),
                     pair([t.get("switch_rate") for t in te], lambda v: sv.fmt(v, 4)) if any(te) else "–",
                     pair([t.get("hamming_rate") for t in te], lambda v: sv.fmt(v, 4)) if any(te) else "–",
                     sv.fmt(fl.get("non_hap_pct"), 2) if fl else "–", sv.fmt(rel.get("reliable_fraction_of_reference"), 3) if rel else "–",
                     sv.fmt(con.get("genotype_concordance"), 4) if con else "–", sv.fmt(con.get("het_loss_rate"), 4) if con else "–",
                     f"{sv.fmt(dip.get('snps'))} / {sv.fmt(dip.get('indels'))}" if dip else "–", f"{sv.fmt(pav.get('snps'))} / {sv.fmt(pav.get('indels'))}" if pav else "–",
                     f"{sv.fmt(kar.get('chrx_fraction'), 3)} / {sv.fmt(kar.get('chry_fraction'), 3)}" if kar.get("chrx_fraction") is not None else "–"])
    cols = ["sample", "mode", "karyotype", "length h1 / h2", "contigs", "NG50", "auNG", "switch error", "Hamming", "Flagger non-Hap %", "reliable fraction",
            "concordance", "het loss", "dipcall SNVs / indels", "PAV SNVs / indels", "chrX / chrY fraction"]
    state = " open" if len(rows) <= open_below else ""
    return (f"<details{state}><summary>{len(rows)} assemblies (one row per member; haplotype values h1 / h2)</summary>"
            + '<div class="wrap">' + _table(cols, rows) + "</div></details>"
            + "<p class='small muted'>Flagger non-Hap: Err + Dup + Col share of the diploid assembly; reliable fraction: dip.bed minus Flagger's regions over the reference; "
              "concordance and het loss against the member's read-based SNVs. The assembly metrics also appear in the distributions above.</p>")


SEVERITY_RANK = {"pass": 0, "info": 1, "warn": 2, "fail": 3}
IDENTITY_SAMPLE_CODES = ("sex_mismatch", "sex_not_inferred")   # the per-sample flags that are identity questions too


def _worst(severities: list[str]) -> str:
    return max(severities, key=lambda s: SEVERITY_RANK.get(s, 0), default="pass")


def identity_flags(pairs: list[dict[str, object]], members_tbl: list[dict[str, object]], digests: dict[str, dict[str, object]],
                   order: list[str], thr: dict[str, float]) -> dict[str, list[dict[str, object]]]:
    """The cohort-level identity flags per member, from cohort_qc's tables (docs/guide/07-results.md, "Identity"):
    `duplicate` (fail, both members of a pair at or above duplicate_relatedness_min), `pedigree_mismatch` (fail, child and
    parent of a declared pair outside pedigree_relatedness_min..max, or missing from the table), `related` (warn, both
    members of an undeclared pair at or above related_relatedness_max), `related_hub` (fail, a member related that way to
    related_hub_min members or more: contamination, whose pairs count against nobody else), `het_outlier` (warn, heterozygosity ratio beyond
    median ± het_ratio_mad_max scaled MADs, from four members), and `sex_mismatch` (warn) for a member whose own run
    predates the somalier extract when cohort_qc's chrX call disagrees with the sheet or the coverage inference."""
    flags: dict[str, list[dict[str, object]]] = {sid: [] for sid in order}

    def flag(sid: str, code: str, severity: str, message: str, value: object = None, threshold: object = None) -> None:
        if sid in flags:
            flags[sid].append({"code": code, "severity": severity, "message": message, "value": value, "threshold": threshold})

    dup_min, ped_min, ped_max, rel_max = (thr["duplicate_relatedness_min"], thr["pedigree_relatedness_min"],
                                          thr["pedigree_relatedness_max"], thr["related_relatedness_max"])
    hub_min = int(thr.get("related_hub_min", 10))
    # a member related to many others is the anomaly, not its partners: contamination inflates the shared heterozygous
    # sites and so somalier's relatedness to everyone (the 864 cohort, 2026-11-09: one sample "related" to 858 of 863).
    # The hubs are taken out one by one, the most connected first, until nobody reaches related_hub_min; their pairs do
    # not count against anyone
    undeclared: list[tuple[str, str, float]] = []
    for p in pairs:
        a, b, r = str(p.get("sample_a") or ""), str(p.get("sample_b") or ""), p.get("relatedness")
        expected = p.get("expected")
        if isinstance(r, (int, float)) and rel_max <= r < dup_min and not (isinstance(expected, (int, float)) and expected > 0):
            undeclared.append((a, b, float(r)))
    hubs: dict[str, int] = {}
    while undeclared:
        degree: dict[str, int] = {}
        for a, b, _r in undeclared:
            degree[a] = degree.get(a, 0) + 1
            degree[b] = degree.get(b, 0) + 1
        sid, n = max(degree.items(), key=lambda kv: (kv[1], kv[0]))
        if n < hub_min:
            break
        hubs[sid] = n
        undeclared = [t for t in undeclared if sid not in t[:2]]
    for sid, n in hubs.items():
        flag(sid, "related_hub", "fail", f"related at or above {rel_max:g} to {n} members with no pedigree link in the sheet (at least "
             f"{hub_min}): a contaminated or mixed sample looks related to everyone; its partners are not flagged on its account", n, hub_min)
    rel: dict[tuple[str, str], float] = {}
    for p in pairs:
        a, b, r = str(p.get("sample_a") or ""), str(p.get("sample_b") or ""), p.get("relatedness")
        if not isinstance(r, (int, float)):
            continue
        rel[(a, b)] = rel[(b, a)] = float(r)
        expected = p.get("expected")
        declared = isinstance(expected, (int, float)) and expected > 0
        if r >= dup_min:
            for x, y in ((a, b), (b, a)):
                flag(x, "duplicate", "fail", f"same individual as {y} (relatedness {r:.3f}, at or above {dup_min:g})", r, dup_min)
        elif r >= rel_max and not declared and a not in hubs and b not in hubs:
            for x, y in ((a, b), (b, a)):
                flag(x, "related", "warn", f"related to {y} (relatedness {r:.3f}) with no pedigree link in the sheet", r, rel_max)
    for sid in order:
        sample: dict[str, object] = (digests.get(sid) or {}).get("sample") or {}  # type: ignore[assignment]
        for role in ("father", "mother"):
            parent = str(sample.get(f"{role}_id") or "")
            if not parent or parent not in flags or sid not in digests:
                continue   # a parent outside the cohort is not checked
            r = rel.get((sid, parent))
            if r is None:
                flag(sid, "pedigree_mismatch", "warn", f"declared {role} {parent} is in the cohort but the pair is not in the relatedness table")
            elif not ped_min <= r <= ped_max:
                side = "below" if r < ped_min else "above"
                flag(sid, "pedigree_mismatch", "fail", f"declared {role} {parent} but relatedness {r:.3f} is {side} the parent-child range "
                     f"{ped_min:g}–{ped_max:g}", r, [ped_min, ped_max])
                flag(parent, "pedigree_mismatch", "fail", f"declared {role} of {sid} but relatedness {r:.3f} is {side} the parent-child range "
                     f"{ped_min:g}–{ped_max:g}", r, [ped_min, ped_max])
    het = {str(r.get("sample_id")): float(r["het_ratio"]) for r in members_tbl if isinstance(r.get("het_ratio"), (int, float))}  # type: ignore[arg-type]
    if len(het) >= 4:
        med = statistics.median(het.values())
        mad = statistics.median(abs(v - med) for v in het.values())
        if mad > 0:
            lim = thr["het_ratio_mad_max"] * 1.4826 * mad
            for sid, v in het.items():
                if abs(v - med) > lim:
                    flag(sid, "het_outlier", "warn", f"heterozygosity ratio {v:.3f} is beyond the cohort's median {med:.3f} ± {thr['het_ratio_mad_max']:g} MAD "
                         f"({med - lim:.3f}–{med + lim:.3f}): contamination or a sample mix", v, [med - lim, med + lim])
    for r in members_tbl:
        sid = str(r.get("sample_id") or "")
        d = digests.get(sid)
        if not d or not ((d.get("identity") or {}).get("legacy")):  # type: ignore[union-attr]
            continue   # members with their own somalier call carry the flag in their digest
        som_raw = str(r.get("somalier_sex") or "")
        som = SEX_TO_GROUP.get(som_raw.upper())
        if not som:
            continue
        sample = d.get("sample") or {}  # type: ignore[assignment]
        cov: dict[str, object] = d.get("coverage") or {}  # type: ignore[assignment]
        others = [(name, raw, SEX_TO_GROUP.get(str(raw or "").upper())) for name, raw in (("the sheet", sample.get("sex_sheet")), ("coverage", cov.get("inferred_sex")))]
        dis = [(name, raw, g) for name, raw, g in others if g and g != som]
        if dis:
            flag(sid, "sex_mismatch", "warn", f"somalier's chrX call in cohort_qc is {som_raw} ({som}) but " + ", ".join(f"{n} says {raw} ({g})" for n, raw, g in dis),
                 som_raw, sample.get("sex_sheet"))
    return flags


def identity_pass(d: dict[str, object], sheet_order: dict[str, int]) -> tuple[list[str], list[tuple[str, str]]]:
    """The members that pass the identity checks of a cohort digest, in sheet order, and the exclusions with their
    reason: a member whose identity status is fail is out; of an undeclared related pair (relatedness at or above
    related_relatedness_max, no expectation from the sheet's pedigree) the later member in the sheet is out. A pair
    the sheet declares (parent and child) stays: the pedigree is the analyst's statement, the surprises are the
    undeclared relatives and the duplicates (2026-11-05, the user's call)."""
    ident: dict[str, object] = d.get("identity", {}) or {}  # type: ignore[assignment]
    flags: dict[str, list[dict[str, object]]] = ident.get("flags", {}) or {}  # type: ignore[assignment]
    rows: list[dict[str, object]] = [r for r in d.get("members", []) if r.get("present")]  # type: ignore[union-attr]
    order = sorted((str(r["id"]) for r in rows), key=lambda s: sheet_order.get(s, len(sheet_order)))
    status = {str(r["id"]): str(r.get("identity_status", "pass")) for r in rows}
    thr: dict[str, float] = d.get("thresholds", {}) or {}  # type: ignore[assignment]
    rel_max = float(thr.get("related_relatedness_max", 0.2))
    pairs: list[dict[str, object]] = ident.get("pairs", []) or []  # type: ignore[assignment]
    hubs = {sid for sid, fl in flags.items() if any(f.get("code") == "related_hub" for f in fl)}   # out as failures; their pairs count against nobody
    pairs = [p for p in pairs if p.get("sample_a") not in hubs and p.get("sample_b") not in hubs]
    kept: list[str] = []
    excluded: list[tuple[str, str]] = []
    for sid in order:
        if status.get(sid) == "fail":
            why = "; ".join(str(f.get("message")) for f in flags.get(sid, []) if f.get("severity") == "fail") or "an identity failure"
            excluded.append((sid, why))
            continue
        partner = None
        for p in pairs:
            r, exp = p.get("relatedness"), p.get("expected")
            if not isinstance(r, (int, float)) or r < rel_max or (isinstance(exp, (int, float)) and exp > 0):
                continue   # below the threshold, or declared in the sheet
            other = p.get("sample_b") if p.get("sample_a") == sid else p.get("sample_a") if p.get("sample_b") == sid else None
            if other is not None and str(other) in kept:
                partner = (str(other), float(r))
                break
        if partner:
            excluded.append((sid, f"related to {partner[0]} (relatedness {partner[1]:.3f}) with no pedigree link in the sheet; {partner[0]} comes first"))
            continue
        kept.append(sid)
    return kept, excluded


def build_cohort_digest(cohort: Cohort, members: dict[str, dict[str, object]], cruns: dict[str, RunRecord], mode: str,
                        thr: dict[str, float], *, genes_of_interest: list[str] | None = None) -> dict[str, object]:
    warnings: list[str] = []
    rows: list[dict[str, object]] = []
    for sid in cohort.members:
        m = members.get(sid)
        if not m:
            rows.append({"id": sid, "present": False})
            continue
        qc = m.get("qc", {}) if isinstance(m.get("qc"), dict) else {}
        sample = m.get("sample", {}) if isinstance(m.get("sample"), dict) else {}
        cov = m.get("coverage", {}) if isinstance(m.get("coverage"), dict) else {}
        ident = m.get("identity", {}) if isinstance(m.get("identity"), dict) else {}
        row: dict[str, object] = {"id": sid, "present": True, "sex_sheet": sample.get("sex_sheet"), "sex_inferred": cov.get("inferred_sex"),
                                  "sex_somalier": ident.get("somalier") or "", "sex_used": ident.get("resolved") or "", "sex_source": ident.get("source") or "",
                                  "sex_legacy": bool(ident.get("legacy")),
                                  "identity_status": _worst([str(f.get("severity")) for f in qc.get("flags", []) if f.get("code") in IDENTITY_SAMPLE_CODES]),  # type: ignore[union-attr]
                                  "qc_status": qc.get("status", "pass"), "flags": sorted({str(f.get("code")) for f in qc.get("flags", [])}),  # type: ignore[union-attr]
                                  "runs": {s: r.get("run_id") for s, r in (m.get("runs") or {}).items() if isinstance(r, dict)},  # type: ignore[union-attr]
                                  "ugc_wgw_version": next((r.get("ugc_wgw_version") for r in (m.get("runs") or {}).values() if isinstance(r, dict)), None),  # type: ignore[union-attr]
                                  "link": f"{sid}.summary.html", "has_reads": _has_read_data(m)}
        for key, _label, _d, _u, path, _t, _dir in MEMBER_METRICS:
            row[key] = _num(_dig(m, path))
        targeted = m.get("targeted", {}) if isinstance(m.get("targeted"), dict) else {}
        kivvi = targeted.get("kivvi", {}) if isinstance(targeted.get("kivvi"), dict) else {}
        row["karyotype"] = ((m.get("assembly") or {}).get("karyotype") or {}).get("karyotype") if isinstance(m.get("assembly"), dict) else None  # type: ignore[union-attr]
        row["kivvi_kiv2_cn"] = _num((kivvi.get("kiv2") or {}).get("allele_cn"))
        row["kivvi_d4z4_cn"] = _num((kivvi.get("d4z4") or {}).get("allele_cn"))
        para = targeted.get("paraphase", {}) if isinstance(targeted.get("paraphase"), dict) else {}
        row["paraphase_failed"] = para.get("n_failed")
        rows.append(row)
    present = [r for r in rows if r.get("present")]
    distributions: dict[str, object] = {}
    for key, label, digits, unit, _path, tkey, direction in MEMBER_METRICS:
        values = {str(r["id"]): float(r[key]) for r in present if r.get(key) is not None}  # type: ignore[arg-type]
        if not values:
            continue
        med, mad, outliers = _outliers(values)
        limit = thr.get(tkey) if tkey else None
        beyond = sorted(k for k, v in values.items() if limit is not None and ((v < limit) if direction == "min" else (v > limit)))
        distributions[key] = {"label": label, "digits": digits, "unit": unit, "values": values, "median": med, "mad": mad,
                              "outliers": outliers, "threshold": limit, "direction": direction, "beyond_threshold": beyond,
                              "below_threshold": beyond if direction == "min" else []}
    # cohort stages' own outputs
    callsets: dict[str, object] = {}
    for stage, run in cruns.items():
        try:
            outs = Outputs.load(run)
        except UgcError as exc:
            warnings.append(f"{stage}: {exc}")
            continue
        entry: dict[str, object] = {"run_id": run.run_id, "attempt": run.attempt, "finished_at": run.finished_at}
        try:
            if stage == "cohort_call":
                p = outs.opt("cohort_small_variant_vcf")
                if isinstance(p, str) and Path(p).exists():
                    entry["small_variant_vcf"] = {"path": p, "bytes": Path(p).stat().st_size, "samples": sp.vcf_samples(Path(p))}
            elif stage == "cohort_merge":
                p = outs.opt("cohort_sv_vcf")
                if isinstance(p, str) and Path(p).exists():
                    entry["sv"] = sp.parse_sv_vcf(sp.read_binary_lines(Path(p)))
                    entry["sv"]["bytes"] = Path(p).stat().st_size  # type: ignore[index]
                p = outs.opt("cohort_trgt_lps")
                if isinstance(p, str) and Path(p).exists():
                    cols, n = sp.tsv_dims(Path(p))
                    entry["trgt_lps"] = {"loci": n, "samples": cols[2:], "bytes": Path(p).stat().st_size}
                p = outs.opt("cohort_trgt_vcf")
                if isinstance(p, str) and Path(p).exists():
                    entry["trgt_vcf"] = {"bytes": Path(p).stat().st_size, "samples": sp.vcf_samples(Path(p))}
                p = outs.opt("cohort_small_variant_vcf")
                if isinstance(p, str) and Path(p).exists():
                    entry["small_variant_vcf"] = {"path": p, "bytes": Path(p).stat().st_size, "samples": sp.vcf_samples(Path(p))}
                p = outs.opt("cohort_mt_small_variant_vcf")
                if isinstance(p, str) and Path(p).exists():
                    entry["mt"] = sp.parse_mt_vcf(sp.read_binary_lines(Path(p)))
                    entry["mt"]["bytes"] = Path(p).stat().st_size  # type: ignore[index]
                p = outs.opt("cohort_mt_sv_vcf")
                if isinstance(p, str) and Path(p).exists():
                    entry["mt_svs"] = {"bytes": Path(p).stat().st_size, "records": sp.parse_mt_vcf(sp.read_binary_lines(Path(p))).get("records")}
                p = outs.opt("cohort_mt_haplogroups")
                if isinstance(p, str) and Path(p).exists():
                    entry["mt_haplogroups"] = sp.parse_haplogroups(sp.read_lines(Path(p)))
            elif stage == "cohort_pangenome":
                entry.update(_pangenome_entry(outs, thr))
            elif stage == "cohort_qc":
                p = outs.opt("qc_summary")
                if isinstance(p, str) and Path(p).exists():
                    entry["summary"] = sp.parse_qc_summary(sp.read_lines(Path(p)))
                p = outs.opt("qc_members")
                if isinstance(p, str) and Path(p).exists():
                    entry["members"] = sp.parse_qc_members(sp.read_lines(Path(p)))
                p = outs.opt("qc_related_pairs")
                if isinstance(p, str) and Path(p).exists():
                    entry["related_pairs"] = sp.parse_somalier_pairs(sp.read_lines(Path(p)))
                p = outs.opt("qc_ancestry")
                if isinstance(p, str) and Path(p).exists():
                    arows = sp.parse_ancestry(sp.read_lines(Path(p)))
                    # the labelled set as a light background (label, PC1, PC2); the members with their probabilities
                    entry["ancestry_background"] = [[r["given"], r["pcs"][0], r["pcs"][1]] for r in arows if r.get("given") and len(r.get("pcs", [])) >= 2]  # type: ignore[index, arg-type]
                    entry["ancestry_members"] = [{"sample_id": r["sample_id"], "predicted": r["predicted"], "probs": r["probs"], "pcs": r["pcs"][:2]}  # type: ignore[index]
                                                 for r in arows if not r.get("given") and r.get("sample_id") in members]
                    entry["ancestry"] = {str(r["sample_id"]): str(r["predicted"]) for r in entry["ancestry_members"] if r.get("predicted")}  # type: ignore[union-attr, index]
                entry["legacy_extracted"] = [str(r.get("sample_id")) for r in entry.get("members", []) if str(r.get("extracted_here")) == "yes"]  # type: ignore[union-attr]
                for name in ("qc_html", "qc_ancestry_html"):
                    p = outs.opt(name)
                    if isinstance(p, str) and Path(p).exists():
                        entry[name] = {"path": p, "bytes": Path(p).stat().st_size}
            elif stage == "cohort_methyl":
                p = outs.opt("methyl_summary")
                if isinstance(p, str) and Path(p).exists():
                    entry["summary"] = sp.parse_methyl_summary(sp.read_lines(Path(p)))
                p = outs.opt("methyl_samples")
                if isinstance(p, str) and Path(p).exists():
                    entry["samples"] = sp.parse_methyl_samples(sp.read_lines(Path(p)))
                p = outs.opt("methyl_imprinting")
                if isinstance(p, str) and Path(p).exists():
                    entry["imprinting"] = sp.table(sp.read_lines(Path(p)))[1]
                p = outs.opt("methyl_regions")
                if isinstance(p, str) and Path(p).exists():
                    entry["regions"] = sp.parse_methyl_regions(sp.read_lines(Path(p)))
                p = outs.opt("methyl_events")
                if isinstance(p, str) and Path(p).exists():
                    entry["events"] = sp.parse_methyl_events(sp.read_lines(Path(p)), genes=genes_of_interest or ())
                    entry["events"]["bytes"] = Path(p).stat().st_size  # type: ignore[index]
                for name in ("cohort_background", "methyl_collection", "joint_segments"):
                    p = outs.opt(name)
                    if isinstance(p, str) and Path(p).exists():
                        entry[name] = {"path": p, "bytes": Path(p).stat().st_size}
                entry["genes_of_interest"] = list(genes_of_interest or [])
            elif stage == "cohort_freq":
                p = outs.opt("freq_summary")
                if isinstance(p, str) and Path(p).exists():
                    entry["summary"] = sp.parse_freq_summary(sp.read_lines(Path(p)))
                p = outs.opt("freq_samples")
                if isinstance(p, str) and Path(p).exists():
                    entry["samples"] = sp.parse_freq_samples(sp.read_lines(Path(p)))
                for name in ("small_variant_freq_vcf", "sv_freq_vcf", "mt_small_variant_freq_vcf", "mt_sv_freq_vcf"):
                    p = outs.opt(name)
                    if isinstance(p, str) and Path(p).exists():
                        entry[name] = {"bytes": Path(p).stat().st_size}
        except Exception as exc:  # noqa: BLE001 - one bad file must not sink the report
            warnings.append(f"{stage}: {type(exc).__name__}: {exc}")
        callsets[stage] = entry

    assembly = {sid: members[sid].get("assembly") for sid in cohort.members if sid in members and members[sid].get("assembly")}
    present_members = {sid: members[sid] for sid in cohort.members if sid in members}

    # identity: cohort_qc's tables, the cohort-level flags, merged into the member rows
    qc_entry: dict[str, object] = callsets.get("cohort_qc") if isinstance(callsets.get("cohort_qc"), dict) else {}  # type: ignore[assignment]
    identity: dict[str, object] = {"available": bool(qc_entry.get("members") is not None and qc_entry.get("related_pairs") is not None)}
    if identity["available"]:
        qc_members: list[dict[str, object]] = qc_entry.get("members", [])  # type: ignore[assignment]
        qc_pairs: list[dict[str, object]] = qc_entry.get("related_pairs", [])  # type: ignore[assignment]
        iflags = identity_flags(qc_pairs, qc_members, present_members, list(cohort.members), thr)
        by_id = {str(r.get("sample_id")): r for r in qc_members}
        for r in present:
            fl = iflags.get(str(r["id"]), [])
            r["flags"] = sorted(set(r["flags"]) | {str(f["code"]) for f in fl})  # type: ignore[arg-type]
            r["qc_status"] = _worst([str(r.get("qc_status", "pass"))] + [str(f["severity"]) for f in fl])
            r["identity_status"] = _worst([str(r.get("identity_status", "pass"))] + [str(f["severity"]) for f in fl])
            qr = by_id.get(str(r["id"]))
            if qr is not None and r.get("sex_legacy") and not r.get("sex_somalier"):
                r["sex_somalier"] = qr.get("somalier_sex") or ""
        counts: dict[str, int] = {}
        for fl in iflags.values():
            for code in {str(f["code"]) for f in fl}:
                counts[code] = counts.get(code, 0) + 1
        identity.update({"flags": iflags, "counts": counts, "pairs": qc_pairs,
                         "members": {sid: dict(by_id.get(sid, {"sample_id": sid}), flags=iflags.get(sid, [])) for sid in cohort.members},
                         "n_fail": sum(1 for r in present if r.get("identity_status") == "fail"),
                         "legacy_extracted": qc_entry.get("legacy_extracted", [])})

    # methylation: cohort_methyl's samples table joined to the member rows (events per member, the MAD outliers)
    me_entry: dict[str, object] = callsets.get("cohort_methyl") if isinstance(callsets.get("cohort_methyl"), dict) else {}  # type: ignore[assignment]
    methylation: dict[str, object] = {"available": bool(me_entry.get("samples") is not None)}
    if methylation["available"]:
        mrows: list[dict[str, object]] = me_entry.get("samples", [])  # type: ignore[assignment]
        by_id_m = {str(r.get("sample_id")): r for r in mrows}
        values = {sid: float(r["n_events"]) for sid, r in by_id_m.items() if isinstance(r.get("n_events"), (int, float)) and sid in present_members}
        mflags: dict[str, list[dict[str, object]]] = {}
        med, mad, _ = _outliers(values) if values else (None, None, [])
        limit = (med + thr.get("methyl_events_mad_max", 3.0) * mad) if (med is not None and mad is not None and len(values) >= 4) else None
        for sid, v in values.items():
            if limit is not None and v > limit and mad:   # a constant cohort (mad 0) flags nobody
                mflags.setdefault(sid, []).append({"code": "methyl_outlier", "severity": "warn",
                                                   "message": f"{int(v)} methylation events against the cohort background, above the cohort's median {sv.fmt(med)} + "
                                                              f"{sv.fmt(thr.get('methyl_events_mad_max', 3.0), 1)} MAD ({sv.fmt(limit, 1)})", "value": v, "threshold": limit})
        for r in present:
            fl = mflags.get(str(r["id"]), [])
            if fl:
                r["flags"] = sorted(set(r["flags"]) | {str(f["code"]) for f in fl})  # type: ignore[arg-type]
                r["qc_status"] = _worst([str(r.get("qc_status", "pass"))] + [str(f["severity"]) for f in fl])
            mr = by_id_m.get(str(r["id"]))
            r["methyl_events"] = _num(mr.get("n_events")) if mr else None
            r["imprinting_anomalous"] = _num(mr.get("imprinting_anomalous")) if mr else None
        if values:
            distributions["methyl_events"] = {"label": "methylation events (cohort_methyl)", "digits": 0, "unit": "", "values": values, "median": med, "mad": mad,
                                              "outliers": sorted(mflags), "threshold": None, "direction": "max", "beyond_threshold": [], "below_threshold": []}
        msumm: dict[str, dict[str, object]] = me_entry.get("summary") if isinstance(me_entry.get("summary"), dict) else {}  # type: ignore[assignment]
        methylation.update({"flags": mflags, "counts": {"methyl_outlier": len(mflags)},
                            "members": {sid: by_id_m.get(sid, {"sample_id": sid}) for sid in cohort.members},
                            "legacy_reported": [sid for sid, r in by_id_m.items() if str(r.get("report_source")) == "stage"],
                            # each comparison's summary_comparison counts and where its differing regions lie (chrX, chrY, autosomal)
                            "comparisons": {pair: {"counts": counts, "contigs": (msumm.get("compare_contigs") or {}).get(pair, {})}
                                            for pair, counts in (msumm.get("compare") or {}).items()}})

    sex: dict[str, object] = {"sheet": {}, "inferred": {}, "somalier": {}, "used": {}, "karyotype": {},
                              "mismatches": [r["id"] for r in present if "sex_mismatch" in r["flags"]],  # type: ignore[operator]
                              "karyotype_mismatches": [r["id"] for r in present if "asm_karyotype_mismatch" in r["flags"]],  # type: ignore[operator]
                              "not_inferred": [r["id"] for r in present if not r.get("sex_inferred") and r.get("has_reads")],
                              "unknown": [r["id"] for r in present if not r.get("sex_used") and r.get("has_reads")],
                              "sources": {}}
    for r in present:
        for kind, key in (("sex_sheet", "sheet"), ("sex_inferred", "inferred"), ("sex_somalier", "somalier"), ("sex_used", "used")):
            k = str(r.get(kind) or "unknown").upper()
            sex[key][k] = sex[key].get(k, 0) + 1  # type: ignore[index, call-overload]
        if r.get("has_reads"):
            src = str(r.get("sex_source") or "none")
            sex["sources"][src] = sex["sources"].get(src, 0) + 1  # type: ignore[index, call-overload]
        if r.get("karyotype"):
            k = str(r["karyotype"])
            sex["karyotype"][k] = sex["karyotype"].get(k, 0) + 1  # type: ignore[index, call-overload]
    ancestry = qc_entry.get("ancestry")
    try:
        cm_entry: dict[str, object] = callsets.get("cohort_merge") if isinstance(callsets.get("cohort_merge"), dict) else {}  # type: ignore[assignment]
        targeted = export.targeted_summary(present_members, export.strata_of(present_members, ("sex", "ancestry"),
                                                                             ancestry=ancestry if isinstance(ancestry, dict) else None),
                                           cohort_haplogroups=cm_entry.get("mt_haplogroups"))  # type: ignore[arg-type]
    except Exception as exc:  # noqa: BLE001 - one bad digest must not sink the page
        warnings.append(f"targeted callers: {type(exc).__name__}: {exc}")
        targeted = {}
    versions = {"ugc_wgw": sorted({str(r.get("ugc_wgw_version")) for r in present if r.get("ugc_wgw_version")}),
                "digest_schema": sorted({str(members[sid].get("schema")) for sid in members})}
    man = next((_manifest_of(cruns.get(s)) for s in ("cohort_freq", "cohort_merge", "cohort_call", "cohort_pangenome", "cohort_qc") if cruns.get(s)), {})
    prov = _provenance(man) if man else next((m.get("provenance", {}) for m in members.values()), {})
    return {"schema": DIGEST_SCHEMA, "generated_at": utc_now(), "mode": mode,
            "cohort": {"id": cohort.cohort_id, "frozen_at": cohort.frozen_at, "sha256": cohort.sample_list_sha256, "n": len(cohort.members)},
            "runs": {s: {"run_id": r.run_id, "ugc_wgw_version": r.ugc_wgw_version, "attempt": r.attempt, "finished_at": r.finished_at} for s, r in cruns.items()},
            "members": rows, "missing_members": [r["id"] for r in rows if not r.get("present")],
            "distributions": distributions, "sex": sex, "identity": identity, "methylation": methylation, "callsets": callsets, "targeted": targeted,
            "assembly": assembly,
            "versions": versions,
            "thresholds": dict(thr), "provenance": prov, "warnings": warnings}


def _identity_html(ident: dict[str, object], cq: dict[str, object], present: list[dict[str, object]], thr: dict[str, float], link) -> str:  # type: ignore[no-untyped-def]
    """The Identity section: cards, the relatedness histogram, the pairs worth a look, the flags, the ancestry projection."""
    counts: dict[str, int] = ident.get("counts", {}) or {}  # type: ignore[assignment]
    summ: dict[str, dict[str, dict[str, object]]] = cq.get("summary", {}) or {}  # type: ignore[assignment]
    pairs_n = (summ.get("pairs") or {}).get("n", {})
    members_n = (summ.get("members") or {}).get("n", {})
    hist_raw: dict[str, object] = (summ.get("pairs") or {}).get("hist", {})
    legacy: list[str] = ident.get("legacy_extracted", []) or []  # type: ignore[assignment]
    dup_min, ped_min, ped_max, rel_max = (thr.get("duplicate_relatedness_min", 0.9), thr.get("pedigree_relatedness_min", 0.35),
                                          thr.get("pedigree_relatedness_max", 0.65), thr.get("related_relatedness_max", 0.2))
    body = _cards([("pairs compared", pairs_n.get("total")), ("pairs worth a look", pairs_n.get("kept")), ("duplicates", counts.get("duplicate", 0)),
                   ("pedigree mismatches", counts.get("pedigree_mismatch", 0)), ("related, undeclared", counts.get("related", 0)),
                   ("heterozygosity outliers", counts.get("het_outlier", 0)), ("extracts made here", len(legacy)),
                   ("members with ancestry", members_n.get("with_ancestry", 0))])
    bins = sorted(((k, float(sv.finite(v))) for k, v in hist_raw.items()), key=lambda kv: float(kv[0]) if kv[0].replace(".", "", 1).replace("-", "", 1).isdigit() else 99)
    pairs: list[dict[str, object]] = ident.get("pairs", []) or []  # type: ignore[assignment]
    flags: dict[str, list[dict[str, object]]] = ident.get("flags", {}) or {}  # type: ignore[assignment]

    def pair_flags(p: dict[str, object]) -> str:
        """The flags this pair itself raises (the same rules as identity_flags), not the union of its members' flags."""
        r, exp = p.get("relatedness"), p.get("expected")
        if not isinstance(r, (int, float)):
            return ""
        declared = isinstance(exp, (int, float)) and exp > 0
        codes = []
        if r >= dup_min:
            codes.append("duplicate")
        elif r >= rel_max and not declared:
            codes.append("related")
        if declared and not ped_min <= r <= ped_max:
            codes.append("pedigree_mismatch")
        return ", ".join(codes)

    pair_rows = [[link(p.get("sample_a")), link(p.get("sample_b")), sv.fmt(p.get("relatedness"), 3), sv.fmt(p.get("ibs0")),
                  sv.fmt(p.get("expected"), 2) if p.get("expected") is not None else "–", pair_flags(p)] for p in pairs]
    body += ('<div class="grid"><div><h3>Relatedness of every pair</h3>'
             + (sv.hist(bins, width=460, height=170, y_label="pairs") if bins else "<p class='muted'>No pairs (one member).</p>")
             + f"<p class='small muted'>somalier relatedness in bins of 0.05 (the bin's lower edge). Thresholds: related at or above {rel_max:g}, "
             f"parent-child {ped_min:g}–{ped_max:g}, duplicate at or above {dup_min:g}; siblings and parent-child both sit near 0.5, "
             f"second degree near 0.25.</p></div>"
             + "<div><h3>Pairs worth a look</h3>"
             + (big_table(["sample a", "sample b", "relatedness", "IBS0", "expected", "flags"], pair_rows,
                          what=f"{len(pair_rows)} pairs with relatedness at or above 0.1 or a pedigree expectation", open_below=40)
                if pair_rows else "<p class='muted'>No pair is related and the sheet declares no pedigree.</p>") + "</div></div>")
    flag_rows = [[link(sid), f.get("code"), Raw(f"<span class='flag-{_esc(f.get('severity'))}'>{_esc(f.get('severity'))}</span>"), f.get("message")]
                 for sid in (str(r["id"]) for r in present) for f in flags.get(sid, [])]
    body += "<h3>Identity flags</h3>" + (big_table(["sample", "flag", "severity", "message"], flag_rows, what=f"{len(flag_rows)} identity flags", open_below=40)
                                        if flag_rows else "<p class='muted'>None: no duplicate, no pedigree mismatch, no undeclared relative, no heterozygosity outlier.</p>")
    body += ("<p class='small muted'>duplicate and pedigree_mismatch fail; related and het_outlier warn. <code>ugc-wgw summary --cohort &lt;id&gt; --pass-list FILE --pass-on identity</code> "
             "writes the members without a failure, keeping the first of each undeclared related pair in sheet order (a pair the sheet declares stays).</p>")
    anc_members: list[dict[str, object]] = cq.get("ancestry_members", []) or []  # type: ignore[assignment]
    if anc_members:
        label_counts: dict[str, int] = {}
        for r in anc_members:
            label_counts[str(r.get("predicted") or "unknown")] = label_counts.get(str(r.get("predicted") or "unknown"), 0) + 1
        background = [(str(b[0]), float(b[1]), float(b[2])) for b in cq.get("ancestry_background", []) or []]  # type: ignore[union-attr]
        points = [(str(r["sample_id"]), float(r["pcs"][0]), float(r["pcs"][1]), str(r.get("predicted") or "")) for r in anc_members if len(r.get("pcs", [])) >= 2]  # type: ignore[index, arg-type]
        anc_rows = [[link(r["sample_id"]), r.get("predicted") or "–", sv.fmt(max((r.get("probs") or {}).values(), default=None), 2),  # type: ignore[union-attr]
                     ", ".join(f"{k} {sv.fmt(v, 2)}" for k, v in sorted((r.get("probs") or {}).items(), key=lambda kv: -float(kv[1] or 0)))] for r in anc_members]  # type: ignore[union-attr]
        body += ("<h3>Ancestry (projection onto the 1000 Genomes superpopulations)</h3>"
                 "<p class='small muted'>somalier ancestry projects each member onto the principal components of the 1000 Genomes samples and assigns the "
                 "nearest superpopulation (AFR, AMR, EAS, EUR, SAS): a description of genetic structure for stratifying frequencies, not an origin.</p>"
                 '<div class="grid"><div>' + sv.bars_h([(k, float(v), "") for k, v in sorted(label_counts.items())], width=300, label_w=60, cls="m3") + "</div><div>"
                 + sv.scatter(background, points, width=460, height=320) + "</div></div>"
                 + big_table(["sample", "superpopulation", "probability", "all probabilities"], anc_rows, what=f"{len(anc_rows)} members with an ancestry label", open_below=40))
    else:
        body += ("<h3>Ancestry</h3><p class='muted'>Not run: the 1000 Genomes reference set (<code>somalier_1kg_tar</code> in the reference map) is not installed, "
                 "or run_ancestry was off.</p>")
    if legacy:
        body += (f"<p class='small muted'>Extracts made in this stage from the aligned BAM for {len(legacy)} member(s) whose run predates ugc-wgw 0.14: "
                 + ", ".join(_esc(x) for x in legacy[:50]) + ("…" if len(legacy) > 50 else "") + ".</p>")
    if cq.get("qc_html"):
        body += f"<p class='small muted'>somalier's own interactive report: <code>{_esc((cq.get('qc_html') or {}).get('path'))}</code>.</p>"  # type: ignore[union-attr]
    return body


def _methyl_html(meth: dict[str, object], cme: dict[str, object], events_dist: dict[str, object] | None, thr: dict[str, float], link) -> str:  # type: ignore[no-untyped-def]
    """The Methylation section: the background and its labels, events per member, the recurrent regions and genes, the
    imprinting screen, the comparisons, deconvolution when run."""
    summ: dict[str, dict[str, dict[str, object]]] = cme.get("summary", {}) or {}  # type: ignore[assignment]
    cohort_s = summ.get("cohort", {})
    val = lambda metric: (cohort_s.get(metric) or {}).get("")  # noqa: E731
    labels: dict[str, object] = (summ.get("labels", {}).get("n_members") or {})
    notes = summ.get("notes", {})
    below = notes.get("label_below_min", {}) or {}
    body = _cards([("members", val("n_members")), ("regions", val("n_regions")), ("events (Hypo/Hyper rows)", val("n_events")),
                   ("labels", len(labels)), ("reports made in the stage", val("reports_made_in_stage")),
                   ("members flagged methyl_outlier", (meth.get("counts") or {}).get("methyl_outlier", 0))])  # type: ignore[union-attr]
    body += ("<h3>Background labels</h3>"
             + _table(["label", "members", "note"], [[lab, sv.fmt(n), ("below the minimum: " + str(below[lab])) if lab in below else ""] for lab, n in labels.items()], cls="")
             + f"<p class='small muted'>Every member is in <code>ALL</code>; its sex group and the project's <code>methyl_label_columns</code> add labels. A label needs "
               f"{sv.fmt(val('background_min_datasets'))} members before a member can be Hypo- or Hyper-methylated against it (MethBat's --min-datasets); "
               "below that every comparison is <code>InsufficientData</code>.</p>")
    if events_dist:
        body += ("<h3>Events per member</h3>"
                 + sv.strip_or_band([{"label": "methylation events", "points": sorted(events_dist["values"].items()), "outliers": events_dist.get("outliers", []),  # type: ignore[union-attr, attr-defined]
                                      "threshold": None, "beyond": [], "digits": 0, "unit": ""}], width=700)
                 + "<p class='small muted'>Rows of a member's scored profile whose compare_label is HypoMethylated, HyperMethylated, HypoASM or HyperASM, over "
                   "every label the member carries; red points lie beyond median + 3 MAD (<code>methyl_outlier</code>, warn). A member with many events "
                   "against a healthy background is the one to look at; a region flagged in many members is a cohort artefact, not many events.</p>")
    genes = summ.get("genes", {}).get("members_with_event", {}) or {}
    types = summ.get("region_type", {}).get("n_events", {}) or {}
    regions_entry: dict[str, object] = cme.get("regions", {}) or {}  # type: ignore[assignment]
    top_regions: list[dict[str, object]] = regions_entry.get("top", []) or []  # type: ignore[assignment]
    if genes or types or top_regions:
        body += '<div class="grid"><div><h3>Genes with the most members carrying an event</h3>'
        body += (sv.bars_h(sorted(((g, float(n), "") for g, n in genes.items()), key=lambda t: -t[1])[:25], width=300, label_w=120, cls="m3") if genes
                 else "<p class='muted'>No event.</p>")
        body += "</div><div><h3>Events by region type</h3>"
        body += (sv.bars_h(sorted(((t, float(n), "") for t, n in types.items()), key=lambda t: -t[1]), width=300, label_w=100, cls="m5") if types else "<p class='muted'>None.</p>")
        body += "</div></div>"
        if top_regions:
            body += ("<h3>Recurrent regions (label ALL)</h3>"
                     + big_table(["region", "position", "genes", "type", "members with an event", "hypo", "hyper", "hypo ASM", "hyper ASM", "datasets"],
                                 [[r.get("region"), f"{r.get('chrom')}:{sv.fmt(r.get('start'))}-{sv.fmt(r.get('end'))}", r.get("gene_symbol"), r.get("region_type"),
                                   sv.fmt(r.get("events")), sv.fmt(r.get("n_hypo")), sv.fmt(r.get("n_hyper")), sv.fmt(r.get("n_hypo_asm")), sv.fmt(r.get("n_hyper_asm")),
                                   sv.fmt(r.get("n_datasets"))] for r in top_regions],
                                 what=f"{len(top_regions)} regions with the most members carrying an event", open_below=30)
                     + "<p class='small muted'>A region flagged in a large share of members is a reference or mappability artefact to exclude, not a shared rare event; "
                       "the full table is <code>methyl_regions</code> and the export's <code>methyl.regions.tsv</code>.</p>")
    ev: dict[str, object] = cme.get("events", {}) or {}  # type: ignore[assignment]
    goi: list[str] = cme.get("genes_of_interest", []) or []  # type: ignore[assignment]
    if goi:
        rows: list[dict[str, object]] = ev.get("rows", []) or []  # type: ignore[assignment]
        body += (f"<h3>Events in the genes of interest</h3><p class='small muted'>methyl_genes_of_interest: {_esc(', '.join(goi))}.</p>"
                 + (big_table(["sample", "region", "genes", "type", "label", "call", "methylation z", "delta %", "ASM z"],
                              [[link(r.get("sample_id")), r.get("region"), r.get("gene_symbol"), r.get("region_type"), r.get("label"), r.get("compare"),
                                sv.fmt(r.get("zscore"), 2), sv.fmt(r.get("delta"), 1), sv.fmt(r.get("asm_zscore"), 2)] for r in rows],
                              what=f"{len(rows)} events in the genes of interest", open_below=50) if rows else "<p class='muted'>None.</p>"))
    imp_rows: list[dict[str, str]] = cme.get("imprinting", []) or []  # type: ignore[assignment]
    if imp_rows:
        per_region: dict[str, dict[str, int]] = {}
        anomalous_members: dict[str, list[str]] = {}
        for r in imp_rows:
            reg = r.get("region_label", "")
            c = per_region.setdefault(reg, {"PASS": 0, "Inconclusive": 0, "Anomalous": 0, "AnomalousQcWarning": 0})
            c[r.get("report_summary", "")] = c.get(r.get("report_summary", ""), 0) + 1
            if r.get("report_summary") in ("Anomalous", "AnomalousQcWarning"):
                anomalous_members.setdefault(reg, []).append(r.get("sample_id", ""))
        imp_summary = summ.get("imprinting", {})
        n_anom = sum(1 for sid, v in (imp_summary.get("Anomalous") or {}).items() if isinstance(v, (int, float)) and v > 0)
        body += ("<h3>Imprinting screen (methbat report)</h3>"
                 + _cards([("regions", len(per_region)), ("members with an Anomalous region", n_anom),
                           ("reports made in the stage", len(meth.get("legacy_reported", []) or []))])  # type: ignore[arg-type]
                 + _table(["region", "PASS", "inconclusive", "anomalous", "anomalous with QC warning", "members anomalous"],
                          [[reg, c["PASS"], c["Inconclusive"], c["Anomalous"], c["AnomalousQcWarning"],
                            Raw(", ".join(str(link(m)) for m in anomalous_members.get(reg, [])[:12]) + (" …" if len(anomalous_members.get(reg, [])) > 12 else ""))]
                           for reg, c in per_region.items()], cls="")
                 + "<p class='small muted'>Expected allele-specific methylation at every region; <code>Inconclusive</code> mostly means low haplotype coverage "
                   "or few phased sites. The member table is <code>methyl_imprinting</code>; each member's page lists its regions.</p>")
    comp = summ.get("compare", {})
    if comp:
        contigs: dict[str, dict[str, object]] = summ.get("compare_contigs", {}) or {}

        def where(pair: str) -> str:
            cls = {k: v for k, v in (contigs.get(pair) or {}).items() if isinstance(v, (int, float))}
            total = sum(cls.values())
            if not total:
                return ""
            return f"{sv.fmt(total)} differing region(s): " + ", ".join(f"{k} {sv.fmt(v)}" for k, v in sorted(cls.items(), key=lambda kv: -kv[1]))

        crow = [[pair, ", ".join(f"{lab} {sv.fmt(n)}" for lab, n in sorted(counts.items())),
                 "; ".join(x for x in (where(pair), str((notes.get("compare") or {}).get(pair, ""))) if x)] for pair, counts in comp.items()]
        body += "<h3>Group comparisons (methbat compare)</h3>" + _table(["pair (baseline:compare)", "summary_comparison counts", "where, note"], crow, cls="")
        body += ("<p class='small muted'>The default <code>XX:XY</code> pair is a check as much as a result: chrX regions should differ and little else. "
                 "Each comparison file is an output of the stage.</p>")
    dec = summ.get("deconvolution", {}).get("top_cell_type", {}) or {}
    if dec:
        counts_c: dict[str, int] = {}
        for cell in dec.values():
            counts_c[str(cell)] = counts_c.get(str(cell), 0) + 1
        body += ("<h3>Cell-type deconvolution (methbat deconvolve, nanomix 39Bisulfite atlas, an example)</h3>"
                 + sv.bars_h(sorted(((c, float(n), "") for c, n in counts_c.items()), key=lambda t: -t[1]), width=300, label_w=140)
                 + "<p class='small muted'>Each member's largest estimated fraction; the atlas is not validated for our tissue.</p>")
    return body


def render_cohort(d: dict[str, object], tools: list[dict[str, object]], *, project_name: str, project_url: str) -> str:
    c: dict[str, object] = d["cohort"]  # type: ignore[assignment]
    cid = str(c["id"])
    rows: list[dict[str, object]] = d.get("members", [])  # type: ignore[assignment]
    present = [r for r in rows if r.get("present")]
    dist: dict[str, dict[str, object]] = d.get("distributions", {})  # type: ignore[assignment]
    sex: dict[str, object] = d.get("sex", {})  # type: ignore[assignment]
    callsets: dict[str, dict[str, object]] = d.get("callsets", {})  # type: ignore[assignment]
    runs: dict[str, dict[str, object]] = d.get("runs", {})  # type: ignore[assignment]
    thr: dict[str, float] = d.get("thresholds", {})  # type: ignore[assignment]
    sections: list[tuple[str, str]] = []
    parts: list[str] = []

    def add(sec_id: str, title: str, body: str) -> None:
        sections.append((sec_id, title))
        parts.append(_sec(sec_id, title, body))

    def link(sid: object) -> Raw:
        return Raw(f'<a href="{_esc(sid)}.summary.html">{_esc(sid)}</a>')

    status_counts = {s: sum(1 for r in present if r.get("qc_status") == s) for s in ("pass", "info", "warn", "fail")}
    add("sec-overview", "Overview",
        f"<p class='muted'>cohort <b>{_esc(cid)}</b> · mode {_esc(d.get('mode'))} · frozen {_esc(c.get('frozen_at'))} · "
        f"sample list sha256 <code>{_esc(str(c.get('sha256'))[:12])}</code> · generated {_esc(d.get('generated_at'))}</p>"
        + _cards([("members", c.get("n")), ("with results", len(present)), ("missing", len(rows) - len(present)),
                  ("QC pass", status_counts["pass"] + status_counts["info"]), ("QC warn", status_counts["warn"]), ("QC fail", status_counts["fail"]),
                  ("sex mismatches", len(sex.get("mismatches", []))), ("cohort stages run", len(runs))])  # type: ignore[arg-type]
        + "<h3>Cohort runs summarised</h3>"
        + (_table(["stage", "run", "ugc-wgw version", "attempt", "finished"], [[s, r.get("run_id"), r.get("ugc_wgw_version"), r.get("attempt"), r.get("finished_at") or ""] for s, r in runs.items()], cls="")
           if runs else "<p class='muted'>No successful cohort stage yet.</p>")
        + (("<h3>Members without a summarised run</h3><p>" + ", ".join(_esc(x) for x in d.get("missing_members", [])) + "</p>") if d.get("missing_members") else "")  # type: ignore[union-attr]
        + (("<h3>Warnings</h3><ul class='plain small'>" + "".join(f"<li>{_esc(w)}</li>" for w in d.get("warnings", [])) + "</ul>") if d.get("warnings") else ""))  # type: ignore[union-attr]

    member_rows = [[link(r["id"]), r.get("sex_sheet") or "–", r.get("sex_inferred") or "–", (f"{r.get('sex_used')} ({r.get('sex_source')})" if r.get("sex_used") else "–"),
                    _status_badge(str(r.get("qc_status", "pass"))),
                    ", ".join(r.get("flags", [])),  # type: ignore[arg-type]
                    sv.fmt(r.get("depth_mean"), 1), sv.fmt(r.get("mapped_pct"), 2), sv.fmt(r.get("read_n50")), sv.fmt(r.get("snv")), sv.fmt(r.get("indel")),
                    sv.fmt(r.get("tstv"), 2), sv.fmt(r.get("sv_total")), sv.fmt_bp(r.get("block_ng50")), sv.fmt(r.get("trgt_genotyped")), sv.fmt_bp(r.get("roh_total")),
                    sv.fmt(r.get("kivvi_kiv2_cn"), 1), sv.fmt(r.get("kivvi_d4z4_cn"), 1)]
                   for r in present] + [[link(r["id"]), "–", "–", "–", "–", "no summarised run"] + [""] * 12 for r in rows if not r.get("present")]
    n_flagged = sum(1 for r in present if r.get("qc_status") in ("warn", "fail"))
    add("sec-members", "Members",
        big_table(["sample", "sex (sheet)", "sex (coverage)", "sex used", "QC", "flags", "depth", "mapped %", "read N50", "SNVs", "indels",
                   "Ts/Tv", "SVs", "block NG50", "TRGT genotyped", "ROH total", "KIV-2 CN", "D4Z4 CN"], member_rows,
                  what=f"{len(member_rows)} members ({n_flagged} with a QC warning or failure)")
        + "<p class='small muted'>Click a header to sort; type in the box to keep the rows that mention the text (an id, a flag, 'fail'); "
        "each sample links to its own page.</p>")

    groups = []
    outlier_rows = []
    for key, label, digits, unit, _path, _tkey, direction in MEMBER_METRICS:
        g = dist.get(key)
        if not g:
            continue
        groups.append({"label": label, "points": sorted(g["values"].items()), "outliers": g.get("outliers", []), "threshold": g.get("threshold"),  # type: ignore[union-attr]
                       "beyond": g.get("beyond_threshold", []), "digits": digits, "unit": unit})
        for sid in g.get("outliers", []):  # type: ignore[union-attr]
            outlier_rows.append([link(sid), label, sv.fmt(g["values"].get(sid), digits, unit), sv.fmt(g.get("median"), digits, unit)])  # type: ignore[union-attr]
        for sid in g.get("beyond_threshold", g.get("below_threshold", [])):  # type: ignore[union-attr]
            side = "above" if g.get("direction", direction) == "max" else "below"
            outlier_rows.append([link(sid), label, sv.fmt(g["values"].get(sid), digits, unit), f"{side} threshold {sv.fmt(g.get('threshold'), digits, unit)}"])  # type: ignore[union-attr]
    n_members = max((len(g["points"]) for g in groups), default=0)
    legend = ("One row per metric on its own scale: points are samples (hover for the id), the box is the interquartile range, "
              "the line the median; red points lie beyond median ± 3 MAD (from four members up), a dashed line marks a QC threshold."
              if n_members < sv.BAND_FROM else
              f"One row per metric on its own scale, drawn for {n_members} members: the thick line spans the 5th to 95th percentile, "
              "the box the interquartile range, the orange line the median, the small bars under the axis the members per bin; "
              "only the outliers (beyond median ± 3 MAD) and the threshold failures are drawn as red points (hover for the id); "
              "a dashed line marks a QC threshold. Every member's value is in the digest JSON.")
    add("sec-distributions", "Distributions across members",
        f"<p class='small muted'>{legend}</p>"
        + sv.strip_or_band(groups, width=700)
        + "<h3>Outliers and threshold failures</h3>"
        + (big_table(["sample", "metric", "value", "reference"], outlier_rows, what=f"{len(outlier_rows)} outliers and threshold failures", open_below=60)
           if outlier_rows else "<p class='muted'>None.</p>"))

    sheet_counts: dict[str, int] = sex.get("sheet", {})  # type: ignore[assignment]
    inf_counts: dict[str, int] = sex.get("inferred", {})  # type: ignore[assignment]
    som_counts: dict[str, int] = sex.get("somalier", {})  # type: ignore[assignment]
    used_counts: dict[str, int] = sex.get("used", {})  # type: ignore[assignment]
    src_counts: dict[str, int] = sex.get("sources", {})  # type: ignore[assignment]
    kar_counts: dict[str, int] = sex.get("karyotype", {})  # type: ignore[assignment]
    kar_of = {"MALE": "XY", "FEMALE": "XX"}
    ratio_groups = [g for g in groups if g["label"] in ("chrX / autosomes", "chrY / autosomes", "assembly: non-PAR chrX fraction", "assembly: non-PAR chrY fraction")]
    sex_keys = sorted(set(sheet_counts) | set(inf_counts) | set(som_counts) | set(used_counts))
    add("sec-sex", "Sex consistency",
        '<div class="grid"><div>' + _table(["sex", "in the sheet", "inferred from coverage", "somalier (chrX)", "used by the run"] + (["karyotype from the assembly"] if kar_counts else []),
                                          [[k, sheet_counts.get(k, 0), inf_counts.get(k, 0), som_counts.get(k, 0), used_counts.get(k, 0)]
                                           + ([kar_counts.get(kar_of.get(k, k), 0)] if kar_counts else []) for k in sex_keys]
                                          + ([["unknown karyotype", "", "", "", "", kar_counts.get("unknown", 0)]] if kar_counts.get("unknown") else []), cls="")
        + "<p class='small muted'>The run's sex (sawfish, TRGT, the frequency strata) is the sheet's when given, else the coverage inference, else somalier's "
        f"chrX call: {', '.join(f'{k} {v}' for k, v in sorted(src_counts.items())) or 'no read-based member'}.</p>"
        + f"<p>mismatches: {', '.join(_esc(x) for x in sex.get('mismatches', [])) or 'none'}; not inferred from coverage: {', '.join(_esc(x) for x in sex.get('not_inferred', [])) or 'none'}"  # type: ignore[union-attr]
        + f"; no sex from any source: {', '.join(_esc(x) for x in sex.get('unknown', [])) or 'none'}"  # type: ignore[union-attr]
        + (f"; karyotype mismatches: {', '.join(_esc(x) for x in sex.get('karyotype_mismatches', []))}" if sex.get("karyotype_mismatches") else "") + "</p>"  # type: ignore[union-attr]
        + "</div><div>" + sv.strip_or_band(ratio_groups, width=420, label_w=170) + "</div></div>")

    # identity (cohort_qc)
    ident: dict[str, object] = d.get("identity", {}) or {}  # type: ignore[assignment]
    cq: dict[str, object] = callsets.get("cohort_qc", {}) or {}
    if d.get("mode") != "assembly":
        add("sec-identity", "Identity (cohort_qc)", _identity_html(ident, cq, present, thr, link) if ident.get("available") else
            "<p class='muted'>cohort_qc has not run for this cohort (<code>ugc-wgw submit --cohort " + _esc(cid) + " --stage cohort_qc</code>).</p>")

    cc = callsets.get("cohort_call", {})
    cm = callsets.get("cohort_merge", {})
    cf = callsets.get("cohort_freq", {})
    body = ""
    if cc:
        jv: dict[str, object] = cc.get("small_variant_vcf", {})  # type: ignore[assignment]
        body += "<h3>Joint small-variant call set (cohort_call, GLnexus)</h3>" + _kv([["samples in the VCF", sv.fmt(len(jv.get("samples", [])))], ["size", _bytes(jv.get("bytes"))]])  # type: ignore[arg-type]
    if cm:
        svm: dict[str, object] = cm.get("sv", {})  # type: ignore[assignment]
        lps: dict[str, object] = cm.get("trgt_lps", {})  # type: ignore[assignment]
        supp: dict[str, object] = svm.get("supp", {})  # type: ignore[assignment]
        n_sv_samples = len(svm.get("samples", []))  # type: ignore[arg-type]
        n_records = sv.finite(svm.get("n")) or 1.0
        singletons = sv.finite(supp.get("1"))
        in_all = sv.finite(supp.get(str(n_sv_samples))) if n_sv_samples else 0.0
        body += ("<h3>Merged structural variants (cohort_merge, svx)</h3>"
                 + _cards([("records", sv.fmt(svm.get("n"))), ("samples", sv.fmt(n_sv_samples)), ("size", _bytes(svm.get("bytes"))),  # type: ignore[arg-type]
                           ("singletons", f"{sv.fmt(singletons)} ({100 * singletons / n_records:.1f} %)"),
                           ("in every member", f"{sv.fmt(in_all)} ({100 * in_all / n_records:.1f} %)")])
                 + '<div class="grid"><div><h4>By type</h4>' + sv.bars_h([(t, sv.finite(v), "") for t, v in sorted(svm.get("counts", {}).items())], width=300, label_w=60)  # type: ignore[union-attr]
                 + "</div><div><h4>Samples supporting a merged variant (SUPP)</h4>"
                 + sv.hist(supp_bins(supp, n_sv_samples), width=420, y_label="records", cumulative=True)
                 + "<p class='small muted'>Bars: records by the number of members carrying the variant, in bins that scale with the cohort; "
                 "line: the cumulative share (right axis).</p></div></div>"
                 + "<h3>Merged tandem repeats (trgt merge, trgt-lps)</h3>"
                 + _kv([["loci in the trgt-lps table", sv.fmt(lps.get("loci"))], ["samples", sv.fmt(len(lps.get("samples", [])))], ["size", _bytes(lps.get("bytes"))]]))  # type: ignore[arg-type]
        if cm.get("small_variant_vcf"):
            jv2: dict[str, object] = cm.get("small_variant_vcf", {})  # type: ignore[assignment]
            body += "<h3>Joint small-variant call set (cohort_merge, GLnexus)</h3>" + _kv([["samples in the VCF", sv.fmt(len(jv2.get("samples", [])))], ["size", _bytes(jv2.get("bytes"))]])  # type: ignore[arg-type]
    add("sec-callsets", "Cohort call sets", body or "<p class='muted'>No cohort call set summarised.</p>")

    fbody = ""
    if cf:
        summ: dict[str, dict[str, dict[str, object]]] = cf.get("summary", {})  # type: ignore[assignment]
        samples_block = summ.get("samples", {})
        fbody += _cards([("members counted", sv.fmt((samples_block.get("n") or {}).get("total"))), ("XX", sv.fmt((samples_block.get("n") or {}).get("XX"))),
                         ("XY", sv.fmt((samples_block.get("n") or {}).get("XY"))), ("unknown sex", sv.fmt((samples_block.get("n") or {}).get("unknown"))),
                         ("sex from the sheet", sv.fmt((samples_block.get("n_sex_source") or {}).get("sheet"))),
                         ("sex inferred", sv.fmt((samples_block.get("n_sex_source") or {}).get("inferred")))])
        for resource, title in (("small_variants", "small variants"), ("structural_variants", "structural variants"),
                                ("mt_small_variants", "mitochondrial small variants (mitorsaw; AC_Het heteroplasmic, AC_Hom homoplasmic members; no sex strata)"),
                                ("mt_svs", "mitochondrial SVs")):
            r = summ.get(resource)
            if r:
                fbody += _freq_resource_html(title, r)  # type: ignore[arg-type]
        fs: list[dict[str, str]] = cf.get("samples", [])  # type: ignore[assignment]
        if fs:
            has_mt = any(r.get("in_mt_vcf") not in (None, "", "none") for r in fs)
            fbody += ("<details><summary>members in the frequency VCFs</summary>"
                      + _table(["sample", "group", "sex source", "in small-variant VCF", "in SV VCF"] + (["in mt VCF"] if has_mt else []),
                               [[link(r.get("sample_id")), r.get("group"), r.get("sex_source"), r.get("in_small_variant_vcf"), r.get("in_sv_vcf")]
                                + ([r.get("in_mt_vcf")] if has_mt else []) for r in fs]) + "</details>")
    add("sec-freq", "Allele frequencies (cohort_freq)", fbody or "<p class='muted'>cohort_freq has not run.</p>")

    # methylation (cohort_methyl)
    meth: dict[str, object] = d.get("methylation", {}) or {}  # type: ignore[assignment]
    cme: dict[str, object] = callsets.get("cohort_methyl", {}) or {}
    if d.get("mode") != "assembly":
        add("sec-methyl", "Methylation (cohort_methyl)", _methyl_html(meth, cme, dist.get("methyl_events"), thr, link) if cme else
            "<p class='muted'>cohort_methyl has not run for this cohort (<code>ugc-wgw submit --cohort " + _esc(cid) + " --stage cohort_methyl</code>).</p>")

    # targeted callers across the cohort (export.targeted_summary: the same tables `ugc-wgw export` writes)
    tg: dict[str, object] = d.get("targeted", {}) or {}  # type: ignore[assignment]
    if tg:
        para_tg: dict[str, object] = tg.get("paraphase", {}) or {}  # type: ignore[assignment]
        pgx_tg: dict[str, object] = tg.get("pgx", {}) or {}  # type: ignore[assignment]
        hla_tg: dict[str, object] = tg.get("hla", {}) or {}  # type: ignore[assignment]
        mito_tg: dict[str, object] = tg.get("mito", {}) or {}  # type: ignore[assignment]
        strata_names: list[str] = tg.get("strata", []) or []  # type: ignore[assignment]
        tbody = (f"<p class='small muted'>From the members' digests, the computation `ugc-wgw export` writes as <code>cn</code>, <code>pgx</code>, "
                 f"<code>hla</code> and <code>mito</code> tables (chapter 07). Strata: {_esc(', '.join(strata_names))}"
                 + ("; ancestry is somalier's projection onto the 1000 Genomes superpopulations." if "ancestry" in strata_names
                    else "; ancestry strata appear once <code>cohort_qc</code> has run with its ancestry table.") + "</p>")
        regions: list[dict[str, object]] = [r for r in (para_tg.get("regions", []) or []) if r.get("n_called")]  # type: ignore[union-attr]
        tbody += "<h3>Copy number (Paraphase)</h3>" + _cards([
            ("regions with a call", para_tg.get("n_regions_called")), ("regions", para_tg.get("n_regions")),
            ("with a gene copy number", para_tg.get("n_regions_gene_cn")), ("members with a failed region", para_tg.get("n_members_with_failed"))])
        if regions:
            panels = []
            for reg in regions:
                counts: dict[str, int] = reg.get("cn_counts", {}) or {}  # type: ignore[assignment]
                gcounts: dict[str, int] = reg.get("gene_cn_counts", {}) or {}  # type: ignore[assignment]
                panel = (f"<div><h4>{_esc(reg.get('region'))} <span class='muted small'>{_esc(reg.get('genes'))} · {reg.get('n_called')} called"
                         f"{', ' + str(reg.get('n_failed')) + ' failed' if reg.get('n_failed') else ''}</span></h4>"
                         + sv.hist([(str(k), float(v)) for k, v in counts.items()], width=300, height=150, y_label="members"))
                if gcounts:
                    panel += "<p class='small muted'>gene copy number</p>" + sv.hist([(str(k), float(v)) for k, v in gcounts.items()], width=300, height=130, cls="m2", y_label="members")
                panels.append(panel + "</div>")
            tbody += (f"<p class='small muted'>Total copy number per member for the {len(regions)} regions with the most variation across the cohort "
                      "(then by members called); every region is in the table and in the export.</p><div class='grid'>" + "".join(panels) + "</div>")
        cn_freq: list[dict[str, object]] = para_tg.get("freq", []) or []  # type: ignore[assignment]
        all_rows = [r for r in cn_freq if r.get("stratum_kind") == "ALL" and r.get("metric") == "total_cn"]
        by_region: dict[str, list[dict[str, object]]] = {}
        for r in all_rows:
            by_region.setdefault(str(r.get("region")), []).append(r)
        region_rows = [[reg, rs[0].get("genes"), rs[0].get("n_called"), rs[0].get("n_members"),
                        ", ".join(f"{sv.fmt(r.get('value'))} × {r.get('count')}" for r in rs)] for reg, rs in sorted(by_region.items())]
        tbody += (big_table(["region", "genes", "called", "members", "copy number × members"], region_rows, what=f"{len(region_rows)} Paraphase regions with a call", open_below=40)
                  if region_rows else "<p class='muted'>No Paraphase region has a copy number in this cohort.</p>")
        genes: dict[str, dict[str, object]] = pgx_tg.get("genes", {}) or {}  # type: ignore[assignment]
        tbody += "<h3>Pharmacogenomics (StarPhase)</h3>"
        if genes:
            rates = [float(g["called_rate"]) for g in genes.values() if g.get("called_rate") is not None]
            tbody += _cards([("genes typed", len(genes)), ("median called rate", sv.fmt(statistics.median(rates), 2) if rates else "–")])
            gene_panels = []
            for gname, g in sorted(genes.items()):
                top: list = g.get("top", []) or []  # type: ignore[assignment]
                if not top:
                    continue
                gene_panels.append(f"<div><h4>{_esc(gname)} <span class='muted small'>{g.get('n_called')} of {g.get('n_members')} called</span></h4>"
                                   + sv.bars_h([(str(dip), float(n), "") for dip, n in top], width=300, label_w=130, row_h=15) + "</div>")
            tbody += "<div class='grid'>" + "".join(gene_panels) + "</div>"
            pgx_freq: list[dict[str, object]] = pgx_tg.get("freq", []) or []  # type: ignore[assignment]
            tbody += big_table(["gene", "diplotype", "stratum", "count", "fraction", "called", "members"],
                               [[r.get("gene"), r.get("simple_diplotype"), f"{r.get('stratum_kind')}:{r.get('stratum')}" if r.get("stratum_kind") != "ALL" else "ALL",
                                 r.get("count"), sv.fmt(r.get("fraction"), 3), r.get("n_called"), r.get("n_members")] for r in pgx_freq],
                               what=f"{len(pgx_freq)} diplotype frequency rows", open_below=1)
        else:
            tbody += "<p class='muted'>No StarPhase results.</p>"
        hla_freq: list[dict[str, object]] = hla_tg.get("freq", []) or []  # type: ignore[assignment]
        tbody += "<h3>HLA (StarPhase)</h3>"
        tbody += (big_table(["gene", "allele", "stratum", "count", "fraction", "alleles", "members"],
                            [[r.get("gene"), r.get("allele"), f"{r.get('stratum_kind')}:{r.get('stratum')}" if r.get("stratum_kind") != "ALL" else "ALL",
                              r.get("count"), sv.fmt(r.get("fraction"), 3), r.get("n_called"), r.get("n_members")] for r in hla_freq],
                            what=f"{len(hla_freq)} HLA allele frequency rows (both haplotypes pooled)", open_below=1)
                  if hla_freq else "<p class='muted'>No HLA alleles.</p>")
        tbody += "<h3>Mitochondria (mitorsaw)</h3>"
        unexplained: dict[str, float] = mito_tg.get("unexplained", {}) or {}  # type: ignore[assignment]
        flagged: list[str] = mito_tg.get("flagged", []) or []  # type: ignore[assignment]
        hg_counts: dict[str, int] = mito_tg.get("haplogroups", {}) or {}  # type: ignore[assignment]
        tbody += _cards([("members with haplotypes", mito_tg.get("n_members")), ("with two or more", len(mito_tg.get("multi_haplotype", []) or [])),  # type: ignore[arg-type]
                         ("minor haplotype abundance (median)", sv.fmt(mito_tg.get("minor_abundance_median"), 3)),
                         ("above the unexplained threshold", len(flagged)), ("members with a haplogroup", mito_tg.get("n_haplogrouped")),
                         ("haplogroups", len(hg_counts))])
        if hg_counts:
            macro: dict[str, int] = {}
            for hg, n in hg_counts.items():
                macro[hg[:1]] = macro.get(hg[:1], 0) + n
            tbody += ("<p class='small muted'>Haplogroup of each member's major haplotype (haplogrep3; from the member's run, or from cohort_merge for a run "
                      "before 0.15); the macro-haplogroup is its first letter. The cohort's heteroplasmic sites are in the frequency section and the export.</p>"
                      '<div class="grid"><div>' + sv.bars_h(sorted(((k, float(v), "") for k, v in macro.items()), key=lambda t: -t[1]), width=300, label_w=60, cls="m3")
                      + "</div><div>" + sv.bars_h(sorted(((k, float(v), "") for k, v in hg_counts.items()), key=lambda t: -t[1])[:30], width=300, label_w=110, cls="m5")
                      + "</div></div>")
        if unexplained:
            tbody += sv.strip_or_band([{"label": "unexplained fingerprint fraction", "points": sorted(unexplained.items()), "outliers": [],
                                        "threshold": thr.get("mito_unexplained_max"), "beyond": flagged, "digits": 3, "unit": ""}], width=700)
            if flagged:
                tbody += "<p class='small'>Above <code>mito_unexplained_max</code>: " + ", ".join(str(link(s)) for s in flagged) + "</p>"
        else:
            tbody += "<p class='muted'>No member has a fingerprint fraction (mitorsaw reports none when no heteroplasmic site is seen).</p>"
        add("sec-targeted", "Targeted callers", tbody)

    cp = callsets.get("cohort_pangenome", {})
    if cp or d.get("mode") == "assembly":
        pbody = ""
        if cp:
            graph: dict[str, object] = cp.get("graph", {})  # type: ignore[assignment]
            per_chrom: dict[str, dict[str, object]] = graph.get("chromosomes", {})  # type: ignore[assignment]
            split: dict[str, object] = cp.get("split", {})  # type: ignore[assignment]
            psumm: dict[str, dict[str, dict[str, object]]] = cp.get("summary", {})  # type: ignore[assignment]
            ps = psumm.get("samples", {})
            name_map: dict[str, dict[str, object]] = cp.get("name_map") or {}  # type: ignore[assignment]
            assigned_total = sv.finite(split.get("assigned_bp")) + sv.finite(split.get("ambiguous_bp"))
            pbody += _cards([("chromosomes built", sv.fmt(len(cp.get("chromosomes", [])))), ("members", sv.fmt((ps.get("n") or {}).get("total"))),  # type: ignore[arg-type]
                             ("in the graph VCF", sv.fmt((ps.get("n_in_vcf") or {}).get("pangenome"))),
                             ("trio-binned (phased)", sv.fmt(sum(1 for v in name_map.values() if v.get("trio"))) if name_map else "–"),
                             ("XX", sv.fmt((ps.get("n") or {}).get("XX"))), ("XY", sv.fmt((ps.get("n") or {}).get("XY"))),
                             ("karyotype unknown", sv.fmt((ps.get("n") or {}).get("unknown"))),
                             ("strata from the karyotype", sv.fmt((ps.get("n_sex_source") or {}).get("karyotype"))),
                             ("graph nodes", sv.fmt(graph.get("nodes"))), ("edges", sv.fmt(graph.get("edges"))), ("graph length", sv.fmt_bp(graph.get("length"))),
                             ("ambiguous sequence", f"{sv.fmt(sv.finite(split.get('ambiguous_bp')) / assigned_total * 100, 2)}%" if assigned_total else "–"),
                             ("graph VCF", _bytes((cp.get("pangenome_vcf") or {}).get("bytes"))), ("vcfwave VCF", _bytes((cp.get("pangenome_wave_vcf") or {}).get("bytes"))),  # type: ignore[union-attr]
                             ("frequency VCF", _bytes((cp.get("freq_vcf") or {}).get("bytes")))])  # type: ignore[union-attr]
            if per_chrom:
                has_sizes = any(g.get("ref_length") for g in per_chrom.values())
                rows_pc = []
                for c, g in per_chrom.items():
                    ref_len = sv.finite(g.get("ref_length"))
                    rows_pc.append([c, sv.fmt(g.get("nodes")), sv.fmt(g.get("edges")), sv.fmt_bp(g.get("length"))]
                                   + ([sv.fmt_bp(g.get("ref_length")), sv.fmt(sv.finite(g.get("length")) / ref_len, 3) if ref_len else "–",
                                       f"{sv.fmt_bp(g.get('hap_bp_min'))} / {sv.fmt_bp(g.get('hap_bp_avg'))} / {sv.fmt_bp(g.get('hap_bp_max'))}"] if has_sizes else [])
                                   + [sv.fmt(sv.finite(g.get("nodes")) / (ref_len / 1e6), 0) if ref_len else "–"])
                pbody += ("<h3>Per chromosome</h3>"
                          + _table(["chromosome", "nodes", "edges", "graph length"] + (["reference length", "graph / reference", "haplotype bp min / avg / max"] if has_sizes else []) + ["nodes per reference Mb"],
                                   rows_pc, cls=""))
                if has_sizes:
                    pbody += ('<div class="grid"><div><h4>Graph length over reference length</h4>'
                              + sv.bars_h([(c, sv.finite(g.get("length")) / sv.finite(g.get("ref_length")), "") for c, g in per_chrom.items() if sv.finite(g.get("ref_length"))],
                                          width=300, label_w=60, digits=3)
                              + "</div><div><h4>Shortest haplotype as a fraction of the reference</h4>"
                              + sv.bars_h([(c, sv.finite(g.get("hap_bp_min")) / sv.finite(g.get("ref_length")), "") for c, g in per_chrom.items() if sv.finite(g.get("ref_length"))],
                                          width=300, label_w=60, digits=3, vmax=1.0)
                              + "</div></div>")
            if split:
                pbody += ("<h3>Contigs binned by the split (cactus-graphmap-split)</h3>"
                          + _cards([("assigned", sv.fmt_bp(split.get("assigned_bp"))), ("ambiguous, in no graph", sv.fmt_bp(split.get("ambiguous_bp"))),
                                    ("ambiguous contigs", sv.fmt(split.get("ambiguous_contigs")))])
                          + "<details><summary>per haplotype</summary>"
                          + _table(["haplotype", "contigs assigned", "bp assigned", "contigs ambiguous", "bp ambiguous"],
                                   [[g.get("genome"), sv.fmt(g.get("assigned")), sv.fmt_bp(g.get("assigned_bp")), sv.fmt(g.get("ambiguous")), sv.fmt_bp(g.get("ambiguous_bp"))]
                                    for g in split.get("genomes", [])]) + "</details>")  # type: ignore[union-attr]
            if psumm.get("pangenome"):
                pbody += _freq_resource_html("Allele frequencies of the vcfwave VCF (missing genotypes kept: no assembled sequence)", psumm["pangenome"])  # type: ignore[arg-type]
            tv: dict[str, object] = cp.get("truvari", {})  # type: ignore[assignment]
            if tv:
                params: dict[str, object] = tv.get("params") or {}  # type: ignore[assignment]
                tp = sv.finite(tv.get("TP-base"))
                gt_ok = sv.finite(tv.get("TP-base_TP-gt"))
                pbody += ("<h3>Graph SVs against the read-based merge (truvari bench)</h3>"
                          + _cards([("recall", sv.fmt(tv.get("recall"), 4)), ("precision", sv.fmt(tv.get("precision"), 4)), ("F1", sv.fmt(tv.get("f1"), 4)),
                                    ("read-based SVs compared", sv.fmt(tv.get("base cnt"))), ("graph variants compared", sv.fmt(tv.get("comp cnt"))),
                                    ("TP (base)", sv.fmt(tv.get("TP-base"))), ("FN", sv.fmt(tv.get("FN"))), ("FP", sv.fmt(tv.get("FP"))),
                                    ("genotype concordance", sv.fmt(tv.get("gt_concordance"), 4)),
                                    ("TP with matching genotype", f"{sv.fmt(gt_ok)} ({sv.fmt(gt_ok / tp * 100, 1)}%)" if tp and tv.get("TP-base_TP-gt") is not None else "–")])
                          + "<p class='small muted'>Base: the cohort_merge SV VCF; comparison: the vcfwave VCF split to biallelic records; variants of "
                          + f"{sv.fmt(params.get('sizemin') or 50)} to {sv.fmt(params.get('sizemax') or 50000)} bp, sequence similarity not required (pctseq {sv.fmt(params.get('pctseq') if params.get('pctseq') is not None else 0, 1)}); "
                            "genotype concordance over the matched pairs.</p>")
            levels = (psumm.get("pangenome") or {}) if psumm else {}
            has_sizes = per_chrom and any(g.get("ref_length") for g in per_chrom.values())
            structure = ""
            if any(g.get("genome_bp_total") for g in per_chrom.values()):
                structure += ("<h4>Haplotype sequence in each chromosome's graph</h4>"
                              + _table(["chromosome", "haplotypes with sequence", "bp min / avg / max", "reference gaps (n, bp)"],
                                       [[c, sv.fmt(g.get("genome_haplotypes")), f"{sv.fmt_bp(g.get('genome_bp_min'))} / {sv.fmt_bp(g.get('genome_bp_avg'))} / {sv.fmt_bp(g.get('genome_bp_max'))}",
                                         f"{sv.fmt((g.get('refgaps') or {}).get('n'))}, {sv.fmt_bp((g.get('refgaps') or {}).get('bp'))}" if g.get("refgaps") else "–"]  # type: ignore[union-attr]
                                        for c, g in per_chrom.items()], cls="")
                              + "<p class='small muted'>From the join's path statistics (every haplotype's path fragments) and its reference gaps: the reference intervals no "
                                "haplotype covers in the graph (centromeres, telomeres, unassembled regions).</p>")
            degree: dict[str, object] = cp.get("degree") or {}  # type: ignore[assignment]
            lv: dict[str, object] = levels.get("records_by_level") or {}  # type: ignore[union-attr, assignment]
            vg_type: dict[str, object] = levels.get("records_by_vg_type") or {}  # type: ignore[union-attr, assignment]
            len_bins: dict[str, object] = levels.get("records_by_len_bin") or {}  # type: ignore[union-attr, assignment]
            grid = ""
            if degree:
                grid += "<div><h4>Node degree (sides)</h4>" + sv.hist([(k, sv.finite(v)) for k, v in sorted(degree.items(), key=lambda kv: (len(kv[0]), kv[0]))], width=420, y_label="node sides") + "</div>"
            if lv:
                grid += "<div><h4>Records by nesting level (vg deconstruct LV)</h4>" + sv.bars_h([(str(k), sv.finite(v), "") for k, v in sorted(lv.items())], width=300, label_w=60) + "</div>"
            if vg_type:
                grid += "<div><h4>Records by vcfwave type</h4>" + sv.bars_h([(str(k), sv.finite(v), "") for k, v in sorted(vg_type.items())], width=300, label_w=70) + "</div>"
            if len_bins:
                order_bins = ["1", "2-9", "10-49", "50-999", "1k-9.9k", ">=10k"]
                grid += "<div><h4>Allele length (vcfwave LEN)</h4>" + sv.hist([(k, sv.finite(len_bins.get(k))) for k in order_bins if k in len_bins], width=420, y_label="records") + "</div>"
            if grid:
                structure += f'<div class="grid">{grid}</div>'
                if lv or vg_type:
                    structure += ("<p class='small muted'>Level 0 records are top-level bubbles, deeper levels sit inside another record's alleles; vcfwave's TYPE and LEN are set on the "
                                  "records it realigned, the rest count as '.'.</p>")
            if structure:
                pbody += "<h3>Graph structure</h3>" + structure
            members: list[dict[str, object]] = cp.get("members") or []  # type: ignore[assignment]
            mdist: dict[str, dict[str, object]] = cp.get("member_distributions") or {}  # type: ignore[assignment]
            if members:
                groups = [{"label": g["label"], "points": sorted(g["values"].items()), "outliers": g.get("outliers", []), "threshold": g.get("threshold"),  # type: ignore[union-attr, index]
                           "digits": g.get("digits", 0), "unit": ""} for g in mdist.values()]
                low = (mdist.get("call_rate") or {}).get("beyond_threshold") or []
                pbody += ("<h3>Members' genotypes in the graph VCF</h3>"
                          + _cards([("members", sv.fmt(len(members))), ("median call rate", sv.fmt((mdist.get("call_rate") or {}).get("median"), 4)),
                                    (f"call rate below {sv.fmt(thr.get('pg_member_call_rate_min'), 2)}", sv.fmt(len(low))),
                                    ("median singletons", sv.fmt((mdist.get("n_singletons") or {}).get("median")))])
                          + sv.strip(groups, width=700)
                          + (("<p>Members below the call-rate threshold: " + ", ".join(str(link(x)) for x in low) + ".</p>") if low else "")  # type: ignore[arg-type]
                          + "<details><summary>per member</summary>"
                          + _table(["sample", "call rate", "het", "hom-alt", "hom-ref", "singletons", "missing", "sites"],
                                   [[link(m.get("sample_id")), sv.fmt(m.get("call_rate"), 4), sv.fmt(m.get("n_hets")), sv.fmt(m.get("n_nonref_hom")), sv.fmt(m.get("n_ref_hom")),
                                     sv.fmt(m.get("n_singletons")), sv.fmt(m.get("n_missing")), sv.fmt(m.get("n_sites"))] for m in members]) + "</details>"
                          + "<p class='small muted'>bcftools stats over the vcfwave VCF: a missing genotype is a site where the member has no assembled sequence, so the call "
                            "rate is the member's assembly coverage of the graph; singletons are alleles only this member carries. The threshold is a signal for the next "
                            "cohort freeze, not a sample flag.</p>")
            pm: list[dict[str, str]] = cp.get("samples", [])  # type: ignore[assignment]
            if pm:
                cs: dict[str, object] = cp.get("contig_sizes") or {}  # type: ignore[assignment]
                assigned: dict[str, object] = cs.get("assigned_bp") or {}  # type: ignore[assignment]
                ambiguous: dict[str, object] = cs.get("ambiguous_bp") or {}  # type: ignore[assignment]

                def genome_bp(sid: object, table: dict[str, object]) -> str:
                    safe = str((name_map.get(str(sid)) or {}).get("safe_name") or sid)
                    vals = [table.get(f"{safe}.{h}") for h in (1, 2)]
                    return " / ".join(sv.fmt_bp(v) for v in vals) if any(v is not None for v in vals) else "–"

                pbody += ("<details><summary>members in the graph VCF</summary>"
                          + _table(["sample", "group", "sex source", "trio-binned", "in the graph VCF", "assigned bp h1 / h2", "ambiguous bp h1 / h2"],
                                   [[link(r.get("sample_id")), r.get("group"), r.get("sex_source"), "yes" if (name_map.get(str(r.get("sample_id"))) or {}).get("trio") else "no",
                                     r.get("in_pangenome_vcf"), genome_bp(r.get("sample_id"), assigned), genome_bp(r.get("sample_id"), ambiguous)] for r in pm]) + "</details>")
            pbody += ("<p class='muted'>Genotypes of members that are not trio-binned are unphased: hifiasm's hap1/hap2 are local phase blocks, so the graph's paths "
                      "of those members are not haplotypes (no LD, imputation or PanGenie-style use without statistical phasing). Contigs the split left ambiguous are in no graph.</p>")
        add("sec-pangenome", "Pangenome graph (cohort_pangenome)", pbody or "<p class='muted'>cohort_pangenome has not run.</p>")

    asm: dict[str, dict[str, object]] = d.get("assembly", {})  # type: ignore[assignment]
    if asm:
        add("sec-assembly", "Assemblies", _assemblies_html(asm, link, open_below=50))

    versions: dict[str, list[str]] = d.get("versions", {})  # type: ignore[assignment]
    note = f" Members were summarised from ugc-wgw versions {', '.join(versions.get('ugc', []))}." if len(versions.get("ugc_wgw", [])) > 1 else ""
    add("sec-software", "Software and workflow", _software_sections(d.get("provenance", {}), tools, project_url, note))  # type: ignore[arg-type]

    nav = _nav(f"{cid} · {project_name}", sections, "<div class='sub'>members</div>" + "".join(f'<a href="{_esc(r["id"])}.summary.html">{_esc(r["id"])}</a>' for r in present[:200])
               + ("<div class='sub'>…</div>" if len(present) > 200 else ""))
    return _page(f"{cid} cohort summary", nav, f"<h1>{_esc(cid)} cohort summary</h1>" + "".join(parts))
