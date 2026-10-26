"""The cohort page of `ugc-wgw summary`: rendered from the members' per-sample digests (never from their files) plus
the cohort stages' own small outputs (merged SV VCF header and INFO, trgt-lps dimensions, cohort_freq summary).
See docs/guide/07-results.md ("Analysis summary")."""
from __future__ import annotations

import statistics
from pathlib import Path

from . import summary_parsers as sp
from . import summary_svg as sv
from .db import Cohort, RunRecord
from .layout import run_files
from .outputs import Outputs
from .report import Raw, _bytes, _esc, _table
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


def build_cohort_digest(cohort: Cohort, members: dict[str, dict[str, object]], cruns: dict[str, RunRecord], mode: str,
                        thr: dict[str, float]) -> dict[str, object]:
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
        row: dict[str, object] = {"id": sid, "present": True, "sex_sheet": sample.get("sex_sheet"), "sex_inferred": cov.get("inferred_sex"),
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
    sex = {"sheet": {}, "inferred": {}, "karyotype": {}, "mismatches": [r["id"] for r in present if "sex_mismatch" in r["flags"]],  # type: ignore[operator]
           "karyotype_mismatches": [r["id"] for r in present if "asm_karyotype_mismatch" in r["flags"]],  # type: ignore[operator]
           "not_inferred": [r["id"] for r in present if not r.get("sex_inferred") and r.get("has_reads")]}
    for r in present:
        for kind in ("sex_sheet", "sex_inferred"):
            k = str(r.get(kind) or "unknown").upper()
            sex["sheet" if kind == "sex_sheet" else "inferred"][k] = sex["sheet" if kind == "sex_sheet" else "inferred"].get(k, 0) + 1  # type: ignore[index]
        if r.get("karyotype"):
            k = str(r["karyotype"])
            sex["karyotype"][k] = sex["karyotype"].get(k, 0) + 1  # type: ignore[index]

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
            elif stage == "cohort_pangenome":
                entry.update(_pangenome_entry(outs, thr))
            elif stage == "cohort_freq":
                p = outs.opt("freq_summary")
                if isinstance(p, str) and Path(p).exists():
                    entry["summary"] = sp.parse_freq_summary(sp.read_lines(Path(p)))
                p = outs.opt("freq_samples")
                if isinstance(p, str) and Path(p).exists():
                    entry["samples"] = sp.parse_freq_samples(sp.read_lines(Path(p)))
                for name in ("small_variant_freq_vcf", "sv_freq_vcf"):
                    p = outs.opt(name)
                    if isinstance(p, str) and Path(p).exists():
                        entry[name] = {"bytes": Path(p).stat().st_size}
        except Exception as exc:  # noqa: BLE001 - one bad file must not sink the report
            warnings.append(f"{stage}: {type(exc).__name__}: {exc}")
        callsets[stage] = entry

    assembly = {sid: members[sid].get("assembly") for sid in cohort.members if sid in members and members[sid].get("assembly")}
    versions = {"ugc_wgw": sorted({str(r.get("ugc_wgw_version")) for r in present if r.get("ugc_wgw_version")}),
                "digest_schema": sorted({str(members[sid].get("schema")) for sid in members})}
    man = next((_manifest_of(cruns.get(s)) for s in ("cohort_freq", "cohort_merge", "cohort_call", "cohort_pangenome") if cruns.get(s)), {})
    prov = _provenance(man) if man else next((m.get("provenance", {}) for m in members.values()), {})
    return {"schema": DIGEST_SCHEMA, "generated_at": utc_now(), "mode": mode,
            "cohort": {"id": cohort.cohort_id, "frozen_at": cohort.frozen_at, "sha256": cohort.sample_list_sha256, "n": len(cohort.members)},
            "runs": {s: {"run_id": r.run_id, "ugc_wgw_version": r.ugc_wgw_version, "attempt": r.attempt, "finished_at": r.finished_at} for s, r in cruns.items()},
            "members": rows, "missing_members": [r["id"] for r in rows if not r.get("present")],
            "distributions": distributions, "sex": sex, "callsets": callsets, "assembly": assembly, "versions": versions,
            "thresholds": dict(thr), "provenance": prov, "warnings": warnings}


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

    member_rows = [[link(r["id"]), r.get("sex_sheet") or "–", r.get("sex_inferred") or "–", _status_badge(str(r.get("qc_status", "pass"))),
                    ", ".join(r.get("flags", [])),  # type: ignore[arg-type]
                    sv.fmt(r.get("depth_mean"), 1), sv.fmt(r.get("mapped_pct"), 2), sv.fmt(r.get("read_n50")), sv.fmt(r.get("snv")), sv.fmt(r.get("indel")),
                    sv.fmt(r.get("tstv"), 2), sv.fmt(r.get("sv_total")), sv.fmt_bp(r.get("block_ng50")), sv.fmt(r.get("trgt_genotyped")), sv.fmt_bp(r.get("roh_total")),
                    sv.fmt(r.get("kivvi_kiv2_cn"), 1), sv.fmt(r.get("kivvi_d4z4_cn"), 1)]
                   for r in present] + [[link(r["id"]), "–", "–", "–", "no summarised run"] + [""] * 12 for r in rows if not r.get("present")]
    add("sec-members", "Members",
        '<div class="wrap">' + _table(["sample", "sex (sheet)", "sex (inferred)", "QC", "flags", "depth", "mapped %", "read N50", "SNVs", "indels",
                                       "Ts/Tv", "SVs", "block NG50", "TRGT genotyped", "ROH total", "KIV-2 CN", "D4Z4 CN"], member_rows) + "</div>"
        + "<p class='small muted'>Click a header to sort; each sample links to its own page.</p>")

    groups = []
    outlier_rows = []
    for key, label, digits, unit, _path, _tkey, direction in MEMBER_METRICS:
        g = dist.get(key)
        if not g:
            continue
        groups.append({"label": label, "points": sorted(g["values"].items()), "outliers": g.get("outliers", []), "threshold": g.get("threshold"),  # type: ignore[union-attr]
                       "digits": digits, "unit": unit})
        for sid in g.get("outliers", []):  # type: ignore[union-attr]
            outlier_rows.append([link(sid), label, sv.fmt(g["values"].get(sid), digits, unit), sv.fmt(g.get("median"), digits, unit)])  # type: ignore[union-attr]
        for sid in g.get("beyond_threshold", g.get("below_threshold", [])):  # type: ignore[union-attr]
            side = "above" if g.get("direction", direction) == "max" else "below"
            outlier_rows.append([link(sid), label, sv.fmt(g["values"].get(sid), digits, unit), f"{side} threshold {sv.fmt(g.get('threshold'), digits, unit)}"])  # type: ignore[union-attr]
    add("sec-distributions", "Distributions across members",
        "<p class='small muted'>One row per metric on its own scale: points are samples (hover for the id), the box is the interquartile range, "
        "the line the median; red points lie beyond median ± 3 MAD (from four members up), a dashed line marks a QC threshold.</p>"
        + sv.strip(groups, width=700)
        + "<h3>Outliers and threshold failures</h3>"
        + (_table(["sample", "metric", "value", "reference"], outlier_rows, cls="") if outlier_rows else "<p class='muted'>None.</p>"))

    sheet_counts: dict[str, int] = sex.get("sheet", {})  # type: ignore[assignment]
    inf_counts: dict[str, int] = sex.get("inferred", {})  # type: ignore[assignment]
    kar_counts: dict[str, int] = sex.get("karyotype", {})  # type: ignore[assignment]
    kar_of = {"MALE": "XY", "FEMALE": "XX"}
    ratio_groups = [g for g in groups if g["label"] in ("chrX / autosomes", "chrY / autosomes", "assembly: non-PAR chrX fraction", "assembly: non-PAR chrY fraction")]
    sex_keys = sorted(set(sheet_counts) | set(inf_counts))
    add("sec-sex", "Sex consistency",
        '<div class="grid"><div>' + _table(["sex", "in the sheet", "inferred from coverage"] + (["karyotype from the assembly"] if kar_counts else []),
                                          [[k, sheet_counts.get(k, 0), inf_counts.get(k, 0)] + ([kar_counts.get(kar_of.get(k, k), 0)] if kar_counts else []) for k in sex_keys]
                                          + ([["unknown karyotype", "", "", kar_counts.get("unknown", 0)]] if kar_counts.get("unknown") else []), cls="")
        + f"<p>mismatches: {', '.join(_esc(x) for x in sex.get('mismatches', [])) or 'none'}; not inferred: {', '.join(_esc(x) for x in sex.get('not_inferred', [])) or 'none'}"  # type: ignore[union-attr]
        + (f"; karyotype mismatches: {', '.join(_esc(x) for x in sex.get('karyotype_mismatches', []))}" if sex.get("karyotype_mismatches") else "") + "</p>"  # type: ignore[union-attr]
        + "</div><div>" + sv.strip(ratio_groups, width=420, label_w=170) + "</div></div>")

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
        body += ("<h3>Merged structural variants (cohort_merge, svx)</h3>"
                 + _cards([("records", sv.fmt(svm.get("n"))), ("samples", sv.fmt(len(svm.get("samples", [])))), ("size", _bytes(svm.get("bytes")))])  # type: ignore[arg-type]
                 + '<div class="grid"><div><h4>By type</h4>' + sv.bars_h([(t, sv.finite(v), "") for t, v in sorted(svm.get("counts", {}).items())], width=300, label_w=60)  # type: ignore[union-attr]
                 + "</div><div><h4>Samples supporting a merged variant (SUPP)</h4>"
                 + sv.hist([(k, sv.finite(v)) for k, v in sorted(supp.items(), key=lambda kv: int(kv[0]))], width=420, y_label="records") + "</div></div>"
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
        for resource in ("small_variants", "structural_variants"):
            r = summ.get(resource)
            if r:
                fbody += _freq_resource_html(resource.replace("_", " "), r)  # type: ignore[arg-type]
        fs: list[dict[str, str]] = cf.get("samples", [])  # type: ignore[assignment]
        if fs:
            fbody += "<details><summary>members in the frequency VCFs</summary>" + _table(["sample", "group", "sex source", "in small-variant VCF", "in SV VCF"],
                                                                                          [[link(r.get("sample_id")), r.get("group"), r.get("sex_source"), r.get("in_small_variant_vcf"), r.get("in_sv_vcf")] for r in fs]) + "</details>"
    add("sec-freq", "Allele frequencies (cohort_freq)", fbody or "<p class='muted'>cohort_freq has not run.</p>")

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
