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
from .summary import DIGEST_SCHEMA, _cards, _kv, _nav, _page, _sec, _software_sections, _status_badge
from .util import UgcError, read_json, utc_now

# (key, label, digits, unit, path into the sample digest, threshold key or None)
MEMBER_METRICS: list[tuple[str, str, int, str, tuple[str, ...], str | None]] = [
    ("depth_mean", "mean depth", 1, "x", ("coverage", "mean"), "depth_mean_min"),
    ("mapped_pct", "mapped reads %", 2, "%", ("reads", "mapped_read_percent"), "mapped_read_percent_min"),
    ("read_q_median", "read quality (median)", 1, "", ("reads", "read_quality_median"), "read_quality_median_min"),
    ("read_n50", "read N50 (bp)", 0, "", ("reads", "read_length_n50"), None),
    ("snv", "SNVs", 0, "", ("small_variants", "snv"), None),
    ("indel", "indels", 0, "", ("small_variants", "indel"), None),
    ("tstv", "Ts/Tv", 2, "", ("small_variants", "tstv"), None),
    ("hethom", "het/hom", 2, "", ("small_variants", "hethom"), None),
    ("sv_total", "SVs", 0, "", ("sv", "_total"), None),
    ("phased_bp", "phased bases", 0, "", ("phasing", "phased_bp"), None),
    ("block_ng50", "phase block NG50", 0, "", ("phasing", "block_ng50"), None),
    ("trgt_genotyped", "TRGT loci genotyped", 0, "", ("tandem_repeats", "genotyped"), None),
    ("cpg", "CpG sites", 0, "", ("methylation", "cpg", "combined"), None),
    ("roh_total", "ROH total (bp)", 0, "", ("small_variants", "roh", "total_bp"), None),
    ("chrx_ratio", "chrX / autosomes", 2, "", ("coverage", "ratios", "chrX_ratio"), None),
    ("chry_ratio", "chrY / autosomes", 3, "", ("coverage", "ratios", "chrY_ratio"), None),
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
                                  "link": f"{sid}.summary.html"}
        for key, _label, _d, _u, path, _t in MEMBER_METRICS:
            row[key] = _num(_dig(m, path))
        targeted = m.get("targeted", {}) if isinstance(m.get("targeted"), dict) else {}
        kivvi = targeted.get("kivvi", {}) if isinstance(targeted.get("kivvi"), dict) else {}
        row["kivvi_kiv2_cn"] = _num((kivvi.get("kiv2") or {}).get("allele_cn"))
        row["kivvi_d4z4_cn"] = _num((kivvi.get("d4z4") or {}).get("allele_cn"))
        para = targeted.get("paraphase", {}) if isinstance(targeted.get("paraphase"), dict) else {}
        row["paraphase_failed"] = para.get("n_failed")
        rows.append(row)
    present = [r for r in rows if r.get("present")]
    distributions: dict[str, object] = {}
    for key, label, digits, unit, _path, tkey in MEMBER_METRICS:
        values = {str(r["id"]): float(r[key]) for r in present if r.get(key) is not None}  # type: ignore[arg-type]
        if not values:
            continue
        med, mad, outliers = _outliers(values)
        below = sorted(k for k, v in values.items() if tkey and v < thr.get(tkey, float("-inf")))
        distributions[key] = {"label": label, "digits": digits, "unit": unit, "values": values, "median": med, "mad": mad,
                              "outliers": outliers, "threshold": thr.get(tkey) if tkey else None, "below_threshold": below}
    sex = {"sheet": {}, "inferred": {}, "mismatches": [r["id"] for r in present if "sex_mismatch" in r["flags"]],  # type: ignore[operator]
           "not_inferred": [r["id"] for r in present if not r.get("sex_inferred")]}
    for r in present:
        for kind in ("sex_sheet", "sex_inferred"):
            k = str(r.get(kind) or "unknown").upper()
            sex["sheet" if kind == "sex_sheet" else "inferred"][k] = sex["sheet" if kind == "sex_sheet" else "inferred"].get(k, 0) + 1  # type: ignore[index]

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
    man = next((_manifest_of(cruns.get(s)) for s in ("cohort_freq", "cohort_merge", "cohort_call") if cruns.get(s)), {})
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
    for key, label, digits, unit, _path, _tkey in MEMBER_METRICS:
        g = dist.get(key)
        if not g:
            continue
        groups.append({"label": label, "points": sorted(g["values"].items()), "outliers": g.get("outliers", []), "threshold": g.get("threshold"),  # type: ignore[union-attr]
                       "digits": digits, "unit": unit})
        for sid in g.get("outliers", []):  # type: ignore[union-attr]
            outlier_rows.append([link(sid), label, sv.fmt(g["values"].get(sid), digits, unit), sv.fmt(g.get("median"), digits, unit)])  # type: ignore[union-attr]
        for sid in g.get("below_threshold", []):  # type: ignore[union-attr]
            outlier_rows.append([link(sid), label, sv.fmt(g["values"].get(sid), digits, unit), f"below threshold {sv.fmt(g.get('threshold'), digits, unit)}"])  # type: ignore[union-attr]
    add("sec-distributions", "Distributions across members",
        "<p class='small muted'>One row per metric on its own scale: points are samples (hover for the id), the box is the interquartile range, "
        "the line the median; red points lie beyond median ± 3 MAD (from four members up), a dashed line marks a QC threshold.</p>"
        + sv.strip(groups, width=700)
        + "<h3>Outliers and threshold failures</h3>"
        + (_table(["sample", "metric", "value", "reference"], outlier_rows, cls="") if outlier_rows else "<p class='muted'>None.</p>"))

    sheet_counts: dict[str, int] = sex.get("sheet", {})  # type: ignore[assignment]
    inf_counts: dict[str, int] = sex.get("inferred", {})  # type: ignore[assignment]
    ratio_groups = [g for g in groups if g["label"] in ("chrX / autosomes", "chrY / autosomes")]
    add("sec-sex", "Sex consistency",
        '<div class="grid"><div>' + _table(["sex", "in the sheet", "inferred from coverage"],
                                          [[k, sheet_counts.get(k, 0), inf_counts.get(k, 0)] for k in sorted(set(sheet_counts) | set(inf_counts))], cls="")
        + f"<p>mismatches: {', '.join(_esc(x) for x in sex.get('mismatches', [])) or 'none'}; not inferred: {', '.join(_esc(x) for x in sex.get('not_inferred', [])) or 'none'}</p>"  # type: ignore[union-attr]
        + "</div><div>" + sv.strip(ratio_groups, width=420, label_w=120) + "</div></div>")

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
            if not r:
                continue
            rec = r.get("records", {})
            fbody += (f"<h3>{_esc(resource.replace('_', ' '))}</h3>"
                      + _cards([("records", sv.fmt(rec.get("total"))), ("with AC > 0", sv.fmt(rec.get("ac_gt0"))), ("singletons", sv.fmt(rec.get("singletons"))),
                                ("mean call rate", sv.fmt((r.get("mean_call_rate") or {}).get("all"), 4)), ("call rate < 0.9", sv.fmt(rec.get("call_rate_below_0.9"))),
                                ("with AC > 0 in XX", sv.fmt((r.get("records_ac_gt0_in") or {}).get("XX"))), ("in XY", sv.fmt((r.get("records_ac_gt0_in") or {}).get("XY")))])
                      + '<div class="grid"><div><h4>Allele-frequency bins</h4>' + sv.hist([(k, sv.finite(v)) for k, v in (r.get("records_by_af_bin") or {}).items()], width=420, y_label="records")
                      + "</div><div><h4>By type</h4>" + sv.bars_h([(k, sv.finite(v), "") for k, v in sorted((r.get("records_by_type") or {}).items())], width=300, label_w=70)
                      + "</div><div><h4>By contig class</h4>" + sv.bars_h([(k, sv.finite(v), "") for k, v in sorted((r.get("records_by_contig_class") or {}).items())], width=300, label_w=70)
                      + "</div><div><h4>By filter</h4>" + sv.bars_h([(k, sv.finite(v), "") for k, v in sorted((r.get("records_by_filter") or {}).items())], width=300, label_w=110) + "</div></div>")
        fs: list[dict[str, str]] = cf.get("samples", [])  # type: ignore[assignment]
        if fs:
            fbody += "<details><summary>members in the frequency VCFs</summary>" + _table(["sample", "group", "sex source", "in small-variant VCF", "in SV VCF"],
                                                                                          [[link(r.get("sample_id")), r.get("group"), r.get("sex_source"), r.get("in_small_variant_vcf"), r.get("in_sv_vcf")] for r in fs]) + "</details>"
    add("sec-freq", "Allele frequencies (cohort_freq)", fbody or "<p class='muted'>cohort_freq has not run.</p>")

    asm: dict[str, dict[str, object]] = d.get("assembly", {})  # type: ignore[assignment]
    if asm:
        asm_rows = []
        for sid, a in asm.items():
            for h in a.get("haplotypes", []):  # type: ignore[union-attr]
                st = h.get("stats", {})
                pf = h.get("paftools", {})
                asm_rows.append([link(sid), "trio" if a.get("trio") else "single", h.get("name"), sv.fmt_bp(st.get("total_length")), sv.fmt(st.get("n_sequences")),
                                 sv.fmt_bp(st.get("n50")) if st.get("n50") else "–", sv.fmt_bp(st.get("aun")), sv.fmt(pf.get("snps")), sv.fmt(pf.get("indels"))])
        add("sec-assembly", "Assemblies", _table(["sample", "mode", "haplotype", "total length", "contigs", "N50", "auN", "SNVs vs reference", "indels vs reference"], asm_rows))

    versions: dict[str, list[str]] = d.get("versions", {})  # type: ignore[assignment]
    note = f" Members were summarised from ugc-wgw versions {', '.join(versions.get('ugc', []))}." if len(versions.get("ugc_wgw", [])) > 1 else ""
    add("sec-software", "Software and workflow", _software_sections(d.get("provenance", {}), tools, project_url, note))  # type: ignore[arg-type]

    nav = _nav(f"{cid} · {project_name}", sections, "<div class='sub'>members</div>" + "".join(f'<a href="{_esc(r["id"])}.summary.html">{_esc(r["id"])}</a>' for r in present[:200])
               + ("<div class='sub'>…</div>" if len(present) > 200 else ""))
    return _page(f"{cid} cohort summary", nav, f"<h1>{_esc(cid)} cohort summary</h1>" + "".join(parts))
