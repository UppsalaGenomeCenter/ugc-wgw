"""The cohort page of `ugc-wgw summary`: rendered from the members' per-sample digests (never from their files) plus
the cohort stages' own small outputs (merged SV VCF header and INFO, trgt-lps dimensions, cohort_freq summary).
See docs/guide/07-results.md ("Analysis summary")."""
from __future__ import annotations

import math
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
from .summary import DIGEST_SCHEMA, _cards, _clip, _has_read_data, _kv, _nav, _page, _sec, _software_sections, _status_badge, member_nav
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


# The atlas's cell types by lineage, for reading a composition at a glance (nanomix's 39 bisulfite cell types; a type the
# map does not know is `other`). The fractions themselves are never regrouped: the flag and the tables use the cell types.
CELL_GROUPS: dict[str, tuple[str, ...]] = {
    "blood": ("granulocyte", "monocyte", "T-cell", "B-cell", "NK-cell", "erythrocyte_progenitor"),
    "epithelial": ("kidney_epithelial", "thyroid_epithelial", "fallopian_epithelial", "ovary_epithelial", "epiderminal_keratinocyte",
                   "head_neck_epithelial", "lung_epithelial", "prostate_epithelial", "bladder_epithelial", "breast_luminal", "breast_basal",
                   "lung_alveolar_epithelial", "gallbladder", "gastric_epithelial", "colon_epithelial", "small_int_epithelial"),
    "mesenchymal": ("adipocyte", "endothelial_cell", "colon_fibroblasts", "heart_fibroblasts", "dermal_fibroblast", "skeletal_muscle",
                    "smooth_muscle", "heart_cardiomyocyte", "bone_osteoblast"),
    "neural": ("oligodendrocyte", "neuron"),
    "liver and pancreas": ("hepatocyte", "pancreas_duct", "pancreas_acinar", "pancreas_delta", "pancreas_beta", "pancreas_alpha"),
}
CELL_GROUP_OF = {cell: group for group, cells in CELL_GROUPS.items() for cell in cells}


def _cell_composition(parsed: dict[str, object] | None, present: dict[str, dict[str, object]], thr: dict[str, float]) -> dict[str, object]:
    """The cohort's cell-type composition from cohort_methyl's deconvolution table: per cell type the median, 5th and
    95th percentile and maximum over the members; per member its fractions (above 0.1%), its lineage sums, and its
    distance from the cohort's median composition (total variation distance: half the summed absolute differences,
    0 = the median, 1 = nothing shared). `cell_composition_outlier` (warn): a distance beyond the cohort's median +
    `cell_composition_mad_max` scaled MADs and at least `cell_composition_distance_min`: another tissue or sample
    type, a cell line, contamination, or an unusual blood picture. Cohort-relative, so a saliva cohort is judged
    against saliva."""
    members: dict[str, dict[str, object]] = {sid: m for sid, m in ((parsed or {}).get("members") or {}).items() if sid in present}  # type: ignore[union-attr]
    if not members:
        return {"available": False}
    cell_types: list[str] = list((parsed or {}).get("cell_types") or [])  # type: ignore[arg-type]
    fr = {sid: m["fractions"] for sid, m in members.items()}  # type: ignore[index]
    stats: dict[str, dict[str, float]] = {}
    for cell in cell_types:
        vals = sorted(float(f.get(cell, 0.0)) for f in fr.values())  # type: ignore[union-attr]
        n = len(vals)
        stats[cell] = {"median": statistics.median(vals), "p5": vals[max(0, math.ceil(0.05 * n) - 1)], "p95": vals[max(0, math.ceil(0.95 * n) - 1)],
                       "max": vals[-1], "n_above_5pct": sum(1 for v in vals if v >= 0.05)}
    distance = {sid: 0.5 * sum(abs(float(f.get(c, 0.0)) - stats[c]["median"]) for c in cell_types) for sid, f in fr.items()}  # type: ignore[union-attr]
    med, mad, _ = _outliers(distance)
    floor = float(thr.get("cell_composition_distance_min", 0.25))
    k = float(thr.get("cell_composition_mad_max", 4.0))
    limit = max(floor, med + k * 1.4826 * mad) if (med is not None and mad is not None and len(distance) >= 4) else None
    flags: dict[str, list[dict[str, object]]] = {}
    out_members: dict[str, dict[str, object]] = {}
    for sid, m in members.items():
        f: dict[str, float] = fr[sid]  # type: ignore[assignment]
        groups: dict[str, float] = {}
        for cell, v in f.items():
            g = CELL_GROUP_OF.get(cell, "other")
            groups[g] = groups.get(g, 0.0) + float(v)
        out_members[sid] = {"fractions": {c: round(float(v), 4) for c, v in f.items() if float(v) >= 0.001},
                            "groups": {g: round(v, 4) for g, v in sorted(groups.items(), key=lambda kv: -kv[1]) if v >= 0.001},
                            "distance": round(distance[sid], 4), "top_cell_type": m.get("top_cell_type"),
                            "residual_error": m.get("residual_error"), "raw_fraction_sum": m.get("raw_fraction_sum")}
        if limit is not None and distance[sid] > limit:
            diffs = sorted(cell_types, key=lambda c: -abs(float(f.get(c, 0.0)) - stats[c]["median"]))[:3]
            shown = ", ".join(f"{c} {float(f.get(c, 0.0)):.0%} against {stats[c]['median']:.0%}" for c in diffs)
            flags[sid] = [{"code": "cell_composition_outlier", "severity": "warn",
                           "message": f"cell-type composition far from the cohort's (distance {distance[sid]:.2f}, limit {limit:.2f}; {shown}): another "
                                      "tissue or sample type, a cell line, contamination or an unusual blood picture; read the member's methylation "
                                      "events with it", "value": round(distance[sid], 4), "threshold": round(limit, 4)}]
    order = sorted(cell_types, key=lambda c: (-stats[c]["median"], -stats[c]["max"], c))
    atlas = next((str(m.get("atlas_md5")) for m in members.values() if m.get("atlas_md5")), "")
    return {"available": True, "n_members": len(members), "cell_types": order, "stats": {c: {k2: round(v, 4) for k2, v in stats[c].items()} for c in order},
            "members": out_members, "distance": {"median": med, "mad": mad, "limit": limit}, "flags": flags, "atlas_md5": atlas,
            "group_of": {c: CELL_GROUP_OF.get(c, "other") for c in order}}


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


def growth_from_hist(hist: dict[int, float], n: int) -> tuple[list[float], list[float]]:
    """From `hist[c]` = the amount carried by exactly c of n members: for k = 1..n, the expected amount carried by at
    least one of k members drawn at random (the growth curve) and by all k of them (the core). An amount with c
    carriers is missed by a draw of k with probability C(n-c, k) / C(n, k) and is in all of them with C(c, k) / C(n, k)."""
    lg = [math.lgamma(i + 1) for i in range(n + 1)]

    def log_choose(a: int, b: int) -> float:
        return lg[a] - lg[b] - lg[a - b]

    levels = [(min(c, n), v) for c, v in hist.items() if c > 0 and v]
    union, core = [], []
    for k in range(1, n + 1):
        u = co = 0.0
        denom = log_choose(n, k)
        for c, v in levels:
            u += v if k > n - c else v * (1.0 - math.exp(log_choose(n - c, k) - denom))
            if c >= k:
                co += v * math.exp(log_choose(c, k) - denom)
        union.append(u)
        core.append(co)
    return union, core


GROWTH_CHECKPOINTS = (1, 2, 5, 10, 20, 50, 100, 200, 400, 800)
MQTL_TOP_ROWS = 500   # rows of cohort_mqtl's top table kept in the digest and on the page; the file has them all


def _growth_html(growth: dict[str, object]) -> str:
    """The cohort page's growth block: the sequence the reference lacks against the number of members (expected for
    members at random, and in the cohort's own order), everything the members carry with its core, and the same at a
    few sizes as a table with what one more member adds there."""
    n = int(growth.get("members") or 0)
    nonref: list[float] = growth.get("nonref_bp") or []  # type: ignore[assignment]
    ordered: list[float] = growth.get("ordered_nonref_bp") or []  # type: ignore[assignment]
    allbp: list[float] = growth.get("all_bp") or []  # type: ignore[assignment]
    core: list[float] = growth.get("core_bp") or []  # type: ignore[assignment]
    nodes: list[float] = growth.get("nonref_node") or []  # type: ignore[assignment]
    series = [("members at random (expected)", "ln2", [(i + 1, v) for i, v in enumerate(nonref)])]
    if ordered and len(ordered) == len(nonref):
        series.append(("in the cohort's order", "ln3", [(i + 1, v) for i, v in enumerate(ordered)]))
    body = ("<h3>Growth: what each further member adds</h3>"
            + sv.lines(series, x_label="members in the graph", y_label="sequence the reference lacks", y_format=sv.fmt_bp))
    if allbp and core:
        body += sv.lines([("carried by at least one", "ln2", [(i + 1, v) for i, v in enumerate(allbp)]),
                          ("carried by all (core)", "ln1", [(i + 1, v) for i, v in enumerate(core)])],
                         x_label="members in the graph", y_label="sequence the members carry", y_format=sv.fmt_bp)
    marks = [k for k in GROWTH_CHECKPOINTS if k < n] + [n]
    rows = []
    for k in marks:
        gain = (nonref[k - 1] - nonref[k - 2]) if k >= 2 else nonref[0]
        rows.append([k, sv.fmt_bp(nonref[k - 1]), sv.fmt_bp(gain), f"{nonref[k - 1] / nonref[-1] * 100:.1f}%" if nonref[-1] else "–",
                     sv.fmt(nodes[k - 1]) if len(nodes) >= k else "–", sv.fmt_bp(core[k - 1]) if len(core) >= k else "–"])
    body += _table(["members", "non-reference sequence", "added by the last one", "share of the cohort's", "non-reference nodes", "core"], rows, cls="")
    body += ("<p class='small muted'>From panacus's coverage histograms of the chromosome graphs, the members' paths grouped per sample and the "
             "reference's path left out: the expected amount for that many members drawn at random from the cohort. Non-reference means "
             "not on the reference's path; the core is what every one of them carries. A curve that still climbs at the right edge "
             "says the next member brings new sequence; satellite arrays and other poorly aligned regions keep it from ever going flat.</p>")
    return body


def pangenome_growth(tables: list[dict[str, object]], members: list[str]) -> dict[str, object]:
    """The chromosomes' growth tables as genome-wide curves over the cohort's members: the histograms add up because
    every chromosome is counted against the whole cohort (a member without sequence in one carries nothing there).
    `nonref_bp` and `all_bp` are expected values for k members drawn at random, `core_bp` what all k share,
    `ordered_nonref_bp` the cumulative amount in the cohort's own order (the ranking's, when frozen from one)."""
    n = max([int(t.get("cohort") or 0) for t in tables] + [0]) or len(members)
    if n <= 0 or not tables:
        return {}
    summed: dict[tuple[str, str], dict[int, float]] = {}
    for t in tables:
        for key, h in (t.get("hist") or {}).items():  # type: ignore[union-attr]
            acc = summed.setdefault(key, {})
            for c, v in h.items():
                acc[c] = acc.get(c, 0.0) + v
    out: dict[str, object] = {"members": n, "chromosomes": len(tables)}
    for name, key in (("nonref_bp", ("nonref", "bp")), ("nonref_node", ("nonref", "node")), ("all_bp", ("all", "bp"))):
        if key in summed:
            union, core = growth_from_hist(summed[key], n)
            out[name] = [round(v) for v in union]
            if name == "all_bp":
                out["core_bp"] = [round(v) for v in core]
    order = [m for m in members] or []
    pos = {m: i for i, m in enumerate(order)}
    for name, which in (("ordered_nonref_bp", "nonref"), ("ordered_all_bp", "all")):
        total = [0.0] * len(order)
        seen = False
        for t in tables:
            rows = [(pos[m], v) for m, v in (t.get("ordered") or {}).get(which, []) if m in pos]  # type: ignore[union-attr]
            if not rows:
                continue
            seen = True
            rows.sort()
            cur, j = 0.0, 0
            for i in range(len(order)):       # a chromosome without this member keeps its last value
                while j < len(rows) and rows[j][0] <= i:
                    cur = rows[j][1]
                    j += 1
                total[i] += cur
        if seen:
            out[name] = [round(v) for v in total]
    return out


def _mqtl_html(entry: dict[str, object]) -> str:
    """The cohort page's cohort_mqtl section: what was tested, how many pairs stand at the false-discovery rate, the
    inflation factor with what moves it, and the top pairs."""
    s: dict[str, dict[str, object]] = entry.get("summary") or {}  # type: ignore[assignment]

    def g(metric: str, key: str) -> object:
        return (s.get(metric) or {}).get(key)

    islands: dict[str, object] = s.get("islands") or {}
    variants: dict[str, object] = s.get("variants_tested") or {}
    notes = list((s.get("note") or {}).keys())
    fdr = g("parameters", "fdr")
    lam = g("inflation", "lambda")
    body = ("<p class='small muted'>A research screen, not a diagnostic result: does carrying a structural variant, or a longer tandem repeat, go "
            "with a different methylation of a CpG island within "
            f"{_esc(sv.fmt_bp(g('parameters', 'window_bp')))}? One tissue, CpG islands only, autosomes only; an association is not a mechanism.</p>")
    body += _cards([("members in the model", sv.fmt(g("members", "used"))), ("islands tested", sv.fmt(islands.get("tested"))),
                    ("variants tested", ", ".join(f"{sv.fmt(v)} {k}" for k, v in sorted(variants.items())) or "0"),
                    ("pairs tested", sv.fmt(g("pairs", "total"))),
                    (f"pairs at q ≤ {fdr}", sv.fmt(g("hits", "total"))), ("islands with one", sv.fmt(g("hits", "islands"))),
                    ("inflation factor", sv.fmt(lam, 2) if lam is not None else "–")])
    for note in notes:
        body += f"<p class='warn'>{_esc(note)}.</p>"
    covs = [k for k in (s.get("covariate") or {}) if k != "intercept"]
    body += ("<p class='small'>Covariates in the model: " + (_esc(", ".join(covs)) if covs else "none besides the intercept")
             + ". Sex is there when both are present; the ancestry components and the blood cell fractions only when every member has them.</p>")
    left = {k: v for k, v in islands.items() if k.startswith("not tested")}
    if left:
        body += ("<details><summary>islands not tested</summary>"
                 + _table(["why", "islands"], [[k.replace("not tested: ", ""), sv.fmt(v)] for k, v in sorted(left.items())], cls="") + "</details>")
    top: list[dict[str, str]] = entry.get("top") or []  # type: ignore[assignment]
    if top:
        rows = [[f"{r.get('chrom')}:{r.get('start')}-{r.get('end')}", r.get("gene_symbol") or "–", r.get("region_type") or "–", r.get("variant_type"),
                 _clip(str(r.get("variant_id")), 40), sv.fmt(sp._num(r.get("variant_len"))), sv.fmt(sp._num(r.get("distance"))), r.get("n_minor"),
                 f"{r.get('effect')} {str(r.get('effect_unit') or '').replace('points ', '')}", r.get("p"), r.get("q")] for r in top]
        n_top = int(entry.get("n_top") or len(top))
        body += (f"<h3>Pairs at q ≤ {_esc(fdr)}" + (f" (the first {len(top)} of {n_top}, by p-value)" if n_top > len(top) else "") + "</h3>"
                 + big_table(["island", "gene", "region", "variant", "id", "length (bp)", "distance (bp)", "members off the commonest genotype",
                              "effect (methylation points)", "p", "q"], rows, what="pairs"))  # type: ignore[arg-type]
    elif not notes:
        body += "<p class='muted'>No pair at that false-discovery rate.</p>"
    body += ("<p class='small muted'>Methylation is rank-normalised across the members before the test; the effect is the slope on the "
             "original scale. The inflation factor is the median test statistic over its expectation: about 1 when most pairs are null, "
             "above it when relatives, population structure or a batch run through the cohort (the exclude_ids input leaves members out). "
             "Every tested pair is in the pairs table of the export.</p>")
    return body


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
    paths = outs.opt("chrom_growth")      # 0.22.0: panacus histograms and the ordered growth per chromosome
    if isinstance(paths, list):
        tables = [sp.parse_growth(sp.read_lines(Path(q))) for q in paths if isinstance(q, str) and Path(q).exists()]
        growth = pangenome_growth(tables, list((entry.get("name_map") or {}).keys()))  # type: ignore[union-attr]
        if growth:
            entry["growth"] = growth
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
    # a member related to many others is the anomaly, not its partners (summary_parsers.related_hubs)
    hubs = sp.related_hubs(pairs, lo=rel_max, hi=dup_min, hub_min=hub_min)
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
                        thr: dict[str, float], *, genes_of_interest: list[str] | None = None,
                        member_tools: dict[str, str] | None = None) -> dict[str, object]:
    """`member_tools`: member -> "aligner + caller" (summary.member_tools), joined here so the digests stay as they are.
    More than one set among the members is the cohort warning `mixed_tools`: nothing is refused, the page says it."""
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
                                  "link": f"{sid}.summary.html", "has_reads": _has_read_data(m),
                                  "tools": (member_tools or {}).get(sid)}
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
                p = outs.opt("methyl_deconvolution")   # 0.20.0; absent on earlier runs
                if isinstance(p, str) and Path(p).exists():
                    entry["deconvolution"] = sp.parse_methyl_deconvolution(sp.read_lines(Path(p)))
                p = outs.opt("methyl_events")
                if isinstance(p, str) and Path(p).exists():
                    entry["events"] = sp.parse_methyl_events(sp.read_lines(Path(p)), genes=genes_of_interest or ())
                    entry["events"]["bytes"] = Path(p).stat().st_size  # type: ignore[index]
                for name in ("cohort_background", "methyl_collection", "joint_segments"):
                    p = outs.opt(name)
                    if isinstance(p, str) and Path(p).exists():
                        entry[name] = {"path": p, "bytes": Path(p).stat().st_size}
                entry["genes_of_interest"] = list(genes_of_interest or [])
            elif stage == "cohort_tr":
                for name, parser, key in (("tr_summary", sp.parse_tr_summary, "summary"), ("tr_samples", sp.parse_tr_samples, "samples"),
                                          ("tr_pathogenic_screen", sp.parse_tr_screen, "screen"), ("tr_screen_loci", sp.parse_tr_screen_loci, "screen_loci"),
                                          ("tr_outliers_top", sp.parse_tr_outliers_top, "outliers_top")):
                    p = outs.opt(name)
                    if isinstance(p, str) and Path(p).exists():
                        entry[key] = parser(sp.read_lines(Path(p)))  # type: ignore[operator]
                for name in ("tr_locus_summary", "tr_allele_spectrum", "tr_outlier_pairs"):   # the big tables: sizes only
                    p = outs.opt(name)
                    if isinstance(p, str) and Path(p).exists():
                        entry[name] = {"bytes": Path(p).stat().st_size}
            elif stage == "cohort_mqtl":
                p = outs.opt("mqtl_summary")
                if isinstance(p, str) and Path(p).exists():
                    entry["summary"] = sp.parse_tr_summary(sp.read_lines(Path(p))).get("mqtl", {})   # the same long format
                p = outs.opt("mqtl_top")
                if isinstance(p, str) and Path(p).exists():
                    _cols, rows_top = sp.table(sp.read_lines(Path(p)))
                    entry["n_top"] = len(rows_top)
                    entry["top"] = rows_top[:MQTL_TOP_ROWS]
                for name in ("mqtl_matrix", "mqtl_pairs"):
                    p = outs.opt(name)
                    if isinstance(p, str) and Path(p).exists():
                        entry[name] = {"bytes": Path(p).stat().st_size}
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
        deconv = _cell_composition(me_entry.pop("deconvolution", None), present_members, thr)  # type: ignore[arg-type]
        for sid, fl in (deconv.get("flags") or {}).items():  # type: ignore[union-attr]
            mflags.setdefault(sid, []).extend(fl)
        for r in present:
            fl = (deconv.get("flags") or {}).get(str(r["id"]), [])  # type: ignore[union-attr]
            if fl:
                r["flags"] = sorted(set(r["flags"]) | {str(f["code"]) for f in fl})  # type: ignore[arg-type]
                r["qc_status"] = _worst([str(r.get("qc_status", "pass"))] + [str(f["severity"]) for f in fl])
        msumm: dict[str, dict[str, object]] = me_entry.get("summary") if isinstance(me_entry.get("summary"), dict) else {}  # type: ignore[assignment]
        methylation.update({"flags": mflags, "counts": {"methyl_outlier": sum(1 for fl in mflags.values() if any(f["code"] == "methyl_outlier" for f in fl)),
                                                        "cell_composition_outlier": len(deconv.get("flags") or {})},  # type: ignore[arg-type]
                            "deconvolution": deconv,
                            "members": {sid: by_id_m.get(sid, {"sample_id": sid}) for sid in cohort.members},
                            "legacy_reported": [sid for sid, r in by_id_m.items() if str(r.get("report_source")) == "stage"],
                            # each comparison's summary_comparison counts and where its differing regions lie (chrX, chrY, autosomal)
                            "comparisons": {pair: {"counts": counts, "contigs": (msumm.get("compare_contigs") or {}).get(pair, {})}
                                            for pair, counts in (msumm.get("compare") or {}).items()}})

    # tandem repeats: cohort_tr's samples table joined to the member rows (outliers per member, the MAD outliers), the
    # screen's findings with each member's own dropout class at the locus (from its digest: the stage reads no per-sample file)
    tr_entry: dict[str, object] = callsets.get("cohort_tr") if isinstance(callsets.get("cohort_tr"), dict) else {}  # type: ignore[assignment]
    tandem: dict[str, object] = {"available": bool(tr_entry.get("samples") is not None)}
    if tandem["available"]:
        trows: list[dict[str, object]] = tr_entry.get("samples", [])  # type: ignore[assignment]
        by_id_t = {str(r.get("sample_id")): r for r in trows}
        tvalues = {sid: float(r["n_outliers"]) for sid, r in by_id_t.items() if isinstance(r.get("n_outliers"), (int, float)) and sid in present_members}
        tflags: dict[str, list[dict[str, object]]] = {}
        tmed, tmad, _ = _outliers(tvalues) if tvalues else (None, None, [])
        tmax = thr.get("tr_outlier_count_mad_max", 3.0)
        # scaled MADs (1.4826 x MAD, a standard deviation for a normal count), as the stage's own locus rule: the plain MAD
        # flagged 62 of 864 members (7%) on the first production cohort, the scaled one 29
        tlimit = (tmed + tmax * 1.4826 * tmad) if (tmed is not None and tmad is not None and len(tvalues) >= 4) else None
        for sid, v in tvalues.items():
            if tlimit is not None and v > tlimit and tmad:   # a constant cohort (mad 0) flags nobody
                tflags.setdefault(sid, []).append({"code": "tr_outlier_count_high", "severity": "warn",
                                                   "message": f"{int(v)} expansion outliers at ordinary tandem-repeat loci, above the cohort's median {sv.fmt(tmed)} + "
                                                              f"{sv.fmt(tmax, 1)} scaled MADs ({sv.fmt(tlimit, 1)}): an ancestry that is rare in this cohort, a lower "
                                                              "depth than the others (the count rises as depth falls), or a sample problem",
                                                   "value": v, "threshold": tlimit})
        for r in present:
            fl = tflags.get(str(r["id"]), [])
            if fl:
                r["flags"] = sorted(set(r["flags"]) | {str(f["code"]) for f in fl})  # type: ignore[arg-type]
                r["qc_status"] = _worst([str(r.get("qc_status", "pass"))] + [str(f["severity"]) for f in fl])
            tr_row = by_id_t.get(str(r["id"]))
            r["tr_outliers"] = _num(tr_row.get("n_outliers")) if tr_row else None
        if tvalues:
            distributions["tr_outliers"] = {"label": "expansion outliers at ordinary TR loci (cohort_tr)", "digits": 0, "unit": "", "values": tvalues, "median": tmed,
                                            "mad": tmad, "outliers": sorted(tflags), "threshold": None, "direction": "max", "beyond_threshold": [],
                                            "below_threshold": []}
        screen: dict[str, object] = tr_entry.get("screen") if isinstance(tr_entry.get("screen"), dict) else {}  # type: ignore[assignment]
        findings: list[dict[str, object]] = [dict(f) for f in screen.get("findings", [])]  # type: ignore[union-attr]
        for f in findings:
            loci_of: list[dict[str, object]] = ((present_members.get(str(f.get("sample_id"))) or {}).get("tandem_repeats") or {}).get("disease_loci") or []  # type: ignore[union-attr, assignment]
            f["dropout"] = next((str(l.get("dropout") or "") for l in loci_of if l.get("trid") == f.get("trid")), "")
        tandem.update({"flags": tflags, "counts": {"tr_outlier_count_high": len(tflags)},
                       "members": {sid: by_id_t.get(sid, {"sample_id": sid}) for sid in cohort.members},
                       "findings": findings, "by_class": screen.get("by_class", {}),
                       "outliers_top": tr_entry.get("outliers_top", {})})
        for dup in ("screen", "outliers_top", "samples"):
            tr_entry.pop(dup, None)

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
    # the aligner and caller sets among the members: one set is the normal case and says nothing
    tool_groups: dict[str, int] = {}
    for r in present:
        if r.get("tools"):
            tool_groups[str(r["tools"])] = tool_groups.get(str(r["tools"]), 0) + 1
    if len(tool_groups) > 1:
        shown = "; ".join(f"{k}: {n}" for k, n in sorted(tool_groups.items(), key=lambda kv: (-kv[1], kv[0])))
        warnings.append(f"mixed_tools: the members were made with different aligners or small-variant callers ({shown}); joint "
                        "calls, frequencies and outlier statistics can show that difference as a batch effect")

    return {"schema": DIGEST_SCHEMA, "generated_at": utc_now(), "mode": mode,
            "cohort": {"id": cohort.cohort_id, "frozen_at": cohort.frozen_at, "sha256": cohort.sample_list_sha256, "n": len(cohort.members)},
            "runs": {s: {"run_id": r.run_id, "ugc_wgw_version": r.ugc_wgw_version, "attempt": r.attempt, "finished_at": r.finished_at} for s, r in cruns.items()},
            "members": rows, "missing_members": [r["id"] for r in rows if not r.get("present")],
            "distributions": distributions, "sex": sex, "identity": identity, "methylation": methylation, "tandem_repeats": tandem, "callsets": callsets,
            "targeted": targeted,
            "assembly": assembly,
            "versions": versions,
            "tools": {"groups": tool_groups, "mixed": len(tool_groups) > 1},
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


def _tr_html(trd: dict[str, object], ctr: dict[str, object], outlier_dist: dict[str, object] | None, thr: dict[str, float], link) -> str:  # type: ignore[no-untyped-def]
    """The Tandem repeats section: the screen of the named disease loci (rows worth a look, the per-locus table), the
    expansion outliers per member and the loci with the most carriers, the resource files, STRchive's attribution."""
    summ: dict[str, dict[str, dict[str, object]]] = ctr.get("summary", {}) or {}  # type: ignore[assignment]
    cohort_s = summ.get("cohort", {})
    val = lambda metric: (cohort_s.get(metric) or {}).get("")  # noqa: E731
    status: dict[str, object] = (summ.get("loci", {}).get("outlier_status") or {})
    params = {k: (v or {}).get("") for k, v in (summ.get("parameters") or {}).items()}
    findings: list[dict[str, object]] = trd.get("findings", []) or []  # type: ignore[assignment]
    by_class: dict[str, int] = trd.get("by_class", {}) or {}  # type: ignore[assignment]
    classed = [f for f in findings if f.get("class") in ("pathogenic", "intermediate")]
    body = _cards([("members", val("n_members")), ("loci in the catalog", val("n_loci")), ("named disease loci", (summ.get("screen", {}).get("n_loci") or {}).get("")),
                   ("pathogenic-class rows", by_class.get("pathogenic", 0)), ("intermediate-class rows", by_class.get("intermediate", 0)),
                   ("expansion-outlier pairs", val("n_outlier_pairs")),
                   ("loci with a listed outlier", status.get("listed")), ("loci with a common long allele", status.get("common")),
                   ("members flagged tr_outlier_count_high", (trd.get("counts") or {}).get("tr_outlier_count_high", 0))])  # type: ignore[union-attr]
    body += ("<h3>Disease-locus screen (STRchive thresholds)</h3>"
             "<p class='small muted'><b>A research screen, not a diagnostic result.</b> Per member and named locus, each allele's copies of the pathogenic "
             "motif (TRGT's motif counts; allele length over motif length where the table says so) are classed against STRchive's ranges: "
             "<code>pathogenic</code>, <code>intermediate</code>, <code>above_benign</code> and <code>below_benign</code> for the gaps around the benign "
             "range, <code>below_pathogenic</code> where the locus has no benign maximum, <code>not_assessable</code> where the thresholds do not count "
             "the repeat unit. TRGT's counts tolerate interruptions and its frames are wider than the clinical ones at some loci, so a class is a "
             "reason to look at the reads (<code>trgt plot</code> on the member's spanning reads), not a finding. Fit: <code>consistent</code> with the "
             "locus's inheritance, or <code>carrier</code> (one allele at a recessive locus).</p>")
    if classed or findings:
        rows = [[link(str(f.get("sample_id"))), f.get("gene"), _clip(str(f.get("disease") or ""), 40), f.get("inheritance"), Raw(f"<b>{_esc(f.get('class'))}</b>"),
                 f.get("mode_match"), f.get("gt"), f"{sv.fmt(f.get('a1_copies'))} / {sv.fmt(f.get('a2_copies'))}" if f.get("a2_copies") is not None else sv.fmt(f.get("a1_copies")),
                 f"{sv.fmt(f.get('a1_al'))} / {sv.fmt(f.get('a2_al'))}" if f.get("a2_al") is not None else sv.fmt(f.get("a1_al")), sv.fmt(f.get("a1_lps")),
                 f.get("a1_am") or "", f.get("borderline") or "", f.get("dropout") or "", _clip(str(f.get("notes") or ""), 60)] for f in findings]
        body += big_table(["member", "gene", "disease", "inheritance", "class", "fit", "GT", "copies (a1 / a2)", "AL bp (a1 / a2)", "LPS a1", "5mC a1", "borderline",
                           "dropout", "notes"], rows, what=f"{len(rows)} member x locus rows worth a look ({len(classed)} pathogenic or intermediate)", open_below=40)
        body += ("<p class='small muted'>a1 is the allele with more copies. LPS: the pathogenic motif's longest uninterrupted run (shorter than the copies means "
                 "interruptions). 5mC: the allele's mean methylation; read it within a sex group at X-linked loci. Dropout: the member's own class at the locus "
                 "from its run. Listed: the pathogenic and intermediate rows, and any row whose call is in doubt (the allele's sequence and length disagree, "
                 "the ploidy does not fit the sex, the catalog's motifs changed). The other classes are counted below and per locus; every row is in the "
                 "stage's screen table.</p>")
    else:
        body += "<p class='muted'>No member has a pathogenic or intermediate class, or a call in doubt, at a called locus.</p>"
    if by_class:
        body += "<p class='small'>Member x locus rows by class: " + ", ".join(f"{_esc(k)} {sv.fmt(v)}" for k, v in sorted(by_class.items(), key=lambda kv: -kv[1])) + ".</p>"
    loci: list[dict[str, object]] = ctr.get("screen_loci", []) or []  # type: ignore[assignment]
    if loci:
        lrows = [[l.get("gene"), l.get("trid"), l.get("inheritance"), f"{l.get('chrom')}:{sv.fmt(l.get('start'))}", _clip(str(l.get("pathogenic_motifs") or ""), 28),
                  f"{sv.fmt(l.get('benign_max'))} / {sv.fmt(l.get('pathogenic_min'))}", f"{sv.fmt(l.get('n_called'))} of {sv.fmt(l.get('n_expected'))}",
                  sv.fmt(l.get("copies_median"), 1), sv.fmt(l.get("copies_max")), sv.fmt(l.get("n_below_benign")), sv.fmt(l.get("n_above_benign")),
                  sv.fmt(l.get("n_intermediate")), sv.fmt(l.get("n_pathogenic")), sv.fmt(l.get("n_expanded_other_motif")), l.get("call_note")] for l in loci]
        body += ("<h3>Named loci</h3>"
                 + big_table(["gene", "locus", "inheritance", "position", "pathogenic motif(s)", "benign max / pathogenic min", "members called", "copies median",
                              "copies max", "below benign", "above benign", "intermediate", "pathogenic", "other motif", "note"], lrows,
                             what=f"{len(lrows)} named loci ({sum(1 for l in loci if l.get('call_note') != 'ok')} with a call note)", open_below=10)
                 + "<p class='small muted'>Members per class at each locus. <code>mostly_uncalled</code>: a catalog or coverage problem at the locus, not a "
                   "property of the members; <code>few_uncalled</code>: the uncalled members are the ones to check for an allele longer than the reads; "
                   "<code>threshold_or_frame_suspect</code>: more than a fifth of the alleles are above or below the benign range, so the thresholds and "
                   "TRGT's frame disagree here and the classes at this locus say little; <code>not_assessable</code>: the thresholds do not count this "
                   "repeat unit (ADTKD_MUC1: one inserted base inside a unit), so nothing is classed.</p>")
    body += "<h3>Expansion outliers</h3>"
    if outlier_dist:
        body += sv.strip_or_band([{"label": "outliers per member", "points": sorted(outlier_dist["values"].items()), "outliers": outlier_dist.get("outliers", []),  # type: ignore[union-attr, attr-defined]
                                   "threshold": None, "beyond": [], "digits": 0, "unit": ""}], width=700)
    body += ("<p class='small muted'>A member and locus pair is listed when the member's longer allele is at least "
             f"{_esc(sv.fmt(params.get('outlier_min_bp')))} bp and {_esc(sv.fmt(params.get('outlier_mad_min'), 1))} scaled MADs above the cohort's median (the scale "
             "never below one repeat unit), at loci with at least " + _esc(sv.fmt(params.get("outlier_min_called"))) + " called members, and only while "
             "the carriers are few: a locus where many members qualify has a common long allele and lists nobody. The count per member is over the ordinary "
             "loci; red points lie beyond median + " + _esc(sv.fmt(thr.get("tr_outlier_count_mad_max", 3.0), 1)) + " scaled MADs "
             "(<code>tr_outlier_count_high</code>, warn). A high count goes with an ancestry that is rare in the cohort and with a lower depth than the "
             "others, more often than with a sample problem; contractions are not scanned.</p>")
    top_loci = (summ.get("outlier_loci", {}).get("n_carriers") or {})
    named_of = (summ.get("outlier_loci", {}).get("named_locus") or {})
    if top_loci:
        trows = [[trid, sv.fmt(n), named_of.get(trid, "")] for trid, n in sorted(top_loci.items(), key=lambda kv: (-float(kv[1]), kv[0]))[:25]]  # type: ignore[arg-type]
        body += "<h3>Loci with the most outlier members</h3>" + big_table(["locus", "members listed", "named locus"], trows, what=f"{len(trows)} loci", open_below=12)
    body += _kv([["allele-length resource (trid, al, n_alleles, n_members)", _bytes((ctr.get("tr_allele_spectrum") or {}).get("bytes"))],  # type: ignore[union-attr]
                 ["per-locus table (median, MAD, quantiles, status)", _bytes((ctr.get("tr_locus_summary") or {}).get("bytes"))],  # type: ignore[union-attr]
                 ["every outlier pair", _bytes((ctr.get("tr_outlier_pairs") or {}).get("bytes"))]])  # type: ignore[union-attr]
    source = ((summ.get("source", {}).get("tr_map_source") or {}).get("values") or [])
    if source:
        body += "<p class='small muted'>Thresholds: " + _esc(str(source[0])) + ".</p>"  # type: ignore[index]
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
    body += _cell_composition_html(meth.get("deconvolution") or {}, summ, thr, link)  # type: ignore[arg-type]
    return body


def _cell_composition_html(dec: dict[str, object], summ: dict[str, dict[str, dict[str, object]]], thr: dict[str, float], link) -> str:  # type: ignore[no-untyped-def]
    """Cell-type composition under Methylation: the cohort's composition per cell type, the members' spread for the
    main types, the members far from the cohort. A run from before the table (0.16 to 0.19 with run_deconvolve) has
    only each member's top cell type."""
    if not dec.get("available"):
        top = summ.get("deconvolution", {}).get("top_cell_type", {}) or {}
        if not top:
            return ""
        counts_c: dict[str, int] = {}
        for cell in top.values():
            counts_c[str(cell)] = counts_c.get(str(cell), 0) + 1
        return ("<h3>Cell-type composition (methbat deconvolve)</h3>"
                + sv.bars_h(sorted(((c, float(n), "") for c, n in counts_c.items()), key=lambda t: -t[1]), width=300, label_w=140)
                + "<p class='small muted'>Each member's largest estimated fraction; the stage's run predates the composition table (rerun cohort_methyl "
                  "for the full picture).</p>")
    stats: dict[str, dict[str, float]] = dec.get("stats", {})  # type: ignore[assignment]
    members: dict[str, dict[str, object]] = dec.get("members", {})  # type: ignore[assignment]
    group_of: dict[str, str] = dec.get("group_of", {})  # type: ignore[assignment]
    flags: dict[str, object] = dec.get("flags", {})  # type: ignore[assignment]
    dist: dict[str, object] = dec.get("distance", {})  # type: ignore[assignment]
    shown = [c for c in dec.get("cell_types", []) if stats[c]["max"] >= 0.01]  # type: ignore[union-attr]
    blood = sorted(float((m.get("groups") or {}).get("blood", 0.0)) for m in members.values())  # type: ignore[union-attr]
    pct = lambda v: f"{100 * float(v):.1f}%"  # noqa: E731
    body = ("<h3>Cell-type composition (methbat deconvolve)</h3>"
            + _cards([("members with an estimate", dec.get("n_members")), ("cell types above 1% in a member", len(shown)),
                      ("blood lineage, cohort median", pct(statistics.median(blood)) if blood else "–"),
                      ("members flagged cell_composition_outlier", len(flags))])
            + "<p class='small muted'><b>An estimate for research, not a blood count.</b> Each member's 5mC pileup is explained as a mixture of the atlas's "
              "39 cell types (nanomix, built from bisulfite data; non-negative least squares). The main lineages of a blood sample come out as expected "
              "(granulocytes, T cells, monocytes); fractions of a few percent in unrelated tissues are the method's noise, and closely related cell types "
              "share their signal (B and NK cells often read 0). Composition is the first thing to check behind a member with many methylation events.</p>")
    # the picture first: every member as one column, sorted by the cohort's main cell type
    main_types = shown[:5]
    if main_types and len(members) >= 2:
        columns = []
        for sid, m in members.items():
            fr: dict[str, float] = m.get("fractions") or {}  # type: ignore[assignment]
            vals = {c: float(fr.get(c, 0.0)) for c in main_types}
            vals["other"] = max(0.0, sum(float(v) for v in fr.values()) - sum(vals.values()))
            columns.append((sid, vals))
        columns.sort(key=lambda cv: (-cv[1][main_types[0]], cv[0]))
        body += (sv.ribbon(columns, main_types + ["other"], width=700, marked=set(flags))
                 + f"<p class='small muted'>One column per member (above {sv.RIBBON_COLUMNS} members, the mean of adjacent ones), sorted by "
                   f"{_esc(main_types[0])}: the five cell types with the largest cohort medians and the rest. The members flagged below stand at the "
                   "right with a red tick, one column each; hover for the member and its fractions.</p>")
    rows = [[c, group_of.get(c, "other"), pct(stats[c]["median"]), f"{pct(stats[c]['p5'])} – {pct(stats[c]['p95'])}", pct(stats[c]["max"]),
             sv.fmt(stats[c]["n_above_5pct"])] for c in shown]
    body += big_table(["cell type", "lineage", "median", "5th – 95th percentile", "maximum", "members at 5% or more"], rows,
                      what=f"{len(rows)} cell types above 1% in at least one member", open_below=20)
    main = [c for c in shown if stats[c]["median"] >= 0.02][:6]
    if main and len(members) >= 2:
        body += sv.strip_or_band([{"label": c, "points": sorted((sid, 100 * float((m.get("fractions") or {}).get(c, 0.0))) for sid, m in members.items()),  # type: ignore[union-attr]
                                   "outliers": sorted(flags), "threshold": None, "beyond": [], "digits": 0, "unit": "%"} for c in main], width=700)
        body += ("<p class='small muted'>The members' fractions for the cell types with a cohort median of 2% or more; red points are the members flagged "
                 "below.</p>")
    if flags:
        frows = []
        for sid in sorted(flags, key=lambda s: -float(members[s].get("distance") or 0.0)):  # type: ignore[arg-type]
            m = members[sid]
            fr = sorted((m.get("fractions") or {}).items(), key=lambda kv: -kv[1])[:4]  # type: ignore[union-attr]
            frows.append([link(sid), sv.fmt(m.get("distance"), 2), ", ".join(f"{c} {pct(v)}" for c, v in fr),
                          ", ".join(f"{g} {pct(v)}" for g, v in list((m.get("groups") or {}).items())[:3])])  # type: ignore[union-attr]
        body += ("<h4>Members far from the cohort's composition</h4>"
                 + big_table(["member", "distance", "largest fractions", "lineages"], frows, what=f"{len(frows)} members", open_below=30))
    lim = dist.get("limit")
    body += ("<p class='small muted'>Distance: half the summed absolute differences between a member's fractions and the cohort's medians (0 = the median "
             "composition, 1 = nothing shared). <code>cell_composition_outlier</code> (warn): beyond the cohort's median + "
             + _esc(sv.fmt(thr.get("cell_composition_mad_max", 4.0), 1)) + " scaled MADs and at least "
             + _esc(sv.fmt(thr.get("cell_composition_distance_min", 0.25), 2))
             + (f" (here {_esc(sv.fmt(lim, 2))})" if isinstance(lim, (int, float)) else " (needs four members)")
             + ": another tissue or sample type, a cell line, contamination, or an unusual blood picture. The stage's table has every fraction "
               "(<code>ugc-wgw export --domains methyl</code>).</p>")
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
                  ("sex mismatches", len(sex.get("mismatches", []))), ("cohort stages run", len(runs)),  # type: ignore[arg-type]
                  *([("aligner and caller sets", len((d.get("tools") or {}).get("groups", {})))] if (d.get("tools") or {}).get("mixed") else [])])  # type: ignore[union-attr]
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
                 + _kv([["rows of the trgt-lps table (locus x motif)", sv.fmt(lps.get("loci"))], ["samples", sv.fmt(len(lps.get("samples", [])))], ["size", _bytes(lps.get("bytes"))]]))  # type: ignore[arg-type]
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

    # tandem repeats (cohort_tr)
    trd: dict[str, object] = d.get("tandem_repeats", {}) or {}  # type: ignore[assignment]
    ctr: dict[str, object] = callsets.get("cohort_tr", {}) or {}
    if d.get("mode") != "assembly":
        add("sec-tr", "Tandem repeats (cohort_tr)", _tr_html(trd, ctr, dist.get("tr_outliers"), thr, link) if ctr and trd.get("available") else
            "<p class='muted'>cohort_tr has not run for this cohort (<code>ugc-wgw submit --cohort " + _esc(cid) + " --stage cohort_tr</code>).</p>")

    # methylation (cohort_methyl)
    meth: dict[str, object] = d.get("methylation", {}) or {}  # type: ignore[assignment]
    cme: dict[str, object] = callsets.get("cohort_methyl", {}) or {}
    if d.get("mode") != "assembly":
        add("sec-methyl", "Methylation (cohort_methyl)", _methyl_html(meth, cme, dist.get("methyl_events"), thr, link) if cme else
            "<p class='muted'>cohort_methyl has not run for this cohort (<code>ugc-wgw submit --cohort " + _esc(cid) + " --stage cohort_methyl</code>).</p>")

    # structural variants and repeats against methylation (cohort_mqtl, optional)
    cmq: dict[str, object] = callsets.get("cohort_mqtl", {}) or {}
    if d.get("mode") != "assembly":
        add("sec-mqtl", "Variants and methylation (cohort_mqtl)", _mqtl_html(cmq) if cmq.get("summary") else
            "<p class='muted'>cohort_mqtl has not run for this cohort (<code>ugc-wgw submit --cohort " + _esc(cid) + " --stage cohort_mqtl</code>).</p>")

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
            growth: dict[str, object] = cp.get("growth") or {}  # type: ignore[assignment]
            if growth.get("nonref_bp"):
                pbody += _growth_html(growth)
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

    nav = _nav(f"{cid} · {project_name}", sections, member_nav(present))
    return _page(f"{cid} cohort summary", nav, f"<h1>{_esc(cid)} cohort summary</h1>" + "".join(parts))
