"""`ugc-wgw export`: the cohort's numbers as files, one TSV per domain, from the members' per-sample digests and the
cohort stages' own outputs. The HTML pages are for reading; the export is for joining. The frequency tables (copy
number per Paraphase region, StarPhase diplotypes, HLA alleles) are computed here and drawn by the cohort page from
the same function, so the two never disagree. See docs/guide/07-results.md ("Cohort tables") and docs/DESIGN.md §8.1."""
from __future__ import annotations

import gzip
import shutil
import statistics
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path

from . import summary
from . import summary_parsers as sp
from .config import Config, read_version
from .db import DB, Cohort, RunRecord
from .inputs import SEX_TO_GROUP
from .layout import stage_dir
from .manifest import CodeInfo
from .outputs import Outputs
from .util import UgcError, diag, sha256_file, utc_now, write_json

EXPORT_SCHEMA = 1
DOMAINS = ("members", "coverage", "variants", "phasing", "cn", "pgx", "hla", "mito", "identity", "tr", "methyl", "assembly")
STRATA = ("none", "sex", "ancestry", "all")
# domain -> (cohort stage, ((output name, file suffix), ...)): the stage's own tables, copied byte for byte when it ran
COPY_DOMAINS: dict[str, tuple[str, tuple[tuple[str, str], ...]]] = {
    "identity": ("cohort_qc", (("qc_pairs", "pairs"), ("qc_samples", "samples"), ("qc_groups", "groups"), ("qc_summary", "summary"),
                               ("qc_members", "members"), ("qc_related_pairs", "related_pairs"), ("qc_ancestry", "ancestry"))),
    "tr": ("cohort_tr", (("tr_allele_spectrum", "spectrum"), ("tr_locus_summary", "loci"), ("tr_outlier_pairs", "outliers"),
                         ("tr_outliers_top", "outliers_top"), ("tr_pathogenic_screen", "screen"), ("tr_screen_loci", "screen_loci"),
                         ("tr_samples", "samples"), ("tr_summary", "summary"))),   # the first three are bgzipped: they keep .gz
    "methyl": ("cohort_methyl", (("cohort_background", "background"), ("methyl_summary", "summary"), ("methyl_regions", "regions"),
                                 ("methyl_imprinting", "imprinting"), ("methyl_samples", "samples"), ("methyl_events", "events"),
                                 ("comparisons", "compare"))),   # an array output: one file per pair, `<prefix>.compare.<baseline>_vs_<compare>.tsv`
}
ANCESTRY_SOURCE = ("cohort_qc", "qc_ancestry")   # somalier's ancestry TSV: `#sample_id`, `predicted_ancestry`, probabilities, PCs
ANCESTRY_NOTE = "ancestry is somalier's projection onto the 1000 Genomes superpopulations (AFR, AMR, EAS, EUR, SAS): a description of structure, not an origin"
STAGE_COLUMNS = ("singleton", "upstream", "downstream", "assembly", "assembly_call")
META_INTERNAL = ("tsv_line", "input_check")   # sample-sheet meta the driver writes, not sheet columns
KIVVI_REGIONS = (("kiv2", "KIV2", "LPA"), ("d4z4", "D4Z4", "DUX4"))
VARIANT_FIXED = ("SNV_count", "INDEL_count", "TSTV_ratio", "HETHOM_ratio", "depth_mean", "inferred_sex", "phased_basepairs",
                 "phase_block_ng50", "trgt_genotyped_count", "trgt_uncalled_count", "cpg_combined_count", "cpg_hap1_count",
                 "cpg_hap2_count", "methbat_methylated_count", "methbat_unmethylated_count", "methbat_asm_count",
                 "sv_DEL_count", "sv_INS_count", "sv_DUP_count", "sv_INV_count", "sv_BND_count", "sv_SWAP_count")


@dataclass
class ExportResult:
    out_dir: Path
    files: list[Path] = field(default_factory=list)
    manifest: dict[str, object] = field(default_factory=dict)
    skipped: list[str] = field(default_factory=list)   # "tr: cohort_tr has not run"
    notes: list[str] = field(default_factory=list)     # "ancestry strata: no cohort_qc run with an ancestry table"


# ---- small helpers --------------------------------------------------------------------------


def _dig(doc: object, *path: str) -> object:
    cur = doc
    for key in path:
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


def _cell(v: object) -> str:
    if v is None:
        return ""
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, float):
        return format(v, ".10g")
    if isinstance(v, (list, tuple)):
        return ",".join(_cell(x) for x in v)
    return str(v)


def parse_domains(items: list[str] | None) -> list[str]:
    """`--domains a,b --domains c` -> ordered list; unknown names are refused."""
    out: list[str] = []
    for item in items or []:
        for d in str(item).split(","):
            d = d.strip()
            if not d:
                continue
            if d not in DOMAINS:
                raise UgcError(f"unknown export domain {d!r} (known: {', '.join(DOMAINS)})")
            if d not in out:
                out.append(d)
    return out or list(DOMAINS)


def sex_group_of(digest: dict[str, object]) -> str:
    """The member's sex group for the strata, by the rule cohort_freq uses: the resolved sex when the run recorded one
    (0.14+), else the sheet, else the coverage inference; "" when none."""
    resolved = _dig(digest, "identity", "resolved")
    if isinstance(resolved, str) and SEX_TO_GROUP.get(resolved.upper()):
        return SEX_TO_GROUP[resolved.upper()]
    sheet = SEX_TO_GROUP.get(str(_dig(digest, "sample", "sex_sheet") or "").upper())
    inferred = SEX_TO_GROUP.get(str(_dig(digest, "coverage", "inferred_sex") or "").upper())
    return sheet or inferred or ""


def read_ancestry(path: Path) -> dict[str, str]:
    """somalier's ancestry TSV: sample id -> predicted superpopulation (rows without a prediction are skipped)."""
    cols, rows = sp.table(sp.read_lines(path))
    out: dict[str, str] = {}
    for r in rows:
        sid, label = r.get("sample_id", ""), r.get("predicted_ancestry", "")
        if sid and label:
            out[sid] = label
    return out


def strata_of(members: dict[str, dict[str, object]], kinds: tuple[str, ...] | list[str],
              ancestry: dict[str, str] | None = None) -> dict[str, dict[str, str]]:
    """{kind: {sample_id: label}}: `ALL` always (every member), `sex` (XX/XY where known), `ancestry` (where labelled)."""
    out: dict[str, dict[str, str]] = {"ALL": {sid: "ALL" for sid in members}}
    if "sex" in kinds:
        out["sex"] = {sid: g for sid, m in members.items() if (g := sex_group_of(m))}
    if "ancestry" in kinds and ancestry:
        out["ancestry"] = {sid: ancestry[sid] for sid in members if sid in ancestry}
    return out


def freq_table(rows: list[dict[str, object]], key: str, strata: dict[str, dict[str, str]], *, by: tuple[str, ...],
               called: str = "called") -> list[dict[str, object]]:
    """Counts and fractions of `key`'s values per group (`by` columns) and stratum, over the members with a call.
    `n_called` counts the group's called rows in the stratum, `n_members` its members with any row (called or not)."""
    out: list[dict[str, object]] = []
    for kind, labels in strata.items():
        counts: dict[tuple, Counter] = defaultdict(Counter)
        n_called: Counter = Counter()
        members: dict[tuple, set] = defaultdict(set)
        for r in rows:
            sid = str(r.get("sample_id"))
            label = labels.get(sid)
            if label is None:
                continue
            group = tuple(r.get(c) for c in by) + (label,)
            members[group].add(sid)
            if not r.get(called):
                continue
            n_called[group] += 1
            counts[group][r.get(key)] += 1
        for group in sorted(counts, key=lambda g: tuple(str(x) for x in g)):
            total = n_called[group]
            for value, n in counts[group].most_common():
                row: dict[str, object] = dict(zip(by, group[:-1]))
                row.update({key: value, "stratum_kind": kind, "stratum": group[-1], "count": n,
                            "fraction": round(n / total, 4) if total else None, "n_called": total, "n_members": len(members[group])})
                out.append(row)
    return out


def write_tsv(path: Path, rows: list[dict[str, object]], columns: list[str], comment: str) -> int:
    """One `#` comment line, the header, the rows; ints as is, floats to 10 significant digits, bools true/false."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fh.write(f"# {comment}\n")
        fh.write("\t".join(columns) + "\n")
        for r in rows:
            fh.write("\t".join(_cell(r.get(c)) for c in columns) + "\n")
    return len(rows)


# ---- per-member rows -------------------------------------------------------------------------


def member_rows(cohort: Cohort, digests: dict[str, dict[str, object]], db: DB, *, ancestry: dict[str, str],
                sex_somalier: dict[str, str]) -> tuple[list[dict[str, object]], list[str]]:
    rows: list[dict[str, object]] = []
    sheet_cols: list[str] = []
    for sid in cohort.members:
        rec = db.get_sample(sid)
        meta = {k: v for k, v in (rec.meta if rec else {}).items() if k not in META_INTERNAL and not isinstance(v, (dict, list))}
        for k in meta:
            if k not in sheet_cols:
                sheet_cols.append(k)
        m = digests.get(sid)
        row: dict[str, object] = {"sample_id": sid, "present": m is not None,
                                  "sex_sheet": (rec.sex if rec else None) or (_dig(m, "sample", "sex_sheet") if m else None),
                                  "sex_inferred": _dig(m, "coverage", "inferred_sex") if m else None,
                                  "sex_group": sex_group_of(m) if m else "",
                                  "sex_somalier": sex_somalier.get(sid, _dig(m, "identity", "somalier") if m else None) or "",
                                  "ancestry": ancestry.get(sid, ""),
                                  "father_id": (rec.father_id if rec else None) or "", "mother_id": (rec.mother_id if rec else None) or ""}
        for k, v in meta.items():
            row[f"sheet_{k}"] = v
        runs = _dig(m, "runs") if m else None
        for stage in STAGE_COLUMNS:
            r = (runs or {}).get(stage) if isinstance(runs, dict) else None
            row[f"v_{stage}"] = r.get("ugc_wgw_version") if isinstance(r, dict) else None
            row[f"run_{stage}"] = r.get("run_id") if isinstance(r, dict) else None
        qc = _dig(m, "qc") if m else None
        row["qc_status"] = qc.get("status") if isinstance(qc, dict) else None
        row["qc_flags"] = sorted({str(f.get("code")) for f in (qc.get("flags", []) if isinstance(qc, dict) else [])})
        row["n_warnings"] = len(m.get("warnings", [])) if m else None  # type: ignore[arg-type]
        row["digest_generated_at"] = m.get("generated_at") if m else None
        rows.append(row)
    return rows, sheet_cols


def coverage_rows(digests: dict[str, dict[str, object]]) -> list[dict[str, object]]:
    rows = []
    for sid, m in digests.items():
        cov: dict[str, object] = m.get("coverage", {}) if isinstance(m.get("coverage"), dict) else {}
        ratios: dict[str, object] = cov.get("ratios", {}) if isinstance(cov.get("ratios"), dict) else {}
        row: dict[str, object] = {"sample_id": sid, "depth_mean": cov.get("mean"),
                                  "haploid_coverage": _dig(cov, "copynum", "haploid_coverage"), "inferred_sex": cov.get("inferred_sex"),
                                  "sex_from_ratio": ratios.get("sex_from_ratio"), "chrX_ratio": ratios.get("chrX_ratio"),
                                  "chrY_ratio": ratios.get("chrY_ratio")}
        reads: dict[str, object] = m.get("reads", {}) if isinstance(m.get("reads"), dict) else {}
        for k in summary.READ_KEYS:
            row[k] = reads.get(k)
        for r in cov.get("chrom", []) if isinstance(cov.get("chrom"), list) else []:
            if isinstance(r, dict) and r.get("chrom"):
                row[f"depth_{r['chrom']}"] = r.get("mean")
        rows.append(row)
    return rows


def variants_rows(digests: dict[str, dict[str, object]]) -> list[dict[str, object]]:
    rows = []
    for sid, m in digests.items():
        smv: dict[str, object] = m.get("small_variants", {}) if isinstance(m.get("small_variants"), dict) else {}
        psc: dict[str, object] = smv.get("psc", {}) if isinstance(smv.get("psc"), dict) else {}
        roh: dict[str, object] = smv.get("roh", {}) if isinstance(smv.get("roh"), dict) else {}
        svb: dict[str, object] = m.get("sv", {}) if isinstance(m.get("sv"), dict) else {}
        counts: dict[str, object] = svb.get("counts", {}) if isinstance(svb.get("counts"), dict) else {}
        tr: dict[str, object] = m.get("tandem_repeats", {}) if isinstance(m.get("tandem_repeats"), dict) else {}
        me: dict[str, object] = m.get("methylation", {}) if isinstance(m.get("methylation"), dict) else {}
        row: dict[str, object] = {"sample_id": sid, "snv": smv.get("snv"), "indel": smv.get("indel"), "tstv": smv.get("tstv"),
                                  "hethom": smv.get("hethom"), "records": _dig(smv, "sn", "number of records"),
                                  "n_hets": psc.get("n_hets"), "n_nonref_hom": psc.get("n_nonref_hom"),
                                  "n_singletons": psc.get("n_singletons"), "n_missing": psc.get("n_missing"),
                                  "roh_n": roh.get("n"), "roh_total_bp": roh.get("total_bp"), "roh_longest_bp": roh.get("longest_bp")}
        sv_values = [counts.get(t) for t in ("DEL", "INS", "DUP", "INV", "BND", "SWAP")]
        for t, v in zip(("DEL", "INS", "DUP", "INV", "BND", "SWAP"), sv_values):
            row[f"sv_{t}"] = v
        row["sv_total"] = sum(_num(v) or 0 for v in sv_values) if any(v is not None for v in sv_values) else None
        row["sv_records"] = svb.get("n")
        row.update({"trgt_genotyped": tr.get("genotyped"), "trgt_uncalled": tr.get("uncalled"), "trgt_dropouts": _dig(tr, "dropouts", "n"),
                    "trgt_disease_called": tr.get("n_named_called"), "trgt_disease_total": tr.get("n_named"),
                    "cpg_combined": _dig(me, "cpg", "combined"), "cpg_hap1": _dig(me, "cpg", "hap1"), "cpg_hap2": _dig(me, "cpg", "hap2"),
                    "methbat_methylated": _dig(me, "methbat", "methylated"), "methbat_unmethylated": _dig(me, "methbat", "unmethylated"),
                    "methbat_asm": _dig(me, "methbat", "asm")})
        stats: dict[str, object] = m.get("stats", {}) if isinstance(m.get("stats"), dict) else {}
        for k in sorted(stats):
            if k in VARIANT_FIXED or k in summary.READ_KEYS or k in row:
                continue
            row[k] = stats[k]
        rows.append(row)
    return rows


def phasing_rows(digests: dict[str, dict[str, object]]) -> list[dict[str, object]]:
    rows = []
    for sid, m in digests.items():
        ph: dict[str, object] = m.get("phasing", {}) if isinstance(m.get("phasing"), dict) else {}
        allv: dict[str, object] = ph.get("all", {}) if isinstance(ph.get("all"), dict) else {}
        het, phased = _num(allv.get("num_heterozygous")), _num(allv.get("num_phased"))
        rows.append({"sample_id": sid, "phased_bp": ph.get("phased_bp"), "block_ng50": ph.get("block_ng50"),
                     "blocks": _dig(ph, "blocks", "n"), "block_n50": _dig(ph, "blocks", "n50"),
                     "num_variants": allv.get("num_variants"), "num_heterozygous": allv.get("num_heterozygous"),
                     "num_phased": allv.get("num_phased"), "num_unphased": allv.get("num_unphased"), "num_blocks": allv.get("num_blocks"),
                     "phased_fraction": round(phased / het, 4) if het and phased is not None else None})
    return rows


def cn_rows(digests: dict[str, dict[str, object]]) -> list[dict[str, object]]:
    """One row per member and Paraphase region, then one per kivvi repeat (KIV2, D4Z4)."""
    rows = []
    for sid, m in digests.items():
        para: dict[str, object] = _dig(m, "targeted", "paraphase") or {}  # type: ignore[assignment]
        for r in para.get("regions", []) if isinstance(para.get("regions"), list) else []:
            if not isinstance(r, dict):
                continue
            failed = bool(r.get("failed"))
            rows.append({"sample_id": sid, "source": "paraphase", "region": r.get("region"), "genes": r.get("genes"),
                         "phase_region": r.get("phase_region"), "total_cn": r.get("total_cn"), "gene_cn": r.get("gene_cn"),
                         "n_haplotypes": r.get("n_haplotypes"), "two_copy": r.get("two_copy"),
                         "haplotypes_by_gene": ";".join(f"{g}:{n}" for g, n in (r.get("haplotypes_by_gene") or {}).items()),  # type: ignore[union-attr]
                         "fusions": ";".join(r.get("fusions") or []), "failed": failed,  # type: ignore[arg-type]
                         "called": not failed and r.get("total_cn") is not None,
                         "depth_median": r.get("depth_median"), "depth_p80": r.get("depth_p80"), "genome_depth": r.get("genome_depth")})
        kivvi: dict[str, object] = _dig(m, "targeted", "kivvi") or {}  # type: ignore[assignment]
        for key, region, gene in KIVVI_REGIONS:
            k: dict[str, object] = kivvi.get(key) or {}  # type: ignore[assignment]
            if not k:
                continue
            cn = _num(k.get("allele_cn"))
            status = str(k.get("call_status") or "")
            failed = "fail" in status.lower()
            rows.append({"sample_id": sid, "source": "kivvi", "region": region, "genes": gene, "phase_region": "", "total_cn": cn,
                         "gene_cn": None, "n_haplotypes": k.get("n_complete_alleles"), "two_copy": None, "haplotypes_by_gene": "",
                         "fusions": "", "failed": failed, "called": cn is not None and not failed,
                         "depth_median": k.get("repeat_depth"), "depth_p80": None, "genome_depth": k.get("genome_depth")})
    return rows


def _starphase_genes(m: dict[str, object]) -> dict[str, dict[str, object]]:
    genes = _dig(m, "targeted", "starphase_meta", "genes")
    if isinstance(genes, dict) and genes:
        return genes  # type: ignore[return-value]
    out: dict[str, dict[str, object]] = {}
    for r in _dig(m, "targeted", "starphase") or []:  # type: ignore[union-attr]
        if isinstance(r, dict) and r.get("gene"):
            dip = str(r.get("diplotype") or "")
            hap1, hap2 = (dip.split("/", 1) + [""])[:2] if "/" in dip else (dip, dip)
            out[str(r["gene"])] = {"hap1": hap1, "hap2": hap2, "diplotype": dip, "simple_diplotype": dip, "n_diplotypes": 1 if dip else 0,
                                   "inexact": False, "called": bool(dip) and sp.starphase_called(hap1, hap2),
                                   "hla": str(r["gene"]).upper().startswith("HLA-")}
    return out


def pgx_rows(digests: dict[str, dict[str, object]]) -> list[dict[str, object]]:
    rows = []
    for sid, m in digests.items():
        for gene, g in _starphase_genes(m).items():
            if g.get("hla"):
                continue
            rows.append({"sample_id": sid, "gene": gene, "diplotype": g.get("diplotype"), "simple_diplotype": g.get("simple_diplotype"),
                         "haplotype_1": g.get("hap1"), "haplotype_2": g.get("hap2"), "n_diplotypes": g.get("n_diplotypes"),
                         "inexact": bool(g.get("inexact")), "called": bool(g.get("called"))})
    return rows


def hla_rows(digests: dict[str, dict[str, object]]) -> list[dict[str, object]]:
    rows = []
    for sid, m in digests.items():
        for gene, g in _starphase_genes(m).items():
            if g.get("hla"):
                rows.append({"sample_id": sid, "gene": gene, "allele_1": g.get("hap1"), "allele_2": g.get("hap2"),
                             "diplotype": g.get("diplotype"), "n_diplotypes": g.get("n_diplotypes"), "called": bool(g.get("called"))})
    return rows


def hla_allele_rows(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    """The HLA rows as one row per allele (both haplotypes), the input of the allele frequency table."""
    out = []
    for r in rows:
        for a in (r.get("allele_1"), r.get("allele_2")):
            out.append({"sample_id": r["sample_id"], "gene": r["gene"], "allele": a, "called": bool(r.get("called"))})
    return out


def mito_rows(digests: dict[str, dict[str, object]], cohort_haplogroups: list[dict[str, object]] | None = None) -> list[dict[str, object]]:
    """One row per member and mitorsaw haplotype, with the haplogroup from the member's own run or, for a member whose run
    predates it, from cohort_merge's gathered table (`haplogroup_source` = run or cohort_merge)."""
    gathered: dict[tuple[str, str], dict[str, object]] = {}
    for r in cohort_haplogroups or []:
        gathered[(str(r.get("sample_id")), str(r.get("haplotype")))] = r
    rows = []
    for sid, m in digests.items():
        mito: dict[str, object] = _dig(m, "targeted", "mitorsaw") or {}  # type: ignore[assignment]
        haps = mito.get("haplotypes") if isinstance(mito.get("haplotypes"), list) else []
        fp: dict[str, object] = mito.get("fingerprint") or {}  # type: ignore[assignment]
        collapsed: dict[str, object] = mito.get("collapsed") or {}  # type: ignore[assignment]
        flagged = any(f.get("code") == "mito_unexplained_high" for f in (_dig(m, "qc", "flags") or []) if isinstance(f, dict))  # type: ignore[union-attr]
        major = str(collapsed.get("major_haplotype") or "")
        if not major:
            # a run before 0.15 has no collapsed summary: the most abundant haplotype is the major one, as the collapse defines it
            # (the 864 cohort, 2026-11-09: 864 haplogroups made in cohort_merge and none counted, for want of this)
            best = max((h for h in haps if isinstance(h, dict)), key=lambda h: _num(h.get("estimated_abundance")) or 0.0, default=None)
            major = str(best.get("label") or "") if best else ""
        for h in haps:
            if not isinstance(h, dict):
                continue
            hg, hq, source = h.get("haplogroup") or "", h.get("haplogroup_quality"), "run"
            if not hg:
                g = gathered.get((sid, str(h.get("label"))))
                if g:
                    hg, hq, source = str(g.get("haplogroup") or ""), g.get("quality"), "cohort_merge"
            rows.append({"sample_id": sid, "haplotype": h.get("label"), "n_haplotypes": len(haps), "seq_len": h.get("seq_len"),
                         "num_ref_variants": h.get("num_ref_variants"), "estimated_abundance": h.get("estimated_abundance"),
                         "unexplained_fraction": fp.get("unexplained_fraction"), "passing_explanation": fp.get("passing_explanation"),
                         "mito_unexplained_high": flagged, "haplogroup": hg, "haplogroup_quality": hq if hg else "",
                         "haplogroup_source": source if hg else "", "major_haplotype": bool(major) and str(h.get("label")) == major,
                         "homoplasmic_sites": collapsed.get("homoplasmic"), "heteroplasmic_sites": collapsed.get("heteroplasmic"),
                         "het_vaf_max": collapsed.get("het_vaf_max")})
    return rows


def mito_freq_rows(hg_major: list[dict[str, object]], strata: dict[str, dict[str, str]]) -> list[dict[str, object]]:
    """The haplogroup frequencies: one row per (kind, value, stratum), kind = haplogroup (the major haplotype's) or
    macro_haplogroup (its first letter)."""
    out: list[dict[str, object]] = []
    for kind in ("haplogroup", "macro_haplogroup"):
        for r in freq_table(hg_major, kind, strata, by=()):
            out.append({"kind": kind, "value": r[kind], **{k: v for k, v in r.items() if k != kind}})
    return out


def haplogroup_rows(mito: list[dict[str, object]]) -> list[dict[str, object]]:
    """One row per member: the major haplotype's haplogroup and its macro-haplogroup (the first letter), for freq_table."""
    out = []
    for r in mito:
        if r.get("major_haplotype") and r.get("haplogroup"):
            hg = str(r["haplogroup"])
            out.append({"sample_id": r["sample_id"], "haplogroup": hg, "macro_haplogroup": hg[:1], "called": True})
    return out


def assembly_rows(digests: dict[str, dict[str, object]]) -> list[dict[str, object]]:
    rows = []
    for sid, m in digests.items():
        asm: dict[str, object] = m.get("assembly") or {}  # type: ignore[assignment]
        if not asm:
            continue
        row: dict[str, object] = {"sample_id": sid, "trio": asm.get("trio"), "assembler": asm.get("assembler"),
                                  "assembler_version": asm.get("assembler_version"), "karyotype": _dig(asm, "karyotype", "karyotype")}
        for k, v in (asm.get("summary") or {}).items():  # type: ignore[union-attr]
            row[k] = v
        rows.append(row)
    return rows


# ---- the targeted callers across the cohort (shared with the cohort page) ------------------------


def targeted_summary(members: dict[str, dict[str, object]], strata: dict[str, dict[str, str]], *, top_regions: int = 12,
                     cohort_haplogroups: list[dict[str, object]] | None = None) -> dict[str, object]:
    """Copy-number counts per Paraphase region, StarPhase diplotype and HLA allele frequencies, the mitorsaw fingerprint
    per member: what the cohort page draws and the export writes, from one computation."""
    cn = cn_rows(members)
    para = [r for r in cn if r.get("source") == "paraphase"]
    regions: dict[str, dict[str, object]] = {}
    for r in para:
        reg = regions.setdefault(str(r["region"]), {"region": r["region"], "genes": r.get("genes"), "n_members": 0, "n_called": 0,
                                                     "n_failed": 0, "cn_counts": Counter(), "gene_cn_counts": Counter(), "values": []})
        reg["n_members"] += 1  # type: ignore[operator]
        if r.get("failed"):
            reg["n_failed"] += 1  # type: ignore[operator]
        if r.get("called"):
            reg["n_called"] += 1  # type: ignore[operator]
            reg["cn_counts"][r.get("total_cn")] += 1  # type: ignore[index]
            reg["values"].append(float(r["total_cn"]))  # type: ignore[union-attr, arg-type]
            if r.get("gene_cn") is not None:
                reg["gene_cn_counts"][r.get("gene_cn")] += 1  # type: ignore[index]
    region_rows = []
    for reg in regions.values():
        vals: list[float] = reg.pop("values")  # type: ignore[assignment]
        reg["variance"] = round(statistics.pvariance(vals), 4) if len(vals) > 1 else 0.0
        reg["cn_counts"] = {str(k): v for k, v in sorted(reg["cn_counts"].items(), key=lambda kv: (kv[0] is None, kv[0]))}  # type: ignore[union-attr]
        reg["gene_cn_counts"] = {str(k): v for k, v in sorted(reg["gene_cn_counts"].items(), key=lambda kv: (kv[0] is None, kv[0]))}  # type: ignore[union-attr]
        region_rows.append(reg)
    region_rows.sort(key=lambda r: (-float(r["variance"]), -int(r["n_called"]), str(r["region"])))  # type: ignore[arg-type]
    cn_freq = freq_table([dict(r, metric="total_cn", value=r.get("total_cn")) for r in para], "value", strata, by=("region", "genes", "metric"))
    cn_freq += freq_table([dict(r, metric="gene_cn", value=r.get("gene_cn"), called=bool(r.get("called")) and r.get("gene_cn") is not None)
                           for r in para], "value", strata, by=("region", "genes", "metric"))
    pgx = pgx_rows(members)
    pgx_freq = freq_table(pgx, "simple_diplotype", strata, by=("gene",))
    genes: dict[str, dict[str, object]] = {}
    for r in pgx:
        g = genes.setdefault(str(r["gene"]), {"gene": r["gene"], "n_members": 0, "n_called": 0})
        g["n_members"] += 1  # type: ignore[operator]
        g["n_called"] += int(bool(r.get("called")))  # type: ignore[operator]
    for g in genes.values():
        g["called_rate"] = round(int(g["n_called"]) / int(g["n_members"]), 4) if g["n_members"] else None  # type: ignore[call-overload]
        g["top"] = [(str(r["simple_diplotype"]), int(r["count"])) for r in pgx_freq  # type: ignore[call-overload]
                    if r["gene"] == g["gene"] and r["stratum_kind"] == "ALL"][:6]
    hla = hla_rows(members)
    hla_freq = freq_table(hla_allele_rows(hla), "allele", strata, by=("gene",))
    hla_genes: dict[str, dict[str, object]] = {}
    for r in hla:
        g = hla_genes.setdefault(str(r["gene"]), {"gene": r["gene"], "n_members": 0, "n_called": 0})
        g["n_members"] += 1  # type: ignore[operator]
        g["n_called"] += int(bool(r.get("called")))  # type: ignore[operator]
    mito = mito_rows(members, cohort_haplogroups)
    hg_major = haplogroup_rows(mito)
    hg_counts: dict[str, int] = {}
    for r in hg_major:
        hg_counts[str(r["haplogroup"])] = hg_counts.get(str(r["haplogroup"]), 0) + 1
    unexplained = {str(r["sample_id"]): _num(r.get("unexplained_fraction")) for r in mito if r.get("haplotype") is not None}
    multi = sorted({str(r["sample_id"]) for r in mito if (r.get("n_haplotypes") or 0) > 1})  # type: ignore[operator]
    minor = [float(r["estimated_abundance"]) for r in mito if r["sample_id"] in multi and _num(r.get("estimated_abundance")) is not None  # type: ignore[arg-type]
             and float(r["estimated_abundance"]) < 0.5]  # type: ignore[arg-type]
    return {"strata": sorted(strata), "n_members": len(members),
            "paraphase": {"regions": region_rows[:top_regions], "n_regions": len(region_rows),
                          "n_regions_called": sum(1 for r in region_rows if r["n_called"]),
                          "n_regions_gene_cn": sum(1 for r in region_rows if r["gene_cn_counts"]),
                          "n_members_with_failed": len({str(r["sample_id"]) for r in para if r.get("failed")}), "freq": cn_freq},
            "pgx": {"genes": genes, "freq": pgx_freq},
            "hla": {"genes": hla_genes, "freq": hla_freq},
            "mito": {"n_members": len({str(r["sample_id"]) for r in mito}), "multi_haplotype": multi,
                     "minor_abundance_median": round(statistics.median(minor), 4) if minor else None,
                     "unexplained": {k: v for k, v in unexplained.items() if v is not None},
                     "flagged": sorted({str(r["sample_id"]) for r in mito if r.get("mito_unexplained_high")}),
                     "n_haplogrouped": len(hg_major), "haplogroups": hg_counts, "freq": mito_freq_rows(hg_major, strata)}}


# ---- copies of the cohort stages' tables --------------------------------------------------------


def copy_outputs(run: RunRecord, specs: tuple[tuple[str, str], ...], out_dir: Path, prefix: str) -> tuple[list[tuple[Path, str]], list[str]]:
    """Copy the named outputs of a cohort run as `<prefix>.<suffix>.tsv`; optional outputs that are absent are listed."""
    outs = Outputs.load(run)
    copied: list[tuple[Path, str]] = []
    absent: list[str] = []
    for name, suffix in specs:
        p = outs.opt(name)
        if isinstance(p, list):
            # an array output (the comparisons): each file keeps what distinguishes it, the part after `.<suffix>.` of its name
            files = [f for f in p if isinstance(f, str) and Path(f).exists()]
            if not files:
                absent.append(name)
                continue
            for f in files:
                tail = Path(f).name.split(f".{suffix}.", 1)[1] if f".{suffix}." in Path(f).name else Path(f).name
                tail = tail[:-4] if tail.endswith(".tsv") else tail
                dest = out_dir / f"{prefix}.{suffix}.{tail}.tsv"
                shutil.copyfile(f, dest)
                copied.append((dest, name))
            continue
        if not isinstance(p, str) or not Path(p).exists():
            absent.append(name)
            continue
        dest = out_dir / (f"{prefix}.{suffix}.tsv" + (".gz" if p.endswith(".gz") else ""))
        shutil.copyfile(p, dest)
        copied.append((dest, name))
    return copied, absent


def _count_rows(path: Path) -> int:
    """The data rows of a written table (comment lines and the header excluded), plain or gzipped."""
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rt") as fh:  # type: ignore[operator]
        return max(0, sum(1 for line in fh if not line.startswith("#")) - 1)


# ---- the export ---------------------------------------------------------------------------------


def build(cfg: Config, db: DB, code: CodeInfo, *, cohort_id: str, domains: list[str], strata: str, out_dir: Path | None,
          any_version: bool, mode: str = "standalone", force: bool = False, thresholds: dict[str, float] | None = None,
          jobs: int = 1, summary_dir: Path | None = None) -> ExportResult:
    """Write the cohort's tables. Digests are shared with `ugc-wgw summary` (reused when fresh, rebuilt and written back
    otherwise); the cohort stages' tables are copied when the stage has a successful run in `mode`."""
    cohort = db.get_cohort(cohort_id)
    if cohort is None:
        raise UgcError(f"unknown cohort {cohort_id}")
    if strata not in STRATA:
        raise UgcError(f"unknown strata {strata!r} (known: {', '.join(STRATA)})")
    version = read_version(cfg.code_dir)
    lookup = None if any_version else version
    thr = thresholds or dict(summary.DEFAULT_THRESHOLDS)
    dest = out_dir or stage_dir(cfg.results_dir, "cohort", cohort_id, version, "export")
    res = ExportResult(out_dir=dest)
    if dest.exists() and any(dest.iterdir()) and not (dest / "export_manifest.json").exists() and not force:
        raise UgcError(f"{dest} holds files but no export_manifest.json; --force to write into it, or --out DIR")
    collected = summary.collect_digests(cfg, db, code, sample_ids=list(cohort.members), mode=mode, version=lookup,
                                        out_dir=summary_dir or cfg.results_dir / "reports" / "summary", force=force, thresholds=thr, jobs=jobs)
    digests = {sid: collected.digests[sid] for sid in cohort.members if sid in collected.digests}
    for sid in collected.skipped:
        res.notes.append(f"{sid}: no successful {mode} run; absent from the per-member tables")
    if not digests:
        res.notes.append(f"no member of {cohort_id} has a successful {mode} run")
        return res
    cruns: dict[str, RunRecord] = {}
    for stage in sorted({v[0] for v in COPY_DOMAINS.values()} | {ANCESTRY_SOURCE[0]}):
        run = summary.latest_success_in_mode(db, "cohort", cohort_id, stage, mode, lookup)
        if run is not None:
            cruns[stage] = run
    ancestry: dict[str, str] = {}
    ancestry_run = cruns.get(ANCESTRY_SOURCE[0])
    ancestry_path: str | None = None
    if ancestry_run is not None:
        p = Outputs.load(ancestry_run).opt(ANCESTRY_SOURCE[1])
        if isinstance(p, str) and Path(p).exists():
            ancestry_path = p
            ancestry = read_ancestry(Path(p))
    kinds: tuple[str, ...] = {"none": (), "sex": ("sex",), "ancestry": ("ancestry",), "all": ("sex", "ancestry")}[strata]
    if "ancestry" in kinds and not ancestry:
        res.notes.append("ancestry strata: no cohort_qc run with an ancestry table for this cohort; ALL" + (" and sex" if "sex" in kinds else "") + " only")
    st = strata_of(digests, kinds, ancestry=ancestry)
    unlabelled = sorted(sid for sid in digests if "ancestry" in st and sid not in st["ancestry"])
    cohort_haplogroups: list[dict[str, object]] = []
    if "cohort_merge" in cruns:
        p = Outputs.load(cruns["cohort_merge"]).opt("cohort_mt_haplogroups")
        if isinstance(p, str) and Path(p).exists():
            cohort_haplogroups = sp.parse_haplogroups(sp.read_lines(Path(p)))
    sex_somalier: dict[str, str] = {}
    if "cohort_qc" in cruns:
        p = Outputs.load(cruns["cohort_qc"]).opt("qc_members")
        if isinstance(p, str) and Path(p).exists():
            _cols, rows = sp.table(sp.read_lines(Path(p)))
            sex_somalier = {r.get("sample_id", ""): r.get("somalier_sex", "") for r in rows if r.get("sample_id")}

    dest.mkdir(parents=True, exist_ok=True)
    stamp = utc_now()
    prefix = cohort_id
    written: list[tuple[str, Path, int]] = []   # domain, path, rows
    domain_files: dict[str, object] = {}

    def comment(domain: str, what: str) -> str:
        return f"ugc-wgw {version} export · cohort {cohort_id} ({mode}) · {domain} · generated {stamp}; {what}"

    def emit(domain: str, suffix: str, rows: list[dict[str, object]], columns: list[str], what: str) -> None:
        path = dest / f"{prefix}.{suffix}.tsv"
        n = write_tsv(path, rows, columns, comment(domain, what))
        written.append((domain, path, n))
        domain_files.setdefault(domain, {"files": []})["files"].append(path.name)  # type: ignore[index]

    strata_note = "strata: ALL" + (", sex (XX/XY by the cohort_freq rule: resolved sex, else the sheet, else coverage)" if "sex" in st else "") \
        + (f"; {ANCESTRY_NOTE}" if "ancestry" in st else "")
    for domain in domains:
        if domain == "members":
            rows, sheet_cols = member_rows(cohort, digests, db, ancestry=ancestry, sex_somalier=sex_somalier)
            cols = ["sample_id", "present", "sex_sheet", "sex_inferred", "sex_group", "sex_somalier", "ancestry", "father_id", "mother_id"] \
                + [f"sheet_{c}" for c in sheet_cols] + [f"{p}_{s}" for s in STAGE_COLUMNS for p in ("v", "run")] \
                + ["qc_status", "qc_flags", "n_warnings", "digest_generated_at"]
            emit(domain, "members", rows, cols, "one row per cohort member in frozen order; present = has a digest; " + strata_note)
        elif domain == "coverage":
            rows = coverage_rows(digests)
            cols = _columns(rows, ["sample_id"])
            emit(domain, "coverage", rows, cols, "mosdepth depth, sawfish haploid coverage, chrX/chrY ratios, read metrics, mean depth per chromosome")
        elif domain == "variants":
            rows = variants_rows(digests)
            emit(domain, "variants", rows, _columns(rows, ["sample_id"]), "per-sample variant counts; the trailing columns are every remaining stat_* value of the stats file")
        elif domain == "phasing":
            rows = phasing_rows(digests)
            emit(domain, "phasing", rows, _columns(rows, ["sample_id"]), "HiPhase statistics; phased_fraction = num_phased / num_heterozygous")
        elif domain == "cn":
            rows = cn_rows(digests)
            cols = ["sample_id", "source", "region", "genes", "phase_region", "total_cn", "gene_cn", "n_haplotypes", "two_copy",
                    "haplotypes_by_gene", "fusions", "failed", "called", "depth_median", "depth_p80", "genome_depth"]
            emit(domain, "cn", rows, cols, "one row per member and Paraphase region, then kivvi KIV2 and D4Z4 (continuous copy number)")
            tg = targeted_summary(digests, st)
            freq = tg["paraphase"]["freq"]  # type: ignore[index]
            emit(domain, "cn.freq", freq, ["region", "genes", "metric", "value", "stratum_kind", "stratum", "count", "fraction", "n_called", "n_members"],  # type: ignore[arg-type]
                 "Paraphase copy-number counts per region (total_cn, and gene_cn where Paraphase separates the gene from its paralogs); kivvi excluded; " + strata_note)
        elif domain == "pgx":
            rows = pgx_rows(digests)
            emit(domain, "pgx", rows, ["sample_id", "gene", "diplotype", "simple_diplotype", "haplotype_1", "haplotype_2", "n_diplotypes", "inexact", "called"],
                 "StarPhase diplotypes per member and non-HLA gene; called = both haplotypes typed")
            emit(domain, "pgx.freq", freq_table(rows, "simple_diplotype", st, by=("gene",)),
                 ["gene", "simple_diplotype", "stratum_kind", "stratum", "count", "fraction", "n_called", "n_members"],
                 "StarPhase core-allele diplotype counts per gene over the called members; " + strata_note)
        elif domain == "hla":
            rows = hla_rows(digests)
            emit(domain, "hla", rows, ["sample_id", "gene", "allele_1", "allele_2", "diplotype", "n_diplotypes", "called"],
                 "StarPhase HLA alleles per member and gene")
            hf = freq_table(hla_allele_rows(rows), "allele", st, by=("gene",))
            for r in hf:
                r["n_alleles"] = r["n_called"]
                r["n_called"] = int(r["n_called"]) // 2  # type: ignore[call-overload]
            emit(domain, "hla.freq", hf, ["gene", "allele", "stratum_kind", "stratum", "count", "fraction", "n_alleles", "n_called", "n_members"],
                 "HLA allele counts per gene, both haplotypes pooled (n_alleles = 2 × n_called); " + strata_note)
        elif domain == "mito":
            rows = mito_rows(digests, cohort_haplogroups)
            emit(domain, "mito", rows, ["sample_id", "haplotype", "n_haplotypes", "seq_len", "num_ref_variants", "estimated_abundance", "major_haplotype",
                                        "haplogroup", "haplogroup_quality", "haplogroup_source", "homoplasmic_sites", "heteroplasmic_sites", "het_vaf_max",
                                        "unexplained_fraction", "passing_explanation", "mito_unexplained_high"],
                 "mitorsaw haplotypes per member with haplogrep3's haplogroup per haplotype (from the member's run, or from cohort_merge for a run "
                 "before 0.15); the collapsed chrM counts are the sample's (major_haplotype marks the most abundant)")
            hg_major = haplogroup_rows(rows)
            emit(domain, "mito.freq", mito_freq_rows(hg_major, st),
                 ["kind", "value", "stratum_kind", "stratum", "count", "fraction", "n_called", "n_members"],
                 "haplogroup of each member's major haplotype, and its macro-haplogroup (first letter), counted per stratum; " + strata_note)
        elif domain == "assembly":
            rows = assembly_rows(digests)
            if not rows:
                res.skipped.append("assembly: no member has an assembly")
                domain_files[domain] = {"skipped": "no member has an assembly"}
                continue
            emit(domain, "assembly", rows, _columns(rows, ["sample_id", "trio", "assembler", "assembler_version", "karyotype"]),
                 "one row per member with an assembly: the flat assembly summary the cohort page draws")
        elif domain in COPY_DOMAINS:
            stage, specs = COPY_DOMAINS[domain]
            run = cruns.get(stage)
            if run is None:
                res.skipped.append(f"{domain}: {stage} has not run for {cohort_id} in {mode} mode")
                domain_files[domain] = {"skipped": f"{stage} has not run"}
                continue
            copied, absent = copy_outputs(run, specs, dest, f"{prefix}.{domain}")
            for path, name in copied:
                written.append((domain, path, -1))
                domain_files.setdefault(domain, {"files": [], "source_run": run.run_id})["files"].append(path.name)  # type: ignore[index]
            if absent:
                domain_files.setdefault(domain, {"files": [], "source_run": run.run_id})["absent"] = absent  # type: ignore[index]
    if force:
        keep = {p.name for _d, p, _n in written}
        for stale in sorted(list(dest.glob(f"{prefix}.*.tsv")) + list(dest.glob(f"{prefix}.*.tsv.gz"))):
            if stale.name not in keep:
                stale.unlink()
                res.notes.append(f"removed stale {stale.name}")
    files_meta = []
    for _domain, path, n in written:
        files_meta.append({"name": path.name, "rows": n if n >= 0 else _count_rows(path),
                           "bytes": path.stat().st_size, "sha256": sha256_file(path)})
    res.files = [p for _d, p, _n in written]
    res.manifest = {
        "schema": EXPORT_SCHEMA, "generated_at": stamp, "ugc_wgw": {"version": version, "git_commit": code.git_commit},
        "cohort": {"id": cohort_id, "frozen_at": cohort.frozen_at, "sha256": cohort.sample_list_sha256, "n": len(cohort.members)},
        "mode": mode, "any_version": any_version, "digest_schema": summary.DIGEST_SCHEMA, "thresholds": dict(thr),
        "members": {"n": len(cohort.members), "with_digest": len(digests), "missing": [s for s in cohort.members if s not in digests]},
        "strata": {"requested": strata, "kinds": sorted(st),
                   "ancestry": {"source_run": ancestry_run.run_id if ancestry_run else None, "path": ancestry_path, "note": ANCESTRY_NOTE,
                                "unlabelled": unlabelled} if "ancestry" in st else None},
        "source_runs": {"members": {sid: {s: r.get("run_id") for s, r in (m.get("runs") or {}).items() if isinstance(r, dict)}  # type: ignore[union-attr]
                                    for sid, m in digests.items()},
                        "cohort": {s: r.run_id for s, r in cruns.items()}},
        "domains": domain_files, "files": files_meta}
    write_json(dest / "export_manifest.json", res.manifest)
    return res


def _columns(rows: list[dict[str, object]], first: list[str]) -> list[str]:
    cols = list(first)
    for r in rows:
        for k in r:
            if k not in cols:
                cols.append(k)
    return cols
