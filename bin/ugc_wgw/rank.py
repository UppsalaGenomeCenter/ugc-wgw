"""`ugc-wgw cohort rank`: order a cohort's members for a pangenome graph by what each adds (docs/DESIGN.md §18,
2026-11-13).

The read-based cohort already says which member carries which structural variant (cohort_merge's SV VCF, one
genotype column per member). The ranking is a greedy cover over it, the idea of SVCollector and of the HPRC's
MaxVar: at every step the member that adds the most SV alleles no earlier member carries, counting the alleles
that are common in the cohort first. A member related to an earlier pick is held back to a second tier, so a
prefix of the ranking is a set of unrelated members. The same pass gives the saturation curves: how much of the
cohort's variation the first n members hold, against the expectation for n members drawn at random, and how much
of a member outside the first n is already there.

Every ladder cohort is a prefix of one saved ranking (`--emit K`), so the cohorts are nested by construction.
Stdlib only; reading the 864-member VCF takes minutes on one core.
"""
from __future__ import annotations

import io
import math
from array import array
from bisect import bisect_left
from dataclasses import dataclass, field
from pathlib import Path

from . import summary_parsers as sp
from . import summary_svg as sv
from .cohorts import read_ids
from .config import Config
from .db import DB
from .manifest import CodeInfo
from .outputs import Outputs
from .report import _esc, _table, big_table
from .util import UgcError, diag

SV_TYPES = ("INS", "DEL", "DUP", "INV")          # BND records come in pairs and CNV calls overlap the DEL/DUP ones
CLASSES = ("ge5", "common", "rare", "singleton")   # by carrier frequency among the ranked members
CLASS_LABELS = {"ge5": "5% and above", "common": "common, below 5%", "rare": "rare, two carriers or more", "singleton": "one carrier"}
CHECKPOINTS = (20, 50, 100, 200, 400)
NON_CARRIER = frozenset({"./.", ".|.", ".", "0/0", "0|0", "0", "0/.", "./0"})
RANK_COLUMNS = ["rank", "sample_id", "tier", "related_to", "relatedness", "new_common", "new_all", "cum_ge5", "cum_common", "cum_rare",
                "cum_singleton", "cum_all", "sv_alleles", "sex", "depth_mean", "read_n50", "qc_status", "ancestry",
                "PC1", "PC2", "PC3", "PC4", "PC5", "reason"]


@dataclass
class Matrix:
    """The carriers of every SV site among the wanted members: `site_carriers[s]` holds member indexes."""
    members: list[str]
    site_carriers: list[array] = field(default_factory=list)
    by_type: dict[str, int] = field(default_factory=dict)
    records: int = 0                                  # records of the wanted types, with or without a carrier
    missing: list[str] = field(default_factory=list)  # wanted members without a column


def _svtype(info: str) -> str:
    i = info.find("SVTYPE=")
    if i < 0 or (i > 0 and info[i - 1] != ";"):
        return ""
    j = info.find(";", i)
    return info[i + 7:j if j >= 0 else len(info)]


def read_sv_matrix(path: Path, members: list[str], types: tuple[str, ...] = SV_TYPES) -> Matrix:
    """Which of `members` carries each SV of `types` in a merged SV VCF. A genotype with any non-reference allele
    is a carrier; `./.` is svx's way of writing a member without the call, so it is not one."""
    out = Matrix(members=list(members))
    cols: list[tuple[int, int]] = []
    with io.TextIOWrapper(sp.open_binary(path), encoding="ascii", errors="replace") as fh:   # gzip or plain, by content
        for line in fh:
            if line.startswith("##"):
                continue
            f = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                where = {name: i for i, name in enumerate(f)}
                cols = [(where[m], k) for k, m in enumerate(members) if m in where and where[m] >= 9]
                out.missing = [m for m in members if m not in where]
                continue
            if len(f) < 10 or not cols:
                continue
            kind = _svtype(f[7])
            if kind not in types:
                continue
            out.records += 1
            if not f[8].startswith("GT"):
                raise UgcError(f"{path}: FORMAT of {f[0]}:{f[1]} does not start with GT")
            carriers = array("H")
            for col, k in cols:
                gt = f[col].split(":", 1)[0]
                if gt in NON_CARRIER:
                    continue
                if any(ch not in "0./|" for ch in gt):
                    carriers.append(k)
            if carriers:
                out.site_carriers.append(carriers)
                out.by_type[kind] = out.by_type.get(kind, 0) + 1
            if out.records % 50000 == 0:
                diag(f"  {out.records} SV records read")
    if not cols:
        raise UgcError(f"{path}: no #CHROM line with a column for any of the {len(members)} members")
    return out


def site_classes(matrix: Matrix, common_min: float) -> bytearray:
    """Per site, the index into CLASSES of its carrier-frequency class among the matrix's members."""
    n = max(1, len(matrix.members))
    out = bytearray(len(matrix.site_carriers))
    for s, carriers in enumerate(matrix.site_carriers):
        c = len(carriers)
        out[s] = 3 if c == 1 else (0 if c / n >= 0.05 else (1 if c / n >= common_min else 2))
    return out


@dataclass
class Ranking:
    order: list[int]                         # member indexes, best first
    tier: dict[int, str]                     # first | unrelated | related
    related_to: dict[int, tuple[int, float]]  # a held-back member -> (the earlier pick, relatedness)
    new_common: list[int]                    # per position: alleles added at carrier frequency >= common_min
    new_all: list[int]
    cum: list[list[int]]                     # per position: sites covered so far, per class
    totals: list[int]                        # sites per class
    first_pos: array                         # per site: the position of its first carrier in the order
    member_sites: list[array]


def greedy_rank(matrix: Matrix, klass: bytearray, *, eligible: list[int], first: list[int],
                related: dict[int, dict[int, float]]) -> Ranking:
    """Order `eligible` (member indexes). `first` are placed at the top as given; then the unrelated tier by greedy
    cover (most not-yet-covered common alleles, ties by all alleles, then by the cohort's own order); a member
    related to a pick is held back and the held-back members are ranked the same way afterwards."""
    n, n_sites = len(matrix.members), len(matrix.site_carriers)
    member_sites = [array("i") for _ in range(n)]
    for s, carriers in enumerate(matrix.site_carriers):
        for m in carriers:
            member_sites[m].append(s)
    is_common = bytearray(1 if k <= 1 else 0 for k in klass)
    gain_all = [len(ms) for ms in member_sites]
    gain_common = [sum(is_common[s] for s in ms) for ms in member_sites]
    covered = bytearray(n_sites)
    first_pos = array("i", [-1]) * n_sites
    totals = [0, 0, 0, 0]
    for k in klass:
        totals[k] += 1
    running = [0, 0, 0, 0]
    rk = Ranking(order=[], tier={}, related_to={}, new_common=[], new_all=[], cum=[], totals=totals, first_pos=first_pos,
                 member_sites=member_sites)
    pool = set(eligible)
    held: set[int] = set()

    def pick(m: int, tier: str) -> None:
        pos = len(rk.order)
        rk.order.append(m)
        rk.tier[m] = tier
        rk.new_common.append(gain_common[m])
        rk.new_all.append(gain_all[m])
        for s in member_sites[m]:
            if covered[s]:
                continue
            covered[s] = 1
            first_pos[s] = pos
            running[klass[s]] += 1
            common = is_common[s]
            for k in matrix.site_carriers[s]:
                gain_all[k] -= 1
                if common:
                    gain_common[k] -= 1
        rk.cum.append(list(running))
        pool.discard(m)
        held.discard(m)
        if tier != "related":
            for r, rel in related.get(m, {}).items():
                if r in pool:
                    pool.discard(r)
                    held.add(r)
                    rk.related_to[r] = (m, rel)

    def best(candidates: set[int]) -> int:
        return max(candidates, key=lambda k: (gain_common[k], gain_all[k], -k))

    for m in first:
        if m in pool or m in held:
            rk.related_to.pop(m, None)
            pick(m, "first")
    while pool:
        pick(best(pool), "unrelated")
    while held:
        pick(best(held), "related")
    return rk


def expected_random(totals_by_count: list[dict[int, int]], n: int) -> list[list[float]]:
    """Per class, the expected number of sites covered by k members drawn at random, k = 1..n: a site with c
    carriers among n members is missed by a draw of k with probability C(n-c, k) / C(n, k)."""
    lg = [math.lgamma(i + 1) for i in range(n + 1)]

    def log_choose(a: int, b: int) -> float:
        return lg[a] - lg[b] - lg[a - b]

    out: list[list[float]] = []
    for k in range(1, n + 1):
        row = []
        for counts in totals_by_count:
            exp = 0.0
            for c, sites in counts.items():
                miss = 0.0 if k > n - c else math.exp(log_choose(n - c, k) - log_choose(n, k))
                exp += sites * (1.0 - miss)
            row.append(exp)
        out.append(row)
    return out


def held_out(rk: Ranking, klass: bytearray, checkpoints: tuple[int, ...] = CHECKPOINTS) -> list[dict[str, object]]:
    """For each checkpoint K below the ranking's length: of every member after the first K, the share of its own
    SV alleles (all, and the common ones) that the first K members already carry."""
    per_member: dict[int, tuple[list[int], list[int]]] = {}
    for m in rk.order:
        fp_all = sorted(rk.first_pos[s] for s in rk.member_sites[m])
        fp_common = sorted(rk.first_pos[s] for s in rk.member_sites[m] if klass[s] <= 1)
        per_member[m] = (fp_all, fp_common)
    rows: list[dict[str, object]] = []
    for cp in checkpoints:
        if cp >= len(rk.order):
            continue
        shares_all, shares_common = [], []
        for m in rk.order[cp:]:
            fp_all, fp_common = per_member[m]
            if fp_all:
                shares_all.append(bisect_left(fp_all, cp) / len(fp_all))
            if fp_common:
                shares_common.append(bisect_left(fp_common, cp) / len(fp_common))
        rows.append({"first": cp, "members_outside": len(rk.order) - cp, **_spread("all", shares_all), **_spread("common", shares_common)})
    return rows


def _quantile(sorted_values: list[float], q: float) -> float:
    if not sorted_values:
        return float("nan")
    i = q * (len(sorted_values) - 1)
    lo = int(math.floor(i))
    hi = min(lo + 1, len(sorted_values) - 1)
    return sorted_values[lo] + (sorted_values[hi] - sorted_values[lo]) * (i - lo)


def _spread(name: str, values: list[float]) -> dict[str, object]:
    v = sorted(values)
    return {f"{name}_median": _quantile(v, 0.5), f"{name}_p05": _quantile(v, 0.05), f"{name}_min": v[0] if v else float("nan")}


# ---- the command ------------------------------------------------------------------------------------------


@dataclass
class Member:
    sample_id: str
    eligible: bool = True
    reason: str = ""
    sex: str = ""
    depth_mean: float | None = None
    read_n50: float | None = None
    qc_status: str = ""
    ancestry: str = ""
    pcs: list[float] = field(default_factory=list)


def _f(v: object, digits: int = 4) -> str:
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return ""
    return f"{v:.{digits}f}" if isinstance(v, float) else str(v)


def member_table(digests: dict[str, dict[str, object]], cohort_members: list[str], *, min_depth: float | None,
                 min_read_n50: float | None) -> list[Member]:
    """One row per cohort member with what decides eligibility: a digest, no read-based QC flag of severity fail
    (the per-sample identity flags are among them; the assembly's `asm_*` flags are not: a ranking is made before
    the assemblies exist and must not change when they do, and `--emit --pass-list` applies them), and the optional
    depth and read-length floors."""
    rows: list[Member] = []
    for sid in cohort_members:
        d = digests.get(sid)
        m = Member(sample_id=sid)
        if d is None:
            m.eligible, m.reason = False, "no successful run to summarise"
            rows.append(m)
            continue
        qc = d.get("qc") if isinstance(d.get("qc"), dict) else {}
        cov = d.get("coverage") if isinstance(d.get("coverage"), dict) else {}
        reads = d.get("reads") if isinstance(d.get("reads"), dict) else {}
        ident = d.get("identity") if isinstance(d.get("identity"), dict) else {}
        m.qc_status = str(qc.get("status") or "pass")  # type: ignore[union-attr]
        m.depth_mean = sp._float(cov.get("mean"))  # type: ignore[union-attr]
        m.read_n50 = sp._float(reads.get("read_length_n50"))  # type: ignore[union-attr]
        m.sex = str(ident.get("resolved") or cov.get("inferred_sex") or "")  # type: ignore[union-attr]
        fails = sorted({str(f.get("code")) for f in qc.get("flags", [])  # type: ignore[union-attr]
                        if f.get("severity") == "fail" and not str(f.get("code")).startswith("asm_")})
        if fails:
            m.eligible, m.reason = False, "QC fail: " + ", ".join(fails)
        elif min_depth is not None and (m.depth_mean is None or m.depth_mean < min_depth):
            m.eligible, m.reason = False, f"mean depth {_f(m.depth_mean, 1) or 'unknown'} below {min_depth:g}"
        elif min_read_n50 is not None and (m.read_n50 is None or m.read_n50 < min_read_n50):
            m.eligible, m.reason = False, f"read N50 {_f(m.read_n50, 0) or 'unknown'} below {min_read_n50:g}"
        rows.append(m)
    return rows


def run(cfg: Config, db: DB, code: CodeInfo, *, cohort_id: str, mode: str, version: str | None, out_dir: Path,
        digest_dir: Path, thresholds: dict[str, float], jobs: int = 1, first_file: Path | None = None,
        min_depth: float | None = None, min_read_n50: float | None = None, common_min: float = 0.01,
        related_min: float = 0.1) -> list[Path]:
    """Rank the cohort's members and write `<cohort>.rank.tsv`, `.rank.curves.tsv`, `.rank.heldout.tsv` and
    `.rank.html` into out_dir; returns the paths."""
    from . import summary  # noqa: PLC0415 - summary imports the page helpers this module shares

    cohort = db.get_cohort(cohort_id)
    if cohort is None:
        raise UgcError(f"unknown cohort {cohort_id}")
    if not 0.0 < common_min < 0.05:
        raise UgcError("--common-min must be above 0 and below 0.05 (the 5% class is fixed)")
    cruns = summary.resolve_cohort_runs(db, cohort_id, mode, version)
    merge = cruns.get("cohort_merge")
    if merge is None:
        raise UgcError(f"cohort {cohort_id} has no successful cohort_merge run in {mode} mode"
                       + ("" if version is None else f" at {version} (--any-version reads earlier ones)"))
    vcf = Outputs.load(merge).opt("cohort_sv_vcf")
    if not isinstance(vcf, str) or not Path(vcf).exists():
        raise UgcError(f"cohort_merge run {merge.run_id} has no cohort_sv_vcf on disk")

    collected = summary.collect_digests(cfg, db, code, sample_ids=list(cohort.members), mode=mode, version=version,
                                        out_dir=digest_dir, thresholds=thresholds, jobs=jobs)
    members = member_table(collected.digests, list(cohort.members), min_depth=min_depth, min_read_n50=min_read_n50)
    by_id = {m.sample_id: m for m in members}

    pairs: list[dict[str, object]] = []
    qc = cruns.get("cohort_qc")
    if qc is not None:
        outs = Outputs.load(qc)
        p = outs.opt("qc_related_pairs")
        if isinstance(p, str) and Path(p).exists():
            pairs = sp.parse_somalier_pairs(sp.read_lines(Path(p)))
        p = outs.opt("qc_ancestry")
        if isinstance(p, str) and Path(p).exists():
            for r in sp.parse_ancestry(sp.read_lines(Path(p))):
                m = by_id.get(str(r.get("sample_id")))
                if m is not None and not r.get("given"):
                    m.ancestry = str(r.get("predicted") or "")
                    m.pcs = [float(x) for x in r.get("pcs", []) if x is not None][:5]  # type: ignore[union-attr]
    else:
        diag(f"cohort {cohort_id}: no cohort_qc run in {mode} mode; relatives are not held back and no ancestry is shown")

    rankable = [m.sample_id for m in members if m.eligible]
    if not rankable:
        why: dict[str, int] = {}
        for m in members:
            why[m.reason] = why.get(m.reason, 0) + 1
        raise UgcError(f"cohort {cohort_id}: no eligible member ("
                       + "; ".join(f"{n} x {r}" for r, n in sorted(why.items(), key=lambda kv: -kv[1])[:4])
                       + "); a QC threshold can be moved for the ranking with --threshold KEY=VALUE")
    diag(f"reading {vcf}")
    matrix = read_sv_matrix(Path(vcf), rankable)
    for sid in matrix.missing:
        by_id[sid].eligible, by_id[sid].reason = False, "no column in the cohort SV VCF"
    if matrix.missing:   # a member without a column carries nothing, so dropping it only renumbers the others
        gone = set(matrix.missing)
        keep_idx = [i for i, sid in enumerate(matrix.members) if sid not in gone]
        remap = {old: new for new, old in enumerate(keep_idx)}
        matrix.site_carriers = [array("H", (remap[k] for k in c)) for c in matrix.site_carriers]
        matrix.members = [matrix.members[i] for i in keep_idx]
    index = {sid: i for i, sid in enumerate(matrix.members)}
    related: dict[int, dict[int, float]] = {}
    for pr in pairs:
        a, b, rel = index.get(str(pr.get("sample_a"))), index.get(str(pr.get("sample_b"))), pr.get("relatedness")
        if a is None or b is None or not isinstance(rel, float) or rel < related_min:
            continue
        related.setdefault(a, {})[b] = rel
        related.setdefault(b, {})[a] = rel
    first: list[int] = []
    if first_file is not None:
        for sid in read_ids(first_file):
            if sid not in index:
                why = by_id[sid].reason if sid in by_id else "not a member of the cohort"
                raise UgcError(f"--first: {sid} cannot be ranked ({why})")
            first.append(index[sid])

    klass = site_classes(matrix, common_min)
    rk = greedy_rank(matrix, klass, eligible=list(range(len(matrix.members))), first=first, related=related)
    n = len(rk.order)
    by_count: list[dict[int, int]] = [{} for _ in CLASSES]
    for s, carriers in enumerate(matrix.site_carriers):
        d = by_count[klass[s]]
        d[len(carriers)] = d.get(len(carriers), 0) + 1
    random_rows = expected_random(by_count, n)
    heldout = held_out(rk, klass)

    out_dir.mkdir(parents=True, exist_ok=True)
    totals = rk.totals
    total_all = sum(totals) or 1

    def frac(covered: int, total: int) -> str:
        return f"{covered / total:.4f}" if total else ""

    rank_path = out_dir / f"{cohort_id}.rank.tsv"
    lines = ["\t".join(RANK_COLUMNS)]
    for pos, mi in enumerate(rk.order):
        m = by_id[matrix.members[mi]]
        rel = rk.related_to.get(mi)
        cum = rk.cum[pos]
        lines.append("\t".join([str(pos + 1), m.sample_id, rk.tier[mi], matrix.members[rel[0]] if rel else "", _f(rel[1], 3) if rel else "",
                                str(rk.new_common[pos]), str(rk.new_all[pos]),
                                *[frac(cum[k], totals[k]) for k in range(4)], frac(sum(cum), total_all),
                                str(len(rk.member_sites[mi])), m.sex, _f(m.depth_mean, 1), _f(m.read_n50, 0), m.qc_status, m.ancestry,
                                *[_f(m.pcs[i], 4) if i < len(m.pcs) else "" for i in range(5)], ""]))
    for m in members:
        if not m.eligible:
            lines.append("\t".join(["", m.sample_id, "excluded", "", "", "", "", "", "", "", "", "", "", m.sex, _f(m.depth_mean, 1),
                                    _f(m.read_n50, 0), m.qc_status, m.ancestry, *[_f(m.pcs[i], 4) if i < len(m.pcs) else "" for i in range(5)],
                                    m.reason]))
    rank_path.write_text("\n".join(lines) + "\n")

    curves_path = out_dir / f"{cohort_id}.rank.curves.tsv"
    head = ["n"] + [f"greedy_{c}" for c in CLASSES] + ["greedy_all"] + [f"random_{c}" for c in CLASSES] + ["random_all"]
    clines = ["\t".join(head)]
    for k in range(n):
        g, r = rk.cum[k], random_rows[k]
        clines.append("\t".join([str(k + 1), *[frac(g[c], totals[c]) for c in range(4)], frac(sum(g), total_all),
                                 *[(f"{r[c] / totals[c]:.4f}" if totals[c] else "") for c in range(4)], f"{sum(r) / total_all:.4f}"]))
    curves_path.write_text("\n".join(clines) + "\n")

    held_path = out_dir / f"{cohort_id}.rank.heldout.tsv"
    hcols = ["first", "members_outside", "all_median", "all_p05", "all_min", "common_median", "common_p05", "common_min"]
    held_path.write_text("\n".join(["\t".join(hcols)] + ["\t".join(_f(r[c]) for c in hcols) for r in heldout]) + "\n")

    html_path = out_dir / f"{cohort_id}.rank.html"
    html_path.write_text(render(cohort_id, mode, matrix, rk, members, by_id, random_rows, heldout, common_min=common_min,
                                related_min=related_min, vcf=str(vcf)))
    tiers = {t: sum(1 for v in rk.tier.values() if v == t) for t in ("first", "unrelated", "related")}
    diag(f"{cohort_id}: ranked {n} of {len(members)} members ({tiers['first']} placed first, {tiers['unrelated']} unrelated, "
         f"{tiers['related']} held back as relatives; {len(members) - n} excluded) on {len(matrix.site_carriers)} SV sites")
    return [rank_path, curves_path, held_path, html_path]


def read_rank(path: Path) -> list[dict[str, str]]:
    _cols, rows = sp.table(sp.read_lines(path))
    if not rows or "sample_id" not in rows[0] or "rank" not in rows[0]:
        raise UgcError(f"{path}: not a ranking written by `cohort rank`")
    return rows


def emit(rank_path: Path, k: int, *, pass_list: Path | None, out: Path, include_related: bool = False) -> list[str]:
    """The first `k` members of a saved ranking, in rank order, that are in `pass_list` (when given): the member
    list of one ladder cohort. Held-back relatives are skipped unless `include_related`; a list of the whole
    ranking (k = 0) always includes them."""
    if k < 0:
        raise UgcError("--emit takes a count of members, or 0 for the whole ranking")
    passing = set(read_ids(pass_list)) if pass_list is not None else None
    ids: list[str] = []
    for r in read_rank(rank_path):
        if not r.get("rank"):
            continue
        if k and r.get("tier") == "related" and not include_related:
            continue
        if passing is not None and r["sample_id"] not in passing:
            continue
        ids.append(r["sample_id"])
        if k and len(ids) == k:
            break
    if k and len(ids) < k:
        raise UgcError(f"{rank_path}: only {len(ids)} ranked member(s) qualify, {k} asked for"
                       + ("" if include_related else " (held-back relatives are left out; --include-related adds them)"))
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("".join(i + "\n" for i in ids))
    return ids


def render(cohort_id: str, mode: str, matrix: Matrix, rk: Ranking, members: list[Member], by_id: dict[str, Member],
           random_rows: list[list[float]], heldout: list[dict[str, object]], *, common_min: float, related_min: float, vcf: str) -> str:
    from . import summary  # noqa: PLC0415

    n = len(rk.order)
    totals = rk.totals
    pct = f"{common_min * 100:g}%"

    def series(k: int) -> tuple[list[tuple[float, float]], list[tuple[float, float]]]:
        if not totals[k]:
            return [], []
        return ([(i + 1, rk.cum[i][k] / totals[k]) for i in range(n)], [(i + 1, random_rows[i][k] / totals[k]) for i in range(n)])

    charts = []
    for k, name in ((0, "ge5"), (1, "common"), (2, "rare")):
        greedy, rand = series(k)
        if not greedy:
            continue
        charts.append(f"<h3>{_esc(CLASS_LABELS[name])}: {totals[k]} SV sites</h3>"
                      + sv.lines([("this ranking", "ln2", greedy), ("members at random (expected)", "ln3", rand)],
                                 x_label="members in the graph", y_label="share of the sites carried by at least one of them", percent=True))
    cards = summary._cards([("members ranked", n), ("excluded", len(members) - n),
                            ("held back as relatives", sum(1 for v in rk.tier.values() if v == "related")),
                            ("SV sites with a carrier", len(matrix.site_carriers)),
                            (f"sites at {pct} carrier frequency or more", totals[0] + totals[1])])
    intro = (f"<p>Members of <b>{_esc(cohort_id)}</b> ordered by the structural variants each adds to the ones before it, from the "
             f"{_esc(mode)} cohort_merge SV VCF ({', '.join(f'{v} {k}' for k, v in sorted(matrix.by_type.items()))}). At every step the member "
             f"with the most not-yet-covered alleles of {pct} carrier frequency or more is taken, ties by all alleles. A member related to "
             f"an earlier pick (relatedness {related_min:g} or more in cohort_qc) is held back to the end, so the first members are "
             "unrelated. Frequencies are carrier frequencies among the ranked members.</p>")
    held_rows = [[r["first"], r["members_outside"], *[f"{float(r[c]) * 100:.1f}%" for c in ("common_median", "common_p05", "common_min",  # type: ignore[arg-type]
                                                                                          "all_median", "all_p05", "all_min")]] for r in heldout]
    held_html = (_table(["first members", "members outside", "common alleles already there: median", "5th percentile", "lowest",
                         "all alleles: median", "5th percentile", "lowest"], held_rows)  # type: ignore[arg-type]
                 if held_rows else "<p class='muted'>The ranking is shorter than the first checkpoint (20 members).</p>")
    rows = []
    for pos, mi in enumerate(rk.order):
        m = by_id[matrix.members[mi]]
        rel = rk.related_to.get(mi)
        cum = rk.cum[pos]
        rows.append([pos + 1, m.sample_id, rk.tier[mi] + (f" to {matrix.members[rel[0]]} ({rel[1]:.2f})" if rel else ""),
                     rk.new_common[pos], rk.new_all[pos],
                     f"{(cum[0] + cum[1]) / ((totals[0] + totals[1]) or 1) * 100:.1f}%", f"{sum(cum) / (sum(totals) or 1) * 100:.1f}%",
                     m.sex, _f(m.depth_mean, 1), _f(m.read_n50, 0), m.ancestry])
    table = big_table(["rank", "member", "tier", "new common alleles", "new alleles", "common sites covered", "all sites covered",
                       "sex", "mean depth", "read N50", "ancestry"], rows, what="members")  # type: ignore[arg-type]
    excluded = [[m.sample_id, m.reason] for m in members if not m.eligible]
    pcs = [(m.sample_id, m.pcs[0], m.pcs[1], m.ancestry) for m in members if len(m.pcs) >= 2]
    top = {matrix.members[mi] for mi in rk.order[:100]}
    scatter = ""
    if pcs:
        scatter = ("<p>The first 100 members (large points) on the cohort's ancestry components (all ranked members, small points): "
                   "a greedy cover tends to reach for the edges of the cloud.</p>"
                   + sv.scatter([(lab or "unlabelled", x, y) for sid, x, y, lab in pcs if sid not in top],
                                [(sid, x, y, lab or "unlabelled") for sid, x, y, lab in pcs if sid in top]))
    body = (f"<h1>{_esc(cohort_id)}: members ranked for a pangenome graph</h1>" + cards + intro
            + summary._sec("sec-curves", "How much the first members hold", "".join(charts) or "<p class='muted'>No data.</p>")
            + summary._sec("sec-heldout", "What a member outside already finds in the graph", held_html)
            + (summary._sec("sec-ancestry", "Where the first picks sit", scatter) if scatter else "")
            + summary._sec("sec-rank", "The ranking", table)
            + (summary._sec("sec-excluded", "Not ranked", _table(["member", "why"], excluded)) if excluded else "")  # type: ignore[arg-type]
            + f"<p class='muted'>Source: {_esc(vcf)}</p>")
    nav = summary._nav(cohort_id, [("sec-curves", "Curves"), ("sec-heldout", "Outside members")]
                       + ([("sec-ancestry", "Ancestry")] if scatter else []) + [("sec-rank", "Ranking")]
                       + ([("sec-excluded", "Not ranked")] if excluded else []))
    return summary._page(f"{cohort_id} ranking", nav, body)


