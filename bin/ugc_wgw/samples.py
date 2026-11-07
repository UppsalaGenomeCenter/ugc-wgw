"""Sample registration from a TSV, input BAM checks and per-stage sample listing (docs/DESIGN.md §8.1)."""
from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Callable, Optional

from . import bam, vcf
from .db import DB, SINGLE_INPUT_KINDS, SampleRecord
from .log import LOGGER, Events
from .outputs import Outputs
from .stages import ENTRY_KINDS, mode_stages, stage_spec
from .util import UgcError, utc_now

ENTRY_COLUMNS = ("aligned_bam", "aligned_bam_index", "small_variant_vcf", "small_variant_vcf_index", "sv_vcf", "sv_vcf_index",
                 "gvcf", "gvcf_index")   # the bam and vcf entry kinds (docs/guide/05-project-setup.md)
COLUMNS = ("sample_id", "sex", "hifi_reads", "fail_reads", "father_id", "mother_id") + ENTRY_COLUMNS
REQUIRED_COLUMNS = ("sample_id",)
INDEX_SUFFIXES = {"bam": (".bai", ".csi"), "vcf": (".tbi", ".csi")}
SEX_VALUES = ("MALE", "FEMALE")
SEX_ALIASES = {"MALE": "MALE", "M": "MALE", "XY": "MALE", "1": "MALE", "FEMALE": "FEMALE", "F": "FEMALE", "XX": "FEMALE", "2": "FEMALE"}
ID_RE = re.compile(r"^[A-Za-z0-9._-]+$")
# Advisory floors for the input check (config.json `input_thresholds`, guide chapter 05); 0 disables one.
INPUT_THRESHOLDS: dict[str, float] = {
    "file_reads_min": 1000.0,     # HiFi reads per BAM below which the file is flagged (a Revio BAM has millions)
    "sample_gbases_min": 30.0,    # HiFi bases per sample (Gb) below which the sample is flagged (~10x of GRCh38)
}
DROPPABLE = ("empty file", "no reads")   # the problems `--drop-empty` may remove a file for; anything else is damage
Progress = Callable[[int, int, str], None]   # (files done, files total, label) after every inspected file


def input_thresholds(overrides: dict[str, object] | None) -> dict[str, float]:
    """INPUT_THRESHOLDS with the project's overrides; unknown keys and non-numbers are errors."""
    out = dict(INPUT_THRESHOLDS)
    for key, value in (overrides or {}).items():
        if key not in INPUT_THRESHOLDS:
            raise UgcError(f"input_thresholds: unknown key {key!r} (known: {', '.join(INPUT_THRESHOLDS)})")
        if isinstance(value, bool) or not isinstance(value, (int, float)) or value < 0:
            raise UgcError(f"input_thresholds.{key}: must be a number >= 0, not {value!r}")
        out[key] = float(value)
    return out


def _split_paths(cell: str | None) -> list[str]:
    if not cell or not cell.strip():
        return []
    return [p.strip() for p in cell.split(",") if p.strip()]


def parse_tsv(path: Path) -> list[SampleRecord]:
    """Header row required. Multi-file cells are comma-separated. Blank cells are null. Extra columns go to meta."""
    added_at = utc_now()
    records: list[SampleRecord] = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        header = reader.fieldnames or []
        missing = [c for c in REQUIRED_COLUMNS if c not in header]
        if missing:
            raise UgcError(f"{path}: header lacks required column(s): {', '.join(missing)}")
        if "hifi_reads" not in header and "aligned_bam" not in header:
            raise UgcError(f"{path}: header lacks an input column: hifi_reads (unaligned reads) or aligned_bam (an aligned BAM, with or without VCFs)")
        for n, row in enumerate(reader, start=2):
            sid = (row.get("sample_id") or "").strip()
            if not sid:
                continue
            sex_raw = (row.get("sex") or "").strip().upper()
            sex_raw = SEX_ALIASES.get(sex_raw, sex_raw)   # M/F, XX/XY and PLINK's 1/2 are accepted spellings
            meta = {k: v for k, v in row.items() if k and k not in COLUMNS and v not in (None, "")}
            rec = SampleRecord(
                sample_id=sid,
                sex=sex_raw or None,
                father_id=(row.get("father_id") or "").strip() or None,
                mother_id=(row.get("mother_id") or "").strip() or None,
                added_at=added_at,
                hifi_reads=_split_paths(row.get("hifi_reads")),
                fail_reads=_split_paths(row.get("fail_reads")),
                meta={"tsv_line": n, **meta},
            )
            for col in ENTRY_COLUMNS:
                setattr(rec, col, (row.get(col) or "").strip() or None)
            rec.entry_kind = derive_kind(rec)[0] or "reads"
            records.append(rec)
    if not records:
        raise UgcError(f"{path}: no samples")
    return records


def derive_kind(rec: SampleRecord) -> tuple[str | None, list[str]]:
    """(entry kind, problems) from which input columns the record fills: reads (hifi_reads), bam (aligned_bam
    alone), vcf (aligned_bam with small_variant_vcf and sv_vcf; gvcf optional). One source of truth per sample:
    reads with an aligned BAM, VCFs without the BAM, one VCF without the other, a gVCF without the VCFs and fail
    reads with a non-reads kind are refused (nallo's entry-point rules made explicit)."""
    problems: list[str] = []
    has_reads, has_bam = bool(rec.hifi_reads), bool(rec.aligned_bam)
    vcfs = [k for k in ("small_variant_vcf", "sv_vcf") if getattr(rec, k)]
    if has_reads and has_bam:
        problems.append("hifi_reads and aligned_bam are both given; one source of truth per sample (reads are aligned here; an aligned BAM is used as is)")
    if (vcfs or rec.gvcf) and not has_bam:
        problems.append("VCFs need the aligned_bam they were called from (phasing, repeats and methylation read it)")
    if len(vcfs) == 1:
        problems.append(f"{vcfs[0]} without {'sv_vcf' if vcfs[0] == 'small_variant_vcf' else 'small_variant_vcf'}: both VCFs or neither")
    if rec.gvcf and len(vcfs) < 2:
        problems.append("gvcf given without small_variant_vcf and sv_vcf")
    if not has_reads and not has_bam:
        problems.append("no input: give hifi_reads, or aligned_bam (with small_variant_vcf and sv_vcf for a pre-called sample)")
    if problems:
        return None, problems
    if has_reads:
        kind = "reads"
    elif len(vcfs) == 2:
        kind = "vcf"
    else:
        kind = "bam"
    if rec.fail_reads and kind != "reads":
        problems.append("fail_reads are only taken with hifi_reads (they are aligned and used at the TRGT fail-reads loci)")
    for col in ENTRY_COLUMNS:
        if col.endswith("_index") and getattr(rec, col) and not getattr(rec, col[:-6]):
            problems.append(f"{col} without {col[:-6]}")
    return (kind if not problems else None), problems


def resolve_index(path: str, given: str | None, kind: str, check_paths: bool) -> tuple[str | None, str | None]:
    """(index path, problem) for a BAM or VCF: the given one, else by convention (<file>.bai/.csi or <file>.tbi/.csi,
    and <stem>.bai for a BAM). Without check_paths the conventional name is taken on trust."""
    candidates = [given] if given else [path + ext for ext in INDEX_SUFFIXES[kind]]
    if kind == "bam" and not given and path.endswith(".bam"):
        candidates.append(path[:-4] + ".bai")
    if not check_paths:
        return candidates[0], None
    for c in candidates:
        if Path(c).expanduser().exists():
            return str(Path(c).expanduser().resolve()), None
    if given:
        return None, f"index not found: {given}"
    return None, f"no index next to {Path(path).name} (looked for {', '.join(Path(c).name for c in candidates)}); give the _index column"


def validate(records: list[SampleRecord], db: DB, check_paths: bool = True, *, inspect: bool = True,
             thresholds: dict[str, float] | None = None, drop_empty: bool = False,
             progress: Optional[Progress] = None) -> list[str]:
    """Return warnings; raise UgcError with every hard problem found. With check_paths and inspect, every read
    file is opened (`inspect_inputs`): a truncated, empty or read-less BAM is a hard problem unless drop_empty
    removes the empty ones from the record."""
    problems: list[str] = []
    warnings: list[str] = []
    seen: set[str] = set()
    known = set(db.list_sample_ids()) | {r.sample_id for r in records}
    for rec in records:
        if not ID_RE.match(rec.sample_id):
            problems.append(f"{rec.sample_id}: sample_id must match [A-Za-z0-9._-]+")
        if rec.sample_id in seen:
            problems.append(f"{rec.sample_id}: duplicated in the TSV")
        seen.add(rec.sample_id)
        if rec.sex is not None and rec.sex not in SEX_VALUES:
            problems.append(f"{rec.sample_id}: sex must be MALE or FEMALE (also M/F, XY/XX, 1/2) or blank, not {rec.sex!r}")
        kind, kind_problems = derive_kind(rec)
        problems.extend(f"{rec.sample_id}: {p}" for p in kind_problems)
        rec.entry_kind = kind or rec.entry_kind
        for col in ("aligned_bam", "small_variant_vcf", "sv_vcf", "gvcf"):
            p = getattr(rec, col)
            if not p:
                continue
            pp = Path(p).expanduser()
            if check_paths and not pp.exists():
                problems.append(f"{rec.sample_id}: {col} path not found: {p}")
            setattr(rec, col, str(pp.resolve()) if pp.exists() else str(pp.absolute()))
            if not (check_paths and not pp.exists()):
                idx, problem = resolve_index(getattr(rec, col), getattr(rec, col + "_index"), "bam" if col == "aligned_bam" else "vcf", check_paths)
                if problem:
                    problems.append(f"{rec.sample_id}: {col}: {problem}")
                setattr(rec, col + "_index", idx)
        for kind in ("hifi_reads", "fail_reads"):
            paths = getattr(rec, kind)
            resolved = []
            for p in paths:
                pp = Path(p).expanduser()
                if check_paths and not pp.exists():
                    problems.append(f"{rec.sample_id}: {kind} path not found: {p}")
                resolved.append(str(pp.resolve()) if pp.exists() else str(pp.absolute()))
            setattr(rec, kind, resolved)
        for parent in (rec.father_id, rec.mother_id):
            if parent and parent not in known:
                warnings.append(f"{rec.sample_id}: parent {parent} is not a registered sample")
    if check_paths and inspect and not problems:
        more, warned = inspect_inputs(records, input_thresholds(thresholds), drop_empty=drop_empty, progress=progress)
        problems.extend(more)
        warnings.extend(warned)
    if problems:
        raise UgcError("invalid samples TSV:\n  " + "\n  ".join(problems))
    return warnings


def _file_findings(sid: str, kind: str, info: bam.BamInfo, thresholds: dict[str, float]) -> list[str]:
    """Advisory findings about one usable file, prefixed for the sample."""
    name = Path(info.path).name
    out = [f"{sid}: {kind} {name}: {w}" for w in info.warnings]
    if kind == "hifi_reads" and thresholds["file_reads_min"] and info.reads < thresholds["file_reads_min"]:
        approx = "" if info.reads_exact else "about "
        out.append(f"{sid}: {kind} {name}: only {approx}{info.reads:,} reads ({bam.fmt_bases(info.bases)}), "
                   f"below file_reads_min {thresholds['file_reads_min']:g}")
    return out


def inspect_inputs(records: list[SampleRecord], thresholds: dict[str, float], *, drop_empty: bool = False,
                   progress: Optional[Progress] = None) -> tuple[list[str], list[str]]:
    """Open every read file of every record (`bam.inspect`). Returns (problems, warnings); fills
    `rec.input_info` per kept path and `rec.meta["input_check"]` per sample. With drop_empty, files that are
    empty or have no reads are removed from the record (listed in the warnings and in meta) instead of being
    problems; a truncated or unreadable file is always a problem."""
    problems: list[str] = []
    warnings: list[str] = []
    n_files = sum(len(r.data_inputs()) for r in records)
    done = 0
    for rec in records:
        total_reads = total_bases = 0
        exact = True
        dropped: list[dict[str, object]] = []
        if rec.entry_kind != "reads":
            done += _inspect_entry_files(rec, problems, warnings, progress, done, n_files)
            continue
        for kind in ("hifi_reads", "fail_reads"):
            kept: list[str] = []
            for p in getattr(rec, kind):
                info = bam.inspect(Path(p))
                done += 1
                if progress:
                    progress(done, n_files, f"{rec.sample_id} {Path(p).name}")
                if info.problems:
                    if drop_empty and all(pr.startswith(DROPPABLE) for pr in info.problems):
                        dropped.append({"kind": kind, "path": p, "reason": info.problems[0]})
                        warnings.append(f"{rec.sample_id}: dropped {kind} {Path(p).name}: {info.problems[0]}")
                        continue
                    problems.extend(f"{rec.sample_id}: {kind} {p}: {pr}" for pr in info.problems)
                    kept.append(p)
                    continue
                kept.append(p)
                rec.input_info[p] = info.to_dict()
                warnings.extend(_file_findings(rec.sample_id, kind, info, thresholds))
                if kind == "hifi_reads":
                    total_reads += info.reads
                    total_bases += info.bases
                    exact = exact and info.reads_exact
            setattr(rec, kind, kept)
        if not rec.hifi_reads:
            problems.append(f"{rec.sample_id}: no hifi_reads left after dropping empty files")
            continue
        floor = thresholds["sample_gbases_min"]
        if floor and total_bases / 1e9 < floor and not any(pr.startswith(rec.sample_id + ":") for pr in problems):
            cov = total_bases / 1e9 / bam.HUMAN_GENOME_GB
            warnings.append(f"{rec.sample_id}: {bam.fmt_bases(total_bases)} of HiFi bases in {len(rec.hifi_reads)} "
                            f"file(s), about {cov:.1f}x of GRCh38, below sample_gbases_min {floor:g}")
        LOGGER.info("inspected %s: %d read file(s), %s of HiFi bases (about %.1fx)%s", rec.sample_id,
                    len(rec.hifi_reads) + len(rec.fail_reads), bam.fmt_bases(total_bases),
                    total_bases / 1e9 / bam.HUMAN_GENOME_GB, "" if exact else ", estimated")
        rec.meta["input_check"] = {
            "checked_at": utc_now(), "hifi_files": len(rec.hifi_reads), "fail_files": len(rec.fail_reads),
            "reads": total_reads, "gbases": round(total_bases / 1e9, 3),
            "coverage": round(total_bases / 1e9 / bam.HUMAN_GENOME_GB, 2), "exact": exact,
            "dropped": dropped,
        }
    return problems, warnings


def _inspect_entry_files(rec: SampleRecord, problems: list[str], warnings: list[str], progress: Optional[Progress],
                         done: int, n_files: int) -> int:
    """The bam and vcf entry kinds: the aligned BAM (must be aligned, i.e. have @SQ lines) and every VCF (BGZF, one
    sample). Fills rec.input_info and rec.meta["input_check"]; returns the number of files looked at."""
    seen = 0
    check: dict[str, object] = {"checked_at": utc_now(), "entry_kind": rec.entry_kind}
    for kind, p in rec.data_inputs():
        seen += 1
        if progress:
            progress(done + seen, n_files, f"{rec.sample_id} {Path(p).name}")
        if kind == "aligned_bam":
            info = bam.inspect(Path(p))
            if info.ok and not info.aligned:
                info.problems.append("not aligned (no @SQ lines in the header): an aligned BAM is expected under aligned_bam; unaligned reads go under hifi_reads")
            info.warnings = [w for w in info.warnings if not w.startswith("aligned input")]
            if info.problems:
                problems.extend(f"{rec.sample_id}: {kind} {p}: {pr}" for pr in info.problems)
                continue
            rec.input_info[p] = info.to_dict()
            warnings.extend(f"{rec.sample_id}: {kind} {Path(p).name}: {w}" for w in info.warnings)
            check["bam"] = {"reads": info.reads, "reads_exact": info.reads_exact, "gbases": round(info.gb, 3), "coverage": round(info.coverage, 2),
                            "n_ref": info.n_ref, "samples": info.samples, "movies": info.movies}
            # the callers name their outputs after the read group (sawfish its sample directory and VCF column, TRGT its
            # VCF column; HiPhase is told the sample id): an aligned BAM's read-group sample name must be the sample id
            if not info.samples:
                problems.append(f"{rec.sample_id}: aligned_bam {p}: no @RG SM in the header; the callers name their outputs after it "
                                "(samtools addreplacerg adds a read group)")
            elif info.samples != [rec.sample_id]:
                problems.append(f"{rec.sample_id}: aligned_bam {p}: @RG SM is {', '.join(info.samples)}, not the sample id: register the "
                                "sample under that name, or rewrite the header first (samtools reheader)")
        else:
            vinfo = vcf.inspect(Path(p))
            if vinfo.problems:
                problems.extend(f"{rec.sample_id}: {kind} {p}: {pr}" for pr in vinfo.problems)
                continue
            rec.input_info[p] = vinfo.to_dict()
            warnings.extend(f"{rec.sample_id}: {kind} {Path(p).name}: {w}" for w in vinfo.warnings)
            check.setdefault("vcfs", {})[kind] = {"samples": vinfo.samples, "n_contigs": vinfo.n_contigs, "fileformat": vinfo.fileformat, "sources": vinfo.sources[:5]}  # type: ignore[index]
            if vinfo.samples and vinfo.samples[0] != rec.sample_id:
                warnings.append(f"{rec.sample_id}: {kind} {Path(p).name}: sample column {vinfo.samples[0]} is renamed to the sample id by the intake stage")
    if not any(pr.startswith(rec.sample_id + ":") for pr in problems):
        LOGGER.info("inspected %s (%s entry): %d file(s)", rec.sample_id, rec.entry_kind, seen)
        rec.meta["input_check"] = check
    return seen


def add_samples(db: DB, events: Events, path: Path, check_paths: bool = True, replace: bool = False, *,
                inspect: bool = True, thresholds: dict[str, object] | None = None, drop_empty: bool = False,
                progress: Optional[Progress] = None) -> int:
    records = parse_tsv(path)
    for w in validate(records, db, check_paths, inspect=inspect, thresholds=thresholds, drop_empty=drop_empty,
                      progress=progress):
        events.emit("sample.warning", level="warning", message=w)
    if not replace:
        dup = [r.sample_id for r in records if db.sample_exists(r.sample_id)]
        if dup:
            raise UgcError("already registered (use --replace): " + ", ".join(dup))
    for rec in records:
        db.upsert_sample(rec, replace=replace)
        check = rec.meta.get("input_check") or {}
        for d in check.get("dropped", []):   # type: ignore[union-attr]
            events.emit("sample.input_dropped", level="warning", sample_id=rec.sample_id, **d)   # type: ignore[arg-type]
        events.emit("sample.added", sample_id=rec.sample_id, entry_kind=rec.entry_kind, hifi_reads=len(rec.hifi_reads),
                    fail_reads=len(rec.fail_reads), replaced=replace,
                    gbases=check.get("gbases"), reads=check.get("reads"))   # type: ignore[union-attr]
    return len(records)


def check_samples(db: DB, ids: list[str], thresholds: dict[str, object] | None, *, stored: bool = False,
                  update: bool = True, progress: Optional[Progress] = None
                  ) -> tuple[list[dict[str, object]], list[str], list[str]]:
    """Inspect the registered read files of `ids` (all samples when empty) again, or show what registration
    recorded (`stored`). Returns (one row per file, problems, warnings); with update the stored info is refreshed."""
    thr = input_thresholds(thresholds)
    rows: list[dict[str, object]] = []
    problems: list[str] = []
    warnings: list[str] = []
    sample_ids = ids or db.list_sample_ids()
    unknown = [sid for sid in sample_ids if not db.sample_exists(sid)]
    if unknown:
        raise UgcError("unknown sample(s): " + ", ".join(unknown))
    recs = {sid: db.get_sample(sid) for sid in sample_ids}
    n_files = sum(len(r.data_inputs()) for r in recs.values() if r is not None)
    done = 0
    for sid in sample_ids:
        rec = recs[sid]
        assert rec is not None
        total_bases = 0
        for kind, p in rec.data_inputs():
            if True:
                done += 1
                if progress and not stored:
                    progress(done, n_files, f"{sid} {Path(p).name}")
                if kind in ("small_variant_vcf", "sv_vcf", "gvcf"):
                    if stored:
                        vd = rec.input_info.get(p)
                        vinfo = vcf.VcfInfo.from_dict(vd) if vd else None
                    else:
                        vinfo = vcf.inspect(Path(p))
                        if update and vinfo.ok:
                            db.update_input_info(sid, p, vinfo.to_dict())
                    if vinfo is None:
                        rows.append({"sample_id": sid, "kind": kind, "file": Path(p).name, "status": "not inspected"})
                        continue
                    if vinfo.problems:
                        problems.extend(f"{sid}: {kind} {p}: {pr}" for pr in vinfo.problems)
                    rows.append({"sample_id": sid, "kind": kind, "file": Path(p).name, "size": bam.fmt_bytes(vinfo.bytes),
                                 "reads": "", "bases": "", "movie": ",".join(vinfo.samples),
                                 "status": ("problem: " + "; ".join(vinfo.problems)) if vinfo.problems else ("warning: " + "; ".join(vinfo.warnings) if vinfo.warnings else "ok")})
                    warnings.extend(f"{sid}: {kind} {Path(p).name}: {w}" for w in vinfo.warnings)
                    continue
                if stored:
                    info_d = rec.input_info.get(p)
                    if not info_d:
                        rows.append({"sample_id": sid, "kind": kind, "file": Path(p).name, "status": "not inspected"})
                        continue
                    info = bam.BamInfo.from_dict(info_d)
                else:
                    info = bam.inspect(Path(p))
                    if kind == "aligned_bam":
                        info.warnings = [w for w in info.warnings if not w.startswith("aligned input")]
                    if update and info.ok:
                        db.update_input_info(sid, p, info.to_dict())
                findings = [] if info.problems else _file_findings(sid, kind, info, thr)
                if info.problems:
                    problems.extend(f"{sid}: {kind} {p}: {pr}" for pr in info.problems)
                    status = "problem: " + "; ".join(info.problems)
                elif findings:
                    warnings.extend(findings)
                    status = "warning: " + "; ".join(f.split(": ", 2)[-1] for f in findings)
                else:
                    status = "ok"
                if kind == "hifi_reads" and not info.problems:
                    total_bases += info.bases
                rows.append({"sample_id": sid, "kind": kind, "file": Path(p).name, "size": bam.fmt_bytes(info.bytes),
                             "reads": f"{'' if info.reads_exact else '~'}{info.reads:,}" if info.reads else "",
                             "bases": bam.fmt_bases(info.bases) if info.bases else "",
                             "movie": ",".join(info.movies), "status": status})
        floor = thr["sample_gbases_min"]
        if not stored and floor and rec.hifi_reads and total_bases / 1e9 < floor \
                and not any(pr.startswith(sid + ":") for pr in problems):
            warnings.append(f"{sid}: {bam.fmt_bases(total_bases)} of HiFi bases, about "
                            f"{total_bases / 1e9 / bam.HUMAN_GENOME_GB:.1f}x of GRCh38, below sample_gbases_min {floor:g}")
    return rows, problems, warnings


def remove_samples(db: DB, events: Events, ids: list[str], *, force: bool, results_dir: Path) -> list[dict[str, object]]:
    """Unregister samples: every id is validated before the first deletion. Results on disk are never touched."""
    problems: list[str] = []
    for sid in ids:
        if not db.sample_exists(sid):
            problems.append(f"{sid}: unknown sample")
            continue
        cohorts = db.sample_cohorts(sid)
        if cohorts:
            problems.append(f"{sid}: member of cohort(s) {', '.join(cohorts)}; cohorts are immutable")
        runs = db.runs_for("sample", sid)
        if any(r.is_active for r in runs):
            problems.append(f"{sid}: has active run(s); stop or reconcile them first")
        elif runs and not force:
            problems.append(f"{sid}: has {len(runs)} recorded run(s); use --force to delete the rows (results on disk are kept)")
    if problems:
        raise UgcError("cannot remove:\n  " + "\n  ".join(problems))
    out = []
    for sid in ids:
        counts = db.delete_sample(sid, with_runs=force)
        kept = str(Path(results_dir) / "samples" / sid)
        events.emit("sample.removed", sample_id=sid, runs_deleted=counts["runs"], results_kept=kept)
        out.append({"sample_id": sid, "runs_deleted": counts["runs"], "results_kept": kept})
    return out


def list_samples(db: DB, mode: str, ugc_wgw_version: str | None, stage: str | None = None,
                 statuses: list[str] | None = None) -> list[dict[str, object]]:
    """Rows: sample_id, sex, and the latest status per sample stage of the mode ('-' when never run)."""
    stages = [s for s in mode_stages(mode) if stage_spec(s).subject_type == "sample"]
    if stage:
        if stage not in stages:
            raise UgcError(f"stage {stage} is not a sample stage of mode {mode}")
        stages = [stage]
    latest = db.latest_per_stage(ugc_wgw_version)
    rows: list[dict[str, object]] = []
    for rec in db.list_samples():
        row: dict[str, object] = {"sample_id": rec.sample_id, "sex": rec.sex or "", "entry_kind": rec.entry_kind}
        for s in stages:
            run = latest.get(("sample", rec.sample_id, s))
            row[s] = run.status if run else "-"
        if statuses and not any(row[s] in statuses for s in stages):
            continue
        rows.append(row)
    return rows


REBASE_STAGES = ("singleton", "upstream", "intake")


def rebase_sample(db: DB, events: Events, sample_id: str, run_id: str, *, with_vcfs: bool = False) -> SampleRecord:
    """Register a member's own earlier outputs as its inputs: the aligned BAM of a successful singleton, upstream or
    intake run (singleton: the haplotagged BAM), and with with_vcfs its small-variant and SV VCFs (and gVCF when the
    run has one). The entry kind becomes bam or vcf; the reads stay recorded in meta["rebase"]. The next stage run
    reuses the files by path, so a tool bump after alignment costs one downstream per member."""
    rec = db.get_sample(sample_id)
    if rec is None:
        raise UgcError(f"unknown sample {sample_id}")
    run = db.get_run(run_id)
    if run is None or run.subject_type != "sample" or run.subject_id != sample_id:
        raise UgcError(f"{run_id} is not a run of sample {sample_id}")
    if run.status != "success" or run.stage not in REBASE_STAGES:
        raise UgcError(f"{run_id} is a {run.status} {run.stage} run; rebase takes a successful {', '.join(REBASE_STAGES)} run")
    outs = Outputs.load(run)
    bam_path = outs.opt("aligned_hifi_reads") or outs.opt("merged_haplotagged_bam")
    bam_index = outs.opt("aligned_hifi_reads_index") or outs.opt("merged_haplotagged_bam_index")
    if not isinstance(bam_path, str) or not isinstance(bam_index, str):
        raise UgcError(f"{run_id} has no aligned BAM output")
    new = SampleRecord(sample_id=rec.sample_id, sex=rec.sex, father_id=rec.father_id, mother_id=rec.mother_id, added_at=utc_now(),
                       meta={k: v for k, v in rec.meta.items() if k != "input_check"})
    new.aligned_bam, new.aligned_bam_index = bam_path, bam_index
    if with_vcfs:
        small = outs.opt("small_variant_vcf") or outs.opt("phased_small_variant_vcf")
        small_i = outs.opt("small_variant_vcf_index") or outs.opt("phased_small_variant_vcf_index")
        sv = outs.opt("sv_vcf") or outs.opt("phased_sv_vcf")
        sv_i = outs.opt("sv_vcf_index") or outs.opt("phased_sv_vcf_index")
        if not all(isinstance(x, str) for x in (small, small_i, sv, sv_i)):
            raise UgcError(f"{run_id} has no small-variant and SV VCF outputs to rebase on")
        new.small_variant_vcf, new.small_variant_vcf_index, new.sv_vcf, new.sv_vcf_index = small, small_i, sv, sv_i  # type: ignore[assignment]
        g, gi = outs.opt("small_variant_gvcf"), outs.opt("small_variant_gvcf_index")
        if isinstance(g, str) and isinstance(gi, str):
            new.gvcf, new.gvcf_index = g, gi
    kind, problems = derive_kind(new)
    if problems or kind is None:
        raise UgcError("rebase: " + "; ".join(problems))
    new.entry_kind = kind
    missing = [p for _k, p in new.data_inputs() if not Path(p).exists()]
    if missing:
        raise UgcError("rebase: output file(s) no longer on disk: " + ", ".join(missing))
    new.meta["rebase"] = {"run_id": run_id, "stage": run.stage, "ugc_wgw_version": run.ugc_wgw_version, "at": utc_now(),
                          "previous_entry_kind": rec.entry_kind, "previous_hifi_reads": list(rec.hifi_reads), "previous_fail_reads": list(rec.fail_reads),
                          "previous_aligned_bam": rec.aligned_bam}
    new.added_at = rec.added_at
    db.update_sample(new)
    events.emit("sample.rebased", sample_id=sample_id, run_id=run_id, stage=run.stage, entry_kind=kind, with_vcfs=with_vcfs)
    return new
