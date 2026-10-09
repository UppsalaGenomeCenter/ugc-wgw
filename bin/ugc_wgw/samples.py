"""Sample registration from a TSV, input BAM checks and per-stage sample listing (docs/DESIGN.md §8.1)."""
from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Callable, Optional

from . import bam
from .db import DB, SampleRecord
from .log import LOGGER, Events
from .stages import mode_stages, stage_spec
from .util import UgcError, utc_now

COLUMNS = ("sample_id", "sex", "hifi_reads", "fail_reads", "father_id", "mother_id")
REQUIRED_COLUMNS = ("sample_id", "hifi_reads")
SEX_VALUES = ("MALE", "FEMALE")
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
        for n, row in enumerate(reader, start=2):
            sid = (row.get("sample_id") or "").strip()
            if not sid:
                continue
            sex_raw = (row.get("sex") or "").strip().upper()
            meta = {k: v for k, v in row.items() if k and k not in COLUMNS and v not in (None, "")}
            records.append(
                SampleRecord(
                    sample_id=sid,
                    sex=sex_raw or None,
                    father_id=(row.get("father_id") or "").strip() or None,
                    mother_id=(row.get("mother_id") or "").strip() or None,
                    added_at=added_at,
                    hifi_reads=_split_paths(row.get("hifi_reads")),
                    fail_reads=_split_paths(row.get("fail_reads")),
                    meta={"tsv_line": n, **meta},
                )
            )
    if not records:
        raise UgcError(f"{path}: no samples")
    return records


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
            problems.append(f"{rec.sample_id}: sex must be MALE, FEMALE or blank, not {rec.sex!r}")
        if not rec.hifi_reads:
            problems.append(f"{rec.sample_id}: at least one hifi_reads path is required")
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
    n_files = sum(len(r.hifi_reads) + len(r.fail_reads) for r in records)
    done = 0
    for rec in records:
        total_reads = total_bases = 0
        exact = True
        dropped: list[dict[str, object]] = []
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
        events.emit("sample.added", sample_id=rec.sample_id, hifi_reads=len(rec.hifi_reads),
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
    n_files = sum(len(r.hifi_reads) + len(r.fail_reads) for r in recs.values() if r is not None)
    done = 0
    for sid in sample_ids:
        rec = recs[sid]
        assert rec is not None
        total_bases = 0
        for kind in ("hifi_reads", "fail_reads"):
            for p in getattr(rec, kind):
                done += 1
                if progress and not stored:
                    progress(done, n_files, f"{sid} {Path(p).name}")
                if stored:
                    info_d = rec.input_info.get(p)
                    if not info_d:
                        rows.append({"sample_id": sid, "kind": kind, "file": Path(p).name, "status": "not inspected"})
                        continue
                    info = bam.BamInfo.from_dict(info_d)
                else:
                    info = bam.inspect(Path(p))
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
        row: dict[str, object] = {"sample_id": rec.sample_id, "sex": rec.sex or ""}
        for s in stages:
            run = latest.get(("sample", rec.sample_id, s))
            row[s] = run.status if run else "-"
        if statuses and not any(row[s] in statuses for s in stages):
            continue
        rows.append(row)
    return rows
