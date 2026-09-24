"""Sample registration from a TSV and per-stage sample listing (docs/DESIGN.md §8.1)."""
from __future__ import annotations

import csv
import re
from pathlib import Path

from .db import DB, SampleRecord
from .log import Events
from .stages import mode_stages, stage_spec
from .util import UgcError, utc_now

COLUMNS = ("sample_id", "sex", "hifi_reads", "fail_reads", "father_id", "mother_id")
REQUIRED_COLUMNS = ("sample_id", "hifi_reads")
SEX_VALUES = ("MALE", "FEMALE")
ID_RE = re.compile(r"^[A-Za-z0-9._-]+$")


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


def validate(records: list[SampleRecord], db: DB, check_paths: bool = True) -> list[str]:
    """Return warnings; raise UgcError with every hard problem found."""
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
    if problems:
        raise UgcError("invalid samples TSV:\n  " + "\n  ".join(problems))
    return warnings


def add_samples(db: DB, events: Events, path: Path, check_paths: bool = True, replace: bool = False) -> int:
    records = parse_tsv(path)
    for w in validate(records, db, check_paths):
        events.emit("sample.warning", level="warning", message=w)
    if not replace:
        dup = [r.sample_id for r in records if db.sample_exists(r.sample_id)]
        if dup:
            raise UgcError("already registered (use --replace): " + ", ".join(dup))
    for rec in records:
        db.upsert_sample(rec, replace=replace)
        events.emit("sample.added", sample_id=rec.sample_id, hifi_reads=len(rec.hifi_reads),
                    fail_reads=len(rec.fail_reads), replaced=replace)
    return len(records)


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
