"""Frozen cohorts: named, immutable, ordered sample lists with a checksum (docs/DESIGN.md §8.1)."""
from __future__ import annotations

from pathlib import Path

from . import summary_parsers as sp
from .db import DB, Cohort
from .log import LOGGER, Events
from .outputs import Outputs
from .samples import ID_RE
from .util import UgcError, sha256_text, utc_now


def read_ids(path: Path) -> list[str]:
    ids: list[str] = []
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            ids.append(s.split()[0])
    if not ids:
        raise UgcError(f"{path}: no sample IDs")
    dup = sorted({i for i in ids if ids.count(i) > 1})
    if dup:
        raise UgcError(f"{path}: duplicated IDs: {', '.join(dup)}")
    return ids


def sample_list_sha256(ids: list[str]) -> str:
    return sha256_text("\n".join(ids) + "\n")


def duplicate_warnings(db: DB, ids: list[str], threshold: float) -> list[str]:
    """Pairs of the new list that a finished cohort_qc run of this project reported at or above `threshold`
    relatedness (a duplicate submission, or monozygotic twins): one line per pair, from the latest such run."""
    want = set(ids)
    seen: dict[tuple[str, str], tuple[str, float, str]] = {}
    for run in db.all_runs():
        if run.stage != "cohort_qc" or run.status != "success":
            continue
        try:
            path = Outputs.load(run).opt("qc_related_pairs")
        except UgcError:
            continue
        if not isinstance(path, str) or not Path(path).exists():
            continue
        try:
            _cols, rows = sp.table(sp.read_lines(Path(path)))
        except (OSError, ValueError):
            continue
        for r in rows:
            a, b = r.get("sample_a", ""), r.get("sample_b", "")
            try:
                rel = float(r.get("relatedness", ""))
            except ValueError:
                continue
            if a in want and b in want and rel >= threshold:
                key = tuple(sorted((a, b)))
                prev = seen.get(key)
                if prev is None or (run.finished_at or "") > prev[2]:
                    seen[key] = (run.subject_id, rel, run.finished_at or "")   # type: ignore[assignment]
    return [f"{a} and {b} look like the same individual (relatedness {rel:.3f} in cohort_qc of {cid})"
            for (a, b), (cid, rel, _t) in sorted(seen.items())]


def freeze(db: DB, events: Events, cohort_id: str, ids_path: Path, *, duplicate_threshold: float | None = None) -> Cohort:
    """Freeze the list; with a threshold, warn (never refuse) about members a cohort_qc run showed to be duplicates."""
    if not ID_RE.match(cohort_id):
        raise UgcError("cohort_id must match [A-Za-z0-9._-]+")
    ids = read_ids(ids_path)
    warnings = duplicate_warnings(db, ids, duplicate_threshold) if duplicate_threshold is not None else []
    sha = sample_list_sha256(ids)
    db.freeze_cohort(cohort_id, ids, sha, utc_now())
    cohort = db.get_cohort(cohort_id)
    assert cohort is not None
    events.emit("cohort.frozen", cohort_id=cohort_id, members=len(ids), sample_list_sha256=sha)
    for w in warnings:
        LOGGER.warning("cohort %s: %s", cohort_id, w)
        events.emit("cohort.warning", cohort_id=cohort_id, message=w)
    return cohort
