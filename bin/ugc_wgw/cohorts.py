"""Frozen cohorts: named, immutable, ordered sample lists with a checksum (docs/DESIGN.md §8.1)."""
from __future__ import annotations

from pathlib import Path

from .db import DB, Cohort
from .log import Events
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


def freeze(db: DB, events: Events, cohort_id: str, ids_path: Path) -> Cohort:
    if not ID_RE.match(cohort_id):
        raise UgcError("cohort_id must match [A-Za-z0-9._-]+")
    ids = read_ids(ids_path)
    sha = sample_list_sha256(ids)
    db.freeze_cohort(cohort_id, ids, sha, utc_now())
    cohort = db.get_cohort(cohort_id)
    assert cohort is not None
    events.emit("cohort.frozen", cohort_id=cohort_id, members=len(ids), sample_list_sha256=sha)
    return cohort
