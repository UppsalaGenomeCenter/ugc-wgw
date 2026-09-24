"""Small helpers shared by the ugc-wgw driver: errors, time, hashing, atomic JSON, diagnostics."""
from __future__ import annotations

import datetime as dt
import hashlib
import json
import math
import os
import socket
import sys
from pathlib import Path


class UgcError(Exception):
    """A user-facing error: printed as `[ugc-wgw] error: ...`, exit status 1."""


UTC_FMT = "%Y-%m-%dT%H:%M:%SZ"  # what utc_now() produces; fixed width, so strings compare chronologically


def utc_now() -> str:
    """ISO-8601 UTC timestamp with a Z suffix and no fractional seconds."""
    return dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def parse_utc(ts: str) -> dt.datetime:
    return dt.datetime.strptime(ts, UTC_FMT).replace(tzinfo=dt.timezone.utc)


def utc_plus(seconds: float, start: str | None = None) -> str:
    """`start` (default now) plus `seconds`, rounded up to a whole second, in UTC_FMT."""
    base = parse_utc(start) if start else dt.datetime.now(dt.timezone.utc).replace(microsecond=0)
    return (base + dt.timedelta(seconds=math.ceil(seconds))).strftime(UTC_FMT)


def seconds_until(ts: str, now: str | None = None) -> float:
    base = parse_utc(now) if now else dt.datetime.now(dt.timezone.utc)
    return max(0.0, (parse_utc(ts) - base).total_seconds())


def utc_stamp() -> str:
    """Compact UTC timestamp for identifiers, e.g. 20260910T120000Z."""
    return dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def iso_from_mtime(path: Path) -> str:
    return (
        dt.datetime.fromtimestamp(path.stat().st_mtime, dt.timezone.utc)
        .replace(microsecond=0)
        .isoformat()
        .replace("+00:00", "Z")
    )


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def write_json(path: Path, doc: object) -> None:
    """Atomically write pretty, key-sorted JSON with a trailing newline."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    with open(tmp, "w") as fh:
        json.dump(doc, fh, indent=2, sort_keys=True)
        fh.write("\n")
    os.replace(tmp, path)


def read_json(path: Path) -> object:
    with open(path) as fh:
        return json.load(fh)


def diag(msg: str) -> None:
    print(f"[ugc-wgw] {msg}", file=sys.stderr, flush=True)


def hostname() -> str:
    return socket.gethostname()


def pid_alive(pid: int | None) -> bool:
    """True if a process with this PID exists (or we cannot tell because of permissions)."""
    if not pid or pid <= 0:
        return False
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True
