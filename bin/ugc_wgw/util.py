"""Small helpers shared by the ugc-wgw driver: errors, time, hashing, atomic JSON, diagnostics."""
from __future__ import annotations

import datetime as dt
import hashlib
import json
import math
import os
import shutil
import socket
import sys
import time
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


def hms(seconds: float) -> str:
    s = int(round(seconds))
    if s < 60:
        return f"{s}s"
    if s < 3600:
        return f"{s // 60}m{s % 60:02d}s"
    return f"{s // 3600}h{(s % 3600) // 60:02d}m"


class Meter:
    """Progress of a long loop on stderr: a line rewritten in place on a terminal, one plain line every
    `interval` seconds otherwise (a log or a batch job), nothing at all for work shorter than `quiet` seconds.
    `update(done, total, label)` after each item; `close(summary)` clears the line and prints the summary."""

    def __init__(self, what: str, *, stream=None, interval: float = 10.0, quiet: float = 0.5,
                 clock=time.monotonic):
        self.what, self.interval, self.quiet, self.clock = what, interval, quiet, clock
        self.stream = stream if stream is not None else sys.stderr
        self.tty = bool(getattr(self.stream, "isatty", lambda: False)())
        self.start = self.clock()
        self.last = float("-inf")
        self.done = self.total = 0
        self.drawn = False

    @property
    def elapsed(self) -> float:
        return self.clock() - self.start

    def update(self, done: int, total: int, label: str = "") -> None:
        self.done, self.total = done, total
        now = self.clock()
        if now - self.start < self.quiet or done >= total:
            return
        if now - self.last < (0.2 if self.tty else self.interval):
            return
        self.last = now
        elapsed = now - self.start
        left = (total - done) * elapsed / done if done else 0.0
        pct = 100.0 * done / total if total else 100.0
        text = f"[ugc-wgw] {self.what} {done}/{total} ({pct:.0f}%), {hms(elapsed)} elapsed, about {hms(left)} left"
        if label:
            text += f": {label}"
        if self.tty:
            width = shutil.get_terminal_size((100, 20)).columns
            self.stream.write("\r\x1b[K" + text[: max(20, width - 1)])
            self.drawn = True
        else:
            self.stream.write(text + "\n")
        self.stream.flush()

    def close(self, summary: str | None = None) -> None:
        if self.drawn:
            self.stream.write("\r\x1b[K")
            self.drawn = False
        if summary:
            self.stream.write(f"[ugc-wgw] {summary}\n")
        self.stream.flush()


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
