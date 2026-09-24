"""The project lease in `<project>/.ugc-wgw/lock`: who drives the project, and whether they are still alive.

flock guards the file on one host; the lease line (`pid= host= since= last_seen=`) is what another host can
read. A driver refreshes `last_seen` every poll; a lease older than `lease_seconds` is expired.
"""
from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

from .util import hostname, parse_utc, utc_now


@dataclass
class Lease:
    pid: int | None
    host: str | None
    since: str | None
    last_seen: str | None

    def expired(self, lease_seconds: float, now: str | None = None) -> bool:
        if not self.last_seen:
            return True
        try:
            age = (parse_utc(now or utc_now()) - parse_utc(self.last_seen)).total_seconds()
        except ValueError:
            return True
        return age > lease_seconds

    def is_mine(self) -> bool:
        return self.pid == os.getpid() and self.host == hostname()

    def describe(self) -> str:
        return f"pid={self.pid} host={self.host} since={self.since} last_seen={self.last_seen or 'never'}"


def parse(text: str) -> Lease | None:
    line = next((ln for ln in text.splitlines() if ln.strip()), "")
    if not line:
        return None
    fields = dict(tok.split("=", 1) for tok in line.split() if "=" in tok)
    pid = fields.get("pid")
    return Lease(int(pid) if pid and pid.isdigit() else None, fields.get("host"), fields.get("since"),
                 fields.get("last_seen"))


def read(path: Path) -> Lease | None:
    try:
        return parse(Path(path).read_text())
    except OSError:
        return None


def render(pid: int, host: str, since: str, last_seen: str) -> str:
    return f"pid={pid} host={host} since={since} last_seen={last_seen}\n"
