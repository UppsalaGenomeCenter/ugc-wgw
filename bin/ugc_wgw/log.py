"""Human log (`.ugc-wgw/logs/ugc-wgw.log`) and the JSONL event stream (`.ugc-wgw/logs/events.jsonl` + `events` table).

The file log is always plain. The stderr copy is coloured by meaning when `--color` allows it (guide chapter 11):
the event token by kind (green success, red failure, cyan submission, dimmed waiting, magenta driver actions),
`subject=`/`stage=` values in bold, failure details in red, run ids and paths dimmed; never in the file.
"""
from __future__ import annotations

import json
import logging
import os
import re
import sys
import time
from pathlib import Path
from typing import IO

from .db import DB
from .util import UgcError, utc_now

LOGGER = logging.getLogger("ugc-wgw")
LOG_FORMAT = "%(asctime)sZ %(levelname)s %(message)s"
DATE_FORMAT = "%Y-%m-%dT%H:%M:%S"
COLOR_MODES = ("auto", "always", "never")

RESET, BOLD, DIM = "\x1b[0m", "\x1b[1m", "\x1b[2m"
RED, GREEN, YELLOW, MAGENTA, CYAN = "\x1b[31m", "\x1b[32m", "\x1b[33m", "\x1b[35m", "\x1b[36m"
LEVEL_COLOR = {"DEBUG": DIM, "INFO": DIM, "WARNING": BOLD + YELLOW, "ERROR": BOLD + RED, "CRITICAL": BOLD + RED}
EVENT_COLOR = {
    "run.success": GREEN,
    "run.failed": RED, "run.cancelled": RED, "run.driver_lost": RED,
    "run.created": CYAN, "run.submitted": CYAN,
    "run.blocked": DIM, "run.backoff": DIM,
    "run.auto_retry": MAGENTA, "run.scancel": MAGENTA, "run.reconciled": MAGENTA,
    "submit.start": BOLD, "submit.stop": BOLD, "submit.progress": DIM,
    "progress": BOLD + CYAN,
}
FAILURE_EVENTS = {"run.failed", "run.cancelled", "run.driver_lost"}
FIELD_RE = re.compile(r"\b(subject|stage|mode|cohort)=(\S+)")
BAD_RE = re.compile(r"\b(error_class|kind|exit_code|exit_status|reason|message)=(\S+)")
DIM_RE = re.compile(r"\b(run_id|pid|cmd|run_dir|task_dir|stderr|engine_message)=(\S+)")


def color_enabled(mode: str, stream: IO[str] | None = None) -> bool:
    """`always`, `never`, or `auto`: only on a terminal, with NO_COLOR unset and TERM not `dumb`."""
    mode = (mode or "auto").lower()
    if mode not in COLOR_MODES:
        raise UgcError(f"--color/UGC_WGW_COLOR must be one of {', '.join(COLOR_MODES)}, not {mode!r}")
    if mode == "always":
        return True
    if mode == "never" or os.environ.get("NO_COLOR") or os.environ.get("TERM", "") == "dumb":
        return False
    stream = stream if stream is not None else sys.stderr
    return bool(getattr(stream, "isatty", lambda: False)())


class ColorFormatter(logging.Formatter):
    """The plain format with ANSI colours; stripping the escapes gives exactly the file log's line."""

    def __init__(self) -> None:
        super().__init__(LOG_FORMAT, datefmt=DATE_FORMAT)
        self.converter = time.gmtime

    def format(self, record: logging.LogRecord) -> str:
        ts = self.formatTime(record, self.datefmt) + "Z"
        msg = record.getMessage()
        first, newline, tail = msg.partition("\n")
        token, space, rest = first.partition(" ")
        colour = EVENT_COLOR.get(token, CYAN if token.startswith("run.") else "")
        if record.levelno >= logging.WARNING:
            colour = colour or LEVEL_COLOR[record.levelname]
        if token in FAILURE_EVENTS or record.levelno >= logging.WARNING:
            rest = BAD_RE.sub(lambda m: f"{m.group(1)}={RED}{m.group(2)}{RESET}", rest)
        rest = FIELD_RE.sub(lambda m: f"{m.group(1)}={BOLD}{m.group(2)}{RESET}", rest)
        rest = DIM_RE.sub(lambda m: f"{DIM}{m.group(1)}={m.group(2)}{RESET}", rest)
        if token == "progress":
            body = f"{colour}{first}{RESET}"
        else:
            body = f"{colour}{token}{RESET}{space}{rest}" if colour else f"{token}{space}{rest}"
        level = f"{LEVEL_COLOR.get(record.levelname, '')}{record.levelname}{RESET}"
        return f"{DIM}{ts}{RESET} {level} {body}{newline}{tail}"


def setup_logging(logs_dir: Path, verbose: bool = False, color: str = "auto") -> None:
    logs_dir.mkdir(parents=True, exist_ok=True)
    fmt = logging.Formatter(LOG_FORMAT, datefmt=DATE_FORMAT)
    fmt.converter = time.gmtime
    LOGGER.setLevel(logging.INFO)
    LOGGER.handlers.clear()
    fh = logging.FileHandler(logs_dir / "ugc-wgw.log")
    fh.setFormatter(fmt)
    fh.setLevel(logging.INFO)
    LOGGER.addHandler(fh)
    sh = logging.StreamHandler()
    sh.setFormatter(ColorFormatter() if color_enabled(color, sh.stream) else fmt)
    sh.setLevel(logging.INFO if verbose else logging.WARNING)
    LOGGER.addHandler(sh)
    LOGGER.propagate = False


class Events:
    """One JSON object per line: {"ts", "run_id", "level", "event", "detail"}; mirrored into the DB and the log."""

    def __init__(self, db: DB, path: Path):
        self.db = db
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)

    def emit(self, event: str, run_id: str | None = None, level: str = "info", **detail: object) -> None:
        ts = utc_now()
        rec = {"ts": ts, "run_id": run_id, "level": level, "event": event, "detail": detail}
        with open(self.path, "a") as fh:
            fh.write(json.dumps(rec, sort_keys=True, default=str) + "\n")
        self.db.add_event(ts, run_id, level, event, detail)
        msg = f"{event} run_id={run_id} " + " ".join(f"{k}={v}" for k, v in detail.items())
        getattr(LOGGER, level if level in ("debug", "info", "warning", "error") else "info")(msg)
