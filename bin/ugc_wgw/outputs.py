"""Reading a finished run's outputs.json (namespaced keys) and miniwdl error documents."""
from __future__ import annotations

from dataclasses import dataclass

from .db import RunRecord
from .layout import run_files
from .stages import namespace
from .util import UgcError, read_json


class MissingOutputError(UgcError):
    def __init__(self, subject_type: str, subject_id: str, stage: str, run_id: str | None, name: str):
        self.subject_type, self.subject_id, self.stage, self.run_id, self.name = subject_type, subject_id, stage, run_id, name
        where = f"run {run_id}" if run_id else "no successful run"
        super().__init__(f"{subject_type} {subject_id}: stage {stage} has no output {name!r} ({where})")


class Outputs:
    def __init__(self, run: RunRecord, doc: dict[str, object]):
        self.run = run
        self.stage = run.stage
        self.doc = doc

    @classmethod
    def load(cls, run: RunRecord) -> "Outputs":
        path = run_files(run.run_path).outputs
        if not path.exists():
            raise MissingOutputError(run.subject_type, run.subject_id, run.stage, run.run_id, "outputs.json")
        doc = read_json(path)
        if not isinstance(doc, dict):
            raise UgcError(f"malformed {path}")
        return cls(run, doc)

    def key(self, name: str) -> str:
        return f"{namespace(self.stage)}.{name}"

    def opt(self, name: str) -> object | None:
        return self.doc.get(self.key(name))

    def req(self, name: str) -> object:
        v = self.opt(name)
        if v is None:
            raise MissingOutputError(self.run.subject_type, self.run.subject_id, self.stage, self.run.run_id, name)
        return v

    def arr(self, name: str) -> list[object]:
        v = self.req(name)
        if not isinstance(v, list):
            raise UgcError(f"{self.run.run_id}: output {name} is not an array")
        return v

    def opt_arr(self, name: str) -> list[object] | None:
        v = self.opt(name)
        if v is None:
            return None
        if not isinstance(v, list):
            raise UgcError(f"{self.run.run_id}: output {name} is not an array")
        return v


@dataclass
class ErrorInfo:
    error_class: str
    message: str
    task_dir: str | None = None   # directory of the failed task (innermost `dir`, else the top-level `from_dir`)
    node: str | None = None       # workflow node id, e.g. call-deepvariant-03-chr3
    exit_status: int | None = None


def error_info(doc: object) -> ErrorInfo:
    """Follow miniwdl's `cause` chain to the triggering error and keep what a classifier needs."""
    if not isinstance(doc, dict):
        return ErrorInfo("Unknown", "")
    cur: dict[str, object] = doc
    while isinstance(cur.get("cause"), dict):
        cur = cur["cause"]  # type: ignore[assignment]
    cls = str(cur.get("error") or "Unknown")
    msg = cur.get("message")
    exit_status = cur.get("exit_status")
    if not msg and exit_status is not None:
        msg = f"exit_status={exit_status}"
        if cur.get("stderr_file"):
            msg += f" stderr={cur['stderr_file']}"
    task_dir = cur.get("dir") if cur is not doc else None
    if not task_dir:
        task_dir = doc.get("from_dir")
    node = cur.get("node") or cur.get("run")
    return ErrorInfo(cls, str(msg or ""), str(task_dir) if task_dir else None, str(node) if node else None,
                     int(exit_status) if isinstance(exit_status, int) else None)


def innermost_error(doc: object) -> tuple[str, str]:
    """(class, message) of the triggering error; see error_info."""
    info = error_info(doc)
    return (info.error_class, info.message)
