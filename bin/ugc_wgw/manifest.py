"""run_manifest.json (docs/DESIGN.md §9.2): written by the driver only, once per finished attempt."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from .config import git_commit, read_version
from .db import RunRecord
from .engine import RunResult
from .layout import run_files
from .stages import namespace
from .util import hostname, read_json, sha256_file, write_json

MANIFEST_SCHEMA = 1


@dataclass
class CodeInfo:
    version: str
    git_commit: str
    upstream: dict[str, object] = field(default_factory=dict)
    references: dict[str, object] = field(default_factory=dict)
    containers: list[dict[str, str]] = field(default_factory=list)


def _lock_section(path: Path, key: str) -> dict[str, object]:
    if not path.exists():
        return {}
    doc = read_json(path)
    return dict(doc.get(key, {})) if isinstance(doc, dict) else {}


def _images(path: Path) -> list[str]:
    if not path.exists():
        return []
    out = []
    for line in path.read_text().splitlines():
        s = line.strip()
        if s and not s.startswith("#"):
            out.append(s)
    return out


def load_code_info(code_dir: Path) -> CodeInfo:
    images = set(_images(code_dir / "image_manifest.ugc-wgw.txt"))
    images |= set(_images(code_dir / "vendor" / "hifi-human-wgs-wdl" / "image_manifest.txt"))
    return CodeInfo(
        version=read_version(code_dir),
        git_commit=git_commit(code_dir),
        upstream=_lock_section(code_dir / "upstream.lock", "sources"),
        references=_lock_section(code_dir / "references.lock", "bundles"),
        containers=[{"task": "", "image": img} for img in sorted(images)],
    )


def build(run: RunRecord, result: RunResult, code: CodeInfo, engine: dict[str, str],
          members: list[dict[str, object]], status: str) -> dict[str, object]:
    files = run_files(run.run_path)
    doc: dict[str, object] = {
        "schema": MANIFEST_SCHEMA,
        "ugc_pacbio_wgw": {"version": code.version, "git_commit": code.git_commit},
        "upstream": code.upstream,
        "references": code.references,
        "engine": dict(engine),
        "stage": run.stage,
        "mode": run.mode,
        "subject": {"type": run.subject_type, "id": run.subject_id},
        "cohort_members": members,
        "inputs_file": str(files.inputs),
        "inputs_sha256": run.inputs_sha256 or (sha256_file(files.inputs) if files.inputs.exists() else None),
        "containers": code.containers,
        "containers_source": "image_manifests",
        "run_id": run.run_id,
        "run_dir": run.run_dir,
        "status": status,
        "attempt": run.attempt,
        "started_at": run.started_at,
        "finished_at": run.finished_at,
        "host": run.host or hostname(),
        "slurm_job_ids": list(run.meta.get("slurm_job_ids", [])),  # type: ignore[arg-type]
        "exit_code": result.exit_code,
        "resource_policy": run.meta.get("resource_policy"),  # {path, sha256} of the site policy at launch, or None
    }
    if status != "success":
        doc["error"] = {
            "class": result.error_class,
            "kind": run.error_kind,
            "message": run.error_message or result.error_message,
            "engine_message": result.error_message,
            "exit_code": result.exit_code,
            "exit_status": result.exit_status,
            "task_dir": result.task_dir,
            "node": result.node,
            "evidence": run.meta.get("error_evidence"),
            "not_before": run.not_before,
        }
    if status == "success" and result.outputs:
        path = result.outputs.get(f"{namespace(run.stage)}.ugc_wgw_manifest")
        if isinstance(path, str) and Path(path).exists():
            try:
                doc["ugc_wgw_manifest"] = read_json(Path(path))
            except (OSError, ValueError):
                doc["ugc_wgw_manifest"] = None
    return doc


def write(run: RunRecord, doc: dict[str, object]) -> Path:
    path = run_files(run.run_path).manifest
    write_json(path, doc)
    return path
