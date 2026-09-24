"""Launching `miniwdl run` and reading what it left behind (verified against miniwdl 1.15.0)."""
from __future__ import annotations

import os
import subprocess
from dataclasses import dataclass
from pathlib import Path

from .config import Config
from .layout import run_files
from .outputs import error_info
from .stages import wdl_path
from .util import read_json


def miniwdl_cmd(cfg: Config, stage: str, attempt_path: Path) -> list[str]:
    files = run_files(attempt_path)
    return [
        str(cfg.miniwdl), "run", str(wdl_path(cfg.code_dir, stage)),
        "-i", str(files.inputs),
        "--dir", f"{attempt_path}/.",      # trailing /. = run exactly in this directory
        "-o", str(files.run_json),          # atomic success-or-error JSON (implies --error-json)
        "--cfg", str(cfg.miniwdl_cfg),
        "--no-color", "--log-json",
    ]


def launch(cmd: list[str], attempt_path: Path) -> subprocess.Popen:
    files = run_files(attempt_path)
    env = dict(os.environ)
    env["NO_COLOR"] = "1"
    out = open(files.stdout, "ab")
    err = open(files.stderr, "ab")
    try:
        return subprocess.Popen(cmd, cwd=str(attempt_path), stdout=out, stderr=err, env=env, start_new_session=True)
    finally:
        out.close()
        err.close()


def _run(cmd: list[str]) -> str:
    try:
        res = subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=60)
    except (OSError, subprocess.SubprocessError):
        return "unknown"
    text = res.stdout.strip() or res.stderr.strip()
    return text.splitlines()[0] if text else "unknown"


def probe(cfg: Config) -> dict[str, str]:
    info = {"miniwdl": _run([str(cfg.miniwdl), "--version"]), "miniwdl_slurm": "unknown", "apptainer": "unknown",
            "ugc_wgw_miniwdl": "unknown"}
    python = (cfg.venv_dir / "bin" / "python") if cfg.venv_dir else cfg.miniwdl.parent / "python"
    if python.exists():
        info["miniwdl_slurm"] = _run([str(python), "-c",
                                      "import importlib.metadata as m; print(m.version('miniwdl-slurm'))"])
        info["ugc_wgw_miniwdl"] = _run([str(python), "-c",
                                    "import importlib.metadata as m; print(m.version('ugc-wgw-miniwdl'))"])
    info["apptainer"] = _run(["apptainer", "--version"])
    return info


@dataclass
class RunResult:
    status: str                      # "success" | "failed"
    exit_code: int | None
    error_class: str | None = None
    error_message: str = ""
    outputs: dict[str, object] | None = None
    task_dir: str | None = None    # failed task's directory, from miniwdl's error document
    node: str | None = None
    exit_status: int | None = None


def _failed(doc: object, exit_code: int | None) -> RunResult:
    info = error_info(doc)
    return RunResult("failed", exit_code, info.error_class, info.message, task_dir=info.task_dir, node=info.node,
                     exit_status=info.exit_status)


def read_result(attempt_path: Path, exit_code: int | None) -> RunResult | None:
    """Decide from the run directory; None if undecidable (no run.json, outputs.json or error.json)."""
    files = run_files(attempt_path)
    run_doc: object = read_json(files.run_json) if files.run_json.exists() else None
    err_doc: object = read_json(files.error) if files.error.exists() else None
    failed_exit = exit_code not in (None, 0)
    if isinstance(run_doc, dict) and "outputs" in run_doc and files.outputs.exists() and not failed_exit:
        outputs = read_json(files.outputs)
        return RunResult("success", exit_code, outputs=outputs if isinstance(outputs, dict) else {})
    if isinstance(run_doc, dict) and "error" in run_doc:
        return _failed(run_doc, exit_code)
    if isinstance(err_doc, dict):
        return _failed(err_doc, exit_code)
    if files.outputs.exists() and not failed_exit:
        outputs = read_json(files.outputs)
        return RunResult("success", exit_code, outputs=outputs if isinstance(outputs, dict) else {})
    if failed_exit:
        return RunResult("failed", exit_code, "NoResult", f"miniwdl exited {exit_code} without run.json")
    return None
