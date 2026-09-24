"""Project configuration: `<project>/.ugc-wgw/config.json` and what `ugc-wgw init` derives from an install prefix.

See docs/DESIGN.md §8 (driver), §9.1 (layout), §14 (install prefix layout).
"""
from __future__ import annotations

import json
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

from .log import LOGGER
from .util import UgcError, read_json, write_json

CONFIG_SCHEMA = 1
DEFAULT_REGISTRY = "ghcr.io/uppsalagenomecenter"
DEEPVARIANT_MODES = ("cpu", "gpu", "parabricks")   # small-variant caller of singleton/upstream (guide chapter 12)


@dataclass
class Config:
    project_dir: Path
    code_dir: Path
    miniwdl: Path
    miniwdl_cfg: Path
    results_dir: Path
    ref_map_file: Path
    venv_dir: Path | None = None
    ugc_wgw_container_registry: str = DEFAULT_REGISTRY
    backend: str = "HPC"
    preemptible: bool = True
    max_inflight: int = 4
    poll_interval: float = 30.0
    assembly_use_parents: bool = True
    auto_retry_max: int = 2
    backoff_seconds: float = 300.0
    cancel_orphans: bool = True
    lease_seconds: float = 900.0
    progress_interval: float = 300.0
    stage_inputs: dict[str, dict[str, object]] = field(default_factory=dict)
    summary_thresholds: dict[str, float] = field(default_factory=dict)   # ugc-wgw summary QC thresholds (docs/guide/07-results.md)
    project_url: str = ""                                                # link printed in the analysis summaries
    deepvariant: str = "cpu"        # cpu | gpu (DeepVariant call_variants on 1 GPU) | parabricks (pbrun deepvariant)
    gpu_type: str = ""              # SLURM gres type (`a100`): --gres gpu:<type>:N; empty = gpu:N
    parabricks_gpus: int = 4        # GPUs per Parabricks task (the WDL's own default is 4)

    @property
    def ugc_wgw_dir(self) -> Path:
        return self.project_dir / ".ugc-wgw"

    @property
    def config_path(self) -> Path:
        return self.ugc_wgw_dir / "config.json"

    @property
    def db_path(self) -> Path:
        return self.ugc_wgw_dir / "state.sqlite"

    @property
    def logs_dir(self) -> Path:
        return self.ugc_wgw_dir / "logs"

    @property
    def lock_path(self) -> Path:
        return self.ugc_wgw_dir / "lock"

    def stage_overrides(self, stage: str) -> dict[str, object]:
        """Per-stage input overrides from config.json, keyed `ugc_wgw_<stage>` (bare stage name also accepted)."""
        merged: dict[str, object] = {}
        merged.update(self.stage_inputs.get(stage, {}))
        merged.update(self.stage_inputs.get(f"ugc_wgw_{stage}", {}))
        return merged

    def to_json(self) -> dict[str, object]:
        return {
            "schema": CONFIG_SCHEMA,
            "code_dir": str(self.code_dir),
            "miniwdl": str(self.miniwdl),
            "miniwdl_cfg": str(self.miniwdl_cfg),
            "venv_dir": str(self.venv_dir) if self.venv_dir else None,
            "results_dir": str(self.results_dir),
            "ref_map_file": str(self.ref_map_file),
            "ugc_wgw_container_registry": self.ugc_wgw_container_registry,
            "backend": self.backend,
            "preemptible": self.preemptible,
            "max_inflight": self.max_inflight,
            "poll_interval": self.poll_interval,
            "assembly_use_parents": self.assembly_use_parents,
            "auto_retry_max": self.auto_retry_max,
            "backoff_seconds": self.backoff_seconds,
            "cancel_orphans": self.cancel_orphans,
            "lease_seconds": self.lease_seconds,
            "progress_interval": self.progress_interval,
            "stage_inputs": self.stage_inputs,
            "summary_thresholds": self.summary_thresholds,
            "project_url": self.project_url,
            "deepvariant": self.deepvariant,
            "gpu_type": self.gpu_type,
            "parabricks_gpus": self.parabricks_gpus,
        }

    @classmethod
    def from_json(cls, project_dir: Path, doc: dict[str, object]) -> "Config":
        if doc.get("schema") != CONFIG_SCHEMA:
            raise UgcError(f"unsupported config schema {doc.get('schema')!r} in {project_dir / '.ugc-wgw' / 'config.json'}")
        for old in ("ugc_wgw_ref_map_file", "tertiary_map_file"):
            if doc.get(old):
                LOGGER.warning("config key %s is obsolete since ugc-pacbio-wgw 0.2.0 (one reference map, no tertiary analysis); ignored", old)

        def p(key: str) -> Path:
            return Path(str(doc[key]))

        def opt(key: str) -> Path | None:
            v = doc.get(key)
            return Path(str(v)) if v else None

        return cls(
            project_dir=project_dir,
            code_dir=p("code_dir"),
            miniwdl=p("miniwdl"),
            miniwdl_cfg=p("miniwdl_cfg"),
            results_dir=p("results_dir"),
            ref_map_file=p("ref_map_file"),
            venv_dir=opt("venv_dir"),
            ugc_wgw_container_registry=str(doc.get("ugc_wgw_container_registry", DEFAULT_REGISTRY)),
            backend=str(doc.get("backend", "HPC")),
            preemptible=bool(doc.get("preemptible", True)),
            max_inflight=int(doc.get("max_inflight", 4)),
            poll_interval=float(doc.get("poll_interval", 30.0)),
            assembly_use_parents=bool(doc.get("assembly_use_parents", True)),
            auto_retry_max=int(doc.get("auto_retry_max", 2)),
            backoff_seconds=float(doc.get("backoff_seconds", 300.0)),
            cancel_orphans=bool(doc.get("cancel_orphans", True)),
            lease_seconds=float(doc.get("lease_seconds", 900.0)),
            progress_interval=float(doc.get("progress_interval", 300.0)),
            stage_inputs=dict(doc.get("stage_inputs", {})),  # type: ignore[arg-type]
            summary_thresholds=dict(doc.get("summary_thresholds", {}) or {}),  # type: ignore[arg-type]
            project_url=str(doc.get("project_url") or ""),
            deepvariant=check_deepvariant(str(doc.get("deepvariant") or "cpu")),
            gpu_type=str(doc.get("gpu_type") or "").strip(),
            parabricks_gpus=check_gpus(doc.get("parabricks_gpus", 4)),
        )


def check_deepvariant(mode: str) -> str:
    mode = mode.strip().lower()
    if mode not in DEEPVARIANT_MODES:
        raise UgcError(f"deepvariant must be one of {', '.join(DEEPVARIANT_MODES)}, not {mode!r}")
    return mode


def check_gpus(value: object) -> int:
    try:
        n = int(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        raise UgcError(f"parabricks_gpus must be a positive integer, not {value!r}") from None
    if n < 1:
        raise UgcError(f"parabricks_gpus must be a positive integer, not {value!r}")
    return n


def derive_from_install(install: Path) -> dict[str, Path]:
    """Paths inside an install prefix version dir (`<root>/versions/<v>` or `<root>/current`), DESIGN §14."""
    install = install.resolve()
    paths = {
        "code": install / "code",
        "miniwdl": install / "venv" / "bin" / "miniwdl",
        "cfg": install / "miniwdl.cfg",
        "venv": install / "venv",
    }
    missing = [f"{k}: {v}" for k, v in paths.items() if not v.exists()]
    if missing:
        raise UgcError("install prefix is incomplete; missing " + ", ".join(missing))
    return paths


def read_version(code_dir: Path) -> str:
    path = code_dir / "VERSION"
    try:
        version = path.read_text().strip()
    except OSError as exc:
        raise UgcError(f"cannot read {path}: {exc}") from None
    if not version:
        raise UgcError(f"{path} is empty")
    return version


def git_commit(code_dir: Path) -> str:
    """The commit of the code: from git for a checkout, else from the bundle's manifest.json next to
    the installed code dir (<prefix>/versions/<v>/manifest.json); "unknown" otherwise."""
    if (code_dir / ".git").exists():
        try:
            res = subprocess.run(
                ["git", "-C", str(code_dir), "rev-parse", "HEAD"],
                check=True, capture_output=True, text=True, timeout=30,
            )
            if res.stdout.strip():
                return res.stdout.strip()
        except (OSError, subprocess.SubprocessError):
            pass
    manifest = code_dir.parent / "manifest.json"
    if manifest.exists():
        try:
            doc = json.loads(manifest.read_text())
            commit = doc.get("ugc_pacbio_wgw", {}).get("git_commit")
            if isinstance(commit, str) and commit:
                return commit
        except (OSError, ValueError, AttributeError):
            pass
    return "unknown"


def init_project(
    project_dir: Path,
    *,
    code: Path,
    miniwdl: Path,
    cfg: Path,
    ref_map: Path,
    results: Path | None = None,
    venv: Path | None = None,
    registry: str = DEFAULT_REGISTRY,
    max_inflight: int = 4,
    poll_interval: float = 30.0,
    assembly_use_parents: bool = True,
    auto_retry_max: int = 2,
    backoff_seconds: float = 300.0,
    cancel_orphans: bool = True,
    lease_seconds: float = 900.0,
    progress_interval: float = 300.0,
    deepvariant: str = "cpu",
    gpu_type: str = "",
    parabricks_gpus: int = 4,
) -> Config:
    project_dir = project_dir.resolve()
    config = Config(
        project_dir=project_dir,
        code_dir=code.resolve(),
        miniwdl=miniwdl.resolve(),
        miniwdl_cfg=cfg.resolve(),
        results_dir=(results or project_dir).resolve(),
        ref_map_file=ref_map.resolve(),
        venv_dir=venv.resolve() if venv else None,
        ugc_wgw_container_registry=registry,
        max_inflight=max_inflight,
        poll_interval=poll_interval,
        assembly_use_parents=assembly_use_parents,
        auto_retry_max=auto_retry_max,
        backoff_seconds=backoff_seconds,
        cancel_orphans=cancel_orphans,
        lease_seconds=lease_seconds,
        progress_interval=progress_interval,
        deepvariant=check_deepvariant(deepvariant),
        gpu_type=gpu_type.strip(),
        parabricks_gpus=check_gpus(parabricks_gpus),
    )
    if config.ugc_wgw_dir.exists():
        raise UgcError(f"{config.ugc_wgw_dir} already exists; refusing to re-initialise")
    for label, path in (("code dir", config.code_dir), ("miniwdl", config.miniwdl), ("ref map", config.ref_map_file)):
        if not path.exists():
            raise UgcError(f"{label} not found: {path}")
    read_version(config.code_dir)
    config.logs_dir.mkdir(parents=True)
    config.results_dir.mkdir(parents=True, exist_ok=True)
    save(config)
    return config


def load(project_dir: Path) -> Config:
    project_dir = project_dir.resolve()
    path = project_dir / ".ugc-wgw" / "config.json"
    if not path.exists():
        raise UgcError(f"not a ugc-wgw project: {project_dir} (run `ugc-wgw init`)")
    doc = read_json(path)
    if not isinstance(doc, dict):
        raise UgcError(f"malformed {path}")
    return Config.from_json(project_dir, doc)


def save(config: Config) -> None:
    write_json(config.config_path, config.to_json())
