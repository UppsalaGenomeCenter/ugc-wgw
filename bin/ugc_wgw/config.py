"""Project configuration: `<project>/.ugc-wgw/config.json` and what `ugc-wgw init` derives from an install prefix.

See docs/DESIGN.md §8 (driver), §9.1 (layout), §14 (install prefix layout).
"""
from __future__ import annotations

import json
import os
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

from .log import LOGGER
from .util import UgcError, read_json, sha256_file, write_json

CONFIG_SCHEMA = 1
DEFAULT_REGISTRY = "ghcr.io/uppsalagenomecenter"
DEEPVARIANT_MODES = ("cpu", "gpu", "parabricks")   # small-variant caller of singleton/upstream (guide chapter 12)
# the project keys a site profile may set (`ugc-wgw init --profile`, guide chapter 05); `description` is allowed and ignored
PROFILE_KEYS = ("max_inflight", "poll_interval", "assembly_use_parents", "deepvariant", "gpu_type", "parabricks_gpus",
                "stage_inputs", "prices", "summary_thresholds", "input_thresholds")


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
    input_thresholds: dict[str, float] = field(default_factory=dict)     # input BAM check floors (samples.INPUT_THRESHOLDS)
    project_url: str = ""                                                # link printed in the analysis summaries
    deepvariant: str = "cpu"        # cpu | gpu (DeepVariant call_variants on 1 GPU) | parabricks (pbrun deepvariant)
    gpu_type: str = ""              # SLURM gres type (`a100`): --gres gpu:<type>:N; empty = gpu:N
    parabricks_gpus: int = 4        # GPUs per Parabricks task (the WDL's own default is 4)
    prices: dict[str, object] = field(default_factory=dict)   # unit prices for `usage`/`report` cost columns (guide chapter 07)
    accounting: bool = True         # read `sacct` for every finished run into accounting.json (no-op without sacct)
    accounting_timeout: float = 120.0   # seconds per sacct call
    profile: dict[str, object] = field(default_factory=dict)   # the site profile init applied: name, path, sha256, applied keys

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
        for key in (stage, f"ugc_wgw_{stage}"):
            value = self.stage_inputs.get(key, {})
            if not isinstance(value, dict):
                raise UgcError(f"config.json stage_inputs.{key} must be an object of input: value, not {type(value).__name__}")
            merged.update(value)
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
            "input_thresholds": self.input_thresholds,
            "project_url": self.project_url,
            "deepvariant": self.deepvariant,
            "gpu_type": self.gpu_type,
            "parabricks_gpus": self.parabricks_gpus,
            "prices": self.prices,
            "accounting": self.accounting,
            "accounting_timeout": self.accounting_timeout,
            "profile": self.profile,
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
            input_thresholds=dict(doc.get("input_thresholds", {}) or {}),  # type: ignore[arg-type]
            project_url=str(doc.get("project_url") or ""),
            deepvariant=check_deepvariant(str(doc.get("deepvariant") or "cpu")),
            gpu_type=str(doc.get("gpu_type") or "").strip(),
            parabricks_gpus=check_gpus(doc.get("parabricks_gpus", 4)),
            prices=dict(doc.get("prices") or {}),  # type: ignore[arg-type]
            accounting=bool(doc.get("accounting", True)),
            accounting_timeout=float(doc.get("accounting_timeout", 120.0)),
            profile=dict(doc.get("profile") or {}),  # type: ignore[arg-type]
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


def install_root(install: Path) -> Path | None:
    """`<root>/versions/<v>` or `<root>/current` -> `<root>`; None when the directory is not laid out that way."""
    real = Path(install).resolve()
    return real.parent.parent if real.parent.name == "versions" else None


def find_profile(spec: str, *, install: Path | None, code: Path) -> Path:
    """A profile by file path, else by name under `<root>/profiles/` of the install, else the code's
    `backends/hpc/profiles/` (the examples shipped with the bundle)."""
    p = Path(spec)
    if p.suffix == ".json" or "/" in spec:
        if p.is_file():
            return p.resolve()
        raise UgcError(f"profile file not found: {spec}")
    candidates: list[Path] = []
    root = install_root(install) if install else None
    if root is not None:
        candidates.append(root / "profiles" / f"{spec}.json")
    candidates.append(Path(code) / "backends" / "hpc" / "profiles" / f"{spec}.json")
    for c in candidates:
        if c.is_file():
            return c.resolve()
    raise UgcError(f"profile {spec!r} not found; looked for " + ", ".join(str(c) for c in candidates))


def load_profile(path: Path) -> dict[str, object]:
    """Read and check a site profile: the project keys that depend on the cluster and the sample set."""
    try:
        doc = read_json(path)
    except (OSError, ValueError) as exc:
        raise UgcError(f"profile {path}: {exc}") from None
    if not isinstance(doc, dict):
        raise UgcError(f"profile {path}: not a JSON object")
    out: dict[str, object] = {}
    for key, value in doc.items():
        if key == "description":
            continue
        if key not in PROFILE_KEYS:
            raise UgcError(f"profile {path}: unknown key {key!r} (known: {', '.join(PROFILE_KEYS)}, description)")
        out[key] = value
    if "deepvariant" in out:
        out["deepvariant"] = check_deepvariant(str(out["deepvariant"]))
    if "parabricks_gpus" in out:
        out["parabricks_gpus"] = check_gpus(out["parabricks_gpus"])
    if "gpu_type" in out:
        out["gpu_type"] = str(out["gpu_type"]).strip()
    if "max_inflight" in out and (not isinstance(out["max_inflight"], int) or out["max_inflight"] < 1):
        raise UgcError(f"profile {path}: max_inflight must be a positive integer")
    if "poll_interval" in out and (not isinstance(out["poll_interval"], (int, float)) or out["poll_interval"] <= 0):
        raise UgcError(f"profile {path}: poll_interval must be a positive number")
    if "assembly_use_parents" in out and not isinstance(out["assembly_use_parents"], bool):
        raise UgcError(f"profile {path}: assembly_use_parents must be true or false")
    for key in ("stage_inputs", "prices", "summary_thresholds", "input_thresholds"):
        if key in out and not isinstance(out[key], dict):
            raise UgcError(f"profile {path}: {key} must be an object")
    for stage, values in (out.get("stage_inputs") or {}).items():  # type: ignore[union-attr]
        if not isinstance(values, dict):
            raise UgcError(f"profile {path}: stage_inputs.{stage} must be an object of input: value")
    return out


def profile_meta(path: Path, profile: dict[str, object]) -> dict[str, object]:
    """What config.json records about the profile init applied: name, path, checksum and the keys it set."""
    applied = []
    for key, value in profile.items():
        if key == "stage_inputs" and isinstance(value, dict):
            applied += [f"stage_inputs.{stage}.{name}" for stage, inputs in value.items() for name in inputs]
        else:
            applied.append(key)
    return {"name": Path(path).stem, "path": str(path), "sha256": sha256_file(Path(path)), "applied": sorted(applied)}


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
    stage_inputs: dict[str, object] | None = None,
    prices: dict[str, object] | None = None,
    summary_thresholds: dict[str, object] | None = None,
    input_thresholds: dict[str, object] | None = None,
    profile: dict[str, object] | None = None,
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
        stage_inputs=dict(stage_inputs or {}),  # type: ignore[arg-type]
        prices=dict(prices or {}),
        summary_thresholds=dict(summary_thresholds or {}),  # type: ignore[arg-type]
        input_thresholds=dict(input_thresholds or {}),  # type: ignore[arg-type]
        profile=dict(profile or {}),
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


def device_of(path: Path) -> int | None:
    """The file system (st_dev) of a path, None when it cannot be stat-ed."""
    try:
        return os.stat(path).st_dev
    except OSError:
        return None


def fs_warnings(cfg: Config, device: Callable[[Path], int | None] = device_of) -> list[str]:
    """Warn when the results live on another file system than the install, where the call cache and the
    references sit by default: miniwdl hardlinks outputs, and one it cannot hardlink across (a cached output
    of a run on the other side, an input passed through) becomes a symlink (guide chapter 04)."""
    out: list[str] = []
    results, install = device(cfg.results_dir), device(cfg.code_dir)
    if results is not None and install is not None and results != install:
        out.append(f"results {cfg.results_dir} and the install {cfg.code_dir} are on different file systems: an output "
                   "miniwdl cannot hardlink across (a call-cache hit of a run on the other side, a passed-through "
                   "input) becomes a symlink into that tree, so keep it, or set [call_cache] get = false for this "
                   "project (guide chapter 04)")
    return out


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
