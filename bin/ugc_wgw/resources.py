"""Task resources as the project's install sees them: declared in the WDL -> site caps -> policy rows.

A read-only view for `ugc-wgw resources` and the submit preflight (docs/guide/12-resources.md). The policy
itself is applied at task launch by the ugc_wgw_resources miniwdl plugin; this module reads the same file
with the same parser (plugins/ugc_wgw_miniwdl/ugc_wgw_miniwdl/policy.py, imported from the code dir, stdlib
only) and the same rendered miniwdl.cfg, so what it prints is what sbatch will be asked for.
"""
from __future__ import annotations

import configparser
import json
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .config import Config
from .util import UgcError, sha256_file

DECLARED_REL = Path("backends") / "hpc" / "resources.declared.tsv"
PARTITION_ARG = re.compile(r"(^|\s)(-p|--partition)(\s|=|$)")
# miniwdl's own byte-size table (WDL/_util.py) for hand-edited memory_max values; the installer renders bytes
MINIWDL_UNITS = {"B": 1, "K": 1000, "KB": 1000, "Ki": 1024, "KiB": 1024, "M": 10**6, "MB": 10**6, "Mi": 1024**2,
                 "MiB": 1024**2, "G": 10**9, "GB": 10**9, "Gi": 1024**3, "GiB": 1024**3, "T": 10**12, "TB": 10**12,
                 "Ti": 1024**4, "TiB": 1024**4}


def policy_module(code_dir: Path) -> Any:
    """The plugin's policy module, loaded from the installed code (the driver runs outside the venv)."""
    pkg = code_dir / "plugins" / "ugc_wgw_miniwdl"
    if str(pkg) not in sys.path:
        sys.path.insert(0, str(pkg))
    try:
        from ugc_wgw_miniwdl import policy  # type: ignore[import-not-found]
    except ImportError as exc:
        raise UgcError(f"resource policy parser not found under {pkg}: {exc}") from None
    return policy


# ---- declared (the generated inventory) -----------------------------------------------

@dataclass
class Declared:
    task: str
    cpu: int | None
    memory: int | None          # bytes
    stages: list[str]
    source: str
    notes: str


def read_declared(code_dir: Path) -> list[Declared]:
    path = code_dir / DECLARED_REL
    if not path.exists():
        raise UgcError(f"task inventory missing: {path} (regenerate with scripts/wdl-resources.py --write)")
    policy = policy_module(code_dir)
    out: list[Declared] = []
    header: list[str] | None = None
    for line in path.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        cells = line.split("\t")
        if header is None:
            header = cells
            continue
        row = dict(zip(header, cells + [""] * (len(header) - len(cells))))
        cpu = int(row["cpu"]) if row.get("cpu", "-").isdigit() else None
        mem = row.get("memory", "-")
        memory = policy.parse_memory(mem) if mem not in ("", "-") else None
        out.append(Declared(row["task"], cpu, memory, [s for s in row.get("stages", "").split(",") if s],
                            row.get("source", ""), row.get("notes", "")))
    return out


# ---- the rendered miniwdl.cfg ----------------------------------------------------------

@dataclass
class EngineCfg:
    path: Path
    cpu_max: int = 0
    memory_max: int = 0             # bytes; 0 or -1 = no cap
    policy_path: str = ""
    extra_args: str = ""
    defaults: dict[str, Any] = field(default_factory=dict)
    backend: str = ""               # [scheduler] container_backend
    run_options: list[str] = field(default_factory=list)   # [singularity] run_options (--nv for the GPU tasks)


def _strip(value: str) -> str:
    v = value.strip()
    if len(v) >= 2 and v[0] == v[-1] and v[0] in "\"'":
        v = v[1:-1]
    return v


def _byte_size(text: str) -> int:
    t = text.strip()
    if t in ("0", "-1", ""):
        return int(t or 0)
    m = re.match(r"^(\d+(?:\.\d+)?)\s*([A-Za-z]*)$", t)
    if not m or (m.group(2) and m.group(2) not in MINIWDL_UNITS):
        raise UgcError(f"memory_max in miniwdl.cfg is not a byte size miniwdl understands: {text!r}")
    return int(float(m.group(1)) * MINIWDL_UNITS.get(m.group(2), 1))


def read_engine_cfg(path: Path) -> EngineCfg:
    eng = EngineCfg(path)
    if not path.exists():
        return eng
    cp = configparser.ConfigParser(interpolation=None, strict=False)
    try:
        cp.read(path)
    except configparser.Error as exc:
        raise UgcError(f"cannot parse {path}: {exc}") from None
    if cp.has_option("task_runtime", "cpu_max"):
        try:
            eng.cpu_max = int(_strip(cp.get("task_runtime", "cpu_max")))
        except ValueError:
            raise UgcError(f"cpu_max in {path} is not an integer") from None
    if cp.has_option("task_runtime", "memory_max"):
        eng.memory_max = _byte_size(_strip(cp.get("task_runtime", "memory_max")))
    if cp.has_option("task_runtime", "defaults"):
        try:
            eng.defaults = json.loads(cp.get("task_runtime", "defaults"))
        except ValueError:
            eng.defaults = {}
    if cp.has_option("ugc_wgw", "resources"):
        eng.policy_path = _strip(cp.get("ugc_wgw", "resources"))
    if cp.has_option("slurm", "extra_args"):
        eng.extra_args = _strip(cp.get("slurm", "extra_args"))
    if cp.has_option("scheduler", "container_backend"):
        eng.backend = _strip(cp.get("scheduler", "container_backend"))
    if cp.has_option("singularity", "run_options"):
        try:
            opts = json.loads(cp.get("singularity", "run_options"))
            eng.run_options = [str(o) for o in opts] if isinstance(opts, list) else []
        except ValueError:
            eng.run_options = []
    return eng


# ---- the GPU tasks -----------------------------------------------------------------------

# The inventory marks the two tasks that ask for GPUs (runtime gpuCount). deepvariant.wdl's call_variants_gpu
# declares one GPU; parabricks_deepvariant.wdl's task defaults to four, set per project by config.json
# `parabricks_gpus` (inputs.PARABRICKS_GPUS_INPUT). Which of them runs is config.json `deepvariant`.
GPU_TASKS = {"deepvariant_call_variants_gpu": "gpu", "run_parabricks_deepvariant": "parabricks"}


@dataclass
class GpuChoice:
    deepvariant: str = "cpu"
    gpu_type: str = ""
    parabricks_gpus: int = 4

    def count(self, task: str) -> int:
        return self.parabricks_gpus if task == "run_parabricks_deepvariant" else 1

    def gres(self, task: str) -> str:
        """The `gpus` cell: what miniwdl-slurm puts after `--gres gpu:`."""
        n = self.count(task)
        return f"{self.gpu_type}:{n}" if self.gpu_type else str(n)

    def selected_task(self) -> str | None:
        return next((task for task, mode in GPU_TASKS.items() if mode == self.deepvariant), None)


# ---- effective values -------------------------------------------------------------------

@dataclass
class Effective:
    task: str
    stages: list[str]
    declared_cpu: int | None
    declared_memory: int | None
    cpu: int | None
    memory: int | None
    time: int | None
    partition: str | None
    constraint: str | None
    gpus: str | None             # `[type:]N` for the GPU tasks, None otherwise
    sources: dict[str, str]      # cpu/memory/time/partition/constraint -> declared | clamp | default | policy <file:line>
    notes: str
    source: str

    @property
    def changed(self) -> bool:
        return any(v.startswith(("clamp", "policy")) for v in self.sources.values())


def effective(dec: Declared, rules: list[Any], eng: EngineCfg, policy: Any, gpu: GpuChoice | None = None) -> Effective:
    gpu = gpu or GpuChoice()
    sources: dict[str, str] = {}
    cpu, memory = dec.cpu, dec.memory
    sources["cpu"] = "declared" if cpu is not None else "dynamic"
    sources["memory"] = "declared" if memory is not None else "dynamic"
    if cpu is not None and eng.cpu_max > 0 and cpu > eng.cpu_max:
        cpu, sources["cpu"] = eng.cpu_max, "clamp"
    if memory is not None and eng.memory_max > 0 and memory > eng.memory_max:
        memory, sources["memory"] = eng.memory_max, "clamp"
    time = eng.defaults.get("time_minutes") if isinstance(eng.defaults.get("time_minutes"), int) else None
    partition = eng.defaults.get("slurm_partition") if isinstance(eng.defaults.get("slurm_partition"), str) else None
    gpus = gpu.gres(dec.task) if dec.task in GPU_TASKS else None
    if gpus is not None and isinstance(eng.defaults.get("slurm_partition_gpu"), str):
        partition = eng.defaults["slurm_partition_gpu"]   # miniwdl-slurm prefers it for tasks with gpuCount
    constraint = None
    sources["time"] = "default" if time is not None else "-"
    sources["partition"] = "default" if partition is not None else "-"
    sources["constraint"] = "-"
    for col, (value, rule) in policy.resolve(rules, dec.task).items():
        src = f"policy {rule.source}"
        if col == "cpu":
            cpu = int(value)
        elif col == "memory":
            memory = int(value)
        elif col == "time":
            time = int(value)
        elif col == "partition":
            partition = str(value)
        elif col == "constraint":
            constraint = str(value)
        sources[col] = src
    return Effective(dec.task, dec.stages, dec.cpu, dec.memory, cpu, memory, time, partition, constraint, gpus, sources,
                     dec.notes, dec.source)


# ---- validation and the table -------------------------------------------------------------

@dataclass
class Resources:
    engine: EngineCfg
    rules: list[Any]
    declared: list[Declared]
    warnings: list[str]
    policy: Any
    gpu: GpuChoice = field(default_factory=GpuChoice)

    def summary(self) -> str:
        pol = f"{self.engine.policy_path} ({len(self.rules)} row{'s' if len(self.rules) != 1 else ''})" if self.engine.policy_path else "none"
        mem = self.policy.format_memory(self.engine.memory_max) if self.engine.memory_max > 0 else "-"
        cpu = self.engine.cpu_max if self.engine.cpu_max > 0 else "-"
        gpu = f" deepvariant={self.gpu.deepvariant}"
        if self.gpu.deepvariant != "cpu":
            gpu += f" gpu_type={self.gpu.gpu_type or '-'}"
        return f"resources: cpu_max={cpu} memory_max={mem} policy={pol}{gpu}"

    def rows(self, *, stage: str | None = None, changed: bool = False) -> list[Effective]:
        out = []
        for dec in self.declared:
            if stage and stage not in dec.stages:
                continue
            eff = effective(dec, self.rules, self.engine, self.policy, self.gpu)
            if changed and not eff.changed:
                continue
            out.append(eff)
        return out


def load(cfg: Config, *, engine_info: dict[str, str] | None = None) -> Resources:
    """Read the inventory, the rendered cfg and the policy; raise UgcError on a malformed policy."""
    policy = policy_module(cfg.code_dir)
    declared = read_declared(cfg.code_dir)
    eng = read_engine_cfg(cfg.miniwdl_cfg)
    rules: list[Any] = []
    warnings: list[str] = []
    if eng.policy_path:
        if Path(eng.policy_path).exists():
            try:
                rules = policy.parse_policy(eng.policy_path)
            except policy.PolicyError as exc:
                raise UgcError(f"resource policy: {exc}") from None
        else:
            warnings.append(f"resource policy file missing: {eng.policy_path} (no rows apply; the installer creates it)")
    elif cfg.miniwdl_cfg.exists():
        warnings.append(f"{cfg.miniwdl_cfg} has no [ugc_wgw] resources entry: no per-task policy applies (older install?)")
    names = {d.task for d in declared}
    for rule in rules:
        if not any(rule.matches(n) for n in names):
            warnings.append(f"policy row {rule.source} ({rule.pattern}) matches no task in the inventory (typo?)")
    if any(r.partition is not None for r in rules) and PARTITION_ARG.search(eng.extra_args):
        warnings.append(f"policy rows set partitions but [slurm] extra_args carries a --partition ({eng.extra_args!r}); "
                        "sbatch keeps the last one, so those rows have no effect: move the partition to SLURM_PARTITION")
    if rules and engine_info is not None and engine_info.get("ugc_wgw_miniwdl", "unknown") == "unknown":
        warnings.append("the engine venv has no ugc_wgw_resources task plugin: the policy rows will not apply "
                        "(install a bundle that ships plugins/ugc_wgw_miniwdl)")
    gpu = GpuChoice(cfg.deepvariant, cfg.gpu_type, cfg.parabricks_gpus)
    task = gpu.selected_task()
    if task is not None and cfg.miniwdl_cfg.exists():
        if "--nv" not in eng.run_options:
            warnings.append(f"deepvariant={cfg.deepvariant} but [singularity] run_options has no --nv, so {task} "
                            "would not see a GPU: set SINGULARITY_NV=1 in site.cfg (or add \"--nv\" to run_options)")
        if eng.backend == "slurm_singularity" and not isinstance(eng.defaults.get("slurm_partition_gpu"), str) \
                and policy.resolve(rules, task).get("partition") is None:
            warnings.append(f"deepvariant={cfg.deepvariant} but no GPU partition: set SLURM_PARTITION_GPU in site.cfg "
                            f"or a policy row with a partition for {task}")
    return Resources(eng, rules, declared, warnings, policy, gpu)


def policy_ref(cfg: Config) -> dict[str, str] | None:
    """{path, sha256} of the policy file the install points at, for the run manifest; None when there is none."""
    try:
        eng = read_engine_cfg(cfg.miniwdl_cfg)
    except UgcError:
        return None
    if not eng.policy_path or not Path(eng.policy_path).exists():
        return None
    return {"path": eng.policy_path, "sha256": sha256_file(Path(eng.policy_path))}


def as_rows(res: Resources, effs: list[Effective]) -> list[dict[str, object]]:
    """Rows for _print_rows / --json."""
    fm, ft = res.policy.format_memory, res.policy.format_time
    rows: list[dict[str, object]] = []
    for e in effs:
        def cell(value: object, src: str, fmt) -> str:  # noqa: E306
            shown = fmt(value) if value is not None else "-"
            return shown if src in ("declared", "default", "-", "dynamic") else f"{shown} ({src.split(' ')[0]})"
        rows.append({
            "task": e.task,
            "stages": ",".join(s[len("ugc_wgw_"):] if s.startswith("ugc_wgw_") else s for s in e.stages),
            "declared": f"{e.declared_cpu if e.declared_cpu is not None else '-'} / "
                        f"{fm(e.declared_memory) if e.declared_memory is not None else '-'}",
            "cpu": cell(e.cpu, e.sources["cpu"], str),
            "memory": cell(e.memory, e.sources["memory"], fm),
            "time": cell(e.time, e.sources["time"], ft),
            "partition": cell(e.partition, e.sources["partition"], str),
            "constraint": cell(e.constraint, e.sources["constraint"], str),
            "gpus": e.gpus if e.gpus is not None else "-",
            "rules": ", ".join(sorted({v.split(" ", 1)[1] for v in e.sources.values() if v.startswith("policy ")})),
            "notes": e.notes,
        })
    return rows
