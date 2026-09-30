"""`ugc-wgw stage-inputs`: every input of every entrypoint, who fills it, and the project's overrides (guide chapter 05).

An entrypoint's own inputs come from its `input {` block and `parameter_meta` (stages.py, stdlib). The
call-qualified inputs of the tasks and subworkflows inside it (`--nested`: `upstream.pbmm2.pbmm2_align_wgs.threads`)
come from miniwdl itself, run in the engine venv by wdl_available_inputs.py, because only the engine knows which
call inputs a workflow leaves unbound. Who fills what is the table below; tests/driver/test_stage_inputs.py pins it
to what the builders in inputs.py actually produce.
"""
from __future__ import annotations

import json
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from .config import Config
from .engine import engine_python
from .inputs import PARABRICKS_GPUS_INPUT
from .stages import STAGES, declared_calls, declared_input_specs, input_descriptions, stage_spec, wdl_path
from .util import UgcError

DRIVER = "driver"            # built from the sheet, the cohort or earlier outputs: leave it alone
CONFIG = "config.json"       # a key of the project configuration
FREE = "stage_inputs"        # nobody: the WDL default applies unless the project sets it

Source = tuple[str, str]     # (who, from where / what a switch does)

_COMMON: dict[str, Source] = {
    "ugc_wgw_version": (DRIVER, "VERSION of the install"),
    "ref_map_file": (CONFIG, "ref_map_file"),
    "backend": (CONFIG, "backend"),
    "preemptible": (CONFIG, "preemptible"),
}
_SAMPLE: dict[str, Source] = {
    "sample_id": (DRIVER, "sample sheet"),
    "hifi_reads": (DRIVER, "sample sheet"),
    "fail_reads": (DRIVER, "sample sheet"),
}
_DEEPVARIANT: dict[str, Source] = {
    "use_gpu": (CONFIG, "deepvariant"),
    "use_parabricks_deepvariant": (CONFIG, "deepvariant"),
    "gpuType": (CONFIG, "gpu_type, GPU flavours only"),
}
_COHORT: dict[str, Source] = {
    "cohort_id": (DRIVER, "frozen cohort"),
    "sample_ids": (DRIVER, "frozen cohort, in its order"),
}
_REGISTRY: dict[str, Source] = {"ugc_wgw_container_registry": (CONFIG, "ugc_wgw_container_registry")}
_MEMBER_VCFS = "singleton (standalone) or downstream (joint) outputs of the members"

SOURCES: dict[str, dict[str, Source]] = {
    "singleton": {**_COMMON, **_SAMPLE, **_DEEPVARIANT},
    "upstream": {**_COMMON, **_SAMPLE, **_DEEPVARIANT},
    "cohort_call": {
        **_COMMON, **_COHORT,
        "gvcfs": (DRIVER, "upstream outputs of the members"),
        "gvcf_indices": (DRIVER, "upstream outputs of the members"),
        "run_sawfish_joint_call": (FREE, "switch: when true the driver also gathers discover_tars and the aligned BAMs"),
        "discover_tars": (DRIVER, "upstream outputs, when run_sawfish_joint_call"),
        "aligned_bams": (DRIVER, "upstream outputs, when run_sawfish_joint_call"),
        "aligned_bam_indices": (DRIVER, "upstream outputs, when run_sawfish_joint_call"),
    },
    "downstream": {
        **_COMMON,
        "sample_id": (DRIVER, "sample sheet"),
        "sex": (DRIVER, "upstream outputs (inferred_sex)"),
        "aligned_hifi_reads": (DRIVER, "upstream outputs"),
        "aligned_hifi_reads_index": (DRIVER, "upstream outputs"),
        "aligned_fail_reads": (DRIVER, "upstream outputs, when the sample has fail reads"),
        "aligned_fail_reads_index": (DRIVER, "upstream outputs, when the sample has fail reads"),
        "stat_depth_mean": (DRIVER, "upstream outputs"),
        "upstream_msg": (DRIVER, "upstream outputs (msg)"),
        "small_variant_vcf": (DRIVER, "cohort_call outputs, the member's slice"),
        "small_variant_vcf_index": (DRIVER, "cohort_call outputs, the member's slice"),
        "sv_vcf": (DRIVER, "cohort_call slice when it joint-called SVs, else upstream outputs"),
        "sv_vcf_index": (DRIVER, "cohort_call slice when it joint-called SVs, else upstream outputs"),
    },
    "cohort_merge": {
        **_COMMON, **_COHORT, **_REGISTRY,
        "sv_vcfs": (DRIVER, _MEMBER_VCFS),
        "sv_vcf_indices": (DRIVER, _MEMBER_VCFS),
        "trgt_vcfs": (DRIVER, _MEMBER_VCFS),
        "trgt_vcf_indices": (DRIVER, _MEMBER_VCFS),
        "run_glnexus": (FREE, "switch: the driver sets true in standalone, false in joint (cohort_call made the "
                              "cohort VCF); when true it also gathers the gVCFs"),
        "gvcfs": (DRIVER, "singleton (standalone) or upstream (joint) outputs, when run_glnexus"),
        "gvcf_indices": (DRIVER, "singleton (standalone) or upstream (joint) outputs, when run_glnexus"),
        "merge_phased_small_variants": (FREE, "switch: when true the driver also gathers the phased small-variant VCFs"),
        "phased_small_variant_vcfs": (DRIVER, "singleton or downstream outputs, when merge_phased_small_variants"),
        "phased_small_variant_vcf_indices": (DRIVER, "singleton or downstream outputs, when merge_phased_small_variants"),
    },
    "cohort_freq": {
        **_COMMON, **_COHORT,
        "sv_vcf": (DRIVER, "cohort_merge outputs"),
        "sv_vcf_index": (DRIVER, "cohort_merge outputs"),
        "small_variant_vcf": (DRIVER, "cohort_call (joint) or cohort_merge (standalone, when it ran GLnexus) outputs"),
        "small_variant_vcf_index": (DRIVER, "cohort_call (joint) or cohort_merge (standalone, when it ran GLnexus) outputs"),
        "sample_sexes": (DRIVER, "sample sheet, else upstream's inferred_sex, else empty"),
        "sample_sex_sources": (DRIVER, "sheet, inferred or empty, per member"),
    },
    "assembly": {
        **_COMMON, **_REGISTRY,
        "sample_id": (DRIVER, "sample sheet"),
        "hifi_reads": (DRIVER, "sample sheet"),
        "father_id": (DRIVER, "sample sheet, trios only (assembly_use_parents)"),
        "mother_id": (DRIVER, "sample sheet, trios only (assembly_use_parents)"),
        "father_hifi_reads": (DRIVER, "sample sheet, trios only (assembly_use_parents)"),
        "mother_hifi_reads": (DRIVER, "sample sheet, trios only (assembly_use_parents)"),
    },
}

# call-qualified inputs the driver sets (inputs.py); every other nested input is the task's own default
NESTED_SOURCES: dict[str, dict[str, Source]] = {
    "singleton": {PARABRICKS_GPUS_INPUT: (CONFIG, "parabricks_gpus, parabricks flavour only")},
    "upstream": {PARABRICKS_GPUS_INPUT: (CONFIG, "parabricks_gpus, parabricks flavour only")},
}

UNSET = "REQUIRED, unset"    # a required input nobody fills: a driver bug, never seen in a release


@dataclass
class Row:
    stage: str
    name: str
    type: str
    default: Optional[str]
    required: bool
    nested: bool
    set_by: str
    note: str
    description: str
    has_override: bool = False
    override: object = None

    def as_dict(self) -> dict[str, object]:
        d: dict[str, object] = {"stage": self.stage, "input": self.name, "type": self.type, "default": self.default,
                                "required": self.required, "nested": self.nested, "set_by": self.set_by,
                                "from": self.note, "description": self.description}
        if self.has_override:
            d["override"] = self.override
        return d

    def as_table_row(self) -> dict[str, object]:
        d: dict[str, object] = {"input": self.name, "type": self.type, "default": "-" if self.default is None else self.default,
                                "set_by": self.set_by, "from": self.note, "description": self.description}
        if self.has_override:
            d["override"] = json.dumps(self.override)
        return d


def source(stage: str, name: str, *, nested: bool) -> Source:
    table = NESTED_SOURCES if nested else SOURCES
    return table.get(stage, {}).get(name, (FREE, ""))


def own_rows(cfg: Config, stage: str) -> list[Row]:
    """The entrypoint's own inputs, in declaration order."""
    wdl = wdl_path(cfg.code_dir, stage)
    descriptions = input_descriptions(wdl)
    rows = []
    for spec in declared_input_specs(wdl):
        who, note = source(stage, spec.name, nested=False)
        if spec.required and who == FREE:
            who = UNSET
        rows.append(Row(stage, spec.name, spec.type, spec.default, spec.required, False, who, note,
                        descriptions.get(spec.name, "")))
    return rows


def available_inputs(cfg: Config, stages: list[str]) -> dict[str, list[dict[str, object]]]:
    """miniwdl's own list per stage (own and call-qualified inputs), from wdl_available_inputs.py in the engine venv."""
    python = engine_python(cfg)
    if not python.exists():
        raise UgcError(f"--nested lists the call inputs with miniwdl's WDL package, which needs the engine's Python; "
                       f"none at {python} (config.json venv_dir or the miniwdl executable's directory)")
    script = Path(__file__).with_name("wdl_available_inputs.py")
    wdls = [str(wdl_path(cfg.code_dir, s)) for s in stages]
    try:
        res = subprocess.run([str(python), str(script), *wdls], capture_output=True, text=True, timeout=600)
    except (OSError, subprocess.SubprocessError) as exc:
        raise UgcError(f"listing the call inputs with {python} failed: {exc}") from None
    if res.returncode != 0:
        tail = "\n".join(res.stderr.strip().splitlines()[-5:])
        raise UgcError(f"listing the call inputs failed ({python} exited {res.returncode}):\n{tail}")
    try:
        doc = json.loads(res.stdout)
    except json.JSONDecodeError as exc:
        raise UgcError(f"wdl_available_inputs.py wrote no JSON: {exc}") from None
    return {stage: list(doc[wdl]) for stage, wdl in zip(stages, wdls)}


def nested_rows(stage: str, listed: list[dict[str, object]], own: list[Row]) -> tuple[list[Row], list[str]]:
    """Rows for the call-qualified inputs, plus warnings where miniwdl and the driver's parser disagree on the
    entrypoint's own inputs (the parser is what `submit` checks with)."""
    rows, warnings = [], []
    by_name = {r.name: r for r in own}
    seen = set()
    for item in listed:
        name = str(item["name"])
        if "." in name:
            who, note = source(stage, name, nested=True)
            rows.append(Row(stage, name, str(item["type"]), item.get("default"), bool(item["required"]), True,  # type: ignore[arg-type]
                            who, note, str(item.get("description") or "")))
            continue
        seen.add(name)
        mine = by_name.get(name)
        if mine is None:
            warnings.append(f"{stage}: miniwdl lists input {name}, the driver's WDL parser does not (inputs.check would refuse it)")
        elif (mine.type, mine.required) != (str(item["type"]), bool(item["required"])):
            warnings.append(f"{stage}: input {name} is {item['type']} required={item['required']} to miniwdl but "
                            f"{mine.type} required={mine.required} to the driver's parser")
    for name in by_name.keys() - seen:
        warnings.append(f"{stage}: the driver's WDL parser lists input {name}, miniwdl does not")
    return rows, warnings


def check_overrides(cfg: Config, stage: str, rows: list[Row], *, nested: bool) -> list[str]:
    """Attach the project's stage_inputs values to the rows; warn about keys submit or miniwdl would refuse
    and about overrides of values the driver fills."""
    warnings: list[str] = []
    by_name = {r.name: r for r in rows}
    calls = declared_calls(wdl_path(cfg.code_dir, stage))
    wf = stage_spec(stage).wdl
    try:
        overrides = cfg.stage_overrides(stage)
    except UgcError:
        return warnings   # not an object: stage_key_warnings said so
    for key, value in overrides.items():
        row = by_name.get(key)
        if row is not None:
            row.has_override, row.override = True, value
            if row.set_by == DRIVER:
                warnings.append(f"stage_inputs.{stage}.{key} overrides a value the driver fills ({row.note})")
            continue
        if "." in key:
            first = key.split(".", 1)[0]
            if first not in calls:
                warnings.append(f"stage_inputs.{stage}.{key}: {first} is not a call of {wf}; submit would refuse it")
            elif nested:
                warnings.append(f"stage_inputs.{stage}.{key}: miniwdl lists no such input of {wf}; the run would fail at start")
            continue
        warnings.append(f"stage_inputs.{stage}.{key} is not an input of {wf}; submit would refuse it")
    return warnings


def stage_key_warnings(cfg: Config) -> list[str]:
    warnings = []
    for key, value in cfg.stage_inputs.items():
        bare = key[len("ugc_wgw_"):] if key.startswith("ugc_wgw_") else key
        if bare not in STAGES:
            warnings.append(f"stage_inputs key {key!r} names no stage (known: {', '.join(STAGES)})")
        elif not isinstance(value, dict):
            warnings.append(f"stage_inputs.{key} must be an object of input: value, not {type(value).__name__}")
    return warnings


def listing(cfg: Config, stages: list[str], *, nested: bool = False) -> tuple[list[Row], list[str]]:
    """Rows for the stages, in stage then declaration order, and the warnings about the project's overrides."""
    rows: list[Row] = []
    warnings = stage_key_warnings(cfg)
    listed = available_inputs(cfg, stages) if nested else {}
    for stage in stages:
        own = own_rows(cfg, stage)
        stage_rows = list(own)
        if nested:
            more, w = nested_rows(stage, listed[stage], own)
            stage_rows += more
            warnings += w
        for r in own:
            if r.set_by == UNSET:
                warnings.append(f"{stage}: required input {r.name} is not filled by the driver (bug: submit would fail)")
        warnings += check_overrides(cfg, stage, stage_rows, nested=nested)
        rows += stage_rows
    return rows, warnings


def summary(rows: list[Row], stages: list[str], *, nested: bool) -> str:
    own = [r for r in rows if not r.nested]
    n = {who: sum(1 for r in own if r.set_by == who) for who in (DRIVER, CONFIG, FREE)}
    text = (f"{len(stages)} stage{'s' if len(stages) != 1 else ''}, {len(own)} workflow inputs: {n[DRIVER]} filled by the "
            f"driver, {n[CONFIG]} from config.json keys, {n[FREE]} yours under stage_inputs")
    if nested:
        text += f"; {sum(1 for r in rows if r.nested)} call-qualified inputs of the tasks inside"
    else:
        text += " (--nested adds the call-qualified inputs of the tasks inside: thread counts, memory, tool options)"
    return text


COLUMNS_HELP = ("columns: default = the WDL default (- = none); set_by = who fills the value: driver (from the sheet, the "
                "cohort or earlier outputs: leave it), config.json (a key of the project config), stage_inputs (nobody: "
                "the default applies unless config.json's stage_inputs sets it); from = where the driver takes it, or "
                "what a switch does; override = the project's current stage_inputs value")
