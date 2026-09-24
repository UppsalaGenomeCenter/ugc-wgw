"""Inputs JSON for any (subject, stage) from config + state DB + prior-stage outputs (docs/ENTRYPOINTS.md §§3-7).

Every builder returns unprefixed input names. `generate()` applies the config's per-stage
overrides, drops None (keeps ""), prefixes `ugc_wgw_<stage>.`, and checks the result against
the entrypoint's declared inputs, because miniwdl rejects unknown keys outright.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from .config import Config
from .db import DB, Cohort, RunRecord, SampleRecord
from .log import LOGGER
from .outputs import MissingOutputError, Outputs
from .stages import declared_calls, declared_inputs, namespace, stage_spec, wdl_path
from .util import UgcError, sha256_file, write_json


@dataclass
class BuildContext:
    cfg: Config
    db: DB
    mode: str
    ugc_wgw_version: str
    subject_type: str
    subject_id: str
    cohort: Cohort | None = None
    any_version: bool = False
    member_runs: list[dict[str, object]] = field(default_factory=list)

    def prior(self, stage: str, subject_type: str, subject_id: str) -> tuple[RunRecord, Outputs]:
        run = self.db.latest_success(subject_type, subject_id, stage, None if self.any_version else self.ugc_wgw_version)
        if run is None:
            raise MissingOutputError(subject_type, subject_id, stage, None, "*")
        outs = Outputs.load(run)
        triple = {"sample_id" if subject_type == "sample" else "cohort_id": subject_id,
                  "stage": stage, "ugc_wgw_version": run.ugc_wgw_version, "run_id": run.run_id}
        if triple not in self.member_runs:
            self.member_runs.append(triple)
        return run, outs

    def sample(self) -> SampleRecord:
        rec = self.db.get_sample(self.subject_id)
        if rec is None:
            raise UgcError(f"unknown sample {self.subject_id}")
        return rec

    def members(self) -> list[str]:
        if self.cohort is None:
            raise UgcError(f"stage needs a frozen cohort (--cohort) for {self.subject_type} {self.subject_id}")
        return list(self.cohort.members)

    def ov(self, stage: str, name: str, default: object) -> object:
        return self.cfg.stage_overrides(stage).get(name, default)


def common_inputs(ctx: BuildContext) -> dict[str, object]:
    return {
        "ugc_wgw_version": ctx.ugc_wgw_version,
        "ref_map_file": str(ctx.cfg.ref_map_file),
        "backend": ctx.cfg.backend,
        "preemptible": ctx.cfg.preemptible,
    }


def _sample_reads(ctx: BuildContext) -> dict[str, object]:
    rec = ctx.sample()
    return {
        "sample_id": rec.sample_id,
        "hifi_reads": list(rec.hifi_reads),
        "fail_reads": list(rec.fail_reads) or None,
    }


# Parabricks' task hard-codes four GPUs (`Int gpuCount = 4` in parabricks_deepvariant.wdl) and no workflow
# exposes it; miniwdl accepts the call-qualified name, which is the same under ugc_wgw_singleton and ugc_wgw_upstream.
PARABRICKS_GPUS_INPUT = "upstream.parabricks_deepvariant.run_parabricks_deepvariant.gpuCount"


def _deepvariant_inputs(ctx: BuildContext) -> dict[str, object]:
    """`use_gpu`, `use_parabricks_deepvariant`, `gpuType` and Parabricks' GPU count from config.json's
    `deepvariant` switch (guide chapter 12); `stage_inputs` still override them."""
    mode = ctx.cfg.deepvariant
    doc: dict[str, object] = {"use_gpu": mode != "cpu", "use_parabricks_deepvariant": mode == "parabricks"}
    if mode != "cpu" and ctx.cfg.gpu_type:
        doc["gpuType"] = ctx.cfg.gpu_type
    if mode == "parabricks":
        doc[PARABRICKS_GPUS_INPUT] = ctx.cfg.parabricks_gpus
    return doc


def build_singleton(ctx: BuildContext) -> dict[str, object]:
    doc = common_inputs(ctx)
    doc.update(_sample_reads(ctx))
    doc.update(_deepvariant_inputs(ctx))
    return doc


def build_upstream(ctx: BuildContext) -> dict[str, object]:
    doc = common_inputs(ctx)
    doc.update(_sample_reads(ctx))
    doc.update(_deepvariant_inputs(ctx))
    return doc


def build_cohort_call(ctx: BuildContext) -> dict[str, object]:
    doc = common_inputs(ctx)
    members = ctx.members()
    doc["cohort_id"] = ctx.subject_id
    doc["sample_ids"] = members
    gvcfs, gvcf_indices, discover, bams, bais = [], [], [], [], []
    joint_sv = bool(ctx.ov("cohort_call", "run_sawfish_joint_call", False))
    for sid in members:
        _, outs = ctx.prior("upstream", "sample", sid)
        gvcfs.append(outs.req("small_variant_gvcf"))
        gvcf_indices.append(outs.req("small_variant_gvcf_index"))
        if joint_sv:
            discover.append(outs.req("discover_tar"))
            bams.append(outs.req("aligned_hifi_reads"))
            bais.append(outs.req("aligned_hifi_reads_index"))
    doc["gvcfs"] = gvcfs
    doc["gvcf_indices"] = gvcf_indices
    if joint_sv:
        doc["run_sawfish_joint_call"] = True
        doc["discover_tars"] = discover
        doc["aligned_bams"] = bams
        doc["aligned_bam_indices"] = bais
    return doc


def build_downstream(ctx: BuildContext) -> dict[str, object]:
    doc = common_inputs(ctx)
    sid = ctx.subject_id
    doc["sample_id"] = sid
    _, up = ctx.prior("upstream", "sample", sid)
    doc["sex"] = up.req("inferred_sex")
    doc["aligned_hifi_reads"] = up.req("aligned_hifi_reads")
    doc["aligned_hifi_reads_index"] = up.req("aligned_hifi_reads_index")
    doc["aligned_fail_reads"] = up.opt("aligned_fail_reads")
    doc["aligned_fail_reads_index"] = up.opt("aligned_fail_reads_index")
    doc["stat_depth_mean"] = up.opt("stat_depth_mean")
    doc["upstream_msg"] = up.opt_arr("msg") or []
    if ctx.mode == "joint":
        if ctx.cohort is None:
            raise UgcError("downstream in joint mode needs --cohort (the cohort_call context)")
        if sid not in ctx.cohort.members:
            raise UgcError(f"sample {sid} is not a member of cohort {ctx.cohort.cohort_id}")
        i = ctx.cohort.members.index(sid)
        _, cc = ctx.prior("cohort_call", "cohort", ctx.cohort.cohort_id)
        doc["small_variant_vcf"] = cc.arr("split_small_variant_vcfs")[i]
        doc["small_variant_vcf_index"] = cc.arr("split_small_variant_vcf_indices")[i]
        split_sv = cc.opt_arr("split_sv_vcfs")
        if split_sv:
            doc["sv_vcf"] = split_sv[i]
            doc["sv_vcf_index"] = cc.arr("split_sv_vcf_indices")[i]
        else:
            doc["sv_vcf"] = up.req("sv_vcf")
            doc["sv_vcf_index"] = up.req("sv_vcf_index")
    else:
        doc["small_variant_vcf"] = up.req("small_variant_vcf")
        doc["small_variant_vcf_index"] = up.req("small_variant_vcf_index")
        doc["sv_vcf"] = up.req("sv_vcf")
        doc["sv_vcf_index"] = up.req("sv_vcf_index")
    return doc


def build_cohort_merge(ctx: BuildContext) -> dict[str, object]:
    doc = common_inputs(ctx)
    members = ctx.members()
    doc["cohort_id"] = ctx.subject_id
    doc["sample_ids"] = members
    doc["ugc_wgw_container_registry"] = ctx.cfg.ugc_wgw_container_registry
    if ctx.mode == "standalone":
        sv_stage, trgt_name, gvcf_stage = "singleton", "phased_trgt_vcf", "singleton"
    else:
        sv_stage, trgt_name, gvcf_stage = "downstream", "trgt_vcf", "upstream"
    # standalone only by default: in joint mode cohort_call already made the cohort GLnexus VCF (DESIGN §7.5)
    run_glnexus = bool(ctx.ov("cohort_merge", "run_glnexus", ctx.mode == "standalone"))
    merge_phased = bool(ctx.ov("cohort_merge", "merge_phased_small_variants", False))
    sv, svi, tr, tri, gv, gvi, ph, phi = [], [], [], [], [], [], [], []
    for sid in members:
        _, outs = ctx.prior(sv_stage, "sample", sid)
        sv.append(outs.req("phased_sv_vcf"))
        svi.append(outs.req("phased_sv_vcf_index"))
        tr.append(outs.req(trgt_name))
        tri.append(outs.req(trgt_name + "_index"))
        if merge_phased:
            ph.append(outs.req("phased_small_variant_vcf"))
            phi.append(outs.req("phased_small_variant_vcf_index"))
        if run_glnexus:
            g = outs if gvcf_stage == sv_stage else ctx.prior(gvcf_stage, "sample", sid)[1]
            gv.append(g.req("small_variant_gvcf"))
            gvi.append(g.req("small_variant_gvcf_index"))
    doc["sv_vcfs"], doc["sv_vcf_indices"] = sv, svi
    doc["trgt_vcfs"], doc["trgt_vcf_indices"] = tr, tri
    doc["run_glnexus"] = run_glnexus
    if run_glnexus:
        doc["gvcfs"], doc["gvcf_indices"] = gv, gvi
    if merge_phased:
        doc["merge_phased_small_variants"] = True
        doc["phased_small_variant_vcfs"], doc["phased_small_variant_vcf_indices"] = ph, phi
    return doc


SEX_TO_GROUP = {"MALE": "XY", "FEMALE": "XX"}


def _sex_group(ctx: BuildContext, sid: str, inferred_stage: str) -> tuple[str, str]:
    """(group, source) for the frequency strata: the sheet's sex first, else the sex upstream inferred from
    coverage (an empty inferred_sex is unknown, as upstream itself treats it), else unknown."""
    rec = ctx.db.get_sample(sid)
    sheet = SEX_TO_GROUP.get((rec.sex or "").upper()) if rec else None
    _, outs = ctx.prior(inferred_stage, "sample", sid)
    inferred = SEX_TO_GROUP.get(str(outs.opt("inferred_sex") or "").upper())
    if sheet and inferred and sheet != inferred:
        LOGGER.warning("cohort_freq: sample %s is %s in the sheet but %s was inferred from coverage; using the sheet",
                       sid, sheet, inferred)
    if sheet:
        return sheet, "sheet"
    if inferred:
        return inferred, "inferred"
    return "", ""


def build_cohort_freq(ctx: BuildContext) -> dict[str, object]:
    doc = common_inputs(ctx)
    members = ctx.members()
    doc["cohort_id"] = ctx.subject_id
    doc["sample_ids"] = members
    _, cm = ctx.prior("cohort_merge", "cohort", ctx.subject_id)
    doc["sv_vcf"] = cm.req("cohort_sv_vcf")
    doc["sv_vcf_index"] = cm.req("cohort_sv_vcf_index")
    if ctx.mode == "joint":
        _, cc = ctx.prior("cohort_call", "cohort", ctx.subject_id)
        doc["small_variant_vcf"] = cc.req("cohort_small_variant_vcf")
        doc["small_variant_vcf_index"] = cc.req("cohort_small_variant_vcf_index")
        inferred_stage = "upstream"
    else:
        # optional: absent when cohort_merge ran with run_glnexus = false; the stage then makes SV frequencies only
        doc["small_variant_vcf"] = cm.opt("cohort_small_variant_vcf")
        doc["small_variant_vcf_index"] = cm.opt("cohort_small_variant_vcf_index")
        inferred_stage = "singleton"
    groups = [_sex_group(ctx, sid, inferred_stage) for sid in members]
    doc["sample_sexes"] = [g for g, _ in groups]
    doc["sample_sex_sources"] = [src for _, src in groups]
    return doc


def _trio_parents(ctx: BuildContext, rec: SampleRecord) -> tuple[SampleRecord, SampleRecord] | None:
    """Both parents registered with reads and distinct from each other and the child -> (father, mother); else None."""
    if not rec.father_id and not rec.mother_id:
        return None
    if not (rec.father_id and rec.mother_id):
        LOGGER.warning("assembly %s: only one parent given (%s/%s); running sample mode",
                       rec.sample_id, rec.father_id, rec.mother_id)
        return None
    if len({rec.sample_id, rec.father_id, rec.mother_id}) != 3:
        LOGGER.warning("assembly %s: pedigree is not a trio (father %s, mother %s); running sample mode",
                       rec.sample_id, rec.father_id, rec.mother_id)
        return None
    father, mother = ctx.db.get_sample(rec.father_id), ctx.db.get_sample(rec.mother_id)
    for label, parent, pid in (("father", father, rec.father_id), ("mother", mother, rec.mother_id)):
        if parent is None or not parent.hifi_reads:
            LOGGER.warning("assembly %s: %s %s is not a registered sample with reads; running sample mode",
                           rec.sample_id, label, pid)
            return None
    assert father is not None and mother is not None
    return father, mother


def build_assembly(ctx: BuildContext) -> dict[str, object]:
    doc = common_inputs(ctx)
    rec = ctx.sample()
    doc["sample_id"] = rec.sample_id
    doc["hifi_reads"] = list(rec.hifi_reads)
    doc["ugc_wgw_container_registry"] = ctx.cfg.ugc_wgw_container_registry
    if ctx.cfg.assembly_use_parents:
        parents = _trio_parents(ctx, rec)
        if parents is not None:
            father, mother = parents
            doc["father_id"] = father.sample_id
            doc["mother_id"] = mother.sample_id
            doc["father_hifi_reads"] = list(father.hifi_reads)
            doc["mother_hifi_reads"] = list(mother.hifi_reads)
    return doc


BUILDERS = {
    "singleton": build_singleton,
    "upstream": build_upstream,
    "cohort_call": build_cohort_call,
    "downstream": build_downstream,
    "cohort_merge": build_cohort_merge,
    "cohort_freq": build_cohort_freq,
    "assembly": build_assembly,
}


def check(code_dir: Path, stage: str, doc: dict[str, object]) -> None:
    """Refuse to hand miniwdl a document it would reject: unknown keys or missing required inputs."""
    declared = declared_inputs(wdl_path(code_dir, stage))
    ns = namespace(stage) + "."
    names = {}
    for k in doc:
        if not k.startswith(ns):
            raise UgcError(f"{stage}: input key {k!r} is not namespaced {ns}")
        names[k[len(ns):]] = k
    calls = declared_calls(wdl_path(code_dir, stage))
    # a dotted name is a nested call input (`<call>.<...>.<input>`), which miniwdl resolves through the
    # call graph and checks at start; here only the first segment is checked against the entrypoint's calls
    unknown = sorted(n for n in names if n not in declared and not ("." in n and n.split(".", 1)[0] in calls))
    if unknown:
        raise UgcError(f"{stage}: inputs not declared by {stage_spec(stage).wdl}: {', '.join(unknown)}")
    missing = sorted(n for n, required in declared.items() if required and n not in names)
    if missing:
        raise UgcError(f"{stage}: required inputs missing: {', '.join(missing)}")


def generate(ctx: BuildContext, stage: str) -> tuple[dict[str, object], list[dict[str, object]]]:
    spec = stage_spec(stage)
    if spec.subject_type != ctx.subject_type:
        raise UgcError(f"stage {stage} runs per {spec.subject_type}, not per {ctx.subject_type}")
    ctx.member_runs = []
    doc = BUILDERS[stage](ctx)
    doc.update(ctx.cfg.stage_overrides(stage))
    doc = {k: v for k, v in doc.items() if v is not None}
    prefixed = {f"{namespace(stage)}.{k}": v for k, v in doc.items()}
    check(ctx.cfg.code_dir, stage, prefixed)
    return prefixed, list(ctx.member_runs)


def write_inputs(path: Path, doc: dict[str, object]) -> str:
    write_json(path, doc)
    return sha256_file(path)
