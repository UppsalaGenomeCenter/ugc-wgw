"""Inputs JSON for any (subject, stage) from config + state DB + prior-stage outputs (docs/ENTRYPOINTS.md §§3-7).

Every builder returns unprefixed input names. `generate()` applies the config's per-stage
overrides, drops None (keeps ""), prefixes `ugc_wgw_<stage>.`, and checks the result against
the entrypoint's declared inputs, because miniwdl rejects unknown keys outright.
"""
from __future__ import annotations

import re

from dataclasses import dataclass, field
from pathlib import Path

from .config import GPU_MODES, Config
from .db import DB, Cohort, RunRecord, SampleRecord
from .log import LOGGER
from .outputs import MissingOutputError, Outputs
from .stages import declared_calls, declared_inputs, namespace, sample_stages, stage_spec, wdl_path, ref_map_keys
from .util import UgcError, read_json, sha256_file, write_json


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
    allow_standalone: bool = False   # joint mode: vcf members without a gVCF are left out of cohort_call

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

    def optional_prior(self, stage: str, subject_type: str, subject_id: str) -> Outputs | None:
        """The latest success of `stage` for the subject at any ugc-wgw version, or None: a result of another mode
        that a builder uses when it exists and never requires. Recorded in member_runs like prior()."""
        run = self.db.latest_success(subject_type, subject_id, stage, None)
        if run is None:
            return None
        outs = Outputs.load(run)
        triple = {"sample_id" if subject_type == "sample" else "cohort_id": subject_id,
                  "stage": stage, "ugc_wgw_version": run.ugc_wgw_version, "run_id": run.run_id}
        if triple not in self.member_runs:
            self.member_runs.append(triple)
        LOGGER.info("%s %s: using %s run %s (ugc-wgw %s)", subject_type, subject_id, stage, run.run_id, run.ugc_wgw_version)
        return outs

    def kind(self, sid: str) -> str:
        rec = self.db.get_sample(sid)
        return (rec.entry_kind if rec else None) or "reads"

    def sample_outputs(self, sid: str) -> "MergedOutputs":
        """The outputs of a member's per-sample runs in this mode, by its entry kind (singleton; upstream and downstream;
        intake and downstream), the stages that have run so far merged (the first stage with a name wins). Raises
        MissingOutputError when none of them has run."""
        stages = sample_stages(self.mode, self.kind(sid))
        parts: list[Outputs] = []
        for stage in stages:
            run = self.db.latest_success("sample", sid, stage, None if self.any_version else self.ugc_wgw_version)
            if run is None:
                continue
            parts.append(self.prior(stage, "sample", sid)[1])
        if not parts:
            raise MissingOutputError("sample", sid, stages[0] if stages else "singleton", None, "*")
        return MergedOutputs(sid, stages, parts)

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


class MergedOutputs:
    """The outputs of a member's per-sample stages seen as one: `opt` takes the first stage that has the name."""

    def __init__(self, sample_id: str, stages: tuple[str, ...], parts: list[Outputs]):
        self.sample_id, self.stages, self.parts = sample_id, stages, parts

    @property
    def run(self) -> RunRecord:
        return self.parts[0].run

    def opt(self, name: str) -> object | None:
        for o in self.parts:
            v = o.opt(name)
            if v is not None:
                return v
        return None

    def req(self, name: str) -> object:
        v = self.opt(name)
        if v is None:
            raise MissingOutputError("sample", self.sample_id, "/".join(o.run.stage for o in self.parts) or self.stages[0], self.parts[0].run.run_id, name)
        return v

    def opt_arr(self, name: str) -> list[object] | None:
        v = self.opt(name)
        return list(v) if isinstance(v, list) else None

    def arr(self, name: str) -> list[object]:
        v = self.req(name)
        return list(v) if isinstance(v, list) else [v]


def entry_record(ctx: BuildContext) -> dict[str, object]:
    """What a sample run's manifest records about its inputs: the entry kind and, for the bam and vcf kinds, every
    external file with its size, modification time and what registration found (no content checksum: the files are
    tens of GB; the registered size and mtime identify them, as for reads)."""
    rec = ctx.sample()
    files: list[dict[str, object]] = []
    if rec.entry_kind != "reads":
        for kind, path in rec.data_inputs():
            info = rec.input_info.get(path) or {}
            entry: dict[str, object] = {"role": kind, "path": path}
            try:
                st = Path(path).stat()
                entry["bytes"], entry["mtime"] = st.st_size, st.st_mtime
            except OSError:
                entry["bytes"], entry["mtime"] = info.get("bytes"), info.get("mtime")
            if info:
                entry["info"] = {k: info.get(k) for k in ("samples", "n_ref", "reads", "reads_exact", "movies", "n_contigs", "fileformat", "sources") if k in info}
            files.append(entry)
            idx = getattr(rec, kind + "_index", None)
            if idx:
                files.append({"role": kind + "_index", "path": idx})
    out: dict[str, object] = {"kind": rec.entry_kind, "inputs": files}
    if rec.meta.get("rebase"):
        out["rebase"] = rec.meta["rebase"]
    return out


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


def _entry_bam(ctx: BuildContext) -> dict[str, object]:
    """The bam entry kind's inputs of upstream: the aligned BAM in place of reads (hifi_reads stays an empty array)."""
    rec = ctx.sample()
    if rec.entry_kind != "bam":
        return {}
    return {"aligned_bam_in": rec.aligned_bam, "aligned_bam_in_index": rec.aligned_bam_index}


# Parabricks' task hard-codes four GPUs (`Int gpuCount = 4` in parabricks_deepvariant.wdl) and no workflow
# exposes it; miniwdl accepts the call-qualified name, which is the same under ugc_wgw_singleton and ugc_wgw_upstream.
PARABRICKS_GPUS_INPUT = "upstream.parabricks_deepvariant.run_parabricks_deepvariant.gpuCount"


def _deepvariant_inputs(ctx: BuildContext) -> dict[str, object]:
    """`use_gpu`, `use_parabricks_deepvariant` and Parabricks' GPU count from config.json's `deepvariant`
    switch (guide chapter 12); `stage_inputs` still override them. The GPU type is deliberately not an
    input: `RuntimeAttributes.gpuType` reaches every task and would be part of every call-cache key, so the
    plugin applies config.json `gpu_type` at submit time instead (engine.run_env)."""
    mode = ctx.cfg.deepvariant
    doc: dict[str, object] = {"use_gpu": mode in GPU_MODES, "use_parabricks_deepvariant": mode == "parabricks"}
    if mode == "parabricks":
        doc[PARABRICKS_GPUS_INPUT] = ctx.cfg.parabricks_gpus
    if mode in ("pbdv", "pbdv-cpu"):
        doc["use_pbdv"] = True   # pblab's caller (preview): the GPU task with use_gpu, else the CPU task
    return doc


def _aligner_inputs(ctx: BuildContext) -> dict[str, object]:
    """`aligner` from config.json's key, only when it is not the default: a project on pbmm2 keeps the inputs
    document it always had."""
    return {"aligner": ctx.cfg.aligner} if ctx.cfg.aligner != "pbmm2" else {}


def tools_record(stage: str, doc: dict[str, object], kind: str) -> dict[str, str] | None:
    """Which aligner and small-variant caller a per-sample run uses, read from its final inputs document (so a
    `stage_inputs` override is told as it ran): the run manifest's `tools` block. pbdv's VCF carries DeepVariant's
    header verbatim, so nothing but this block says which of them made a call set. `external`: the files came
    with the sample (the bam and vcf entry kinds). None for the stages that neither align nor call."""
    if stage == "intake":
        return {"aligner": "external", "small_variant_caller": "external"}
    if stage not in ("singleton", "upstream"):
        return None
    def get(name: str, default: object) -> object:
        return doc.get(f"ugc_wgw_{stage}.{name}", default)
    gpu = bool(get("use_gpu", False))
    if get("use_pbdv", False):
        caller, device = "pbdv", "gpu" if gpu else "cpu"
    elif gpu and get("use_parabricks_deepvariant", False):
        caller, device = "parabricks_deepvariant", "gpu"
    else:
        caller, device = "deepvariant", "gpu" if gpu else "cpu"
    aligner = "external" if kind == "bam" else str(get("aligner", "pbmm2"))
    return {"aligner": aligner, "small_variant_caller": caller, "small_variant_device": device}


def _ref_map(ctx: BuildContext) -> dict[str, str]:
    """The rendered reference map as a dict (two-column TSV, comments and blank lines ignored)."""
    out: dict[str, str] = {}
    try:
        with open(ctx.cfg.ref_map_file) as fh:
            for line in fh:
                s = line.rstrip("\n")
                if not s or s.startswith("#") or "\t" not in s:
                    continue
                k, v = s.split("\t", 1)
                out[k] = v
    except OSError:
        pass
    return out


def _sex_inputs(ctx: BuildContext) -> dict[str, object]:
    """The sheet's sex as the first source of the resolved sex (omitted when blank), and the registry of the
    somalier image the per-sample stages run."""
    rec = ctx.sample()
    return {"declared_sex": rec.sex or None, "ugc_wgw_container_registry": ctx.cfg.ugc_wgw_container_registry}


def build_singleton(ctx: BuildContext) -> dict[str, object]:
    rec = ctx.sample()
    if rec.entry_kind != "reads":
        raise UgcError(f"singleton is for samples registered from reads; {rec.sample_id} is the {rec.entry_kind} kind "
                       f"(its standalone stages: {' -> '.join(sample_stages('standalone', rec.entry_kind))})")
    doc = common_inputs(ctx)
    doc.update(_sample_reads(ctx))
    doc.update(_deepvariant_inputs(ctx))
    doc.update(_aligner_inputs(ctx))
    doc.update(_sex_inputs(ctx))
    return doc


def build_upstream(ctx: BuildContext) -> dict[str, object]:
    rec = ctx.sample()
    if rec.entry_kind == "vcf":
        raise UgcError(f"{rec.sample_id} is the vcf kind: its calls are given, so intake runs in place of upstream")
    doc = common_inputs(ctx)
    doc.update(_sample_reads(ctx))
    doc.update(_entry_bam(ctx))
    doc.update(_deepvariant_inputs(ctx))
    doc.update(_aligner_inputs(ctx))
    doc.update(_sex_inputs(ctx))
    return doc


def build_intake(ctx: BuildContext) -> dict[str, object]:
    """The vcf entry kind: the sheet's aligned BAM and VCFs (indices resolved at registration), the declared sex."""
    rec = ctx.sample()
    if rec.entry_kind != "vcf":
        raise UgcError(f"intake is for samples registered with VCFs (the vcf kind); {rec.sample_id} is the {rec.entry_kind} kind")
    doc = common_inputs(ctx)
    doc["sample_id"] = rec.sample_id
    doc["aligned_bam"], doc["aligned_bam_index"] = rec.aligned_bam, rec.aligned_bam_index
    doc["small_variant_vcf_in"], doc["small_variant_vcf_in_index"] = rec.small_variant_vcf, rec.small_variant_vcf_index
    doc["sv_vcf_in"], doc["sv_vcf_in_index"] = rec.sv_vcf, rec.sv_vcf_index
    doc["gvcf_in"], doc["gvcf_in_index"] = rec.gvcf, rec.gvcf_index
    doc.update(_sex_inputs(ctx))
    return doc


def build_cohort_call(ctx: BuildContext) -> dict[str, object]:
    doc = common_inputs(ctx)
    members = ctx.members()
    doc["cohort_id"] = ctx.subject_id
    doc["sample_ids"] = members
    gvcfs, gvcf_indices, discover, bams, bais = [], [], [], [], []
    joint_sv = bool(ctx.ov("cohort_call", "run_sawfish_joint_call", False))
    kept: list[str] = []
    left_out: list[str] = []
    for sid in members:
        outs = ctx.sample_outputs(sid)
        gvcf = outs.opt("small_variant_gvcf")
        if gvcf is None and ctx.kind(sid) == "vcf":
            if not ctx.allow_standalone:
                raise UgcError(f"cohort_call {ctx.subject_id}: member {sid} (vcf kind) has no gVCF; register its gvcf column or submit with "
                               "--allow-standalone-members")
            left_out.append(sid)
            continue
        kept.append(sid)
        gvcfs.append(gvcf if gvcf is not None else outs.req("small_variant_gvcf"))
        gvcf_indices.append(outs.req("small_variant_gvcf_index"))
        if joint_sv:
            tar = outs.opt("discover_tar")
            if tar is None:
                raise UgcError(f"cohort_call {ctx.subject_id}: run_sawfish_joint_call needs every member's sawfish discover output; "
                               f"{sid} ({ctx.kind(sid)} kind) has none")
            discover.append(tar)
            bams.append(outs.req("aligned_hifi_reads"))
            bais.append(outs.req("aligned_hifi_reads_index"))
    if left_out:
        LOGGER.warning("cohort_call %s: %d vcf member(s) without a gVCF left out of the joint call (--allow-standalone-members): %s",
                       ctx.subject_id, len(left_out), ", ".join(left_out[:5]) + (" ..." if len(left_out) > 5 else ""))
        doc["sample_ids"] = kept
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
    head_stage = "intake" if ctx.kind(sid) == "vcf" else "upstream"
    _, up = ctx.prior(head_stage, "sample", sid)
    sex = up.opt("resolved_sex")
    if sex is None:   # an upstream run from before 0.14.0: mosdepth's inference is all there is
        sex = up.req("inferred_sex")
        LOGGER.info("downstream %s: upstream run %s predates resolved_sex; using inferred_sex %r", sid, up.run.run_id, sex)
    doc["sex"] = sex
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
        _, cc = ctx.prior("cohort_call", "cohort", ctx.cohort.cohort_id)
        called = _run_input(cc.run, "cohort_call", "sample_ids") or list(ctx.cohort.members)
        if sid in called:
            i = called.index(sid)
            doc["small_variant_vcf"] = cc.arr("split_small_variant_vcfs")[i]
            doc["small_variant_vcf_index"] = cc.arr("split_small_variant_vcf_indices")[i]
            split_sv = cc.opt_arr("split_sv_vcfs")
            if split_sv:
                doc["sv_vcf"] = split_sv[i]
                doc["sv_vcf_index"] = cc.arr("split_sv_vcf_indices")[i]
            else:
                doc["sv_vcf"] = up.req("sv_vcf")
                doc["sv_vcf_index"] = up.req("sv_vcf_index")
        else:   # a vcf member left out of the joint call (--allow-standalone-members): phased from its own VCFs
            LOGGER.info("downstream %s: not in cohort_call %s (no gVCF); phasing its own VCFs", sid, cc.run.run_id)
            doc["small_variant_vcf"] = up.req("small_variant_vcf")
            doc["small_variant_vcf_index"] = up.req("small_variant_vcf_index")
            doc["sv_vcf"] = up.req("sv_vcf")
            doc["sv_vcf_index"] = up.req("sv_vcf_index")
    else:
        doc["small_variant_vcf"] = up.req("small_variant_vcf")
        doc["small_variant_vcf_index"] = up.req("small_variant_vcf_index")
        doc["sv_vcf"] = up.req("sv_vcf")
        doc["sv_vcf_index"] = up.req("sv_vcf_index")
    return doc


def _run_input(run: RunRecord, stage: str, name: str) -> list[object] | None:
    """An array input a run was started with (its inputs.json), or None."""
    try:
        doc = read_json(Path(run.inputs_path)) if run.inputs_path else None
    except (OSError, ValueError):
        return None
    v = doc.get(f"{namespace(stage)}.{name}") if isinstance(doc, dict) else None
    return list(v) if isinstance(v, list) else None


def build_cohort_merge(ctx: BuildContext) -> dict[str, object]:
    doc = common_inputs(ctx)
    members = ctx.members()
    doc["cohort_id"] = ctx.subject_id
    doc["sample_ids"] = members
    doc["ugc_wgw_container_registry"] = ctx.cfg.ugc_wgw_container_registry
    # standalone only by default: in joint mode cohort_call already made the cohort GLnexus VCF (DESIGN §7.5)
    run_glnexus = bool(ctx.ov("cohort_merge", "run_glnexus", ctx.mode == "standalone"))
    merge_phased = bool(ctx.ov("cohort_merge", "merge_phased_small_variants", False))
    sv, svi, tr, tri, gv, gvi, ph, phi = [], [], [], [], [], [], [], []
    for sid in members:
        outs = ctx.sample_outputs(sid)   # singleton, or upstream/intake + downstream (downstream: phased_*, trgt_vcf)
        sv.append(outs.req("phased_sv_vcf"))
        svi.append(outs.req("phased_sv_vcf_index"))
        trgt = outs.opt("phased_trgt_vcf")
        tr.append(trgt if trgt is not None else outs.req("trgt_vcf"))
        tri.append(outs.req("phased_trgt_vcf_index") if trgt is not None else outs.req("trgt_vcf_index"))
        if merge_phased:
            ph.append(outs.req("phased_small_variant_vcf"))
            phi.append(outs.req("phased_small_variant_vcf_index"))
        if run_glnexus:
            g = outs.opt("small_variant_gvcf")
            if g is None:
                raise UgcError(f"cohort_merge {ctx.subject_id}: GLnexus needs a gVCF from every member; {sid} ({ctx.kind(sid)} kind) has none "
                               "(register its gvcf column, or set stage_inputs cohort_merge.run_glnexus false)")
            gv.append(g)
            gvi.append(outs.req("small_variant_gvcf_index"))
    doc["sv_vcfs"], doc["sv_vcf_indices"] = sv, svi
    doc["trgt_vcfs"], doc["trgt_vcf_indices"] = tr, tri
    doc["run_glnexus"] = run_glnexus
    if run_glnexus:
        doc["gvcfs"], doc["gvcf_indices"] = gv, gvi
    if merge_phased:
        doc["merge_phased_small_variants"] = True
        doc["phased_small_variant_vcfs"], doc["phased_small_variant_vcf_indices"] = ph, phi
    doc.update(_mito_inputs(ctx, members, ""))
    return doc


def _mito_inputs(ctx: BuildContext, members: list[str], per_stage: str) -> dict[str, object]:
    """cohort_merge's mitochondrial inputs: each member's collapsed chrM sets and haplogroups (upstream/singleton since
    0.15.0), or null plus its mitorsaw outputs for a member whose run predates the collapse (the stage collapses and
    classifies it); merge_mitochondria is on when at least one member brings something. `merge_mitochondria` set to false
    in config.json stage_inputs leaves the whole branch out."""
    if not bool(ctx.ov("cohort_merge", "merge_mitochondria", True)):
        return {"merge_mitochondria": False}
    doc: dict[str, object] = {}
    small, smalli, svs, svsi, hg, mv, mvi, mj = ([] for _ in range(8))
    legacy: list[str] = []
    without: list[str] = []
    for sid in members:
        outs = ctx.sample_outputs(sid)
        vcf = outs.opt("mt_small_variant_vcf")
        if vcf is None and outs.opt("mitorsaw_vcf") is None:   # no mitorsaw call at all: a sample entered with its own VCFs
            without.append(sid)
            continue
        if vcf is not None:
            small.append(vcf)
            smalli.append(outs.req("mt_small_variant_vcf_index"))
            svs.append(outs.req("mt_sv_vcf"))
            svsi.append(outs.req("mt_sv_vcf_index"))
            hg.append(outs.opt("mt_haplogroups"))
            mv.append(None)
            mvi.append(None)
            mj.append(None)
        else:
            legacy.append(sid)
            small.append(None)
            smalli.append(None)
            svs.append(None)
            svsi.append(None)
            hg.append(None)
            mv.append(outs.req("mitorsaw_vcf"))
            mvi.append(outs.req("mitorsaw_vcf_index"))
            mj.append(outs.req("mitorsaw_hap_stats"))
    if without:
        LOGGER.warning("cohort_merge %s: %d member(s) have no mitorsaw call (the vcf entry kind: calls made outside this project): the "
                       "mitochondrial merge is left out for this cohort: %s", ctx.subject_id, len(without), ", ".join(without[:5]) + (" ..." if len(without) > 5 else ""))
        return {"merge_mitochondria": False}
    if legacy:
        LOGGER.warning("cohort_merge %s: %d member(s) have no collapsed mitochondrial set (per-sample runs before 0.15.0): the stage "
                       "collapses and classifies their mitorsaw calls: %s", ctx.subject_id, len(legacy), ", ".join(legacy[:5]) + (" ..." if len(legacy) > 5 else ""))
    doc["merge_mitochondria"] = True
    doc["mt_small_variant_vcfs"], doc["mt_small_variant_vcf_indices"] = small, smalli
    doc["mt_sv_vcfs"], doc["mt_sv_vcf_indices"] = svs, svsi
    doc["mt_haplogroups"] = hg
    if legacy:
        doc["mitorsaw_vcfs"], doc["mitorsaw_vcf_indices"], doc["mitorsaw_hap_stats"] = mv, mvi, mj
    return doc


SEX_TO_GROUP = {"MALE": "XY", "FEMALE": "XX"}


SEX_SOURCES = ("declared", "coverage", "somalier")   # the sources of a resolved sex, in precedence order


def _sex_group(ctx: BuildContext, sid: str, inferred_stage: str) -> tuple[str, str]:
    """(group, source) for the strata: the resolved sex of the member's per-sample run (0.14+: declared, else
    coverage, else somalier, with its source); for a run from before that, the sheet's sex (`declared`), else the sex
    upstream inferred from coverage (`coverage`; an empty inferred_sex is unknown, as upstream itself treats it)."""
    rec = ctx.db.get_sample(sid)
    sheet = SEX_TO_GROUP.get((rec.sex or "").upper()) if rec else None
    del inferred_stage   # the member's own per-sample stages decide (singleton, upstream or intake)
    outs = ctx.sample_outputs(sid)   # recorded as provenance even when the sheet decides
    if sheet:   # the sheet is the first source of the resolved sex as well
        return sheet, "declared"
    resolved = outs.opt("resolved_sex")
    if resolved is not None:
        group = SEX_TO_GROUP.get(str(resolved).upper())
        source = str(outs.opt("resolved_sex_source") or "")
        return (group, source if source in SEX_SOURCES else "") if group else ("", "")
    inferred = SEX_TO_GROUP.get(str(outs.opt("inferred_sex") or "").upper())
    if inferred:
        return inferred, "coverage"
    return "", ""


def build_cohort_qc(ctx: BuildContext) -> dict[str, object]:
    """Identity QC: every member's somalier extract (null, plus the aligned BAM, for a member whose per-sample run
    predates the extract: the stage makes it), the sheet's pedigree, the resolved sex groups, and the ancestry
    switch (on when the reference map names an installed 1000 Genomes set)."""
    doc = common_inputs(ctx)
    members = ctx.members()
    doc["cohort_id"] = ctx.subject_id
    doc["sample_ids"] = members
    extracts: list[object] = []
    bams: list[object] = []
    bais: list[object] = []
    legacy: list[str] = []
    for sid in members:
        outs = ctx.sample_outputs(sid)
        som = outs.opt("somalier_extract")
        extracts.append(som)
        if som is None:
            legacy.append(sid)
            bam_name = "aligned_hifi_reads" if outs.opt("aligned_hifi_reads") is not None else "merged_haplotagged_bam"
            bams.append(outs.req(bam_name))
            bais.append(outs.req(bam_name + "_index"))
        else:
            bams.append(None)
            bais.append(None)
    doc["somalier_files"] = extracts
    if legacy:
        LOGGER.warning("cohort_qc %s: %d member(s) have no somalier_extract (per-sample runs before 0.14.0): the stage extracts "
                       "from their aligned BAMs: %s", ctx.subject_id, len(legacy), ", ".join(legacy[:5]) + (" ..." if len(legacy) > 5 else ""))
        doc["aligned_bams"] = bams
        doc["aligned_bam_indices"] = bais
    groups = [_sex_group(ctx, sid, "") for sid in members]
    doc["sample_sexes"] = [g for g, _ in groups]
    doc["sample_sex_sources"] = [src for _, src in groups]
    recs = {sid: ctx.db.get_sample(sid) for sid in members}
    doc["father_ids"] = [(recs[sid].father_id if recs[sid] else "") or "" for sid in members]
    doc["mother_ids"] = [(recs[sid].mother_id if recs[sid] else "") or "" for sid in members]
    tar = _ref_map(ctx).get("somalier_1kg_tar")
    doc["run_ancestry"] = bool(ctx.ov("cohort_qc", "run_ancestry", bool(tar) and Path(str(tar)).exists()))
    doc["ugc_wgw_container_registry"] = ctx.cfg.ugc_wgw_container_registry
    return doc


def build_cohort_freq(ctx: BuildContext) -> dict[str, object]:
    doc = common_inputs(ctx)
    members = ctx.members()
    doc["cohort_id"] = ctx.subject_id
    doc["sample_ids"] = members
    _, cm = ctx.prior("cohort_merge", "cohort", ctx.subject_id)
    doc["sv_vcf"] = cm.req("cohort_sv_vcf")
    doc["sv_vcf_index"] = cm.req("cohort_sv_vcf_index")
    # the mitochondrial sets: absent when cohort_merge ran without merge_mitochondria (or at 0.14)
    for name in ("mt_small_variant_vcf", "mt_sv_vcf"):
        doc[name] = cm.opt("cohort_" + name)
        doc[name + "_index"] = cm.opt("cohort_" + name + "_index")
    if doc["mt_small_variant_vcf"] is None:
        LOGGER.info("cohort_freq %s: cohort_merge has no mitochondrial set; the mt frequencies are skipped", ctx.subject_id)
    if ctx.mode == "joint":
        _, cc = ctx.prior("cohort_call", "cohort", ctx.subject_id)
        doc["small_variant_vcf"] = cc.req("cohort_small_variant_vcf")
        doc["small_variant_vcf_index"] = cc.req("cohort_small_variant_vcf_index")
    else:
        # optional: absent when cohort_merge ran with run_glnexus = false; the stage then makes SV frequencies only
        doc["small_variant_vcf"] = cm.opt("cohort_small_variant_vcf")
        doc["small_variant_vcf_index"] = cm.opt("cohort_small_variant_vcf_index")
    groups = [_sex_group(ctx, sid, "") for sid in members]
    doc["sample_sexes"] = [g for g, _ in groups]
    doc["sample_sex_sources"] = [src for _, src in groups]
    doc["sample_entry_kinds"] = [ctx.kind(sid) for sid in members]
    return doc


def build_cohort_tr(ctx: BuildContext) -> dict[str, object]:
    """The cohort's tandem repeats: cohort_merge's merged TRGT VCF and trgt-lps table, the members and their sex groups
    (which decide how many members a chrY locus expects and the ploidy check on chrX). No per-sample file is read."""
    doc = common_inputs(ctx)
    members = ctx.members()
    doc["cohort_id"] = ctx.subject_id
    doc["sample_ids"] = members
    _, cm = ctx.prior("cohort_merge", "cohort", ctx.subject_id)
    doc["trgt_vcf"] = cm.req("cohort_trgt_vcf")
    doc["trgt_vcf_index"] = cm.req("cohort_trgt_vcf_index")
    doc["trgt_lps"] = cm.req("cohort_trgt_lps")
    doc["sample_sexes"] = [_sex_group(ctx, sid, "")[0] for sid in members]
    return doc


LABEL_SAFE = re.compile(r"[^A-Za-z0-9_.-]+")


def methyl_label(column: str, value: object) -> str:
    """A background label from a sheet column and its value: `column=value`, the value reduced to [A-Za-z0-9_.-] (MethBat
    splits labels on ';'; the prefix keeps two columns' values apart)."""
    v = LABEL_SAFE.sub("_", str(value).strip()).strip("_")
    return f"{column}={v}" if v else ""


def build_cohort_methyl(ctx: BuildContext) -> dict[str, object]:
    """The methylation layer: every member's 5mC pileup and MethBat profile (singleton in standalone mode, downstream in
    joint), its methbat report (null for a run from before the task: the stage makes it), and the labels for the
    background: the resolved sex group as cohort_freq derives it, then the project's methyl_label_columns as
    column=value. A member whose run has no pileup (no methylation tags) is left out with a warning."""
    doc = common_inputs(ctx)
    members = ctx.members()
    doc["cohort_id"] = ctx.subject_id
    columns = list(ctx.cfg.methyl_label_columns)
    known_columns: set[str] = set()
    kept: list[str] = []
    pileups: list[object] = []
    profiles: list[object] = []
    reports: list[object] = []
    labels: list[str] = []
    legacy: list[str] = []
    skipped: list[str] = []
    for sid in members:
        outs = ctx.sample_outputs(sid)
        pileup, profile = outs.opt("cpg_pileup_bed"), outs.opt("methbat_profile")
        if pileup is None or profile is None:
            skipped.append(sid)
            continue
        kept.append(sid)
        pileups.append(pileup)
        profiles.append(profile)
        report = outs.opt("methbat_report")
        reports.append(report)
        if report is None:
            legacy.append(sid)
        parts: list[str] = []
        group, _src = _sex_group(ctx, sid, "")
        if group:
            parts.append(group)
        rec = ctx.db.get_sample(sid)
        meta = rec.meta if rec else {}
        for col in columns:
            if col in meta:
                known_columns.add(col)
                lab = methyl_label(col, meta.get(col))
                if lab:
                    parts.append(lab)
        labels.append(";".join(parts))
    unknown = [c for c in columns if c not in known_columns]
    if unknown:
        raise UgcError(f"methyl_label_columns names {', '.join(repr(c) for c in unknown)}, which no member of {ctx.subject_id} has in the "
                       "sample sheet (the sheet's extra columns are the labels available)")
    if skipped:
        LOGGER.warning("cohort_methyl %s: %d member(s) have no 5mC pileup or MethBat profile (no methylation tags in the reads): left out "
                       "of the methylation layer: %s", ctx.subject_id, len(skipped), ", ".join(skipped[:5]) + (" ..." if len(skipped) > 5 else ""))
    if not kept:
        raise UgcError(f"cohort_methyl {ctx.subject_id}: no member has a 5mC pileup and a MethBat profile")
    if legacy:
        LOGGER.warning("cohort_methyl %s: %d member(s) have no methbat_report (per-sample runs before 0.16.0): the stage makes it from "
                       "their pileups: %s", ctx.subject_id, len(legacy), ", ".join(legacy[:5]) + (" ..." if len(legacy) > 5 else ""))
    doc["sample_ids"] = kept
    doc["sample_labels"] = labels
    doc["cpg_pileup_beds"] = pileups
    doc["methbat_profiles"] = profiles
    doc["methbat_reports"] = reports
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


def build_assembly_call(ctx: BuildContext) -> dict[str, object]:
    doc = common_inputs(ctx)
    rec = ctx.sample()
    sid = rec.sample_id
    doc["sample_id"] = sid
    doc["ugc_wgw_container_registry"] = ctx.cfg.ugc_wgw_container_registry
    _, asm = ctx.prior("assembly", "sample", sid)
    doc["hap_fastas"] = asm.arr("hap_fastas")
    doc["haplotype_labels"] = asm.arr("haplotype_labels")
    doc["trio"] = bool(asm.req("trio"))
    doc["flagger_bed"] = asm.opt("flagger_bed")
    # the sample's read-based small-variant calls when it has them: singleton first, else downstream (joint mode),
    # at any version; another mode's result, used for concordance and never a prerequisite
    for stage in ("singleton", "downstream"):
        reads = ctx.optional_prior(stage, "sample", sid)
        if reads is not None and reads.opt("phased_small_variant_vcf") is not None:
            doc["read_small_variant_vcf"] = reads.req("phased_small_variant_vcf")
            doc["read_small_variant_vcf_index"] = reads.req("phased_small_variant_vcf_index")
            break
    return doc


KARYOTYPE_TO_GROUP = {"XX": "XX", "XY": "XY"}   # assembly_call's inferred_karyotype -> frequency stratum; unknown -> ""


def build_cohort_pangenome(ctx: BuildContext) -> dict[str, object]:
    doc = common_inputs(ctx)
    members = ctx.members()
    doc["cohort_id"] = ctx.subject_id
    doc["sample_ids"] = members
    hap1, hap2, trios, sexes, sources = [], [], [], [], []
    for sid in members:
        _, asm = ctx.prior("assembly", "sample", sid)
        haps = asm.arr("hap_fastas")
        if len(haps) != 2:
            raise UgcError(f"cohort_pangenome: assembly of {sid} has {len(haps)} haplotype FASTAs, not 2")
        hap1.append(haps[0])
        hap2.append(haps[1])
        trios.append(bool(asm.req("trio")))
        # the strata come from the karyotype assembly_call inferred from the haplotypes themselves, never the sheet:
        # a member of unknown karyotype is counted in the unstratified tags only
        _, call = ctx.prior("assembly_call", "sample", sid)
        group = KARYOTYPE_TO_GROUP.get(str(call.opt("inferred_karyotype") or "").strip().upper(), "")
        sexes.append(group)
        sources.append("karyotype" if group else "")
    doc["hap1_fastas"] = hap1
    doc["hap2_fastas"] = hap2
    doc["trios"] = trios
    doc["sample_sexes"] = sexes
    doc["sample_sex_sources"] = sources
    # the read-based SV merge to compare the graph with (truvari): the cohort config.json names, else this cohort's own,
    # at any version and in either read mode; used when it exists, never a prerequisite
    sv_cohort = ctx.cfg.pangenome_sv_cohort or ctx.subject_id
    cm = ctx.optional_prior("cohort_merge", "cohort", sv_cohort)
    if cm is not None and cm.opt("cohort_sv_vcf") is not None:
        doc["sv_vcf"] = cm.req("cohort_sv_vcf")
        doc["sv_vcf_index"] = cm.req("cohort_sv_vcf_index")
    elif ctx.cfg.pangenome_sv_cohort:
        LOGGER.warning("cohort_pangenome %s: no successful cohort_merge of cohort %s (config.json pangenome_sv_cohort); "
                       "the graph's SVs are not compared with a read-based merge", ctx.subject_id, sv_cohort)
    return doc


BUILDERS = {
    "singleton": build_singleton,
    "upstream": build_upstream,
    "intake": build_intake,
    "cohort_qc": build_cohort_qc,
    "cohort_call": build_cohort_call,
    "downstream": build_downstream,
    "cohort_merge": build_cohort_merge,
    "cohort_freq": build_cohort_freq,
    "cohort_tr": build_cohort_tr,
    "cohort_methyl": build_cohort_methyl,
    "assembly": build_assembly,
    "assembly_call": build_assembly_call,
    "cohort_pangenome": build_cohort_pangenome,
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
    check_ref_map(ctx, stage)
    return prefixed, list(ctx.member_runs)


def check_ref_map(ctx: BuildContext, stage: str) -> None:
    """Refuse a reference map that lacks a key the entrypoint reads (a map rendered by an older version, kept by a
    project that moved: the engine would fail late with `Map key not found`)."""
    wdl = wdl_path(ctx.cfg.code_dir, stage)
    needed = ref_map_keys(wdl)
    if not needed:
        return
    if not Path(ctx.cfg.ref_map_file).exists():
        raise UgcError(f"{stage}: the reference map {ctx.cfg.ref_map_file} does not exist (config.json ref_map_file)")
    missing = sorted(needed - set(_ref_map(ctx)))
    if missing:
        raise UgcError(f"{stage}: the reference map {ctx.cfg.ref_map_file} lacks {', '.join(missing)}, read by "
                       f"{stage_spec(stage).wdl}: the map was rendered by an older version; `ugc-wgw --project "
                       f"{ctx.cfg.project_dir} upgrade --install <root>/current` points the project at the installed "
                       f"version's map (<root>/versions/<v>/references/), or set ref_map_file in config.json to it")


def write_inputs(path: Path, doc: dict[str, object]) -> str:
    write_json(path, doc)
    return sha256_file(path)
