#!/usr/bin/env python3
"""Fake `miniwdl` for driver tests: reproduces miniwdl 1.15.0's run-directory conventions without running WDL.

Accepts `miniwdl run <wdl> -i INPUTS --dir DIR/. -o RUN.json [--cfg F] [--no-color] [--log-json]` and
`miniwdl --version`. Validates the inputs the way miniwdl does (namespaced keys, only declared inputs,
required inputs present; exit 2 otherwise), writes workflow.log, then either succeeds (out/, outputs.json,
run.json with "outputs", exit 0) or fails (error.json, run.json with "error", exit 1).

Environment knobs:
  UGC_WGW_FAKE_FAIL=id1,id2   subjects that fail (ids ending in _FAIL always fail)
  UGC_WGW_FAKE_SLEEP=SECONDS  wait before finishing
  UGC_WGW_FAKE_HANG=1         ignore completion; on SIGTERM write an Interrupted error and exit 143
  UGC_WGW_FAKE_VERSION=STR    echo this as ugc_wgw_workflow_version instead of the input
  UGC_WGW_FAKE_TRACE=FILE     append "start <id> <t>" / "end <id> <t>" lines
  UGC_WGW_FAKE_FAIL_CLASS=<class>[:<exit_status>]   error class of a failure (default CommandFailed:1); CommandFailed
                          mimics miniwdl 1.15 (no message, exit_status/stderr_file/node/dir, top-level from_dir)
  UGC_WGW_FAKE_FAIL_ATTEMPTS=1,2   fail only when the attempt directory is attempt-<n> for n in the list
  UGC_WGW_FAKE_SLURM_LOG=TEXT write TEXT (\n decoded) to <attempt>/call-fake/slurm_singularity.log.txt on every run
  UGC_WGW_FAKE_STDERR=TEXT    write TEXT to <attempt>/call-fake/stderr.txt
  UGC_WGW_FAKE_IGNORE_TERM=1  with UGC_WGW_FAKE_HANG: ignore SIGTERM so the driver has to SIGKILL
  UGC_WGW_FAKE_INFERRED_SEX=S1=MALE,S2=FEMALE   inferred_sex per sample for singleton/upstream (default ""); a seed
                          may also carry a non-WDL "inferred_sex" key (tests/driver/helpers.py seed_success)
"""
from __future__ import annotations

import argparse
import gzip
import json
import math
import os
import pathlib
import signal
import sys
import time

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2] / "bin"))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from ugc_wgw.stages import declared_calls, declared_inputs  # noqa: E402
from fixtures import fixture_bytes, stats_values  # noqa: E402

REF = "GRCh38_GIABv3"


def _trace(line: str) -> None:
    path = os.environ.get("UGC_WGW_FAKE_TRACE")
    if path:
        with open(path, "a") as fh:
            fh.write(f"{line} {time.time():.3f}\n")


def _inferred_sex(sid: str, v: dict[str, object]) -> str:
    """mosdepth's inference, faked: the seed's knob, else the UGC_WGW_FAKE_INFERRED_SEX table, else unknown."""
    if "inferred_sex" in v:
        return str(v["inferred_sex"] or "")
    table = dict(item.split("=", 1) for item in os.environ.get("UGC_WGW_FAKE_INFERRED_SEX", "").split(",") if "=" in item)
    return table.get(sid, "")


def _write_json(path: pathlib.Path, doc: object) -> None:
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(doc, indent=2, sort_keys=True) + "\n")
    os.replace(tmp, path)


def make_outputs(stage: str, inputs: dict[str, object], run_dir: pathlib.Path,
                 ugc_wgw_version: str | None = None) -> dict[str, object]:
    """Create out/<name>/... files and return the namespaced outputs dict (absolute paths), like miniwdl."""
    ns = f"ugc_wgw_{stage}."
    v = {k[len(ns):]: val for k, val in inputs.items() if k.startswith(ns)}
    sid = str(v.get("sample_id") or v.get("cohort_id"))
    files: dict[str, object] = {}
    strings: dict[str, object] = {}

    def f(name: str, filename: str) -> None:
        files[name] = filename
        files[name + "_index"] = filename + (".bai" if filename.endswith(".bam") else ".tbi")

    if stage == "singleton":
        f("merged_haplotagged_bam", f"{sid}.{REF}.haplotagged.bam")
        f("phased_small_variant_vcf", f"{sid}.{REF}.small_variants.phased.vcf.gz")
        f("phased_sv_vcf", f"{sid}.{REF}.structural_variants.phased.vcf.gz")
        f("phased_trgt_vcf", f"{sid}.{REF}.trgt.sorted.vcf.gz")
        f("small_variant_gvcf", f"{sid}.{REF}.small_variants.g.vcf.gz")
        f("cpg_pileup_bed", f"{sid}.{REF}.haplotagged.cpg.bed.gz")
        f("kivvi_kiv2_vcf", f"{sid}.{REF}.kivvi.kiv2.vcf.gz")
        f("mosdepth_region_bed", f"{sid}.{REF}.mosdepth.regions.bed.gz")
        files["read_quality_plot"] = f"{sid}.{REF}.read_quality.png"
        files["sv_stats_plot"] = f"{sid}.{REF}.sv_stats.png"
        files["pbstarphase_tsv"] = f"{sid}.{REF}.pbstarphase.pharmcat.tsv"
        files["stats_file"] = f"{sid}.stats.txt"
        files["msg_file"] = f"{sid}.messages.txt"
        for name, fname in (("mosdepth_summary", f"{sid}.{REF}.mosdepth.summary.txt"), ("sv_copynum_summary", f"{sid}.{REF}.copynum.summary.json"),
                            ("sv_copynum_bedgraph", f"{sid}.{REF}.copynum.bedgraph"), ("small_variant_stats", f"{sid}.{REF}.small_variants.vcf.stats.txt"),
                            ("bcftools_roh_out", f"{sid}.{REF}.bcftools_roh.out.gz"), ("phase_stats", f"{sid}.{REF}.hiphase.stats.tsv"),
                            ("phase_blocks", f"{sid}.{REF}.hiphase.blocks.tsv"), ("trgt_coverage_dropouts", f"{sid}.{REF}.trgt.dropouts.txt"),
                            ("methbat_profile", f"{sid}.{REF}.methbat.profile.tsv"), ("paraphase_summary", f"{sid}.paraphase.json"),
                            ("mitorsaw_hap_stats", f"{sid}.{REF}.mitorsaw.json"), ("pbstarphase_summary", f"{sid}.pbstarphase.json"),
                            ("kivvi_kiv2_json", f"{sid}.{REF}.kivvi.kiv2.json"), ("kivvi_d4z4_json", f"{sid}.{REF}.kivvi.d4z4.json")):
            files[name] = fname
        inferred = _inferred_sex(sid, v)
        strings.update({"inferred_sex": inferred, "msg": []})
        strings.update({"stat_" + k: val for k, val in stats_values(sid, inferred).items() if k not in ("sample_id", "inferred_sex")})
        strings["stat_depth_mean"] = "30.1"
    elif stage == "upstream":
        f("aligned_hifi_reads", f"{sid}.{REF}.hifi_reads.bam")
        if v.get("fail_reads"):
            f("aligned_fail_reads", f"{sid}.{REF}.fail_reads.bam")
        f("small_variant_vcf", f"{sid}.{REF}.small_variants.vcf.gz")
        f("small_variant_gvcf", f"{sid}.{REF}.small_variants.g.vcf.gz")
        if v.get("run_sawfish_call", True):
            f("sv_vcf", f"{sid}.{REF}.structural_variants.vcf.gz")
        files["discover_tar"] = f"{sid}.tar"
        f("kivvi_kiv2_vcf", f"{sid}.{REF}.kivvi.kiv2.vcf.gz")
        f("mosdepth_region_bed", f"{sid}.{REF}.mosdepth.regions.bed.gz")
        for name, fname in (("mosdepth_summary", f"{sid}.{REF}.mosdepth.summary.txt"), ("sv_copynum_summary", f"{sid}.{REF}.copynum.summary.json"),
                            ("sv_copynum_bedgraph", f"{sid}.{REF}.copynum.bedgraph"), ("paraphase_output_json", f"{sid}.paraphase.json"),
                            ("mitorsaw_hap_stats", f"{sid}.{REF}.mitorsaw.json"), ("kivvi_kiv2_json", f"{sid}.{REF}.kivvi.kiv2.json"),
                            ("kivvi_d4z4_json", f"{sid}.{REF}.kivvi.d4z4.json")):
            files[name] = fname
        inferred = _inferred_sex(sid, v)
        strings.update({"inferred_sex": inferred, "stat_depth_mean": "30.1", "msg": [f"{sid}: fake upstream message"]})
    elif stage == "cohort_call":
        members = list(v.get("sample_ids", []))  # type: ignore[arg-type]
        f("cohort_small_variant_vcf", f"{sid}.joint.{REF}.small_variants.vcf.gz")
        files["split_small_variant_vcfs"] = [f"{m}.{sid}.joint.{REF}.small_variants.vcf.gz" for m in members]
        files["split_small_variant_vcf_indices"] = [x + ".tbi" for x in files["split_small_variant_vcfs"]]  # type: ignore[union-attr]
        if v.get("run_sawfish_joint_call"):
            f("cohort_sv_vcf", f"{sid}.joint.{REF}.structural_variants.vcf.gz")
            files["split_sv_vcfs"] = [f"{m}.{sid}.joint.{REF}.structural_variants.vcf.gz" for m in members]
            files["split_sv_vcf_indices"] = [x + ".tbi" for x in files["split_sv_vcfs"]]  # type: ignore[union-attr]
    elif stage == "downstream":
        f("merged_haplotagged_bam", f"{sid}.{REF}.haplotagged.bam")
        f("phased_small_variant_vcf", f"{sid}.{REF}.small_variants.phased.vcf.gz")
        f("phased_sv_vcf", f"{sid}.{REF}.structural_variants.phased.vcf.gz")
        f("trgt_vcf", f"{sid}.{REF}.trgt.sorted.vcf.gz")
        f("cpg_pileup_bed", f"{sid}.{REF}.haplotagged.cpg.bed.gz")
        files["read_quality_plot"] = f"{sid}.{REF}.read_quality.png"
        files["sv_stats_plot"] = f"{sid}.{REF}.sv_stats.png"
        files["pbstarphase_tsv"] = f"{sid}.{REF}.pbstarphase.pharmcat.tsv"
        files["stats_file"] = f"{sid}.stats.txt"
        files["msg_file"] = f"{sid}.messages.txt"
        for name, fname in (("small_variant_stats", f"{sid}.{REF}.small_variants.vcf.stats.txt"), ("bcftools_roh_out", f"{sid}.{REF}.bcftools_roh.out.gz"),
                            ("phase_stats", f"{sid}.{REF}.hiphase.stats.tsv"), ("phase_blocks", f"{sid}.{REF}.hiphase.blocks.tsv"),
                            ("trgt_coverage_dropouts", f"{sid}.{REF}.trgt.dropouts.txt"), ("methbat_profile", f"{sid}.{REF}.methbat.profile.tsv"),
                            ("pbstarphase_json", f"{sid}.pbstarphase.json")):
            files[name] = fname
        strings.update({"msg": list(v.get("upstream_msg", []))})  # type: ignore[arg-type]
        strings.update({"stat_" + k: val for k, val in stats_values(sid, str(v.get("sex") or "")).items() if k not in ("sample_id", "inferred_sex")})
    elif stage == "assembly":
        trio = bool(v.get("father_hifi_reads")) and bool(v.get("mother_hifi_reads"))
        tag = "dip" if trio else "bp"
        haps = ["hap1", "hap2"]
        files["assembly_noseq_gfas"] = [f"{sid}.asm.{tag}.{h}.p_ctg.noseq.gfa" for h in haps] + \
            [f"{sid}.asm.{tag}.{u}_utg.noseq.gfa" for u in ("p", "r")]
        files["assembly_lowQ_beds"] = [f"{sid}.asm.{tag}.{h}.p_ctg.lowQ.bed" for h in haps] + \
            [f"{sid}.asm.{tag}.{u}_utg.lowQ.bed" for u in ("p", "r")]
        files["zipped_assembly_fastas"] = [f"{sid}.asm.{tag}.{h}.p_ctg.fasta.gz" for h in haps]
        files["assembly_stats"] = [f"{sid}.asm.{tag}.{h}.p_ctg.fasta.stats.txt" for h in haps]
        files["asm_bams"] = [f"{sid}.{h}.asm.{REF}.bam" for h in haps]
        files["asm_bam_indices"] = [x + ".bai" for x in files["asm_bams"]]  # type: ignore[union-attr]
        f("merged_asm_bam", f"{sid}.asm.{REF}.bam")
        files["paftools_vcfs"] = [f"{sid}.{h}.asm.{REF}.paftools.vcf.gz" for h in haps]
        files["paftools_vcf_indices"] = [x + ".tbi" for x in files["paftools_vcfs"]]  # type: ignore[union-attr]
        files["paftools_vcf_stats"] = [x + ".stats.txt" for x in files["paftools_vcfs"]]  # type: ignore[union-attr]
        strings.update({"trio": trio, "haplotypes": haps})
        if trio:
            strings.update({"hap1_parent": v.get("father_id"), "hap2_parent": v.get("mother_id")})
    elif stage == "cohort_merge":
        f("cohort_sv_vcf", f"{sid}.merged.{REF}.structural_variants.vcf.gz")
        f("cohort_trgt_vcf", f"{sid}.merged.{REF}.trgt.vcf.gz")
        files["cohort_trgt_lps"] = f"{sid}.merged.{REF}.trgt.lps.tsv"
        if v.get("run_glnexus", True):
            f("cohort_small_variant_vcf", f"{sid}.merged.{REF}.small_variants.vcf.gz")
    elif stage == "cohort_freq":
        if v.get("small_variant_vcf"):
            f("small_variant_freq_vcf", f"{sid}.freq.{REF}.small_variants.vcf.gz")
        f("sv_freq_vcf", f"{sid}.freq.{REF}.structural_variants.vcf.gz")
        files["freq_summary"] = f"{sid}.freq.{REF}.summary.tsv"
        files["freq_samples"] = f"{sid}.freq.{REF}.samples.tsv"
    else:
        raise SystemExit(f"fake miniwdl: unknown stage {stage}")

    manifest_name = f"{sid}.{stage}.ugc_wgw_manifest.json"
    files["ugc_wgw_manifest"] = manifest_name
    out = run_dir / "out"
    outputs: dict[str, object] = {}
    members_of = [str(m) for m in v.get("sample_ids", [])] if isinstance(v.get("sample_ids"), list) else []  # type: ignore[union-attr]
    inferred_for_fixture = str(strings.get("inferred_sex") or "")

    def write_fixture(p: pathlib.Path, name: str, index: int, placeholder: str) -> None:
        content = fixture_bytes(name, sid, stage, inferred_for_fixture, index=index, members=members_of)
        if content is None:
            p.write_text(placeholder)
        elif p.name.endswith(".gz"):
            p.write_bytes(gzip.compress(content))
        else:
            p.write_bytes(content)

    for name, val in files.items():
        if isinstance(val, list):
            width = max(1, math.ceil(math.log10(len(val)))) if len(val) > 1 else 1
            paths = []
            for i, fname in enumerate(val):
                p = out / name / str(i).rjust(width, "0") / fname
                p.parent.mkdir(parents=True, exist_ok=True)
                write_fixture(p, name, i, f"fake {name}[{i}] for {sid}\n")
                paths.append(str(p))
            outputs[ns + name] = paths
        else:
            p = out / name / str(val)
            p.parent.mkdir(parents=True, exist_ok=True)
            if name == "ugc_wgw_manifest":
                members = list(v.get("sample_ids", []))  # type: ignore[arg-type]
                if stage == "assembly" and strings.get("trio"):
                    members = [str(v.get("father_id")), str(v.get("mother_id"))]
                _write_json(p, {"schema": 1, "ugc_wgw_version": ugc_wgw_version, "stage": stage,
                                "subject": {"type": "cohort" if "cohort_id" in v else "sample", "id": sid},
                                "cohort_members": members, "host": "fake", "written_at": "now"})
            else:
                write_fixture(p, name, 0, f"fake {name} for {sid}\n")
            outputs[ns + name] = str(p)
    for name, val in strings.items():
        outputs[ns + name] = val
    outputs[ns + "ugc_wgw_workflow_version"] = os.environ.get("UGC_WGW_FAKE_VERSION", ugc_wgw_version)
    return outputs


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] == "--version":
        print("miniwdl v0.0.0-fake")
        return 0
    ap = argparse.ArgumentParser()
    ap.add_argument("command")
    ap.add_argument("wdl")
    ap.add_argument("-i", "--input", required=True)
    ap.add_argument("-d", "--dir", required=True)
    ap.add_argument("-o", required=True)
    ap.add_argument("--cfg")
    ap.add_argument("--no-color", action="store_true")
    ap.add_argument("--log-json", action="store_true")
    a = ap.parse_args(argv)
    if a.command != "run":
        print(f"fake miniwdl: unsupported command {a.command}", file=sys.stderr)
        return 2
    stage = pathlib.Path(a.wdl).stem
    if not stage.startswith("ugc_wgw_"):
        print(f"fake miniwdl: unexpected wdl {a.wdl}", file=sys.stderr)
        return 2
    stage = stage[len("ugc_wgw_"):]
    run_dir = pathlib.Path(a.dir[:-2] if a.dir.endswith("/.") else a.dir).resolve()
    run_dir.mkdir(parents=True, exist_ok=True)
    inputs = json.loads(pathlib.Path(a.input).read_text())

    # input validation as miniwdl does it
    ns = f"ugc_wgw_{stage}."
    declared = declared_inputs(pathlib.Path(a.wdl))
    calls = declared_calls(pathlib.Path(a.wdl))
    bad = [k for k in inputs if not k.startswith(ns)
           or (k[len(ns):] not in declared and not ("." in k[len(ns):] and k[len(ns):].split(".", 1)[0] in calls))]
    missing = [n for n, req in declared.items() if req and ns + n not in inputs]
    if bad or missing:
        msg = f"check JSON input; unknown: {bad}; missing required: {missing}"
        print(msg, file=sys.stderr)
        _write_json(pathlib.Path(a.o), {"error": "InputError", "message": msg})
        return 2
    sid = str(inputs.get(ns + "sample_id") or inputs.get(ns + "cohort_id"))
    ugc_wgw_version = str(inputs.get(ns + "ugc_wgw_version"))

    (run_dir / "workflow.log").write_text(f"fake miniwdl run {stage} {sid}\n")
    (run_dir / "workflow.log.json").write_text(json.dumps({"level": "NOTICE", "message": f"fake miniwdl run {stage} {sid}"}) + "\n")
    (run_dir / "inputs.json").write_text(json.dumps(inputs, indent=2, sort_keys=True) + "\n")
    task_dir = run_dir / "call-fake"
    if os.environ.get("UGC_WGW_FAKE_SLURM_LOG") or os.environ.get("UGC_WGW_FAKE_STDERR"):
        task_dir.mkdir(exist_ok=True)
        if os.environ.get("UGC_WGW_FAKE_SLURM_LOG"):
            (task_dir / "slurm_singularity.log.txt").write_text(os.environ["UGC_WGW_FAKE_SLURM_LOG"].replace("\\n", "\n") + "\n")
        if os.environ.get("UGC_WGW_FAKE_STDERR"):
            (task_dir / "stderr.txt").write_text(os.environ["UGC_WGW_FAKE_STDERR"].replace("\\n", "\n") + "\n")
    _trace(f"start {sid}")

    def fail(cls: str, message: str, code: int) -> int:
        cause: dict[str, object] = {"error": cls, "exit_status": code, "run": "call-fake", "node": "call-fake",
                                    "dir": str(task_dir)}
        if cls == "CommandFailed":  # miniwdl 1.15: no message key, but the stream files
            cause["stderr_file"] = str(task_dir / "stderr.txt")
            cause["stdout_file"] = str(task_dir / "stdout.txt")
        else:
            cause["message"] = message
        err = {"error": "RunFailed", "workflow": f"ugc_wgw_{stage}", "dir": str(run_dir), "from_dir": str(task_dir),
               "cause": cause}
        _write_json(run_dir / "error.json", err)
        _write_json(pathlib.Path(a.o), err)
        _trace(f"end {sid}")
        return code if code else 2

    if os.environ.get("UGC_WGW_FAKE_HANG") == "1":
        if os.environ.get("UGC_WGW_FAKE_IGNORE_TERM") == "1":
            signal.signal(signal.SIGTERM, signal.SIG_IGN)
        else:
            def on_term(signum: int, _frame: object) -> None:
                raise SystemExit(fail("Interrupted", f"terminated by signal {signum}", 143))
            signal.signal(signal.SIGTERM, on_term)
        while True:
            time.sleep(0.1)

    time.sleep(float(os.environ.get("UGC_WGW_FAKE_SLEEP", "0")))
    failing = set(filter(None, os.environ.get("UGC_WGW_FAKE_FAIL", "").split(",")))
    attempt_no = run_dir.name.split("-")[-1]
    attempts = set(filter(None, os.environ.get("UGC_WGW_FAKE_FAIL_ATTEMPTS", "").split(",")))
    if (sid.endswith("_FAIL") or sid in failing) and (not attempts or attempt_no in attempts):
        spec = os.environ.get("UGC_WGW_FAKE_FAIL_CLASS", "CommandFailed:1")
        cls, _, code = spec.partition(":")
        return fail(cls, f"fake task failed for {sid}", int(code) if code else 1)

    outputs = make_outputs(stage, inputs, run_dir, ugc_wgw_version)
    _write_json(run_dir / "outputs.json", outputs)
    _write_json(pathlib.Path(a.o), {"dir": str(run_dir), "outputs": outputs})
    _trace(f"end {sid}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
