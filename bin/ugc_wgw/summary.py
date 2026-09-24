"""`ugc-wgw summary`: self-contained HTML analysis reports, one per sample and one per cohort (docs/guide/07-results.md,
"Analysis summary"; docs/DESIGN.md §8.1).

Data flow: the latest successful runs of a sample (`singleton`, or `upstream` + `downstream` in joint mode, plus
`assembly` when present) are read once into a JSON digest `<sample>.summary.json` (small numbers only: the big
files are streamed, never copied), the sample page is rendered from the digest, and the cohort page is rendered
from its members' digests, so a large cohort never re-reads sample files. Charts are inline SVG
(`summary_svg.py`), the page needs no network. QC flags are advisory driver-side checks with thresholds from
`config.json` (`summary_thresholds`) or `--threshold`; upstream's own messages are shown as they are.
"""
from __future__ import annotations

import concurrent.futures
import json
import os
from dataclasses import dataclass, field
from pathlib import Path

from . import summary_parsers as sp
from . import summary_svg as sv
from .config import Config
from .db import DB, Cohort, RunRecord
from .inputs import SEX_TO_GROUP
from .layout import run_files
from .log import LOGGER
from .manifest import CodeInfo
from .outputs import Outputs
from .report import CSS as REPORT_CSS
from .report import JS as REPORT_JS
from .report import Raw, _bytes, _esc, _table
from .stages import namespace
from .util import UgcError, diag, read_json, utc_now, write_json

DIGEST_SCHEMA = 1
# The public repository (scripts/export-public.sh fills the owner); `project_url` in config.json overrides it.
DEFAULT_PROJECT_URL = "https://github.com/UppsalaGenomeCenter/ugc-wgw"

DEFAULT_THRESHOLDS: dict[str, float] = {
    "depth_mean_min": 20.0,           # mean depth (x) below which the sample fails
    "mapped_read_percent_min": 95.0,  # mapped reads (%) below which the sample fails
    "read_quality_median_min": 25.0,  # median read quality (Phred) below which the sample fails
}

STAGE_ORDER = ("singleton", "downstream", "upstream")  # which stage's file wins when several carry it
SOURCE_MAP: dict[str, dict[str, str]] = {  # logical name -> {stage: output name}
    "stats_file": {"singleton": "stats_file", "downstream": "stats_file"},
    "msg_file": {"singleton": "msg_file", "downstream": "msg_file"},
    "mosdepth_summary": {"singleton": "mosdepth_summary", "upstream": "mosdepth_summary"},
    "mosdepth_region_bed": {"singleton": "mosdepth_region_bed", "upstream": "mosdepth_region_bed"},
    "sv_copynum_summary": {"singleton": "sv_copynum_summary", "upstream": "sv_copynum_summary"},
    "sv_copynum_bedgraph": {"singleton": "sv_copynum_bedgraph", "upstream": "sv_copynum_bedgraph"},
    "small_variant_stats": {"singleton": "small_variant_stats", "downstream": "small_variant_stats"},
    "bcftools_roh_out": {"singleton": "bcftools_roh_out", "downstream": "bcftools_roh_out"},
    "sv_vcf": {"singleton": "phased_sv_vcf", "downstream": "phased_sv_vcf", "upstream": "sv_vcf"},
    "phase_stats": {"singleton": "phase_stats", "downstream": "phase_stats"},
    "phase_blocks": {"singleton": "phase_blocks", "downstream": "phase_blocks"},
    "trgt_vcf": {"singleton": "phased_trgt_vcf", "downstream": "trgt_vcf"},
    "trgt_coverage_dropouts": {"singleton": "trgt_coverage_dropouts", "downstream": "trgt_coverage_dropouts"},
    "methbat_profile": {"singleton": "methbat_profile", "downstream": "methbat_profile"},
    "paraphase_json": {"singleton": "paraphase_summary", "upstream": "paraphase_output_json"},
    "mitorsaw_hap_stats": {"singleton": "mitorsaw_hap_stats", "upstream": "mitorsaw_hap_stats"},
    "pbstarphase_tsv": {"singleton": "pbstarphase_tsv", "downstream": "pbstarphase_tsv"},
    "pbstarphase_json": {"singleton": "pbstarphase_summary", "downstream": "pbstarphase_json"},
    "kivvi_kiv2_json": {"singleton": "kivvi_kiv2_json", "upstream": "kivvi_kiv2_json"},
    "kivvi_d4z4_json": {"singleton": "kivvi_d4z4_json", "upstream": "kivvi_d4z4_json"},
    "haplotagged_bam": {"singleton": "merged_haplotagged_bam", "downstream": "merged_haplotagged_bam"},
    "aligned_bam": {"upstream": "aligned_hifi_reads"},
    "small_variant_vcf": {"singleton": "phased_small_variant_vcf", "downstream": "phased_small_variant_vcf", "upstream": "small_variant_vcf"},
    "small_variant_gvcf": {"singleton": "small_variant_gvcf", "upstream": "small_variant_gvcf"},
    "cpg_pileup_bed": {"singleton": "cpg_pileup_bed", "downstream": "cpg_pileup_bed"},
}
STAT_ALIASES = {"small_variant_SNV_count": "SNV_count", "small_variant_INDEL_count": "INDEL_count",
                "small_variant_TSTV_ratio": "TSTV_ratio", "small_variant_HETHOM_ratio": "HETHOM_ratio",
                "phased_basepairs": "stat_phased_basepairs"}
READ_KEYS = ["read_count", "read_length_mean", "read_length_median", "read_length_n50", "read_quality_mean",
             "read_quality_median", "mapped_read_count", "mapped_read_percent", "gap_compressed_identity_mean",
             "gap_compressed_identity_median"]
CITATIONS: list[tuple[str, str]] = [
    ("HiFi-human-WGS-WDL (PacBio)", "https://github.com/PacificBiosciences/HiFi-human-WGS-WDL"),
    ("pbmm2 / minimap2", "https://github.com/PacificBiosciences/pbmm2; Li H. Bioinformatics 2018, doi:10.1093/bioinformatics/bty191"),
    ("DeepVariant", "https://github.com/google/deepvariant; Poplin R et al. Nat Biotechnol 2018, doi:10.1038/nbt.4235"),
    ("sawfish", "https://github.com/PacificBiosciences/sawfish"),
    ("TRGT", "https://github.com/PacificBiosciences/trgt; Dolzhenko E et al. Nat Biotechnol 2024, doi:10.1038/s41587-023-02057-3"),
    ("HiPhase", "https://github.com/PacificBiosciences/HiPhase; Holt JM et al. Bioinformatics 2024, doi:10.1093/bioinformatics/btae042"),
    ("mosdepth", "https://github.com/brentp/mosdepth; Pedersen BS, Quinlan AR. Bioinformatics 2018, doi:10.1093/bioinformatics/btx699"),
    ("bcftools / samtools", "https://github.com/samtools/bcftools; Danecek P et al. GigaScience 2021, doi:10.1093/gigascience/giab008"),
    ("Paraphase", "https://github.com/PacificBiosciences/paraphase"),
    ("mitorsaw", "https://github.com/PacificBiosciences/mitorsaw"),
    ("MethBat", "https://github.com/PacificBiosciences/MethBat"),
    ("StarPhase", "https://github.com/PacificBiosciences/pb-StarPhase"),
    ("kivvi", "https://github.com/PacificBiosciences/kivvi"),
    ("pbjam", "https://github.com/PacificBiosciences/pbjam"),
    ("GLnexus", "https://github.com/dnanexus-rnd/GLnexus; Yun T et al. Bioinformatics 2021, doi:10.1093/bioinformatics/btaa1081"),
    ("svx", "https://github.com/PacificBiosciences/svx"),
    ("trgt-lps", "https://github.com/PacificBiosciences/trgt-lps"),
    ("hifiasm", "https://github.com/chhylp123/hifiasm; Cheng H et al. Nat Methods 2021, doi:10.1038/s41592-020-01056-5"),
    ("yak, gfatools", "https://github.com/lh3/yak; https://github.com/lh3/gfatools"),
]


# ---- sources ------------------------------------------------------------------------------


@dataclass
class Sources:
    """Everything the digest builder needs, as plain values (picklable for --jobs)."""
    sample_id: str
    mode: str
    runs: dict[str, RunRecord]                 # stage -> run
    docs: dict[str, dict[str, object]]         # stage -> outputs.json
    sheet_sex: str | None = None
    father_id: str | None = None
    mother_id: str | None = None

    def _stages(self) -> list[str]:
        return [s for s in STAGE_ORDER if s in self.docs]

    def output(self, stage: str, name: str) -> object | None:
        return self.docs.get(stage, {}).get(f"{namespace(stage)}.{name}")

    def path(self, logical: str) -> Path | None:
        for stage in self._stages():
            name = SOURCE_MAP.get(logical, {}).get(stage)
            if name:
                v = self.output(stage, name)
                if isinstance(v, str) and v and Path(v).exists():
                    return Path(v)
        return None

    def scalar(self, name: str) -> str | None:
        for stage in self._stages():
            v = self.output(stage, name)
            if v is not None and not isinstance(v, (list, dict)):
                return str(v)
        return None

    def stats(self) -> dict[str, str]:
        """The stats table: the stats_file columns, then the `stat_*` scalar outputs (aliases applied)."""
        out: dict[str, str] = {}
        p = self.path("stats_file")
        if p:
            out.update(sp.parse_two_line_tsv(sp.read_lines(p)))
        for stage in reversed(self._stages()):
            ns = namespace(stage) + ".stat_"
            for k, v in self.docs[stage].items():
                if k.startswith(ns) and v is not None and not isinstance(v, (list, dict)):
                    name = k[len(ns):]
                    name = STAT_ALIASES.get(name, name)
                    out.setdefault(name, str(v))
        sex = self.scalar("inferred_sex")
        if sex is not None:
            out["inferred_sex"] = sex
        return out

    def messages(self) -> list[str]:
        seen: list[str] = []
        p = self.path("msg_file")
        if p:
            for line in sp.read_lines(p):
                if line.strip() and line not in seen:
                    seen.append(line)
        for stage in ("upstream", "downstream", "singleton"):
            v = self.output(stage, "msg")
            if isinstance(v, list):
                for m in v:
                    if str(m).strip() and str(m) not in seen:
                        seen.append(str(m))
        return seen

    def manifest(self) -> dict[str, object]:
        for stage in self._stages():
            run = self.runs.get(stage)
            if run:
                p = run_files(run.run_path).manifest
                if p.exists():
                    doc = read_json(p)
                    if isinstance(doc, dict):
                        return doc
        return {}


def latest_success_in_mode(db: DB, subject_type: str, subject_id: str, stage: str, mode: str, version: str | None) -> RunRecord | None:
    """The latest successful run of a stage that was submitted in `mode` (the stage names cohort_merge and
    cohort_freq exist in both modes and mean different inputs)."""
    runs = [r for r in db.runs_for(subject_type, subject_id, stage)
            if r.status == "success" and r.mode == mode and (version is None or r.ugc_wgw_version == version)]
    if not runs:
        return None
    return max(runs, key=lambda r: (r.finished_at or "", r.attempt))


def resolve_sources(db: DB, sample_id: str, mode: str, version: str | None) -> Sources | None:
    """The latest successful runs of a sample in `mode` (plus assembly); None when it has none."""
    stages = ("singleton",) if mode == "standalone" else ("upstream", "downstream")
    runs: dict[str, RunRecord] = {}
    docs: dict[str, dict[str, object]] = {}
    for stage in stages + ("assembly",):
        run = latest_success_in_mode(db, "sample", sample_id, stage, "assembly" if stage == "assembly" else mode, version)
        if run is None:
            continue
        try:
            docs[stage] = Outputs.load(run).doc
        except UgcError as exc:
            LOGGER.warning("summary: %s %s: %s", sample_id, stage, exc)
            continue
        runs[stage] = run
    if not any(s in runs for s in stages):
        return None
    rec = db.get_sample(sample_id)
    return Sources(sample_id=sample_id, mode=mode, runs=runs, docs=docs,
                   sheet_sex=(rec.sex or None) if rec else None,
                   father_id=(rec.father_id or None) if rec else None,
                   mother_id=(rec.mother_id or None) if rec else None)


def resolve_cohort_runs(db: DB, cohort_id: str, mode: str, version: str | None) -> dict[str, RunRecord]:
    stages = ("cohort_call", "cohort_merge", "cohort_freq") if mode == "joint" else ("cohort_merge", "cohort_freq")
    out: dict[str, RunRecord] = {}
    for stage in stages:
        run = latest_success_in_mode(db, "cohort", cohort_id, stage, mode, version)
        if run is not None:
            out[stage] = run
    return out


# ---- thresholds and QC flags ------------------------------------------------------------------


def thresholds_for(cfg: Config, overrides: list[str] | None = None) -> dict[str, float]:
    """Defaults, then config.json `summary_thresholds`, then `--threshold KEY=VALUE`; unknown keys are errors."""
    thr = dict(DEFAULT_THRESHOLDS)
    for source, items in (("config.json summary_thresholds", cfg.summary_thresholds.items()),
                          ("--threshold", [tuple(o.split("=", 1)) if "=" in o else (o, "") for o in (overrides or [])])):
        for key, value in items:
            if key not in DEFAULT_THRESHOLDS:
                raise UgcError(f"{source}: unknown threshold {key!r} (known: {', '.join(sorted(DEFAULT_THRESHOLDS))})")
            try:
                thr[key] = float(value)  # type: ignore[arg-type]
            except (TypeError, ValueError):
                raise UgcError(f"{source}: {key} needs a number, got {value!r}") from None
    return thr


def _stat_num(stats: dict[str, object], key: str) -> float | None:
    v = stats.get(key)
    if v in (None, "", "NA"):
        return None
    try:
        return float(str(v))
    except ValueError:
        return None


def qc_flags(digest: dict[str, object], thr: dict[str, float]) -> dict[str, object]:
    """Advisory flags; status is the worst severity (pass < info < warn < fail)."""
    stats: dict[str, object] = digest.get("stats", {})  # type: ignore[assignment]
    flags: list[dict[str, object]] = []

    def flag(code: str, severity: str, message: str, value: object = None, threshold: object = None) -> None:
        flags.append({"code": code, "severity": severity, "message": message, "value": value, "threshold": threshold})

    for code, key, tkey, label, unit in (("low_depth", "depth_mean", "depth_mean_min", "mean depth", "x"),
                                         ("low_mapped", "mapped_read_percent", "mapped_read_percent_min", "mapped reads", "%"),
                                         ("low_read_quality", "read_quality_median", "read_quality_median_min", "median read quality", "")):
        v = _stat_num(stats, key)
        if v is not None and v < thr[tkey]:
            flag(code, "fail", f"{label} {sv.fmt(v, 1)}{unit} is below {sv.fmt(thr[tkey], 1)}{unit}", v, thr[tkey])
    sample: dict[str, object] = digest.get("sample", {})  # type: ignore[assignment]
    sheet = SEX_TO_GROUP.get(str(sample.get("sex_sheet") or "").upper())
    inferred_raw = str(stats.get("inferred_sex") or "")
    inferred = SEX_TO_GROUP.get(inferred_raw.upper())
    if sheet and inferred and sheet != inferred:
        flag("sex_mismatch", "warn", f"sex in the sheet is {sample.get('sex_sheet')} ({sheet}) but coverage infers "
                                    f"{inferred_raw} ({inferred})", inferred_raw, sample.get("sex_sheet"))
    elif not inferred_raw and digest.get("coverage"):
        flag("sex_not_inferred", "warn", "no sex inferred from coverage (autosome depth is zero or mosdepth did not run)")
    targeted: dict[str, object] = digest.get("targeted", {})  # type: ignore[assignment]
    kivvi: dict[str, object] = targeted.get("kivvi", {})  # type: ignore[assignment]
    for region, d in kivvi.items():
        status = str((d or {}).get("call_status") or "")
        if "fail" in status.lower():
            flag("tool_status", "warn", f"kivvi {region}: {status}", status)
    para: dict[str, object] = targeted.get("paraphase", {})  # type: ignore[assignment]
    if para and para.get("n_failed"):
        flag("tool_status", "warn", f"Paraphase: {para['n_failed']} of {para.get('n_regions')} regions failed for coverage",
             para["n_failed"], para.get("n_regions"))
    mito: dict[str, object] = targeted.get("mitorsaw", {})  # type: ignore[assignment]
    if mito and not mito.get("haplotypes"):
        flag("tool_status", "warn", "mitorsaw reported no haplotypes")
    if digest.get("messages"):
        flag("messages", "info", f"{len(digest['messages'])} workflow message(s)")  # type: ignore[arg-type]
    if digest.get("warnings"):
        flag("parser_warnings", "info", f"{len(digest['warnings'])} file(s) could not be summarised")  # type: ignore[arg-type]
    rank = {"pass": 0, "info": 1, "warn": 2, "fail": 3}
    status = max((str(f["severity"]) for f in flags), key=lambda s: rank.get(s, 0), default="pass")
    return {"status": status, "flags": flags, "thresholds": dict(thr)}


# ---- the digest ------------------------------------------------------------------------------


def _bin_series(pairs: list[list[float]], edges: list[float], labels: list[str], value_index: int = 1) -> list[list[object]]:
    out = {label: 0.0 for label in labels}
    for row in pairs:
        x, v = float(row[0]), float(row[value_index])
        out[sp._bin_label(x, edges, labels)] += v
    return [[label, out[label]] for label in labels]


def _dp_bins(dp: list[list[object]]) -> list[list[object]]:
    labels = ["<10", "10-20", "20-30", "30-40", "40-60", "60-100", ">=100"]
    edges = [10, 20, 30, 40, 60, 100]
    out = {label: 0 for label in labels}
    for bin_label, genotypes, _sites in dp:
        s = str(bin_label)
        x = float(s.lstrip("<>")) if s.lstrip("<>").replace(".", "", 1).isdigit() else None
        if x is None:
            continue
        if s.startswith(">"):
            x = x + 1
        out[sp._bin_label(x, edges, labels)] += int(genotypes)
    return [[label, out[label]] for label in labels]


def build_digest(src: Sources, thr: dict[str, float], y_threshold: float = 0.1) -> dict[str, object]:
    """Read the sample's result files once into plain numbers. A file that cannot be read becomes a warning."""
    d: dict[str, object] = {"schema": DIGEST_SCHEMA, "generated_at": utc_now(), "mode": src.mode,
                            "sample": {"id": src.sample_id, "sex_sheet": src.sheet_sex, "father_id": src.father_id,
                                       "mother_id": src.mother_id},
                            "runs": {s: {"run_id": r.run_id, "ugc_wgw_version": r.ugc_wgw_version, "attempt": r.attempt,
                                         "finished_at": r.finished_at, "run_dir": r.run_dir} for s, r in src.runs.items()},
                            "warnings": []}
    warnings: list[str] = d["warnings"]  # type: ignore[assignment]

    def safe(name: str, fn, *args):  # type: ignore[no-untyped-def]
        try:
            return fn(*args)
        except Exception as exc:  # noqa: BLE001 - one bad file must not sink the report
            warnings.append(f"{name}: {type(exc).__name__}: {exc}")
            return None

    def text(logical: str, parser, *args):  # type: ignore[no-untyped-def]
        p = src.path(logical)
        return safe(logical, lambda: parser(sp.read_lines(p), *args)) if p else None

    def binary(logical: str, parser, *args):  # type: ignore[no-untyped-def]
        p = src.path(logical)
        return safe(logical, lambda: parser(sp.read_binary_lines(p), *args)) if p else None

    def jsonf(logical: str, parser):  # type: ignore[no-untyped-def]
        p = src.path(logical)
        return safe(logical, lambda: parser(sp.load_json(p))) if p else None

    stats = safe("stats", src.stats) or {}
    d["stats"] = stats
    d["reads"] = {k: _stat_num(stats, k) for k in READ_KEYS}

    # coverage
    md = text("mosdepth_summary", sp.parse_mosdepth_summary) or {}
    chrom_rows = [r for r in md.get("chrom", []) if r["chrom"] in sp.PRIMARY_CONTIGS or r["chrom"] == "chrM"]
    ratios = sp.coverage_ratios(md.get("chrom", []), y_threshold) if md else {}
    track = binary("mosdepth_region_bed", sp.aggregate_region_bed) or {}
    cn = jsonf("sv_copynum_summary", sp.parse_copynum_summary) or {}
    segs = text("sv_copynum_bedgraph", sp.parse_bedgraph) or {}
    d["coverage"] = {"mean": _stat_num(stats, "depth_mean"), "inferred_sex": stats.get("inferred_sex") or "",
                     "total": md.get("total", {}), "chrom": chrom_rows, "ratios": ratios, "track": track,
                     "copynum": {"haploid_coverage": cn.get("haploid_coverage"),
                                 "chrom": [r for r in cn.get("chrom", []) if r["chrom"] in sp.PRIMARY_CONTIGS],
                                 "segments": segs}}

    # small variants
    bs = text("small_variant_stats", sp.parse_bcftools_stats) or {}
    roh = binary("bcftools_roh_out", sp.parse_roh_out) or {}
    d["small_variants"] = {
        "snv": _stat_num(stats, "SNV_count"), "indel": _stat_num(stats, "INDEL_count"),
        "tstv": _stat_num(stats, "TSTV_ratio"), "hethom": _stat_num(stats, "HETHOM_ratio"),
        "sn": bs.get("sn", {}), "psc": bs.get("psc", {}), "st": bs.get("st", {}),
        "idd": _bin_series(bs.get("idd", []), [-49, -9, 1, 10, 50], ["<=-50", "-49..-10", "-9..-1", "+1..+9", "+10..+49", ">=+50"]) if bs.get("idd") else [],
        "qual": _bin_series([[q, s + i] for q, s, _t, _v, i in bs.get("qual", [])], [10, 20, 30, 40, 50, 60],
                            ["<10", "10-20", "20-30", "30-40", "40-50", "50-60", ">=60"]) if bs.get("qual") else [],
        "dp": _dp_bins(bs.get("dp", [])) if bs.get("dp") else [],
        "roh": roh}

    # structural variants
    svd = binary("sv_vcf", sp.parse_sv_vcf) or {}
    d["sv"] = {"counts": {t: _stat_num(stats, f"sv_{t}_count") for t in ("DEL", "INS", "DUP", "INV", "BND", "SWAP")},
               "vcf_counts": svd.get("counts", {}), "sizes": svd.get("sizes", {}), "size_bins": svd.get("size_bins", sp.SV_SIZE_BINS),
               "filters": svd.get("filters", {}), "n": svd.get("n")}

    # phasing
    ps = text("phase_stats", sp.parse_phase_stats) or {}
    pb = text("phase_blocks", sp.parse_phase_blocks) or {}
    keep = ("chromosome", "num_variants", "num_heterozygous", "num_phased", "num_unphased", "num_blocks",
            "basepairs_per_block_sum", "block_ng50")
    d["phasing"] = {"phased_bp": _stat_num(stats, "stat_phased_basepairs"), "block_ng50": _stat_num(stats, "phase_block_ng50"),
                    "all": {k: ps.get("all", {}).get(k) for k in keep} if ps.get("all") else {},
                    "chrom": [{k: r.get(k) for k in keep} for r in ps.get("chrom", []) if str(r.get("chromosome")) in sp.PRIMARY_CONTIGS],
                    "blocks": pb}

    # tandem repeats
    messages = safe("messages", src.messages) or []
    fr = sp.fail_reads_trids(messages)
    tv = binary("trgt_vcf", sp.parse_trgt_vcf, fr) or {}
    dr = binary("trgt_coverage_dropouts", sp.parse_trgt_dropouts) or {}
    loci = tv.get("loci", [])
    for locus in loci:
        locus["dropout"] = dr.get("named", {}).get(locus["trid"], "")
    d["tandem_repeats"] = {"genotyped": _stat_num(stats, "trgt_genotyped_count"), "uncalled": _stat_num(stats, "trgt_uncalled_count"),
                           "n_records": tv.get("n_records"), "dropouts": {"n": dr.get("n"), "by_class": dr.get("by_class", {})},
                           "n_named": tv.get("n_named"), "n_named_called": tv.get("n_named_called"),
                           "disease_loci": loci, "fail_reads_trids": fr}

    # methylation
    mb = text("methbat_profile", sp.parse_methbat_profile) or {}
    d["methylation"] = {"cpg": {k: _stat_num(stats, f"cpg_{k}_count") for k in ("combined", "hap1", "hap2")},
                        "methbat": {k: _stat_num(stats, f"methbat_{k}_count") for k in ("methylated", "unmethylated", "asm")},
                        "profile": mb}

    # targeted callers
    d["targeted"] = {"paraphase": jsonf("paraphase_json", sp.parse_paraphase) or {},
                     "mitorsaw": jsonf("mitorsaw_hap_stats", sp.parse_mitorsaw) or {},
                     "starphase": text("pbstarphase_tsv", sp.parse_starphase_tsv) or [],
                     "starphase_meta": jsonf("pbstarphase_json", sp.parse_starphase_json) or {},
                     "kivvi": {k: (jsonf(f"kivvi_{k}_json", sp.parse_kivvi) or {}) for k in ("kiv2", "d4z4")}}

    d["messages"] = messages
    files: dict[str, object] = {}
    for logical in SOURCE_MAP:
        p = src.path(logical)
        if p:
            try:
                files[logical] = {"path": str(p), "bytes": p.stat().st_size}
            except OSError:
                files[logical] = {"path": str(p), "bytes": None}
    d["files"] = files

    # assembly (when the sample has one)
    asm: dict[str, object] = {}
    if "assembly" in src.docs:
        doc = src.docs["assembly"]
        ns = namespace("assembly") + "."
        stats_paths = doc.get(ns + "assembly_stats") if isinstance(doc.get(ns + "assembly_stats"), list) else []
        paf_paths = doc.get(ns + "paftools_vcf_stats") if isinstance(doc.get(ns + "paftools_vcf_stats"), list) else []
        haps = doc.get(ns + "haplotypes") if isinstance(doc.get(ns + "haplotypes"), list) else []
        asm = {"trio": bool(doc.get(ns + "trio")), "haplotypes": []}
        for i, p in enumerate(stats_paths):  # type: ignore[arg-type]
            g = safe("assembly_stats", lambda p=p: sp.parse_gfatools_stats(sp.read_lines(Path(str(p))))) or {}
            paf = {}
            if i < len(paf_paths):  # type: ignore[arg-type]
                paf = (safe("paftools_vcf_stats", lambda q=paf_paths[i]: sp.parse_bcftools_stats(sp.read_lines(Path(str(q))))) or {}).get("sn", {})  # type: ignore[index]
            asm["haplotypes"].append({"name": str(haps[i]) if i < len(haps) else f"hap{i + 1}", "stats": g,  # type: ignore[attr-defined, arg-type]
                                      "paftools": {"snps": paf.get("number of SNPs"), "indels": paf.get("number of indels"),
                                                   "records": paf.get("number of records")}})
    d["assembly"] = asm

    # provenance from the primary run's manifest
    man = safe("run_manifest", src.manifest) or {}
    upstream = man.get("upstream") if isinstance(man.get("upstream"), dict) else {}
    references = man.get("references") if isinstance(man.get("references"), dict) else {}
    ugc_wgw_man = man.get("ugc_wgw_manifest") if isinstance(man.get("ugc_wgw_manifest"), dict) else {}
    d["provenance"] = {
        "ugc_wgw": man.get("ugc_pacbio_wgw", {}), "engine": man.get("engine", {}), "host": man.get("host"),
        "finished_at": man.get("finished_at"),
        "upstream": {k: {"tag": v.get("tag"), "commit": v.get("commit")} for k, v in upstream.items() if isinstance(v, dict)},
        "workflow": ugc_wgw_man.get("upstream", {}) if isinstance(ugc_wgw_man, dict) else {},
        "references": {k: {"version": v.get("version"), "build": v.get("build"), "source": v.get("source"),
                           "files": v.get("files", [])} for k, v in references.items() if isinstance(v, dict)},
        "containers": sorted({str(c.get("image")) for c in man.get("containers", []) if isinstance(c, dict)}),
    }
    d["qc"] = qc_flags(d, thr)
    return d


def digest_is_fresh(doc: object, src: Sources) -> bool:
    if not isinstance(doc, dict) or doc.get("schema") != DIGEST_SCHEMA:
        return False
    runs = doc.get("runs")
    if not isinstance(runs, dict) or set(runs) != set(src.runs):
        return False
    return all(isinstance(runs.get(s), dict) and runs[s].get("run_id") == r.run_id for s, r in src.runs.items())


# ---- the pages ---------------------------------------------------------------------------------

SUMMARY_CSS = """
:root{--m1:#2f6f8f;--m2:#c98a2b;--m3:#5b8c5a;--m4:#a24b6f;--m5:#7b6bb0;--m6:#8a8f95;--ok:#2e7d32;--warn:#ef6c00;--bad:#c62828;--lane:#e9edf1;--iqr:#d6dbe0}
@media (prefers-color-scheme: dark){:root:not([data-theme="light"]){--m1:#6fb7c7;--m2:#e0b060;--m3:#8fc48d;--m4:#d07a9c;--m5:#a596d8;--m6:#9aa5b1;--ok:#66bb6a;--warn:#ffa040;--bad:#ef5350;--lane:#232a31;--iqr:#333b44}}
:root[data-theme="dark"]{--m1:#6fb7c7;--m2:#e0b060;--m3:#8fc48d;--m4:#d07a9c;--m5:#a596d8;--m6:#9aa5b1;--ok:#66bb6a;--warn:#ffa040;--bad:#ef5350;--lane:#232a31;--iqr:#333b44}
body{max-width:none;padding:0}
nav.side{position:fixed;left:0;top:0;bottom:0;width:200px;overflow:auto;border-right:1px solid var(--line);background:var(--panel);padding:16px 12px;box-sizing:border-box;font-size:12px}
nav.side a{display:block;padding:4px 6px;color:var(--fg);text-decoration:none;border-radius:3px}nav.side a:hover{background:var(--head)}
nav.side .brand{font-weight:600;margin-bottom:10px;color:var(--accent);word-break:break-all}nav.side .sub{color:var(--muted);margin:10px 6px 4px;text-transform:uppercase;letter-spacing:0.06em;font-size:10px}
main{margin-left:200px;padding:24px 28px 48px;max-width:1240px;box-sizing:border-box}
section{scroll-margin-top:12px}h2{margin-top:40px}h3{font-size:13px;margin:18px 0 6px}h4{font-size:12px;margin:12px 0 4px;color:var(--muted)}
.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(380px,1fr));gap:12px 24px;align-items:start}
.m1{fill:var(--m1)}.m2{fill:var(--m2)}.m3{fill:var(--m3)}.m4{fill:var(--m4)}.m5{fill:var(--m5)}.m6{fill:var(--m6)}
.ln1{stroke:var(--m1);fill:none;stroke-width:1.5}.ln2{stroke:var(--m2);fill:none;stroke-width:2}.ln3{stroke:var(--m4);fill:none;stroke-width:1.5}
.ok{fill:var(--ok)}.warn{fill:var(--warn)}.bad{fill:var(--bad)}.pt{fill:var(--m1);opacity:.8}
.lane{fill:var(--lane)}.iqr{fill:var(--iqr)}.base{stroke:var(--warn)}.thr{stroke:var(--bad)}.axis{stroke:var(--grid)}
svg.chart{max-width:100%;height:auto}svg.legend{display:block;margin:4px 0}
.flag-fail{color:var(--bad);font-weight:600}.flag-warn{color:var(--warn);font-weight:600}.flag-pass{color:var(--ok);font-weight:600}.flag-info{color:var(--muted)}
.status{display:inline-block;padding:2px 8px;border-radius:3px;font-weight:600;color:#fff}.status-pass{background:var(--ok)}.status-info{background:var(--m6)}.status-warn{background:var(--warn)}.status-fail{background:var(--bad)}
ul.plain{padding-left:18px}ul.plain li{margin:2px 0}.small{font-size:11px}
@media print{nav.side{display:none}main{margin-left:0}}
"""


def _sec(sec_id: str, title: str, body: str) -> str:
    return f'<section id="{sec_id}"><h2>{_esc(title)}</h2>{body}</section>'


def _cards(items: list[tuple[str, object]]) -> str:
    return '<div class="cards">' + "".join(
        f'<div class="card"><span class="muted">{_esc(k)}</span><b>{v if isinstance(v, Raw) else _esc(v)}</b></div>' for k, v in items) + "</div>"


def _kv(rows: list[list[object]]) -> str:
    return _table(["what", "value"], rows, cls="")


def _nav(brand: str, items: list[tuple[str, str]], extra: str = "") -> str:
    links = "".join(f'<a href="#{sec_id}">{_esc(title)}</a>' for sec_id, title in items)
    return f'<nav class="side"><div class="brand">{_esc(brand)}</div>{links}{extra}</nav>'


def _page(title: str, nav: str, body: str) -> str:
    return (f"<!DOCTYPE html><html lang='en'><head><meta charset='utf-8'><title>{_esc(title)}</title>"
            f"<style>{REPORT_CSS}{SUMMARY_CSS}</style></head><body>{nav}<main>{body}</main>"
            f"<script>{REPORT_JS}</script></body></html>")


def _clip(text: str, n: int) -> Raw:
    """A cell that shows at most n characters, the whole text as a tooltip."""
    if len(text) <= n:
        return Raw(_esc(text))
    return Raw(f'<span title="{_esc(text)}">{_esc(text[: n - 1])}…</span>')


def _status_badge(status: str) -> Raw:
    return Raw(f'<span class="status status-{_esc(status)}">{_esc(status)}</span>')


def _flags_table(qc: dict[str, object]) -> str:
    flags: list[dict[str, object]] = qc.get("flags", [])  # type: ignore[assignment]
    if not flags:
        return "<p>No flags: every check passed.</p>"
    rows = [[Raw(f'<span class="flag-{_esc(f["severity"])}">{_esc(f["severity"])}</span>'), f["code"], f["message"],
             sv.fmt(f.get("value"), 2) if isinstance(f.get("value"), (int, float)) else (f.get("value") or ""),
             sv.fmt(f.get("threshold"), 2) if isinstance(f.get("threshold"), (int, float)) else (f.get("threshold") or "")] for f in flags]
    return _table(["severity", "flag", "message", "value", "threshold"], rows, cls="")


def _software_sections(prov: dict[str, object], tools: list[dict[str, object]], project_url: str, versions_note: str = "") -> str:
    containers = set(prov.get("containers", []))  # type: ignore[arg-type]
    tool_rows = []
    for t in tools:
        digest = str(t.get("digest") or "")
        tag = f"{t['tool']}:{t['version']}"
        in_run = any(digest in c for c in containers) if digest else (any(tag in c for c in containers) or None)
        tool_rows.append([t["tool"], t["version"], t.get("build") or "", Raw(f"<code>{_esc(digest[:12])}</code>") if digest else "",
                          "yes" if in_run else ("–" if in_run is None else "no"), Raw("<br>".join(_esc(f) for f in t.get("files", [])))])  # type: ignore[union-attr]
    engine: dict[str, object] = prov.get("engine", {})  # type: ignore[assignment]
    ugc: dict[str, object] = prov.get("ugc_wgw", {})  # type: ignore[assignment]
    wf = prov.get("workflow", {})
    wf_rows: list[list[object]] = [["ugc-pacbio-wgw", f"{ugc.get('version', '')} ({ugc.get('git_commit', '')})"]]
    for name, u in sorted(prov.get("upstream", {}).items()):  # type: ignore[union-attr]
        wf_rows.append([f"upstream {name}", f"{u.get('tag', '')} ({str(u.get('commit', ''))[:12]})"])
    if isinstance(wf, dict) and wf:
        wf_rows.append(["upstream workflow", f"{wf.get('workflow_name', '')} {wf.get('workflow_version', '')}"])
    for name, r in sorted(prov.get("references", {}).items()):  # type: ignore[union-attr]
        wf_rows.append([f"reference bundle {name}", f"{r.get('version', '')} {r.get('build') or ''} · {r.get('source') or ''}"])
    for k, v in sorted(engine.items()):
        wf_rows.append([f"engine {k}", str(v)])
    wf_rows.append(["run host", str(prov.get("host") or "")])
    wf_rows.append(["project", Raw(f'<a href="{_esc(project_url)}">{_esc(project_url)}</a>')])
    ref_files = []
    for name, r in sorted(prov.get("references", {}).items()):  # type: ignore[union-attr]
        for f in r.get("files", []) or []:
            if isinstance(f, dict):
                ref_files.append([name, f.get("name", ""), sv.fmt(f.get("bytes")), Raw(f"<code>{_esc(f.get('md5') or f.get('sha256') or '')}</code>")])
    body = ("<h3>Tools</h3><p class='small muted'>Versions as pinned in the WDL task definitions of this code version; "
            "the digest column is the container image; 'in run' says whether that image is in the run's manifest." + _esc(versions_note) + "</p>"
            + '<div class="wrap">' + _table(["tool", "version", "build", "image digest", "in run", "defined in"], tool_rows) + "</div>"
            + "<h3>Workflow and references</h3>" + _kv(wf_rows)
            + ("<details><summary>reference files</summary>" + _table(["bundle", "file", "bytes", "checksum"], ref_files, cls="") + "</details>" if ref_files else "")
            + "<h3>Citations</h3><ul class='plain small'>" + "".join(f"<li><b>{_esc(n)}</b>: {_esc(c)}</li>" for n, c in CITATIONS) + "</ul>"
            + "<p class='small muted'>The workflow itself: HiFi-human-WGS-WDL and this composition layer (above). Please cite the tools whose results you use.</p>")
    return body


def render_sample(d: dict[str, object], tools: list[dict[str, object]], *, project_name: str, project_url: str,
                  cohort_links: list[tuple[str, str]] | None = None) -> str:
    s: dict[str, object] = d["sample"]  # type: ignore[assignment]
    sid = str(s["id"])
    stats: dict[str, object] = d.get("stats", {})  # type: ignore[assignment]
    qc: dict[str, object] = d.get("qc", {})  # type: ignore[assignment]
    reads: dict[str, object] = d.get("reads", {})  # type: ignore[assignment]
    cov: dict[str, object] = d.get("coverage", {})  # type: ignore[assignment]
    smv: dict[str, object] = d.get("small_variants", {})  # type: ignore[assignment]
    svd: dict[str, object] = d.get("sv", {})  # type: ignore[assignment]
    ph: dict[str, object] = d.get("phasing", {})  # type: ignore[assignment]
    tr: dict[str, object] = d.get("tandem_repeats", {})  # type: ignore[assignment]
    me: dict[str, object] = d.get("methylation", {})  # type: ignore[assignment]
    tg: dict[str, object] = d.get("targeted", {})  # type: ignore[assignment]
    prov: dict[str, object] = d.get("provenance", {})  # type: ignore[assignment]
    runs: dict[str, dict[str, object]] = d.get("runs", {})  # type: ignore[assignment]
    inferred = str(cov.get("inferred_sex") or "") or "not inferred"
    sections: list[tuple[str, str]] = []
    parts: list[str] = []

    def add(sec_id: str, title: str, body: str) -> None:
        sections.append((sec_id, title))
        parts.append(_sec(sec_id, title, body))

    # overview
    cards = [("QC", _status_badge(str(qc.get("status", "pass")))), ("sex (sheet)", s.get("sex_sheet") or "–"),
             ("sex (inferred)", inferred), ("mean depth", sv.fmt(cov.get("mean"), 2, "x")),
             ("reads", sv.fmt(reads.get("read_count"))), ("read N50", sv.fmt(reads.get("read_length_n50"), 0, " bp")),
             ("SNVs", sv.fmt(smv.get("snv"))), ("indels", sv.fmt(smv.get("indel"))),
             ("SVs", sv.fmt(sum(sv.finite(v) for v in svd.get("counts", {}).values()))),  # type: ignore[union-attr]
             ("phase block NG50", sv.fmt_bp(ph.get("block_ng50")))]
    run_rows = [[stage, r.get("run_id"), r.get("ugc_wgw_version"), r.get("attempt"), r.get("finished_at") or ""] for stage, r in runs.items()]
    add("sec-overview", "Overview",
        f"<p class='muted'>sample <b>{_esc(sid)}</b> · mode {_esc(d.get('mode'))} · father {_esc(s.get('father_id') or '–')} · "
        f"mother {_esc(s.get('mother_id') or '–')} · generated {_esc(d.get('generated_at'))}</p>"
        + _cards(cards) + "<h3>QC flags</h3>" + _flags_table(qc)
        + "<p class='small muted'>Thresholds: " + _esc(", ".join(f"{k} = {v:g}" for k, v in sorted(qc.get("thresholds", {}).items()))) + "</p>"  # type: ignore[union-attr]
        + "<h3>Runs summarised</h3>" + _table(["stage", "run", "ugc-wgw version", "attempt", "finished"], run_rows, cls=""))

    # reads
    add("sec-reads", "Reads", '<div class="grid"><div>' + _kv([
        ["reads", sv.fmt(reads.get("read_count"))], ["mapped reads", f"{sv.fmt(reads.get('mapped_read_count'))} ({sv.fmt(reads.get('mapped_read_percent'), 2)}%)"],
        ["read length mean / median / N50", f"{sv.fmt(reads.get('read_length_mean'), 0)} / {sv.fmt(reads.get('read_length_median'), 0)} / {sv.fmt(reads.get('read_length_n50'), 0)} bp"],
        ["read quality mean / median", f"{sv.fmt(reads.get('read_quality_mean'), 1)} / {sv.fmt(reads.get('read_quality_median'), 1)}"],
        ["gap-compressed identity mean / median", f"{sv.fmt(reads.get('gap_compressed_identity_mean'), 2)} / {sv.fmt(reads.get('gap_compressed_identity_median'), 2)} %"]])
        + "</div></div>")

    # coverage
    ratios: dict[str, object] = cov.get("ratios", {})  # type: ignore[assignment]
    chrom_rows: list[dict[str, object]] = cov.get("chrom", [])  # type: ignore[assignment]
    track: dict[str, object] = cov.get("track", {})  # type: ignore[assignment]
    cn: dict[str, object] = cov.get("copynum", {})  # type: ignore[assignment]
    segs: dict[str, list[list[float]]] = cn.get("segments", {})  # type: ignore[assignment]
    order = {c: i for i, c in enumerate(sp.PRIMARY_CONTIGS)}
    track_rows = [{"chrom": c, "len": t["len"], "bins": t["bins"], "bin_bp": track.get("bin_bp"), "steps": segs.get(c, [])}
                  for c, t in sorted(track.get("chroms", {}).items(), key=lambda kv: (order.get(kv[0], 99), kv[0]))]  # type: ignore[union-attr]
    cn_rows = [[r["chrom"], r.get("most_common_cn"), sv.fmt_bp(r.get("total_bases"))] for r in cn.get("chrom", [])]  # type: ignore[union-attr]
    cov_cards = [("mean depth (mosdepth)", sv.fmt(cov.get("mean"), 2, "x")),
                 ("haploid coverage (sawfish, GC-corrected)", sv.fmt(cn.get("haploid_coverage"), 1, "x")),
                 ("chrX / autosomes", sv.fmt(ratios.get("chrX_ratio"), 2)), ("chrY / autosomes", sv.fmt(ratios.get("chrY_ratio"), 3)),
                 ("chrY threshold (female max)", sv.fmt(ratios.get("y_threshold"), 2)),
                 ("sex from ratio", ratios.get("sex_from_ratio") or "–"), ("sex reported", inferred)]
    add("sec-coverage", "Coverage",
        _cards(cov_cards)
        + '<div class="grid"><div><h3>Mean depth per chromosome</h3>'
        + sv.bars_h([(str(r["chrom"]), sv.finite(r.get("mean")), "") for r in chrom_rows], digits=2, unit="x", width=300, label_w=60)
        + "</div><div><h3>Most common copy number (sawfish)</h3>"
        + (_table(["chromosome", "copy number", "bases"], cn_rows, cls="") if cn_rows else "<p class='muted'>Not available.</p>")
        + "</div></div>"
        + f"<h3>Depth along the genome ({sv.fmt_bp(track.get('bin_bp')) if track.get('bin_bp') else '–'} bins) with copy-number segments</h3>"
        + sv.legend([("depth per bin", "m1"), ("copy number (0–4 scale)", "m4"), ("mean depth, dashed", "warn")])
        + sv.track(track_rows, baseline=sv.finite(cov.get("mean")) or None))

    # small variants
    roh: dict[str, object] = smv.get("roh", {})  # type: ignore[assignment]
    sn: dict[str, object] = smv.get("sn", {})  # type: ignore[assignment]
    psc: dict[str, object] = smv.get("psc", {})  # type: ignore[assignment]
    st: dict[str, object] = smv.get("st", {})  # type: ignore[assignment]
    by_chrom_roh = sorted(roh.get("by_chrom", {}).items(), key=lambda kv: -sv.finite(kv[1]))[:12]  # type: ignore[union-attr]
    add("sec-small-variants", "Small variants",
        _cards([("SNVs", sv.fmt(smv.get("snv"))), ("indels", sv.fmt(smv.get("indel"))), ("Ts/Tv", sv.fmt(smv.get("tstv"), 2)),
                ("het/hom", sv.fmt(smv.get("hethom"), 2)), ("records (bcftools)", sv.fmt(sn.get("number of records"))),
                ("het genotypes", sv.fmt(psc.get("n_hets"))), ("hom-alt genotypes", sv.fmt(psc.get("n_nonref_hom"))),
                ("singletons", sv.fmt(psc.get("n_singletons"))), ("missing genotypes", sv.fmt(psc.get("n_missing")))])
        + '<div class="grid"><div><h3>Substitution spectrum</h3>' + sv.bars_h([(k, sv.finite(v), "") for k, v in sorted(st.items())], width=300, label_w=60)
        + "</div><div><h3>Indel length (sites)</h3>" + sv.hist([(str(l), sv.finite(v)) for l, v in smv.get("idd", [])], width=420)  # type: ignore[union-attr]
        + "</div><div><h3>Site quality (SNVs and indels)</h3>" + sv.hist([(str(l), sv.finite(v)) for l, v in smv.get("qual", [])], width=420)  # type: ignore[union-attr]
        + "</div><div><h3>Depth at called sites (genotypes)</h3>" + sv.hist([(str(l), sv.finite(v)) for l, v in smv.get("dp", [])], width=420)  # type: ignore[union-attr]
        + "</div></div><h3>Runs of homozygosity (bcftools roh)</h3>"
        + _cards([("runs", sv.fmt(roh.get("n"))), ("total", sv.fmt_bp(roh.get("total_bp"))), ("longest", sv.fmt_bp(roh.get("longest_bp")))])
        + '<div class="grid"><div>' + sv.hist([(str(l), sv.finite(v)) for l, v in roh.get("hist", [])], width=420, y_label="runs")  # type: ignore[union-attr]
        + "</div><div><h3>Longest total per chromosome</h3>" + sv.bars_h([(str(c), sv.finite(v) / 1e6, "Mb") for c, v in by_chrom_roh], digits=1, width=300, label_w=60) + "</div></div>")

    # structural variants
    sizes: dict[str, dict[str, float]] = svd.get("sizes", {})  # type: ignore[assignment]
    size_bins: list[str] = svd.get("size_bins", sp.SV_SIZE_BINS)  # type: ignore[assignment]
    counts: dict[str, object] = svd.get("counts", {})  # type: ignore[assignment]
    add("sec-sv", "Structural variants",
        _cards([(t, sv.fmt(v)) for t, v in counts.items()] + [("records in the VCF", sv.fmt(svd.get("n")))])
        + '<div class="grid"><div><h3>Calls by type</h3>' + sv.bars_h([(t, sv.finite(v), "") for t, v in counts.items()], width=300, label_w=60)
        + "</div><div><h3>Size by type</h3>" + sv.stacked([(t, {b: sv.finite(v) for b, v in sizes.get(t, {}).items()}) for t in sorted(sizes)], size_bins, width=380)
        + "</div></div><h3>Filters</h3>" + (_table(["filter", "records"], [[k, sv.fmt(v)] for k, v in sorted(svd.get("filters", {}).items())], cls="") if svd.get("filters") else "<p class='muted'>Not available.</p>"))  # type: ignore[union-attr]

    # phasing
    pall: dict[str, object] = ph.get("all", {})  # type: ignore[assignment]
    blocks: dict[str, object] = ph.get("blocks", {})  # type: ignore[assignment]
    het = sv.finite(pall.get("num_heterozygous"))
    add("sec-phasing", "Phasing",
        _cards([("phased bases", sv.fmt_bp(ph.get("phased_bp"))), ("block NG50", sv.fmt_bp(ph.get("block_ng50"))),
                ("block N50", sv.fmt_bp(blocks.get("n50"))), ("blocks", sv.fmt(blocks.get("n") or pall.get("num_blocks"))),
                ("heterozygous variants", sv.fmt(pall.get("num_heterozygous"))),
                ("phased", f"{sv.fmt(pall.get('num_phased'))} ({sv.finite(pall.get('num_phased')) / het * 100:.1f}%)" if het else sv.fmt(pall.get("num_phased")))])
        + '<div class="grid"><div><h3>Block length</h3>' + sv.hist([(str(l), sv.finite(v)) for l, v in blocks.get("hist", [])], width=420, y_label="blocks")  # type: ignore[union-attr]
        + "</div><div><h3>Phased heterozygous variants per chromosome</h3>"
        + sv.bars_h([(str(r.get("chromosome")), sv.finite(r.get("num_phased")) / sv.finite(r.get("num_heterozygous")) * 100 if sv.finite(r.get("num_heterozygous")) else 0.0,
                      f"of {sv.fmt(r.get('num_heterozygous'))}") for r in ph.get("chrom", [])], vmax=100, digits=1, unit="%", width=300, label_w=60)  # type: ignore[union-attr]
        + "</div></div><details><summary>per chromosome</summary>"
        + _table(["chromosome", "variants", "heterozygous", "phased", "unphased", "blocks", "phased bp", "block NG50"],
                 [[r.get("chromosome"), sv.fmt(r.get("num_variants")), sv.fmt(r.get("num_heterozygous")), sv.fmt(r.get("num_phased")),
                   sv.fmt(r.get("num_unphased")), sv.fmt(r.get("num_blocks")), sv.fmt_bp(r.get("basepairs_per_block_sum")), sv.fmt_bp(r.get("block_ng50"))]
                  for r in ph.get("chrom", [])]) + "</details>")  # type: ignore[union-attr]

    # tandem repeats
    dro: dict[str, object] = tr.get("dropouts", {})  # type: ignore[assignment]
    loci: list[dict[str, object]] = tr.get("disease_loci", [])  # type: ignore[assignment]
    locus_rows = [[l.get("gene"), l.get("disease"), f"{l.get('chrom')}:{sv.fmt(l.get('pos'))}-{sv.fmt(l.get('end'))}", _clip(str(l.get("motifs") or ""), 36),
                   Raw(f"<b>{_esc(l.get('gt'))}</b>") if l.get("called") else l.get("gt"), _clip(str(l.get("al") or ""), 24), _clip(str(l.get("mc") or ""), 24),
                   l.get("sd"), l.get("dropout") or "", "yes" if l.get("fail_reads") else ""] for l in loci]
    add("sec-tr", "Tandem repeats",
        _cards([("loci genotyped", sv.fmt(tr.get("genotyped"))), ("loci uncalled", sv.fmt(tr.get("uncalled"))),
                ("catalog records", sv.fmt(tr.get("n_records"))), ("dropout loci", sv.fmt(dro.get("n"))),
                ("disease loci", f"{sv.fmt(tr.get('n_named_called'))} called of {sv.fmt(tr.get('n_named'))}")])
        + '<div class="grid"><div><h3>Coverage dropouts by class</h3>' + sv.bars_h([(k, sv.finite(v), "") for k, v in sorted(dro.get("by_class", {}).items())], width=300, label_w=120)  # type: ignore[union-attr]
        + "</div></div><h3>Known disease loci (STRchive set in the catalog)</h3>"
        + "<p class='small muted'>GT, allele lengths (AL, bp), motif counts (MC) and spanning reads (SD) as TRGT reports them; "
        + "fail reads: locus genotyped with the fail reads included (" + _esc(", ".join(tr.get("fail_reads_trids", []))) + ").</p>"  # type: ignore[arg-type]
        + ('<div class="wrap">' + _table(["gene", "disease", "locus", "motifs", "GT", "AL", "MC", "SD", "dropout", "fail reads"], locus_rows) + "</div>" if locus_rows else "<p class='muted'>Not available.</p>"))

    # methylation
    cpg: dict[str, object] = me.get("cpg", {})  # type: ignore[assignment]
    mbt: dict[str, object] = me.get("methbat", {})  # type: ignore[assignment]
    prof: dict[str, object] = me.get("profile", {})  # type: ignore[assignment]
    add("sec-methylation", "Methylation",
        _cards([("CpG sites (combined)", sv.fmt(cpg.get("combined"))), ("hap1", sv.fmt(cpg.get("hap1"))), ("hap2", sv.fmt(cpg.get("hap2"))),
                ("regions methylated", sv.fmt(mbt.get("methylated"))), ("unmethylated", sv.fmt(mbt.get("unmethylated"))),
                ("allele-specific", sv.fmt(mbt.get("asm"))), ("profile regions", sv.fmt(prof.get("n_regions")))])
        + '<div class="grid"><div><h3>Region categories (MethBat profile)</h3>' + sv.bars_h([(k, sv.finite(v), "") for k, v in sorted(prof.get("labels", {}).items())], width=250, label_w=170)  # type: ignore[union-attr]
        + "</div><div><h3>Combined methylation per region (%)</h3>" + sv.hist([(str(l), sv.finite(v)) for l, v in prof.get("hist", [])], width=420, y_label="regions")  # type: ignore[union-attr]
        + "</div></div>")

    # targeted callers
    para: dict[str, object] = tg.get("paraphase", {})  # type: ignore[assignment]
    mito: dict[str, object] = tg.get("mitorsaw", {})  # type: ignore[assignment]
    star: list[dict[str, str]] = tg.get("starphase", [])  # type: ignore[assignment]
    star_meta: dict[str, object] = tg.get("starphase_meta", {})  # type: ignore[assignment]
    kivvi: dict[str, dict[str, object]] = tg.get("kivvi", {})  # type: ignore[assignment]
    para_rows = [[r["region"], r.get("genes"), sv.fmt(r.get("depth_median"), 1), sv.fmt(r.get("total_cn")), sv.fmt(r.get("n_haplotypes")),
                  Raw('<span class="flag-warn">failed for coverage</span>') if r.get("failed") else "ok"] for r in para.get("regions", [])]  # type: ignore[union-attr]
    mito_rows = [[h.get("label"), sv.fmt(h.get("seq_len")), sv.fmt(h.get("num_ref_variants")), sv.fmt(h.get("estimated_abundance"), 3)] for h in mito.get("haplotypes", [])]  # type: ignore[union-attr]
    kivvi_rows = [[{"kiv2": "LPA KIV-2", "d4z4": "DUX4 D4Z4"}.get(k, k), sv.fmt(v.get("allele_cn"), 1) if isinstance(v.get("allele_cn"), (int, float)) else (v.get("allele_cn") or "–"),
                   sv.fmt(v.get("n_complete_alleles")), sv.fmt(v.get("repeat_depth"), 1), sv.fmt(v.get("genome_depth"), 1), v.get("call_status") or "ok"] for k, v in kivvi.items()]
    db_meta = star_meta.get("database", {})
    add("sec-targeted", "Targeted callers",
        "<h3>Paraphase (paralogous regions)</h3>"
        + _cards([("regions", sv.fmt(para.get("n_regions"))), ("with a copy number", sv.fmt(para.get("n_with_cn"))), ("failed for coverage", sv.fmt(para.get("n_failed")))])
        + ("<details><summary>per region</summary>" + _table(["region", "genes", "depth median", "total copy number", "haplotypes", "status"], para_rows) + "</details>" if para_rows else "<p class='muted'>Not available.</p>")
        + "<h3>mitorsaw (mitochondrial haplotypes)</h3>"
        + (_table(["haplotype", "length", "variants vs reference", "estimated abundance"], mito_rows, cls="") if mito_rows else "<p class='muted'>Not available.</p>")
        + (f"<p class='small muted'>fingerprint explanation passing: {_esc(mito.get('fingerprint', {}).get('passing_explanation'))}</p>" if mito.get("fingerprint") else "")  # type: ignore[union-attr]
        + "<h3>StarPhase (pharmacogenomic diplotypes)</h3>"
        + (f"<p class='small muted'>StarPhase {_esc(star_meta.get('version') or '')}; database: " + _esc(", ".join(f"{k} {v}" for k, v in sorted(db_meta.items()))) + "</p>" if star_meta else "")  # type: ignore[union-attr]
        + (_table(["gene", "diplotype"], [[r["gene"], r["diplotype"]] for r in star]) if star else "<p class='muted'>Not available.</p>")
        + "<h3>kivvi (KIV-2 and D4Z4 repeats)</h3>"
        + (_table(["repeat", "copy number", "complete alleles", "repeat depth", "genome depth", "status"], kivvi_rows, cls="") if kivvi_rows else "<p class='muted'>Not available.</p>"))

    # assembly
    asm: dict[str, object] = d.get("assembly", {})  # type: ignore[assignment]
    if asm:
        hap_rows = [[h["name"], sv.fmt_bp(h["stats"].get("total_length")), sv.fmt(h["stats"].get("n_sequences")), sv.fmt_bp(h["stats"].get("n50")) if h["stats"].get("n50") else "–",
                     sv.fmt_bp(h["stats"].get("aun")), sv.fmt(h["paftools"].get("snps")), sv.fmt(h["paftools"].get("indels"))] for h in asm.get("haplotypes", [])]  # type: ignore[union-attr]
        add("sec-assembly", "Assembly", f"<p class='muted'>{'trio-binned' if asm.get('trio') else 'single-sample'} hifiasm assembly</p>"
            + _table(["haplotype", "total length", "contigs", "N50", "auN", "SNVs vs reference", "indels vs reference"], hap_rows, cls=""))

    # messages and files
    msgs: list[str] = d.get("messages", [])  # type: ignore[assignment]
    warns: list[str] = d.get("warnings", [])  # type: ignore[assignment]
    files: dict[str, dict[str, object]] = d.get("files", {})  # type: ignore[assignment]
    file_rows = [[k, Raw(f'<a href="file://{_esc(v.get("path"))}"><code>{_esc(Path(str(v.get("path"))).name)}</code></a>'), _bytes(v.get("bytes"))]  # type: ignore[arg-type]
                 for k, v in sorted(files.items())]
    add("sec-messages", "Messages and files",
        "<h3>Workflow messages</h3>" + ("<ul class='plain'>" + "".join(f"<li>{_esc(m)}</li>" for m in msgs) + "</ul>" if msgs else "<p class='muted'>None.</p>")
        + "<h3>Summary warnings</h3>" + ("<ul class='plain small'>" + "".join(f"<li>{_esc(w)}</li>" for w in warns) + "</ul>" if warns else "<p class='muted'>None: every file was read.</p>")
        + "<h3>Files summarised</h3>" + _table(["what", "file", "size"], file_rows))

    add("sec-software", "Software and workflow", _software_sections(prov, tools, project_url))

    extra = ""
    if cohort_links:
        extra = "<div class='sub'>cohorts</div>" + "".join(f'<a href="{_esc(href)}">{_esc(name)}</a>' for name, href in cohort_links)
    nav = _nav(f"{sid} · {project_name}", sections, extra)
    header = f"<h1>{_esc(sid)} analysis summary</h1>"
    return _page(f"{sid} analysis summary", nav, header + "".join(parts))


# ---- orchestration -----------------------------------------------------------------------------


def _digest_job(args: tuple[Sources, dict[str, float], float]) -> dict[str, object]:
    src, thr, y = args
    return build_digest(src, thr, y)


def y_threshold_of(code: CodeInfo) -> float:
    b = code.references.get("hifi-wdl-resources")
    if isinstance(b, dict) and isinstance(b.get("scalars"), dict):
        try:
            return float(b["scalars"].get("max_norm_female_chrY_depth", 0.1))  # type: ignore[union-attr]
        except (TypeError, ValueError):
            pass
    return 0.1


def run(cfg: Config, db: DB, code: CodeInfo, *, sample_ids: list[str] | None, cohort_id: str | None, mode: str,
        version: str | None, out_dir: Path, force: bool = False, thresholds: dict[str, float] | None = None,
        jobs: int = 1) -> list[Path]:
    """Write the per-sample reports (and the cohort report when asked); returns the report paths written."""
    thr = thresholds or dict(DEFAULT_THRESHOLDS)
    cohort: Cohort | None = None
    if cohort_id:
        cohort = db.get_cohort(cohort_id)
        if cohort is None:
            raise UgcError(f"unknown cohort {cohort_id}")
    targets: list[str] = []
    for sid in (sample_ids or []) + (cohort.members if cohort else []) + ([] if (sample_ids or cohort) else db.list_sample_ids()):
        if sid not in targets:
            targets.append(sid)
    for sid in targets:
        if not db.sample_exists(sid):
            raise UgcError(f"unknown sample {sid}")
    out_dir.mkdir(parents=True, exist_ok=True)
    tools = sp.scan_tool_versions(cfg.code_dir)
    y = y_threshold_of(code)
    project_name = cfg.project_dir.name or "project"
    written: list[Path] = []
    digests: dict[str, dict[str, object]] = {}
    sources: dict[str, Sources] = {}
    to_build: list[str] = []
    for sid in targets:
        src = resolve_sources(db, sid, mode, version)
        if src is None:
            diag(f"{sid}: no successful {mode} run; skipped")
            continue
        sources[sid] = src
        dpath = out_dir / f"{sid}.summary.json"
        if not force and dpath.exists():
            doc = read_json(dpath)
            if digest_is_fresh(doc, src):
                digests[sid] = doc  # type: ignore[assignment]
                continue
        to_build.append(sid)
    if to_build:
        work = [(sources[sid], thr, y) for sid in to_build]
        if jobs > 1 and len(work) > 1:
            with concurrent.futures.ProcessPoolExecutor(max_workers=jobs) as pool:
                results = list(pool.map(_digest_job, work))
        else:
            results = [_digest_job(w) for w in work]
        for sid, digest in zip(to_build, results):
            digests[sid] = digest
    cohort_links = [(c, f"{c}.summary.html") for c in ([cohort.cohort_id] if cohort else [])]
    for i, sid in enumerate(targets):
        if sid not in digests:
            continue
        digest = digests[sid]
        digest["qc"] = qc_flags(digest, thr)  # thresholds may have changed since the digest was written
        dpath = out_dir / f"{sid}.summary.json"
        write_json(dpath, digest)
        hpath = out_dir / f"{sid}.summary.html"
        links = cohort_links or [(c, f"{c}.summary.html") for c in db.sample_cohorts(sid) if (out_dir / f"{c}.summary.html").exists()]
        hpath.write_text(render_sample(digest, tools, project_name=project_name, project_url=cfg.project_url or DEFAULT_PROJECT_URL,
                                       cohort_links=links))
        written.append(hpath)
        diag(f"{sid}: digest {'rebuilt' if sid in to_build else 'reused'}, wrote {hpath.name}")
        if (i + 1) % 50 == 0:
            diag(f"{i + 1}/{len(targets)} samples")
    if cohort:
        from .summary_cohort import build_cohort_digest, render_cohort  # noqa: PLC0415 - keeps the cohort page optional
        cruns = resolve_cohort_runs(db, cohort.cohort_id, mode, version)
        cdigest = build_cohort_digest(cohort, {sid: digests[sid] for sid in cohort.members if sid in digests}, cruns, mode, thr)
        write_json(out_dir / f"{cohort.cohort_id}.summary.json", cdigest)
        cpath = out_dir / f"{cohort.cohort_id}.summary.html"
        cpath.write_text(render_cohort(cdigest, tools, project_name=project_name, project_url=cfg.project_url or DEFAULT_PROJECT_URL))
        written.append(cpath)
        diag(f"{cohort.cohort_id}: wrote {cpath.name} ({len(cdigest.get('members', []))} members)")  # type: ignore[arg-type]
    return written
