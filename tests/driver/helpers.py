"""Shared fixtures for the driver tests."""
from __future__ import annotations

import contextlib
import io
import json
import os
import pathlib
import struct
import zlib
import sys
from html.parser import HTMLParser

from ugc_wgw import cli, config, layout
from ugc_wgw.db import DB, RunRecord
from ugc_wgw.util import utc_now

REPO = pathlib.Path(__file__).resolve().parents[2]
FAKE = pathlib.Path(__file__).resolve().parent / "fake_miniwdl.py"
sys.path.insert(0, str(FAKE.parent))
from fake_miniwdl import make_outputs  # noqa: E402

VERSION = (REPO / "VERSION").read_text().strip()


class HtmlChecker(HTMLParser):
    """Balanced tags (SVG shapes are void) and the number of document-level <title> elements."""
    VOID = {"meta", "br", "hr", "img", "input", "link", "line", "rect", "polyline", "path", "circle", "polygon"}

    def __init__(self):
        super().__init__()
        self.stack, self.errors, self.titles = [], [], 0

    def handle_starttag(self, tag, attrs):
        if tag in self.VOID:
            return
        if tag == "title" and "svg" not in self.stack:  # SVG <title> tooltips do not count
            self.titles += 1
        self.stack.append(tag)

    def handle_startendtag(self, tag, attrs):
        pass

    def handle_endtag(self, tag):
        if tag in self.VOID:
            return
        if not self.stack or self.stack[-1] != tag:
            self.errors.append(f"unexpected </{tag}> (open: {self.stack[-3:]})")
        else:
            self.stack.pop()


def make_project(tmp: pathlib.Path, **kw: object) -> config.Config:
    tmp = pathlib.Path(tmp)
    ref = tmp / "ugc_wgw_ref_map.GRCh38_GIABv3.tsv"
    ref.write_text("name\tGRCh38_GIABv3\nscatter_regions\t/dev/null\n")
    cfg_file = tmp / "miniwdl.cfg"
    cfg_file.write_text("[scheduler]\n")
    args = dict(code=REPO, miniwdl=FAKE, cfg=cfg_file, ref_map=ref, poll_interval=0.05, max_inflight=4)
    args.update(kw)
    return config.init_project(tmp / "proj", **args)  # type: ignore[arg-type]


BGZF_EOF = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")


def bgzf_block(payload: bytes) -> bytes:
    """One BGZF block (gzip member with the BC extra field) holding `payload` (at most 64 KB)."""
    comp = zlib.compressobj(6, zlib.DEFLATED, -15)
    cdata = comp.compress(payload) + comp.flush()
    bsize = 18 + len(cdata) + 8 - 1
    return (b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00" + struct.pack("<H", bsize) + cdata
            + struct.pack("<II", zlib.crc32(payload) & 0xFFFFFFFF, len(payload)))


def bam_bytes(*, reads: int = 10, read_len: int = 1000, movie: str = "m84000_240101_000000_s1", sample: str = "S",
              aligned: bool = False, header_only: bool = False, read_group: bool = True) -> bytes:
    """A small, valid uBAM (BGZF blocks + EOF marker) with `reads` unmapped records of `read_len` bases."""
    header = "@HD\tVN:1.6\tSO:unknown\n"
    if read_group:
        header += f"@RG\tID:rg1\tPL:PACBIO\tPU:{movie}\tSM:{sample}\n"
    refs = [("chr20", 64444167)] if aligned else []
    for name, length in refs:
        header += f"@SQ\tSN:{name}\tLN:{length}\n"
    parts = [b"BAM\x01" + struct.pack("<i", len(header)) + header.encode() + struct.pack("<i", len(refs))]
    for name, length in refs:
        parts.append(struct.pack("<i", len(name) + 1) + name.encode() + b"\0" + struct.pack("<i", length))
    if not header_only:
        seq = bytes([0x11] * ((read_len + 1) // 2))   # all-A, 4 bits per base
        qual = bytes([30]) * read_len
        for i in range(reads):
            qname = f"{movie}/{i}/ccs".encode() + b"\0"
            # refID pos l_read_name mapq bin n_cigar_op flag l_seq next_refID next_pos tlen
            core = struct.pack("<iiBBHHHIiii", -1, -1, len(qname), 255, 4680, 0, 4, read_len, -1, -1, 0)
            rec = core + qname + seq + qual
            parts.append(struct.pack("<i", len(rec)) + rec)
    data = b"".join(parts)
    out = b"".join(bgzf_block(data[i:i + 65280]) for i in range(0, len(data), 65280))
    return out + BGZF_EOF


def write_bam(path: pathlib.Path, **kw: object) -> pathlib.Path:
    """Write `bam_bytes(**kw)` to path; `truncated=True` cuts the file short (no EOF marker), `empty=True`
    writes zero bytes, `pbi_reads=N` adds a `.pbi` index declaring N reads."""
    path = pathlib.Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    truncated = bool(kw.pop("truncated", False))
    empty = bool(kw.pop("empty", False))
    pbi_reads = kw.pop("pbi_reads", None)
    data = b"" if empty else bam_bytes(**kw)  # type: ignore[arg-type]
    if truncated:
        data = data[: max(1, len(data) * 2 // 3)]
    path.write_bytes(data)
    if pbi_reads is not None:
        pbi = b"PBI\x01" + struct.pack("<I", 0x00000301) + struct.pack("<H", 0) + struct.pack("<I", int(pbi_reads)) + bytes(18)
        pathlib.Path(str(path) + ".pbi").write_bytes(bgzf_block(pbi) + BGZF_EOF)
    return path


def write_tsv(tmp: pathlib.Path, rows: list[dict[str, str]], name: str = "samples.tsv") -> pathlib.Path:
    """Write a samples TSV; every relative path mentioned is created as a small valid uBAM under tmp/data
    (an existing file is left as it is, so a test may write its own BAM there first)."""
    tmp = pathlib.Path(tmp)
    cols = ["sample_id", "sex", "hifi_reads", "fail_reads", "father_id", "mother_id"]
    extra = sorted({k for r in rows for k in r} - set(cols))
    cols += extra
    data = tmp / "data"
    data.mkdir(exist_ok=True)
    lines = ["\t".join(cols)]
    for r in rows:
        row = dict(r)
        for kind in ("hifi_reads", "fail_reads"):
            paths = []
            for p in (row.get(kind) or "").split(","):
                p = p.strip()
                if not p:
                    continue
                full = p if p.startswith("/") else str(data / p)
                if not pathlib.Path(full).exists():
                    write_bam(pathlib.Path(full), sample=str(row.get("sample_id", "S")))
                paths.append(full)
            row[kind] = ",".join(paths)
        lines.append("\t".join(row.get(c, "") for c in cols))
    path = tmp / name
    path.write_text("\n".join(lines) + "\n")
    return path


def write_ids(tmp: pathlib.Path, ids: list[str], name: str = "ids.txt") -> pathlib.Path:
    path = pathlib.Path(tmp) / name
    path.write_text("\n".join(ids) + "\n")
    return path


def fake_scancel(tmp: pathlib.Path) -> tuple[str, pathlib.Path]:
    """A fake `scancel` that appends its arguments to $UGC_WGW_FAKE_SCANCEL_LOG; returns (bin dir, log path)."""
    bin_dir = pathlib.Path(tmp) / "bin"
    bin_dir.mkdir(exist_ok=True)
    log = pathlib.Path(tmp) / "scancel.log"
    exe = bin_dir / "scancel"
    exe.write_text('#!/bin/sh\necho "$@" >> "$UGC_WGW_FAKE_SCANCEL_LOG"\n')
    exe.chmod(0o755)
    return str(bin_dir), log


def fake_sacct(tmp: pathlib.Path, text: str) -> tuple[str, pathlib.Path]:
    """A fake `sacct` that appends its arguments to $UGC_WGW_FAKE_SACCT_LOG and prints the canned `text`; when
    $UGC_WGW_FAKE_SACCT_REJECT names a field present in the --format, it fails like an old sacct. Returns (bin dir, args log)."""
    bin_dir = pathlib.Path(tmp) / "bin"
    bin_dir.mkdir(exist_ok=True)
    out = pathlib.Path(tmp) / "sacct.out"
    out.write_text(text)
    log = pathlib.Path(tmp) / "sacct.log"
    exe = bin_dir / "sacct"
    exe.write_text('#!/bin/sh\necho "$@" >> "$UGC_WGW_FAKE_SACCT_LOG"\n'
                   'if [ -n "$UGC_WGW_FAKE_SACCT_REJECT" ]; then case " $* " in *"$UGC_WGW_FAKE_SACCT_REJECT"*)\n'
                   '  echo "sacct: error: Invalid field requested: \\"$UGC_WGW_FAKE_SACCT_REJECT\\"" >&2; exit 1;; esac; fi\n'
                   f'cat "{out}"\n')
    exe.chmod(0o755)
    return str(bin_dir), log


def sacct_env(tmp: pathlib.Path, text: str) -> tuple[dict[str, str], pathlib.Path]:
    bin_dir, log = fake_sacct(tmp, text)
    return {"PATH": f"{bin_dir}:{os.environ['PATH']}", "UGC_WGW_FAKE_SACCT_LOG": str(log)}, log


def scancel_env(tmp: pathlib.Path) -> tuple[dict[str, str], pathlib.Path]:
    bin_dir, log = fake_scancel(tmp)
    return {"PATH": f"{bin_dir}:{os.environ['PATH']}", "UGC_WGW_FAKE_SCANCEL_LOG": str(log)}, log


def write_lease(cfg: config.Config, *, host: str, age_seconds: float = 0, pid: int = 1) -> None:
    """Write a lease line as another driver would have, without holding the flock."""
    from ugc_wgw import lease as lease_mod
    from ugc_wgw.util import utc_plus
    seen = utc_plus(-age_seconds)
    cfg.lock_path.write_text(lease_mod.render(pid, host, seen, seen))


def run_cli(args: list[str], env: dict[str, str] | None = None) -> tuple[int, str, str]:
    """Run the CLI in-process, capturing stdout/stderr. Env overrides are applied for the call."""
    old_env = dict(os.environ)
    if env:
        os.environ.update(env)
    out, err = io.StringIO(), io.StringIO()
    try:
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            code = cli.main(args)
    finally:
        os.environ.clear()
        os.environ.update(old_env)
    return code, out.getvalue(), err.getvalue()


def seed_success(cfg: config.Config, db: DB, stage: str, subject_type: str, subject_id: str,
                 inputs: dict[str, object], *, mode: str = "standalone", version: str = VERSION,
                 attempt: int = 1) -> RunRecord:
    """Insert a successful run with realistic outputs.json/out/ produced by the fake's output generator."""
    ns = f"ugc_wgw_{stage}."
    doc = {ns + k: v for k, v in inputs.items()}
    doc.setdefault(ns + "ugc_wgw_version", version)
    stage_path = layout.stage_dir(cfg.results_dir, subject_type, subject_id, version, stage)
    attempt_path = layout.make_attempt_dir(stage_path / f"attempt-{attempt}")
    files = layout.run_files(attempt_path)
    files.inputs.write_text(json.dumps(doc, indent=2, sort_keys=True) + "\n")
    outputs = make_outputs(stage, doc, attempt_path, version)
    files.outputs.write_text(json.dumps(outputs, indent=2, sort_keys=True) + "\n")
    files.run_json.write_text(json.dumps({"dir": str(attempt_path), "outputs": outputs}) + "\n")
    run = RunRecord(run_id=f"{subject_id}-{stage}-a{attempt}-seed", subject_type=subject_type, subject_id=subject_id,
                    stage=stage, mode=mode, ugc_wgw_version=version, run_dir=str(attempt_path), status="pending",
                    attempt=attempt, inputs_path=str(files.inputs), started_at=utc_now())
    db.insert_run(run)
    db.finalize_run(run.run_id, "success", utc_now(), exit_code=0)
    layout.point_current(stage_path, attempt)
    run.status = "success"
    return run


def seed_files(cfg: config.Config, db: DB, stage: str, subject_type: str, subject_id: str, files: dict[str, str | None], *,
               mode: str = "standalone", version: str = VERSION, attempt: int | None = None,
               scalars: dict[str, object] | None = None) -> RunRecord:
    """Insert a successful run whose outputs are hand-written files (for stages the fake engine does not know, or
    outputs it does not fake): `files` maps an output name to the file's text (None = the optional output is absent);
    `scalars` adds non-File outputs; the attempt is the next free one unless given."""
    ns = f"ugc_wgw_{stage}."
    stage_path = layout.stage_dir(cfg.results_dir, subject_type, subject_id, version, stage)
    if attempt is None:
        attempt = max(layout.existing_attempts(stage_path), default=0) + 1
    attempt_path = layout.make_attempt_dir(stage_path / f"attempt-{attempt}")
    run_files = layout.run_files(attempt_path)
    run_files.inputs.write_text(json.dumps({ns + "ugc_wgw_version": version}, indent=2) + "\n")
    outputs: dict[str, object] = {ns + "ugc_wgw_version": version}
    for name, text in files.items():
        if text is None:
            outputs[ns + name] = None
            continue
        p = attempt_path / "out" / name / f"{subject_id}.{stage}.{name}.tsv"
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(text)
        outputs[ns + name] = str(p)
    for k, v in (scalars or {}).items():
        outputs[ns + k] = v
    run_files.outputs.write_text(json.dumps(outputs, indent=2, sort_keys=True) + "\n")
    run_files.run_json.write_text(json.dumps({"dir": str(attempt_path), "outputs": outputs}) + "\n")
    run = RunRecord(run_id=f"{subject_id}-{stage}-a{attempt}-seedfiles", subject_type=subject_type, subject_id=subject_id,
                    stage=stage, mode=mode, ugc_wgw_version=version, run_dir=str(attempt_path), status="pending",
                    attempt=attempt, inputs_path=str(run_files.inputs), started_at=utc_now())
    db.insert_run(run)
    db.finalize_run(run.run_id, "success", utc_now(), exit_code=0)
    layout.point_current(stage_path, attempt)
    run.status = "success"
    return run


def _deep_merge(base: dict, over: dict) -> dict:
    """Dicts merge recursively; an empty dict or list replaces (the way to blank a block)."""
    for k, v in over.items():
        if isinstance(v, dict) and v and isinstance(base.get(k), dict):
            _deep_merge(base[k], v)
        else:
            base[k] = v
    return base


def make_digest(sid: str, **overrides: object) -> dict[str, object]:
    """A minimal per-sample digest of the current schema with the fixtures' targeted callers; `overrides` deep-merge."""
    from ugc_wgw import summary
    from ugc_wgw import summary_parsers as sp
    from fixtures import MITORSAW, PARAPHASE, STARPHASE_JSON, STARPHASE_TSV, KIVVI_NA, KIVVI_OK  # noqa: E402
    d: dict[str, object] = {
        "schema": summary.DIGEST_SCHEMA, "generated_at": "2026-11-05T00:00:00Z", "mode": "standalone",
        "sample": {"id": sid, "sex_sheet": "MALE", "father_id": None, "mother_id": None},
        "runs": {"singleton": {"run_id": f"{sid}-singleton-a1", "ugc_wgw_version": VERSION, "attempt": 1, "finished_at": "2026-11-05T00:00:00Z"}},
        "warnings": [], "messages": [],
        "stats": {"depth_mean": "30.1", "mapped_read_percent": "99.9", "read_quality_median": "31", "inferred_sex": "MALE",
                  "SNV_count": "3800000", "INDEL_count": "700000", "TSTV_ratio": "1.96", "HETHOM_ratio": "1.3", "foo_extra": "7"},
        "reads": {k: 1.0 for k in summary.READ_KEYS},
        "coverage": {"mean": 30.1, "inferred_sex": "MALE", "total": {}, "chrom": [{"chrom": "chr1", "mean": 30.0}, {"chrom": "chrX", "mean": 15.0}],
                     "ratios": {"chrX_ratio": 0.5, "chrY_ratio": 0.4, "sex_from_ratio": "MALE"}, "track": {},
                     "copynum": {"haploid_coverage": 15.0, "chrom": [], "segments": {}}},
        "small_variants": {"snv": 3800000.0, "indel": 700000.0, "tstv": 1.96, "hethom": 1.3, "sn": {"number of records": 4500000},
                           "psc": {"n_hets": 2000000, "n_nonref_hom": 1500000, "n_singletons": 100, "n_missing": 10}, "st": {}, "idd": [],
                           "qual": [], "dp": [], "roh": {"n": 3, "total_bp": 10000000, "longest_bp": 5000000}},
        "sv": {"counts": {"DEL": 10.0, "INS": 12.0, "DUP": 1.0, "INV": 0.0, "BND": 2.0, "SWAP": 0.0}, "vcf_counts": {}, "sizes": {},
               "size_bins": [], "filters": {}, "n": 25},
        "phasing": {"phased_bp": 2.8e9, "block_ng50": 300000.0, "all": {"num_variants": 100, "num_heterozygous": 60, "num_phased": 54,
                    "num_unphased": 6, "num_blocks": 3}, "chrom": [], "blocks": {"n": 3, "n50": 400000, "total_bp": 900000}},
        "tandem_repeats": {"genotyped": 930000.0, "uncalled": 5000.0, "n_records": 935000, "dropouts": {"n": 2, "by_class": {}},
                           "n_named": 70, "n_named_called": 69, "disease_loci": [], "fail_reads_trids": []},
        "methylation": {"cpg": {"combined": 28000000.0, "hap1": 1.0, "hap2": 1.0}, "methbat": {"methylated": 1.0, "unmethylated": 2.0, "asm": 3.0}, "profile": {}},
        "targeted": {"paraphase": sp.parse_paraphase(PARAPHASE), "mitorsaw": sp.parse_mitorsaw(MITORSAW),
                     "starphase": sp.parse_starphase_tsv(STARPHASE_TSV.splitlines()), "starphase_meta": sp.parse_starphase_json(STARPHASE_JSON),
                     "kivvi": {"kiv2": sp.parse_kivvi(KIVVI_NA), "d4z4": sp.parse_kivvi(KIVVI_OK)}},
        "files": {}, "assembly": {}, "provenance": {"ugc_wgw": {"version": VERSION}},
    }
    _deep_merge(d, overrides)  # type: ignore[arg-type]
    d["qc"] = summary.qc_flags(d, dict(summary.DEFAULT_THRESHOLDS))
    return d
