"""Shared fixtures for the driver tests."""
from __future__ import annotations

import contextlib
import io
import json
import os
import pathlib
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


def write_tsv(tmp: pathlib.Path, rows: list[dict[str, str]], name: str = "samples.tsv") -> pathlib.Path:
    """Write a samples TSV; every path mentioned is created as an empty file under tmp/data."""
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
                pathlib.Path(full).parent.mkdir(parents=True, exist_ok=True)
                pathlib.Path(full).touch()
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
