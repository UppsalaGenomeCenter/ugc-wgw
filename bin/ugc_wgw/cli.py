"""`ugc-wgw` command line (docs/DESIGN.md §8.2)."""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

from . import cohorts, config, engine, manifest, progress, report, resources, samples, stage_inputs, summary, usage
from .log import LOGGER
from .db import DB
from .log import Events, setup_logging
from .plan import Selection, compute
from .reconcile import reconcile
from .stages import MODES, STAGES, mode_stages, stage_spec
from .submit import GRACE_SECONDS, Submitter, acquire_lock, dry_run
from .util import UgcError, diag, utc_stamp
from . import inputs as inputs_mod
from . import layout


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="ugc-wgw", description="ugc-pacbio-wgw driver: samples, cohorts, stages, miniwdl runs.")
    p.add_argument("--project", default=".", metavar="DIR", help="project directory holding .ugc-wgw/ (default: .)")
    p.add_argument("-v", "--verbose", action="store_true", help="log INFO to stderr as well as to .ugc-wgw/logs/ugc-wgw.log")
    p.add_argument("--color", default=os.environ.get("UGC_WGW_COLOR", "auto"), metavar="auto|always|never",
                   help="colour the stderr log by meaning (auto: only on a terminal, honouring NO_COLOR; default from UGC_WGW_COLOR)")
    sub = p.add_subparsers(dest="command", required=True)

    s = sub.add_parser("init", help="create <project_dir>/.ugc-wgw/ (config.json, state.sqlite, logs/)")
    s.add_argument("project_dir")
    s.add_argument("--install", metavar="DIR", help="install prefix version dir (<root>/versions/<v> or <root>/current)")
    s.add_argument("--code", metavar="DIR", help="code dir (repo checkout) when not using --install")
    s.add_argument("--miniwdl", metavar="EXE", help="miniwdl executable when not using --install")
    s.add_argument("--cfg", metavar="FILE", help="miniwdl.cfg when not using --install")
    s.add_argument("--venv", metavar="DIR", help="engine venv (for version probing)")
    s.add_argument("--results", metavar="DIR", help="results root (default: the project dir)")
    s.add_argument("--ref-map", required=True, metavar="FILE",
                   help="ugc-wgw reference map TSV of the build (rendered by install-bundle.sh, ref_map= in its report)")
    s.add_argument("--registry", default=config.DEFAULT_REGISTRY, help="registry of the ugc-built images")
    s.add_argument("--profile", metavar="NAME|FILE",
                   help="site profile applied before the flags below: <prefix>/profiles/NAME.json (with --install), else "
                        "backends/hpc/profiles/NAME.json of the code, or a file; keys: " + ", ".join(config.PROFILE_KEYS))
    s.add_argument("--max-inflight", type=int, default=None, help="default concurrent miniwdl runs for submit (default 4)")
    s.add_argument("--poll-interval", type=float, default=None, help="seconds between polls in submit (default 30)")
    s.add_argument("--no-assembly-parents", action="store_true",
                   help="assembly mode never trio-bins, even when both parents are registered")
    s.add_argument("--deepvariant", choices=config.DEEPVARIANT_MODES, default=None,
                   help="small-variant caller of singleton/upstream: cpu (default), gpu (DeepVariant call_variants on "
                        "one GPU) or parabricks (pbrun deepvariant); needs --nv and a GPU partition in the install")
    s.add_argument("--gpu-type", default=None, metavar="TYPE",
                   help="SLURM gres type of the GPU tasks (sbatch --gres gpu:TYPE:N); empty = any GPU (gpu:N)")
    s.add_argument("--parabricks-gpus", type=int, default=None, metavar="N", help="GPUs per Parabricks task (default 4)")

    s = sub.add_parser("samples", help="register and list samples")
    ss = s.add_subparsers(dest="samples_command", required=True)
    a = ss.add_parser("add", help="register samples from a TSV (sample_id, sex, hifi_reads, fail_reads, father_id, mother_id; "
                                  "father_id/mother_id drive trio binning in assembly mode)")
    a.add_argument("tsv")
    a.add_argument("--no-check", action="store_true", help="do not require read paths to exist (nothing is inspected)")
    a.add_argument("--no-inspect", action="store_true",
                   help="require the paths to exist but do not open the BAMs (default: header, first reads, "
                        "BGZF end marker and an estimate of the yield per file)")
    a.add_argument("--drop-empty", action="store_true",
                   help="register a sample without its empty or read-less BAMs instead of refusing the sheet "
                        "(truncated files are always refused)")
    a.add_argument("--replace", action="store_true", help="overwrite already registered samples")
    a = ss.add_parser("check", help="inspect the registered read files again (all samples or the given ids): size, "
                                    "reads, bases, movie and problems per file; exit 1 on a problem")
    a.add_argument("sample_ids", nargs="*", metavar="ID")
    a.add_argument("--stored", action="store_true", help="show what registration recorded instead of reading the files")
    a.add_argument("--json", action="store_true")
    a.add_argument("--tsv", action="store_true")
    a = ss.add_parser("remove", help="unregister samples; refused for cohort members and, without --force, "
                                     "for samples with recorded runs (results on disk are never touched)")
    a.add_argument("sample_ids", nargs="+", metavar="ID")
    a.add_argument("--force", action="store_true", help="also delete the samples' run rows")
    a = ss.add_parser("list", help="list samples with their latest status per stage")
    a.add_argument("--mode", default="standalone", choices=sorted(MODES))
    a.add_argument("--stage", choices=sorted(STAGES))
    a.add_argument("--status", help="comma-separated statuses to keep (e.g. failed,running)")
    a.add_argument("--any-version", action="store_true", help="consider runs of any ugc-wgw version")
    a.add_argument("--json", action="store_true")

    s = sub.add_parser("cohort", help="freeze cohorts")
    cs = s.add_subparsers(dest="cohort_command", required=True)
    a = cs.add_parser("freeze", help="freeze an immutable, ordered sample list under a cohort ID")
    a.add_argument("cohort_id")
    a.add_argument("--samples", required=True, metavar="FILE", help="one sample ID per line")

    for name in ("submit", "retry"):
        s = sub.add_parser(name, help="run stages" if name == "submit" else "re-run failed/cancelled runs of one stage")
        s.add_argument("--mode", default="standalone", choices=sorted(MODES))
        s.add_argument("--stage", choices=sorted(STAGES), required=(name == "retry"))
        g = s.add_mutually_exclusive_group()
        g.add_argument("--samples", nargs="+", metavar="ID")
        g.add_argument("--samples-file", metavar="FILE")
        s.add_argument("--cohort", metavar="ID")
        s.add_argument("--max-inflight", type=int)
        s.add_argument("--poll-interval", type=float)
        s.add_argument("--any-version", action="store_true", help="accept prerequisite outputs from any ugc-wgw version")
        s.add_argument("--dry-run", action="store_true", help="print inputs and miniwdl command lines, change nothing")
        s.add_argument("--takeover", action="store_true",
                       help="take the project lease from a driver on another host that is known to be dead")
        s.add_argument("--grace", type=float, default=GRACE_SECONDS, help=argparse.SUPPRESS)
        if name == "retry":
            s.add_argument("--max-attempts", type=int, default=3)

    s = sub.add_parser("status", help="latest status per subject and stage")
    s.add_argument("--mode", default="standalone", choices=sorted(MODES))
    s.add_argument("--cohort", metavar="ID")
    s.add_argument("--failed", action="store_true")
    s.add_argument("--any-version", action="store_true")
    s.add_argument("--json", action="store_true")

    s = sub.add_parser("logs", help="show the miniwdl workflow log of a run")
    s.add_argument("subject_id")
    s.add_argument("--stage", required=True, choices=sorted(STAGES))
    s.add_argument("--attempt", type=int)
    s.add_argument("--tail", type=int, nargs="?", const=50, metavar="N")
    s.add_argument("-f", "--follow", action="store_true", help="keep printing until the run ends, then its final status")
    s.add_argument("--json", action="store_true", help="show workflow.log.json (JSON lines) instead of workflow.log")

    s = sub.add_parser("progress", help="done/total per stage, what the active runs are doing, and an estimated time left")
    s.add_argument("--mode", default="standalone", choices=sorted(MODES))
    s.add_argument("--stage", choices=sorted(STAGES))
    g = s.add_mutually_exclusive_group()
    g.add_argument("--samples", nargs="+", metavar="ID")
    g.add_argument("--samples-file", metavar="FILE")
    s.add_argument("--cohort", metavar="ID")
    s.add_argument("--any-version", action="store_true")
    s.add_argument("--json", action="store_true")

    s = sub.add_parser("report", help="write a self-contained HTML summary of the project's runs (technical metrics)")
    s.add_argument("--mode", choices=sorted(MODES), help="only runs of this mode (default: all)")
    s.add_argument("--cohort", metavar="ID", help="only the cohort and its members")
    s.add_argument("--any-version", action="store_true", help="runs of every ugc-wgw version (default: the current one)")
    s.add_argument("--sizes", action="store_true", help="also measure each run's out/ (walks the results tree)")
    s.add_argument("--out", metavar="FILE", help="output file (default: <results>/reports/ugc-wgw-report-<stamp>.html)")
    _price_args(s)

    s = sub.add_parser("summary", help="write self-contained HTML analysis summaries: one per sample, one per cohort "
                                       "(what the analysis found, with QC flags)")
    g = s.add_mutually_exclusive_group()
    g.add_argument("--samples", nargs="+", metavar="ID", help="samples to summarise (default: every registered sample)")
    g.add_argument("--samples-file", metavar="FILE")
    s.add_argument("--cohort", metavar="ID", help="also write the cohort summary (covers the cohort's members)")
    s.add_argument("--mode", default="standalone", choices=("standalone", "joint"),
                   help="which runs to summarise: singleton, or upstream + downstream (default: standalone)")
    s.add_argument("--any-version", action="store_true", help="runs of every ugc-wgw version (default: the current one)")
    s.add_argument("--out-dir", metavar="DIR", help="output directory (default: <results>/reports/summary)")
    s.add_argument("--force", action="store_true", help="rebuild the per-sample digests even when they are current")
    s.add_argument("--threshold", action="append", default=[], metavar="KEY=VALUE",
                   help="QC threshold override (keys: " + ", ".join(sorted(summary.DEFAULT_THRESHOLDS)) + "); repeatable")
    s.add_argument("--jobs", type=int, default=1, metavar="N", help="build per-sample digests in N processes")

    s = sub.add_parser("resources", help="what every task will ask SLURM for: declared in the WDL, capped by the site "
                                         "(cpu_max, memory_max), overridden by the resource policy")
    s.add_argument("--stage", choices=sorted(STAGES), help="only tasks of this stage")
    s.add_argument("--changed", action="store_true", help="only tasks a cap or a policy row changes")
    s.add_argument("--json", action="store_true")

    s = sub.add_parser("inputs", help="print the generated inputs JSON for a subject and stage")
    s.add_argument("subject_id")
    s.add_argument("--stage", required=True, choices=sorted(STAGES))
    s.add_argument("--mode", default="standalone", choices=sorted(MODES))
    s.add_argument("--cohort", metavar="ID")
    s.add_argument("--any-version", action="store_true")

    s = sub.add_parser("usage", help="core-hours, GPU-hours, memory, queue wait, disk and a cost estimate per stage, task "
                                     "and sample, from each run's SLURM accounting (accounting.json, read with sacct when "
                                     "the run finished) or, without it, estimated from the workflow log")
    s.add_argument("--mode", choices=sorted(MODES), help="only runs of this mode (default: all)")
    s.add_argument("--cohort", metavar="ID", help="only the cohort and its members")
    g = s.add_mutually_exclusive_group()
    g.add_argument("--samples", nargs="+", metavar="ID", help="only these samples' runs")
    g.add_argument("--samples-file", metavar="FILE")
    s.add_argument("--any-version", action="store_true", help="runs of every ugc-wgw version (default: the current one)")
    s.add_argument("--by", choices=("stage", "task", "subject", "attempt", "all"), default="all",
                   help="which table(s) to print (default: all = stage, task and subject)")
    s.add_argument("--top", type=int, default=20, metavar="N", help="tasks shown, by core-hours (default 20)")
    s.add_argument("--sizes", action="store_true", help="also measure each attempt's out/ and whole directory (walks the results tree)")
    _price_args(s)
    s.add_argument("--collect", action="store_true",
                   help="first read sacct for the finished runs that have no accounting.json, or a partial one (needs sacct)")
    s.add_argument("--refresh", action="store_true", help="re-read sacct for every finished run (implies --collect)")
    f = s.add_mutually_exclusive_group()
    f.add_argument("--json", action="store_true")
    f.add_argument("--tsv", action="store_true", help="tab-separated tables")

    s = sub.add_parser("stage-inputs", help="every input of every stage: type, default, who fills it (the driver, a "
                                            "config.json key, or nobody: yours under stage_inputs), its description and "
                                            "the project's override; warns about overrides submit would refuse")
    g = s.add_mutually_exclusive_group()
    g.add_argument("--stage", choices=sorted(STAGES), help="only this stage")
    g.add_argument("--mode", choices=sorted(MODES), help="only the stages of this mode")
    s.add_argument("--nested", action="store_true",
                   help="also the call-qualified inputs of the tasks and subworkflows inside each entrypoint "
                        "(<call>.<task>.<input>: thread counts, memory, tool options); runs miniwdl in the engine venv")
    s.add_argument("--json", action="store_true")
    return p


def _price_args(s: argparse.ArgumentParser) -> None:
    s.add_argument("--price", action="append", default=[], metavar="KEY=VALUE",
                   help="unit price for the cost columns (keys: " + ", ".join(usage.PRICE_KEYS) + ", currency); repeatable; "
                        "config.json `prices` gives the defaults")
    s.add_argument("--basis", choices=usage.BASES, default="allocated",
                   help="core-hours charged: allocated CPUs (what SLURM gave) or requested (what a cloud would be sized to)")


def _open(a: argparse.Namespace) -> tuple[config.Config, DB, Events, str]:
    cfg = config.load(Path(a.project))
    setup_logging(cfg.logs_dir, a.verbose, a.color)
    db = DB(cfg.db_path)
    for v in db.migrated:
        diag(f"migrated {cfg.db_path} to schema {v}")
    events = Events(db, cfg.logs_dir / "events.jsonl")
    return cfg, db, events, config.read_version(cfg.code_dir)


def _selection(a: argparse.Namespace) -> Selection:
    samples_sel = None
    if getattr(a, "samples", None):
        samples_sel = list(a.samples)
    elif getattr(a, "samples_file", None):
        samples_sel = cohorts.read_ids(Path(a.samples_file))
    return Selection(mode=a.mode, stage=getattr(a, "stage", None), samples=samples_sel, cohort=getattr(a, "cohort", None))


def cmd_init(a: argparse.Namespace) -> int:
    if a.install:
        paths = config.derive_from_install(Path(a.install))
        for override, key in ((a.code, "code"), (a.miniwdl, "miniwdl"), (a.cfg, "cfg"), (a.venv, "venv")):
            if override:
                paths[key] = Path(override)
        code, miniwdl, cfg_file, venv = paths["code"], paths["miniwdl"], paths["cfg"], paths["venv"]
    else:
        if not (a.code and a.miniwdl and a.cfg):
            raise UgcError("init needs --install DIR, or all of --code, --miniwdl and --cfg")
        code, miniwdl, cfg_file = Path(a.code), Path(a.miniwdl), Path(a.cfg)
        venv = Path(a.venv) if a.venv else None
    profile: dict[str, object] = {}
    meta: dict[str, object] = {}
    if a.profile:
        path = config.find_profile(a.profile, install=Path(a.install) if a.install else None, code=code)
        profile = config.load_profile(path)
        meta = config.profile_meta(path, profile)

    def pick(flag: object, key: str, default: object) -> object:
        return flag if flag is not None else profile.get(key, default)

    cfg = config.init_project(
        Path(a.project_dir), code=code, miniwdl=miniwdl, cfg=cfg_file, venv=venv,
        results=Path(a.results) if a.results else None, ref_map=Path(a.ref_map), registry=a.registry,
        max_inflight=int(pick(a.max_inflight, "max_inflight", 4)),  # type: ignore[arg-type]
        poll_interval=float(pick(a.poll_interval, "poll_interval", 30.0)),  # type: ignore[arg-type]
        assembly_use_parents=False if a.no_assembly_parents else bool(profile.get("assembly_use_parents", True)),
        deepvariant=str(pick(a.deepvariant, "deepvariant", "cpu")), gpu_type=str(pick(a.gpu_type, "gpu_type", "")),
        parabricks_gpus=int(pick(a.parabricks_gpus, "parabricks_gpus", 4)),  # type: ignore[arg-type]
        stage_inputs=profile.get("stage_inputs"), prices=profile.get("prices"),  # type: ignore[arg-type]
        summary_thresholds=profile.get("summary_thresholds"), profile=meta,  # type: ignore[arg-type]
    )
    DB(cfg.db_path).close()
    diag(f"initialised {cfg.ugc_wgw_dir} (ugc-wgw version {config.read_version(cfg.code_dir)}, results {cfg.results_dir})")
    for w in config.fs_warnings(cfg):
        diag(f"warning: {w}")
    if meta:
        diag(f"profile {meta['name']} ({meta['path']}) set {', '.join(meta['applied'])}")  # type: ignore[arg-type]
    if not cfg.miniwdl_cfg.exists():
        diag(f"warning: miniwdl cfg does not exist yet: {cfg.miniwdl_cfg}")
    return 0


def cmd_samples(a: argparse.Namespace) -> int:
    cfg, db, events, version = _open(a)
    try:
        if a.samples_command == "add":
            n = samples.add_samples(db, events, Path(a.tsv), check_paths=not a.no_check, replace=a.replace,
                                    inspect=not a.no_inspect, thresholds=cfg.input_thresholds, drop_empty=a.drop_empty)
            diag(f"registered {n} sample(s)")
            return 0
        if a.samples_command == "check":
            rows, problems, warnings = samples.check_samples(db, list(a.sample_ids), cfg.input_thresholds, stored=a.stored)
            _print_rows(rows, a.json, tsv=a.tsv)
            for w in warnings:
                diag(f"warning: {w}")
            for pr in problems:
                diag(f"problem: {pr}")
            diag(f"{len(rows)} file(s) of {len({r['sample_id'] for r in rows})} sample(s): {len(problems)} problem(s), "
                 f"{len(warnings)} warning(s)")
            return 1 if problems else 0
        if a.samples_command == "remove":
            for r in samples.remove_samples(db, events, list(a.sample_ids), force=a.force, results_dir=cfg.results_dir):
                diag(f"removed {r['sample_id']} ({r['runs_deleted']} run row(s) deleted; results on disk kept: {r['results_kept']})")
            return 0
        code = manifest.load_code_info(cfg.code_dir)
        reconcile(cfg, db, events, code, {})
        statuses = [s.strip() for s in a.status.split(",")] if a.status else None
        rows = samples.list_samples(db, a.mode, None if a.any_version else version, a.stage, statuses)
        _print_rows(rows, a.json)
        return 0
    finally:
        db.close()


def cmd_cohort(a: argparse.Namespace) -> int:
    cfg, db, events, _ = _open(a)
    try:
        cohort = cohorts.freeze(db, events, a.cohort_id, Path(a.samples))
        diag(f"froze cohort {cohort.cohort_id}: {len(cohort.members)} samples, sha256 {cohort.sample_list_sha256}")
        return 0
    finally:
        db.close()


def cmd_submit(a: argparse.Namespace, retry: bool = False) -> int:
    cfg, db, events, version = _open(a)
    try:
        sel = _selection(a)
        code = manifest.load_code_info(cfg.code_dir)
        max_attempts = getattr(a, "max_attempts", None)
        res = resources.load(cfg)  # a malformed policy is refused here, before anything is launched
        for w in res.warnings:
            diag(f"warning: {w}")
        if a.dry_run:
            reconcile(cfg, db, events, code, {}, dry=True)
            print(f"# {res.summary()}")
            return dry_run(cfg, db, version, sel, any_version=a.any_version, retry_failed=retry, max_attempts=max_attempts)
        if not cfg.miniwdl_cfg.exists():
            raise UgcError(f"miniwdl cfg not found: {cfg.miniwdl_cfg}")
        lock = acquire_lock(cfg, takeover=a.takeover)
        try:
            engine_info = engine.probe(cfg)
            if res.rules and engine_info.get("ugc_wgw_miniwdl", "unknown") == "unknown":
                diag("warning: the engine venv has no ugc_wgw_resources task plugin: the policy rows will not apply")
            LOGGER.info(res.summary())
            diag(res.summary())
            reconcile(cfg, db, events, code, engine_info)
            sub = Submitter(cfg, db, events, code, engine_info, version, sel,
                            max_inflight=a.max_inflight or cfg.max_inflight,
                            poll=a.poll_interval or cfg.poll_interval, any_version=a.any_version,
                            retry_failed=retry, max_attempts=max_attempts, grace=a.grace, lock=lock)
            return sub.run_loop()
        finally:
            lock.close()
    finally:
        db.close()


def cmd_status(a: argparse.Namespace) -> int:
    cfg, db, events, version = _open(a)
    try:
        code = manifest.load_code_info(cfg.code_dir)
        reconcile(cfg, db, events, code, {})
        lookup = None if a.any_version else version
        latest = db.latest_per_stage(lookup)
        stages = list(mode_stages(a.mode))
        cohort = db.get_cohort(a.cohort) if a.cohort else None
        if a.cohort and cohort is None:
            raise UgcError(f"unknown cohort {a.cohort}")
        rows: list[dict[str, object]] = []
        sample_ids = cohort.members if cohort else db.list_sample_ids()
        for sid in sample_ids:
            row: dict[str, object] = {"subject": sid, "type": "sample"}
            for s in stages:
                if stage_spec(s).subject_type != "sample":
                    continue
                run = latest.get(("sample", sid, s))
                row[s] = f"{run.status}" + (f"({run.attempt})" if run and run.attempt > 1 else "") if run else "-"
            rows.append(row)
        cohort_stages = [s for s in stages if stage_spec(s).subject_type == "cohort"]
        for c in ([cohort] if cohort else db.list_cohorts()) if cohort_stages else []:
            row = {"subject": c.cohort_id, "type": "cohort"}
            for s in stages:
                if stage_spec(s).subject_type != "cohort":
                    continue
                run = latest.get(("cohort", c.cohort_id, s))
                row[s] = f"{run.status}" + (f"({run.attempt})" if run and run.attempt > 1 else "") if run else "-"
            rows.append(row)
        if a.failed or a.json:
            for row in rows:
                kind, message = "", ""
                for s in stages:
                    run = latest.get((str(row["type"]), str(row["subject"]), s))
                    if run is not None and run.status in ("failed", "cancelled"):
                        kind, message = run.error_kind or "", run.error_message or ""
                        break
                row["kind"], row["message"] = kind, message
        if a.failed:
            rows = [r for r in rows if any(str(v).startswith(("failed", "cancelled")) for k, v in r.items() if k in stages)]
        _print_rows(rows, a.json)
        return 0
    finally:
        db.close()


def _follow(cfg: config.Config, db: DB, events: Events, run_id: str, path: Path, tail: int | None) -> int:
    """Print `path` as it grows until the run is terminal; reconcile every 60 s so a dead driver's run ends too."""
    code = manifest.load_code_info(cfg.code_dir)
    pos = 0
    printed_tail = tail is None
    ticks = 0
    while True:
        if path.exists():
            with open(path, errors="replace") as fh:
                if not printed_tail:
                    lines = fh.read().splitlines(keepends=True)
                    for line in lines[-tail:]:  # type: ignore[operator]
                        print(line, end="", flush=True)
                    pos = fh.tell()
                    printed_tail = True
                else:
                    fh.seek(pos)
                    while True:
                        line = fh.readline()
                        if not line or not line.endswith("\n"):
                            break  # wait for the rest of a partial line
                        pos = fh.tell()
                        print(line, end="", flush=True)
        current = db.get_run(run_id)
        if current is None or not current.is_active:
            if current is not None:
                print(f"# {run_id} finished: status={current.status} class={current.error_class or ''} "
                      f"kind={current.error_kind or ''} message={current.error_message or ''}", flush=True)
            return 0
        ticks += 1
        if ticks % 60 == 0:
            reconcile(cfg, db, events, code, {})
        time.sleep(1.0)


def cmd_logs(a: argparse.Namespace) -> int:
    cfg, db, events, _ = _open(a)
    try:
        subject_type = stage_spec(a.stage).subject_type
        runs = db.runs_for(subject_type, a.subject_id, a.stage)
        if a.attempt:
            runs = [r for r in runs if r.attempt == a.attempt]
        if not runs:
            raise UgcError(f"no {a.stage} run recorded for {subject_type} {a.subject_id}")
        run = runs[-1]
        files = layout.run_files(run.run_path)
        if a.json:
            path = files.workflow_log_json
        else:
            path = files.workflow_log if files.workflow_log.exists() or a.follow else files.stderr
        print(f"# {run.run_id} status={run.status} attempt={run.attempt} dir={run.run_dir}")
        print(f"# {path}")
        if a.follow:
            return _follow(cfg, db, events, run.run_id, path, a.tail)
        if not path.exists():
            diag("no log written yet" if run.is_active else f"no {path.name} for this run")
            return 0
        lines = path.read_text(errors="replace").splitlines()
        if a.tail:
            lines = lines[-a.tail:]
        for line in lines:
            print(line)
        return 0
    finally:
        db.close()


def cmd_progress(a: argparse.Namespace) -> int:
    cfg, db, events, version = _open(a)
    try:
        code = manifest.load_code_info(cfg.code_dir)
        reconcile(cfg, db, events, code, {}, dry=True)
        sel = _selection(a)
        plan = compute(db, cfg, sel, version, any_version=a.any_version)
        summary = progress.summarize(db, cfg, plan, version, mode=sel.mode, max_inflight=cfg.max_inflight)
        if a.json:
            print(json.dumps(summary.as_detail(), indent=2, sort_keys=True))
        else:
            print(summary.format())
        return 0
    finally:
        db.close()


def cmd_report(a: argparse.Namespace) -> int:
    cfg, db, events, version = _open(a)
    try:
        code = manifest.load_code_info(cfg.code_dir)
        reconcile(cfg, db, events, code, {})
        if a.cohort and db.get_cohort(a.cohort) is None:
            raise UgcError(f"unknown cohort {a.cohort}")
        text = report.build(cfg, db, code, ugc_wgw_version=None if a.any_version else version, mode=a.mode,
                            cohort=a.cohort, sizes=a.sizes, prices=usage.parse_prices(cfg, a.price), basis=a.basis)
        out = Path(a.out) if a.out else cfg.results_dir / "reports" / f"ugc-wgw-report-{utc_stamp()}.html"
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(text)
        diag(f"wrote {out}")
        print(out)
        return 0
    finally:
        db.close()


def cmd_summary(a: argparse.Namespace) -> int:
    cfg, db, events, version = _open(a)
    try:
        code = manifest.load_code_info(cfg.code_dir)
        reconcile(cfg, db, events, code, {})
        sel = _selection(a)
        thr = summary.thresholds_for(cfg, a.threshold)
        out_dir = Path(a.out_dir) if a.out_dir else cfg.results_dir / "reports" / "summary"
        written = summary.run(cfg, db, code, sample_ids=sel.samples, cohort_id=a.cohort, mode=a.mode,
                              version=None if a.any_version else version, out_dir=out_dir, force=a.force,
                              thresholds=thr, jobs=max(1, a.jobs))
        for p in written:
            print(p)
        if not written:
            diag("nothing written: no sample has a successful run to summarise")
            return 1
        return 0
    finally:
        db.close()


def cmd_resources(a: argparse.Namespace) -> int:
    cfg = config.load(Path(a.project))
    res = resources.load(cfg)
    for w in res.warnings:
        diag(f"warning: {w}")
    effs = res.rows(stage=a.stage, changed=a.changed)
    if not a.json:
        print(f"# {res.summary()}")
        print("# columns: declared = cpu / memory in the WDL; cpu, memory, time, partition, constraint = what sbatch is asked "
              "for, with (clamp) for the site caps and (policy) for a policy row; gpus = --gres gpu:[type:]N of the GPU "
              "tasks (config.json deepvariant, gpu_type, parabricks_gpus); rules = the rows that applied")
    _print_rows(resources.as_rows(res, effs), a.json)
    return 0


def cmd_inputs(a: argparse.Namespace) -> int:
    cfg, db, _, version = _open(a)
    try:
        spec = stage_spec(a.stage)
        cohort = db.get_cohort(a.cohort) if a.cohort else None
        if a.cohort and cohort is None:
            raise UgcError(f"unknown cohort {a.cohort}")
        if spec.subject_type == "cohort":
            if cohort is None:
                cohort = db.get_cohort(a.subject_id)
            if cohort is None:
                raise UgcError(f"unknown cohort {a.subject_id}")
        ctx = inputs_mod.BuildContext(cfg, db, a.mode, version, spec.subject_type, a.subject_id, cohort=cohort,
                                      any_version=a.any_version)
        doc, _ = inputs_mod.generate(ctx, a.stage)
        print(json.dumps(doc, indent=2, sort_keys=True))
        return 0
    finally:
        db.close()


def cmd_usage(a: argparse.Namespace) -> int:
    cfg, db, events, version = _open(a)
    try:
        code = manifest.load_code_info(cfg.code_dir)
        reconcile(cfg, db, events, code, {})
        prices = usage.parse_prices(cfg, a.price)
        samples_sel = list(a.samples) if a.samples else cohorts.read_ids(Path(a.samples_file)) if a.samples_file else None
        ver = None if a.any_version else version
        runs = usage.select_runs(db, mode=a.mode, cohort=a.cohort, samples=samples_sel, ugc_wgw_version=ver)
        if a.collect or a.refresh:
            statuses = usage.collect_missing(cfg, db, events, runs, refresh=a.refresh)
            diag("accounting: " + (", ".join(f"{k} {v}" for k, v in sorted(statuses.items())) or "no finished runs"))
        summary_, attempts = usage.build(cfg, db, mode=a.mode, cohort=a.cohort, samples=samples_sel, ugc_wgw_version=ver,
                                         sizes=a.sizes)
        for note in summary_.notes:
            diag(f"note: {note}")
        if a.json:
            print(json.dumps(usage.as_json(summary_, attempts, prices, a.basis), indent=2))
            return 0
        tables: list[tuple[str, list[dict[str, object]]]] = []
        if a.by in ("stage", "all"):
            tables.append(("stage", usage.rows_stage(summary_, prices, a.basis, mode=a.mode, sizes=a.sizes)))
        if a.by in ("task", "all"):
            tables.append((f"task, top {a.top} by core-hours", usage.rows_task(summary_, prices, a.basis, top=a.top)))
        if a.by in ("subject", "all"):
            tables.append(("subject", usage.rows_subject(summary_, prices, a.basis, sizes=a.sizes)))
        if a.by == "attempt":
            tables.append(("attempt", usage.rows_attempts(attempts, prices, a.basis, sizes=a.sizes)))
        if not a.tsv:
            print(f"# {usage.summary_line(summary_, prices, a.basis)}")
        for i, (title, rows) in enumerate(tables):
            if not a.tsv or len(tables) > 1:
                print(("\n" if i or not a.tsv else "") + f"# by {title}")
            _print_rows(rows, False, tsv=a.tsv)
        return 0
    finally:
        db.close()


def cmd_stage_inputs(a: argparse.Namespace) -> int:
    cfg = config.load(Path(a.project))
    if a.stage:
        stages = [a.stage]
    elif a.mode:
        stages = list(mode_stages(a.mode))
    else:
        stages = list(STAGES)
    rows, warnings = stage_inputs.listing(cfg, stages, nested=a.nested)
    for w in warnings:
        diag(f"warning: {w}")
    if a.json:
        print(json.dumps([r.as_dict() for r in rows], indent=2))
        return 0
    print(f"# {stage_inputs.summary(rows, stages, nested=a.nested)}")
    print(f"# {stage_inputs.COLUMNS_HELP}")
    for stage in stages:
        mine = [r for r in rows if r.stage == stage]
        own = sum(1 for r in mine if not r.nested)
        extra = f", {len(mine) - own} call-qualified" if a.nested else ""
        print(f"\n# {stage} ({stage_spec(stage).wdl}): {own} inputs{extra}")
        _print_rows([r.as_table_row() for r in mine], False)
    return 0


def _print_rows(rows: list[dict[str, object]], as_json: bool, *, tsv: bool = False) -> None:
    if as_json:
        print(json.dumps(rows, indent=2, sort_keys=True))
        return
    if not rows:
        print("(no rows)")
        return
    cols: list[str] = []
    for r in rows:
        for c in r:
            if c not in cols:
                cols.append(c)
    if tsv:
        print("\t".join(cols))
        for r in rows:
            print("\t".join(str(r.get(c, "")) for c in cols))
        return
    widths = {c: max(len(c), *(len(str(r.get(c, ""))) for r in rows)) for c in cols}
    print("  ".join(c.ljust(widths[c]) for c in cols))
    for r in rows:
        print("  ".join(str(r.get(c, "")).ljust(widths[c]) for c in cols))


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    a = parser.parse_args(argv)
    try:
        if a.command == "init":
            return cmd_init(a)
        if a.command == "samples":
            return cmd_samples(a)
        if a.command == "cohort":
            return cmd_cohort(a)
        if a.command == "submit":
            return cmd_submit(a)
        if a.command == "retry":
            return cmd_submit(a, retry=True)
        if a.command == "status":
            return cmd_status(a)
        if a.command == "logs":
            return cmd_logs(a)
        if a.command == "inputs":
            return cmd_inputs(a)
        if a.command == "report":
            return cmd_report(a)
        if a.command == "summary":
            return cmd_summary(a)
        if a.command == "progress":
            return cmd_progress(a)
        if a.command == "resources":
            return cmd_resources(a)
        if a.command == "stage-inputs":
            return cmd_stage_inputs(a)
        if a.command == "usage":
            return cmd_usage(a)
        parser.error(f"unknown command {a.command}")
    except UgcError as exc:
        diag(f"error: {exc}")
        return 1
    except KeyboardInterrupt:
        diag("interrupted")
        return 130
    return 2
