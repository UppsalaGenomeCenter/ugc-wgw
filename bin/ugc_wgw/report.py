"""`ugc-wgw report`: a self-contained HTML summary of a project's runs (technical metrics, not analysis results).

Sources: the state DB (runs, events), each attempt's run_manifest.json, and miniwdl's workflow.log.json,
whose NOTICE lines pair `task setup` with `done` / `done (cached)` per task logger and carry the
`runtime.cpu adjusted to host limit` warnings. No external assets: inline CSS, inline SVG, a few lines
of inline JavaScript for sortable tables. See docs/guide/07-results.md ("Run report").
"""
from __future__ import annotations

import html
import json
import os
import statistics
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path

from .config import Config
from .db import DB, RunRecord
from .layout import run_files
from .manifest import CodeInfo
from .stages import mode_stages
from .util import parse_utc, read_json, utc_now

STATUS_COLOR = {"success": "#2e7d32", "failed": "#c62828", "cancelled": "#ef6c00",
                "running": "#1565c0", "submitted": "#5c6bc0", "pending": "#9e9e9e"}


# ---- data ----------------------------------------------------------------------

@dataclass
class TaskStat:
    name: str
    call_id: str
    start: float
    end: float | None = None
    cached: bool = False
    cpu_requested: int | None = None
    cpu_granted: int | None = None
    policy: str = ""              # the ugc_wgw_resources plugin's summary, e.g. "cpu 64→48, partition -→fat"
    retries: int = 0
    failed: bool = False

    @property
    def duration(self) -> float | None:
        return None if self.end is None else max(0.0, self.end - self.start)


@dataclass
class RunReport:
    run: RunRecord
    tasks: list[TaskStat] = field(default_factory=list)
    ignored_keys: Counter = field(default_factory=Counter)
    manifest: dict[str, object] = field(default_factory=dict)
    out_bytes: int | None = None

    @property
    def duration(self) -> float | None:
        if not self.run.started_at or not self.run.finished_at:
            return None
        return (parse_utc(self.run.finished_at) - parse_utc(self.run.started_at)).total_seconds()


def parse_workflow_log(path: Path) -> tuple[list[TaskStat], Counter]:
    """Per-task timings from miniwdl's JSON-lines workflow log; the task logger name is the key."""
    tasks: dict[str, TaskStat] = {}
    ignored: Counter = Counter()
    try:
        fh = open(path, errors="replace")
    except OSError:
        return [], ignored
    with fh:
        for line in fh:
            try:
                d = json.loads(line)
            except ValueError:
                continue
            src = str(d.get("source", ""))
            msg = str(d.get("message", ""))
            ts = d.get("timestamp")
            if not isinstance(ts, (int, float)) or ".t:" not in src:
                continue
            if msg == "task setup":
                call_id = src.rsplit(".t:", 1)[-1]
                tasks[src] = TaskStat(str(d.get("name", call_id)), call_id, float(ts))
            elif src in tasks:
                t = tasks[src]
                if msg.startswith("done"):
                    t.end = float(ts)
                    t.cached = "cached" in msg
                elif msg == "runtime.cpu adjusted to host limit":
                    t.cpu_requested, t.cpu_granted = _int(d.get("original")), _int(d.get("adjusted"))
                elif msg == "ugc-wgw resource policy applied":
                    t.policy = str(d.get("summary") or "")
                elif msg == "failed task will be retried":
                    t.retries += 1
                elif msg == "ignored runtime settings":
                    for k in d.get("keys") or []:
                        ignored[str(k)] += 1
                elif d.get("level") == "ERROR" and "failed" in msg:
                    t.failed = True
    return list(tasks.values()), ignored


def _int(v: object) -> int | None:
    try:
        return int(v)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None


def _dir_bytes(path: Path) -> int:
    total = 0
    for root, _dirs, files in os.walk(path):
        for f in files:
            try:
                total += os.lstat(os.path.join(root, f)).st_size
            except OSError:
                pass
    return total


def collect(cfg: Config, db: DB, *, mode: str | None, cohort: str | None, ugc_wgw_version: str | None,
            sizes: bool = False) -> list[RunReport]:
    members: set[str] | None = None
    if cohort:
        c = db.get_cohort(cohort)
        members = set(c.members) if c else set()
    out: list[RunReport] = []
    for run in db.all_runs(ugc_wgw_version):
        if mode and run.mode != mode:
            continue
        if cohort and not (run.cohort_id == cohort or run.subject_id == cohort or run.subject_id in (members or set())):
            continue
        rr = RunReport(run)
        files = run_files(run.run_path)
        rr.tasks, rr.ignored_keys = parse_workflow_log(files.workflow_log_json)
        if files.manifest.exists():
            try:
                doc = read_json(files.manifest)
                rr.manifest = doc if isinstance(doc, dict) else {}
            except (OSError, ValueError):
                rr.manifest = {}
        if sizes and (run.run_path / "out").exists():
            rr.out_bytes = _dir_bytes(run.run_path / "out")
        out.append(rr)
    return out


# ---- formatting helpers ----------------------------------------------------------

def _esc(v: object) -> str:
    return html.escape("" if v is None else str(v))


def _dur(seconds: float | None) -> str:
    if seconds is None:
        return "–"
    s = int(round(seconds))
    if s < 60:
        return f"{s}s"
    if s < 3600:
        return f"{s // 60}m {s % 60:02d}s"
    return f"{s // 3600}h {(s % 3600) // 60:02d}m"


def _bytes(n: int | None) -> str:
    if n is None:
        return "–"
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024 or unit == "TB":
            return f"{n:.0f} {unit}" if unit == "B" else f"{n:.1f} {unit}"
        n /= 1024  # type: ignore[assignment]
    return str(n)


def _ids(ids: list[object]) -> str:
    if not ids:
        return ""
    shown = ", ".join(str(i) for i in ids[:4])
    return shown if len(ids) <= 4 else f"{len(ids)}: {shown}, …"


def _epoch(ts: str | None) -> float | None:
    return parse_utc(ts).timestamp() if ts else None


def _stats(values: list[float]) -> tuple[str, str, str, str]:
    if not values:
        return ("–", "–", "–", "–")
    return (_dur(statistics.median(values)), _dur(statistics.fmean(values)), _dur(max(values)), _dur(sum(values)))


def _table(headers: list[str], rows: list[list[object]], cls: str = "sortable") -> str:
    th = "".join(f"<th>{_esc(h)}</th>" for h in headers)
    body = "".join("<tr>" + "".join(f"<td>{c if isinstance(c, Raw) else _esc(c)}</td>" for c in r) + "</tr>" for r in rows)
    return f'<table class="{cls}"><thead><tr>{th}</tr></thead><tbody>{body}</tbody></table>'


class Raw(str):
    """A cell whose HTML is already escaped/built."""


def _nice_tick(span: float) -> float:
    for step in (60, 300, 900, 1800, 3600, 7200, 21600, 43200, 86400, 172800, 604800):
        if span / step <= 12:
            return float(step)
    return 604800.0


def _time_axis(t0: float, t1: float, x0: float, width: float, y: float) -> str:
    span = max(1.0, t1 - t0)
    step = _nice_tick(span)
    parts = []
    t = t0 - (t0 % step) + step
    while t <= t1:
        x = x0 + (t - t0) / span * width
        label = _dur(t - t0)
        parts.append(f'<line x1="{x:.1f}" y1="0" x2="{x:.1f}" y2="{y:.1f}" stroke="#ddd"/>'
                     f'<text x="{x:.1f}" y="{y + 12:.1f}" font-size="10" text-anchor="middle">+{label}</text>')
        t += step
    return "".join(parts)


# ---- charts --------------------------------------------------------------------------

def timeline_svg(reports: list[RunReport]) -> str:
    rows = [r for r in reports if r.run.started_at]
    if not rows:
        return "<p>No started runs.</p>"
    now = utc_now()
    rows.sort(key=lambda r: (r.run.started_at or ""))
    t0 = min(_epoch(r.run.started_at) or 0 for r in rows)
    t1 = max((_epoch(r.run.finished_at or now) or 0) for r in rows)
    span = max(1.0, t1 - t0)
    label_w, width, row_h = 260, 700, 14
    height = len(rows) * row_h + 24
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{label_w + width + 10}" height="{height}" '
             f'font-family="sans-serif">', _time_axis(t0, t1, label_w, width, len(rows) * row_h)]
    for i, r in enumerate(rows):
        y = i * row_h
        s = _epoch(r.run.started_at) or t0
        e = _epoch(r.run.finished_at or now) or s
        x = label_w + (s - t0) / span * width
        w = max(2.0, (e - s) / span * width)
        color = STATUS_COLOR.get(r.run.status, "#999")
        label = f"{r.run.subject_id} / {r.run.stage} #{r.run.attempt}"
        title = f"{label}: {r.run.status}, {_dur(r.duration)}, {r.run.started_at} → {r.run.finished_at or 'running'}"
        parts.append(f'<text x="{label_w - 6}" y="{y + 11}" font-size="10" text-anchor="end">{_esc(label)}</text>'
                     f'<rect x="{x:.1f}" y="{y + 2}" width="{w:.1f}" height="{row_h - 4}" fill="{color}"><title>{_esc(title)}</title></rect>')
    parts.append("</svg>")
    return "".join(parts)


def concurrency_svg(reports: list[RunReport]) -> str:
    """Runs and tasks in flight over time (step charts)."""
    now = utc_now()
    run_iv = [(_epoch(r.run.started_at), _epoch(r.run.finished_at or now)) for r in reports if r.run.started_at]
    task_iv = [(t.start, t.end if t.end is not None else (_epoch(now) or t.start)) for r in reports for t in r.tasks]
    if not run_iv:
        return "<p>No started runs.</p>"
    t0 = min(s for s, _ in run_iv if s is not None)
    t1 = max(e for _, e in run_iv if e is not None)
    span = max(1.0, t1 - t0)
    label_w, width, height = 40, 920, 120

    def series(iv: list[tuple[float | None, float | None]], color: str) -> tuple[str, int]:
        events: list[tuple[float, int]] = []
        for s, e in iv:
            if s is None or e is None:
                continue
            events.append((s, 1))
            events.append((e, -1))
        events.sort(key=lambda e: (e[0], -e[1]))  # at equal times, starts before ends
        pts, level, peak = [(t0, 0)], 0, 0
        for t, d in events:
            pts.append((t, level))
            level += d
            peak = max(peak, level)
            pts.append((t, level))
        pts.append((t1, level))
        return (pts, color), peak  # type: ignore[return-value]

    (run_pts, _), run_peak = series(run_iv, STATUS_COLOR["running"])
    (task_pts, _), task_peak = series(task_iv, "#8e24aa")
    peak = max(1, run_peak, task_peak)

    def poly(pts: list[tuple[float, int]], color: str) -> str:
        coords = " ".join(f"{label_w + (t - t0) / span * width:.1f},{height - v / peak * (height - 20):.1f}" for t, v in pts)
        return f'<polyline points="{coords}" fill="none" stroke="{color}" stroke-width="1.5"/>'

    return (f'<svg xmlns="http://www.w3.org/2000/svg" width="{label_w + width + 10}" height="{height + 16}" font-family="sans-serif">'
            f'{_time_axis(t0, t1, label_w, width, height)}'
            f'<text x="2" y="12" font-size="10">{peak}</text><text x="2" y="{height}" font-size="10">0</text>'
            f'{poly(run_pts, STATUS_COLOR["running"])}{poly(task_pts, "#8e24aa")}'
            f'<text x="{label_w + 8}" y="12" font-size="10" fill="{STATUS_COLOR["running"]}">runs in flight (peak {run_peak})</text>'
            f'<text x="{label_w + 180}" y="12" font-size="10" fill="#8e24aa">tasks in flight (peak {task_peak})</text></svg>')


def bars_svg(items: list[tuple[str, float, str]]) -> str:
    """Horizontal bars: (label, value in seconds, note)."""
    if not items:
        return "<p>No finished tasks.</p>"
    top = max(v for _, v, _ in items) or 1.0
    label_w, width, row_h = 260, 600, 16
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{label_w + width + 120}" height="{len(items) * row_h + 4}" font-family="sans-serif">']
    for i, (label, v, note) in enumerate(items):
        y = i * row_h
        w = v / top * width
        parts.append(f'<text x="{label_w - 6}" y="{y + 12}" font-size="10" text-anchor="end">{_esc(label)}</text>'
                     f'<rect x="{label_w}" y="{y + 2}" width="{w:.1f}" height="{row_h - 4}" fill="#607d8b"/>'
                     f'<text x="{label_w + w + 4:.1f}" y="{y + 12}" font-size="10">{_esc(_dur(v))} {_esc(note)}</text>')
    parts.append("</svg>")
    return "".join(parts)


def progress_svg(events: list[dict[str, object]]) -> str:
    """done/total over time from submit.progress events."""
    pts = []
    for e in events:
        if e["event"] != "submit.progress":
            continue
        d = e["detail"]
        if isinstance(d, dict) and d.get("total"):
            pts.append((_epoch(str(e["ts"])) or 0.0, 100.0 * float(d["done"]) / float(d["total"]), d.get("eta_seconds")))
    if len(pts) < 2:
        return "<p>No progress samples (the driver emits submit.progress while it runs).</p>"
    t0, t1 = pts[0][0], pts[-1][0]
    span = max(1.0, t1 - t0)
    label_w, width, height = 40, 920, 120
    coords = " ".join(f"{label_w + (t - t0) / span * width:.1f},{height - v / 100.0 * (height - 20):.1f}" for t, v, _ in pts)
    last = pts[-1]
    return (f'<svg xmlns="http://www.w3.org/2000/svg" width="{label_w + width + 10}" height="{height + 16}" font-family="sans-serif">'
            f'{_time_axis(t0, t1, label_w, width, height)}'
            f'<text x="2" y="12" font-size="10">100%</text><text x="2" y="{height}" font-size="10">0%</text>'
            f'<polyline points="{coords}" fill="none" stroke="#2e7d32" stroke-width="1.5"/>'
            f'<text x="{label_w + 8}" y="12" font-size="10" fill="#2e7d32">done, last {last[1]:.0f}% (eta then {_dur(last[2]) if last[2] is not None else "–"})</text></svg>')


# ---- the page --------------------------------------------------------------------------

CSS = """
:root{--bg:#f7f8f9;--fg:#1d2329;--muted:#5a6470;--line:#d6dbe0;--panel:#ffffff;--head:#eaeef2;--zebra:#f2f4f6;--accent:#1d5f6f;--grid:#d6dbe0}
@media (prefers-color-scheme: dark){:root:not([data-theme="light"]){--bg:#15191d;--fg:#e4e8ec;--muted:#9aa5b1;--line:#333b44;--panel:#1c2126;--head:#232a31;--zebra:#1a1f24;--accent:#6fb7c7;--grid:#333b44}}
:root[data-theme="dark"]{--bg:#15191d;--fg:#e4e8ec;--muted:#9aa5b1;--line:#333b44;--panel:#1c2126;--head:#232a31;--zebra:#1a1f24;--accent:#6fb7c7;--grid:#333b44}
body{font-family:-apple-system,"Segoe UI",Helvetica,Arial,sans-serif;margin:0;padding:24px 16px 48px;color:var(--fg);background:var(--bg);max-width:1200px;font-variant-numeric:tabular-nums}
h1{font-size:22px;margin:0 0 4px;letter-spacing:-0.01em}h2{font-size:15px;margin-top:36px;border-bottom:1px solid var(--line);padding-bottom:4px;color:var(--accent);text-transform:uppercase;letter-spacing:0.06em}
p{max-width:70ch}table{border-collapse:collapse;font-size:12px;margin:8px 0;background:var(--panel)}th,td{border:1px solid var(--line);padding:3px 6px;text-align:left;vertical-align:top}
th{background:var(--head);cursor:pointer;user-select:none;font-weight:600}tr:nth-child(even) td{background:var(--zebra)}
.cards{display:flex;flex-wrap:wrap;gap:12px}.card{border:1px solid var(--line);background:var(--panel);padding:8px 14px;min-width:120px}
.card b{display:block;font-size:20px}.muted{color:var(--muted)}code{font-size:11px}.wrap{overflow-x:auto}
.kind-transient{color:#1565c0}.kind-resource{color:#ef6c00}.kind-input{color:#8e24aa}.kind-tool{color:#c62828}
details summary{cursor:pointer}a{color:var(--accent)}
svg text{fill:var(--muted)}svg line{stroke:var(--grid)}
"""

JS = """
document.querySelectorAll('table.sortable th').forEach(function(th){th.addEventListener('click',function(){
var t=th.closest('table'),i=Array.from(th.parentNode.children).indexOf(th),rows=Array.from(t.tBodies[0].rows);
var asc=!(th.dataset.asc==='1');th.dataset.asc=asc?'1':'0';
rows.sort(function(a,b){var x=a.cells[i].textContent,y=b.cells[i].textContent,nx=parseFloat(x),ny=parseFloat(y);
var c=(!isNaN(nx)&&!isNaN(ny))?nx-ny:x.localeCompare(y);return asc?c:-c;});
rows.forEach(function(r){t.tBodies[0].appendChild(r);});});});
"""


def render(cfg: Config, code: CodeInfo, reports: list[RunReport], *, ugc_wgw_version: str, mode: str | None,
           cohort: str | None, events: list[dict[str, object]], n_samples: int, n_cohorts: int) -> str:
    runs = [r.run for r in reports]
    by_status = Counter(r.status for r in runs)
    kinds = Counter(r.error_kind for r in runs if r.status in ("failed", "cancelled"))
    all_tasks = [t for r in reports for t in r.tasks]
    finished_tasks = [t for t in all_tasks if t.duration is not None]
    cached = sum(1 for t in all_tasks if t.cached)
    auto_retries = sum(1 for e in events if e["event"] == "run.auto_retry")
    durations = [r.duration for r in reports if r.duration is not None]
    starts = [r.started_at for r in runs if r.started_at]
    ends = [r.finished_at for r in runs if r.finished_at]
    span = None
    if starts and ends:
        span = (parse_utc(max(ends)) - parse_utc(min(starts))).total_seconds()
    engine = next((r.manifest.get("engine") for r in reversed(reports) if r.manifest.get("engine")), {}) or {}
    hosts = sorted({r.host for r in runs if r.host})

    # stage table
    stage_rows = []
    stages = [s for s in (mode_stages(mode) if mode else sorted({r.stage for r in runs}))]
    for stage in stages:
        rs = [r for r in reports if r.run.stage == stage]
        if not rs:
            continue
        st = Counter(r.run.status for r in rs)
        subjects_done = len({r.run.subject_id for r in rs if r.run.status == "success"})
        d = [r.duration for r in rs if r.duration is not None and r.run.status == "success"]
        med, mean, mx, total = _stats(d)  # type: ignore[arg-type]
        n_tasks = sum(len(r.tasks) for r in rs)
        n_cached = sum(1 for r in rs for t in r.tasks if t.cached)
        stage_rows.append([stage, len(rs), subjects_done, st["success"], st["failed"], st["cancelled"],
                           st["running"] + st["submitted"] + st["pending"], sum(1 for r in rs if r.run.attempt > 1),
                           med, mean, mx, total, n_tasks, n_cached])

    # task table
    per_task: dict[str, list[TaskStat]] = defaultdict(list)
    for t in all_tasks:
        per_task[t.name].append(t)
    task_rows = []
    for name, ts in sorted(per_task.items(), key=lambda kv: -sum(t.duration or 0 for t in kv[1] if not t.cached)):
        d = [t.duration for t in ts if t.duration is not None and not t.cached]
        med, mean, mx, total = _stats(d)  # type: ignore[arg-type]
        adj = [t for t in ts if t.cpu_granted is not None]
        cpu = f"{adj[0].cpu_requested}→{adj[0].cpu_granted} ({len(adj)})" if adj else ""
        pol = [t for t in ts if t.policy]
        policy = f"{pol[0].policy} ({len(pol)})" if pol else ""
        task_rows.append([name, len(ts), sum(1 for t in ts if t.cached), sum(t.retries for t in ts),
                          sum(1 for t in ts if t.failed), med, mean, mx, total, cpu, policy])
    slowest = []
    for name, ts in per_task.items():
        d = [t.duration for t in ts if t.duration is not None and not t.cached]
        if d:
            slowest.append((name, statistics.median(d), f"median of {len(d)}"))
    slowest.sort(key=lambda x: -x[1])

    # runs table
    run_rows = []
    for r in sorted(reports, key=lambda r: (r.run.subject_id, r.run.stage, r.run.attempt)):
        manifest_link = ""
        files = run_files(r.run.run_path)
        if files.manifest.exists():
            manifest_link = f'<a href="file://{_esc(files.manifest)}">manifest</a>'
        kind = r.run.error_kind or ""
        run_rows.append([r.run.subject_id, r.run.subject_type, r.run.stage, r.run.mode, r.run.attempt,
                         Raw(f'<span style="color:{STATUS_COLOR.get(r.run.status, "#333")}">{_esc(r.run.status)}</span>'),
                         Raw(f'<span class="kind-{_esc(kind)}">{_esc(kind)}</span>'), r.run.started_at or "", _dur(r.duration),
                         "" if r.run.exit_code is None else r.run.exit_code,
                         f"{len(r.tasks)} ({sum(1 for t in r.tasks if t.cached)} cached)",
                         _ids(list(r.manifest.get("slurm_job_ids", []) or [])), _bytes(r.out_bytes) if r.out_bytes is not None else "",
                         Raw(manifest_link), r.run.error_message or ""])

    # failures
    fail_rows = []
    for r in reports:
        if r.run.status not in ("failed", "cancelled"):
            continue
        err = r.manifest.get("error") or {}
        ev = err.get("evidence") or {}
        fail_rows.append([r.run.subject_id, r.run.stage, r.run.attempt, r.run.status, r.run.error_kind or "",
                          r.run.error_class or "", err.get("node") or "", err.get("exit_status") if isinstance(err, dict) else "",
                          Raw(f"<code>{_esc(ev.get('line') or '')}</code>"), err.get("task_dir") or "", r.run.not_before or ""])

    # resource adjustments: miniwdl's caps and the ugc_wgw_resources policy
    cpu_adj = Counter((t.name, t.cpu_requested, t.cpu_granted) for t in all_tasks if t.cpu_granted is not None)
    adjust_rows = [[n, "cap (cpu_max)", f"{a}→{b}", c] for (n, a, b), c in sorted(cpu_adj.items())]
    adjust_rows += [[n, "policy", summ, c] for (n, summ), c in sorted(Counter((t.name, t.policy) for t in all_tasks if t.policy).items())]
    ignored: Counter = Counter()
    for r in reports:
        ignored.update(r.ignored_keys)

    # driver sessions and notable events
    session_rows = []
    for e in events:
        if e["event"] in ("submit.start", "submit.stop"):
            d = e["detail"]
            session_rows.append([e["ts"], e["event"], d.get("mode", ""), d.get("cohort", "") or "", d.get("max_inflight", ""),
                                 d.get("launched", ""), d.get("failures", ""), d.get("interrupted", "")])
    notable = [e for e in events if e["event"] in ("run.reconciled", "run.scancel", "run.auto_retry", "sample.removed")
               or e["level"] in ("warning", "error")]

    def section(title: str, body: str) -> str:
        return f"<h2>{_esc(title)}</h2>{body}"

    cards = [("samples", n_samples), ("cohorts", n_cohorts), ("runs", len(runs)), ("success", by_status["success"]),
             ("failed", by_status["failed"]), ("cancelled", by_status["cancelled"]),
             ("active", by_status["running"] + by_status["submitted"] + by_status["pending"]),
             ("automatic retries", auto_retries), ("tasks", len(all_tasks)),
             ("cache hits", f"{cached} ({(100 * cached // len(all_tasks)) if all_tasks else 0}%)"),
             ("run time, summed", _dur(sum(durations)) if durations else "–"), ("campaign span", _dur(span))]
    cards_html = '<div class="cards">' + "".join(f'<div class="card"><span class="muted">{_esc(k)}</span><b>{_esc(v)}</b></div>' for k, v in cards) + "</div>"
    scope = f"mode {mode or 'all'}, cohort {cohort or 'all'}, version {ugc_wgw_version or 'any'}"
    name = cfg.project_dir.name or "project"
    header = (f"<h1>{_esc(name)} run report</h1><p class='muted'>{_esc(cfg.project_dir)} · results {_esc(cfg.results_dir)} · "
              f"{_esc(scope)} · generated {utc_now()}</p>")
    prov_rows = [["ugc-pacbio-wgw", f"{code.version} ({code.git_commit})"]]
    for name, src in sorted(code.upstream.items()):
        if isinstance(src, dict):
            prov_rows.append([f"upstream {name}", f"{src.get('tag', '')} ({str(src.get('commit', ''))[:12]})"])
    for name, b in sorted(code.references.items()):
        if isinstance(b, dict):
            prov_rows.append([f"reference bundle {name}", str(b.get("version", ""))])
    for k, v in sorted(engine.items()):
        prov_rows.append([f"engine {k}", str(v)])
    prov_rows.append(["hosts", ", ".join(hosts)])
    prov_rows.append(["containers", f"{len(code.containers)} image digests (see the manifests)"])

    parts = [
        f"<!DOCTYPE html><html lang='en'><head><meta charset='utf-8'><title>{_esc(name)} run report</title>",
        f"<style>{CSS}</style></head><body>", header, cards_html,
        section("Provenance", _table(["what", "value"], prov_rows, cls="")),
        section("Stages", '<div class="wrap">' + _table(["stage", "runs", "subjects done", "success", "failed", "cancelled", "active",
                                                          "re-attempts", "median", "mean", "max", "total", "tasks", "cached"], stage_rows) + "</div>"),
        section("Timeline", '<div class="wrap">' + timeline_svg(reports) + "</div>"),
        section("Concurrency", '<div class="wrap">' + concurrency_svg(reports) + "</div>"),
        section("Progress over time", '<div class="wrap">' + progress_svg(events) + "</div>"),
        section("Slowest tasks (median wall time, cache hits excluded)", '<div class="wrap">' + bars_svg(slowest[:20]) + "</div>"),
        section("Tasks", '<div class="wrap">' + _table(["task", "calls", "cached", "miniwdl retries", "failed", "median", "mean",
                                                         "max", "total", "cpu requested→granted (n)", "policy (n)"], task_rows) + "</div>"),
        section("Runs", '<div class="wrap">' + _table(["subject", "type", "stage", "mode", "attempt", "status", "kind", "started",
                                                        "duration", "exit", "tasks", "slurm jobs", "out size", "manifest", "message"], run_rows) + "</div>"),
        section("Failures", '<div class="wrap">' + (_table(["subject", "stage", "attempt", "status", "kind", "class", "node",
                                                             "exit status", "evidence", "task dir", "not before"], fail_rows) if fail_rows else "<p>None.</p>") + "</div>"),
        section("Resource adjustments",
                _table(["task", "by", "change", "occurrences"], adjust_rows)
                if adjust_rows else "<p>None: no request exceeded the site caps and no policy row matched.</p>"
                if all_tasks else "<p>No task logs.</p>"),
        (f"<p class='muted'>ignored runtime keys: {_esc(', '.join(f'{k} ({v})' for k, v in ignored.most_common()))}</p>" if ignored else ""),
        section("Driver sessions", _table(["time", "event", "mode", "cohort", "max_inflight", "launched", "failures", "interrupted"], session_rows)
                if session_rows else "<p>None.</p>"),
        section("Notable events", "<details><summary>" + f"{len(notable)} events (reconciliations, cancellations, automatic retries, warnings, errors)" + "</summary>" +
                _table(["time", "level", "event", "run", "detail"], [[e["ts"], e["level"], e["event"], e["run_id"] or "",
                                                                     Raw(f"<code>{_esc(json.dumps(e['detail'], sort_keys=True))[:400]}</code>")] for e in notable]) + "</details>"),
        f"<script>{JS}</script></body></html>",
    ]
    return "\n".join(parts)


def build(cfg: Config, db: DB, code: CodeInfo, *, ugc_wgw_version: str | None, mode: str | None = None,
          cohort: str | None = None, sizes: bool = False) -> str:
    reports = collect(cfg, db, mode=mode, cohort=cohort, ugc_wgw_version=ugc_wgw_version, sizes=sizes)
    events = db.all_events()
    return render(cfg, code, reports, ugc_wgw_version=ugc_wgw_version or "", mode=mode, cohort=cohort, events=events,
                  n_samples=len(db.list_sample_ids()), n_cohorts=len(db.list_cohorts()))
