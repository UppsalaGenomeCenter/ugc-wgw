"""Inline-SVG chart primitives for `ugc-wgw summary` (docs/guide/07-results.md, "Analysis summary").

Every function returns one `<svg>` element as a string, or a `<p>` when there is nothing to draw. Colours come
from CSS classes defined by the page (`m1`…`m6` for series, `ok`/`warn`/`bad` for states, `ln1`… for lines), so
the charts follow the light and dark themes; tooltips are SVG `<title>` children. Pixel sizes only, no fonts
to load, no scripts.
"""
from __future__ import annotations

import hashlib
import html
import math

SERIES = ["m1", "m2", "m3", "m4", "m5", "m6"]


def esc(v: object) -> str:
    return html.escape("" if v is None else str(v), quote=True)


def finite(v: object) -> float:
    try:
        x = float(v)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return 0.0
    return x if math.isfinite(x) else 0.0


def fmt(v: object, digits: int = 0, unit: str = "") -> str:
    """Thousands-separated number, `–` for None, with an optional unit."""
    if v is None:
        return "–"
    try:
        x = float(v)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return str(v)
    if not math.isfinite(x):
        return "–"
    if digits == 0 and abs(x - round(x)) < 1e-9:
        s = f"{int(round(x)):,}"
    else:
        s = f"{x:,.{digits}f}"
    return f"{s}{unit}"


def fmt_bp(v: object) -> str:
    x = finite(v)
    if x >= 1e9:
        return f"{x / 1e9:.2f} Gb"
    if x >= 1e6:
        return f"{x / 1e6:.1f} Mb"
    if x >= 1e3:
        return f"{x / 1e3:.0f} kb"
    return f"{x:.0f} bp"


def nice_ticks(lo: float, hi: float, n: int = 5) -> list[float]:
    """Tick positions on a 1-2-5 sequence covering [lo, hi]."""
    lo, hi = finite(lo), finite(hi)
    if hi <= lo:
        hi = lo + 1.0
    raw = (hi - lo) / max(n, 1)
    mag = 10 ** math.floor(math.log10(raw)) if raw > 0 else 1.0
    step = next((m * mag for m in (1, 2, 5, 10) if m * mag >= raw), 10 * mag)
    t = math.floor(lo / step) * step
    ticks = []
    while len(ticks) < 50:
        ticks.append(round(t, 10))
        if t >= hi - step * 1e-9:
            break
        t += step
    return ticks


def _svg(width: float, height: float, body: str, cls: str = "chart") -> str:
    return (f'<svg xmlns="http://www.w3.org/2000/svg" class="{cls}" width="{width:.0f}" height="{height:.0f}" '
            f'viewBox="0 0 {width:.0f} {height:.0f}" font-family="sans-serif" font-size="10">{body}</svg>')


def _jitter(key: str, span: float) -> float:
    h = int(hashlib.sha1(key.encode()).hexdigest()[:8], 16)
    return (h / 0xFFFFFFFF - 0.5) * span


# ---- charts ----------------------------------------------------------------------


def bars_h(items: list[tuple[str, float, str]], *, width: int = 560, label_w: int = 150, row_h: int = 16,
           cls: str = "m1", vmax: float | None = None, digits: int = 0, unit: str = "") -> str:
    """Labelled horizontal bars: (label, value, note); the value is printed after the bar."""
    if not items:
        return "<p class='muted'>No data.</p>"
    top = vmax if vmax else max(finite(v) for _, v, _ in items) or 1.0
    parts = []
    for i, (label, value, note) in enumerate(items):
        y = i * row_h
        w = max(0.0, min(finite(value) / top, 1.0)) * width
        text = f"{fmt(value, digits, unit)} {note}".strip()
        parts.append(f'<text x="{label_w - 6}" y="{y + 12}" text-anchor="end">{esc(label)}</text>'
                     f'<rect class="{cls}" x="{label_w}" y="{y + 2}" width="{w:.1f}" height="{row_h - 4}">'
                     f'<title>{esc(label)}: {esc(text)}</title></rect>'
                     f'<text x="{label_w + w + 4:.1f}" y="{y + 12}">{esc(text)}</text>')
    return _svg(label_w + width + 160, len(items) * row_h + 4, "".join(parts))


def hist(bins: list[tuple[str, float]], *, width: int = 560, height: int = 170, cls: str = "m1",
         y_label: str = "", digits: int = 0) -> str:
    """A vertical histogram of pre-binned counts: (label, value)."""
    if not bins or all(finite(v) <= 0 for _, v in bins):
        return "<p class='muted'>No data.</p>"
    left, bottom, top_pad = 52, 44, (18 if y_label else 8)
    plot_w, plot_h = width - left - 8, height - bottom - top_pad
    vmax = max(finite(v) for _, v in bins) or 1.0
    ticks = nice_ticks(0, vmax, 4)
    vmax = max(vmax, ticks[-1])
    slot = plot_w / len(bins)
    bw = slot * 0.8
    rotate = len(bins) > 12 or max(len(str(label)) for label, _ in bins) > 6
    parts = []
    for t in ticks:
        y = top_pad + plot_h - t / vmax * plot_h
        parts.append(f'<line x1="{left}" y1="{y:.1f}" x2="{left + plot_w}" y2="{y:.1f}"/>'
                     f'<text x="{left - 4}" y="{y + 3:.1f}" text-anchor="end">{esc(fmt(t))}</text>')
    for i, (label, value) in enumerate(bins):
        v = finite(value)
        h = v / vmax * plot_h
        x = left + i * slot + (slot - bw) / 2
        y = top_pad + plot_h - h
        parts.append(f'<rect class="{cls}" x="{x:.1f}" y="{y:.1f}" width="{bw:.1f}" height="{h:.1f}">'
                     f'<title>{esc(label)}: {esc(fmt(v, digits))}</title></rect>')
        lx, ly = x + bw / 2, top_pad + plot_h + 12
        if rotate:
            parts.append(f'<text x="{lx:.1f}" y="{ly}" text-anchor="end" transform="rotate(-35 {lx:.1f} {ly})">{esc(label)}</text>')
        else:
            parts.append(f'<text x="{lx:.1f}" y="{ly}" text-anchor="middle">{esc(label)}</text>')
    if y_label:
        parts.append(f'<text x="{left}" y="10" font-size="9">{esc(y_label)}</text>')
    return _svg(width, height, "".join(parts))


def stacked(rows: list[tuple[str, dict[str, float]]], categories: list[str], *, width: int = 560, label_w: int = 150,
            row_h: int = 18, legend: bool = True, percent: bool = False) -> str:
    """Horizontal stacked bars: (label, {category: value}); categories take the series colours in order."""
    rows = [(label, {k: finite(v) for k, v in vals.items()}) for label, vals in rows]
    if not rows or not categories or all(sum(vals.values()) <= 0 for _, vals in rows):
        return "<p class='muted'>No data.</p>"
    totals = [sum(vals.get(c, 0.0) for c in categories) for _, vals in rows]
    top = 1.0 if percent else (max(totals) or 1.0)
    y0 = 18 if legend else 0
    parts = []
    if legend:
        x = label_w
        for i, c in enumerate(categories):
            cls = SERIES[i % len(SERIES)]
            parts.append(f'<rect class="{cls}" x="{x}" y="3" width="10" height="10"/>'
                         f'<text x="{x + 13}" y="12">{esc(c)}</text>')
            x += 13 + 7 * len(str(c)) + 14
    for r, ((label, vals), total) in enumerate(zip(rows, totals)):
        y = y0 + r * row_h
        parts.append(f'<text x="{label_w - 6}" y="{y + 12}" text-anchor="end">{esc(label)}</text>')
        x = float(label_w)
        for i, c in enumerate(categories):
            v = vals.get(c, 0.0)
            if v <= 0:
                continue
            frac = (v / total) if percent and total > 0 else (v / top)
            w = frac * width
            cls = SERIES[i % len(SERIES)]
            share = f" ({v / total * 100:.1f}%)" if total > 0 else ""
            parts.append(f'<rect class="{cls}" x="{x:.1f}" y="{y + 2}" width="{w:.1f}" height="{row_h - 4}">'
                         f'<title>{esc(label)} · {esc(c)}: {esc(fmt(v))}{esc(share)}</title></rect>')
            x += w
        parts.append(f'<text x="{x + 4:.1f}" y="{y + 12}">{esc(fmt(total))}</text>')
    return _svg(label_w + width + 90, y0 + len(rows) * row_h + 4, "".join(parts))


def track(rows: list[dict[str, object]], *, width: int = 900, row_h: int = 24, label_w: int = 48,
          baseline: float | None = None, vmax: float | None = None, steps_max: float = 4.0,
          area_cls: str = "m1", step_cls: str = "ln3", unit: str = "x") -> str:
    """Karyogram rows: per chromosome an area of binned values (depth), a dashed baseline, and optional
    step segments (copy number) drawn on their own 0..steps_max scale. rows: {chrom, len, bins, steps?}."""
    rows = [r for r in rows if r.get("len")]
    if not rows:
        return "<p class='muted'>No data.</p>"
    plot_w = width - label_w - 8
    max_len = max(finite(r["len"]) for r in rows) or 1.0
    values = [finite(v) for r in rows for v in (r.get("bins") or []) if v is not None]  # type: ignore[union-attr]
    top = vmax if vmax else (max(values) if values else 1.0) or 1.0
    if baseline and not vmax:
        top = max(top, baseline * 2.5)
    inner = row_h - 4
    parts = []
    for i, r in enumerate(rows):
        y_base = i * row_h + row_h - 2
        length = finite(r["len"])
        w = length / max_len * plot_w
        parts.append(f'<text x="{label_w - 6}" y="{y_base - 6}" text-anchor="end">{esc(r["chrom"])}</text>'
                     f'<rect class="lane" x="{label_w}" y="{y_base - inner}" width="{w:.1f}" height="{inner}"/>')
        bins = r.get("bins") or []
        bin_bp = finite(r.get("bin_bp")) or (length / max(len(bins), 1))  # type: ignore[arg-type]
        if bins:
            pts = []
            for j, v in enumerate(bins):  # type: ignore[arg-type]
                if v is None:
                    continue
                x = label_w + min(j * bin_bp, length) / max_len * plot_w
                h = min(finite(v) / top, 1.0) * inner
                pts.append((x, y_base - h))
            if pts:
                d = f"M{pts[0][0]:.1f},{y_base:.1f} " + " ".join(f"L{x:.1f},{y:.1f}" for x, y in pts) + f" L{pts[-1][0]:.1f},{y_base:.1f} Z"
                parts.append(f'<path class="{area_cls}" d="{d}"><title>{esc(r["chrom"])}: mean per bin, max {esc(fmt(top, 1))}{esc(unit)}</title></path>')
        if baseline:
            yb = y_base - min(baseline / top, 1.0) * inner
            parts.append(f'<line class="base" x1="{label_w}" y1="{yb:.1f}" x2="{label_w + w:.1f}" y2="{yb:.1f}" stroke-dasharray="3,3"/>')
        steps = r.get("steps") or []
        if steps:
            seg = []
            for s, e, v in steps:  # type: ignore[misc]
                ys = y_base - min(finite(v) / steps_max, 1.0) * inner
                x1 = label_w + finite(s) / max_len * plot_w
                x2 = label_w + finite(e) / max_len * plot_w
                seg.append(f"M{x1:.1f},{ys:.1f} L{x2:.1f},{ys:.1f}")
            parts.append(f'<path class="{step_cls}" d="{" ".join(seg)}"><title>{esc(r["chrom"])}: copy number segments (0–{steps_max:g})</title></path>')
    return _svg(width, len(rows) * row_h + 4, "".join(parts))


def strip(groups: list[dict[str, object]], *, width: int = 560, label_w: int = 170, row_h: int = 30) -> str:
    """One row per metric: points (id, value) with deterministic jitter, the median line, the IQR box, the
    outliers in the `bad` class and an optional threshold line. groups: {label, points, outliers?, threshold?,
    digits?, unit?}."""
    groups = [g for g in groups if g.get("points")]
    if not groups:
        return "<p class='muted'>No data.</p>"
    plot_w = width - 8
    parts = []
    for i, g in enumerate(groups):
        y0 = i * row_h
        yc = y0 + row_h / 2
        pts = [(str(k), finite(v)) for k, v in g["points"] if v is not None]  # type: ignore[union-attr]
        if not pts:
            continue
        vals = sorted(v for _, v in pts)
        thr = g.get("threshold")
        lo = min(vals + ([finite(thr)] if thr is not None else []))
        hi = max(vals + ([finite(thr)] if thr is not None else []))
        if hi <= lo:
            lo, hi = lo - 1.0, hi + 1.0
        pad = (hi - lo) * 0.06
        lo, hi = lo - pad, hi + pad

        def sx(v: float) -> float:
            return label_w + (v - lo) / (hi - lo) * plot_w

        digits = int(finite(g.get("digits")))
        unit = str(g.get("unit") or "")
        q1, med, q3 = vals[len(vals) // 4], vals[len(vals) // 2], vals[(3 * len(vals)) // 4]
        parts.append(f'<text x="{label_w - 6}" y="{yc + 4:.1f}" text-anchor="end">{esc(g["label"])}</text>'
                     f'<line class="axis" x1="{label_w}" y1="{yc:.1f}" x2="{label_w + plot_w}" y2="{yc:.1f}"/>'
                     f'<rect class="iqr" x="{sx(q1):.1f}" y="{yc - 7:.1f}" width="{max(sx(q3) - sx(q1), 1):.1f}" height="14">'
                     f'<title>IQR {esc(fmt(q1, digits, unit))} – {esc(fmt(q3, digits, unit))}</title></rect>'
                     f'<line class="ln2" x1="{sx(med):.1f}" y1="{yc - 9:.1f}" x2="{sx(med):.1f}" y2="{yc + 9:.1f}">'
                     f'<title>median {esc(fmt(med, digits, unit))}</title></line>'
                     f'<text x="{label_w}" y="{y0 + 9}" font-size="8">{esc(fmt(vals[0], digits, unit))}</text>'
                     f'<text x="{label_w + plot_w}" y="{y0 + 9}" font-size="8" text-anchor="end">{esc(fmt(vals[-1], digits, unit))}</text>')
        if thr is not None:
            parts.append(f'<line class="thr" x1="{sx(finite(thr)):.1f}" y1="{y0 + 3}" x2="{sx(finite(thr)):.1f}" y2="{y0 + row_h - 3}" stroke-dasharray="3,2">'
                         f'<title>threshold {esc(fmt(thr, digits, unit))}</title></line>')
        outliers = set(str(x) for x in (g.get("outliers") or []))  # type: ignore[union-attr]
        for key, v in pts:
            cls = "bad" if key in outliers else "pt"
            parts.append(f'<circle class="{cls}" cx="{sx(v):.1f}" cy="{yc + _jitter(key, row_h - 12):.1f}" r="3">'
                         f'<title>{esc(key)}: {esc(fmt(v, digits, unit))}</title></circle>')
    return _svg(label_w + width, len(groups) * row_h + 4, "".join(parts))


def legend(items: list[tuple[str, str]]) -> str:
    """A small inline legend: (label, class)."""
    if not items:
        return ""
    parts = []
    x = 0
    for label, cls in items:
        parts.append(f'<rect class="{cls}" x="{x}" y="3" width="10" height="10"/><text x="{x + 13}" y="12">{esc(label)}</text>')
        x += 13 + 7 * len(label) + 14
    return _svg(x, 16, "".join(parts), cls="legend")
