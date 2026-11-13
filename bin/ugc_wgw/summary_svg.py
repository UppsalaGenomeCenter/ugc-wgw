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
         y_label: str = "", digits: int = 0, cumulative: bool = False) -> str:
    """A vertical histogram of pre-binned counts: (label, value). `cumulative` adds the running share of the
    total as a line on a right-hand 0-100 % axis."""
    if not bins or all(finite(v) <= 0 for _, v in bins):
        return "<p class='muted'>No data.</p>"
    left, bottom, top_pad = 52, 44, (18 if y_label else 8)
    right = 40 if cumulative else 8
    plot_w, plot_h = width - left - right, height - bottom - top_pad
    vmax = max(finite(v) for _, v in bins) or 1.0
    if vmax <= 8 and all(float(finite(v)).is_integer() for _, v in bins):
        ticks = [float(t) for t in range(0, int(vmax) + 1)]   # counts: one tick per integer, no repeated labels
    else:
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
    if cumulative:
        total = sum(finite(v) for _, v in bins) or 1.0
        run = 0.0
        pts = []
        for i, (_, value) in enumerate(bins):
            run += finite(value)
            pts.append(f"{left + i * slot + slot / 2:.1f},{top_pad + plot_h - run / total * plot_h:.1f}")
        parts.append(f'<polyline class="ln2" points="{" ".join(pts)}"><title>cumulative share</title></polyline>')
        for pct in (0, 50, 100):
            y = top_pad + plot_h - pct / 100 * plot_h
            parts.append(f'<text x="{left + plot_w + 4}" y="{y + 3:.1f}" font-size="9">{pct}%</text>')
    if y_label:
        parts.append(f'<text x="{left}" y="10" font-size="9">{esc(y_label)}</text>')
    return _svg(width, height, "".join(parts))


BAND_FROM = 150   # strip() draws every member below this many points; band() above


def band(groups: list[dict[str, object]], *, width: int = 560, label_w: int = 170, row_h: int = 38, nbins: int = 60) -> str:
    """strip() for large cohorts: one row per metric with the 5-95 % whisker, the interquartile box, the median
    line, a density of all members in `nbins` bins under the axis, and only the outliers and threshold failures
    as points (hover for the id). Same `groups` as strip(); the full values stay in the digest."""
    groups = [g for g in groups if g.get("points")]
    if not groups:
        return "<p class='muted'>No data.</p>"
    plot_w = width - 8
    parts = []
    for i, g in enumerate(groups):
        y0 = i * row_h
        yc = y0 + 13
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
        n = len(vals)
        q = lambda f: vals[min(n - 1, int(f * n))]  # noqa: E731
        p5, q1, med, q3, p95 = q(0.05), q(0.25), q(0.5), q(0.75), q(0.95)
        parts.append(f'<text x="{label_w - 6}" y="{yc + 4:.1f}" text-anchor="end">{esc(g["label"])}</text>'
                     f'<line class="axis" x1="{label_w}" y1="{yc:.1f}" x2="{label_w + plot_w}" y2="{yc:.1f}"/>'
                     f'<line class="axis" x1="{sx(p5):.1f}" y1="{yc:.1f}" x2="{sx(p95):.1f}" y2="{yc:.1f}" stroke-width="2">'
                     f'<title>5-95 % {esc(fmt(p5, digits, unit))} – {esc(fmt(p95, digits, unit))}</title></line>'
                     f'<rect class="iqr" x="{sx(q1):.1f}" y="{yc - 7:.1f}" width="{max(sx(q3) - sx(q1), 1):.1f}" height="14">'
                     f'<title>IQR {esc(fmt(q1, digits, unit))} – {esc(fmt(q3, digits, unit))} (n={n})</title></rect>'
                     f'<line class="ln2" x1="{sx(med):.1f}" y1="{yc - 9:.1f}" x2="{sx(med):.1f}" y2="{yc + 9:.1f}">'
                     f'<title>median {esc(fmt(med, digits, unit))}</title></line>'
                     f'<text x="{label_w}" y="{y0 + 8}" font-size="8">{esc(fmt(vals[0], digits, unit))}</text>'
                     f'<text x="{label_w + plot_w}" y="{y0 + 8}" font-size="8" text-anchor="end">{esc(fmt(vals[-1], digits, unit))}</text>')
        # density: members per bin as bars under the axis, 12 px tall at the fullest bin
        counts = [0] * nbins
        for v in vals:
            counts[min(nbins - 1, int((v - lo) / (hi - lo) * nbins))] += 1
        cmax = max(counts) or 1
        bw = plot_w / nbins
        for b, c in enumerate(counts):
            if c:
                h = max(1.0, c / cmax * 12)
                parts.append(f'<rect class="pt" x="{label_w + b * bw:.1f}" y="{yc + 9:.1f}" width="{max(bw - 0.5, 0.5):.1f}" height="{h:.1f}">'
                             f'<title>{c} member{"s" if c != 1 else ""} in {esc(fmt(lo + b * (hi - lo) / nbins, digits, unit))} – '
                             f'{esc(fmt(lo + (b + 1) * (hi - lo) / nbins, digits, unit))}</title></rect>')
        if thr is not None:
            parts.append(f'<line class="thr" x1="{sx(finite(thr)):.1f}" y1="{y0 + 3}" x2="{sx(finite(thr)):.1f}" y2="{y0 + row_h - 3}" stroke-dasharray="3,2">'
                         f'<title>threshold {esc(fmt(thr, digits, unit))}</title></line>')
        flagged = set(str(x) for x in (g.get("outliers") or [])) | set(str(x) for x in (g.get("beyond") or []))  # type: ignore[union-attr]
        for key, v in pts:
            if key in flagged:
                parts.append(f'<circle class="bad" cx="{sx(v):.1f}" cy="{yc + _jitter(key, 10):.1f}" r="3">'
                             f'<title>{esc(key)}: {esc(fmt(v, digits, unit))}</title></circle>')
    return _svg(label_w + width, len(groups) * row_h + 4, "".join(parts))


def strip_or_band(groups: list[dict[str, object]], **kw: object) -> str:
    """strip() below BAND_FROM points in the largest group, band() from there."""
    n = max((len(g.get("points") or []) for g in groups), default=0)  # type: ignore[arg-type]
    if n >= BAND_FROM:
        kw.pop("row_h", None)
        return band(groups, **kw)  # type: ignore[arg-type]
    return strip(groups, **kw)  # type: ignore[arg-type]


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


RIBBON_COLUMNS = 200   # ribbon() draws one column per member up to here, then the mean of adjacent members per column


def ribbon(columns: list[tuple[str, dict[str, float]]], categories: list[str], *, width: int = 700, height: int = 150,
           marked: set[str] | None = None, max_columns: int = RIBBON_COLUMNS) -> str:
    """A composition ribbon: one thin stacked column per member, in the order given (the caller sorts), each column
    its categories as shares of 100%, the categories in the series colours. Beyond `max_columns` members a column is
    the mean of adjacent members, so the picture and the page stay the same size at any cohort size. `marked`
    members are never averaged away: they stand at the right, after a gap, one wider column each with a red tick.
    (label, {category: fraction}); fractions need not sum to 1."""
    cols = [(label, {c: finite(vals.get(c)) for c in categories}) for label, vals in columns]
    cols = [(label, vals) for label, vals in cols if sum(vals.values()) > 0]
    if not cols or not categories:
        return "<p class='muted'>No data.</p>"
    marked = marked or set()
    plain = [c for c in cols if c[0] not in marked]
    odd = [c for c in cols if c[0] in marked]
    per = max(1, math.ceil(len(plain) / max_columns))
    groups = [plain[i:i + per] for i in range(0, len(plain), per)]
    left, top, legend_h = 34, 6, 18
    odd_w = min(12.0, max(4.0, width * 0.25 / len(odd))) if odd else 0.0
    gap = 10.0 if odd and groups else 0.0
    w = (width - gap - odd_w * len(odd)) / max(1, len(groups))
    parts = []
    x = float(left)
    for i, c in enumerate(categories):
        parts.append(f'<rect class="{SERIES[i % len(SERIES)]}" x="{x:.0f}" y="3" width="10" height="10"/><text x="{x + 13:.0f}" y="12">{esc(c)}</text>')
        x += 13 + 7 * len(str(c)) + 14
    y0 = top + legend_h
    for tick in (0, 50, 100):
        y = y0 + height * (1 - tick / 100)
        parts.append(f'<text x="{left - 4}" y="{y + 3:.1f}" text-anchor="end">{tick}%</text><line x1="{left}" x2="{left + width}" y1="{y:.1f}" y2="{y:.1f}"/>')
    for g, members in enumerate(groups):
        mean = {c: sum(vals[c] for _, vals in members) / len(members) for c in categories}
        total = sum(mean.values()) or 1.0
        name = members[0][0] if len(members) == 1 else f"{members[0][0]} … {members[-1][0]} ({len(members)} members, mean)"
        tip = name + ": " + ", ".join(f"{c} {mean[c] / total * 100:.0f}%" for c in categories if mean[c] / total >= 0.005)
        gx = left + g * w
        y = float(y0 + height)
        segs = []
        for i, c in enumerate(categories):
            h = mean[c] / total * height
            if h <= 0:
                continue
            y -= h
            segs.append(f'<rect class="{SERIES[i % len(SERIES)]}" x="{gx:.2f}" y="{y:.1f}" width="{max(w - (0.5 if w >= 3 else 0), 0.5):.2f}" height="{h:.1f}"/>')
        parts.append(f"<g><title>{esc(tip)}</title>{''.join(segs)}</g>")
    for k, (label, vals) in enumerate(odd):
        total = sum(vals.values()) or 1.0
        gx = left + len(groups) * w + gap + k * odd_w
        tip = label + ": " + ", ".join(f"{c} {vals[c] / total * 100:.0f}%" for c in categories if vals[c] / total >= 0.005)
        y = float(y0 + height)
        segs = []
        for i, c in enumerate(categories):
            h = vals[c] / total * height
            if h <= 0:
                continue
            y -= h
            segs.append(f'<rect class="{SERIES[i % len(SERIES)]}" x="{gx:.2f}" y="{y:.1f}" width="{odd_w - 1:.2f}" height="{h:.1f}"/>')
        parts.append(f"<g><title>{esc(tip)}</title>{''.join(segs)}"
                     f'<rect class="bad" x="{gx:.2f}" y="{y0 + height + 2}" width="{odd_w - 1:.2f}" height="5"/></g>')
    return _svg(left + width + 10, y0 + height + 12, "".join(parts))


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


def lines(series: list[tuple[str, str, list[tuple[float, float]]]], *, width: int = 620, height: int = 240, x_label: str = "",
          y_label: str = "", percent: bool = False) -> str:
    """Line chart of one or more series over a shared numeric x: (label, line class, [(x, y)]). With `percent`
    the y values are shares of 1 and the axis runs from 0 to 100 %. The legend is drawn above the plot."""
    pts = [(finite(x), finite(y)) for _, _, s in series for x, y in s]
    if not pts:
        return "<p class='muted'>No data.</p>"
    left, bottom, top, right = 46, 34, (40 if y_label else 22), 12   # the y label, then the legend, then the plot
    plot_w, plot_h = width - left - right, height - top - bottom
    lo_x, hi_x = min(x for x, _ in pts), max(x for x, _ in pts)
    if hi_x <= lo_x:
        hi_x = lo_x + 1.0
    hi_y = 1.0 if percent else max(y for _, y in pts) or 1.0
    y_ticks = [0.0, 0.25, 0.5, 0.75, 1.0] if percent else nice_ticks(0, hi_y, 4)
    hi_y = max(hi_y, y_ticks[-1])

    def sx(v: float) -> float:
        return left + (v - lo_x) / (hi_x - lo_x) * plot_w

    def sy(v: float) -> float:
        return top + plot_h - v / hi_y * plot_h

    parts = []
    for t in y_ticks:
        parts.append(f'<line x1="{left}" y1="{sy(t):.1f}" x2="{left + plot_w}" y2="{sy(t):.1f}"/>'
                     f'<text x="{left - 4}" y="{sy(t) + 3:.1f}" text-anchor="end">{esc(f"{t * 100:.0f}%" if percent else fmt(t))}</text>')
    whole = all(float(x).is_integer() for x, _ in pts)   # counts on the x axis: no tick between two of them
    seen: set[str] = set()
    for t in nice_ticks(lo_x, hi_x, 6):
        if lo_x <= t <= hi_x and not (whole and not float(t).is_integer()) and fmt(t) not in seen:
            seen.add(fmt(t))
            parts.append(f'<text x="{sx(t):.1f}" y="{top + plot_h + 12}" text-anchor="middle">{esc(fmt(t))}</text>')
    x = float(left)
    for label, cls, s in series:
        if not s:
            continue
        step = max(1, len(s) // 600)   # at most about 600 points per line; the last one always
        keep = s[::step] + ([s[-1]] if (len(s) - 1) % step else [])
        parts.append(f'<polyline class="{cls}" points="{" ".join(f"{sx(finite(a)):.1f},{sy(finite(b)):.1f}" for a, b in keep)}">'
                     f'<title>{esc(label)}</title></polyline>')
        parts.append(f'<line class="{cls}" x1="{x:.0f}" y1="{top - 10}" x2="{x + 18:.0f}" y2="{top - 10}"/>'
                     f'<text x="{x + 22:.0f}" y="{top - 7}">{esc(label)}</text>')
        x += 22 + 6.2 * len(label) + 16
    if x_label:
        parts.append(f'<text x="{left + plot_w / 2:.1f}" y="{height - 6}" text-anchor="middle">{esc(x_label)}</text>')
    if y_label:
        parts.append(f'<text x="{left}" y="10" font-size="9">{esc(y_label)}</text>')
    return _svg(width, height, "".join(parts))


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


def scatter(background: list[tuple[str, float, float]], points: list[tuple[str, float, float, str]], *,
            width: int = 460, height: int = 320, x_label: str = "PC1", y_label: str = "PC2") -> str:
    """A PCA plane: the background set as small points coloured by label (SERIES classes), the members as larger
    points with their id and label on hover. background: (label, x, y); points: (id, x, y, label)."""
    allx = [x for _, x, _ in background] + [x for _, x, _, _ in points]
    ally = [y for _, _, y in background] + [y for _, _, y, _ in points]
    if not allx:
        return "<p class='muted'>No data.</p>"
    left, bottom, top, right = 44, 28, 8, 8
    plot_w, plot_h = width - left - right, height - top - bottom
    lo_x, hi_x, lo_y, hi_y = min(allx), max(allx), min(ally), max(ally)
    if hi_x <= lo_x:
        lo_x, hi_x = lo_x - 1.0, hi_x + 1.0
    if hi_y <= lo_y:
        lo_y, hi_y = lo_y - 1.0, hi_y + 1.0
    pad_x, pad_y = (hi_x - lo_x) * 0.05, (hi_y - lo_y) * 0.05
    lo_x, hi_x, lo_y, hi_y = lo_x - pad_x, hi_x + pad_x, lo_y - pad_y, hi_y + pad_y

    def sx(v: float) -> float:
        return left + (v - lo_x) / (hi_x - lo_x) * plot_w

    def sy(v: float) -> float:
        return top + (hi_y - v) / (hi_y - lo_y) * plot_h

    labels = sorted({lab for lab, _, _ in background} | {lab for _, _, _, lab in points if lab})
    cls_of = {lab: SERIES[i % len(SERIES)] for i, lab in enumerate(labels)}
    parts = [f'<rect class="plot" x="{left}" y="{top}" width="{plot_w}" height="{plot_h}" fill="none"/>',
             f'<line class="axis" x1="{left}" y1="{top + plot_h}" x2="{left + plot_w}" y2="{top + plot_h}"/>',
             f'<line class="axis" x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_h}"/>',
             f'<text x="{left + plot_w / 2:.1f}" y="{height - 8}" text-anchor="middle">{esc(x_label)}</text>',
             f'<text x="12" y="{top + plot_h / 2:.1f}" text-anchor="middle" transform="rotate(-90 12 {top + plot_h / 2:.1f})">{esc(y_label)}</text>']
    for lab, x, y in background:
        parts.append(f'<circle class="{cls_of.get(lab, "m1")} bg" cx="{sx(finite(x)):.1f}" cy="{sy(finite(y)):.1f}" r="2" opacity="0.35"/>')
    for sid, x, y, lab in points:
        parts.append(f'<circle class="pt {cls_of.get(lab, "")}" cx="{sx(finite(x)):.1f}" cy="{sy(finite(y)):.1f}" r="4.5" stroke="#000" stroke-width="0.8">'
                     f'<title>{esc(sid)}{(": " + esc(lab)) if lab else ""}</title></circle>')
    return _svg(width, height, "".join(parts)) + legend([(lab, cls_of[lab]) for lab in labels])
