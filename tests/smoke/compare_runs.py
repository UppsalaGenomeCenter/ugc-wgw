#!/usr/bin/env python3
"""Two routes over the same members, side by side: what the per-sample digests say after the tools changed.

  compare_runs.py --a-label production --b-label pblab --digests-a DIR --digests-b DIR [--pairs FILE] --out compare.tsv ID...

The digests are `<id>.summary.json` of `ugc-wgw summary` for each project. --pairs is a table the harness computes
from the call sets themselves (member, what, shared, only_a, only_b). Writes one row per member and metric (a, b,
delta) and prints the metrics as rows with one column pair per member. Standard library only: the same script reads
two digest directories brought back from a cluster. Used by tests/smoke/run.sh eval (docs/guide/12-resources.md)."""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path


def dig(d: object, *path: object) -> object:
    for key in path:
        if isinstance(key, tuple):   # ("chrom", "chr20") / ("chromosome", "chr20"): the list entry with that value
            field, value = key
            d = next((x for x in d if isinstance(x, dict) and x.get(field) == value), None) if isinstance(d, list) else None
        else:
            d = d.get(key) if isinstance(d, dict) else None
        if d is None:
            return None
    return d


def named_loci(d: dict) -> str:
    rows = [x for x in dig(d, "tandem_repeats", "disease_loci") or [] if x.get("called")]
    return ";".join(f"{x['trid']}={x.get('al')}" for x in sorted(rows, key=lambda x: x["trid"]))


def starphase(d: dict) -> str:
    return ";".join(sorted(f"{x.get('gene')}={x.get('diplotype')}" for x in dig(d, "targeted", "starphase") or []))


def flags(d: dict) -> str:
    return ",".join(sorted(str(f.get("code")) for f in dig(d, "qc", "flags") or []))


def sv_total(d: dict) -> object:
    counts = dig(d, "sv", "counts")
    return sum(counts.values()) if isinstance(counts, dict) else None


# (metric, getter, digits or None for text)
METRICS = [
    ("mapped reads", lambda d: dig(d, "reads", "mapped_read_count"), 0),
    ("mapped %", lambda d: dig(d, "reads", "mapped_read_percent"), 2),
    ("identity mean %", lambda d: dig(d, "reads", "gap_compressed_identity_mean"), 2),
    ("depth chr20", lambda d: dig(d, "coverage", "chrom", ("chrom", "chr20"), "mean"), 2),
    ("inferred sex", lambda d: dig(d, "coverage", "inferred_sex"), None),
    ("SNVs", lambda d: dig(d, "small_variants", "snv"), 0),
    ("indels", lambda d: dig(d, "small_variants", "indel"), 0),
    ("Ts/Tv", lambda d: dig(d, "small_variants", "tstv"), 3),
    ("het/hom", lambda d: dig(d, "small_variants", "hethom"), 3),
    ("SVs", sv_total, 0),
    ("phased bp", lambda d: dig(d, "phasing", "phased_bp"), 0),
    ("phase blocks chr20", lambda d: dig(d, "phasing", "chrom", ("chromosome", "chr20"), "num_blocks"), 0),
    ("block NG50 chr20", lambda d: dig(d, "phasing", "chrom", ("chromosome", "chr20"), "block_ng50"), 0),
    ("phased hets chr20", lambda d: dig(d, "phasing", "chrom", ("chromosome", "chr20"), "num_phased"), 0),
    ("TRGT loci genotyped", lambda d: dig(d, "tandem_repeats", "genotyped"), 0),
    ("named loci AL", named_loci, None),
    ("CpG sites combined", lambda d: dig(d, "methylation", "cpg", "combined"), 0),
    ("CpG sites hap1", lambda d: dig(d, "methylation", "cpg", "hap1"), 0),
    ("CpG sites hap2", lambda d: dig(d, "methylation", "cpg", "hap2"), 0),
    ("MethBat methylated", lambda d: dig(d, "methylation", "methbat", "methylated"), 0),
    ("MethBat unmethylated", lambda d: dig(d, "methylation", "methbat", "unmethylated"), 0),
    ("MethBat ASM", lambda d: dig(d, "methylation", "methbat", "asm"), 0),
    ("Paraphase regions with CN", lambda d: dig(d, "targeted", "paraphase", "n_with_cn"), 0),
    ("StarPhase diplotypes", starphase, None),
    ("mt haplogroup", lambda d: dig(d, "targeted", "mitorsaw", "haplogroup"), None),
    ("QC status", lambda d: dig(d, "qc", "status"), None),
    ("QC flags", flags, None),
]


def show(value: object, digits: int | None) -> str:
    if value is None or value == "":
        return "NA"
    if digits is None or not isinstance(value, (int, float)):
        return str(value)
    return f"{value:.{digits}f}"


def compare(a: dict, b: dict) -> list[tuple[str, str, str, str]]:
    """(metric, a, b, delta): delta is b - a for numbers, `same` or `differs` for text, NA when either is missing."""
    rows = []
    for name, getter, digits in METRICS:
        va, vb = getter(a), getter(b)
        if digits is None or not isinstance(va, (int, float)) or not isinstance(vb, (int, float)):
            delta = "NA" if va in (None, "") or vb in (None, "") else ("same" if va == vb else "differs")
        else:
            delta = f"{vb - va:+.{digits}f}"
        rows.append((name, show(va, digits), show(vb, digits), delta))
    return rows


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--a-label", default="a")
    ap.add_argument("--b-label", default="b")
    ap.add_argument("--digests-a", required=True, type=Path)
    ap.add_argument("--digests-b", required=True, type=Path)
    ap.add_argument("--pairs", type=Path, help="member, what, shared, only_a, only_b (tab-separated, with a header)")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("members", nargs="+")
    args = ap.parse_args(argv)
    table: dict[str, dict[str, tuple[str, str, str]]] = {}
    order: list[str] = []
    done: list[str] = []
    with args.out.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["member", "metric", args.a_label, args.b_label, "delta"])
        for sid in args.members:
            pa, pb = args.digests_a / f"{sid}.summary.json", args.digests_b / f"{sid}.summary.json"
            if not pa.exists() or not pb.exists():
                print(f"compare: {sid} has no digest in both ({pa.exists()=}, {pb.exists()=}); skipped", file=sys.stderr)
                continue
            done.append(sid)
            for name, va, vb, delta in compare(json.loads(pa.read_text()), json.loads(pb.read_text())):
                w.writerow([sid, name, va, vb, delta])
                if name not in table:
                    table[name] = {}
                    order.append(name)
                table[name][sid] = (va, vb, delta)
        if args.pairs and args.pairs.exists():
            for r in csv.DictReader(args.pairs.open(), delimiter="\t"):
                if r["member"] not in done:
                    continue
                shared, oa, ob = int(r["shared"]), int(r["only_a"]), int(r["only_b"])
                total = shared + oa + ob
                name = f"{r['what']}: shared / only {args.a_label} / only {args.b_label}"
                cell = (f"{shared}", f"{oa} / {ob}", f"{100 * shared / total:.2f}% shared" if total else "NA")
                w.writerow([r["member"], name, shared, f"{oa}/{ob}", cell[2]])
                if name not in table:
                    table[name] = {}
                    order.append(name)
                table[name][r["member"]] = cell
    if not done:
        print("compare: no member with a digest in both directories", file=sys.stderr)
        return 1
    width = max(len(n) for n in order)
    print(f"{'metric':{width}}  " + "  ".join(f"{sid + ' ' + args.a_label + ' > ' + args.b_label:>34}" for sid in done))
    for name in order:
        cells = []
        for sid in done:
            va, vb, delta = table[name].get(sid, ("NA", "NA", "NA"))
            text = delta if len(va) + len(vb) > 26 else f"{va} > {vb}" + ("" if delta in ("same", "NA") or va == vb else f" ({delta})")
            cells.append(f"{text:>34}")
        print(f"{name:{width}}  " + "  ".join(cells))
    return 0


if __name__ == "__main__":
    sys.exit(main())
