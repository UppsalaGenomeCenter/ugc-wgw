#!/usr/bin/env python3
"""Generate references/methbat_map.<ref>.tsv: every CpG island of MethBat's region file with the genes near it
and a region type, so that cohort_methyl's per-region rows and `<cohort>.methyl.events.tsv` can be filtered by
gene (docs/DESIGN.md §13; the methbat-annotation design note).

Inputs
  --regions     MethBat's region TSV (the data container's cpgIslandExt.sorted.hg38.tsv): chrom, start, end,
                cpg_label|region_label; 0-based half-open coordinates, unique rows.
  --gtf         a GENCODE GTF (plain or .gz) at a pinned release; only `gene` rows are read.
  --imprinting  MethBat's report regions (hg38_imprinting_targets.tsv): chrom, start, end, ...

Rules (the generator's version string changes when they do)
  genes        --genes hgnc (default): every gene row with an `hgnc_id` attribute (41,987 of GENCODE 50's 78,733;
               the ids analysts filter by); --genes named: every row with a `gene_name` whose gene_type is not a
               pseudogene (lncRNAs without an HGNC name included).
  TSS          the gene's start on + or its end on - (gene-level; transcripts are not read).
  genes near   a gene whose TSS lies within --tss-window bp (10,000: PacBio's rule for reviewing rare events) of
               the island, or whose body overlaps it; `;`-joined, nearest TSS first, then symbol. Islands
               without a gene keep their row with empty gene columns.
  region_type  the first that applies: imprinted (the island overlaps an imprinting region), promoter (a TSS
               within --promoter bp, 1,000, of the island), gene_body (the island overlaps a gene), intergenic.

Output (stdout): `##` header lines (generator, rules, inputs with sha256), then
  chrom  start  end  region_label  hgnc_id  gene_symbol  region_type

Usage:
  scripts/make-methbat-map.py --regions cpgIslandExt.sorted.hg38.tsv --gtf gencode.v50.annotation.gtf.gz \\
      --imprinting references/methbat_report_regions.hg38.tsv > references/methbat_map.GRCh38_GIABv3.tsv
"""
from __future__ import annotations

import argparse
import bisect
import gzip
import hashlib
import io
import re
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

GENERATOR_VERSION = "1.0.0"
ATTR = re.compile(r'(\S+) "([^"]*)"')


@dataclass
class Gene:
    start: int   # 0-based
    end: int     # exclusive
    strand: str
    hgnc_id: str
    symbol: str

    @property
    def tss(self) -> int:
        return self.start if self.strand == "+" else self.end - 1


def _open(path: Path) -> io.TextIOBase:
    with open(path, "rb") as probe:
        magic = probe.read(2)
    if magic == b"\x1f\x8b":
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8")   # type: ignore[arg-type]
    return open(path, encoding="utf-8")


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_genes(path: Path, selection: str) -> dict[str, list[Gene]]:
    """Gene rows of the GTF per chromosome, sorted by start."""
    genes: dict[str, list[Gene]] = defaultdict(list)
    with _open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "gene":
                continue
            attrs = dict(ATTR.findall(cols[8]))
            hgnc = attrs.get("hgnc_id", "")
            symbol = attrs.get("gene_name", "")
            if selection == "hgnc":
                if not hgnc:
                    continue
            else:   # named
                if not symbol or "pseudogene" in attrs.get("gene_type", ""):
                    continue
            genes[cols[0]].append(Gene(int(cols[3]) - 1, int(cols[4]), cols[6], hgnc, symbol or hgnc))
    for rows in genes.values():
        rows.sort(key=lambda g: (g.start, g.end))
    return genes


def read_regions(path: Path) -> tuple[list[tuple[str, int, int, str]], str]:
    """(rows, label column name) of a MethBat region TSV."""
    rows: list[tuple[str, int, int, str]] = []
    with open(path, encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for required in ("chrom", "start", "end"):
            if required not in idx:
                sys.exit(f"{path}: column {required!r} missing (header: {', '.join(header)})")
        label_col = "region_label" if "region_label" in idx else "cpg_label" if "cpg_label" in idx else None
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            rows.append((c[idx["chrom"]], int(c[idx["start"]]), int(c[idx["end"]]), c[idx[label_col]] if label_col else ""))
    return rows, label_col or ""


class ChromIndex:
    """Genes of one chromosome: TSS bisection for the window query, a prefix maximum of ends for the overlap query."""

    def __init__(self, genes: list[Gene]) -> None:
        self.by_start = genes
        self.starts = [g.start for g in genes]
        self.end_prefix_max: list[int] = []
        m = -1
        for g in genes:
            m = max(m, g.end)
            self.end_prefix_max.append(m)
        self.by_tss = sorted(genes, key=lambda g: g.tss)
        self.tss = [g.tss for g in self.by_tss]

    def tss_within(self, start: int, end: int, window: int) -> list[Gene]:
        lo = bisect.bisect_left(self.tss, start - window)
        hi = bisect.bisect_left(self.tss, end + window)
        return self.by_tss[lo:hi]

    def overlapping(self, start: int, end: int) -> list[Gene]:
        i = bisect.bisect_left(self.starts, end) - 1   # genes starting before the island ends
        hits: list[Gene] = []
        while i >= 0 and self.end_prefix_max[i] > start:
            if self.by_start[i].end > start:
                hits.append(self.by_start[i])
            i -= 1
        return hits


def tss_distance(g: Gene, start: int, end: int) -> int:
    t = g.tss
    if t < start:
        return start - t
    if t >= end:
        return t - end + 1
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--regions", required=True, type=Path, help="MethBat region TSV (CpG islands)")
    ap.add_argument("--gtf", required=True, type=Path, help="GENCODE GTF, plain or gzipped")
    ap.add_argument("--imprinting", required=True, type=Path, help="MethBat report regions TSV (imprinting targets)")
    ap.add_argument("--tss-window", type=int, default=10000, help="a gene is near an island when its TSS is within this many bp (default 10000)")
    ap.add_argument("--promoter", type=int, default=1000, help="region_type promoter when a TSS is within this many bp (default 1000)")
    ap.add_argument("--genes", choices=("hgnc", "named"), default="hgnc", help="which gene rows count (default hgnc: rows with an hgnc_id)")
    ap.add_argument("--gtf-label", default=None, help="how to name the GTF in the header (default: the file name)")
    args = ap.parse_args()

    genes = read_genes(args.gtf, args.genes)
    index = {chrom: ChromIndex(rows) for chrom, rows in genes.items()}
    regions, label_col = read_regions(args.regions)
    imprinting, _ = read_regions(args.imprinting)
    imprinted_by_chrom: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for chrom, start, end, _label in imprinting:
        imprinted_by_chrom[chrom].append((start, end))

    out = sys.stdout
    out.write(f"##methbat_map_generator=scripts/make-methbat-map.py {GENERATOR_VERSION}\n")
    out.write(f"##methbat_map_rules=genes={args.genes}; tss_window={args.tss_window}; promoter={args.promoter}; "
              "region_type=imprinted>promoter>gene_body>intergenic; TSS=gene start (+) or end (-); near=TSS within the window or body overlap\n")
    for role, path in (("regions", args.regions), ("gtf", args.gtf), ("imprinting", args.imprinting)):
        name = args.gtf_label if (role == "gtf" and args.gtf_label) else path.name
        out.write(f"##methbat_map_input={role}\t{name}\tsha256:{sha256_file(path)}\n")
    out.write("chrom\tstart\tend\tregion_label\thgnc_id\tgene_symbol\tregion_type\n")

    types: Counter[str] = Counter()
    with_genes = 0
    for chrom, start, end, label in regions:
        ci = index.get(chrom)
        near: dict[tuple[str, str], Gene] = {}
        if ci is not None:
            for g in ci.tss_within(start, end, args.tss_window):
                near[(g.hgnc_id, g.symbol)] = g
            overlapping = ci.overlapping(start, end)
            for g in overlapping:
                near[(g.hgnc_id, g.symbol)] = g
        else:
            overlapping = []
        ordered = sorted(near.values(), key=lambda g: (tss_distance(g, start, end), g.symbol))
        if any(s < end and e > start for s, e in imprinted_by_chrom.get(chrom, ())):
            rtype = "imprinted"
        elif any(tss_distance(g, start, end) <= args.promoter for g in ordered):
            rtype = "promoter"
        elif overlapping:
            rtype = "gene_body"
        else:
            rtype = "intergenic"
        types[rtype] += 1
        if ordered:
            with_genes += 1
        out.write("\t".join([chrom, str(start), str(end), label,
                             ";".join(g.hgnc_id for g in ordered), ";".join(g.symbol for g in ordered), rtype]) + "\n")
    n_genes = sum(len(v) for v in genes.values())
    print(f"[make-methbat-map] {len(regions)} islands (label column {label_col or 'none'}), {n_genes} genes ({args.genes}), "
          f"{with_genes} islands with a gene; region_type: " + ", ".join(f"{k} {types[k]}" for k in ("imprinted", "promoter", "gene_body", "intergenic")),
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
