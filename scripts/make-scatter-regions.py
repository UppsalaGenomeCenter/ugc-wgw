#!/usr/bin/env python3
"""Generate references/scatter_regions.<ref>.tsv from a FASTA index (.fai).

Cohort-stage tasks (GLnexus, svx, bcftools merge) scatter over these shards
unconditionally — see docs/DESIGN.md §7.1. Output is BED-like and tab-separated:

    shard   chrom   start   end          (start 0-based, end exclusive)

Default: one shard per primary contig (chr1..22, X, Y, M / 1..22, X, Y, MT).
--rest-shard NAME groups every other contig (unplaced, unlocalized, decoys) into
one extra shard with one row per contig, so that region-scoped cohort steps see
the same contigs as upstream's unscoped ones. A shard name may therefore span
several rows; consumers group rows by the first column.
--chunk-size splits contigs into fixed-size shards, e.g. for finer GLnexus jobs.

Usage:
  scripts/make-scatter-regions.py GRCh38.fa.fai --rest-shard rest > references/scatter_regions.GRCh38.tsv
  scripts/make-scatter-regions.py GRCh38.fa.fai --chunk-size 50000000 > references/scatter_regions.GRCh38.50Mb.tsv
"""
from __future__ import annotations

import argparse
import re
import sys

PRIMARY = re.compile(r"^(chr)?([1-9]|1[0-9]|2[0-2]|X|Y|M|MT)$")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fai", help="FASTA index (.fai)")
    ap.add_argument("--chunk-size", type=int, default=0, help="split contigs into shards of at most this many bp (0 = whole contig)")
    ap.add_argument("--include", default=None, help="regex of contig names to include (default: primary chromosomes)")
    ap.add_argument("--all-contigs", action="store_true", help="include every contig in the .fai (alts, decoys, unplaced)")
    ap.add_argument("--rest-shard", default=None, metavar="NAME", help="collect contigs not matched by --include into one shard called NAME")
    args = ap.parse_args()

    include = re.compile(args.include) if args.include else PRIMARY

    print("#shard\tchrom\tstart\tend")
    n = 0
    rest: list[tuple[str, int]] = []
    with open(args.fai) as fh:
        for line in fh:
            if not line.strip():
                continue
            chrom, length = line.split("\t")[:2]
            length = int(length)
            if not args.all_contigs and not include.match(chrom):
                if args.rest_shard:
                    rest.append((chrom, length))
                continue
            if args.chunk_size and args.chunk_size < length:
                start = 0
                idx = 1
                while start < length:
                    end = min(start + args.chunk_size, length)
                    print(f"{chrom}_{idx:03d}\t{chrom}\t{start}\t{end}")
                    start, idx, n = end, idx + 1, n + 1
            else:
                print(f"{chrom}\t{chrom}\t0\t{length}")
                n += 1
    if rest:
        for chrom, length in rest:
            print(f"{args.rest_shard}\t{chrom}\t0\t{length}")
        n += 1
        print(f"[scatter-regions] shard {args.rest_shard!r} holds {len(rest)} contigs", file=sys.stderr)
    print(f"[scatter-regions] wrote {n} shards", file=sys.stderr)
    return 0 if n else 1


if __name__ == "__main__":
    sys.exit(main())
