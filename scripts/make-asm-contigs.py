#!/usr/bin/env python3
"""Generate the assembly stages' contig lists for one reference build (docs/DESIGN.md §10.4, §13).

  keep   every contig of a FASTA index except chrEBV and the GIABv3 decoys matching
         ^(KMT2C|KCNJ18|MAP2K3)_ : the contigs of the decoy-free `asm` reference that the installer
         derives with `samtools faidx -r` (references.lock ugc-wgw-asm.built). PAV's authors advise
         against ALT, PATCH and decoy contigs; every kept contig is byte-identical to the shipped build.
  graph  chr1..22, chrX, chrY: the chromosomes the cohort pangenome is built per (cactus-graphmap-split
         --refContigs); everything else goes to its "other" contig.

Usage:
  scripts/make-asm-contigs.py keep <ref.fa.fai> --expect 194 > references/asm_keep_contigs.GRCh38_GIABv3.txt
  scripts/make-asm-contigs.py graph > references/asm_graph_chromosomes.GRCh38.txt
"""
from __future__ import annotations

import argparse
import re
import sys

DROP = re.compile(r"^(chrEBV|(KMT2C|KCNJ18|MAP2K3)_.*)$")
GRAPH = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]


def keep(fai: str, expect: int | None) -> list[str]:
    names = [line.split("\t", 1)[0] for line in open(fai) if line.strip()]
    kept = [n for n in names if not DROP.match(n)]
    if expect is not None and len(kept) != expect:
        sys.exit(f"{fai}: kept {len(kept)} of {len(names)} contigs, expected {expect}")
    return kept


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    k = sub.add_parser("keep", help="contigs of the decoy-free assembly reference")
    k.add_argument("fai")
    k.add_argument("--expect", type=int, default=None, help="fail unless exactly this many contigs are kept")
    sub.add_parser("graph", help="the pangenome's chromosomes")
    a = ap.parse_args()
    names = keep(a.fai, a.expect) if a.cmd == "keep" else GRAPH
    sys.stdout.write("".join(n + "\n" for n in names))


if __name__ == "__main__":
    main()
