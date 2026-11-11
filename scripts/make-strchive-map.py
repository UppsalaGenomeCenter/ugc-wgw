#!/usr/bin/env python3
"""Generate the two reference tables of cohort_tr (docs/DESIGN.md §13; the cohort-tr design note):

  references/tr_disease_loci.strchive-<tag>.tsv   one row per named locus of the TRGT catalog, with STRchive's
                                                  thresholds and the catalog's motifs that count as pathogenic
  references/tr_named_overlaps.tsv                the ordinary catalog loci that overlap a named locus ("twins":
                                                  the same repeat genotyped a second time under another motif or frame)

Inputs
  --strchive   STRchive's `data/STRchive-loci.json` at a pinned tag (the tag upstream's catalog was built from:
               the data container's trgt/README.md names it). CC BY 4.0, https://strchive.org.
  --catalog    the TRGT catalog BED of the reference data container (`adotto_strchive_<date>.hg38.bed.gz`):
               chrom, start, end, `ID=…;MOTIFS=…;STRUC=…[;INCLUDE_FAIL_READS]`; 0-based half-open.

Rules (the generator's version string changes when they do)
  join         a catalog locus is named when its ID is not `chr<N>_<start>_<end>_<motif>`; every named locus must
               have a STRchive entry with the same id and the reverse, else the generator stops (the catalog and
               the JSON are of different tags).
  thresholds   STRchive's, in repeat copies: benign_min/max, intermediate_min/max, pathogenic_min/max; a null
               stays empty. A locus whose benign_max reaches its pathogenic_min stops the generator unless it is
               in NESTED_RANGES (looked at by hand: tiny counts where the benign range ends where the pathogenic
               one starts).
  motifs       STRchive's pathogenic motifs (reference orientation) are matched to the catalog's motifs on the
               reference strand: exact, a rotation, or a rotation compatible under IUPAC codes (the catalog writes
               GCN). `pathogenic_catalog_motifs` lists the catalog's spellings that matched: the tasks sum TRGT's
               MC over exactly these. `motif_match` is the weakest kind used, `partial:` in front when a
               pathogenic motif matched nothing, `none` when none did (the tasks then fall back to length).
  frame        `frame_differs` = yes when the catalog's interval is not STRchive's GRCh38 one (TRGT's frame is
               wider at composite loci).
  twins        an ordinary locus overlapping a named one by at least one base.
  assessable   `no` for the loci in NOT_COPY_COUNTED, with the reason in `assessable_note`: STRchive's thresholds
               there are not counts of the repeat unit, so the screen reports the alleles and classes nothing. A
               listed id that STRchive no longer has stops the generator.

Output: `##` header lines (generator, rules, attribution, inputs with sha256), then the header row and the rows,
sorted by chromosome order of the catalog and start.

Usage:
  scripts/make-strchive-map.py --strchive STRchive-loci.json --strchive-tag v2.4.4 \\
      --catalog adotto_strchive_20250827.hg38.bed.gz \\
      --out-loci references/tr_disease_loci.strchive-2.4.4.tsv --out-overlaps references/tr_named_overlaps.tsv
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import re
import sys
from pathlib import Path

GENERATOR_VERSION = "1.1.0"
ORDINARY_ID = re.compile(r"^chr[^_]+_\d+_\d+_")
# benign_max >= pathogenic_min in STRchive itself: one or two copies apart, read by hand (2026-11-11, v2.4.4)
NESTED_RANGES = {"HMNR7_VWA1", "pre-MIR7-2_CHNG3"}
# thresholds that do not count the catalog's repeat unit, read by hand from STRchive's `details` (2026-11-11, v2.4.4). On an
# 864-member cohort ADTKD_MUC1 classed 88 members pathogenic and 776 below benign by TRGT's motif counts against these.
NOT_COPY_COUNTED = {
    "ADTKD_MUC1": "the disease variant is one inserted base in a 7-cytosine tract inside one unit of a 60 bp VNTR; STRchive's "
                  "thresholds (benign 7, pathogenic 8) count that tract, not copies of the unit, and TRGT's motif counts cannot "
                  "tell the two units apart",
}
IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "N": "ACGT", "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
         "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG"}
MATCH_ORDER = ("exact", "rotation", "iupac")
LOCI_COLUMNS = ["trid", "gene", "disease_id", "disease", "inheritance", "chrom", "start", "end", "motif_len", "benign_min", "benign_max",
                "intermediate_min", "intermediate_max", "pathogenic_min", "pathogenic_max", "reference_motifs", "pathogenic_motifs",
                "catalog_motifs", "pathogenic_catalog_motifs", "motif_match", "frame_differs", "twin_trids", "locus_structure",
                "assessable", "assessable_note"]


def _open(path: Path) -> io.TextIOBase:
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")  # type: ignore[return-value]
    return open(path)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def rotations(motif: str) -> list[str]:
    return [motif[i:] + motif[:i] for i in range(len(motif))]


def compatible(a: str, b: str) -> bool:
    """Same length and every position shares a base under the IUPAC codes."""
    return len(a) == len(b) and all(set(IUPAC.get(x, x)) & set(IUPAC.get(y, y)) for x, y in zip(a, b))


def match_kind(pathogenic: str, catalog: str) -> str | None:
    """How a STRchive pathogenic motif matches a catalog motif on the reference strand, or None."""
    p, c = pathogenic.upper(), catalog.upper()
    if p == c:
        return "exact"
    if len(p) != len(c):
        return None
    rots = rotations(p)
    if c in rots:
        return "rotation"
    if any(compatible(r, c) for r in rots):
        return "iupac"
    return None


def match_motifs(pathogenic: list[str], catalog: list[str]) -> tuple[list[str], str]:
    """The catalog motifs that count as pathogenic (catalog order) and the locus's motif_match value."""
    matched: list[str] = []
    kinds: set[str] = set()
    unmatched = 0
    for p in pathogenic:
        hits = [(match_kind(p, c), c) for c in catalog]
        hits = [(k, c) for k, c in hits if k]
        if not hits:
            unmatched += 1
            continue
        for k, c in hits:
            kinds.add(k)
            if c not in matched:
                matched.append(c)
    if not matched:
        return [], "none"
    weakest = max(kinds, key=MATCH_ORDER.index)
    return [c for c in catalog if c in matched], ("partial:" if unmatched else "") + weakest


def read_catalog(path: Path) -> tuple[dict[str, tuple[str, int, int, list[str]]], list[tuple[str, int, int, str]], list[str]]:
    """Named loci {trid: (chrom, start, end, motifs)}, the ordinary loci [(chrom, start, end, trid)], the contig order."""
    named: dict[str, tuple[str, int, int, list[str]]] = {}
    ordinary: list[tuple[str, int, int, str]] = []
    contigs: list[str] = []
    with _open(path) as fh:
        for n, line in enumerate(fh, 1):
            if not line.strip() or line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 4:
                raise SystemExit(f"{path}:{n}: fewer than four columns")
            fields = dict(f.split("=", 1) for f in c[3].split(";") if "=" in f)
            trid, motifs = fields.get("ID", ""), fields.get("MOTIFS", "")
            if not trid or not motifs:
                raise SystemExit(f"{path}:{n}: no ID or MOTIFS in {c[3][:80]!r}")
            if c[0] not in contigs:
                contigs.append(c[0])
            if ORDINARY_ID.match(trid):
                ordinary.append((c[0], int(c[1]), int(c[2]), trid))
            else:
                if trid in named:
                    raise SystemExit(f"{path}:{n}: named locus {trid} twice")
                named[trid] = (c[0], int(c[1]), int(c[2]), motifs.split(","))
    return named, ordinary, contigs


def num(v: object) -> str:
    """STRchive's copies as written in the table: empty for null, an integer when integral."""
    if v is None:
        return ""
    f = float(v)  # type: ignore[arg-type]
    return str(int(f)) if f == int(f) else repr(f)


def overlaps(named: dict[str, tuple[str, int, int, list[str]]], ordinary: list[tuple[str, int, int, str]]) -> list[tuple[str, str, int]]:
    """(ordinary_trid, named_trid, overlap_bp) for every ordinary locus sharing at least one base with a named one."""
    by_chrom: dict[str, list[tuple[int, int, str]]] = {}
    for trid, (chrom, start, end, _m) in named.items():
        by_chrom.setdefault(chrom, []).append((start, end, trid))
    out = []
    for chrom, start, end, trid in ordinary:
        for ns, ne, ntrid in by_chrom.get(chrom, ()):
            bp = min(end, ne) - max(start, ns)
            if bp > 0:
                out.append((trid, ntrid, bp))
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--strchive", required=True, type=Path, help="STRchive-loci.json at a pinned tag")
    ap.add_argument("--strchive-tag", required=True, help="the tag of the JSON, for the header (e.g. v2.4.4)")
    ap.add_argument("--catalog", required=True, type=Path, help="the TRGT catalog BED (plain or gzipped)")
    ap.add_argument("--catalog-label", default=None, help="how to name the catalog in the header (default: the file name)")
    ap.add_argument("--out-loci", required=True, type=Path)
    ap.add_argument("--out-overlaps", required=True, type=Path)
    args = ap.parse_args()

    loci = json.loads(args.strchive.read_text())
    by_id = {str(entry["id"]): entry for entry in loci}
    if len(by_id) != len(loci):
        raise SystemExit(f"{args.strchive}: duplicate ids")
    named, ordinary, contigs = read_catalog(args.catalog)
    missing, extra = sorted(set(named) - set(by_id)), sorted(set(by_id) - set(named))
    if missing or extra:
        raise SystemExit(f"the catalog's named loci and STRchive {args.strchive_tag} differ: in the catalog only {missing[:10]}; "
                         f"in STRchive only {extra[:10]} (is this the tag the catalog was built from?)")
    pairs = overlaps(named, ordinary)
    twins: dict[str, list[str]] = {}
    for otrid, ntrid, _bp in pairs:
        twins.setdefault(ntrid, []).append(otrid)

    rows = []
    notes = []
    for trid, (chrom, start, end, motifs) in named.items():
        e = by_id[trid]
        bmax, pmin = e.get("benign_max"), e.get("pathogenic_min")
        if bmax is not None and pmin is not None and float(bmax) >= float(pmin) and trid not in NESTED_RANGES:
            raise SystemExit(f"{trid}: benign_max {bmax} reaches pathogenic_min {pmin}; look at the locus and add it to NESTED_RANGES "
                             f"or fix the classification rule (docs/guide/08-stages.md, cohort_tr)")
        pathogenic = [str(m) for m in e.get("pathogenic_motif_reference_orientation") or []]
        pcm, kind = match_motifs(pathogenic, motifs)
        if kind != "exact":
            notes.append(f"{trid}: motif_match {kind}")
        frame = "no" if (e.get("chrom"), e.get("start_hg38"), e.get("stop_hg38")) == (chrom, start, end) else "yes"
        structure = e.get("locus_structure")
        rows.append({
            "trid": trid, "gene": e.get("gene") or "", "disease_id": e.get("disease_id") or "", "disease": e.get("disease") or "",
            "inheritance": "/".join(e.get("inheritance") or []), "chrom": chrom, "start": start, "end": end,
            "motif_len": num(e.get("motif_len")), "benign_min": num(e.get("benign_min")), "benign_max": num(bmax),
            "intermediate_min": num(e.get("intermediate_min")), "intermediate_max": num(e.get("intermediate_max")),
            "pathogenic_min": num(pmin), "pathogenic_max": num(e.get("pathogenic_max")),
            "reference_motifs": ",".join(e.get("reference_motif_reference_orientation") or []),
            "pathogenic_motifs": ",".join(pathogenic), "catalog_motifs": ",".join(motifs),
            "pathogenic_catalog_motifs": ",".join(pcm), "motif_match": kind, "frame_differs": frame,
            "twin_trids": ",".join(twins.get(trid, [])),
            "locus_structure": structure if isinstance(structure, str) else "".join(structure or []),
            "assessable": "no" if trid in NOT_COPY_COUNTED else "yes", "assessable_note": NOT_COPY_COUNTED.get(trid, ""),
        })
    stale = sorted(set(NOT_COPY_COUNTED) - {str(r["trid"]) for r in rows})
    if stale and len(rows) > 20:   # the real table; the unit tests' five-locus fixtures need not carry the listed loci
        raise SystemExit(f"NOT_COPY_COUNTED names loci the inputs do not have: {', '.join(stale)}; look at the list")
    order = {c: i for i, c in enumerate(contigs)}
    rows.sort(key=lambda r: (order[str(r["chrom"])], int(r["start"])))  # type: ignore[call-overload]

    catalog_label = args.catalog_label or args.catalog.name
    header = [
        f"##tr_map_generator=scripts/make-strchive-map.py {GENERATOR_VERSION}",
        f"##tr_map_source=STRchive {args.strchive_tag} (https://strchive.org, https://github.com/dashnowlab/STRchive), licensed CC BY 4.0; "
        "thresholds are STRchive's, in repeat copies; this table is a derived subset",
        "##tr_map_rules=join by locus id; motifs matched on the reference strand (exact, rotation, IUPAC-compatible rotation); "
        f"nested ranges accepted: {', '.join(sorted(NESTED_RANGES))}; not assessable by copy count: {', '.join(sorted(NOT_COPY_COUNTED))}",
        f"##tr_map_input=strchive\t{args.strchive.name}@{args.strchive_tag}\tsha256:{sha256_file(args.strchive)}",
        f"##tr_map_input=catalog\t{catalog_label}\tsha256:{sha256_file(args.catalog)}",
    ]
    with open(args.out_loci, "w") as out:
        out.write("\n".join(header) + "\n")
        out.write("\t".join(LOCI_COLUMNS) + "\n")
        for r in rows:
            out.write("\t".join(str(r[c]).replace("\t", " ").replace("\n", " ") for c in LOCI_COLUMNS) + "\n")
    start_of = {t: (order[c], s) for c, s, _e, t in ordinary}
    with open(args.out_overlaps, "w") as out:
        out.write("\n".join(header) + "\n")
        out.write("ordinary_trid\tnamed_trid\toverlap_bp\n")
        for otrid, ntrid, bp in sorted(pairs, key=lambda p: (start_of[p[0]], p[1])):
            out.write(f"{otrid}\t{ntrid}\t{bp}\n")
    print(f"{len(rows)} named loci ({sum(1 for r in rows if r['twin_trids'])} with a twin), {len(pairs)} overlap pairs; "
          f"{len(notes)} motif note(s): {'; '.join(notes) or 'none'}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
