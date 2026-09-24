#!/usr/bin/env python3
"""Convert a per-file sequencing manifest into the driver's sample sheet (guide chapter 05).

Input: one or more CSV or TSV files (delimiter detected from the header line, header required) with
one row per BAM file, as a sequencing facility or a PLINK-style pedigree export writes them:

    project,sample,file,family_id,paternal_id,maternal_id,sex,phenotype
    study_plates1-8,study_0001,/proj/data/.../m84045_250210_121343_s1.hifi_reads.bc2001.bam,study_0001,0,0,2,0

Output: one row per sample, tab-separated, with the driver's columns first (`sample_id`, `sex`,
`hifi_reads`, `fail_reads`, `father_id`, `mother_id`) and every other manifest column kept as
metadata (`family_id`, `phenotype`, `project`, ...). Rules:

  - column names are matched case-insensitively; `sample`/`sample_id`, `file`/`path`/`bam`/
    `hifi_reads`, `paternal_id`/`father_id`, `maternal_id`/`mother_id` are accepted
  - the files of a sample are joined with commas in manifest order; a file whose name contains
    `.fail_reads.` goes to `fail_reads`, everything else to `hifi_reads`; duplicates are dropped
  - sex: PLINK coding 1 = MALE, 2 = FEMALE, 0 = unknown (blank); M/F/MALE/FEMALE also accepted
  - parents: `0`, `NA` and empty mean none; a parent that is not a sample of the manifest is a
    warning (the driver warns again at `samples add`)
  - sex, parents and every metadata column must agree across the rows of one sample
  - `--check` verifies that every file exists (off by default: manifests usually name HPC paths)

Exit 1 on any error, with every problem listed. Then: ugc-wgw samples add samples.tsv (chapter 05).
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import OrderedDict
from pathlib import Path

DRIVER_COLUMNS = ("sample_id", "sex", "hifi_reads", "fail_reads", "father_id", "mother_id")
ALIASES = {
    "sample_id": ("sample_id", "sample", "sampleid", "sample_name", "iid"),
    "file": ("file", "path", "bam", "hifi_reads", "reads", "filename"),
    "sex": ("sex", "gender"),
    "father_id": ("father_id", "paternal_id", "father", "pat", "pid"),
    "mother_id": ("mother_id", "maternal_id", "mother", "mat", "mid"),
}
SEX = {"1": "MALE", "2": "FEMALE", "0": "", "": "", "M": "MALE", "F": "FEMALE", "MALE": "MALE", "FEMALE": "FEMALE",
       "NA": "", "UNKNOWN": "", "U": ""}
NO_PARENT = {"", "0", "NA", "NONE", "-"}


def detect_delimiter(path: Path) -> str:
    with open(path, newline="") as fh:
        head = fh.readline()
    return "\t" if head.count("\t") >= head.count(",") and "\t" in head else ","


def resolve_columns(header: list[str], path: Path) -> dict[str, str]:
    """Manifest column name per role; the sample id and the file are required."""
    lower = {h.strip().lower(): h for h in header if h}
    found: dict[str, str] = {}
    for role, names in ALIASES.items():
        for name in names:
            if name in lower:
                found[role] = lower[name]
                break
    for role in ("sample_id", "file"):
        if role not in found:
            sys.exit(f"{path}: no column for the {role.replace('_', ' ')} (accepted: {', '.join(ALIASES[role])}); "
                     f"header is {', '.join(header)}")
    return found


def convert(paths: list[Path], project: str | None, check: bool) -> tuple[list[str], list[dict[str, str]], list[str]]:
    """(output columns, one row per sample, warnings); exits with every hard problem listed."""
    samples: OrderedDict[str, dict[str, object]] = OrderedDict()
    problems: list[str] = []
    warnings: list[str] = []
    meta_columns: list[str] = []
    for path in paths:
        delim = detect_delimiter(path)
        with open(path, newline="") as fh:
            reader = csv.DictReader(fh, delimiter=delim)
            header = [h for h in (reader.fieldnames or []) if h and h.strip()]
            if not header:
                sys.exit(f"{path}: empty file or missing header")
            cols = resolve_columns(header, path)
            role_of = {v: k for k, v in cols.items()}
            for h in header:
                if h not in role_of and h not in meta_columns:
                    meta_columns.append(h)
            for n, row in enumerate(reader, start=2):
                where = f"{path}:{n}"
                sid = (row.get(cols["sample_id"]) or "").strip()
                file = (row.get(cols["file"]) or "").strip()
                if not sid and not file:
                    continue
                if not sid or not file:
                    problems.append(f"{where}: sample id and file are both required")
                    continue
                if project is not None and (row.get("project") or "").strip() != project:
                    continue
                sex_raw = (row.get(cols["sex"]) or "").strip().upper() if "sex" in cols else ""
                if sex_raw not in SEX:
                    problems.append(f"{where}: sex {sex_raw!r} is not 1/2/0 or M/F (sample {sid})")
                    sex = None
                else:
                    sex = SEX[sex_raw]
                parents = {}
                for role in ("father_id", "mother_id"):
                    v = (row.get(cols[role]) or "").strip() if role in cols else ""
                    parents[role] = "" if v.upper() in NO_PARENT else v
                meta = {h: (row.get(h) or "").strip() for h in header if h not in role_of}
                rec = samples.get(sid)
                if rec is None:
                    rec = samples[sid] = {"sex": sex, "father_id": parents["father_id"], "mother_id": parents["mother_id"],
                                          "meta": meta, "hifi_reads": [], "fail_reads": [], "first": where}
                else:
                    for key, val in (("sex", sex), ("father_id", parents["father_id"]), ("mother_id", parents["mother_id"])):
                        if rec[key] != val:
                            problems.append(f"{where}: {key} {val!r} differs from {rec[key]!r} at {rec['first']} (sample {sid})")
                    for k, v in meta.items():
                        if rec["meta"].get(k, "") != v:
                            problems.append(f"{where}: {k} {v!r} differs from {rec['meta'].get(k, '')!r} at {rec['first']} (sample {sid})")
                kind = "fail_reads" if ".fail_reads." in Path(file).name else "hifi_reads"
                if file in rec["hifi_reads"] or file in rec["fail_reads"]:
                    warnings.append(f"{where}: duplicate file dropped: {file}")
                    continue
                if check and not Path(file).exists():
                    problems.append(f"{where}: file not found: {file}")
                rec[kind].append(file)
    for sid, rec in samples.items():
        if not rec["hifi_reads"]:
            problems.append(f"{sid}: no hifi_reads file (only fail reads?)")
        for role in ("father_id", "mother_id"):
            parent = rec[role]
            if parent and parent not in samples:
                warnings.append(f"{sid}: {role} {parent} is not a sample of the manifest")
            if parent == sid:
                problems.append(f"{sid}: is its own {role}")
    if not samples:
        problems.append("no samples" + (f" of project {project!r}" if project else ""))
    if problems:
        sys.exit("manifest problems:\n  " + "\n  ".join(problems))
    columns = list(DRIVER_COLUMNS) + meta_columns
    rows = []
    for sid, rec in samples.items():
        row = {"sample_id": sid, "sex": rec["sex"] or "", "hifi_reads": ",".join(rec["hifi_reads"]),
               "fail_reads": ",".join(rec["fail_reads"]), "father_id": rec["father_id"], "mother_id": rec["mother_id"]}
        row.update({k: rec["meta"].get(k, "") for k in meta_columns})
        rows.append(row)
    return columns, rows, warnings


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n", 1)[0], epilog="Details: the module docstring.")
    ap.add_argument("manifest", nargs="+", type=Path, help="per-file manifest(s), CSV or TSV with a header")
    ap.add_argument("-o", "--output", type=Path, help="sample sheet to write (default: stdout)")
    ap.add_argument("--project", help="keep only rows whose `project` column equals this")
    ap.add_argument("--check", action="store_true", help="fail when a file does not exist on this machine")
    a = ap.parse_args(argv)
    columns, rows, warnings = convert(a.manifest, a.project, a.check)
    for w in warnings:
        print(f"warning: {w}", file=sys.stderr)
    out = open(a.output, "w", newline="") if a.output else sys.stdout
    try:
        writer = csv.DictWriter(out, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    finally:
        if a.output:
            out.close()
    n_files = sum(len(r["hifi_reads"].split(",")) + (len(r["fail_reads"].split(",")) if r["fail_reads"] else 0) for r in rows)
    n_sex = sum(1 for r in rows if r["sex"])
    n_par = sum(1 for r in rows if r["father_id"] or r["mother_id"])
    print(f"{len(rows)} sample(s), {n_files} file(s), sex known for {n_sex}, parents for {n_par}; "
          f"{len(warnings)} warning(s)" + (f"; wrote {a.output}" if a.output else ""), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
