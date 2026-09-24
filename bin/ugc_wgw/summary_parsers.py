"""Parsers for the result files `ugc-wgw summary` reads (docs/guide/07-results.md, "Analysis summary").

Every parser takes an iterable of lines (text, or bytes for the big files), skips malformed rows and returns
plain JSON-able dicts; JSON files are passed in already loaded. Nothing here opens the state DB or renders HTML:
`summary.py` decides which files to read and turns the results into the per-sample digest. Formats were taken
from the real files (the smoke results of 2026-09-22).
"""
from __future__ import annotations

import gzip
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import IO, Iterable, Iterator

PRIMARY_CONTIGS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
AUTOSOME_RE = re.compile(r"^(chr)?\d{1,2}$")
CHRX_RE = re.compile(r"^(chr)?X$")
CHRY_RE = re.compile(r"^(chr)?Y$")

SV_SIZE_BINS = ["<100", "100-1k", "1k-10k", "10k-100k", "100k-1M", ">1M"]
SV_SIZE_EDGES = [100, 1_000, 10_000, 100_000, 1_000_000]
ROH_BINS = ["<100k", "100k-500k", "500k-1M", "1M-5M", ">5M"]
ROH_EDGES = [100_000, 500_000, 1_000_000, 5_000_000]
BLOCK_BINS = ["<10k", "10k-100k", "100k-1M", "1M-10M", ">10M"]
BLOCK_EDGES = [10_000, 100_000, 1_000_000, 10_000_000]


# ---- reading -----------------------------------------------------------------------


def _is_gzip(path: Path) -> bool:
    with open(path, "rb") as fh:
        return fh.read(2) == b"\x1f\x8b"


def open_text(path: Path) -> IO[str]:
    """Open a text file, gzip-transparently (by magic bytes), decoding errors replaced."""
    if _is_gzip(path):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def open_binary(path: Path) -> IO[bytes]:
    if _is_gzip(path):
        return gzip.open(path, "rb")  # type: ignore[return-value]
    return open(path, "rb")


def read_lines(path: Path) -> Iterator[str]:
    with open_text(path) as fh:
        for line in fh:
            yield line.rstrip("\n")


def read_binary_lines(path: Path) -> Iterator[bytes]:
    with open_binary(path) as fh:
        for line in fh:
            yield line


def load_json(path: Path) -> object:
    with open_text(path) as fh:
        return json.load(fh)


def _float(v: object) -> float | None:
    if isinstance(v, bytes):
        v = v.decode("ascii", "replace")
    try:
        x = float(str(v).strip())
    except (TypeError, ValueError):
        return None
    if x != x or x in (float("inf"), float("-inf")):
        return None
    return x


def _int(v: object) -> int | None:
    if isinstance(v, bytes):
        v = v.decode("ascii", "replace")
    try:
        return int(str(v).strip())
    except (TypeError, ValueError):
        x = _float(v)
        return int(x) if x is not None else None


def _num(v: object) -> object:
    """int, float or the stripped string, for header-driven tables."""
    s = str(v).strip()
    i = _int(s) if re.fullmatch(r"-?\d+", s) else None
    if i is not None:
        return i
    f = _float(s)
    return f if f is not None and re.fullmatch(r"-?\d*\.?\d+(e[-+]?\d+)?", s, re.I) else s


def _bin_label(value: float, edges: list[float], labels: list[str]) -> str:
    for edge, label in zip(edges, labels):
        if value < edge:
            return label
    return labels[-1]


def contig_class(chrom: str) -> str:
    """autosome, chrX, chrY, chrM or other (the rule of workflows/ugc_wgw/cohort/freq.wdl)."""
    c = chrom[3:] if chrom.startswith("chr") else chrom
    if c.isdigit():
        return "autosome"
    if c == "X":
        return "chrX"
    if c == "Y":
        return "chrY"
    if c in ("M", "MT"):
        return "chrM"
    return "other"


def _rows(lines: Iterable[str]) -> Iterator[list[str]]:
    for line in lines:
        s = line.rstrip("\n")
        if not s or s.startswith("#"):
            continue
        yield s.split("\t")


def table(lines: Iterable[str]) -> tuple[list[str], list[dict[str, str]]]:
    """A header-driven TSV: (columns, rows as dicts). `##` lines are skipped, a `#` on the header is stripped."""
    cols: list[str] = []
    rows: list[dict[str, str]] = []
    for line in lines:
        s = line.rstrip("\n")
        if not s or s.startswith("##"):
            continue
        parts = s.split("\t")
        if not cols:
            cols = [c.lstrip("#") for c in parts]
            continue
        if s.startswith("#"):
            continue
        rows.append({c: parts[i] if i < len(parts) else "" for i, c in enumerate(cols)})
    return cols, rows


# ---- per-sample files ----------------------------------------------------------------


def parse_two_line_tsv(lines: Iterable[str]) -> dict[str, str]:
    """stats_file: one header line and one value line (consolidate_stats writes it transposed)."""
    cols, rows = table(lines)
    return dict(rows[0]) if rows else {}


def parse_bcftools_stats(lines: Iterable[str]) -> dict[str, object]:
    """bcftools stats: SN, TSTV, ST, IDD, QUAL, DP and PSC sections (the `# CODE` header lines are comments)."""
    sn: dict[str, int] = {}
    tstv: dict[str, float] = {}
    st: dict[str, int] = {}
    idd: list[list[int]] = []
    qual: list[list[float]] = []
    dp: list[list[object]] = []
    psc: dict[str, object] = {}
    for p in _rows(lines):
        code = p[0]
        try:
            if code == "SN" and len(p) >= 4:
                sn[p[2].rstrip(":").strip()] = int(p[3])
            elif code == "TSTV" and len(p) >= 5:
                tstv = {"ts": int(p[2]), "tv": int(p[3]), "ratio": float(p[4])}
            elif code == "ST" and len(p) >= 4:
                st[p[2]] = int(p[3])
            elif code == "IDD" and len(p) >= 4:
                idd.append([int(p[2]), int(p[3])])
            elif code == "QUAL" and len(p) >= 7:
                qual.append([float(p[2]), int(p[3]), int(p[4]), int(p[5]), int(p[6])])
            elif code == "DP" and len(p) >= 6:
                dp.append([p[2], int(p[3]), int(p[5])])
            elif code == "PSC" and len(p) >= 14:
                psc = {"sample": p[2], "n_ref_hom": int(p[3]), "n_nonref_hom": int(p[4]), "n_hets": int(p[5]),
                       "n_transitions": int(p[6]), "n_transversions": int(p[7]), "n_indels": int(p[8]),
                       "average_depth": float(p[9]), "n_singletons": int(p[10]), "n_missing": int(p[13])}
        except ValueError:
            continue
    return {"sn": sn, "tstv": tstv, "st": st, "idd": idd, "qual": qual, "dp": dp, "psc": psc}


def parse_mosdepth_summary(lines: Iterable[str]) -> dict[str, object]:
    """mosdepth summary: per-contig rows, their `<contig>_region` twins and the `total` rows."""
    chrom: list[dict[str, object]] = []
    region: dict[str, float] = {}
    total: dict[str, object] = {}
    cols, rows = table(lines)
    for r in rows:
        name = r.get("chrom", "")
        mean = _float(r.get("mean"))
        if name.endswith("_region"):
            base = name[: -len("_region")]
            if base != "total" and mean is not None:
                region[base] = mean
            continue
        rec = {"chrom": name, "length": _int(r.get("length")), "bases": _int(r.get("bases")), "mean": mean,
               "min": _float(r.get("min")), "max": _float(r.get("max"))}
        if name == "total":
            total = rec
        else:
            chrom.append(rec)
    return {"chrom": chrom, "region": region, "total": total}


def coverage_ratios(chrom_rows: list[dict[str, object]], y_threshold: float) -> dict[str, object]:
    """chrX and chrY mean depth over the mean of the autosome means, and the sex upstream's rule would infer."""
    auto = [float(r["mean"]) for r in chrom_rows if AUTOSOME_RE.match(str(r.get("chrom", ""))) and r.get("mean") is not None]
    x = [float(r["mean"]) for r in chrom_rows if CHRX_RE.match(str(r.get("chrom", ""))) and r.get("mean") is not None]
    y = [float(r["mean"]) for r in chrom_rows if CHRY_RE.match(str(r.get("chrom", ""))) and r.get("mean") is not None]
    auto_mean = sum(auto) / len(auto) if auto else 0.0
    out: dict[str, object] = {"auto_mean": auto_mean, "chrX_ratio": None, "chrY_ratio": None,
                              "sex_from_ratio": None, "y_threshold": y_threshold}
    if auto_mean > 0:
        if x:
            out["chrX_ratio"] = (sum(x) / len(x)) / auto_mean
        if y:
            yr = (sum(y) / len(y)) / auto_mean
            out["chrY_ratio"] = yr
            out["sex_from_ratio"] = "MALE" if yr > y_threshold else "FEMALE"
    return out


def aggregate_region_bed(lines: Iterable[bytes], bin_bp: int = 1_000_000,
                         contigs: Iterable[str] | None = PRIMARY_CONTIGS) -> dict[str, object]:
    """mosdepth regions bed (500 bp windows): mean depth per `bin_bp` bin per contig, one pass, exact."""
    contig_list = list(contigs) if contigs is not None else None
    keep = set(contig_list) if contig_list is not None else None
    sums: dict[str, dict[int, float]] = defaultdict(dict)
    counts: dict[str, dict[int, int]] = defaultdict(dict)
    ends: dict[str, int] = {}
    for line in lines:
        parts = line.split(b"\t", 3)
        if len(parts) < 4:
            continue
        chrom = parts[0].decode("ascii", "replace")
        if keep is not None and chrom not in keep:
            continue
        try:
            start, end, depth = int(parts[1]), int(parts[2]), float(parts[3])
        except ValueError:
            continue
        b = start // bin_bp
        s = sums[chrom]
        s[b] = s.get(b, 0.0) + depth
        c = counts[chrom]
        c[b] = c.get(b, 0) + 1
        if end > ends.get(chrom, 0):
            ends[chrom] = end
    order = [c for c in contig_list if c in ends] if contig_list is not None else sorted(ends)
    chroms: dict[str, object] = {}
    for chrom in order:
        n = (ends[chrom] + bin_bp - 1) // bin_bp
        s, c = sums[chrom], counts[chrom]
        chroms[chrom] = {"len": ends[chrom],
                         "bins": [round(s[i] / c[i], 3) if c.get(i) else None for i in range(n)]}
    return {"bin_bp": bin_bp, "chroms": chroms}


def parse_copynum_summary(doc: object) -> dict[str, object]:
    """sawfish copy-number summary: haploid coverage and per-contig most common copy number."""
    if not isinstance(doc, dict):
        return {}
    chrom_doc = doc.get("chromosomes")
    rows: list[dict[str, object]] = []
    if isinstance(chrom_doc, dict):
        names = [c for c in PRIMARY_CONTIGS if c in chrom_doc] + sorted(c for c in chrom_doc if c not in PRIMARY_CONTIGS)
        for name in names:
            d = chrom_doc.get(name)
            if not isinstance(d, dict):
                continue
            per_cn = d.get("bases_per_copy_number") if isinstance(d.get("bases_per_copy_number"), dict) else {}
            rows.append({"chrom": name, "most_common_cn": _int(d.get("most_common_copy_number")),
                         "total_bases": _int(d.get("total_copy_number_bases")),
                         "cn_bases": {str(k): _int(v) for k, v in per_cn.items()}})
    return {"haploid_coverage": _float(doc.get("gc_bias_corrected_haploid_coverage")),
            "sample": doc.get("sample_name"), "chrom": rows}


def parse_bedgraph(lines: Iterable[str], contigs: Iterable[str] | None = PRIMARY_CONTIGS) -> dict[str, list[list[float]]]:
    """A 4-column bedgraph (the sawfish copy-number segments) as {contig: [[start, end, value], ...]}."""
    keep = set(contigs) if contigs is not None else None
    out: dict[str, list[list[float]]] = defaultdict(list)
    for p in _rows(lines):
        if len(p) < 4 or (keep is not None and p[0] not in keep):
            continue
        try:
            out[p[0]].append([int(p[1]), int(p[2]), float(p[3])])
        except ValueError:
            continue
    return dict(out)


def parse_roh_out(lines: Iterable[bytes]) -> dict[str, object]:
    """bcftools roh output: the RG (region) lines only; ST (per-site state) lines are skipped."""
    n = 0
    total = 0
    longest = 0
    by_chrom: Counter = Counter()
    hist: Counter = Counter()
    for line in lines:
        if not line.startswith(b"RG\t"):
            continue
        p = line.rstrip(b"\n").split(b"\t")
        if len(p) < 6:
            continue
        try:
            length = int(p[5])
        except ValueError:
            continue
        n += 1
        total += length
        longest = max(longest, length)
        by_chrom[p[2].decode("ascii", "replace")] += length
        hist[_bin_label(length, ROH_EDGES, ROH_BINS)] += 1
    return {"n": n, "total_bp": total, "longest_bp": longest, "by_chrom": dict(by_chrom),
            "hist": [[label, hist.get(label, 0)] for label in ROH_BINS]}


def _info(field: bytes) -> dict[bytes, bytes]:
    out: dict[bytes, bytes] = {}
    for item in field.split(b";"):
        k, _, v = item.partition(b"=")
        out[k] = v
    return out


def parse_sv_vcf(lines: Iterable[bytes]) -> dict[str, object]:
    """A sawfish (or svx-merged) SV VCF: records by SVTYPE, SVLEN size bins, filters, SUPP and samples."""
    counts: Counter = Counter()
    sizes: dict[str, Counter] = defaultdict(Counter)
    filters: Counter = Counter()
    supp: Counter = Counter()
    samples: list[str] = []
    n = 0
    for line in lines:
        if line.startswith(b"#"):
            if line.startswith(b"#CHROM"):
                samples = [s.decode("ascii", "replace") for s in line.rstrip(b"\n").split(b"\t")[9:]]
            continue
        p = line.split(b"\t", 8)
        if len(p) < 8:
            continue
        n += 1
        info = _info(p[7].rstrip(b"\n"))
        svtype = info.get(b"SVTYPE", b"?").decode("ascii", "replace") or "?"
        counts[svtype] += 1
        filters[p[6].decode("ascii", "replace")] += 1
        svlen = info.get(b"SVLEN")
        if svlen is not None:
            try:
                size = abs(int(svlen.split(b",")[0]))
                sizes[svtype][_bin_label(size, SV_SIZE_EDGES, SV_SIZE_BINS)] += 1
            except ValueError:
                pass
        s = info.get(b"SUPP")
        if s is not None:
            try:
                supp[int(s)] += 1
            except ValueError:
                pass
    return {"n": n, "counts": dict(counts), "sizes": {t: dict(c) for t, c in sizes.items()},
            "size_bins": SV_SIZE_BINS, "filters": dict(filters),
            "supp": {str(k): v for k, v in sorted(supp.items())}, "samples": samples}


def parse_phase_stats(lines: Iterable[str]) -> dict[str, object]:
    """HiPhase stats: per-chromosome rows and the `all` roll-up."""
    cols, rows = table(lines)
    out_rows: list[dict[str, object]] = []
    all_row: dict[str, object] = {}
    for r in rows:
        rec = {c: _num(v) for c, v in r.items()}
        if r.get("chromosome") == "all":
            all_row = rec
        else:
            out_rows.append(rec)
    return {"all": all_row, "chrom": out_rows}


def parse_phase_blocks(lines: Iterable[str]) -> dict[str, object]:
    """HiPhase blocks: count, N50 by block length and a length histogram."""
    cols, rows = table(lines)
    lengths: list[int] = []
    for r in rows:
        s, e = _int(r.get("start")), _int(r.get("end"))
        if s is not None and e is not None and e >= s:
            lengths.append(e - s)
    hist: Counter = Counter(_bin_label(x, BLOCK_EDGES, BLOCK_BINS) for x in lengths)
    n50 = 0
    if lengths:
        half = sum(lengths) / 2
        acc = 0
        for x in sorted(lengths, reverse=True):
            acc += x
            if acc >= half:
                n50 = x
                break
    return {"n": len(lengths), "n50": n50, "total_bp": sum(lengths),
            "hist": [[label, hist.get(label, 0)] for label in BLOCK_BINS]}


POSITIONAL_TRID_RE = re.compile(rb"TRID=[^;_\t]+_\d+_\d+_")


def parse_trgt_vcf(lines: Iterable[bytes], fail_reads_trids: Iterable[str] = ()) -> dict[str, object]:
    """TRGT VCF: the named (disease) loci with their genotype fields; positional loci are counted only."""
    fail = set(fail_reads_trids)
    loci: list[dict[str, object]] = []
    n = 0
    for line in lines:
        if line.startswith(b"#"):
            continue
        n += 1
        if b"TRID=chr" in line and POSITIONAL_TRID_RE.search(line):
            continue
        p = line.rstrip(b"\n").split(b"\t")
        if len(p) < 10:
            continue
        info = _info(p[7])
        trid = info.get(b"TRID", b"").decode("ascii", "replace")
        if not trid or POSITIONAL_TRID_RE.match(b"TRID=" + trid.encode()):
            continue
        keys = p[8].decode("ascii", "replace").split(":")
        vals = p[9].decode("ascii", "replace").split(":")
        fmt = {k: (vals[i] if i < len(vals) else ".") for i, k in enumerate(keys)}
        disease, _, gene = trid.partition("_")
        loci.append({"trid": trid, "disease": disease, "gene": gene or trid,
                     "chrom": p[0].decode("ascii", "replace"), "pos": _int(p[1]),
                     "end": _int(info.get(b"END", b"")), "motifs": info.get(b"MOTIFS", b"").decode("ascii", "replace"),
                     "gt": fmt.get("GT", "."), "al": fmt.get("AL", "."), "mc": fmt.get("MC", "."),
                     "sd": fmt.get("SD", "."), "ap": fmt.get("AP", "."), "am": fmt.get("AM", "."),
                     "called": fmt.get("GT", ".") not in (".", "./.", ".|."),
                     "fail_reads": trid in fail})
    loci.sort(key=lambda d: (str(d["gene"]).upper(), str(d["trid"])))
    return {"n_records": n, "n_named": len(loci), "n_named_called": sum(1 for d in loci if d["called"]), "loci": loci}


def parse_trgt_dropouts(lines: Iterable[bytes]) -> dict[str, object]:
    """TRGT coverage dropouts: counts per dropout class, and the class of every named locus."""
    by_class: Counter = Counter()
    named: dict[str, str] = {}
    col = -1
    trid_col = 3
    first = True
    n = 0
    for line in lines:
        p = line.rstrip(b"\n").split(b"\t")
        if first:
            first = False
            if b"dropout" in p:
                col = p.index(b"dropout")
                trid_col = p.index(b"trid") if b"trid" in p else 3
                continue
        if len(p) <= max(col, trid_col, 0) and col >= 0:
            continue
        n += 1
        label = p[col].decode("ascii", "replace")
        by_class[label] += 1
        t = p[trid_col]
        if b"STRUC=<TR>" in t or (t.startswith(b"ID=") and not t.startswith(b"ID=chr")
                                  and not POSITIONAL_TRID_RE.match(b"TRID=" + t[3:])):
            trid = t.split(b";", 1)[0]
            if trid.startswith(b"ID="):
                named[trid[3:].decode("ascii", "replace")] = label
    return {"n": n, "by_class": dict(by_class), "named": named}


def parse_methbat_profile(lines: Iterable[str]) -> dict[str, object]:
    """MethBat profile: summary_label counts and a histogram of the combined methylation per region."""
    cols, rows = table(lines)
    labels: Counter = Counter()
    hist = [0] * 10
    with_value = 0
    for r in rows:
        labels[r.get("summary_label", "") or "?"] += 1
        m = _float(r.get("mean_combined_methyl"))
        if m is not None:
            with_value += 1
            hist[min(int(m // 10), 9)] += 1
    return {"n_regions": len(rows), "n_with_value": with_value, "labels": dict(labels),
            "hist": [[f"{i * 10}-{i * 10 + 10}", hist[i]] for i in range(10)]}


def parse_paraphase(doc: object) -> dict[str, object]:
    """Paraphase summary: one row per region (copy number, haplotypes, coverage failure)."""
    regions: list[dict[str, object]] = []
    if isinstance(doc, dict):
        for name in sorted(doc, key=lambda s: str(s).upper()):
            d = doc[name]
            if not isinstance(d, dict):
                continue
            depth = d.get("region_depth") if isinstance(d.get("region_depth"), dict) else {}
            haps = d.get("final_haplotypes")
            regions.append({"region": str(name), "genes": str(d.get("genes_in_region") or ""),
                            "depth_median": _float(depth.get("median")) if depth else None,
                            "total_cn": _int(d.get("total_cn")) if d.get("total_cn") is not None else None,
                            "n_haplotypes": len(haps) if isinstance(haps, (dict, list)) else 0,
                            "failed": bool(d.get("failed_for_coverage"))})
    return {"n_regions": len(regions), "n_failed": sum(1 for r in regions if r["failed"]),
            "n_with_cn": sum(1 for r in regions if r["total_cn"] is not None), "regions": regions}


def parse_mitorsaw(doc: object) -> dict[str, object]:
    """mitorsaw haplotype statistics."""
    haps: list[dict[str, object]] = []
    stats: dict[str, object] = {}
    if isinstance(doc, dict):
        for h in doc.get("haplotypes") or []:
            if isinstance(h, dict):
                haps.append({"label": str(h.get("label", "")), "seq_len": _int(h.get("seq_len")),
                             "num_ref_variants": _int(h.get("num_ref_variants")),
                             "estimated_abundance": _float(h.get("estimated_abundance"))})
        fs = doc.get("fingerprint_stats")
        if isinstance(fs, dict):
            stats = {"passing_explanation": fs.get("passing_explanation"),
                     "unexplained_fraction": _float(fs.get("unexplained_fraction"))}
    return {"haplotypes": haps, "fingerprint": stats}


def parse_starphase_tsv(lines: Iterable[str]) -> list[dict[str, str]]:
    """StarPhase TSV (`#gene\\tdiplotype`)."""
    cols, rows = table(lines)
    return [{"gene": r.get("gene", ""), "diplotype": r.get("diplotype", "")} for r in rows if r.get("gene")]


def parse_starphase_json(doc: object) -> dict[str, object]:
    """StarPhase JSON: the tool and database versions, and how many genes have a diplotype."""
    if not isinstance(doc, dict):
        return {}
    genes = doc.get("gene_details") if isinstance(doc.get("gene_details"), dict) else {}
    n_called = 0
    for d in genes.values():
        dips = d.get("diplotypes") if isinstance(d, dict) else None
        if dips and all("NO_READS" not in str(x.get("diplotype", "")) for x in dips if isinstance(x, dict)):
            n_called += 1
    meta = doc.get("database_metadata") if isinstance(doc.get("database_metadata"), dict) else {}
    return {"version": doc.get("pbstarphase_version"), "database": {k: str(v) for k, v in meta.items()},
            "n_genes": len(genes), "n_called": n_called}


def parse_kivvi(doc: object) -> dict[str, object]:
    """kivvi JSON: copy number (a number or "NA") and the call status."""
    if not isinstance(doc, dict):
        return {}
    cn = doc.get("allele_cn")
    depth = doc.get("depth_summary") if isinstance(doc.get("depth_summary"), dict) else {}
    extra = doc.get("additional") if isinstance(doc.get("additional"), dict) else {}
    alleles = doc.get("complete_alleles")
    return {"allele_cn": _float(cn) if _float(cn) is not None else (str(cn) if cn is not None else None),
            "genome_depth": _float(depth.get("genome_depth")), "repeat_depth": _float(depth.get("repeat_depth")),
            "n_complete_alleles": len(alleles) if isinstance(alleles, list) else 0,
            "call_status": str(extra.get("call_status")) if extra.get("call_status") is not None else None}


def fail_reads_trids(messages: Iterable[str]) -> list[str]:
    """The loci genotyped with fail reads, from the `INCLUDE_FAIL_READS regions: A,B,C` workflow message."""
    out: list[str] = []
    for m in messages:
        if m.startswith("INCLUDE_FAIL_READS regions:"):
            out.extend(x.strip() for x in m.split(":", 1)[1].split(",") if x.strip())
    return out


# ---- cohort and assembly files ---------------------------------------------------------


def parse_freq_summary(lines: Iterable[str]) -> dict[str, dict[str, dict[str, object]]]:
    """cohort_freq summary (resource, metric, key, value) as nested dicts; numbers converted."""
    out: dict[str, dict[str, dict[str, object]]] = {}
    for p in _rows(lines):
        if len(p) < 4 or p[0] == "resource":
            continue
        out.setdefault(p[0], {}).setdefault(p[1], {})[p[2]] = _num(p[3])
    return out


def parse_freq_samples(lines: Iterable[str]) -> list[dict[str, str]]:
    cols, rows = table(lines)
    return rows


def parse_gfatools_stats(lines: Iterable[str]) -> dict[str, object]:
    """gfatools stat: SZ (total length), NN (sequences), NL rows (x, Nx, Lx), AU (auN)."""
    out: dict[str, object] = {"total_length": None, "n_sequences": None, "aun": None, "nl": [], "n50": None, "l50": None}
    nl: list[list[int]] = []
    for p in _rows(lines):
        if p[0] == "CC":
            continue
        try:
            if p[0] == "SZ":
                out["total_length"] = int(p[1])
            elif p[0] == "NN":
                out["n_sequences"] = int(p[1])
            elif p[0] == "AU":
                out["aun"] = int(float(p[1]))
            elif p[0] == "NL" and len(p) >= 4:
                nl.append([int(p[1]), int(p[2]), int(p[3])])
        except (ValueError, IndexError):
            continue
    out["nl"] = nl
    for x, nx, lx in nl:
        if x == 50:
            out["n50"], out["l50"] = nx, lx
    return out


def tsv_dims(path: Path) -> tuple[list[str], int]:
    """(header columns, data rows) of a TSV, counting newlines in binary chunks; fine for multi-GB files."""
    cols: list[str] = []
    n = 0
    with open_binary(path) as fh:
        first = fh.readline()
        cols = first.rstrip(b"\n").decode("utf-8", "replace").split("\t") if first else []
        while True:
            chunk = fh.read(1 << 20)
            if not chunk:
                break
            n += chunk.count(b"\n")
    return cols, n


def vcf_samples(path: Path) -> list[str]:
    """The sample columns of a VCF, reading only its header."""
    with open_binary(path) as fh:
        for line in fh:
            if line.startswith(b"#CHROM"):
                return [s.decode("ascii", "replace") for s in line.rstrip(b"\n").split(b"\t")[9:]]
            if not line.startswith(b"#"):
                break
    return []


# ---- software versions from the WDL ---------------------------------------------------------

DOCKER_RE = re.compile(r'docker:\s*"(?P<image>[^"]+)"\s*(?:#\s*(?P<comment>.*?))?\s*$')
VERSION_STRING_RE = re.compile(r'String\s+(?P<tool>\w+?)_version\s*=\s*"(?P<version>[^"]+)"')
DIGEST_RE = re.compile(r"@sha256:([0-9a-f]{64})")


def _split_comment(comment: str) -> tuple[str, str]:
    """`1.7.0_build2` -> (1.7.0, build2); `pb_wdl_base:build4` -> (build4, ""); a trailing parenthetical is dropped."""
    c = re.sub(r"\s*\(.*\)\s*$", "", comment).strip()
    if ":" in c and " " not in c:
        c = c.split(":", 1)[1]
    m = re.match(r"^(?P<version>.*?)_(?P<build>build\d+)$", c)
    if m:
        return m.group("version"), m.group("build")
    return c, ""


def scan_tool_versions(code_dir: Path, roots: Iterable[str] = ("vendor/hifi-human-wgs-wdl/workflows", "workflows")) -> list[dict[str, object]]:
    """Tools, versions and image digests from the `docker:` pins (and `String <tool>_version`) in the WDL files."""
    found: dict[tuple[str, str], dict[str, object]] = {}
    for root in roots:
        base = code_dir / root
        if not base.is_dir():
            continue
        for path in sorted(base.rglob("*.wdl")):
            rel = str(path.relative_to(code_dir))
            try:
                text = path.read_text(errors="replace")
            except OSError:
                continue
            for line in text.splitlines():
                m = DOCKER_RE.search(line)
                if m:
                    image = m.group("image")
                    if "~{" in image and "@sha256" not in image and ":" not in image.split("/")[-1]:
                        continue  # docker: "~{docker_image}" or docker: docker_image: the version is elsewhere
                    name = image.split("/")[-1]
                    tool = re.split(r"[@:]", name)[0]
                    d = DIGEST_RE.search(image)
                    digest = d.group(1) if d else ""
                    version, build = _split_comment(m.group("comment") or "")
                    if not version and not digest:
                        version = name.split(":", 1)[1] if ":" in name else ""
                    key = (tool, digest or version)
                    rec = found.setdefault(key, {"tool": tool, "version": version, "build": build, "digest": digest,
                                                 "image": image.split("}/")[-1] if "~{" in image else image, "files": []})
                    if rel not in rec["files"]:  # type: ignore[operator]
                        rec["files"].append(rel)  # type: ignore[union-attr]
                    continue
                m = VERSION_STRING_RE.search(line)
                if m and m.group("tool") not in ("workflow", "ugc_wgw"):
                    tool, version = m.group("tool"), m.group("version")
                    key = (tool, version)
                    rec = found.setdefault(key, {"tool": tool, "version": version, "build": "", "digest": "",
                                                 "image": "", "files": []})
                    if rel not in rec["files"]:  # type: ignore[operator]
                        rec["files"].append(rel)  # type: ignore[union-attr]
    return sorted(found.values(), key=lambda r: (str(r["tool"]).lower(), str(r["version"])))
