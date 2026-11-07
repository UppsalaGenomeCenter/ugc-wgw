"""Cheap inspection of a given VCF with the standard library (the vcf entry kind, docs/DESIGN.md §8.1).

A BGZF-compressed VCF is read block by block with the BAM module's stream reader until the `#CHROM` line: the
header gives the sample columns (exactly one is required), the contigs and the producing programs. Nothing
here decides; `samples.py` turns the findings into errors and warnings, and the `intake` stage repeats the
checks against the installed reference.
"""
from __future__ import annotations

import zlib
from dataclasses import asdict, dataclass, field
from pathlib import Path

from . import bam

HEADER_LIMIT = 64 << 20   # bytes of inflated header before giving up (a VCF header is kilobytes to a few MB)
SOURCE_PREFIXES = ("##fileformat", "##source", "##DeepVariant", "##sawfish", "##GLnexus", "##bcftools_", "##reference")


@dataclass
class VcfInfo:
    """What one look at a VCF found. `problems` make the file unusable; `warnings` are advisory."""
    path: str
    bytes: int = 0
    mtime: float = 0.0
    eof_ok: bool = False
    bgzf: bool = False
    fileformat: str = ""
    samples: list[str] = field(default_factory=list)
    n_contigs: int = 0
    contigs: list[str] = field(default_factory=list)   # the first few contig names
    sources: list[str] = field(default_factory=list)   # ##source and tool version lines, as written
    problems: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not self.problems

    def to_dict(self) -> dict[str, object]:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict[str, object]) -> "VcfInfo":
        return cls(**{k: v for k, v in d.items() if k in cls.__dataclass_fields__})  # type: ignore[arg-type]


def inspect(path: Path) -> VcfInfo:
    """Look at one VCF: existence, size, BGZF, EOF marker, header (samples, contigs, sources); never raises."""
    info = VcfInfo(path=str(path))
    try:
        st = path.stat()
    except OSError as exc:
        info.problems.append(f"not readable: {exc.strerror or exc}")
        return info
    info.bytes, info.mtime = st.st_size, st.st_mtime
    if st.st_size == 0:
        info.problems.append("empty file (0 bytes)")
        return info
    info.eof_ok = bam.has_eof_marker(path)
    text = b""
    try:
        with open(path, "rb") as raw:
            stream = bam._Stream(raw)
            info.bgzf = True
            while len(text) < HEADER_LIMIT:
                chunk = stream.read(1 << 16)
                if not chunk:
                    break
                text += chunk
                nl = text.find(b"\n#CHROM")
                if nl >= 0 or (text.startswith(b"#CHROM")):
                    end = text.find(b"\n", nl + 1) if nl >= 0 else text.find(b"\n")
                    if end >= 0:
                        text = text[:end]
                        break
    except (OSError, EOFError, zlib.error, ValueError) as exc:
        if not info.bgzf or "BGZF" in str(exc):
            info.problems.append("not BGZF-compressed (a plain-text or gzip VCF): bgzip and index it")
        else:
            info.problems.append(f"unreadable VCF: {exc}")
        return info
    if not info.eof_ok:
        info.problems.append("truncated: the BGZF end-of-file marker is missing (incomplete copy?)")
    lines = text.decode("utf-8", errors="replace").splitlines()
    if not lines or not lines[0].startswith("##fileformat=VCF"):
        info.problems.append("not a VCF (no ##fileformat=VCF header line)")
        return info
    info.fileformat = lines[0][len("##fileformat="):]
    header_line = None
    for line in lines:
        if line.startswith("##contig=<"):
            info.n_contigs += 1
            if len(info.contigs) < 10:
                body = line[len("##contig=<"):].rstrip(">")
                for part in body.split(","):
                    if part.startswith("ID="):
                        info.contigs.append(part[3:])
                        break
        elif line.startswith(SOURCE_PREFIXES) and not line.startswith("##fileformat"):
            info.sources.append(line[2:][:300])
        elif line.startswith("#CHROM"):
            header_line = line
    if header_line is None:
        info.problems.append("no #CHROM header line (header only, or the header is longer than 64 MB)")
        return info
    cols = header_line.split("\t")
    info.samples = cols[9:] if len(cols) > 9 else []
    if len(info.samples) != 1:
        info.problems.append(f"{len(info.samples)} sample column(s); exactly one is required" if info.samples else "a sites-only VCF (no sample column)")
    if info.n_contigs == 0:
        info.warnings.append("no ##contig lines in the header: the intake stage checks the records' contigs instead")
    return info
