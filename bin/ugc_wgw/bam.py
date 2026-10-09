"""Cheap inspection of input HiFi BAMs with the standard library (docs/DESIGN.md §8.1).

A BAM is a series of BGZF blocks (gzip members of at most 64 KB with a `BC` extra field), which zlib
inflates, so the driver can look at a file without samtools: the end-of-file marker (a truncated copy
lacks it), the header (read groups, reference sequences), the first records (none at all is the
failure that otherwise surfaces hours later in pbmm2 or the merge) and, from a sample of records and
the compressed bytes they occupied, an estimate of the read count and the bases. A `.pbi` next to
the BAM gives the exact read count for free. Nothing here decides; `samples.py` turns the findings
into errors and warnings, and `preflight.py` repeats the cheap part before a run starts.
"""
from __future__ import annotations

import gzip
import struct
import zlib
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import BinaryIO

# BGZF end-of-file block (SAM spec §4.1.2): an empty deflate block with the BC extra field.
EOF_BLOCK = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")
BAM_MAGIC = b"BAM\x01"
PBI_MAGIC = b"PBI\x01"
SAMPLE_RECORDS = 1000          # records read for the yield estimate
HUMAN_GENOME_GB = 3.1          # haploid bases of GRCh38, for the coverage figure


@dataclass
class BamInfo:
    """What one look at a BAM found. `problems` make the file unusable; `warnings` are advisory."""
    path: str
    bytes: int = 0
    mtime: float = 0.0
    eof_ok: bool = False
    aligned: bool = False
    n_ref: int = 0
    read_groups: int = 0
    movies: list[str] = field(default_factory=list)
    samples: list[str] = field(default_factory=list)
    records_sampled: int = 0
    bases_sampled: int = 0
    reads: int = 0                 # exact (pbi, or the whole file was read) or estimated
    bases: int = 0                 # exact only when the whole file was read; else estimated
    reads_exact: bool = False
    pbi: bool = False
    problems: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not self.problems

    @property
    def gb(self) -> float:
        return self.bases / 1e9

    @property
    def coverage(self) -> float:
        return self.gb / HUMAN_GENOME_GB

    @property
    def mean_read_length(self) -> int:
        return self.bases_sampled // self.records_sampled if self.records_sampled else 0

    def to_dict(self) -> dict[str, object]:
        d = asdict(self)
        d["gb"] = round(self.gb, 3)
        d["mean_read_length"] = self.mean_read_length
        return d

    @classmethod
    def from_dict(cls, d: dict[str, object]) -> "BamInfo":
        return cls(**{k: v for k, v in d.items() if k in cls.__dataclass_fields__})  # type: ignore[arg-type]


def has_eof_marker(path: Path) -> bool:
    """True when the file ends with the BGZF end-of-file block (a complete BGZF file always does)."""
    try:
        size = path.stat().st_size
        if size < len(EOF_BLOCK):
            return False
        with open(path, "rb") as fh:
            fh.seek(size - len(EOF_BLOCK))
            return fh.read(len(EOF_BLOCK)) == EOF_BLOCK
    except OSError:
        return False


def pbi_reads(path: Path) -> int | None:
    """Read count from a PacBio index next to the BAM, or None when there is none or it is unreadable."""
    pbi = Path(str(path) + ".pbi")
    if not pbi.is_file():
        return None
    try:
        with gzip.open(pbi, "rb") as fh:
            head = fh.read(14)
    except (OSError, EOFError, zlib.error):
        return None
    if len(head) < 14 or head[:4] != PBI_MAGIC:
        return None
    return struct.unpack("<I", head[10:14])[0]   # magic(4) version(4) flags(2) n_reads(4)


class _Stream:
    """Sequential reader over the inflated payload of a BGZF file that knows how many compressed bytes the
    bytes handed out so far occupied (pro rata within the current block)."""

    def __init__(self, raw: BinaryIO):
        self.raw = raw
        self.buf = b""
        self.pos = 0
        self.block_bytes = 0        # compressed size of the block in buf
        self.before = 0             # compressed bytes of the blocks before buf
        self.ended = False

    def _next_block(self) -> bool:
        head = self.raw.read(12)
        if not head:
            self.ended = True
            return False
        if len(head) < 12 or head[:4] != b"\x1f\x8b\x08\x04":
            raise ValueError("not a BGZF file" if self.before == 0 else "malformed BGZF block")
        xlen, = struct.unpack("<H", head[10:12])
        extra = self.raw.read(xlen)
        i = extra.find(b"BC\x02\x00")
        if i < 0 or len(extra) < i + 6:
            raise ValueError("not a BGZF file (gzip member without the BC field)")
        bsize, = struct.unpack("<H", extra[i + 4:i + 6])
        rest = self.raw.read(bsize + 1 - 12 - xlen)
        if len(rest) < 8:
            raise ValueError("truncated BGZF block")
        payload = zlib.decompress(rest[:-8], -15)
        self.before += self.block_bytes
        self.buf, self.pos, self.block_bytes = payload, 0, 12 + xlen + len(rest)
        return True

    def read(self, n: int) -> bytes:
        out = bytearray()
        while len(out) < n:
            if self.pos >= len(self.buf):
                if not self._next_block():
                    break
                continue
            take = self.buf[self.pos:self.pos + n - len(out)]
            out += take
            self.pos += len(take)
        return bytes(out)

    def at_end(self) -> bool:
        if self.pos < len(self.buf):
            return False
        while self._next_block():
            if self.buf:
                return False
        return True

    @property
    def consumed(self) -> float:
        frac = self.pos / len(self.buf) if self.buf else 1.0
        return self.before + frac * self.block_bytes


def _header(fh: _Stream) -> tuple[str, int]:
    """(SAM header text, number of reference sequences); raises ValueError on a non-BAM stream."""
    if fh.read(4) != BAM_MAGIC:
        raise ValueError("not a BAM file (bad magic)")
    l_text, = struct.unpack("<i", fh.read(4))
    text = fh.read(l_text).decode("utf-8", errors="replace")
    n_ref, = struct.unpack("<i", fh.read(4))
    for _ in range(n_ref):
        l_name, = struct.unpack("<i", fh.read(4))
        fh.read(l_name + 4)
    return text.rstrip("\x00"), n_ref


def _records(fh: _Stream, limit: int) -> tuple[int, int]:
    """Read up to `limit` alignment records: (records, bases)."""
    n = bases = 0
    while n < limit:
        head = fh.read(4)
        if len(head) < 4:
            if head:
                raise ValueError("record cut short")
            break
        block_size, = struct.unpack("<i", head)
        if block_size < 32:
            raise ValueError(f"corrupt record (block size {block_size})")
        rec = fh.read(block_size)
        if len(rec) < block_size:
            raise ValueError("record cut short")
        l_seq, = struct.unpack("<i", rec[16:20])
        n += 1
        bases += max(l_seq, 0)
    return n, bases


def inspect(path: Path, *, sample_records: int = SAMPLE_RECORDS) -> BamInfo:
    """Look at one BAM: existence, size, EOF marker, header, first records, pbi; never raises."""
    info = BamInfo(path=str(path))
    try:
        st = path.stat()
    except OSError as exc:
        info.problems.append(f"not readable: {exc.strerror or exc}")
        return info
    info.bytes, info.mtime = st.st_size, st.st_mtime
    if st.st_size == 0:
        info.problems.append("empty file (0 bytes)")
        return info
    info.eof_ok = has_eof_marker(path)
    if not info.eof_ok:
        info.problems.append("truncated: the BGZF end-of-file marker is missing (incomplete copy?)")
    rgs: list[list[str]] = []
    n = bases = 0
    at_end = False
    consumed = 0.0
    try:
        with open(path, "rb") as raw:
            fh = _Stream(raw)
            text, info.n_ref = _header(fh)
            info.aligned = info.n_ref > 0
            rgs = [line.split("\t") for line in text.splitlines() if line.startswith("@RG")]
            info.read_groups = len(rgs)
            for fields in rgs:
                for f in fields[1:]:
                    if f.startswith("PU:") and f[3:] not in info.movies:
                        info.movies.append(f[3:])
                    elif f.startswith("SM:") and f[3:] not in info.samples:
                        info.samples.append(f[3:])
            n, bases = _records(fh, sample_records)
            consumed = fh.consumed
            at_end = n < sample_records or fh.at_end()
    except (OSError, EOFError, zlib.error, ValueError, struct.error) as exc:
        if info.eof_ok:   # a truncated file is already reported; anything else is damage of its own
            info.problems.append(f"unreadable BAM: {exc}")
        elif n == 0 and not rgs:
            return info
    info.records_sampled, info.bases_sampled = n, bases
    if info.problems:
        return info
    if n == 0:
        info.problems.append("no reads (header only)")
        return info
    exact = pbi_reads(path)
    if exact is not None:
        info.pbi, info.reads, info.reads_exact = True, exact, True
    elif at_end:
        info.reads, info.reads_exact = n, True
    else:
        info.reads = int(st.st_size / (consumed / n)) if consumed else n
    info.bases = int(info.reads * (bases / n))
    if info.aligned:
        info.warnings.append(f"aligned input ({info.n_ref} reference sequences): upstream strips the alignments "
                             "and disables alignment chunking")
    if not rgs:
        info.warnings.append("no @RG read group in the header")
    return info


def fmt_bases(bases: float) -> str:
    if bases >= 1e9:
        return f"{bases / 1e9:.1f} Gb"
    if bases >= 1e6:
        return f"{bases / 1e6:.1f} Mb"
    return f"{bases / 1e3:.0f} kb"


def fmt_bytes(n: float) -> str:
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024:
            return f"{n:.0f} B" if unit == "B" else f"{n:.1f} {unit}"
        n /= 1024.0
    return f"{n:.1f} TB"


def quick_check(path: Path, expected_bytes: int | None = None) -> str | None:
    """The cheap part for the moment before a run starts: present, not empty, same size as registered,
    EOF marker in place. Returns the problem or None."""
    try:
        st = path.stat()
    except OSError as exc:
        return f"not readable: {exc.strerror or exc}"
    if st.st_size == 0:
        return "empty file (0 bytes)"
    if expected_bytes is not None and st.st_size != expected_bytes:
        return f"size changed since registration ({fmt_bytes(expected_bytes)} -> {fmt_bytes(st.st_size)})"
    if not has_eof_marker(path):
        return "truncated: the BGZF end-of-file marker is missing"
    return None
