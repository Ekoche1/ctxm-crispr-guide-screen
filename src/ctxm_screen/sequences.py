"""Reading and validating nucleotide sequences.

We keep this deliberately small. A sequence is "clean" only if it contains
A, C, G and T. Anything else (N, IUPAC ambiguity codes, gaps) is reported and
the record is excluded, because a guide-matching step that silently treats an
N as a match would overstate conservation.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

_COMPLEMENT = str.maketrans("ACGT", "TGCA")
_VALID = set("ACGT")


@dataclass(frozen=True)
class Record:
    """One FASTA record: the first word of the header is the identifier."""

    id: str
    description: str
    seq: str


@dataclass(frozen=True)
class Exclusion:
    """A record we refused to use, and why. Nothing is dropped silently."""

    id: str
    reason: str


def parse_fasta(path: str | Path) -> Iterator[Record]:
    """Yield records from a FASTA file. Sequences are upper-cased."""
    header: str | None = None
    chunks: list[str] = []
    with open(path, encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield _make_record(header, chunks)
                header, chunks = line[1:], []
            else:
                if header is None:
                    raise ValueError(f"{path}: sequence data before first '>' header")
                chunks.append(line)
    if header is not None:
        yield _make_record(header, chunks)


def _make_record(header: str, chunks: list[str]) -> Record:
    parts = header.split(None, 1)
    ident = parts[0] if parts else ""
    desc = parts[1] if len(parts) > 1 else ""
    return Record(id=ident, description=desc, seq="".join(chunks).upper().replace("U", "T"))


def reverse_complement(seq: str) -> str:
    """Reverse complement of an A/C/G/T string."""
    return seq.translate(_COMPLEMENT)[::-1]


def validate_records(
    records: list[Record], min_length: int = 1
) -> tuple[list[Record], list[Exclusion]]:
    """Split records into usable ones and explained exclusions.

    Excluded: empty id, duplicate id, too short, or any non-ACGT character.
    """
    kept: list[Record] = []
    excluded: list[Exclusion] = []
    seen: set[str] = set()
    for rec in records:
        if not rec.id:
            excluded.append(Exclusion("(blank)", "empty identifier"))
        elif rec.id in seen:
            excluded.append(Exclusion(rec.id, "duplicate identifier"))
        elif len(rec.seq) < min_length:
            excluded.append(Exclusion(rec.id, f"shorter than {min_length} bp"))
        elif not set(rec.seq) <= _VALID:
            bad = "".join(sorted(set(rec.seq) - _VALID))
            excluded.append(Exclusion(rec.id, f"non-ACGT characters: {bad}"))
        else:
            kept.append(rec)
        seen.add(rec.id)
    return kept, excluded


def gc_fraction(seq: str) -> float:
    """Fraction of G or C. Returns 0.0 for an empty string."""
    return (seq.count("G") + seq.count("C")) / len(seq) if seq else 0.0
