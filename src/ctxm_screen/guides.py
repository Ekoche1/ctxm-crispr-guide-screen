"""Enumerate candidate Cas12a protospacers on a reference sequence.

Cas12a (LbCas12a / AsCas12a) needs a 5'-TTTV-3' PAM (V = A, C or G)
immediately 5' of the 20-nt protospacer on the strand that carries it. The
crRNA spacer has the same sequence as that protospacer (with U for T).

We look on both strands. Positions are reported in the *forward* coordinates of
the reference (0-based, start of the protospacer, ``end`` exclusive) so they
can be plotted along the gene.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

from .sequences import gc_fraction, reverse_complement

PAM_REGEX = re.compile(r"TTT[ACG]")
SPACER_LEN = 20
PAM_LEN = 4


@dataclass(frozen=True)
class Candidate:
    id: str
    strand: str  # "+" or "-"
    pam: str
    spacer: str  # 5'->3' on the strand carrying the PAM
    start: int  # 0-based forward coordinate, inclusive
    end: int  # exclusive
    gc: float


def has_pam(upstream4: str) -> bool:
    """True if the 4 nt immediately 5' of a protospacer form TTTV."""
    return bool(PAM_REGEX.fullmatch(upstream4))


def enumerate_candidates(
    ref: str, ref_name: str = "ref", spacer_len: int = SPACER_LEN
) -> list[Candidate]:
    """All PAM-valid protospacers of ``spacer_len`` on both strands of ``ref``."""
    out: list[Candidate] = []
    n = len(ref)
    for strand, seq in (("+", ref), ("-", reverse_complement(ref))):
        for i in range(PAM_LEN, n - spacer_len + 1):
            pam = seq[i - PAM_LEN : i]
            if not has_pam(pam):
                continue
            spacer = seq[i : i + spacer_len]
            if strand == "+":
                start, end = i, i + spacer_len
            else:  # map reverse-strand coordinates back onto the forward strand
                start, end = n - (i + spacer_len), n - i
            out.append(
                Candidate(
                    id=f"{ref_name}_{strand}{start + 1}",
                    strand=strand,
                    pam=pam,
                    spacer=spacer,
                    start=start,
                    end=end,
                    gc=gc_fraction(spacer),
                )
            )
    return sorted(out, key=lambda c: (c.start, c.strand))
