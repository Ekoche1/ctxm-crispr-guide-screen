"""Per-position conservation of a panel relative to a reference sequence.

Each panel sequence is globally aligned to the reference (substitutions and
gaps both allowed). For every reference position we record the fraction of
panel sequences that carry the same base there. A reference position that falls
in a deletion in a panel sequence counts as *not* conserved for that sequence.
"""

from __future__ import annotations

import numpy as np
from Bio import Align


def _aligner() -> Align.PairwiseAligner:
    al = Align.PairwiseAligner()
    al.mode = "global"
    al.match_score, al.mismatch_score = 1, -1
    al.open_gap_score, al.extend_gap_score = -3, -1
    return al


def match_vector(ref: str, other: str, aligner: Align.PairwiseAligner | None = None) -> np.ndarray:
    """Boolean array over reference positions: True where ``other`` has the same base."""
    aligner = aligner or _aligner()
    aln = aligner.align(ref, other)[0]
    vec = np.zeros(len(ref), dtype=bool)
    for (r0, r1), (o0, o1) in zip(*aln.aligned):
        vec[r0:r1] = np.frombuffer(ref[r0:r1].encode(), dtype=np.uint8) == np.frombuffer(
            other[o0:o1].encode(), dtype=np.uint8
        )
    return vec


def conservation_profile(ref: str, panel_seqs: list[str]) -> np.ndarray:
    """Fraction of panel sequences matching the reference at each reference position."""
    if not panel_seqs:
        raise ValueError("empty panel")
    al = _aligner()
    total = np.zeros(len(ref))
    for s in panel_seqs:
        total += match_vector(ref, s, al)
    return total / len(panel_seqs)
