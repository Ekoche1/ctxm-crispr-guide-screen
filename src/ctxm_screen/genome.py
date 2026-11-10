"""Locate a reference gene inside a genome assembly.

We do not trust record annotations alone ("CTX-M-15" appearing in a GenBank
title or feature). Instead we find the gene by sequence:

1. every 16-mer of the reference is indexed (``SEED_K``);
2. every genome position on both strands is looked up; each seed hit votes for
   a "diagonal" (genome position minus reference position);
3. the best-supported diagonal gives a candidate region, which is verified by
   a local alignment. Identity is reported over the *whole* reference length,
   so a partial or divergent hit gets a low score rather than a flattering one.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass

from Bio import Align

from .sequences import reverse_complement

SEED_K = 16
MIN_SEED_VOTES = 5


@dataclass(frozen=True)
class GeneHit:
    found: bool
    strand: str  # "+" / "-" / "."
    start: int  # 0-based start on the *forward* genome strand; -1 if not found
    end: int
    seed_votes: int
    identity_over_ref: float  # matching bases / reference length
    coverage: float  # aligned reference bases / reference length


NOT_FOUND = GeneHit(False, ".", -1, -1, 0, 0.0, 0.0)


def _seed_index(ref: str, k: int) -> dict[str, list[int]]:
    idx: dict[str, list[int]] = {}
    for i in range(len(ref) - k + 1):
        idx.setdefault(ref[i : i + k], []).append(i)
    return idx


def locate_gene(genome: str, ref: str, k: int = SEED_K) -> GeneHit:
    index = _seed_index(ref, k)
    n = len(genome)
    best: tuple[int, str, int] | None = None  # (votes, strand, diagonal)
    for strand, s in (("+", genome), ("-", reverse_complement(genome))):
        votes: Counter[int] = Counter()
        for i in range(len(s) - k + 1):
            hit = index.get(s[i : i + k])
            if hit:
                for r in hit:
                    votes[i - r] += 1
        if votes:
            diag, v = votes.most_common(1)[0]
            if best is None or v > best[0]:
                best = (v, strand, diag)
    if best is None or best[0] < MIN_SEED_VOTES:
        return GeneHit(False, ".", -1, -1, best[0] if best else 0, 0.0, 0.0)

    v, strand, diag = best
    s = genome if strand == "+" else reverse_complement(genome)
    pad = 60
    lo, hi = max(0, diag - pad), min(len(s), diag + len(ref) + pad)
    region = s[lo:hi]
    al = Align.PairwiseAligner()
    al.mode = "local"
    al.match_score, al.mismatch_score = 1, -1
    al.open_gap_score, al.extend_gap_score = -3, -1
    aln = al.align(ref, region)[0]
    matches = aligned = 0
    for (r0, r1), (g0, g1) in zip(*aln.aligned):
        aligned += r1 - r0
        matches += sum(a == b for a, b in zip(ref[r0:r1], region[g0:g1]))
    g_start, g_end = lo + int(aln.aligned[1][0][0]), lo + int(aln.aligned[1][-1][1])
    if strand == "-":  # map back to forward-strand coordinates
        g_start, g_end = n - g_end, n - g_start
    return GeneHit(True, strand, g_start, g_end, v, matches / len(ref), aligned / len(ref))
