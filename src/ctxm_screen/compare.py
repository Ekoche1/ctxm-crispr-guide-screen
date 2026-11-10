"""Find the best match of a guide in each panel sequence, allowing mismatches.

For every sequence we look at every 20-nt window on both strands and count
substitutions against the guide (no insertions/deletions are considered).
Two numbers are reported per guide x sequence:

``min_mm_any``  fewest mismatches over all windows, ignoring the PAM.
``min_mm_pam``  fewest mismatches over windows that have a TTTV PAM directly
                5' of them. This is the case that matters for Cas12a cutting.

``NO_HIT`` (99) is returned when no window exists (sequence too short, or no
PAM-valid window for ``min_mm_pam``).

Limitations (also in docs/METHODS.md): a plain mismatch count treats all
positions equally. Real Cas12a is much less tolerant of mismatches in the
PAM-proximal "seed" than elsewhere, and bulges are ignored.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.lib.stride_tricks import sliding_window_view

from .guides import PAM_LEN
from .sequences import reverse_complement

NO_HIT = 99


def _encode(seq: str) -> np.ndarray:
    return np.frombuffer(seq.encode("ascii"), dtype=np.uint8)


@dataclass
class IndexedSequence:
    """A panel sequence prepared for fast window comparison on both strands."""

    id: str
    seq: str
    spacer_len: int = 20

    def __post_init__(self) -> None:
        self._strands = []
        for s in (self.seq, reverse_complement(self.seq)):
            n = len(s)
            if n < self.spacer_len:
                continue
            enc = _encode(s)
            windows = sliding_window_view(enc, self.spacer_len)
            self._strands.append((windows, _pam_mask(enc, windows.shape[0])))

    def best_matches(self, spacer: str) -> tuple[int, int]:
        """Return (min_mm_any, min_mm_pam) for one guide spacer."""
        g = _encode(spacer)
        best_any, best_pam = NO_HIT, NO_HIT
        for windows, pam_ok in self._strands:
            mm = (windows != g).sum(axis=1)
            best_any = min(best_any, int(mm.min()))
            if pam_ok.any():
                best_pam = min(best_pam, int(mm[pam_ok].min()))
        return best_any, best_pam

    def count_exact_with_pam(self, spacer: str) -> int:
        """Number of exact, PAM-valid matches on both strands (copy number)."""
        g = _encode(spacer)
        total = 0
        for windows, pam_ok in self._strands:
            total += int(((windows != g).sum(axis=1) == 0)[pam_ok].sum())
        return total


_A, _C, _G, _T = (ord(c) for c in "ACGT")


def _pam_mask(enc: np.ndarray, n_windows: int) -> np.ndarray:
    """Vectorised TTTV test: True for window i if enc[i-4:i] is T,T,T,[A/C/G].

    Equivalent to ``has_pam`` on every window, but fast enough for whole
    genomes. Any non-ACGT symbol (e.g. N) never counts as V.
    """
    mask = np.zeros(n_windows, dtype=bool)
    if n_windows > PAM_LEN:
        i = np.arange(PAM_LEN, n_windows)
        v = enc[i - 1]
        mask[PAM_LEN:] = (
            (enc[i - 4] == _T) & (enc[i - 3] == _T) & (enc[i - 2] == _T)
            & ((v == _A) | (v == _C) | (v == _G))
        )
    return mask


def index_panel(records, spacer_len: int = 20) -> list[IndexedSequence]:
    return [IndexedSequence(r.id, r.seq, spacer_len) for r in records]


def screen_guide(spacer: str, panel: list[IndexedSequence]) -> list[tuple[str, int, int]]:
    """(sequence id, min_mm_any, min_mm_pam) for every sequence in the panel."""
    return [(s.id, *s.best_matches(spacer)) for s in panel]
