import random

import regex

from ctxm_screen.compare import NO_HIT, IndexedSequence
from ctxm_screen.guides import enumerate_candidates, has_pam
from ctxm_screen.sequences import reverse_complement

GUIDE = "ACGTACGTAGCTAGCTAGCA"  # 20 nt


def test_has_pam_accepts_tttv_rejects_tttt():
    assert all(has_pam(p) for p in ("TTTA", "TTTC", "TTTG"))
    assert not has_pam("TTTT") and not has_pam("ATTA") and not has_pam("TTT")


def test_enumerate_finds_forward_candidate_with_coordinates():
    ref = "GG" + "TTTA" + GUIDE + "CC"
    cands = [c for c in enumerate_candidates(ref) if c.strand == "+"]
    assert len(cands) == 1
    c = cands[0]
    assert c.spacer == GUIDE and c.pam == "TTTA"
    assert (c.start, c.end) == (6, 26) and ref[c.start : c.end] == GUIDE


def test_enumerate_finds_reverse_strand_candidate_and_maps_back():
    # A guide on the minus strand: the plus strand carries the reverse complement.
    plus = reverse_complement("TTTC" + GUIDE) + "AAAA"
    cands = [c for c in enumerate_candidates(plus) if c.strand == "-"]
    hit = [c for c in cands if c.spacer == GUIDE]
    assert len(hit) == 1
    c = hit[0]
    assert reverse_complement(plus[c.start : c.end]) == GUIDE


def test_exact_match_with_pam_is_zero_mismatches():
    s = IndexedSequence("t", "AA" + "TTTG" + GUIDE + "GG")
    assert s.best_matches(GUIDE) == (0, 0)


def test_missing_pam_gives_zero_any_but_no_pam_hit():
    s = IndexedSequence("t", "AA" + "TTTT" + GUIDE + "GG")
    any_mm, pam_mm = s.best_matches(GUIDE)
    assert any_mm == 0 and pam_mm > 0  # the TTTT is not a valid PAM


def test_mismatch_counted_and_reverse_strand_found():
    mutated = GUIDE[:5] + ("A" if GUIDE[5] != "A" else "C") + GUIDE[6:]
    plus = reverse_complement("TTTA" + mutated)
    s = IndexedSequence("t", "CCCC" + plus + "CCCC")
    assert s.best_matches(GUIDE) == (1, 1)


def test_short_sequence_returns_no_hit():
    assert IndexedSequence("t", "ACGT").best_matches(GUIDE) == (NO_HIT, NO_HIT)


def test_vectorised_pam_mask_matches_per_window_has_pam():
    import numpy as np
    from ctxm_screen.compare import _encode, _pam_mask

    rng = random.Random(3)
    for _ in range(30):
        # include N so we also check that ambiguous bases never count as V
        s = "".join(rng.choice("ACGTTTN") for _ in range(200))
        n_windows = len(s) - 20 + 1
        expected = np.array([i >= 4 and has_pam(s[i - 4 : i]) for i in range(n_windows)])
        assert (_pam_mask(_encode(s), n_windows) == expected).all()


def test_count_exact_with_pam_counts_copies_on_both_strands():
    plus = "AA" + "TTTA" + GUIDE + "GG"
    minus = reverse_complement("CC" + "TTTC" + GUIDE + "AA")
    s = IndexedSequence("t", plus + "CCCCCC" + minus)
    assert s.count_exact_with_pam(GUIDE) == 2


def test_numpy_search_agrees_with_independent_regex_engine():
    """Cross-check the numpy search against the `regex` fuzzy matcher."""
    rng = random.Random(7)
    for _ in range(40):
        seq = "".join(rng.choice("ACGT") for _ in range(300))
        guide = "".join(rng.choice("ACGT") for _ in range(20))
        mine, _ = IndexedSequence("t", seq).best_matches(guide)
        best = NO_HIT
        for strand in (seq, reverse_complement(seq)):
            for k in range(0, 21):
                if regex.search(f"(?:{guide}){{s<={k}}}", strand):
                    best = min(best, k)
                    break
        assert mine == best
