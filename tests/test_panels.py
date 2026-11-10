from ctxm_screen.panels import (
    MIN_IDENTITY,
    allele_name,
    accession,
    assign_groups,
    sample_evenly,
)
from ctxm_screen.sequences import Record

# A 120 bp random-looking sequence and mutated copies, to test the rule.
BASE = (
    "ATGGTTAAAAAATCACTGCGCCAGTTCACGCTGATGGCGACGGCAACCGTCACGCTGTTGTTAGGAAGT"
    "GTGCCGCTGTATGCGCAAACGGCGGACGTACAGCAAAAACTTGCCGAACTTGAACGGCAGTCGGGAGGC"
)


def mutate(seq: str, n: int) -> str:
    """Flip the first n bases to a different base (A<->C, G<->T)."""
    swap = {"A": "C", "C": "A", "G": "T", "T": "G"}
    return "".join(swap[c] if i < n else c for i, c in enumerate(seq))


def rec(allele: str, seq: str) -> Record:
    return Record(f"WP_1|NG_{allele}.1|1|1|{allele}|fam|prod", "", seq)


def test_header_parsing():
    r = rec("blaX-1", "ACGT")
    assert allele_name(r) == "blaX-1" and accession(r) == "NG_blaX-1.1"


def test_close_allele_assigned_to_nearest_anchor():
    anchors = {"G1": rec("a1", BASE), "G2": rec("a2", mutate(BASE, 60))}
    close = rec("q", mutate(BASE, 3))
    (a,) = assign_groups([close], anchors)
    assert a.group == "G1" and a.best_identity > MIN_IDENTITY


def test_equidistant_allele_is_unassigned_despite_high_identity():
    # G2 differs from G1 at the first 6 bases; the query differs from G1 at the
    # first 3 only, so it is 3 bases from each anchor: identity is high (~0.97)
    # for both, but the lead is zero, so it must not be forced into a group.
    anchors = {"G1": rec("a1", BASE), "G2": rec("a2", mutate(BASE, 6))}
    query = rec("q", mutate(BASE, 3))
    (a,) = assign_groups([query], anchors)
    assert a.best_identity > MIN_IDENTITY
    assert a.group == "unassigned" and "lead" in a.note


def test_distant_allele_is_unassigned():
    anchors = {"G1": rec("a1", BASE), "G2": rec("a2", mutate(BASE, 30))}
    far = rec("f", mutate(BASE, 100))
    (a,) = assign_groups([far], anchors)
    assert a.group == "unassigned"


def test_sample_evenly_is_deterministic_and_capped():
    names = [f"x{i:03d}" for i in range(100)]
    s1, s2 = sample_evenly(names, 10), sample_evenly(list(reversed(names)), 10)
    assert s1 == s2 and len(s1) == 10 and s1[0] == "x000"
    assert sample_evenly(["b", "a"], 10) == ["a", "b"]
