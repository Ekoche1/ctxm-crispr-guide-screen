import random

from ctxm_screen.genome import locate_gene
from ctxm_screen.sequences import reverse_complement

rng = random.Random(11)
REF = "".join(rng.choice("ACGT") for _ in range(300))


def junk(n: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(n))


def test_finds_forward_copy_with_exact_coordinates():
    left = junk(2000)
    genome = left + REF + junk(1500)
    hit = locate_gene(genome, REF)
    assert hit.found and hit.strand == "+"
    assert (hit.start, hit.end) == (len(left), len(left) + len(REF))
    assert hit.identity_over_ref == 1.0 and hit.coverage == 1.0


def test_finds_reverse_strand_copy_and_maps_coordinates_back():
    left = junk(1800)
    genome = left + reverse_complement(REF) + junk(900)
    hit = locate_gene(genome, REF)
    assert hit.found and hit.strand == "-"
    assert genome[hit.start : hit.end] == reverse_complement(REF)


def test_reports_substitutions_in_identity():
    mutated = REF[:100] + ("A" if REF[100] != "A" else "C") + REF[101:]
    mutated = mutated[:200] + ("A" if mutated[200] != "A" else "C") + mutated[201:]
    hit = locate_gene(junk(1000) + mutated + junk(1000), REF)
    assert hit.found
    assert abs(hit.identity_over_ref - (len(REF) - 2) / len(REF)) < 1e-9


def test_absent_gene_is_not_found():
    assert not locate_gene(junk(5000), REF).found
