import numpy as np
import pytest

from ctxm_screen.conservation import conservation_profile, match_vector

REF = "ATGGTTAAAAAATCACTGCGCCAGTTCACGCTGATGGCGACGGCAACCGTCACGCTG"


def test_identical_sequence_is_fully_conserved():
    assert match_vector(REF, REF).all()


def test_single_substitution_marks_only_that_position():
    other = REF[:10] + ("C" if REF[10] != "C" else "G") + REF[11:]
    v = match_vector(REF, other)
    assert not v[10] and v.sum() == len(REF) - 1


def test_deletion_counts_as_not_conserved():
    other = REF[:20] + REF[23:]  # delete 3 bases
    v = match_vector(REF, other)
    assert v.sum() <= len(REF) - 3 and v[:20].all()


def test_profile_is_fraction_across_panel():
    mut = REF[:5] + ("C" if REF[5] != "C" else "G") + REF[6:]
    prof = conservation_profile(REF, [REF, REF, REF, mut])
    assert prof[5] == pytest.approx(0.75)
    assert prof[0] == pytest.approx(1.0)
    assert np.all((prof >= 0) & (prof <= 1))


def test_empty_panel_raises():
    with pytest.raises(ValueError):
        conservation_profile(REF, [])
