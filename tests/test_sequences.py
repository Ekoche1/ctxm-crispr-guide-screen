import pytest

from ctxm_screen.sequences import (
    Record,
    gc_fraction,
    parse_fasta,
    reverse_complement,
    validate_records,
)


def test_reverse_complement_known_value():
    assert reverse_complement("AACG") == "CGTT"


def test_reverse_complement_is_involution():
    s = "ATGCCGTTAGC"
    assert reverse_complement(reverse_complement(s)) == s


def test_parse_fasta_multiline_and_description(tmp_path):
    f = tmp_path / "x.fa"
    f.write_text(">a first record\nacgt\nTTGG\n\n>b\nCC\n")
    recs = list(parse_fasta(f))
    assert recs == [Record("a", "first record", "ACGTTTGG"), Record("b", "", "CC")]


def test_parse_fasta_rejects_data_before_header(tmp_path):
    f = tmp_path / "bad.fa"
    f.write_text("ACGT\n>a\nAC\n")
    with pytest.raises(ValueError):
        list(parse_fasta(f))


def test_validate_flags_ambiguous_duplicate_short():
    recs = [
        Record("ok", "", "ACGTACGT"),
        Record("amb", "", "ACGNACGT"),
        Record("ok", "", "ACGTACGT"),
        Record("tiny", "", "AC"),
    ]
    kept, excl = validate_records(recs, min_length=5)
    assert [r.id for r in kept] == ["ok"]
    reasons = {e.id: e.reason for e in excl}
    assert "non-ACGT" in reasons["amb"] and "N" in reasons["amb"]
    assert reasons["tiny"].startswith("shorter")
    assert any(e.reason == "duplicate identifier" for e in excl)


def test_gc_fraction():
    assert gc_fraction("GGCC") == 1.0
    assert gc_fraction("ATAT") == 0.0
    assert gc_fraction("") == 0.0
