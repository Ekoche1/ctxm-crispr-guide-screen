"""The README's numbers must match the result files (run scripts 01-07 first)."""

import json
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
R = ROOT / "results"
pytestmark = pytest.mark.skipif(not (R / "run_info.json").exists(), reason="analysis outputs not generated")


def readme() -> str:
    return (ROOT / "README.md").read_text(encoding="utf-8")


def test_readme_table_is_exactly_the_generated_table():
    table = (R / "passing_guides.md").read_text(encoding="utf-8").strip()
    assert table in readme()


def test_readme_panel_and_candidate_counts_match_run_info():
    info = json.loads((R / "run_info.json").read_text())
    text = readme()
    for n in info["panel_sizes"].values():
        assert f"**{n}**" in text or f"{n} " in text
    assert f"{info['candidates']} PAM-valid" in text
    assert f"**{info['candidates_pass_both']} passed both**" in text


def test_passing_guides_agree_between_tables():
    sm = pd.read_csv(R / "guide_summary.tsv", sep="\t")
    gs = pd.read_csv(R / "genome_pilot_guide_summary.tsv", sep="\t")
    assert set(sm[sm.passes_both].guide) <= set(gs.guide)
    assert (gs.positives_with_exact_pam_hit <= gs.positives_total).all()


def test_independent_check_passed():
    chk = json.loads((R / "independent_check.json").read_text())
    assert chk["disagreements"] == 0 and not chk["coordinate_or_pam_errors"]
