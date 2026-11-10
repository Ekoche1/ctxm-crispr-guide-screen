"""Step 2: screen candidate Cas12a guides from blaCTX-M-15 against all panels.

Reads  data/panels/*.fasta
Writes results/guide_vs_sequence.tsv   (every guide x every sequence)
       results/guide_summary.tsv       (one row per candidate guide)
       results/conservation_profile.tsv (per-position conservation)
       results/run_info.json           (versions, parameters, counts)

Pre-specified screening criteria (decided before looking at results; see
docs/METHODS.md):
  * conserved:    exact match WITH a valid PAM in >= 95% of the target panel
  * specific:     >= 3 mismatches (with valid PAM) to EVERY non-target sequence
These mirror PathoGD's "at least 3 mismatches" rule. They are screening
filters, not measures of diagnostic performance.
"""

import csv
import json
import platform
import sys
from datetime import date
from pathlib import Path

import Bio
import numpy as np
import pandas as pd
import regex

from ctxm_screen import __version__
from ctxm_screen.compare import index_panel, screen_guide
from ctxm_screen.conservation import conservation_profile
from ctxm_screen.guides import enumerate_candidates
from ctxm_screen.sequences import parse_fasta

ROOT = Path(__file__).resolve().parents[1]
PANELS, RESULTS = ROOT / "data" / "panels", ROOT / "results"
RESULTS.mkdir(exist_ok=True)

CONSERVED_MIN = 0.95
SPECIFIC_MIN_MM = 3
REFERENCE_ALLELE = "blaCTX-M-15"
PANEL_FILES = {
    "target_g1": "target_g1.fasta",
    "nontarget_other_ctxm_groups": "nontarget_other_ctxm_groups.fasta",
    "nontarget_other_beta_lactamases": "nontarget_other_beta_lactamases.fasta",
    "unassigned_ctxm": "unassigned_ctxm.fasta",
}


def main() -> None:
    panels = {name: list(parse_fasta(PANELS / f)) for name, f in PANEL_FILES.items()}
    ref = next(r for r in panels["target_g1"] if r.id == REFERENCE_ALLELE)
    print(f"reference {ref.id}: {len(ref.seq)} bp")

    candidates = enumerate_candidates(ref.seq, ref.id)
    print(f"candidate guides (TTTV PAM, 20 nt, both strands): {len(candidates)}")

    indexed = {name: index_panel(recs) for name, recs in panels.items()}
    rows = []
    for c in candidates:
        for name, idx in indexed.items():
            for sid, mm_any, mm_pam in screen_guide(c.spacer, idx):
                rows.append((c.id, name, sid, mm_any, mm_pam))
    long = pd.DataFrame(rows, columns=["guide", "panel", "sequence", "min_mm_any_pam", "min_mm_with_pam"])
    long.to_csv(RESULTS / "guide_vs_sequence.tsv", sep="\t", index=False, lineterminator="\n")

    summary = []
    for c in candidates:
        g = long[long.guide == c.id]
        tgt = g[g.panel == "target_g1"]
        nt = g[g.panel.isin(["nontarget_other_ctxm_groups", "nontarget_other_beta_lactamases"])]
        nt_ctxm = g[g.panel == "nontarget_other_ctxm_groups"]
        nt_bla = g[g.panel == "nontarget_other_beta_lactamases"]
        una = g[g.panel == "unassigned_ctxm"]
        n_exact = int((tgt.min_mm_with_pam == 0).sum())
        cons = n_exact / len(tgt)
        min_nt = int(nt.min_mm_with_pam.min())
        summary.append({
            "guide": c.id, "strand": c.strand, "start_1based": c.start + 1, "end": c.end,
            "pam": c.pam, "spacer_5to3": c.spacer, "gc": round(c.gc, 2),
            "target_n": len(tgt), "target_exact_with_pam": n_exact,
            "target_conservation": round(cons, 4),
            "target_max_mm": int(tgt.min_mm_with_pam.max()),
            "nontarget_n": len(nt),
            "nontarget_min_mm_with_pam": min_nt,
            "nontarget_min_mm_any_pam": int(nt.min_mm_any_pam.min()),
            "nontarget_ctxm_min_mm": int(nt_ctxm.min_mm_with_pam.min()),
            "nontarget_bla_min_mm": int(nt_bla.min_mm_with_pam.min()),
            "n_nontarget_le2mm": int((nt.min_mm_with_pam <= 2).sum()),
            "unassigned_exact_with_pam": int((una.min_mm_with_pam == 0).sum()),
            "unassigned_n": len(una),
            "conserved": cons >= CONSERVED_MIN,
            "specific": min_nt >= SPECIFIC_MIN_MM,
        })
    sm = pd.DataFrame(summary)
    sm["passes_both"] = sm.conserved & sm.specific
    sm.to_csv(RESULTS / "guide_summary.tsv", sep="\t", index=False, lineterminator="\n")

    prof = conservation_profile(ref.seq, [r.seq for r in panels["target_g1"]])
    prof_other = conservation_profile(ref.seq, [r.seq for r in panels["nontarget_other_ctxm_groups"]])
    pd.DataFrame({
        "position_1based": np.arange(1, len(ref.seq) + 1), "ref_base": list(ref.seq),
        "conservation_target_g1": prof.round(4),
        "identity_other_ctxm_groups": prof_other.round(4),
    }).to_csv(RESULTS / "conservation_profile.tsv", sep="\t", index=False, lineterminator="\n")

    info = {
        "run_date": date.today().isoformat(),
        "python": sys.version.split()[0], "platform": platform.platform(),
        "ctxm_screen": __version__, "biopython": Bio.__version__, "numpy": np.__version__,
        "pandas": pd.__version__, "regex": regex.__version__,
        "reference": ref.id, "reference_length_bp": len(ref.seq),
        "spacer_len": 20, "pam": "TTTV (5' of protospacer)",
        "conserved_min": CONSERVED_MIN, "specific_min_mm": SPECIFIC_MIN_MM,
        "panel_sizes": {k: len(v) for k, v in panels.items()},
        "candidates": len(candidates),
        "candidates_conserved": int(sm.conserved.sum()),
        "candidates_specific": int(sm.specific.sum()),
        "candidates_pass_both": int(sm.passes_both.sum()),
        "target_g1_conservation_mean": round(float(prof.mean()), 4),
        "target_g1_fully_conserved_positions": int((prof == 1.0).sum()),
    }
    (RESULTS / "run_info.json").write_text(json.dumps(info, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(info, indent=2))

    show = ["guide", "pam", "gc", "target_exact_with_pam", "target_conservation",
            "nontarget_min_mm_with_pam", "nontarget_ctxm_min_mm", "nontarget_bla_min_mm",
            "n_nontarget_le2mm", "conserved", "specific"]
    with pd.option_context("display.width", 250, "display.max_columns", 30):
        print(sm[show].to_string(index=False))


if __name__ == "__main__":
    main()
