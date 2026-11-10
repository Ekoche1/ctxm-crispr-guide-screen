"""Step 3: independently re-derive the screening results with a different engine.

The main analysis (compare.py) uses numpy sliding windows. Here we recompute
every guide x sequence minimum-mismatch value with the `regex` library's fuzzy
matching, writing the PAM *into the pattern*, and compare. We also re-derive
guide coordinates and PAM counts with plain string operations.

Mismatch values are capped at CAP (anything larger is reported as ">=CAP")
so the slow fuzzy search does not run to 20 mismatches.

Writes results/independent_check.json and exits non-zero on any disagreement.
"""

import json
import re
import sys
from pathlib import Path

import pandas as pd
import regex

from ctxm_screen.sequences import parse_fasta, reverse_complement

ROOT = Path(__file__).resolve().parents[1]
RESULTS, PANELS = ROOT / "results", ROOT / "data" / "panels"
CAP = 12


def fuzzy_min(guide: str, seq: str, with_pam: bool) -> int:
    """Smallest k <= CAP with a substitution-only match, on either strand."""
    pat_core = f"(?:{guide})"
    for k in range(CAP + 1):
        pattern = f"(?<=TTT[ACG])(?:{guide}){{s<={k}}}" if with_pam else f"(?:{guide}){{s<={k}}}"
        for strand in (seq, reverse_complement(seq)):
            if regex.search(pattern, strand, flags=regex.BESTMATCH):
                return k
    return CAP + 1


def main() -> int:
    summary = pd.read_csv(RESULTS / "guide_summary.tsv", sep="\t")
    long = pd.read_csv(RESULTS / "guide_vs_sequence.tsv", sep="\t")
    seqs = {}
    for f in PANELS.glob("*.fasta"):
        for r in parse_fasta(f):
            seqs[(f.stem, r.id)] = r.seq
    guide_seq = dict(zip(summary.guide, summary.spacer_5to3))

    # 1) every guide x sequence value, both with and without the PAM requirement
    disagreements = []
    for row in long.itertuples(index=False):
        s = seqs[(row.panel, row.sequence)]
        g = guide_seq[row.guide]
        for col, with_pam in (("min_mm_any_pam", False), ("min_mm_with_pam", True)):
            mine = min(getattr(row, col), CAP + 1)
            other = fuzzy_min(g, s, with_pam)
            if mine != other:
                disagreements.append((row.guide, row.panel, row.sequence, col, mine, other))
    n_pairs = len(long) * 2

    # 2) coordinates: spacer must equal the reference slice (or its reverse complement)
    ref = seqs[("target_g1", "blaCTX-M-15")]
    coord_bad = []
    for r in summary.itertuples(index=False):
        s0, e0 = r.start_1based - 1, r.end
        piece = ref[s0:e0]
        spacer_on_plus = piece if r.strand == "+" else reverse_complement(piece)
        if spacer_on_plus != r.spacer_5to3:
            coord_bad.append(r.guide)
        # the 4 nt before the protospacer on its own strand must be the stated PAM
        pam_region = ref[s0 - 4 : s0] if r.strand == "+" else reverse_complement(ref[e0 : e0 + 4])
        if pam_region != r.pam:
            coord_bad.append(r.guide + " (PAM)")

    # 3) number of TTTV sites: independent lookahead counts on each strand
    n_plus = len(re.findall(r"(?=TTT[ACG].{20})", ref))
    n_minus = len(re.findall(r"(?=TTT[ACG].{20})", reverse_complement(ref)))

    out = {
        "pairs_compared": n_pairs, "cap": CAP,
        "disagreements": len(disagreements),
        "first_disagreements": disagreements[:10],
        "coordinate_or_pam_errors": coord_bad,
        "independent_candidate_count_plus": n_plus,
        "independent_candidate_count_minus": n_minus,
        "independent_candidate_total": n_plus + n_minus,
        "pipeline_candidate_total": int(len(summary)),
    }
    (RESULTS / "independent_check.json").write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(out, indent=2))
    ok = not disagreements and not coord_bad and out["independent_candidate_total"] == out["pipeline_candidate_total"]
    print("INDEPENDENT CHECK:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
