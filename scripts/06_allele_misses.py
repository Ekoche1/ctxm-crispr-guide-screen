"""Step 6: which target alleles does each passing guide miss?

Reads  results/guide_summary.tsv, results/guide_vs_sequence.tsv
Writes results/passing_guides_misses.tsv     (one row per guide x missed group-1 allele)
       results/passing_guides_hybrids.tsv    (exact hits among the 11 unassigned alleles)
A "miss" means no exact match with a valid PAM. A single-mismatch miss may or
may not still be cleaved in the lab; this study cannot say.
"""

from pathlib import Path

import pandas as pd

R = Path(__file__).resolve().parents[1] / "results"
long = pd.read_csv(R / "guide_vs_sequence.tsv", sep="\t")
sm = pd.read_csv(R / "guide_summary.tsv", sep="\t")
passing = sm[sm.passes_both].guide.tolist()

target = long[(long.panel == "target_g1") & long.guide.isin(passing)]
misses = target[target.min_mm_with_pam > 0][["guide", "sequence", "min_mm_with_pam"]]
misses.sort_values(["guide", "min_mm_with_pam", "sequence"]).to_csv(
    R / "passing_guides_misses.tsv", sep="\t", index=False, lineterminator="\n")

hit = target.assign(hit=target.min_mm_with_pam == 0).pivot(index="sequence", columns="guide", values="hit")
print("passing guides:", passing)
print("group-1 alleles missed by ALL passing guides:", list(hit.index[~hit.any(axis=1)]))
print("group-1 alleles hit by <= 2 passing guides:",
      {i: int(hit.loc[i].sum()) for i in hit.index[hit.sum(axis=1) <= 2]})

un = long[(long.panel == "unassigned_ctxm") & long.guide.isin(passing)]
hyb = un.assign(exact=un.min_mm_with_pam == 0).groupby("guide").agg(
    exact_hits=("exact", "sum"), unassigned_n=("exact", "size")).reset_index()
hyb.to_csv(R / "passing_guides_hybrids.tsv", sep="\t", index=False, lineterminator="\n")
print(hyb.to_string(index=False))
print("\nmisses per guide:\n", misses.groupby("guide").size().to_string())
