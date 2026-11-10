"""Step 7: build the summary tables used in the README, from result files only.

Reads  results/guide_summary.tsv, genome_pilot_gene.tsv, genome_pilot_guides.tsv,
       passing_guides_misses.tsv, passing_guides_hybrids.tsv
Writes results/genome_pilot_guide_summary.tsv, results/passing_guides.md
and replaces the block between the PASSING_TABLE markers in README.md, so no
number in that table is ever typed by hand.
"""

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
R = ROOT / "results"
START, END = "<!-- PASSING_TABLE_START -->", "<!-- PASSING_TABLE_END -->"
MIN_IDENTITY_POSITIVE = 0.99  # a genome counts as gene-positive only if the located gene is >=99% identical to CTX-M-15


def main() -> None:
    sm = pd.read_csv(R / "guide_summary.tsv", sep="\t")
    gene = pd.read_csv(R / "genome_pilot_gene.tsv", sep="\t")
    g = pd.read_csv(R / "genome_pilot_guides.tsv", sep="\t")

    gene["verified_positive"] = gene.found & (gene.identity_over_ref >= MIN_IDENTITY_POSITIVE)
    g = g.merge(gene[["accession", "verified_positive"]], on="accession")
    pos, neg = g[g.verified_positive], g[~g.verified_positive]
    n_pos, n_neg = pos.accession.nunique(), neg.accession.nunique()

    rows = []
    for guide, p in pos.groupby("guide"):
        n = neg[neg.guide == guide]
        rows.append({
            "guide": guide,
            "positives_with_exact_pam_hit": int((p.exact_copies_with_pam >= 1).sum()),
            "positives_total": n_pos,
            "total_exact_copies_in_positives": int(p.exact_copies_with_pam.sum()),
            "gene_free_min_mm_with_pam": int(n.min_mm_with_pam.min()),
            "gene_free_exact_hits": int(n.exact_copies_with_pam.sum()),
            "gene_free_total": n_neg,
        })
    gs = pd.DataFrame(rows)
    gs.to_csv(R / "genome_pilot_guide_summary.tsv", sep="\t", index=False, lineterminator="\n")

    misses = pd.read_csv(R / "passing_guides_misses.tsv", sep="\t").groupby("guide").size().rename("misses")
    hyb = pd.read_csv(R / "passing_guides_hybrids.tsv", sep="\t").set_index("guide")[["exact_hits"]]
    d = (sm[sm.passes_both].merge(gs, on="guide").set_index("guide").join(misses).join(hyb).reset_index())
    d["misses"] = d.misses.fillna(0).astype(int)

    lines = [
        "| Guide (strand, forward position) | PAM | Spacer 5'→3' (DNA) | GC | Group-1 alleles matched | Fewest mismatches to any non-target allele | "
        f"Fewest mismatches in {n_neg} gene-free chromosomes | Unassigned/hybrid alleles matched |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for r in d.itertuples():
        lines.append(
            f"| {r.guide.split('_')[1]} | {r.pam} | `{r.spacer_5to3}` | {r.gc:.2f} | "
            f"{r.target_exact_with_pam}/{r.target_n} ({r.target_conservation * 100:.1f}%) | "
            f"{r.nontarget_min_mm_with_pam} | {r.gene_free_min_mm_with_pam} | {r.exact_hits}/{r.unassigned_n} |")
    table = "\n".join(lines)
    (R / "passing_guides.md").write_text(table + "\n", encoding="utf-8")

    readme = ROOT / "README.md"
    text = readme.read_text(encoding="utf-8")
    if START not in text or END not in text:
        raise SystemExit("README.md is missing the PASSING_TABLE markers")
    head, rest = text.split(START, 1)
    _, tail = rest.split(END, 1)
    readme.write_text(f"{head}{START}\n{table}\n{END}{tail}", encoding="utf-8", newline="\n")
    print(table)
    print(f"\ngene-positive sequences: {n_pos}; gene-free chromosomes: {n_neg}")


if __name__ == "__main__":
    main()
