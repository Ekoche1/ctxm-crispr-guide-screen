"""Step 4: figure summarising the guide screen.

Panel A  conservation along blaCTX-M-15 (target group 1 vs other CTX-M groups)
         with each candidate guide drawn as a bar at its position.
Panel B  per-guide trade-off: target conservation vs fewest non-target mismatches.

Colours: Okabe-Ito blue / vermilion (colourblind-safe) + grey. Pass/fail is
also encoded by marker shape and labels, not colour alone.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
R = ROOT / "results"
BLUE, VERM, GREY, INK = "#0072B2", "#D55E00", "#8C8C8C", "#222222"

sm = pd.read_csv(R / "guide_summary.tsv", sep="\t")
prof = pd.read_csv(R / "conservation_profile.tsv", sep="\t")

plt.rcParams.update({"font.size": 10, "axes.edgecolor": GREY, "axes.labelcolor": INK,
                     "xtick.color": INK, "ytick.color": INK, "axes.spines.top": False,
                     "axes.spines.right": False})
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 9), gridspec_kw={"height_ratios": [1.1, 1]})

# ---- Panel A
WIN = 20  # same length as a guide spacer
smooth = lambda s: s.rolling(WIN, center=True, min_periods=1).mean()
ax1.plot(prof.position_1based, smooth(prof.conservation_target_g1), color=BLUE, lw=1.8,
         label=f"Target: CTX-M group 1 alleles (n=133)")
ax1.plot(prof.position_1based, smooth(prof.identity_other_ctxm_groups), color=GREY, lw=1.8,
         label=f"Non-target: other CTX-M groups (n=123)")
for i, r in sm.iterrows():
    row = i % 3  # stagger rows so neighbouring guides do not overlap
    y = 0.44 - 0.05 * row
    col = BLUE if r.passes_both else GREY
    ax1.plot([r.start_1based, r.end], [y, y], color=col, lw=5, solid_capstyle="butt")
ax1.set_ylim(0.26, 1.03)
ax1.set_yticks([0.5, 0.75, 1.0])
ax1.set_xlim(0, 876)
ax1.set_xlabel("Position along blaCTX-M-15 coding sequence (bp)")
ax1.set_ylabel(f"Fraction of panel with same base\n({WIN}-bp rolling mean)")
ax1.set_title("A. Conservation along blaCTX-M-15 and positions of the 20 candidate guides",
              loc="left", pad=30)
ax1.text(876, 0.27,"bars below: candidate guides (blue = passes both screens, grey = does not)",
         ha="right", va="bottom", fontsize=8.5, color=INK)
ax1.legend(loc="lower left", bbox_to_anchor=(0.0, 1.0), ncol=2, frameon=False, fontsize=9)

# ---- Panel B
for _, r in sm.iterrows():
    ok = r.passes_both
    ax2.scatter(r.target_conservation * 100, r.nontarget_min_mm_with_pam,
                s=70, color=BLUE if ok else VERM, marker="o" if ok else "X",
                edgecolor="white", linewidth=1.2, zorder=3)
    label = r.guide.split("_")[1]
    # -205 and +209 are neighbours: put each label on its own side of its own point
    offset = {"-205": (-34, 5), "+209": (6, -12)}.get(label, (6, 4))
    ax2.annotate(label, (r.target_conservation * 100, r.nontarget_min_mm_with_pam),
                 textcoords="offset points", xytext=offset, fontsize=8, color=INK)
ax2.axvline(95, color=GREY, ls="--", lw=1)
ax2.axhline(3, color=GREY, ls="--", lw=1)
ax2.text(95.2, ax2.get_ylim()[1], " >=95% conserved", va="top", fontsize=8.5, color=INK)
ax2.set_xlabel("Target conservation: % of group-1 alleles with exact match and valid PAM")
ax2.set_ylabel("Fewest mismatches to any non-target\n(valid PAM)")
ax2.set_title("B. Conservation vs specificity (circle = passes both; cross = fails at least one)", loc="left")
ax2.text(ax2.get_xlim()[0] + 0.5, 3.15, "3 mismatches", fontsize=8.5, color=INK)

fig.text(0.01, 0.005, "In-silico screen of reference alleles (AMRFinderPlus DB 2026-08-07.1). "
         "Not a measure of diagnostic sensitivity or specificity.", fontsize=8, color=INK)
fig.tight_layout(rect=(0, 0.02, 1, 1))
out = R / "figure1_guide_screen.png"
fig.savefig(out, dpi=170)
print("wrote", out)
