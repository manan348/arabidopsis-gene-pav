#!/usr/bin/env python3
"""Figure 22: where are the absent genes on the chromosomes, and where are the likely errors? (Phase 6), env `map`.

Run after scripts/compare_liftoff_pav.py:
    python scripts/plot_absent_map.py

Input (read only): results/assembly/liftoff_vs_pav_genes.tsv.gz
Output: results/figures/22_validation_absent_genes__chromosome_ideogram.png (300 dpi)

For the 4 accessions with an assembly, every gene the reads call absent is drawn at its Col-0 position:
  confirmed = the assembly does not have it intact (missing, partial or only a look-alike elsewhere)
  likely read error = Liftoff finds it intact at its own place
Centromere positions are approximate TAIR10 values (used only as landmarks, +-2 Mb band shaded).
Descriptive only: no statistical test.
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ACCS = ["An-1", "Cvi-0", "Kyoto", "Ler-1"]
CHR_LEN = {"1": 30427671, "2": 19698289, "3": 23459830, "4": 18585056, "5": 26975502}
CEN = {"1": 15086045, "2": 3607929, "3": 13799418, "4": 3956021, "5": 11725024}  # approximate (TAIR10)
BLUE, ORANGE, GREY, LIGHT = "#2a78d6", "#eb6834", "#898781", "#c3c2b7"
INK, INK2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e1e0d9", "#fcfcfb"

plt.rcParams.update({
    "font.size": 9, "axes.edgecolor": LIGHT, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "axes.titlesize": 10, "axes.titleweight": "bold", "axes.titlelocation": "left",
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": False,
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
})


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--liftoff", default="results/assembly/liftoff_vs_pav_genes.tsv.gz")
    ap.add_argument("--out", default="results/figures/22_validation_absent_genes__chromosome_ideogram.png")
    ap.add_argument("--cen-band", type=float, default=2e6)
    a = ap.parse_args()

    d = pd.read_csv(a.liftoff, sep="\t", dtype={"ref_chr": str})
    ra = d[d["call"] == "absent"].copy()
    ra["error"] = ra["asm_class"] == "found"
    ra["near_cen"] = (ra["ref_start"] - ra["ref_chr"].map(CEN)).abs() <= a.cen_band
    allg = d[d["accession"] == ACCS[0]]
    base = ((allg["ref_start"] - allg["ref_chr"].map(CEN)).abs() <= a.cen_band).mean()
    print(f"share of ALL genes within {a.cen_band / 1e6:.0f} Mb of a centromere: {100 * base:.1f}%")
    for acc in ACCS:
        x = ra[ra["accession"] == acc]
        c, e = x[~x["error"]], x[x["error"]]
        print(f"{acc:<6} confirmed {len(c):3d} ({100 * c['near_cen'].mean():.0f}% near a centromere) | "
              f"likely errors {len(e):2d} ({100 * e['near_cen'].mean():.0f}% near a centromere)")

    fig, axes = plt.subplots(5, 1, figsize=(12, 8.2), sharex=True, gridspec_kw={"hspace": 0.35, "top": 0.94})
    for ax, ch in zip(axes, ["1", "2", "3", "4", "5"]):
        L = CHR_LEN[ch] / 1e6
        ax.axvspan((CEN[ch] - a.cen_band) / 1e6, (CEN[ch] + a.cen_band) / 1e6, color=GRID, lw=0, zorder=0)
        for i, acc in enumerate(ACCS):
            y = i
            ax.plot([0, L], [y, y], color=LIGHT, lw=1, zorder=1)
            x = ra[(ra["accession"] == acc) & (ra["ref_chr"] == ch)]
            c, e = x[~x["error"]], x[x["error"]]
            ax.scatter(c["ref_start"] / 1e6, np.full(len(c), y), marker="|", s=60, color="#1c4f8f", lw=1, zorder=2,
                       label="absent, confirmed by the assembly" if (ch == "1" and i == 0) else None)
            ax.scatter(e["ref_start"] / 1e6, np.full(len(e), y), marker="o", s=22, color=ORANGE, zorder=3,
                       edgecolor="white", lw=0.5,
                       label="absent by reads, but found intact (likely error)" if (ch == "1" and i == 0) else None)
        ax.set_yticks(range(len(ACCS)))
        ax.set_yticklabels(ACCS, fontsize=8)
        ax.set_ylim(len(ACCS) - 0.4, -0.6)
        ax.set_xlim(-0.3, 31)
        ax.set_title(f"Chromosome {ch}", fontsize=9, pad=3)
        ax.spines["left"].set_visible(False)
        ax.tick_params(axis="y", length=0)
    axes[-1].set_xlabel("position on the Col-0 chromosome (Mb); grey band = approximate centromere +-2 Mb")
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, frameon=False, ncol=2, loc="lower center", bbox_to_anchor=(0.5, -0.02), fontsize=8,
               markerscale=1.3)
    fig.suptitle("Genes the reads call absent, by position (4 accessions with an assembly)",
                 x=0.01, ha="left", fontsize=11, y=0.995)
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    fig.savefig(a.out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
