#!/usr/bin/env python3
"""Figure 18: do the assemblies confirm the read-based calls? (Phase 6), run in env `map`.

Run after scripts/compare_liftoff_pav.py:
    python scripts/plot_liftoff.py

Input (read only): results/assembly/liftoff_vs_pav_summary.tsv
Output: results/figures/18_assembly_liftoff__100pct_stacked_bar_chart.png (300 dpi)
    (a) what the assembly shows for the genes the reads call ABSENT (found / elsewhere / partial / missing)
    (b) what the reads called the genes that are MISSING from the assembly (absent / uncertain / ambiguous / present)
Assembly classes are defined in scripts/compare_liftoff_pav.py. The Ler assembly is "Ler" on NCBI.
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
BLUE, ORANGE, GREY, LIGHT = "#2a78d6", "#eb6834", "#898781", "#c3c2b7"
INK, INK2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e1e0d9", "#fcfcfb"
# (a) assembly classes, from "confirms the absent call" (dark) to "contradicts it" (orange)
A_CLASSES = [("missing", "not placed (missing)", "#1c4f8f"), ("partial", "under half of the gene aligns", BLUE),
             ("elsewhere", "placed elsewhere (look-alike?)", "#8fb8ea"), ("found", "found intact in place", ORANGE)]
# (b) read calls, from "caught" (dark) to "missed" (orange)
B_CALLS = [("absent", "called absent", "#1c4f8f"), ("uncertain", "uncertain (low-MAPQ reads)", BLUE),
           ("ambiguous", "ambiguous", "#8fb8ea"), ("present", "called present", ORANGE)]

plt.rcParams.update({
    "font.size": 9, "axes.edgecolor": LIGHT, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "axes.titlesize": 10, "axes.titleweight": "bold", "axes.titlelocation": "left",
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": False,
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
})


def stacked(ax, table, parts, title, xlabel):
    ys = np.arange(len(ACCS))
    tot = table.sum(axis=1).to_numpy()
    left = np.zeros(len(ACCS))
    for key, lab, col in parts:
        v = table[key].to_numpy()
        pct = 100 * v / tot
        ax.barh(ys, pct, left=left, color=col, height=0.62, label=lab)
        for yi, l, p, n in zip(ys, left, pct, v):
            if p >= 6:
                ax.text(l + p / 2, yi, f"{n}", ha="center", va="center", fontsize=8,
                        color="white" if col in ("#1c4f8f", BLUE, ORANGE) else INK)
        left += pct
    ax.set_yticks(ys)
    ax.set_yticklabels([f"{a}\n(n={t})" for a, t in zip(ACCS, tot)])
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=2, fontsize=8)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--summary", default="results/assembly/liftoff_vs_pav_summary.tsv")
    ap.add_argument("--out", default="results/figures/18_assembly_liftoff__100pct_stacked_bar_chart.png")
    a = ap.parse_args()

    s = pd.read_csv(a.summary, sep="\t")
    t = s.pivot_table(index=["accession", "call"], columns="asm_class", values="genes", fill_value=0, aggfunc="sum")
    ta = t.xs("absent", level="call").reindex(index=ACCS, columns=[k for k, _, _ in A_CLASSES], fill_value=0)
    tb = t["missing"].unstack("call").reindex(index=ACCS, columns=[k for k, _, _ in B_CALLS], fill_value=0)
    for acc in ACCS:
        print(f"{acc:<7} read-absent {int(ta.loc[acc].sum())}: confirmed (not found intact) "
              f"{100 * (1 - ta.loc[acc, 'found'] / ta.loc[acc].sum()):.0f}% | missing in assembly "
              f"{int(tb.loc[acc].sum())}: called present {int(tb.loc[acc, 'present'])} "
              f"({100 * tb.loc[acc, 'present'] / tb.loc[acc].sum():.0f}%)")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 4.6), gridspec_kw={"wspace": 0.32})
    stacked(ax1, ta, A_CLASSES, "(a) Genes the reads call absent: what the assembly shows",
            "% of read-absent genes")
    stacked(ax2, tb, B_CALLS, "(b) Genes missing from the assembly: what the reads called them",
            "% of genes not placed by Liftoff")
    fig.suptitle("Checking the read-based calls against 4 PacBio assemblies (Liftoff, gene spans)",
                 x=0.01, ha="left", fontsize=11)
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    fig.savefig(a.out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
