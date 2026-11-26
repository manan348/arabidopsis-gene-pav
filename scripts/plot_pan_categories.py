#!/usr/bin/env python3
"""Figure 27: pan-genome categories (core / soft-core / shell / private) of the 27,444 genes, env `map`.

Run:
    python scripts/plot_pan_categories.py

Input (read only): results/pav/pav_long.tsv.gz.  Output: results/figures/27_pangenome_categories__donut_chart_histogram.png (300 dpi)

A gene counts as present in an accession unless the reads call it `absent` (uncertain and ambiguous calls count as
present, the same rule as the 26,533-gene core in fig 29). Categories by the number of accessions with the gene:
  core 8/8, soft-core 7/8, shell 2-6/8, private 1/8 (with Col-0 as the reference, private = Col-0 only).
Only genes annotated in Col-0 can be counted: genes found only in other accessions are invisible to this method, so
the shell and private shares are lower bounds. Descriptive only: no statistical test.
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

LIGHT, INK, INK2, GRID, SURFACE = "#c3c2b7", "#0b0b0b", "#52514e", "#e1e0d9", "#fcfcfb"
CAT = [("core", "8 of 8", "#1c4f8f"), ("soft-core", "7 of 8", "#2a78d6"), ("shell", "2-6 of 8", "#eb6834"),
       ("private", "1 of 8", "#8e2c0c")]
plt.rcParams.update({
    "font.size": 9, "axes.edgecolor": LIGHT, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "axes.titlesize": 10, "axes.titleweight": "bold", "axes.titlelocation": "left",
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
})


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pav", default="results/pav/pav_long.tsv.gz")
    ap.add_argument("--out", default="results/figures/27_pangenome_categories__donut_chart_histogram.png")
    a = ap.parse_args()

    p = pd.read_csv(a.pav, sep="\t", usecols=["gene", "accession", "call"])
    n_acc = p["accession"].nunique()
    assert n_acc == 8
    n_with = n_acc - (p["call"] == "absent").groupby(p["gene"]).sum()
    hist = n_with.value_counts().reindex(range(1, 9), fill_value=0)
    cats = {"core": hist[8], "soft-core": hist[7], "shell": hist[2:6].sum(), "private": hist[1]}
    total = int(hist.sum())
    assert total == 27444
    for k, v in cats.items():
        print(f"{k:<10} {v:6,d}  {100 * v / total:6.2f}%")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.8), gridspec_kw={"width_ratios": [1, 1.35], "wspace": 0.3})
    vals = [cats[c] for c, _, _ in CAT]
    cols = [c for _, _, c in CAT]
    ax1.pie(vals, colors=cols, startangle=90, counterclock=False, wedgeprops={"width": 0.36, "edgecolor": SURFACE})
    ax1.text(0, 0.08, f"{total:,}", ha="center", va="center", fontsize=15, fontweight="bold")
    ax1.text(0, -0.14, "genes in Col-0", ha="center", va="center", fontsize=8.5, color=INK2)
    labels = [f"{c} ({s}): {v:,} ({100 * v / total:.2f}%)" for (c, s, _), v in zip(CAT, vals)]
    ax1.legend([plt.Rectangle((0, 0), 1, 1, color=c) for c in cols], labels, frameon=False, loc="upper center",
               bbox_to_anchor=(0.5, 0.02), fontsize=8.5)
    ax1.set_title("(a) Gene categories")
    ax1.set_aspect("equal")

    x = np.arange(1, 9)
    col_of = {1: CAT[3][2], 7: CAT[1][2], 8: CAT[0][2]}
    bars = ax2.bar(x, hist.values, color=[col_of.get(i, CAT[2][2]) for i in x], width=0.7)
    ax2.set_yscale("log")
    ax2.set_ylim(5, hist.max() * 4)
    for b, v in zip(bars, hist.values):
        ax2.text(b.get_x() + b.get_width() / 2, v * 1.15, f"{v:,}", ha="center", fontsize=8, color=INK2)
    ax2.set_xticks(x)
    ax2.set_xlabel("number of accessions with the gene (of 8)")
    ax2.set_ylabel("number of genes (log scale)")
    ax2.set_title("(b) How many accessions carry each gene")
    ax2.grid(axis="y", color=GRID, lw=0.6)
    ax2.set_axisbelow(True)
    fig.suptitle("Pan-genome categories of the Col-0 gene set (present = not called absent by reads)",
                 x=0.01, ha="left", fontsize=11)
    fig.text(0.01, -0.04, "Genes found only in other accessions cannot be seen with a Col-0 reference, so shell and "
             "private counts are lower bounds.", fontsize=7.5, color=INK2)
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    fig.savefig(a.out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
