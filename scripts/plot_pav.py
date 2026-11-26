#!/usr/bin/env python3
"""First PAV figures (Phase 4). Run from the project root after scripts/build_pav.py:
    python scripts/plot_pav.py [--pavdir results/pav] [--outdir results/figures]

Writes three PNG files:
  06_pav_calls_per_accession__grouped_bar_chart.png   absent / ambiguous / uncertain gene counts per accession
  04_pav_gene_depth__histograms.png  per-gene depth (MAPQ >= 10) / sample median, with the 0.2 cut-off
  08_pav_variable_genes__binary_heatmap.png      genes with >= 1 absence (rows) x accessions (columns)
No clustering is done here (that needs scipy; planned for Phase 7).
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

ORDER = ["Col-0", "An-1", "Kyoto", "Kondara", "Cvi-0", "Ler-1", "Se-0", "Ale-Stenar-44-4"]
C_PRESENT, C_AMBIG, C_UNC, C_ABSENT = "#e6e6e6", "#fdb863", "#5e8fc7", "#b2182b"


def fig_calls(long, out):
    counts = long.groupby(["accession", "call"]).size().unstack(fill_value=0).reindex(ORDER)
    cols = [("absent", C_ABSENT), ("ambiguous", C_AMBIG), ("uncertain", C_UNC)]
    x = np.arange(len(ORDER))
    w = 0.27
    fig, ax = plt.subplots(figsize=(10, 4.8))
    for i, (name, colour) in enumerate(cols):
        bars = ax.bar(x + (i - 1) * w, counts[name], w, label=name, color=colour)
        ax.bar_label(bars, fontsize=7, padding=2)
    ax.set_xticks(x)
    ax.set_xticklabels(ORDER, rotation=30, ha="right")
    ax.set_ylabel("number of genes (of 27,444)")
    ax.set_title("Genes not called present, per accession (Col-0 = control)")
    ax.legend(frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def fig_depth(long, out):
    fig, axes = plt.subplots(2, 4, figsize=(13, 6), sharex=True, sharey=True)
    bins = np.linspace(0, 3, 61)
    for ax, s in zip(axes.ravel(), ORDER):
        d = long.loc[long["accession"] == s, "norm_depth"].clip(upper=3)
        ax.hist(d, bins=bins, color="#6a8caf")
        ax.axvline(0.2, color=C_ABSENT, lw=1.2, ls="--")
        ax.axvline(1.0, color="#333333", lw=0.8)
        ax.set_yscale("log")
        ax.set_title(f"{s} (median {int(long.loc[long['accession'] == s, 'median_depth'].iloc[0])}x)", fontsize=9)
        ax.spines[["top", "right"]].set_visible(False)
    for ax in axes[1]:
        ax.set_xlabel("gene depth / sample median (capped at 3)")
    for ax in axes[:, 0]:
        ax.set_ylabel("genes (log)")
    fig.suptitle("Per-gene depth relative to each sample's median (dashed red = 0.2 x median cut-off)", fontsize=11)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def fig_heatmap(long, out):
    m = long.pivot(index="gene", columns="accession", values="call")[ORDER]
    var = m[(m == "absent").any(axis=1)].copy()
    code = {"present": 0, "ambiguous": 1, "uncertain": 2, "absent": 3}
    num = var.apply(lambda col: col.map(code)).astype(int)
    # sort: most absences first, then by the pattern of absences so identical patterns sit together
    key = pd.DataFrame({"n": -(num == 3).sum(axis=1)})
    for c in ORDER:
        key[c] = -(num[c] == 3).astype(int)
    num = num.loc[key.sort_values(["n"] + ORDER).index]
    cmap = ListedColormap([C_PRESENT, C_AMBIG, C_UNC, C_ABSENT])
    fig, ax = plt.subplots(figsize=(6.5, 9))
    ax.imshow(num.values, aspect="auto", cmap=cmap, vmin=-0.5, vmax=3.5, interpolation="nearest")
    ax.set_xticks(range(len(ORDER)))
    ax.set_xticklabels(ORDER, rotation=40, ha="right")
    ax.set_yticks([])
    ax.set_ylabel(f"{len(num)} genes absent in at least one accession")
    ax.set_title("Variable genes (sorted by number of absences)")
    ax.legend(handles=[Patch(color=C_PRESENT, label="present"), Patch(color=C_AMBIG, label="ambiguous"),
                       Patch(color=C_UNC, label="uncertain"), Patch(color=C_ABSENT, label="absent")],
              loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pavdir", default="results/pav")
    ap.add_argument("--outdir", default="results/figures")
    a = ap.parse_args()
    long = pd.read_csv(os.path.join(a.pavdir, "pav_long.tsv.gz"), sep="\t")
    os.makedirs(a.outdir, exist_ok=True)
    fig_calls(long, os.path.join(a.outdir, "06_pav_calls_per_accession__grouped_bar_chart.png"))
    fig_depth(long, os.path.join(a.outdir, "04_pav_gene_depth__histograms.png"))
    fig_heatmap(long, os.path.join(a.outdir, "08_pav_variable_genes__binary_heatmap.png"))
    print("wrote 3 figures to", a.outdir)


if __name__ == "__main__":
    main()
