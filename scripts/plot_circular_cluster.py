#!/usr/bin/env python3
"""Figure 30: circular cluster plot of the 911 variable genes, env `map`.

Run (after scripts/go_enrichment.R):
    python scripts/plot_circular_cluster.py

Inputs (read only): results/pav/pav_long.tsv.gz, results/go/go_enrichment.tsv
Output: results/figures/30_pangenome_gene_patterns__circular_dendrogram_heatmap.png (300 dpi)

Centre: genes clustered by their absence pattern across the 7 non-reference accessions (Hamming distance on
absent yes/no, average linkage), so genes missing in the same accessions sit together.
Rings (inside to outside): the read call in each of the 7 accessions, the number of accessions missing the gene,
and the GO group (from the enriched terms, as in the PAV explorer). Col-0 is left out: it has no absent genes.
Descriptive only: no statistical test.
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap
from scipy.cluster.hierarchy import dendrogram, linkage

ACCS = ["An-1", "Kyoto", "Kondara", "Cvi-0", "Ler-1", "Se-0", "Ale-Stenar-44-4"]
CALL_COL = {"present": "#e6edf6", "ambiguous": "#c3c2b7", "uncertain": "#f3c9a8", "absent": "#eb6834"}
GROUPS = [("ADP binding", "NLR resistance gene", "#2a78d6"),
          ("NAD+ nucleosidase activity, cyclic ADP-ribose generating", "NLR resistance gene", "#2a78d6"),
          ("defense response to other organism", "other defense gene", "#9ec5f0"),
          ("sesquiterpene biosynthetic process", "terpene synthase", "#e3b23c"),
          ("RNA-DNA hybrid ribonuclease activity", "RNase H (TE-like)", "#8e6bbf"),
          ("protein localization to cell surface", "small unnamed cluster", "#2a9d8f")]
OTHER = "#ecebe6"
INK, INK2, SURFACE = "#0b0b0b", "#52514e", "#fcfcfb"
GAP = np.deg2rad(40)          # empty wedge at the top for ring labels
R_LEAF, RING_W = 1.0, 0.075   # dendrogram leaves at r = 1; each ring this thick


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pav", default="results/pav/pav_long.tsv.gz")
    ap.add_argument("--go", default="results/go/go_enrichment.tsv")
    ap.add_argument("--out", default="results/figures/30_pangenome_gene_patterns__circular_dendrogram_heatmap.png")
    a = ap.parse_args()

    p = pd.read_csv(a.pav, sep="\t", usecols=["gene", "accession", "call"])
    calls = p.pivot(index="gene", columns="accession", values="call")
    calls = calls[(calls == "absent").any(axis=1)][ACCS]
    assert len(calls) == 911
    X = (calls == "absent").to_numpy(dtype=float)
    Z = linkage(X, method="average", metric="hamming")
    dn = dendrogram(Z, no_plot=True)
    order = dn["leaves"]
    n = len(order)
    genes = calls.index[order]

    go = pd.read_csv(a.go, sep="\t")
    go = go[go["set"] == "absent"]
    grp = {}
    for term, label, col in GROUPS:
        row = go[go["Description"] == term]
        if len(row):
            for g in row.iloc[0]["geneID"].split("/"):
                grp.setdefault(g, (label, col))

    span = 2 * np.pi - GAP
    theta = lambda x: np.pi / 2 - GAP / 2 - (np.asarray(x) - 5) / (10 * n) * span  # clockwise from the gap
    leaf_t = theta(np.arange(n) * 10 + 5)
    width = span / n
    dmax = max(max(d) for d in dn["dcoord"])
    rad = lambda d: R_LEAF * (1 - np.asarray(d) / dmax) * 0.97 + 0.03

    fig = plt.figure(figsize=(10, 10.6))
    ax = fig.add_axes([0.05, 0.1, 0.9, 0.85], projection="polar")
    ax.set_axis_off()
    for xs, ds in zip(dn["icoord"], dn["dcoord"]):
        t, r = theta(xs), rad(ds)
        ax.plot([t[0], t[1]], [r[0], r[1]], color=INK2, lw=0.35)
        ax.plot([t[3], t[2]], [r[3], r[2]], color=INK2, lw=0.35)
        arc = np.linspace(t[1], t[2], 30)
        ax.plot(arc, np.full_like(arc, r[1]), color=INK2, lw=0.35)

    r0 = R_LEAF + 0.04
    for i, acc in enumerate(ACCS):
        cols = [CALL_COL[c] for c in calls.loc[genes, acc]]
        ax.bar(leaf_t, RING_W * 0.92, width=width, bottom=r0 + i * RING_W, color=cols, linewidth=0)
        ax.text(np.pi / 2, r0 + (i + 0.45) * RING_W, acc.replace("Ale-Stenar-44-4", "Ale-S.-44-4"), ha="center",
                va="center", fontsize=6.5, color=INK2)
    rn = r0 + len(ACCS) * RING_W + 0.03
    nmiss = X[order].sum(axis=1).astype(int)
    cmap = ListedColormap(plt.cm.YlOrRd(np.linspace(0.15, 0.95, 7)))
    ax.bar(leaf_t, RING_W * 0.92, width=width, bottom=rn, color=cmap((nmiss - 1) / 6), linewidth=0)
    ax.text(np.pi / 2, rn + 0.45 * RING_W, "missing in", ha="center", va="center", fontsize=6.5, color=INK2)
    rg = rn + RING_W + 0.03
    gcols = [grp.get(g, ("other", OTHER))[1] for g in genes]
    ax.bar(leaf_t, RING_W * 1.3, width=width, bottom=rg, color=gcols, linewidth=0)
    ax.text(np.pi / 2, rg + 0.65 * RING_W, "GO group", ha="center", va="center", fontsize=6.5, color=INK2)
    ax.set_ylim(0, rg + RING_W * 1.5)

    h1 = [plt.Rectangle((0, 0), 1, 1, color=c) for c in CALL_COL.values()]
    l1 = list(CALL_COL)
    seen = dict.fromkeys((l, c) for _, l, c in GROUPS)
    h2 = [plt.Rectangle((0, 0), 1, 1, color=c) for _, c in seen] + [plt.Rectangle((0, 0), 1, 1, color=OTHER)]
    l2 = [f"{l} ({sum(1 for g in genes if grp.get(g, ('',))[0] == l)})" for l, _ in seen] + \
         [f"other ({sum(1 for g in genes if g not in grp)})"]
    fig.legend(h1, l1, title="read call", loc="lower left", bbox_to_anchor=(0.04, 0.0), ncol=2, frameon=False,
               fontsize=8, title_fontsize=8.5)
    fig.legend(h2, l2, title="GO group (enriched terms)", loc="lower center", bbox_to_anchor=(0.52, 0.0), ncol=2,
               frameon=False, fontsize=8, title_fontsize=8.5)
    cax = fig.add_axes([0.80, 0.045, 0.15, 0.015])
    cb = plt.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0.5, 7.5)), cax=cax,
                      orientation="horizontal", ticks=range(1, 8))
    cb.ax.tick_params(labelsize=7)
    cb.set_label("accessions missing the gene", fontsize=8)
    cb.outline.set_visible(False)
    fig.suptitle("911 variable genes clustered by which accessions lack them", x=0.03, y=0.985, ha="left",
                 fontsize=12, fontweight="bold")
    fig.text(0.03, 0.952, "Centre: average-linkage tree of absence patterns (7 accessions; Col-0 has none). "
             "Rings: read call per accession, number of accessions missing the gene, GO group.",
             fontsize=8, color=INK2)
    print("genes per number of accessions missing:", np.bincount(nmiss, minlength=8)[1:].tolist())
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    fig.savefig(a.out, dpi=300, facecolor=SURFACE)
    plt.close(fig)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
