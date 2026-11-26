#!/usr/bin/env python3
"""Figure 33: do gene losses group the accessions the same way SNPs do? Tanglegram + Mantel test, env `map`.

Run:
    python scripts/plot_pav_vs_snp_tree.py

Inputs (read only): results/pav/pav_long.tsv.gz, results/analysis/snp_distance_matrix.tsv
Outputs: results/figures/33_relatedness_pav_vs_snp__tanglegram_mantel_test.png (300 dpi), results/analysis/pav_distance_matrix.tsv

PAV distance between two accessions = number of genes called absent in exactly one of them.
SNP distance = share of differing sites (fig 31). Both trees: UPGMA (average linkage).
Mantel test: Pearson and Spearman correlation of the 28 pairwise distances; the p-value is exact, from all
8! = 40,320 relabellings of one matrix (one-sided, correlation >= observed). Repeated without Col-0 (7! = 5,040)
because Col-0 is the reference: it has 0 absent genes by construction and its SNP distances are biased.
With 8 accessions the test has little power; a high correlation is suggestive, not proof.
"""
import argparse
import itertools
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform
from scipy.stats import rankdata

ORDER = ["Col-0", "An-1", "Kyoto", "Kondara", "Cvi-0", "Ler-1", "Se-0", "Ale-Stenar-44-4"]
COL = {"Col-0": "#52514e", "An-1": "#2a78d6", "Ler-1": "#1c4f8f", "Ale-Stenar-44-4": "#7fb3ea",
       "Se-0": "#e3b23c", "Cvi-0": "#eb6834", "Kondara": "#8e6bbf", "Kyoto": "#2a9d8f"}
LIGHT, INK, INK2, GRID, SURFACE = "#c3c2b7", "#0b0b0b", "#52514e", "#e1e0d9", "#fcfcfb"
plt.rcParams.update({
    "font.size": 9, "axes.edgecolor": LIGHT, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "axes.titlesize": 10, "axes.titleweight": "bold", "axes.titlelocation": "left",
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
})


def mantel_exact(A, B, spearman=False):
    n = len(A)
    iu = np.triu_indices(n, 1)
    a, b = A[iu], B[iu]
    if spearman:
        a, b = rankdata(a), rankdata(b)
    r_obs = np.corrcoef(a, b)[0, 1]
    perms = np.array(list(itertools.permutations(range(n))))
    # B relabelled by each permutation, upper triangle only
    Bp = B[perms[:, iu[0]], perms[:, iu[1]]]
    if spearman:
        Bp = np.apply_along_axis(rankdata, 1, Bp)
    a0 = (a - a.mean()) / a.std()
    bz = (Bp - Bp.mean(axis=1, keepdims=True)) / Bp.std(axis=1, keepdims=True)
    r = (bz * a0).mean(axis=1)
    return r_obs, float((r >= r_obs - 1e-12).mean()), len(perms)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pav", default="results/pav/pav_long.tsv.gz")
    ap.add_argument("--snp", default="results/analysis/snp_distance_matrix.tsv")
    ap.add_argument("--matrix", default="results/analysis/pav_distance_matrix.tsv")
    ap.add_argument("--out", default="results/figures/33_relatedness_pav_vs_snp__tanglegram_mantel_test.png")
    a = ap.parse_args()

    p = pd.read_csv(a.pav, sep="\t", usecols=["gene", "accession", "call"])
    X = (p.pivot(index="gene", columns="accession", values="call") == "absent")[ORDER].to_numpy()
    P = np.array([[float((X[:, i] != X[:, j]).sum()) for j in range(8)] for i in range(8)])
    S = pd.read_csv(a.snp, sep="\t", index_col=0).loc[ORDER, ORDER].to_numpy()
    pd.DataFrame(P.astype(int), index=ORDER, columns=ORDER).to_csv(a.matrix, sep="\t")

    res = {}
    for label, idx in [("all 8", list(range(8))), ("without Col-0", list(range(1, 8)))]:
        A, B = S[np.ix_(idx, idx)], P[np.ix_(idx, idx)]
        rp, pp, n = mantel_exact(A, B)
        rs, ps, _ = mantel_exact(A, B, spearman=True)
        res[label] = (rp, pp, rs, ps, n)
        print(f"Mantel {label:<14} Pearson r = {rp:.3f} (p = {pp:.4f}) | Spearman rho = {rs:.3f} (p = {ps:.4f}) "
              f"| {n:,} permutations")

    fig = plt.figure(figsize=(12, 5.6))
    axl = fig.add_axes([0.03, 0.14, 0.3, 0.7])
    axm = fig.add_axes([0.33, 0.14, 0.3, 0.7])
    axr = fig.add_axes([0.63, 0.14, 0.3, 0.7])
    dl = dendrogram(linkage(squareform(S, checks=False), "average"), orientation="left", ax=axl, no_labels=True,
                    color_threshold=0, above_threshold_color=INK2)
    dr = dendrogram(linkage(squareform(P, checks=False), "average"), orientation="right", ax=axr, no_labels=True,
                    color_threshold=0, above_threshold_color=INK2)
    for ax in (axl, axr):
        ax.set_ylim(0, 80)
        for s in ["top", "right", "left"]:
            ax.spines[s].set_visible(False)
        ax.tick_params(axis="y", length=0)
        ax.grid(axis="x", color=GRID, lw=0.6)
        ax.set_axisbelow(True)
    axr.spines["right"].set_visible(False)
    axl.set_xlabel("SNP distance (share of differing sites)")
    axr.set_xlabel("PAV distance (genes absent in only one of the two)")
    axl.set_title("(a) SNP tree")
    axr.set_title("(b) gene presence/absence tree", loc="right")

    yl = {ORDER[i]: 5 + 10 * k for k, i in enumerate(dl["leaves"])}
    yr = {ORDER[i]: 5 + 10 * k for k, i in enumerate(dr["leaves"])}
    axm.set_xlim(0, 1)
    axm.set_ylim(0, 80)
    axm.axis("off")
    for acc in ORDER:
        c = COL[acc]
        short = acc.replace("Ale-Stenar-44-4", "Ale-S.-44-4")
        axm.text(0.02, yl[acc], short, ha="left", va="center", fontsize=9, color=c, fontweight="bold")
        axm.text(0.98, yr[acc], short, ha="right", va="center", fontsize=9, color=c, fontweight="bold")
        axm.plot([0.3, 0.7], [yl[acc], yr[acc]], color=c, lw=1.6, alpha=0.85)
    rp, pp, rs, ps, n = res["all 8"]
    rp7, pp7, rs7, ps7, n7 = res["without Col-0"]
    fig.suptitle("Do gene losses group the accessions the same way SNPs do?", x=0.03, y=0.985, ha="left",
                 fontsize=12, fontweight="bold")
    fig.text(0.03, 0.915, f"Mantel test (exact): all 8 accessions r = {rp:.2f} (p = {pp:.3f}), Spearman rho = "
             f"{rs:.2f} (p = {ps:.3f});  without Col-0 r = {rp7:.2f} (p = {pp7:.3f}), rho = {rs7:.2f} "
             f"(p = {ps7:.3f}).  UPGMA trees; lines join the same accession.", fontsize=8.5, color=INK2)
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    fig.savefig(a.out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("wrote", a.out, "and", a.matrix)


if __name__ == "__main__":
    main()
