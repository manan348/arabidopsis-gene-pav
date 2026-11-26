#!/usr/bin/env python3
"""Figure 31: how close are the accessions to each other in their SNPs? (Phase 7), run in env `map`.

Run:
    python scripts/plot_snp_pca.py

Input (read only): results/variants/filtered_gt.tsv.gz
    columns: CHROM POS TYPE QUAL DP and then one genotype per accession (Col-0 ... Ale-Stenar-44-4)
Output: results/figures/31_relatedness_snps__pcoa_plot.png (300 dpi); results/variants/snp_distance_matrix.tsv

Method: only biallelic SNPs (genotype 0 or 1; everything else is treated as missing, because low depth
leaves many missing genotypes). For every pair of accessions, distance = share of sites at which the two
differ, among sites genotyped in both. The 8 x 8 distance matrix is turned into a PCA-like map by classical
multidimensional scaling (principal coordinates). With only 8 accessions this is a rough picture of
relatedness, not a population-genetic analysis. Col-0 is the reference, so every distance to it is the
divergence from the reference. Descriptive only: no statistical test.
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ORDER = ["Col-0", "An-1", "Kyoto", "Kondara", "Cvi-0", "Ler-1", "Se-0", "Ale-Stenar-44-4"]

BLUE, ORANGE, GREY, LIGHT = "#2a78d6", "#eb6834", "#898781", "#c3c2b7"
INK, INK2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e1e0d9", "#fcfcfb"

plt.rcParams.update({
    "font.size": 9, "axes.edgecolor": LIGHT, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "axes.titlesize": 10, "axes.titleweight": "bold", "axes.titlelocation": "left",
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": False,
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
})


def pair_distances(path, chunk=250_000):
    n = len(ORDER)
    diff = np.zeros((n, n))
    both = np.zeros((n, n))
    n_snp = 0
    for ch in pd.read_csv(path, sep="\t", header=None, dtype=str, chunksize=chunk):
        ch = ch[ch[2] == "SNP"]
        g = ch.iloc[:, 5:5 + n].to_numpy()
        m = ((g == "0") | (g == "1")).astype(np.float64)  # genotyped as ref or alt
        x = (g == "1").astype(np.float64)                  # 1 = alt allele (0 where missing)
        n_snp += len(ch)
        diff += x.T @ m + m.T @ x - 2 * x.T @ x            # sum over sites of m_i m_j (x_i - x_j)^2
        both += m.T @ m
    return diff / both, both, n_snp


def pcoa(D):
    """Classical multidimensional scaling of a distance matrix."""
    n = len(D)
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J @ (D ** 2) @ J
    w, v = np.linalg.eigh(B)
    idx = np.argsort(w)[::-1]
    w, v = w[idx], v[:, idx]
    pos = w > 1e-12
    coords = v[:, pos] * np.sqrt(w[pos])
    return coords, w[pos] / w[pos].sum()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gt", default="results/variants/filtered_gt.tsv.gz")
    ap.add_argument("--out", default="results/figures/31_relatedness_snps__pcoa_plot.png")
    ap.add_argument("--matrix", default="results/variants/snp_distance_matrix.tsv")
    a = ap.parse_args()

    D, both, n_snp = pair_distances(a.gt)
    print(f"SNP rows read: {n_snp:,}; sites genotyped in both accessions of a pair: "
          f"{int(both[np.triu_indices(len(ORDER), 1)].min()):,} to {int(both[np.triu_indices(len(ORDER), 1)].max()):,}")
    assert np.allclose(D, D.T) and np.allclose(np.diag(D), 0)
    os.makedirs(os.path.dirname(a.matrix) or ".", exist_ok=True)
    pd.DataFrame(D, index=ORDER, columns=ORDER).round(4).to_csv(a.matrix, sep="\t")
    print(pd.DataFrame(D, index=ORDER, columns=ORDER).round(3).to_string())

    coords, share = pcoa(D)
    print("share of variation: " + ", ".join(f"PC{i + 1} {100 * s:.1f}%" for i, s in enumerate(share[:4])))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 5), gridspec_kw={"width_ratios": [1, 1.05], "wspace": 0.28})
    for i, acc in enumerate(ORDER):
        col = GREY if acc == "Col-0" else BLUE
        ax1.scatter(coords[i, 0], coords[i, 1], s=70, color=col, zorder=3)
        dx, dy, ha = (-8, -14, "right") if acc == "Se-0" else (6, 5, "left")  # Se-0 and An-1 lie close together
        ax1.annotate(acc, (coords[i, 0], coords[i, 1]), xytext=(dx, dy), textcoords="offset points", fontsize=8.5,
                     color=INK2, ha=ha)
    span0 = np.ptp(coords[:, 0])
    span1 = np.ptp(coords[:, 1])
    ax1.set_xlim(coords[:, 0].min() - 0.12 * span0, coords[:, 0].max() + 0.35 * span0)
    ax1.set_ylim(coords[:, 1].min() - 0.15 * span1, coords[:, 1].max() + 0.15 * span1)
    ax1.set_xlabel(f"PC1 ({100 * share[0]:.0f}% of variation)")
    ax1.set_ylabel(f"PC2 ({100 * share[1]:.0f}% of variation)")
    ax1.set_title("(a) Accession map from SNP differences")
    ax1.grid(color=GRID, lw=0.6)
    ax1.set_axisbelow(True)

    im = ax2.imshow(D, cmap="Blues", vmin=0, vmax=D.max())
    for i in range(len(ORDER)):
        for j in range(len(ORDER)):
            ax2.text(j, i, f"{D[i, j]:.2f}", ha="center", va="center", fontsize=7.5,
                     color="white" if D[i, j] > 0.6 * D.max() else INK)
    ax2.set_xticks(range(len(ORDER)))
    ax2.set_xticklabels(ORDER, rotation=40, ha="right")
    ax2.set_yticks(range(len(ORDER)))
    ax2.set_yticklabels(ORDER)
    for sp in ax2.spines.values():
        sp.set_visible(False)
    ax2.set_title("(b) Share of SNP sites where two accessions differ")
    cb = fig.colorbar(im, ax=ax2, fraction=0.046, pad=0.03)
    cb.outline.set_visible(False)

    fig.suptitle("How close are the accessions? (biallelic SNPs, missing genotypes skipped pair by pair)",
                 x=0.01, ha="left", fontsize=11)
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    fig.savefig(a.out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
