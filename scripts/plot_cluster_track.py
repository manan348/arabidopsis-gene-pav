#!/usr/bin/env python3
"""Figure 39: read depth over the Chr5 small-gene cluster (~14.3 Mb) in all 8 accessions, env `map`.

First export read depth (env `map`, samtools; MAPQ >= 10 as in the PAV calls):
    samtools depth -a -Q 10 -r NC_003076.8:14300001-14410000 \
        data/bam/{Col-0,An-1,Kyoto,Kondara,Cvi-0,Ler-1,Se-0,Ale-Stenar-44-4}/*.markdup.bam > results/assembly/fig29_depth.tsv
Then:
    python scripts/plot_cluster_track.py

Inputs (read only): results/assembly/fig29_depth.tsv, data/reference/genes.bed, results/pav/pav_long.tsv.gz,
    results/go/go_enrichment.tsv
Output: results/figures/39_function_chr5_cluster__coverage_track.png (300 dpi)
Depth is divided by each accession's median nuclear depth (the PAV normaliser). Orange bars under a track mark
genes the reads call absent in that accession. Highlighted genes = the "protein localization to cell surface"
GO term (the small unnamed cluster from fig 34/fig 35).
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

ORDER = ["Col-0", "An-1", "Kyoto", "Kondara", "Cvi-0", "Ler-1", "Se-0", "Ale-Stenar-44-4"]
MEDIAN = {"Col-0": 11, "An-1": 7, "Kyoto": 8, "Kondara": 9, "Cvi-0": 7, "Ler-1": 12, "Se-0": 9, "Ale-Stenar-44-4": 8}
CHROM, START, END = "NC_003076.8", 14300000, 14410000  # 0-based start, as in the BED
TERM = "protein localization to cell surface"
BLUE, ORANGE, TEAL, LIGHT = "#2a78d6", "#eb6834", "#2a9d8f", "#c3c2b7"
INK, INK2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e1e0d9", "#fcfcfb"
plt.rcParams.update({
    "font.size": 9, "axes.edgecolor": LIGHT, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "axes.titlesize": 10, "axes.titleweight": "bold", "axes.titlelocation": "left",
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
})


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--depth", default="results/assembly/fig29_depth.tsv")
    ap.add_argument("--bed", default="data/reference/genes.bed")
    ap.add_argument("--pav", default="results/pav/pav_long.tsv.gz")
    ap.add_argument("--go", default="results/go/go_enrichment.tsv")
    ap.add_argument("--bin", type=int, default=250)
    ap.add_argument("--out", default="results/figures/39_function_chr5_cluster__coverage_track.png")
    a = ap.parse_args()

    dp = pd.read_csv(a.depth, sep="\t", header=None)
    assert dp.shape[1] == 2 + len(ORDER) and len(dp) == END - START, "depth file does not match the region"
    dp.columns = ["chrom", "pos"] + ORDER
    binned = dp.groupby((dp["pos"] - 1 - START) // a.bin)[ORDER].mean()
    x = (START + binned.index * a.bin + a.bin / 2) / 1e6

    bed = pd.read_csv(a.bed, sep="\t", header=None, usecols=[0, 1, 2, 3], names=["chrom", "s", "e", "gene"])
    genes = bed[(bed["chrom"] == CHROM) & (bed["e"] > START) & (bed["s"] < END)].drop_duplicates("gene")
    go = pd.read_csv(a.go, sep="\t")
    cl = set(go[(go["set"] == "absent") & (go["Description"] == TERM)].iloc[0]["geneID"].split("/"))
    p = pd.read_csv(a.pav, sep="\t", usecols=["gene", "accession", "call"])
    p = p[p["gene"].isin(genes["gene"])]
    absent = set(zip(p.loc[p["call"] == "absent", "gene"], p.loc[p["call"] == "absent", "accession"]))
    hl = genes[genes["gene"].isin(cl)]
    print(f"{len(genes)} genes in the window, {len(hl)} from the cluster term: {', '.join(hl['gene'])}")
    for acc in ORDER:
        n = sum((g, acc) in absent for g in hl["gene"])
        print(f"  {acc:<16} cluster genes called absent: {n} of {len(hl)}")

    n = len(ORDER)
    fig, axes = plt.subplots(n + 1, 1, figsize=(11, 9.5), sharex=True,
                             gridspec_kw={"height_ratios": [1] * n + [1.1], "hspace": 0.18})
    for ax, acc in zip(axes, ORDER):
        ax.fill_between(x, 0, (binned[acc] / MEDIAN[acc]).clip(upper=3), step="mid", color=BLUE, alpha=0.8, lw=0)
        for _, g in hl.iterrows():
            ax.axvspan(g["s"] / 1e6, g["e"] / 1e6, color=TEAL, alpha=0.18, lw=0, zorder=0)
        for _, g in genes.iterrows():
            if (g["gene"], acc) in absent:
                ax.plot([g["s"] / 1e6, g["e"] / 1e6], [-0.25, -0.25], color=ORANGE, lw=3.5, solid_capstyle="butt")
        ax.axhline(1, color=INK2, lw=0.6, ls=":")
        ax.set_ylim(-0.45, 3)
        ax.set_yticks([0, 1, 2, 3])
        ax.tick_params(axis="y", labelsize=7)
        ax.set_ylabel(acc.replace("Ale-Stenar-44-4", "Ale-S.-44-4"), rotation=0, ha="right", va="center",
                      fontsize=8.5)
        ax.tick_params(axis="x", length=0)
    axes[0].set_title("Read depth relative to each accession's median (MAPQ >= 10, capped at 3x); "
                      "orange bar = gene called absent; teal band = cluster gene")

    ax = axes[-1]
    k = 0
    for _, g in genes.sort_values("s").iterrows():
        is_cl = g["gene"] in cl
        s, e = max(g["s"], START) / 1e6, min(g["e"], END) / 1e6
        ax.add_patch(plt.Rectangle((s, 0.62), e - s, 0.28, color=TEAL if is_cl else LIGHT, lw=0))
        if is_cl:  # label only the cluster genes, on alternating rows
            ax.text((s + e) / 2, [0.02, 0.21, 0.40][k % 3], g["gene"], ha="center", va="bottom", fontsize=6.8,
                    color=TEAL, fontweight="bold")
            k += 1
    ax.set_ylim(0, 1.0)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xlim(START / 1e6, END / 1e6)
    ax.set_xlabel("position on Col-0 chromosome 5 (Mb); grey = other genes, teal = cluster genes")
    fig.suptitle("The Chr5 small-gene cluster (~14.3 Mb) is a multi-copy region: even Col-0 has few uniquely "
                 "mapped reads", x=0.01, ha="left", fontsize=12, fontweight="bold", y=0.99)
    fig.text(0.01, 0.955,"In Col-0 all-read depth is normal (MAPQ 0: 100% covered) but MAPQ >= 10 depth is low, so "
             "the calls there are 'uncertain/ambiguous'; 'absent' (orange) means no reads at all, even at MAPQ 0.",
             fontsize=8.5, color=INK2)
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    fig.savefig(a.out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
