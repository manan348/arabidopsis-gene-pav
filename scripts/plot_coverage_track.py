#!/usr/bin/env python3
"""Figure 23: one deleted gene seen by three methods (FAR4, AT3G44540, absent in Cvi-0), run in env `map`.

First export read depth for the region from all 8 BAMs (env `map`, samtools; MAPQ >= 10 as in the PAV calls):
    samtools depth -a -Q 10 -r NC_003074.8:16110001-16150000 data/bam/{Col-0,An-1,Kyoto,Kondara,Cvi-0,Ler-1,Se-0,Ale-Stenar-44-4}/*.markdup.bam \
        > results/assembly/fig16_depth.tsv
Then:
    python scripts/plot_coverage_track.py

Inputs (read only): results/assembly/fig16_depth.tsv (chrom, pos, then one depth column per BAM in the order above),
                    data/reference/genes.bed, results/sv/pass_deletions.tsv, results/assembly/liftoff_vs_pav_genes.tsv.gz
Output: results/figures/23_validation_FAR4_example__coverage_track.png (300 dpi)
Depth is divided by each accession's median nuclear depth (the PAV normaliser) so the tracks are comparable.
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ORDER = ["Col-0", "An-1", "Kyoto", "Kondara", "Cvi-0", "Ler-1", "Se-0", "Ale-Stenar-44-4"]
MEDIAN = {"Col-0": 11, "An-1": 7, "Kyoto": 8, "Kondara": 9, "Cvi-0": 7, "Ler-1": 12, "Se-0": 9, "Ale-Stenar-44-4": 8}
CHROM, START, END = "NC_003074.8", 16110000, 16150000  # 0-based start, as in the BED
GENE, ACC = "AT3G44540", "Cvi-0"
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
    ap.add_argument("--depth", default="results/assembly/fig16_depth.tsv")
    ap.add_argument("--bed", default="data/reference/genes.bed")
    ap.add_argument("--dels", default="results/sv/pass_deletions.tsv")
    ap.add_argument("--liftoff", default="results/assembly/liftoff_vs_pav_genes.tsv.gz")
    ap.add_argument("--bin", type=int, default=200, help="bp per point (depth averaged)")
    ap.add_argument("--out", default="results/figures/23_validation_FAR4_example__coverage_track.png")
    a = ap.parse_args()

    dp = pd.read_csv(a.depth, sep="\t", header=None)
    assert dp.shape[1] == 2 + len(ORDER), f"expected {2 + len(ORDER)} columns, got {dp.shape[1]}"
    dp.columns = ["chrom", "pos"] + ORDER
    assert (dp["chrom"] == CHROM).all() and len(dp) == END - START, "depth file does not match the region"
    b = (dp["pos"] - 1 - START) // a.bin
    binned = dp.groupby(b)[ORDER].mean()
    x = (START + binned.index * a.bin + a.bin / 2) / 1e6

    bed = pd.read_csv(a.bed, sep="\t", header=None, usecols=[0, 1, 2, 3, 5], names=["chrom", "s", "e", "gene", "strand"])
    genes = bed[(bed["chrom"] == CHROM) & (bed["e"] > START) & (bed["s"] < END)]
    g = genes[genes["gene"] == GENE].iloc[0]
    dels = pd.read_csv(a.dels, sep="\t", header=None, names=["acc", "chrom", "pos", "end"])
    dl = dels[(dels["acc"] == ACC) & (dels["chrom"] == CHROM) & (dels["pos"] < g["e"]) & (dels["end"] > g["s"])
              & (dels["end"] - dels["pos"] <= 100_000)]
    lo = pd.read_csv(a.liftoff, sep="\t")
    lo_row = lo[(lo["gene"] == GENE) & (lo["accession"] == ACC)].iloc[0]
    inside = dp[(dp["pos"] > g["s"]) & (dp["pos"] <= g["e"])][ORDER]
    print(f"{GENE} {CHROM}:{g['s'] + 1}-{g['e']} ({g['e'] - g['s']} bp); mean depth inside the gene:")
    print("  " + ", ".join(f"{k} {v:.1f}" for k, v in inside.mean().items()))
    print(f"Delly deletion(s) in {ACC} over the gene: {[(int(p), int(e)) for p, e in zip(dl['pos'], dl['end'])]}")
    print(f"Liftoff in the {ACC} assembly: {lo_row['asm_class']}; read call: {lo_row['call']}")

    n = len(ORDER)
    fig, axes = plt.subplots(n + 1, 1, figsize=(10, 9), sharex=True,
                             gridspec_kw={"height_ratios": [1] * n + [0.9], "hspace": 0.15})
    for ax, acc in zip(axes, ORDER):
        y = binned[acc] / MEDIAN[acc]
        col = ORANGE if acc == ACC else BLUE
        ax.fill_between(x, 0, y, step="mid", color=col, alpha=0.8, lw=0)
        ax.axvspan(g["s"] / 1e6, g["e"] / 1e6, color=GRID, lw=0, zorder=0)
        ax.axhline(1, color=GREY, lw=0.6, ls=":")
        ax.set_ylim(0, 2.6)
        ax.set_yticks([0, 1, 2])
        ax.set_ylabel(acc, rotation=0, ha="right", va="center", fontsize=8.5,
                      color=ORANGE if acc == ACC else INK2, fontweight="bold" if acc == ACC else "normal")
        ax.tick_params(axis="x", length=0)
    axes[0].set_title("Read depth (relative to each accession's median; MAPQ >= 10)  -  grey band = FAR4 (AT3G44540)")

    ax = axes[-1]
    for _, r in genes.iterrows():
        hl = r["gene"] == GENE
        ax.add_patch(plt.Rectangle((max(r["s"], START) / 1e6, 0.55), (min(r["e"], END) - max(r["s"], START)) / 1e6,
                                   0.3, color=INK if hl else GREY, lw=0))
        mid = (max(r["s"], START) + min(r["e"], END)) / 2e6
        ax.text(mid, 1.0, f"{r['gene']}{' (FAR4)' if hl else ''} {r['strand']}", ha="center", va="bottom",
                fontsize=7.5, color=INK if hl else INK2, fontweight="bold" if hl else "normal")
    for _, r in dl.iterrows():
        ax.plot([r["pos"] / 1e6, r["end"] / 1e6], [0.2, 0.2], color=ORANGE, lw=4, solid_capstyle="butt")
        ax.text(r["end"] / 1e6 + 0.0004, 0.2, f"Delly deletion in {ACC} ({(r['end'] - r['pos']) / 1000:.1f} kb)",
                va="center", fontsize=7.5, color=ORANGE)
    ax.set_ylim(0, 1.6)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xlim(START / 1e6, END / 1e6)
    ax.set_xlabel(f"position on Col-0 chromosome 3 (Mb)")
    fig.suptitle(f"One gene, three methods: FAR4 is missing in {ACC} (reads: no coverage; Delly: deletion; "
                 f"assembly: Liftoff cannot place it)", x=0.01, ha="left", fontsize=10.5)
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    fig.savefig(a.out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
