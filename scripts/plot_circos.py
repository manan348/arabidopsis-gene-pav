#!/usr/bin/env python3
"""Genome overview as a circos plot: the main data layers of the project on one circle. Env `map`.

Run (after scripts/plot_diversity.py and scripts/plot_nlr_map.py):
    python scripts/plot_circos.py

Inputs (read only): data/reference/genes.bed, results/variants/pi_windows.tsv, results/pav/pav_long.tsv.gz,
    results/sv/pass_deletions.tsv, results/analysis/nlr_genes.tsv
Output: results/figures/02_overview_genome__circos_plot.png

Rings, outside to inside (500 kb windows, Col-0 coordinates): chromosomes (ticks every 5 Mb, centromere dark),
NLR genes (ticks; red = absent in >= 1 accession), gene density, SNPs per kb, nucleotide diversity pi,
absent-gene calls (summed over the 7 accessions), Delly PASS deletions (<= 100 kb, 7 accessions).
Descriptive only; each ring is scaled to its own maximum.
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

CHR = {"NC_003070.9": "1", "NC_003071.7": "2", "NC_003074.8": "3", "NC_003075.7": "4", "NC_003076.8": "5"}
CHR_LEN = {"1": 30427671, "2": 19698289, "3": 23459830, "4": 18585056, "5": 26975502}
CEN = {"1": 15086045, "2": 3607929, "3": 13799418, "4": 3956021, "5": 11725024}  # approximate (TAIR10)
INK, INK2, LIGHT, SURFACE = "#0b0b0b", "#52514e", "#c3c2b7", "#fcfcfb"
GAP = 0.03  # radians between chromosomes
RINGS = [  # (key, label, colour, inner radius, height)
    ("genes", "gene density", "#7fb3ea", 0.80, 0.10),
    ("snp", "SNPs per kb", "#2a78d6", 0.68, 0.10),
    ("pi", "diversity (pi)", "#1c4f8f", 0.56, 0.10),
    ("absent", "absent genes", "#eb6834", 0.44, 0.10),
    ("dels", "Delly deletions", "#8e6bbf", 0.32, 0.10),
]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bed", default="data/reference/genes.bed")
    ap.add_argument("--pi", default="results/variants/pi_windows.tsv")
    ap.add_argument("--pav", default="results/pav/pav_long.tsv.gz")
    ap.add_argument("--dels", default="results/sv/pass_deletions.tsv")
    ap.add_argument("--nlr", default="results/analysis/nlr_genes.tsv")
    ap.add_argument("--win", type=int, default=500_000)
    ap.add_argument("--out", default="results/figures/02_overview_genome__circos_plot.png")
    a = ap.parse_args()

    total = sum(CHR_LEN.values())
    scale = (2 * np.pi - GAP * len(CHR_LEN)) / total
    start, cur = {}, np.pi / 2  # clockwise from the top
    for ch, L in CHR_LEN.items():
        start[ch] = cur
        cur -= L * scale + GAP
    ang = lambda ch, pos: start[ch] - np.asarray(pos) * scale

    bed = pd.read_csv(a.bed, sep="\t", header=None, usecols=[0, 1, 2, 3], names=["chrom", "s", "e", "gene"])
    bed = bed[bed["chrom"].isin(CHR)].drop_duplicates("gene")
    bed["chr"] = bed["chrom"].map(CHR)
    bed["win"] = bed["s"] // a.win
    p = pd.read_csv(a.pav, sep="\t", usecols=["gene", "call"])
    nmiss = (p["call"] == "absent").groupby(p["gene"]).sum()
    bed["n_missing"] = bed["gene"].map(nmiss).fillna(0)
    pi = pd.read_csv(a.pi, sep="\t", dtype={"chr": str})
    pi["win"] = pi["start"] // a.win
    dl = pd.read_csv(a.dels, sep="\t", header=None, names=["acc", "chrom", "pos", "end"])
    dl = dl[(dl["acc"] != "Col-0") & dl["chrom"].isin(CHR) & (dl["end"] - dl["pos"] <= 100_000)]
    dl["chr"] = dl["chrom"].map(CHR)
    dl["win"] = dl["pos"] // a.win
    nlr = pd.read_csv(a.nlr, sep="\t", dtype={"chr": str})

    tracks = {}
    for ch, L in CHR_LEN.items():
        nw = (L - 1) // a.win + 1
        idx = range(nw)
        b = bed[bed["chr"] == ch]
        q = pi[pi["chr"] == ch].groupby("win")
        size = np.minimum(a.win, L - np.arange(nw) * a.win)
        tracks[ch] = {
            "genes": b.groupby("win").size().reindex(idx, fill_value=0).to_numpy() / size * 1e6,
            "snp": q["n_snp"].sum().reindex(idx, fill_value=0).to_numpy() / size * 1e3,
            "pi": q["pi"].mean().reindex(idx, fill_value=0).to_numpy(),
            "absent": b.groupby("win")["n_missing"].sum().reindex(idx, fill_value=0).to_numpy(),
            "dels": dl[dl["chr"] == ch].groupby("win").size().reindex(idx, fill_value=0).to_numpy(),
            "size": size,
        }
    vmax = {k: max(np.percentile(t[k], 99.5) for t in tracks.values()) for k, *_ in RINGS}

    fig = plt.figure(figsize=(10, 10.8))
    ax = fig.add_axes([0.04, 0.1, 0.92, 0.82], projection="polar")
    ax.set_axis_off()
    ax.set_ylim(0, 1.13)
    for ch, L in CHR_LEN.items():
        th0, th1 = ang(ch, 0), ang(ch, L)
        t = np.linspace(th0, th1, 200)
        ax.fill_between(t, 0.93, 0.97, color="#ecebe6", lw=0)
        ax.plot(t, np.full_like(t, 0.97), color=INK2, lw=0.8)
        ax.bar(ang(ch, CEN[ch]), 0.04, width=1.2e6 * scale, bottom=0.93, color=INK2, lw=0)
        for mb in range(0, int(L / 1e6) + 1, 5):
            th = ang(ch, mb * 1e6)
            ax.plot([th, th], [0.97, 0.985], color=INK2, lw=0.6)
            if mb and mb % 10 == 0:
                ax.text(th, 1.02, f"{mb}", ha="center", va="center", fontsize=6.5, color=INK2)
        mid = ang(ch, L / 2)
        ax.text(mid, 1.09, f"Chr{ch}", ha="center", va="center", fontsize=11, fontweight="bold")
        nl = nlr[nlr["chr"] == ch]
        for _, r in nl.iterrows():
            th = ang(ch, r["start"])
            var = r["n_missing"] > 0
            ax.plot([th, th], [0.905, 0.925], color="#c4501f" if var else INK2, lw=1.1 if var else 0.6)
        tr = tracks[ch]
        nw = len(tr["size"])
        th_c = ang(ch, np.arange(nw) * a.win + tr["size"] / 2)
        width = tr["size"] * scale
        for key, _, col, r0, hgt in RINGS:
            ax.bar(th_c, hgt, width=width, bottom=r0, color="#f3f2ee", lw=0)
            vals = np.clip(tr[key] / vmax[key], 0, 1) * hgt
            ax.bar(th_c, vals, width=width, bottom=r0, color=col, lw=0)
    for key, lab, col, r0, hgt in RINGS:
        ax.text(np.pi / 2 + GAP * 0.0, r0 + hgt / 2, "", fontsize=1)
    h = [plt.Rectangle((0, 0), 1, 1, color=c) for _, _, c, _, _ in RINGS]
    labels = [f"{lab} (ring {i + 1} from outside)" for i, (_, lab, _, _, _) in enumerate(RINGS)]
    h = [plt.Line2D([], [], color="#c4501f", lw=2), plt.Line2D([], [], color=INK2, lw=1)] + h
    labels = ["NLR gene absent in >= 1 accession", "NLR gene present in all"] + labels
    fig.legend(h, labels, loc="lower center", ncol=3, frameon=False, fontsize=8.5, bbox_to_anchor=(0.5, 0.01))
    fig.suptitle("Genome overview: the Arabidopsis genome with all main data layers", x=0.03, y=0.975, ha="left",
                 fontsize=13, fontweight="bold")
    fig.text(0.03, 0.95, "500 kb windows on the Col-0 genome; each ring scaled to its own maximum. Near the centromeres "
             "(dark marks) gene density drops,\nwhile diversity and absent genes rise. NLR genes (outer ticks) sit "
             "mostly on the chromosome arms, often in clusters.", fontsize=8.5, color=INK2, va="top")
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    fig.savefig(a.out, dpi=300, facecolor=SURFACE)
    plt.close(fig)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
