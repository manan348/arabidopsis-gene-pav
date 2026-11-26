#!/usr/bin/env python3
"""Figure 37: NLR disease-resistance genes on the chromosomes, which ones are absent and in how many accessions. Env `map`.

Run:
    python scripts/plot_nlr_map.py

Inputs (read only): data/reference/TAIR10.1_Col-0/GCF_000001735.4/genomic.gff (RefSeq product names),
    data/reference/genes.bed, results/pav/pav_long.tsv.gz
Outputs: results/figures/37_function_nlr_genes__chromosome_ideogram_fisher_test.png (300 dpi), results/pav/nlr_genes.tsv

NLR = protein-coding gene whose RefSeq CDS product name contains NB-ARC or NBS (TIR-NBS-LRR, CC-NBS-LRR,
TIR-NBS, NB-ARC ...). A name-based list, so it is an approximation of the full NLR set.
Test: Fisher's exact test, NLR vs other genes, absent in >= 1 accession or not (two-sided).
Clusters: NLRs within 200 kb of another NLR count as clustered (a common working definition).
"""
import argparse
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap
from scipy.stats import fisher_exact

CHR = {"NC_003070.9": "1", "NC_003071.7": "2", "NC_003074.8": "3", "NC_003075.7": "4", "NC_003076.8": "5"}
CHR_LEN = {"1": 30427671, "2": 19698289, "3": 23459830, "4": 18585056, "5": 26975502}
CEN = {"1": 15086045, "2": 3607929, "3": 13799418, "4": 3956021, "5": 11725024}  # approximate (TAIR10)
PAT = re.compile(r"NB-ARC|NBS", re.I)
LIGHT, INK, INK2, GRID, SURFACE = "#c3c2b7", "#0b0b0b", "#52514e", "#e1e0d9", "#fcfcfb"
plt.rcParams.update({
    "font.size": 9, "axes.edgecolor": LIGHT, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "axes.titlesize": 10, "axes.titleweight": "bold", "axes.titlelocation": "left",
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
})


def nlr_genes(gff):
    out = {}
    with open(gff) as f:
        for line in f:
            if "\tCDS\t" not in line:
                continue
            m = re.search(r"TAIR:(AT[1-5]G\d{5})", line)
            pr = re.search(r";product=([^;\n]*)", line)
            if m and pr and PAT.search(pr.group(1)):
                out.setdefault(m.group(1), pr.group(1).replace("%2C", ","))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gff", default="data/reference/TAIR10.1_Col-0/GCF_000001735.4/genomic.gff")
    ap.add_argument("--bed", default="data/reference/genes.bed")
    ap.add_argument("--pav", default="results/pav/pav_long.tsv.gz")
    ap.add_argument("--table", default="results/pav/nlr_genes.tsv")
    ap.add_argument("--out", default="results/figures/37_function_nlr_genes__chromosome_ideogram_fisher_test.png")
    a = ap.parse_args()

    nlr = nlr_genes(a.gff)
    bed = pd.read_csv(a.bed, sep="\t", header=None, usecols=[0, 1, 2, 3], names=["chrom", "s", "e", "gene"])
    bed = bed.drop_duplicates("gene").set_index("gene")
    p = pd.read_csv(a.pav, sep="\t", usecols=["gene", "call"])
    nmiss = (p["call"] == "absent").groupby(p["gene"]).sum()
    nlr = {g: v for g, v in nlr.items() if g in nmiss.index}
    t = pd.DataFrame({"gene": list(nlr), "product": list(nlr.values())})
    t["chr"] = t["gene"].map(bed["chrom"]).map(CHR)
    t["start"] = t["gene"].map(bed["s"]) + 1
    t["n_missing"] = t["gene"].map(nmiss).astype(int)
    t = t.sort_values(["chr", "start"]).reset_index(drop=True)
    t["clustered"] = False
    for _, g in t.groupby("chr"):
        s = g["start"].to_numpy()
        gap = np.diff(s)
        near = np.zeros(len(s), bool)
        near[:-1] |= gap <= 200_000
        near[1:] |= gap <= 200_000
        t.loc[g.index, "clustered"] = near
    t.to_csv(a.table, sep="\t", index=False)

    var = nmiss > 0
    is_nlr = nmiss.index.isin(t["gene"])
    tab = [[int((is_nlr & var).sum()), int((is_nlr & ~var).sum())],
           [int((~is_nlr & var).sum()), int((~is_nlr & ~var).sum())]]
    odds, pval = fisher_exact(tab)
    nv = tab[0][0]
    print(f"NLR genes: {len(t)} ({100 * t['clustered'].mean():.0f}% in clusters); absent in >= 1 accession: {nv} "
          f"({100 * nv / len(t):.1f}%) vs other genes {100 * tab[1][0] / sum(tab[1]):.1f}%; "
          f"Fisher odds ratio {odds:.1f}, p = {pval:.1e}")
    vt = t[t["n_missing"] > 0]
    print(f"variable NLRs in clusters: {int(vt['clustered'].sum())} of {len(vt)}")

    cmap = ListedColormap(plt.cm.YlOrRd(np.linspace(0.15, 0.95, 7)))
    fig, ax = plt.subplots(figsize=(12, 5.2))
    for i, ch in enumerate(["1", "2", "3", "4", "5"]):
        y = 4 - i
        L = CHR_LEN[ch] / 1e6
        ax.add_patch(plt.Rectangle((0, y - 0.09), L, 0.18, fc="#ecebe6", ec=LIGHT, lw=0.8))
        ax.plot([CEN[ch] / 1e6] * 2, [y - 0.16, y + 0.16], color=INK2, lw=2.5, solid_capstyle="butt")
        g = t[t["chr"] == ch]
        c = g[g["n_missing"] == 0]
        ax.vlines(c["start"] / 1e6, y - 0.09, y + 0.09, color="#2a78d6", lw=0.9)
        v = g[g["n_missing"] > 0]
        ax.scatter(v["start"] / 1e6, np.full(len(v), y + 0.27), s=14 + 9 * v["n_missing"],
                   c=cmap((v["n_missing"] - 1) / 6), edgecolor=INK2, lw=0.4, zorder=3)
        ax.text(-0.5, y, f"Chr{ch}", ha="right", va="center", fontsize=9.5)
    ax.set_xlim(-2.6, 31)
    ax.set_ylim(-0.6, 4.7)
    ax.set_yticks([])
    for s in ["top", "right", "left"]:
        ax.spines[s].set_visible(False)
    ax.set_xlabel("position on the Col-0 chromosome (Mb); black bar = approximate centromere")
    h1 = plt.Line2D([], [], color="#2a78d6", lw=1.5, label="NLR present in all 8")
    h2 = plt.Line2D([], [], marker="o", ls="", color=cmap(0.5), markeredgecolor=INK2,
                    label="NLR absent in some accessions (dot above)")
    fig.legend(handles=[h1, h2], frameon=False, loc="upper left", bbox_to_anchor=(0.12, -0.0), ncol=2,
               fontsize=8.5)
    cax = fig.add_axes([0.72, -0.05, 0.18, 0.025])
    cb = plt.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0.5, 7.5)), cax=cax,
                      orientation="horizontal", ticks=range(1, 8))
    cb.set_label("accessions missing the gene", fontsize=8)
    cb.ax.tick_params(labelsize=7.5)
    cb.outline.set_visible(False)
    fig.suptitle(f"NLR disease-resistance genes: {nv} of {len(t)} are absent in at least one accession "
                 f"({100 * nv / len(t):.0f}% vs {100 * tab[1][0] / sum(tab[1]):.1f}% of other genes; "
                 f"Fisher p = {pval:.0e})", x=0.01, ha="left", fontsize=11.5, fontweight="bold")
    fig.text(0.01, 0.905, "NLRs from RefSeq product names (NB-ARC / NBS), an approximation of the full set. "
             "Ticks on the bar = NLRs present in all accessions; dots = NLRs missing in some.",
             fontsize=8.5, color=INK2)
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    fig.savefig(a.out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("wrote", a.out, "and", a.table)


if __name__ == "__main__":
    main()
