#!/usr/bin/env python3
"""Figure 21: are short genes called absent more often, and are those calls less reliable? Env `map`.

Run:
    python scripts/plot_length_bias.py

Inputs (read only): results/pav/pav_long.tsv.gz, data/reference/genes.bed, results/assembly/liftoff_vs_pav_genes.tsv.gz
Output: results/figures/21_validation_length_bias__line_plot_confidence_intervals.png (300 dpi)

(a) Share of genes the reads call absent, by gene length, for each of the 7 accessions.
(b) For the 4 accessions with an assembly: share of read-absent genes that the assembly confirms (Liftoff does not
    find the gene intact), by gene length. If short-gene calls were mostly method errors, (b) would drop for short
    genes. Confidence intervals: Wilson 95% for a proportion. Length = gene span in Col-0 (incl. introns).
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ORDER = ["An-1", "Kyoto", "Kondara", "Cvi-0", "Ler-1", "Se-0", "Ale-Stenar-44-4"]
ASM = ["An-1", "Kyoto", "Cvi-0", "Ler-1"]
COL = {"An-1": "#2a78d6", "Ler-1": "#1c4f8f", "Ale-Stenar-44-4": "#7fb3ea", "Se-0": "#e3b23c", "Cvi-0": "#eb6834",
       "Kondara": "#8e6bbf", "Kyoto": "#2a9d8f"}
EDGES = [0, 500, 1000, 2000, 3000, 5000, 10000, 1e9]
LABELS = ["<0.5", "0.5-1", "1-2", "2-3", "3-5", "5-10", ">10"]
LIGHT, INK, INK2, GRID, SURFACE = "#c3c2b7", "#0b0b0b", "#52514e", "#e1e0d9", "#fcfcfb"
plt.rcParams.update({
    "font.size": 9, "axes.edgecolor": LIGHT, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "axes.titlesize": 10, "axes.titleweight": "bold", "axes.titlelocation": "left",
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
})


def wilson(k, n, z=1.96):
    k, n = np.asarray(k, float), np.asarray(n, float)
    with np.errstate(invalid="ignore", divide="ignore"):
        p = k / n
        den = 1 + z ** 2 / n
        c = (p + z ** 2 / (2 * n)) / den
        h = z * np.sqrt(p * (1 - p) / n + z ** 2 / (4 * n ** 2)) / den
    return p, c - h, c + h


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pav", default="results/pav/pav_long.tsv.gz")
    ap.add_argument("--bed", default="data/reference/genes.bed")
    ap.add_argument("--liftoff", default="results/assembly/liftoff_vs_pav_genes.tsv.gz")
    ap.add_argument("--out", default="results/figures/21_validation_length_bias__line_plot_confidence_intervals.png")
    a = ap.parse_args()

    bed = pd.read_csv(a.bed, sep="\t", header=None, usecols=[1, 2, 3], names=["s", "e", "gene"])
    length = (bed["e"] - bed["s"]).groupby(bed["gene"]).max()
    bins = pd.Series(pd.cut(length, EDGES, labels=LABELS, right=False), index=length.index)
    n_bin = bins.value_counts().reindex(LABELS)
    print("genes per length bin (kb):", dict(n_bin))

    p = pd.read_csv(a.pav, sep="\t", usecols=["gene", "accession", "call"])
    p["bin"] = p["gene"].map(bins)
    lo = pd.read_csv(a.liftoff, sep="\t", usecols=["gene", "accession", "call", "asm_class"])
    ra = lo[lo["call"] == "absent"].copy()
    ra["bin"] = ra["gene"].map(bins)
    ra["confirmed"] = ra["asm_class"] != "found"

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 4.9), gridspec_kw={"wspace": 0.25})
    x = np.arange(len(LABELS))
    for acc in ORDER:
        d = p[p["accession"] == acc]
        k = (d["call"] == "absent").groupby(d["bin"], observed=False).sum().reindex(LABELS)
        n = d.groupby("bin", observed=False).size().reindex(LABELS)
        ax1.plot(x, 100 * k / n, marker="o", ms=4, color=COL[acc], lw=1.5, label=acc)
    ax1.set_xticks(x)
    ax1.set_xticklabels([f"{l}\n(n={n_bin[l]:,})" for l in LABELS], fontsize=7.5)
    ax1.set_xlabel("gene length (kb)")
    ax1.set_ylabel("% of genes called absent by reads")
    ax1.set_title("(a) Short genes are called absent more often")
    ax1.legend(frameon=False, fontsize=7.5, ncol=2)
    ax1.grid(axis="y", color=GRID, lw=0.6)
    ax1.set_axisbelow(True)

    off = np.linspace(-0.18, 0.18, len(ASM))
    print("\nshare of read-absent genes confirmed by the assembly, by length bin:")
    for o, acc in zip(off, ASM):
        d = ra[ra["accession"] == acc]
        k = d["confirmed"].groupby(d["bin"], observed=False).sum().reindex(LABELS)
        n = d.groupby("bin", observed=False).size().reindex(LABELS)
        pr, l, h = wilson(k, n)
        ok = n.to_numpy() >= 5
        ax2.errorbar(x[ok] + o, 100 * pr[ok], yerr=[100 * (pr - l)[ok], 100 * (h - pr)[ok]], fmt="o", ms=4.5,
                     color=COL[acc], elinewidth=1, capsize=2, label=acc)
        print(f"  {acc:<6} " + "  ".join(f"{lab}: {int(kk)}/{int(nn)}" for lab, kk, nn in zip(LABELS, k, n)))
    pooled = ra.groupby("bin", observed=False)["confirmed"].agg(["sum", "size"]).reindex(LABELS)
    pr, _, _ = wilson(pooled["sum"], pooled["size"])
    okp = pooled["size"].to_numpy() >= 5
    ax2.plot(x[okp], 100 * pr[okp], color=INK2, lw=1, ls="--", label="4 pooled")
    ax2.set_xticks(x)
    ax2.set_xticklabels([f"{l}\n(n={int(m)})" for l, m in zip(LABELS, pooled["size"])], fontsize=7.5)
    ax2.set_xlabel("gene length (kb); n = read-absent calls in the 4 accessions")
    ax2.set_ylabel("% confirmed by the assembly (Wilson 95% CI)")
    ax2.set_ylim(0, 105)
    ax2.set_title("(b) ...but are confirmed by the assemblies about as often")
    ax2.legend(frameon=False, fontsize=7.5, loc="lower right", ncol=3)
    ax2.grid(axis="y", color=GRID, lw=0.6)
    ax2.set_axisbelow(True)
    fig.suptitle("Gene length and absent calls: real biology or a method bias?", x=0.01, ha="left", fontsize=12,
                 fontweight="bold")
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    fig.savefig(a.out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
