#!/usr/bin/env python3
"""Nucleotide diversity (pi) along the chromosomes, with NLR genes and absent-gene density. Env `map`.

Run:
    python scripts/plot_diversity.py

Inputs (read only): results/variants/filtered_gt.tsv.gz (chrom, pos, type, qual, depth, 8 haploid genotypes),
    data/reference/genes.bed, results/pav/pav_long.tsv.gz, results/pav/nlr_genes.tsv
Outputs: results/variants/pi_windows.tsv (100 kb windows), results/variants/pi_genes.tsv (pi per gene),
    results/figures/15_variants_diversity__sliding_window_plot_violin_plot.png (renumbered later into the sequence)

pi per SNP site = 2k(n-k) / (n(n-1)), k = accessions with the alternative allele, n = accessions with a call
(inbred lines, so one allele each). Window pi = sum over SNP sites / window length, i.e. per bp of reference;
sites that are not callable are not removed, so pi is a lower bound in poorly mapped regions (centromeres).
Gene pi = sum over SNP sites inside the gene / gene length. NLR vs other genes: Mann-Whitney U (two-sided),
on genes present in all 8 accessions so that deleted genes do not lower the values.
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

ACC = ["Col-0", "An-1", "Kyoto", "Kondara", "Cvi-0", "Ler-1", "Se-0", "Ale-Stenar-44-4"]
CHR = {"NC_003070.9": "1", "NC_003071.7": "2", "NC_003074.8": "3", "NC_003075.7": "4", "NC_003076.8": "5"}
CHR_LEN = {"1": 30427671, "2": 19698289, "3": 23459830, "4": 18585056, "5": 26975502}
CEN = {"1": 15086045, "2": 3607929, "3": 13799418, "4": 3956021, "5": 11725024}  # approximate (TAIR10)
BLUE, ORANGE, LIGHT, INK, INK2, GRID, SURFACE = "#2a78d6", "#eb6834", "#c3c2b7", "#0b0b0b", "#52514e", "#e1e0d9", "#fcfcfb"
plt.rcParams.update({
    "font.size": 9, "axes.edgecolor": LIGHT, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "axes.titlesize": 10, "axes.titleweight": "bold", "axes.titlelocation": "left",
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
})


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gt", default="results/variants/filtered_gt.tsv.gz")
    ap.add_argument("--bed", default="data/reference/genes.bed")
    ap.add_argument("--pav", default="results/pav/pav_long.tsv.gz")
    ap.add_argument("--nlr", default="results/pav/nlr_genes.tsv")
    ap.add_argument("--win", type=int, default=100_000)
    ap.add_argument("--windows", default="results/variants/pi_windows.tsv")
    ap.add_argument("--genes", default="results/variants/pi_genes.tsv")
    ap.add_argument("--out", default="results/figures/15_variants_diversity__sliding_window_plot_violin_plot.png")
    a = ap.parse_args()

    g = pd.read_csv(a.gt, sep="\t", header=None, names=["chrom", "pos", "type", "qual", "dp"] + ACC,
                    na_values=".", dtype={k: "float32" for k in ACC})
    g = g[(g["type"] == "SNP") & g["chrom"].isin(CHR)]
    G = g[ACC].to_numpy()
    n = np.sum(~np.isnan(G), axis=1)
    k = np.nansum(G, axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        pi_site = np.where(n >= 2, 2 * k * (n - k) / (n * (n - 1)), 0.0)
    s = pd.DataFrame({"chr": g["chrom"].map(CHR).to_numpy(), "pos": g["pos"].to_numpy(), "pi": pi_site})
    s = s[s["pi"] > 0]
    print(f"{len(g):,} SNP sites, {len(s):,} segregating with >= 2 calls")

    s["win"] = (s["pos"] - 1) // a.win
    rows = []
    for ch, L in CHR_LEN.items():
        nw = (L - 1) // a.win + 1
        sub = s[s["chr"] == ch].groupby("win")["pi"].agg(["sum", "size"]).reindex(range(nw), fill_value=0)
        size = np.minimum(a.win, L - np.arange(nw) * a.win)
        rows.append(pd.DataFrame({"chr": ch, "start": np.arange(nw) * a.win, "end": np.arange(nw) * a.win + size,
                                  "n_snp": sub["size"].to_numpy(), "pi": sub["sum"].to_numpy() / size}))
    w = pd.concat(rows, ignore_index=True)
    w.to_csv(a.windows, sep="\t", index=False, float_format="%.6g")
    print(f"genome-wide pi (per bp, lower bound): {s['pi'].sum() / sum(CHR_LEN.values()):.5f}")

    bed = pd.read_csv(a.bed, sep="\t", header=None, usecols=[0, 1, 2, 3], names=["chrom", "s", "e", "gene"])
    bed = bed[bed["chrom"].isin(CHR)].drop_duplicates("gene")
    bed["chr"] = bed["chrom"].map(CHR)
    pis = []
    for ch, b in bed.groupby("chr"):
        ss = s[s["chr"] == ch].sort_values("pos")
        pos, cum = ss["pos"].to_numpy(), np.concatenate([[0], np.cumsum(ss["pi"].to_numpy())])
        lo = np.searchsorted(pos, b["s"].to_numpy() + 1)
        hi = np.searchsorted(pos, b["e"].to_numpy(), side="right")
        pis.append(pd.DataFrame({"gene": b["gene"].to_numpy(), "chr": ch, "start": b["s"].to_numpy() + 1,
                                 "pi": (cum[hi] - cum[lo]) / (b["e"] - b["s"]).to_numpy()}))
    gp = pd.concat(pis, ignore_index=True)
    p = pd.read_csv(a.pav, sep="\t", usecols=["gene", "call"])
    nmiss = (p["call"] == "absent").groupby(p["gene"]).sum()
    gp["n_missing"] = gp["gene"].map(nmiss).fillna(0).astype(int)
    nlr = set(pd.read_csv(a.nlr, sep="\t")["gene"])
    gp["nlr"] = gp["gene"].isin(nlr)
    gp.to_csv(a.genes, sep="\t", index=False, float_format="%.6g")
    core = gp[gp["n_missing"] == 0]
    a_n, a_o = core.loc[core["nlr"], "pi"], core.loc[~core["nlr"], "pi"]
    u = mannwhitneyu(a_n, a_o, alternative="two-sided")
    print(f"gene pi, genes present in all 8: NLR median {a_n.median():.5f} (n={len(a_n)}) vs other "
          f"{a_o.median():.5f} (n={len(a_o)}); ratio {a_n.median() / a_o.median():.1f}x; Mann-Whitney p = {u.pvalue:.1e}")

    # absent genes per window (sum over accessions), for the second track
    bed["win"] = bed["s"] // a.win
    bed["n_missing"] = bed["gene"].map(nmiss).fillna(0)
    absw = bed.groupby(["chr", "win"])["n_missing"].sum()

    fig = plt.figure(figsize=(13, 10))
    gs = fig.add_gridspec(6, 2, width_ratios=[3.2, 1], hspace=0.55, wspace=0.18)
    ymax = np.percentile(w["pi"], 99.7) * 1.1
    for i, (ch, L) in enumerate(CHR_LEN.items()):
        ax = fig.add_subplot(gs[i, 0])
        ww = w[w["chr"] == ch]
        x = (ww["start"] + ww["end"]) / 2e6
        ax.fill_between(x, 0, ww["pi"].clip(upper=ymax), color=BLUE, alpha=0.75, lw=0, step="mid")
        ax.axvspan((CEN[ch] - 1e6) / 1e6, (CEN[ch] + 1e6) / 1e6, color=GRID, lw=0, zorder=0)
        ab = absw.loc[ch] if ch in absw.index.get_level_values(0) else pd.Series(dtype=float)
        ax2 = ax.twinx()
        ax2.bar((ab.index * a.win + a.win / 2) / 1e6, ab.values, width=a.win / 1e6 * 1.6, color=ORANGE, alpha=0.9,
                lw=0)
        ax2.set_ylim(0, max(absw.max() * 1.15, 1))
        ax2.set_zorder(ax.get_zorder() + 1)
        ax2.set_yticks([])
        ax2.spines["right"].set_visible(False)
        nl = gp[(gp["chr"] == ch) & gp["nlr"]]
        ax.scatter(nl["start"] / 1e6, np.full(len(nl), ymax * 0.97), marker="v", s=14, color=INK, zorder=4)
        ax.set_xlim(0, 30.6)
        ax.set_ylim(0, ymax)
        ax.set_ylabel(f"Chr{ch}\npi", rotation=0, ha="right", va="center", fontsize=8.5)
        ax.tick_params(axis="y", labelsize=7)
        if i == 0:
            ax.set_title("(a) Diversity along the 5 chromosomes")
    ax.set_xlabel("position on the Col-0 chromosome (Mb)")

    axb = fig.add_subplot(gs[:, 1])
    groups = [("other\ngenes", a_o, LIGHT), ("NLR\ngenes", a_n, BLUE),
              ("variable\ngenes", gp.loc[gp["n_missing"] > 0, "pi"], ORANGE)]
    data = [np.log10(v[v > 0]) for _, v, _ in groups]
    vp = axb.violinplot(data, showmedians=True, widths=0.8)
    for body, (_, _, c) in zip(vp["bodies"], groups):
        body.set_facecolor(c)
        body.set_alpha(0.7)
    for part in ("cmedians", "cmins", "cmaxes", "cbars"):
        vp[part].set_color(INK2)
    axb.set_xticks([1, 2, 3])
    axb.set_xticklabels([f"{l}\nn={len(v):,}" for l, v, _ in groups], fontsize=8)
    axb.set_ylabel("gene pi (log10; genes with >= 1 SNP)")
    axb.set_title(f"(b) NLR genes: {a_n.median() / a_o.median():.1f}x more diverse\n(Mann-Whitney p = {u.pvalue:.0e})")
    axb.grid(axis="y", color=GRID, lw=0.5)
    fig.suptitle("Genetic diversity: highest around the centromeres, and twice as high in NLR resistance genes",
                 x=0.01, ha="left", fontsize=12, fontweight="bold", y=0.965)
    fig.text(0.01, 0.94, "(a) blue = nucleotide diversity pi in 100 kb windows; orange = absent-gene calls per window "
             "(summed over accessions); triangles = NLR genes; grey band = centromere +-1 Mb.", fontsize=8.5,
             color=INK2)
    fig.text(0.01, 0.925, "pi from 1.6 M filtered SNPs of 8 inbred accessions, per bp of reference (a lower bound where "
             "reads map poorly). Gene pi includes introns; (b) NLR and other genes present in all 8 accessions.",
             fontsize=8.5, color=INK2)
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    fig.savefig(a.out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("wrote", a.out, a.windows, a.genes)


if __name__ == "__main__":
    main()
