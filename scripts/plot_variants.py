#!/usr/bin/env python3
"""Figures from the filtered SNP/indel calls (Phase 5), run in env `map` after exporting the genotypes.

Step 1 (env `variants`): export the genotypes of the filtered VCF to a table
    bcftools query -f '%CHROM\\t%POS\\t%TYPE\\t%QUAL\\t%INFO/DP[\\t%GT]\\n' results/variants/filtered.vcf.gz | gzip > results/variants/filtered_gt.tsv.gz
Step 2 (env `map`):
    python scripts/plot_variants.py

Inputs (read only): results/variants/filtered_gt.tsv.gz, results/pav/pav_long.tsv.gz,
                    data/reference/genes.bed, the reference .fai
Outputs (PNG, 300 dpi) in results/figures/:
    13_variants__bar_chart_pca_distance_heatmap.png            variants per accession, SNP density, PCA, pairwise differences
    14_variants_snp_density__ecdf_dumbbell_plot.png     SNP density around absent vs present genes

Genotypes come from haploid calling (--ploidy 1): 0 = reference allele, 1 or more = alternative allele,
missing = no call. Only SNPs are used for density, PCA and distances. Everything is descriptive.
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.gridspec import GridSpec

ORDER = ["Col-0", "An-1", "Kyoto", "Kondara", "Cvi-0", "Ler-1", "Se-0", "Ale-Stenar-44-4"]
OTHERS = ORDER[1:]
CHROMS = {"NC_003070.9": "Chr1", "NC_003071.7": "Chr2", "NC_003074.8": "Chr3",
          "NC_003075.7": "Chr4", "NC_003076.8": "Chr5"}

BLUE, ORANGE, GREY, LIGHT = "#2a78d6", "#eb6834", "#898781", "#c3c2b7"
INK, INK2, GRID, SURFACE, CTRL = "#0b0b0b", "#52514e", "#e1e0d9", "#fcfcfb", "#52514e"
SEQ = LinearSegmentedColormap.from_list("blue_seq", ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"])

plt.rcParams.update({
    "font.size": 9, "axes.edgecolor": LIGHT, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "axes.titlesize": 10, "axes.titleweight": "bold", "axes.titlelocation": "left",
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": False,
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
})


# ----------------------------------------------------------------------------- data
def load_genotypes(path):
    names = ["chrom", "pos", "vtype", "qual", "dp"] + ORDER
    dtypes = {g: "float32" for g in ORDER}
    df = pd.read_csv(path, sep="\t", header=None, names=names, na_values=["."], dtype=dtypes)
    assert set(df["chrom"]) <= set(CHROMS), f"unexpected sequences: {set(df['chrom']) - set(CHROMS)}"
    assert df["pos"].is_monotonic_increasing or df["chrom"].nunique() > 1
    return df


def load_chrom_lengths(fai_path):
    fai = pd.read_csv(fai_path, sep="\t", header=None, usecols=[0, 1], names=["chrom", "length"]).set_index("chrom")
    assert set(CHROMS) <= set(fai.index), "reference index lacks the nuclear chromosomes"
    return fai.loc[list(CHROMS), "length"]


def save(fig, path):
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


# ----------------------------------------------------------------------------- figure 4
def fig4(df, chrom_len, out, bin_bp=100_000):
    alt = df[ORDER].to_numpy() >= 1  # NaN (no call) compares False
    is_snp = (df["vtype"] == "SNP").to_numpy()
    is_indel = df["vtype"].str.contains("INDEL").to_numpy()
    print(f"sites: {len(df)}  SNP: {is_snp.sum()}  with an INDEL: {is_indel.sum()}  other types: {(~is_snp & ~is_indel).sum()}")

    fig = plt.figure(figsize=(15, 9.5))
    gs = GridSpec(2, 3, height_ratios=[1, 0.8], hspace=0.38, wspace=0.28)

    # (a) non-reference calls per accession
    ax = fig.add_subplot(gs[0, 0])
    snp_n = alt[is_snp].sum(axis=0)
    ind_n = alt[is_indel].sum(axis=0)
    y = np.arange(len(ORDER))
    ax.barh(y - 0.19, snp_n / 1e3, height=0.36, color=BLUE, label="SNPs")
    ax.barh(y + 0.19, ind_n / 1e3, height=0.36, color=ORANGE, label="indels")
    for yi, v in zip(y, snp_n):
        ax.text(v / 1e3 + 6, yi - 0.19, f"{int(v):,}", va="center", fontsize=7.5, color=INK2)
    ax.set_yticks(y)
    ax.set_yticklabels(ORDER)
    ax.invert_yaxis()
    ax.set_xlim(0, snp_n.max() / 1e3 * 1.25)
    ax.set_xlabel("sites differing from Col-0 (thousands)")
    ax.set_title("(a) Differences from the reference")
    ax.legend(frameon=False, loc="upper right")  # top row (Col-0) is nearly empty, so nothing is covered

    # (c) PCA on sites genotyped in all 8 accessions
    complete = ~np.isnan(df[ORDER].to_numpy()).any(axis=1) & is_snp
    X = alt[complete].astype(np.float64)  # sites x accessions, 1 = alternative allele
    M = X.T - X.T.mean(axis=0, keepdims=True)
    U, S, _ = np.linalg.svd(M, full_matrices=False)
    coords, var = U * S, S ** 2 / np.sum(S ** 2)
    ax = fig.add_subplot(gs[0, 1])
    colours = [CTRL if a == "Col-0" else BLUE for a in ORDER]
    ax.scatter(coords[:, 0], coords[:, 1], s=60, color=colours, edgecolor=SURFACE, linewidth=1.2, zorder=3)
    for a, (px, py) in zip(ORDER, coords[:, :2]):
        ax.annotate(a, (px, py), xytext=(5, 4), textcoords="offset points", fontsize=8, color=INK2)
    ax.set_xlabel(f"PC1 ({100 * var[0]:.0f}% of the variation)")
    ax.set_ylabel(f"PC2 ({100 * var[1]:.0f}%)")
    ax.set_title(f"(b) PCA of the accessions ({int(complete.sum()):,} SNPs)")
    ax.grid(color=GRID, lw=0.6)
    ax.set_axisbelow(True)

    # (d) pairwise differences, per Mb of genome
    genome_mb = chrom_len.sum() / 1e6
    n = len(ORDER)
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            D[i, j] = np.abs(X[:, i] - X[:, j]).sum() / genome_mb
    ax = fig.add_subplot(gs[0, 2])
    im = ax.imshow(D, cmap=SEQ)
    for i in range(n):
        for j in range(n):
            ax.text(j, i, f"{D[i, j]:.0f}", ha="center", va="center", fontsize=7.5,
                    color="white" if D[i, j] > 0.6 * D.max() else INK)
    ax.set_xticks(range(n))
    ax.set_xticklabels(ORDER, rotation=45, ha="right")
    ax.set_yticks(range(n))
    ax.set_yticklabels(ORDER)
    ax.set_title("(c) SNP differences between accessions")
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("differing SNP sites per Mb\n(sites called in all 8 accessions)")

    # (e) SNP density along the chromosomes
    sub = gs[1, :].subgridspec(1, len(CHROMS) + 1, width_ratios=[chrom_len[c] for c in CHROMS] + [3e6], wspace=0.04)
    mats = {}
    for c in CHROMS:
        edges = np.arange(0, chrom_len[c] + bin_bp, bin_bp)
        rows = []
        in_c = (df["chrom"] == c).to_numpy() & is_snp
        pos = df["pos"].to_numpy()[in_c]
        alt_c = alt[in_c]
        for k in range(len(ORDER)):
            rows.append(np.histogram(pos[alt_c[:, k]], bins=edges)[0])
        mats[c] = np.array(rows)
    vmax = np.percentile(np.concatenate([m.ravel() for m in mats.values()]), 99.5)
    first = None
    for k, c in enumerate(CHROMS):
        ax = fig.add_subplot(sub[0, k], sharey=first)
        first = first or ax
        im = ax.imshow(mats[c], cmap=SEQ, vmin=0, vmax=vmax, aspect="auto", interpolation="nearest",
                       extent=(0, chrom_len[c] / 1e6, len(ORDER), 0))
        ax.set_title(CHROMS[c], fontsize=9)
        if k == 0:  # panel heading on its own line above the chromosome names, so it cannot run into "Chr2"
            ax.text(0, 1.16, f"(d) SNPs per {bin_bp // 1000} kb along the genome", transform=ax.transAxes,
                    fontsize=10, fontweight="bold", ha="left", va="bottom")
        ax.set_xlabel("Mb")
        if k == 0:
            ax.set_yticks(np.arange(len(ORDER)) + 0.5)
            ax.set_yticklabels(ORDER)
        else:
            ax.set_xticks(np.arange(5, chrom_len[c] / 1e6, 5))  # no "0" tick, so it cannot touch the previous chromosome's last label
            plt.setp(ax.get_yticklabels(), visible=False)
            ax.tick_params(axis="y", length=0)
    cax = fig.add_subplot(sub[0, -1])
    cax.axis("off")
    cb = fig.colorbar(im, ax=cax, fraction=0.4, aspect=10)
    cb.set_label(f"SNPs per {bin_bp // 1000} kb (capped at the 99.5th percentile)")
    fig.suptitle("Filtered variant calls of the 8 accessions against the Col-0 reference", x=0.01, ha="left", fontsize=11)
    save(fig, out)


# ----------------------------------------------------------------------------- figure 5
def flank_density(pos_sorted, start, end, chrom_length, w):
    """SNPs per kb in the w bp on each side of a gene (gene itself excluded), clipped at the chromosome ends."""
    lo = np.maximum(start - w, 0)
    hi = np.minimum(end + w, chrom_length)
    left = np.searchsorted(pos_sorted, start, side="left") - np.searchsorted(pos_sorted, lo, side="left")
    right = np.searchsorted(pos_sorted, hi, side="right") - np.searchsorted(pos_sorted, end, side="right")
    length = (start - lo) + (hi - end)
    assert (length > 0).all()
    return (left + right) / length * 1000.0


def fig5(df, pav_path, bed_path, chrom_len, out, flank_bp):
    bed = pd.read_csv(bed_path, sep="\t", header=None, usecols=[0, 1, 2, 3], names=["chrom", "start", "end", "gene"])
    long = pd.read_csv(pav_path, sep="\t", usecols=["gene", "accession", "call"])
    calls = long.pivot(index="gene", columns="accession", values="call")
    assert set(bed["gene"]) == set(calls.index), "gene lists of the BED and the PAV table differ"
    calls = calls.loc[bed["gene"]]
    is_snp = (df["vtype"] == "SNP").to_numpy()
    alt = df[ORDER].to_numpy() >= 1

    dens = {}
    for k, a in enumerate(ORDER):
        arr = np.full(len(bed), np.nan)
        for c in CHROMS:
            gi = (bed["chrom"] == c).to_numpy()
            sel = (df["chrom"] == c).to_numpy() & is_snp & alt[:, k]
            pos = np.sort(df["pos"].to_numpy()[sel])
            arr[gi] = flank_density(pos, bed.loc[gi, "start"].to_numpy(), bed.loc[gi, "end"].to_numpy(),
                                    chrom_len[c], flank_bp)
        dens[a] = arr
    assert not any(np.isnan(v).any() for v in dens.values())

    groups = {"present": GREY, "ambiguous / uncertain": BLUE, "absent": ORANGE}

    def group_mask(a, g):
        c = calls[a].to_numpy()
        return (c == "present") if g == "present" else (c == "absent") if g == "absent" else np.isin(c, ["ambiguous", "uncertain"])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.2), gridspec_kw={"width_ratios": [1, 1.1], "wspace": 0.3})
    cap = np.percentile(np.concatenate([dens[a] for a in OTHERS]), 99)
    for g, colour in groups.items():
        vals = np.concatenate([dens[a][group_mask(a, g)] for a in OTHERS])
        v = np.sort(np.minimum(vals, cap))
        ax1.step(v, np.arange(1, len(v) + 1) / len(v), where="post", color=colour, lw=1.9,
                 label=f"{g} (n = {len(vals):,}, median {np.median(vals):.1f})")
        print(f"{g}: n = {len(vals)}, median flank density = {np.median(vals):.2f} SNPs/kb")
    ax1.set_xlabel(f"SNPs per kb in the {flank_bp // 1000} kb on each side of the gene (capped at the 99th percentile)")
    ax1.set_ylabel("fraction of genes (cumulative)")
    ax1.set_title("(a) Neighbourhood of present and absent genes")
    ax1.legend(frameon=False, loc="lower right", fontsize=8)
    ax1.grid(color=GRID, lw=0.6)
    ax1.set_axisbelow(True)

    ys = np.arange(len(OTHERS))
    med_p = np.array([np.median(dens[a][group_mask(a, "present")]) for a in OTHERS])
    med_a = np.array([np.median(dens[a][group_mask(a, "absent")]) for a in OTHERS])
    ax2.hlines(ys, med_p, med_a, color=LIGHT, lw=1.6)
    ax2.scatter(med_p, ys, s=55, color=GREY, zorder=3, label="present genes")
    ax2.scatter(med_a, ys, s=55, color=ORANGE, zorder=3, label="absent genes")
    for yi, p, a_ in zip(ys, med_p, med_a):
        ax2.text(max(p, a_) + 0.4, yi, f"{p:.1f} vs {a_:.1f}", va="center", fontsize=8, color=INK2)
    ax2.set_yticks(ys)
    ax2.set_yticklabels(OTHERS)
    ax2.invert_yaxis()
    ax2.set_xlim(0, max(med_p.max(), med_a.max()) * 1.3)
    ax2.set_xlabel("median SNPs per kb around the gene")
    ax2.set_title("(b) Per accession")
    ax2.legend(frameon=False, loc="upper left")  # left of the An-1 dots is empty
    ax2.grid(axis="x", color=GRID, lw=0.6)
    ax2.set_axisbelow(True)
    fig.suptitle("Are absent genes in more diverged regions? (descriptive; a large gap would mean some absent calls may be mapping failures)",
                 x=0.01, ha="left", fontsize=11)
    save(fig, out)


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gt", default="results/variants/filtered_gt.tsv.gz")
    ap.add_argument("--pav", default="results/pav/pav_long.tsv.gz")
    ap.add_argument("--bed", default="data/reference/genes.bed")
    ap.add_argument("--fai", default="data/reference/TAIR10.1_Col-0/GCF_000001735.4/GCF_000001735.4_TAIR10.1_genomic.fna.fai")
    ap.add_argument("--outdir", default="results/figures")
    ap.add_argument("--flank-bp", type=int, default=10_000, help="size of the window on each side of a gene")
    a = ap.parse_args()

    chrom_len = load_chrom_lengths(a.fai)
    df = load_genotypes(a.gt)
    os.makedirs(a.outdir, exist_ok=True)
    fig4(df, chrom_len, os.path.join(a.outdir, "13_variants__bar_chart_pca_distance_heatmap.png"))
    fig5(df, a.pav, a.bed, chrom_len, os.path.join(a.outdir, "14_variants_snp_density__ecdf_dumbbell_plot.png"), a.flank_bp)
    print("wrote 2 figures to", a.outdir)


if __name__ == "__main__":
    main()
