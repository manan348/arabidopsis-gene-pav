#!/usr/bin/env python3
"""More PAV heatmaps (Phase 4). Run from the project root after scripts/build_pav.py:
    python scripts/plot_pav_heatmaps.py [--pavdir results/pav] [--outdir results/figures]

Writes (PNG, 200 dpi):
  11_pav_shared_absences__jaccard_similarity_matrix.png  7 x 7: shared absent genes between accessions (Col-0 left out, it has none)
  12_pav_along_chromosomes__density_heatmap.png        absent genes per 1 Mb bin along chr1-5, one row per accession
  10_pav_relative_depth__diverging_heatmap.png         depth / sample median for the variable genes (red = low, white = 1, blue = high)
  09_pav_clustered__clustered_heatmap_dendrogram.png     variable genes clustered with dendrograms (only if scipy is installed)
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap, TwoSlopeNorm
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch

ORDER = ["Col-0", "An-1", "Kyoto", "Kondara", "Cvi-0", "Ler-1", "Se-0", "Ale-Stenar-44-4"]
CHROMS = {"NC_003070.9": "Chr1", "NC_003071.7": "Chr2", "NC_003074.8": "Chr3",
          "NC_003075.7": "Chr4", "NC_003076.8": "Chr5"}
CODE = {"present": 0, "ambiguous": 1, "uncertain": 2, "absent": 3}
CALL_CMAP = ListedColormap(["#e6e6e6", "#fdb863", "#5e8fc7", "#b2182b"])
CALL_LEGEND = [Patch(color="#e6e6e6", label="present"), Patch(color="#fdb863", label="ambiguous"),
               Patch(color="#5e8fc7", label="uncertain"), Patch(color="#b2182b", label="absent")]
SAVE = dict(dpi=200, bbox_inches="tight")


def variable_genes(long):
    """Return (call codes, norm depth) for genes with >= 1 absence, sorted by number/pattern of absences."""
    calls = long.pivot(index="gene", columns="accession", values="call")[ORDER]
    depth = long.pivot(index="gene", columns="accession", values="norm_depth")[ORDER]
    keep = (calls == "absent").any(axis=1)
    num = calls[keep].apply(lambda col: col.map(CODE)).astype(int)
    key = pd.DataFrame({"n": -(num == 3).sum(axis=1)})
    for c in ORDER:
        key[c] = -(num[c] == 3).astype(int)
    idx = key.sort_values(["n"] + ORDER).index
    return num.loc[idx], depth.loc[idx]


def fig_similarity(long, out):
    accs = ORDER[1:]
    sets = {a: set(long.loc[(long["accession"] == a) & (long["call"] == "absent"), "gene"]) for a in accs}
    n = len(accs)
    shared = np.zeros((n, n), dtype=int)
    jac = np.zeros((n, n))
    for i, a in enumerate(accs):
        for j, b in enumerate(accs):
            inter = len(sets[a] & sets[b])
            union = len(sets[a] | sets[b])
            shared[i, j] = inter
            jac[i, j] = inter / union if union else 0
    fig, ax = plt.subplots(figsize=(6.8, 5.8))
    im = ax.imshow(jac, cmap="YlOrRd", vmin=0, vmax=1)
    for i in range(n):
        for j in range(n):
            ax.text(j, i, shared[i, j], ha="center", va="center", fontsize=9,
                    color="white" if jac[i, j] > 0.55 else "#222222")
    ax.set_xticks(range(n))
    ax.set_xticklabels(accs, rotation=40, ha="right")
    ax.set_yticks(range(n))
    ax.set_yticklabels(accs)
    ax.set_title("Absent genes shared between accessions\n(number = shared genes, diagonal = total; colour = Jaccard)", fontsize=10)
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("Jaccard similarity of absent-gene sets")
    fig.savefig(out, **SAVE)
    plt.close(fig)


def fig_genome(long, bed_path, out, bin_mb=1.0):
    bed = pd.read_csv(bed_path, sep="\t", header=None, usecols=[0, 1, 2, 3],
                      names=["chrom", "start", "end", "gene"])
    bed["mid"] = (bed["start"] + bed["end"]) / 2
    absent = long.loc[long["call"] == "absent", ["gene", "accession"]].merge(bed, on="gene")
    nbins = {c: int(bed.loc[bed["chrom"] == c, "end"].max() // (bin_mb * 1e6)) + 1 for c in CHROMS}
    mats = {}
    for c in CHROMS:
        m = np.zeros((len(ORDER), nbins[c]))
        sub = absent[absent["chrom"] == c]
        for i, a in enumerate(ORDER):
            b = (sub.loc[sub["accession"] == a, "mid"] // (bin_mb * 1e6)).astype(int)
            for k, v in b.value_counts().items():
                m[i, k] = v
        mats[c] = m
    vmax = max(m.max() for m in mats.values())
    fig = plt.figure(figsize=(14, 3.6))
    gs = GridSpec(1, 6, width_ratios=[nbins[c] for c in CHROMS] + [1.2], wspace=0.05)
    first = None
    for k, c in enumerate(CHROMS):
        ax = fig.add_subplot(gs[0, k], sharey=first)
        first = first or ax
        im = ax.imshow(mats[c], cmap="Reds", vmin=0, vmax=vmax, aspect="auto", interpolation="nearest")
        ax.set_title(CHROMS[c], fontsize=10)
        ax.set_xticks(range(0, nbins[c], 5))
        ax.set_xlabel("Mb")
        if k == 0:
            ax.set_yticks(range(len(ORDER)))
            ax.set_yticklabels(ORDER)
        else:
            plt.setp(ax.get_yticklabels(), visible=False)
            ax.tick_params(axis="y", length=0)
    cax = fig.add_subplot(gs[0, 5])
    cax.axis("off")
    cb = fig.colorbar(im, ax=cax, fraction=0.5, aspect=12)
    cb.set_label(f"absent genes per {bin_mb:g} Mb")
    fig.suptitle("Where the absent genes sit along the chromosomes (gene midpoint, Col-0 coordinates)", fontsize=11, y=1.04)
    fig.savefig(out, **SAVE)
    plt.close(fig)


def fig_depth(depth, out):
    fig, ax = plt.subplots(figsize=(6.5, 9))
    norm = TwoSlopeNorm(vmin=0, vcenter=1, vmax=2)
    im = ax.imshow(depth.clip(0, 2).values, aspect="auto", cmap="RdBu", norm=norm, interpolation="nearest")
    ax.set_xticks(range(len(ORDER)))
    ax.set_xticklabels(ORDER, rotation=40, ha="right")
    ax.set_yticks([])
    ax.set_ylabel(f"{len(depth)} genes absent in at least one accession")
    ax.set_title("Depth of the variable genes relative to each sample's median\n(red = low, white = normal, blue = extra copies?)", fontsize=10)
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.03)
    cb.set_label("gene depth / sample median (capped at 2)")
    fig.savefig(out, **SAVE)
    plt.close(fig)


def fig_clustered(num, out):
    try:
        from scipy.cluster import hierarchy as sch
    except ImportError:
        print("scipy not installed: skipping the clustered heatmap")
        return
    X = (num.values == 3).astype(float)  # 1 = absent
    Zr = sch.linkage(X, method="average", metric="cityblock")
    Zc = sch.linkage(X.T, method="average", metric="cityblock")
    ridx = sch.dendrogram(Zr, no_plot=True)["leaves"]
    cidx = sch.dendrogram(Zc, no_plot=True)["leaves"]
    nr, nc = X.shape
    fig = plt.figure(figsize=(7.5, 9.5))
    gs = GridSpec(2, 2, width_ratios=[1.2, 6], height_ratios=[1.2, 12], wspace=0.01, hspace=0.01)
    ax_top = fig.add_subplot(gs[0, 1])
    ax_left = fig.add_subplot(gs[1, 0])
    ax_h = fig.add_subplot(gs[1, 1])
    grey = lambda k: "#555555"
    sch.dendrogram(Zc, ax=ax_top, orientation="top", no_labels=True, link_color_func=grey)
    sch.dendrogram(Zr, ax=ax_left, orientation="left", no_labels=True, link_color_func=grey)
    for a in (ax_top, ax_left):
        a.axis("off")
    ax_top.set_xlim(0, 10 * nc)
    ax_left.set_ylim(0, 10 * nr)
    data = num.values[np.ix_(ridx, cidx)]
    ax_h.imshow(data, aspect="auto", cmap=CALL_CMAP, vmin=-0.5, vmax=3.5,
                interpolation="nearest", origin="lower", extent=(0, nc, 0, nr))
    ax_h.set_xticks(np.arange(nc) + 0.5)
    ax_h.set_xticklabels([ORDER[i] for i in cidx], rotation=40, ha="right")
    ax_h.set_yticks([])
    ax_h.yaxis.tick_right()
    ax_h.legend(handles=CALL_LEGEND, loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False)
    fig.suptitle(f"Clustered PAV of the {nr} variable genes", fontsize=11, y=0.94)
    fig.savefig(out, **SAVE)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pavdir", default="results/pav")
    ap.add_argument("--outdir", default="results/figures")
    ap.add_argument("--bed", default="data/reference/genes.bed")
    a = ap.parse_args()
    long = pd.read_csv(os.path.join(a.pavdir, "pav_long.tsv.gz"), sep="\t")
    os.makedirs(a.outdir, exist_ok=True)
    num, depth = variable_genes(long)
    fig_similarity(long, os.path.join(a.outdir, "11_pav_shared_absences__jaccard_similarity_matrix.png"))
    fig_genome(long, a.bed, os.path.join(a.outdir, "12_pav_along_chromosomes__density_heatmap.png"))
    fig_depth(depth, os.path.join(a.outdir, "10_pav_relative_depth__diverging_heatmap.png"))
    fig_clustered(num, os.path.join(a.outdir, "09_pav_clustered__clustered_heatmap_dendrogram.png"))
    print("done, files in", a.outdir)


if __name__ == "__main__":
    main()
