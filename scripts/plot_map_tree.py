#!/usr/bin/env python3
"""Figure 32: where the 8 accessions come from, and how related they are (SNP tree), env `map`.

Run:
    python scripts/plot_map_tree.py

Inputs (read only): data/meta/accession_coords.tsv (1001 Genomes coordinates, from the AraPheno API),
    results/analysis/snp_distance_matrix.tsv (fig 31), results/pav/pav_long.tsv.gz,
    data/reference/ne_110m_land/ne_110m_land.shp (Natural Earth 1:110m land, public domain; for the outline only)
Output: results/figures/32_relatedness_origins__geographic_map_upgma_dendrogram.png (300 dpi)

Tree: UPGMA (average linkage) on the pairwise SNP distance (share of differing sites). Col-0 is the reference, so
its distances are slightly lower than they would be for a non-reference line. With 8 accessions this shows
relatedness only, not population structure. Col-0 coordinates are the 1001 Genomes record (lab line, Columbia, US).
Point size on the map = number of absent genes. Descriptive only: no statistical test.
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform

ORDER = ["Col-0", "An-1", "Kyoto", "Kondara", "Cvi-0", "Ler-1", "Se-0", "Ale-Stenar-44-4"]
COL = {"Col-0": "#52514e", "An-1": "#2a78d6", "Ler-1": "#1c4f8f", "Ale-Stenar-44-4": "#7fb3ea",
       "Se-0": "#e3b23c", "Cvi-0": "#eb6834", "Kondara": "#8e6bbf", "Kyoto": "#2a9d8f"}
LIGHT, INK, INK2, GRID, SURFACE = "#c3c2b7", "#0b0b0b", "#52514e", "#e1e0d9", "#fcfcfb"
plt.rcParams.update({
    "font.size": 9, "axes.edgecolor": LIGHT, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "axes.titlesize": 10, "axes.titleweight": "bold", "axes.titlelocation": "left",
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
})


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--coords", default="data/meta/accession_coords.tsv")
    ap.add_argument("--dist", default="results/analysis/snp_distance_matrix.tsv")
    ap.add_argument("--pav", default="results/pav/pav_long.tsv.gz")
    ap.add_argument("--land", default="data/reference/ne_110m_land/ne_110m_land.shp")
    ap.add_argument("--out", default="results/figures/32_relatedness_origins__geographic_map_upgma_dendrogram.png")
    a = ap.parse_args()

    c = pd.read_csv(a.coords, sep="\t").set_index("accession").loc[ORDER]
    D = pd.read_csv(a.dist, sep="\t", index_col=0).loc[ORDER, ORDER]
    p = pd.read_csv(a.pav, sep="\t", usecols=["accession", "call"])
    nabs = (p["call"] == "absent").groupby(p["accession"]).sum().reindex(ORDER)

    fig = plt.figure(figsize=(14, 5.6))
    axm = fig.add_axes([0.01, 0.08, 0.64, 0.8])
    axt = fig.add_axes([0.71, 0.12, 0.27, 0.72])

    land = gpd.read_file(a.land)
    land.plot(ax=axm, color="#ecebe6", edgecolor=LIGHT, linewidth=0.4)
    for acc, r in c.iterrows():
        s = 30 + 0.45 * nabs[acc]
        axm.scatter(r["longitude"], r["latitude"], s=s, color=COL[acc], edgecolor="white", lw=0.8, zorder=3)
    off = {"Col-0": (6, -12), "An-1": (-8, 8), "Ler-1": (6, -14), "Ale-Stenar-44-4": (6, 6), "Se-0": (-8, -14),
           "Cvi-0": (8, -4), "Kondara": (6, -16), "Kyoto": (-6, 10)}
    for acc, r in c.iterrows():
        dx, dy = off[acc]
        axm.annotate(f"{acc} ({r['country'].replace('GER', 'Germany')}, {nabs[acc]} absent)",
                     (r["longitude"], r["latitude"]), xytext=(dx, dy), textcoords="offset points", fontsize=7.8,
                     color=INK, ha="right" if dx < 0 else "left")
    axm.set_xlim(-125, 150)
    axm.set_ylim(0, 72)
    axm.set_aspect("equal")
    axm.set_xticks([])
    axm.set_yticks([])
    for s in axm.spines.values():
        s.set_visible(False)
    axm.set_title("(a) Collection sites (point size = number of absent genes)")

    Z = linkage(squareform(D.to_numpy(), checks=False), method="average")
    dn = dendrogram(Z, labels=ORDER, orientation="left", ax=axt, color_threshold=0, above_threshold_color=INK2,
                    leaf_font_size=9)
    for t in axt.get_yticklabels():
        t.set_color(COL[t.get_text()])
        t.set_fontweight("bold")
    axt.set_xlabel("SNP distance (share of differing sites)")
    axt.set_title("(b) UPGMA tree of SNP distance")
    for s in ["top", "right", "left"]:
        axt.spines[s].set_visible(False)
    axt.tick_params(axis="y", length=0)
    axt.grid(axis="x", color=GRID, lw=0.6)
    axt.set_axisbelow(True)
    print("leaf order (top to bottom):", list(reversed(dn["ivl"])))
    fig.suptitle("Eight accessions: origin and relatedness", x=0.01, ha="left", fontsize=12, fontweight="bold")
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    fig.savefig(a.out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
