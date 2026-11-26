#!/usr/bin/env python3
"""Figure 36: GO chord plot, enriched GO terms and the absent genes behind them, env `map`.

Run (after scripts/go_enrichment.R):
    python scripts/plot_go_chord.py            # genes missing in >= 3 accessions (default)
    python scripts/plot_go_chord.py --min 2    # more genes

Inputs (read only): results/go/go_enrichment.tsv, results/pav/pav_long.tsv.gz
Output: results/figures/36_function_go_genes__chord_diagram.png (300 dpi)

Left: absent-set genes in 6 enriched terms, filtered to genes missing in >= --min accessions so the labels stay
readable; the box next to each gene shows how many of the 7 accessions miss it (stands in for the logFC scale of a
classic GOChord). Right: the GO terms, arc length = number of shown genes. A ribbon joins a gene to each of its terms.
Gene label = TAIR symbol, or AT ID when there is none. Same drawing idea as GOplot::GOChord, drawn with matplotlib.
Descriptive only: no new test (term FDRs are from go_enrichment.R).
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap
from matplotlib.patches import PathPatch, Wedge
from matplotlib.path import Path

TERMS = [("ADP binding", "ADP binding (NLR resistance genes)", "#2a78d6"),
         ("NAD+ nucleosidase activity, cyclic ADP-ribose generating", "NAD+ nucleosidase (TIR-domain NLRs)", "#1c4f8f"),
         ("defense response to other organism", "defense response to other organism", "#7fb3ea"),
         ("sesquiterpene biosynthetic process", "sesquiterpene biosynthesis (terpene synthases)", "#e3b23c"),
         ("RNA-DNA hybrid ribonuclease activity", "RNA-DNA hybrid ribonuclease (RNase H, TE-like)", "#8e6bbf"),
         ("protein localization to cell surface", "protein localization to cell surface (small cluster)", "#2a9d8f")]
INK, INK2, SURFACE = "#0b0b0b", "#52514e", "#fcfcfb"
R = 1.0


def arc_pts(t0, t1, r=R, n=40):
    t = np.linspace(t0, t1, n)
    return np.c_[r * np.cos(t), r * np.sin(t)]


def ribbon(ax, g0, g1, t0, t1, color):
    """Ribbon from gene slot [g0, g1] to term slot [t0, t1] (angles, radians), bent through the centre."""
    verts, codes = [], []
    a = arc_pts(g0, g1)
    verts += a.tolist()
    codes += [Path.MOVETO] + [Path.LINETO] * (len(a) - 1)
    verts += [[0, 0], [R * np.cos(t0), R * np.sin(t0)]]
    codes += [Path.CURVE3, Path.CURVE3]
    b = arc_pts(t0, t1)
    verts += b[1:].tolist()
    codes += [Path.LINETO] * (len(b) - 1)
    verts += [[0, 0], [R * np.cos(g0), R * np.sin(g0)]]
    codes += [Path.CURVE3, Path.CURVE3]
    ax.add_patch(PathPatch(Path(verts, codes), fc=color, ec="white", lw=0.25, alpha=0.78, zorder=2))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--go", default="results/go/go_enrichment.tsv")
    ap.add_argument("--pav", default="results/pav/pav_long.tsv.gz")
    ap.add_argument("--min", type=int, default=3, help="keep genes missing in at least this many accessions")
    ap.add_argument("--out", default="results/figures/36_function_go_genes__chord_diagram.png")
    a = ap.parse_args()

    go = pd.read_csv(a.go, sep="\t")
    go = go[go["set"] == "absent"]
    p = pd.read_csv(a.pav, sep="\t", usecols=["gene", "call"])
    nmiss = (p["call"] == "absent").groupby(p["gene"]).sum()

    label, member = {}, {}
    for j, (term, _, _) in enumerate(TERMS):
        row = go[go["Description"] == term].iloc[0]
        for gid, sym in zip(row["geneID"].split("/"), row["symbols"].split("/")):
            if nmiss[gid] >= a.min:
                label[gid] = sym
                member.setdefault(gid, []).append(j)
    genes = sorted(member, key=lambda g: (member[g][0], -nmiss[g], label[g]))
    n = len(genes)
    tcount = [sum(j in member[g] for g in genes) for j in range(len(TERMS))]
    print(f"{n} genes (missing in >= {a.min} accessions); genes per term: "
          + ", ".join(f"{s.split(' (')[0]} {c}" for (_, s, _), c in zip(TERMS, tcount)))

    # angles: genes on the left (top to bottom), terms on the right (bottom to top), gaps between
    g_start, g_end = np.deg2rad(97), np.deg2rad(263)
    gw = (g_end - g_start) / n
    t_start, t_end = np.deg2rad(-80), np.deg2rad(80)
    gap = np.deg2rad(2.5)
    tw_unit = (t_end - t_start - gap * (len(TERMS) - 1)) / sum(tcount)
    # terms from the top of the right side going down, in the same order as the gene blocks on the left,
    # so ribbons run across instead of crossing; each term's slots are filled from its top
    t_rng, cur = [], t_end
    for c in tcount:
        t_rng.append((cur - c * tw_unit, cur))
        cur -= c * tw_unit + gap
    t_next = [r[1] for r in t_rng]

    fig = plt.figure(figsize=(11, 11.4))
    ax = fig.add_axes([0.02, 0.13, 0.96, 0.8])
    ax.set_xlim(-1.75, 1.75)
    ax.set_ylim(-1.55, 1.55)
    ax.set_aspect("equal")
    ax.axis("off")
    cmap = ListedColormap(plt.cm.YlOrRd(np.linspace(0.15, 0.95, 7)))

    for i, g in enumerate(genes):
        a0 = g_start + i * gw
        a1 = a0 + gw * 0.88
        k = len(member[g])
        for s, j in enumerate(member[g]):
            sub0 = a0 + s * (a1 - a0) / k
            sub1 = sub0 + (a1 - a0) / k
            t_next[j] -= tw_unit
            tt0 = t_next[j]
            ribbon(ax, sub0, sub1, tt0 + tw_unit * 0.06, tt0 + tw_unit * 0.94, TERMS[j][2])
        mid = (a0 + a1) / 2
        ax.add_patch(Wedge((0, 0), R + 0.07, np.rad2deg(a0), np.rad2deg(a1), width=0.06,
                           fc=cmap((nmiss[g] - 1) / 6), ec="none", zorder=3))
        deg = np.rad2deg(mid)
        ax.text((R + 0.1) * np.cos(mid), (R + 0.1) * np.sin(mid), label[g], rotation=deg - 180, rotation_mode="anchor",
                ha="right", va="center", fontsize=7.3, color=INK)
    for j, (t0, t1) in enumerate(t_rng):
        ax.add_patch(Wedge((0, 0), R + 0.07, np.rad2deg(t0), np.rad2deg(t1), width=0.07, fc=TERMS[j][2], ec="none",
                           zorder=3))
        mid = (t0 + t1) / 2
        ax.text((R + 0.11) * np.cos(mid), (R + 0.11) * np.sin(mid), str(tcount[j]), ha="left" if np.cos(mid) > 0
                else "right", va="center", fontsize=8.5, fontweight="bold", color=TERMS[j][2])

    h = [plt.Rectangle((0, 0), 1, 1, color=c) for _, _, c in TERMS]
    fig.legend(h, [f"{s} - {c} genes" for (_, s, _), c in zip(TERMS, tcount)], title="GO term (FDR < 0.05)",
               loc="lower left", bbox_to_anchor=(0.04, 0.03), ncol=1, frameon=False, fontsize=8.5,
               title_fontsize=9)
    cax = fig.add_axes([0.70, 0.08, 0.24, 0.016])
    cb = plt.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0.5, 7.5)), cax=cax,
                      orientation="horizontal", ticks=range(1, 8))
    cb.set_label("accessions missing the gene (of 7)", fontsize=8.5)
    cb.ax.tick_params(labelsize=8)
    cb.outline.set_visible(False)
    fig.suptitle("Enriched GO terms and the absent genes behind them", x=0.03, y=0.975, ha="left", fontsize=13,
                 fontweight="bold")
    fig.text(0.03, 0.948, f"Genes missing in at least {a.min} of 7 accessions; a ribbon joins each gene to every "
             "term it carries. Gene names: TAIR symbol, or AT ID if none.", fontsize=8.5, color=INK2)
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    fig.savefig(a.out, dpi=300, facecolor=SURFACE)
    plt.close(fig)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
