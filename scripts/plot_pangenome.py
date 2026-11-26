#!/usr/bin/env python3
"""Figures 8 and 9: shared absences (UpSet) and gene accumulation curve (Phase 7), run in env `map`.

Run:
    python scripts/plot_pangenome.py

Input (read only): results/pav/pav_long.tsv.gz  (gene, accession, call)
Outputs: results/figures/28_pangenome_shared_absences__upset_plot.png, results/figures/29_pangenome_core_genes__accumulation_rarefaction_curve.png (300 dpi)

A gene counts as ABSENT in an accession only when its call is `absent`. Ambiguous and uncertain calls are
not evidence of absence. The accumulation curve shows two lines: "not absent in any" (lenient) and
"called present in all" (strict; this one ends at the 24,259 core genes reported earlier).
The study is reference-based, so genes that exist in other accessions but not in Col-0 cannot be seen.
Eight accessions are a small set: the curve shows a rough shape, not a precise size of the core genome.
Descriptive only: no statistical test.
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

ORDER = ["Col-0", "An-1", "Kyoto", "Kondara", "Cvi-0", "Ler-1", "Se-0", "Ale-Stenar-44-4"]
OTHERS = ORDER[1:]

BLUE, ORANGE, GREY, LIGHT = "#2a78d6", "#eb6834", "#898781", "#c3c2b7"
INK, INK2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e1e0d9", "#fcfcfb"

plt.rcParams.update({
    "font.size": 9, "axes.edgecolor": LIGHT, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "axes.titlesize": 10, "axes.titleweight": "bold", "axes.titlelocation": "left",
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": False,
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
})


def upset(call, out, top_n):
    absent = (call[OTHERS] == "absent")
    n_genes = len(call)
    var = absent[absent.any(axis=1)]
    sizes = var.sum().sort_values(ascending=False)
    sets = list(sizes.index)  # rows of the dot matrix, largest set on top
    pat = var[sets].apply(tuple, axis=1).value_counts()
    deg = pd.Series({p: sum(p) for p in pat.index})
    tab = pd.DataFrame({"n": pat, "deg": deg}).sort_values(["n", "deg"], ascending=[False, True])
    print(f"variable genes (absent in >= 1 of the 7): {len(var)} of {n_genes}")
    print(f"distinct absence patterns: {len(tab)}; private to one accession: "
          f"{int(tab.loc[tab['deg'] == 1, 'n'].sum())}; shared by 2+: {int(tab.loc[tab['deg'] >= 2, 'n'].sum())}; "
          f"absent in all 7: {int(tab.loc[tab['deg'] == 7, 'n'].sum())}")
    show = tab.head(top_n)
    print(f"showing the {len(show)} largest patterns ({int(show['n'].sum())} of {len(var)} variable genes)")

    fig = plt.figure(figsize=(11, 6.2))
    gs = fig.add_gridspec(2, 3, width_ratios=[1, 0.6, 3.4], height_ratios=[2.3, 1.6], wspace=0.03, hspace=0.04)
    ax_b = fig.add_subplot(gs[0, 2])
    ax_m = fig.add_subplot(gs[1, 2], sharex=ax_b)
    ax_l = fig.add_subplot(gs[1, 1], sharey=ax_m)  # column with the accession names
    ax_s = fig.add_subplot(gs[1, 0], sharey=ax_m)

    xs = np.arange(len(show))
    cols = [BLUE if d == 1 else ORANGE for d in show["deg"]]
    ax_b.bar(xs, show["n"], color=cols, width=0.7)
    for xi, v in zip(xs, show["n"]):
        ax_b.text(xi, v + show["n"].max() * 0.015, f"{int(v)}", ha="center", va="bottom", fontsize=8, color=INK2)
    ax_b.set_ylim(0, show["n"].max() * 1.12)
    ax_b.set_ylabel("genes absent in exactly\nthis set of accessions")
    ax_b.set_title("Shared absences: which accessions lack the same genes")
    ax_b.tick_params(axis="x", bottom=False, labelbottom=False)
    ax_b.grid(axis="y", color=GRID, lw=0.6)
    ax_b.set_axisbelow(True)
    ax_b.legend(handles=[Patch(color=BLUE, label="private to one accession"),
                         Patch(color=ORANGE, label="shared by 2 or more")],
                frameon=False, loc="upper right", fontsize=8)

    ys = np.arange(len(sets))
    for xi, p in zip(xs, show.index):
        on = [i for i, v in enumerate(p) if v]
        ax_m.scatter([xi] * len(sets), ys, s=34, color=LIGHT, zorder=2)
        ax_m.scatter([xi] * len(on), on, s=40, color=INK, zorder=3)
        if len(on) > 1:
            ax_m.plot([xi, xi], [min(on), max(on)], color=INK, lw=2, zorder=3)
    ax_m.set_yticks(ys)
    ax_m.set_yticklabels([])
    ax_m.set_ylim(len(sets) - 0.5, -0.5)
    ax_m.set_xlim(-0.6, len(show) - 0.4)
    ax_m.set_xticks([])
    for sp in ["left", "bottom"]:
        ax_m.spines[sp].set_visible(False)
    ax_m.tick_params(left=False)

    ax_s.barh(ys, sizes.to_numpy(), color=GREY, height=0.6)
    for yi, v in zip(ys, sizes):
        ax_s.text(v, yi, f" {int(v)}", va="center", ha="right", fontsize=8, color=INK2)
    ax_s.invert_xaxis()
    ax_s.set_yticks([])
    ax_s.set_xlabel("genes absent in the accession")
    ax_s.spines["left"].set_visible(False)
    ax_s.spines["right"].set_visible(False)
    ax_l.set_xlim(0, 1)
    for yi, name in zip(ys, sets):
        ax_l.text(0.97, yi, name, ha="right", va="center", fontsize=9, color=INK)
    ax_l.axis("off")

    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out)
    return absent


def accumulation(call, absent7, out, n_perm, seed):
    A = (call[ORDER] == "absent").to_numpy()
    P = (call[ORDER] == "present").to_numpy()
    n_genes = len(call)
    rng = np.random.default_rng(seed)
    lenient = np.zeros((n_perm, len(ORDER)))
    strict = np.zeros((n_perm, len(ORDER)))
    for i in range(n_perm):
        perm = rng.permutation(len(ORDER))
        lenient[i] = (~np.logical_or.accumulate(A[:, perm], axis=1)).sum(axis=0)
        strict[i] = (~np.logical_or.accumulate(~P[:, perm], axis=1)).sum(axis=0)
    ks = np.arange(1, len(ORDER) + 1)
    print(f"accumulation, {n_perm} random orders; with all 8 accessions: not absent in any {int(lenient[0, -1])}, "
          f"called present in all {int(strict[0, -1])} (of {n_genes})")
    for k in ks:
        print(f"  k={k}: lenient median {np.median(lenient[:, k - 1]):.0f}, strict median {np.median(strict[:, k - 1]):.0f}")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.6), gridspec_kw={"width_ratios": [1.2, 1], "wspace": 0.28})
    for arr, col, lab in [(lenient, BLUE, "not called absent in any"), (strict, ORANGE, "called present in all")]:
        med = np.median(arr, axis=0)
        ax1.fill_between(ks, np.percentile(arr, 5, axis=0), np.percentile(arr, 95, axis=0), color=col, alpha=0.2, lw=0)
        ax1.plot(ks, med, color=col, lw=2, marker="o", ms=5, label=lab)
        ax1.text(ks[-1] + 0.12, med[-1], f"{med[-1]:,.0f}", va="center", fontsize=8, color=INK2)
    ax1.axhline(n_genes, color=LIGHT, lw=1, ls="--")
    ax1.text(1.0, n_genes + 80, f"all {n_genes:,} genes", fontsize=8, color=INK2)
    ax1.set_xlim(0.8, len(ORDER) + 0.9)
    ax1.set_xticks(ks)
    ax1.set_ylim(22000, n_genes + 500)
    ax1.set_xlabel("number of accessions added (random order, line = median, band = 5-95%)")
    ax1.set_ylabel("core genes (y axis starts at 22,000)")
    ax1.set_title("(a) Core genes shrink as accessions are added")
    ax1.legend(frameon=False, loc="lower left", fontsize=8)
    ax1.grid(axis="y", color=GRID, lw=0.6)
    ax1.set_axisbelow(True)

    n_abs = absent7.sum(axis=1)
    spec = n_abs[n_abs > 0].value_counts().reindex(range(1, 8), fill_value=0)
    ax2.bar(spec.index, spec.to_numpy(), color=BLUE, width=0.7)
    for k, v in spec.items():
        ax2.text(k, v + spec.max() * 0.015, f"{int(v)}", ha="center", va="bottom", fontsize=8, color=INK2)
    ax2.set_ylim(0, spec.max() * 1.12)
    ax2.set_xticks(range(1, 8))
    ax2.set_xlabel("number of the 7 non-Col-0 accessions in which the gene is absent")
    ax2.set_ylabel("number of genes")
    ax2.set_title("(b) Most variable genes are absent in few accessions")
    ax2.grid(axis="y", color=GRID, lw=0.6)
    ax2.set_axisbelow(True)

    fig.suptitle("Reference-based gene presence/absence across 8 accessions (genes missing from Col-0 cannot be seen)",
                 x=0.01, ha="left", fontsize=11)
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pav", default="results/pav/pav_long.tsv.gz")
    ap.add_argument("--out8", default="results/figures/28_pangenome_shared_absences__upset_plot.png")
    ap.add_argument("--out9", default="results/figures/29_pangenome_core_genes__accumulation_rarefaction_curve.png")
    ap.add_argument("--top", type=int, default=20, help="number of absence patterns shown in the UpSet plot")
    ap.add_argument("--perms", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()

    long_ = pd.read_csv(a.pav, sep="\t", usecols=["gene", "accession", "call"])
    call = long_.pivot(index="gene", columns="accession", values="call")[ORDER]
    assert call.notna().all().all(), "missing calls in the PAV table"
    assert (call["Col-0"] != "absent").all(), "the control should have no absent genes"

    absent7 = upset(call, a.out8, a.top)
    accumulation(call, absent7, a.out9, a.perms, a.seed)


if __name__ == "__main__":
    main()
