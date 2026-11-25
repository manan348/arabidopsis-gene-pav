#!/usr/bin/env python3
"""Figure 7: two ways a gene can be lost in an accession (Phase 5), run in env `map`.

Run after scripts/lof_genes.py:
    python scripts/plot_lof.py

Inputs (read only): results/pav/pav_long.tsv.gz, results/variants/lof_genes.tsv.gz
Output: results/figures/fig7_lost_vs_broken_genes.png (300 dpi)
    (a) genes called absent (deleted, no reads) vs genes called present but carrying a HIGH-impact variant
    (b) which kinds of HIGH effect the present-but-broken genes carry (a gene with two kinds counts in both)
Descriptive only: no statistical test. HIGH effects are predictions, not proof that the protein is lost.
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ORDER = ["Col-0", "An-1", "Kyoto", "Kondara", "Cvi-0", "Ler-1", "Se-0", "Ale-Stenar-44-4"]
EFFECTS = ["frameshift_variant", "stop_gained", "splice_donor_variant", "splice_acceptor_variant",
           "stop_lost", "start_lost"]
LABELS = ["frameshift", "stop gained", "splice donor", "splice acceptor", "stop lost", "start lost"]
EFF_COL = ["#2a78d6", "#eb6834", "#2a9d8f", "#8e6bbf", "#c9a227", "#898781"]

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
    ap.add_argument("--pav", default="results/pav/pav_long.tsv.gz")
    ap.add_argument("--lof", default="results/variants/lof_genes.tsv.gz")
    ap.add_argument("--out", default="results/figures/fig7_lost_vs_broken_genes.png")
    a = ap.parse_args()

    pav = pd.read_csv(a.pav, sep="\t", usecols=["gene", "accession", "call"])
    lof = pd.read_csv(a.lof, sep="\t")
    n_absent = pav[pav["call"] == "absent"].groupby("accession").size().reindex(ORDER).fillna(0).astype(int)
    broken = lof[lof["pav_call"] == "present"]
    n_broken = broken.groupby("accession").size().reindex(ORDER).fillna(0).astype(int)
    n_present = pav[pav["call"] == "present"].groupby("accession").size().reindex(ORDER)
    for x in ORDER:
        print(f"{x:<16} absent {n_absent[x]:4d}  present-but-broken {n_broken[x]:5d} "
              f"({100 * n_broken[x] / n_present[x]:.1f}% of present genes)")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 4.8), gridspec_kw={"width_ratios": [1.1, 1], "wspace": 0.3})

    # (a) lost by deletion vs broken by a variant
    xs = np.arange(len(ORDER))
    w = 0.38
    ax1.bar(xs - w / 2, n_absent.to_numpy(), w, color=ORANGE, label="absent (no reads, deleted)")
    ax1.bar(xs + w / 2, n_broken.to_numpy(), w, color=BLUE, label="present but with a HIGH-impact variant")
    for xi, v in zip(xs - w / 2, n_absent):
        ax1.text(xi, v + 40, f"{v}", ha="center", fontsize=7.5, color=INK2)
    for xi, v in zip(xs + w / 2, n_broken):
        ax1.text(xi, v + 40, f"{v:,}", ha="center", fontsize=7.5, color=INK2)
    ax1.set_xticks(xs)
    ax1.set_xticklabels([x.replace("Ale-Stenar-44-4", "Ale-S.-44-4") for x in ORDER], rotation=30, ha="right")
    ax1.set_ylabel("number of genes")
    ax1.set_ylim(0, n_broken.max() * 1.3)
    ax1.set_title("(a) Two ways to lose a gene")
    ax1.legend(frameon=False, loc="upper left", fontsize=8)
    ax1.grid(axis="y", color=GRID, lw=0.6)
    ax1.set_axisbelow(True)

    # (b) kinds of HIGH effect among the present-but-broken genes
    ex = broken.assign(effect=broken["effects"].str.split(";")).explode("effect")
    cnt = ex.groupby(["accession", "effect"]).size().unstack(fill_value=0).reindex(index=ORDER, columns=EFFECTS,
                                                                                    fill_value=0)
    share = 100 * cnt.div(cnt.sum(axis=1), axis=0)
    left = np.zeros(len(ORDER))
    for eff, lab, col in zip(EFFECTS, LABELS, EFF_COL):
        ax2.barh(xs, share[eff].to_numpy(), left=left, color=col, label=lab, height=0.7)
        left += share[eff].to_numpy()
    ax2.set_yticks(xs)
    ax2.set_yticklabels(ORDER)
    ax2.invert_yaxis()
    ax2.set_xlim(0, 100)
    ax2.set_xlabel("% of effect types among present-but-broken genes")
    ax2.set_title("(b) Kind of HIGH-impact effect")
    ax2.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=3, fontsize=8)

    fig.suptitle("Genes lost by deletion are few; genes predicted broken by a variant are 6-13 times more common",
                 x=0.01, ha="left", fontsize=11)
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    fig.savefig(a.out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
