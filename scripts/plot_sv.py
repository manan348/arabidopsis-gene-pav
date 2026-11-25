#!/usr/bin/env python3
"""Figure 6: Delly deletions versus the PAV absent calls (Phase 5), run in env `map`.

Run after scripts/sv_vs_pav.py:
    python scripts/plot_sv.py

Inputs (read only): results/sv/sv_vs_pav_summary.tsv, results/sv/pass_deletions.tsv
Output: results/figures/fig6_deletions_vs_absent_genes.png (300 dpi)
    (a) % of absent and of present genes that lie >= 50% inside a PASS deletion (log scale), fold difference
    (b) number of PASS deletions per accession vs the % of absent genes confirmed (Delly sensitivity)
    (c) deletion lengths (cumulative), with the 100 kb cut used in sv_vs_pav.py
Descriptive only: no statistical test.
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

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


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--summary", default="results/sv/sv_vs_pav_summary.tsv")
    ap.add_argument("--dels", default="results/sv/pass_deletions.tsv")
    ap.add_argument("--max-len", type=int, default=100_000, help="must match the value used in sv_vs_pav.py")
    ap.add_argument("--out", default="results/figures/fig6_deletions_vs_absent_genes.png")
    a = ap.parse_args()

    s = pd.read_csv(a.summary, sep="\t").set_index(["accession", "call"])
    d = pd.read_csv(a.dels, sep="\t", header=None, names=["accession", "chrom", "pos", "end"])
    d["len"] = d["end"] - d["pos"]
    assert (d["len"] > 0).all()

    pct_abs = np.array([s.loc[(x, "absent"), "pct_in_deletion"] for x in OTHERS])
    pct_pre = np.array([100 * s.loc[(x, "present"), "in_deletion"] / s.loc[(x, "present"), "genes"] for x in OTHERS])
    assert (pct_pre > 0).all(), "log scale needs a non-zero present share"
    fold = pct_abs / pct_pre
    n_del = d[d["len"] <= a.max_len].groupby("accession").size().reindex(OTHERS)
    for x, pa, pp, f, n in zip(OTHERS, pct_abs, pct_pre, fold, n_del):
        print(f"{x:<16} absent {pa:5.1f}%  present {pp:5.2f}%  fold {f:5.1f}  PASS deletions <= max-len {n}")

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 4.8), gridspec_kw={"width_ratios": [1.25, 1, 1], "wspace": 0.38})

    # (a) absent vs present genes inside deletions
    ys = np.arange(len(OTHERS))
    ax1.hlines(ys, pct_pre, pct_abs, color=LIGHT, lw=1.6)
    ax1.scatter(pct_pre, ys, s=55, color=GREY, zorder=3, label="present genes")
    ax1.scatter(pct_abs, ys, s=55, color=ORANGE, zorder=3, label="absent genes")
    for yi, pa, f in zip(ys, pct_abs, fold):
        ax1.text(pa * 1.35, yi, f"{f:.0f}x", va="center", fontsize=8, color=INK2)
    ax1.set_xscale("log")
    ax1.set_xlim(0.01, 200)
    ax1.set_xticks([0.01, 0.1, 1, 10, 100])
    ax1.set_xticklabels(["0.01", "0.1", "1", "10", "100"])
    ax1.set_yticks(ys)
    ax1.set_yticklabels(OTHERS)
    ax1.invert_yaxis()
    ax1.set_xlabel("% of genes >= 50% inside a Delly PASS deletion (log scale)")
    ax1.set_title("(a) Absent genes sit in deletions far more often")
    ax1.legend(frameon=False, loc="lower left", fontsize=8)
    ax1.grid(axis="x", color=GRID, lw=0.6)
    ax1.set_axisbelow(True)

    # (b) confirmation rate follows how many deletions Delly found
    ax2.scatter(n_del, pct_abs, s=55, color=BLUE, zorder=3)
    for x, n, pa in zip(OTHERS, n_del, pct_abs):
        ax2.annotate(x, (n, pa), xytext=(5, 4), textcoords="offset points", fontsize=8, color=INK2)
    ax2.set_xlim(0, n_del.max() * 1.3)
    ax2.set_ylim(0, pct_abs.max() * 1.2)
    ax2.set_xlabel(f"PASS deletions found by Delly (<= {a.max_len // 1000} kb)")
    ax2.set_ylabel("% of absent genes inside a deletion")
    ax2.set_title("(b) Confirmation depends on Delly's yield")
    ax2.grid(color=GRID, lw=0.6)
    ax2.set_axisbelow(True)

    # (c) deletion lengths
    for x in OTHERS:
        v = np.sort(d.loc[d["accession"] == x, "len"].to_numpy())
        ax3.step(v, np.arange(1, len(v) + 1) / len(v), where="post", lw=1.2, color=BLUE, alpha=0.55)
    ax3.axvline(a.max_len, color=ORANGE, lw=1.2, ls="--")
    ax3.text(a.max_len * 1.15, 0.05, f"{a.max_len // 1000} kb cut\n({(d['len'] > a.max_len).sum()} dropped)",
             fontsize=8, color=ORANGE)
    ax3.set_xscale("log")
    ax3.set_xlabel("deletion length (bp, log scale)")
    ax3.set_ylabel("fraction of deletions (cumulative)")
    ax3.set_title(f"(c) Deletion lengths (one line per accession, median {d['len'].median():,.0f} bp)")
    ax3.grid(color=GRID, lw=0.6)
    ax3.set_axisbelow(True)

    fig.suptitle("Do the read-depth absent calls agree with Delly deletions (split and discordant reads)?",
                 x=0.01, ha="left", fontsize=11)
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    fig.savefig(a.out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
