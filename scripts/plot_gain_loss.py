#!/usr/bin/env python3
"""Figure 26: genes lost, broken or (possibly) gained per accession, run in env `map`.

Run:
    python scripts/plot_gain_loss.py

Inputs (read only): results/pav/pav_long.tsv.gz, results/variants/lof_genes.tsv.gz
Outputs: results/figures/26_geneloss_lost_broken_gained__bar_chart_error_bars_scatter.png (300 dpi), results/pav/gene_gains.tsv (gain candidates)

  lost    = read call `absent` (Phase 4)
  broken  = read call `present` and a HIGH-impact variant (SnpEff, Phase 5)
  gained  = possible extra copies: normalised read depth (MAPQ >= 10, divided by the accession's median) >= --gain
            (default 2.0) in the accession, breadth >= 0.8, while the gene looks single-copy in Col-0
            (normalised depth 0.7-1.3, breadth >= 0.95). Whiskers show the count at 2.5x (low) and 1.8x (high).
CAUTION: read-depth gains at 7-12x depth are noisy and depend on each sample's depth spread (short reads and
many duplicates widen it); they are candidates, not confirmed duplications. Descriptive only: no statistical test.
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ORDER = ["Col-0", "An-1", "Kyoto", "Kondara", "Cvi-0", "Ler-1", "Se-0", "Ale-Stenar-44-4"]
BLUE, ORANGE, GREY, LIGHT = "#2a78d6", "#eb6834", "#898781", "#c3c2b7"
INK, INK2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e1e0d9", "#fcfcfb"

plt.rcParams.update({
    "font.size": 9, "axes.edgecolor": LIGHT, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "axes.titlesize": 10, "axes.titleweight": "bold", "axes.titlelocation": "left",
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": False,
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
})


def gains(w, b, t):
    single = w["Col-0"].between(0.7, 1.3) & (b["Col-0"] >= 0.95)
    return {a: ((w[a] >= t) & single & (b[a] >= 0.8)) for a in ORDER if a != "Col-0"}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pav", default="results/pav/pav_long.tsv.gz")
    ap.add_argument("--lof", default="results/variants/lof_genes.tsv.gz")
    ap.add_argument("--gain", type=float, default=2.0)
    ap.add_argument("--out", default="results/figures/26_geneloss_lost_broken_gained__bar_chart_error_bars_scatter.png")
    ap.add_argument("--table", default="results/pav/gene_gains.tsv")
    a = ap.parse_args()

    p = pd.read_csv(a.pav, sep="\t", usecols=["gene", "accession", "norm_depth", "breadth_q10", "call"])
    w = p.pivot(index="gene", columns="accession", values="norm_depth")
    b = p.pivot(index="gene", columns="accession", values="breadth_q10")
    lost = p[p["call"] == "absent"].groupby("accession").size().reindex(ORDER, fill_value=0)
    lof = pd.read_csv(a.lof, sep="\t")
    broken = lof[lof["pav_call"] == "present"].groupby("accession").size().reindex(ORDER, fill_value=0)
    g_mid, g_lo, g_hi = gains(w, b, a.gain), gains(w, b, 2.5), gains(w, b, 1.8)
    gain = pd.Series({k: int(v.sum()) for k, v in g_mid.items()}).reindex(ORDER)
    lo = pd.Series({k: int(v.sum()) for k, v in g_lo.items()}).reindex(ORDER)
    hi = pd.Series({k: int(v.sum()) for k, v in g_hi.items()}).reindex(ORDER)
    spread = {acc: float(np.subtract(*np.percentile(w.loc[w["Col-0"].between(0.7, 1.3), acc], [90, 10])))
              for acc in ORDER}
    for acc in ORDER:
        gtxt = "n/a (reference)" if acc == "Col-0" else f"{int(gain[acc])} (range {int(lo[acc])}-{int(hi[acc])})"
        print(f"{acc:<16} lost {int(lost[acc]):4d} | broken {int(broken[acc]):5d} | gained {gtxt:<22} | "
              f"depth spread (90th-10th pct) {spread[acc]:.2f}")

    rows = [(gene, acc, round(float(w.loc[gene, acc]), 3), round(float(w.loc[gene, "Col-0"]), 3))
            for acc, m in g_mid.items() for gene in m[m].index]
    os.makedirs(os.path.dirname(a.table) or ".", exist_ok=True)
    pd.DataFrame(rows, columns=["gene", "accession", "norm_depth", "norm_depth_col0"]).to_csv(a.table, sep="\t",
                                                                                              index=False)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 4.8), gridspec_kw={"width_ratios": [2.2, 1], "wspace": 0.25})
    xs = np.arange(len(ORDER))
    wd = 0.27
    ax1.bar(xs - wd, lost, wd, color=ORANGE, label="lost (absent by read depth)")
    ax1.bar(xs, broken, wd, color=BLUE, label="broken (present, HIGH-impact variant)")
    gx = xs[1:] + wd
    ax1.bar(gx, gain.iloc[1:], wd, color="#2a9d8f", label=f"gained? (depth >= {a.gain:g}x; whisker 2.5x-1.8x)")
    ax1.errorbar(gx, gain.iloc[1:], yerr=[gain.iloc[1:] - lo.iloc[1:], hi.iloc[1:] - gain.iloc[1:]], fmt="none",
                 ecolor=INK2, elinewidth=1, capsize=3)
    ax1.text(xs[0] + wd, 40, "n/a", ha="center", fontsize=7.5, color=INK2)
    for x0, v in zip(xs - wd, lost):
        ax1.text(x0, v + 40, f"{int(v)}", ha="center", fontsize=7, color=INK2)
    for x0, v in zip(xs, broken):
        ax1.text(x0, v + 40, f"{int(v):,}", ha="center", fontsize=7, color=INK2)
    ax1.set_xticks(xs)
    ax1.set_xticklabels([x.replace("Ale-Stenar-44-4", "Ale-S.-44-4") for x in ORDER], rotation=25, ha="right")
    ax1.set_ylabel("number of genes")
    ax1.set_ylim(0, max(broken.max(), hi.max()) * 1.18)
    ax1.set_title("(a) Three kinds of gene change per accession")
    ax1.legend(frameon=False, loc="upper left", fontsize=8)
    ax1.grid(axis="y", color=GRID, lw=0.6)
    ax1.set_axisbelow(True)

    others = ORDER[1:]
    sp = np.array([spread[x] for x in others])
    ax2.scatter(sp, hi.iloc[1:], s=50, color="#2a9d8f", zorder=3)
    for x, s0, h in zip(others, sp, hi.iloc[1:]):
        dx, dy, ha = (-7, -12, "right") if x == "An-1" else (5, 4, "left")  # An-1 and Se-0 lie close together
        ax2.annotate(x.replace("Ale-Stenar-44-4", "Ale-S.-44-4"), (s0, h), xytext=(dx, dy), textcoords="offset points",
                     fontsize=8, color=INK2, ha=ha)
    ax2.set_xlabel("spread of normalised depth (90th - 10th percentile)")
    ax2.set_ylabel("gain candidates at 1.8x")
    ax2.set_title("(b) Gain counts follow data noise")
    ax2.set_xlim(sp.min() - 0.05, sp.max() + 0.12)
    ax2.grid(color=GRID, lw=0.6)
    ax2.set_axisbelow(True)
    fig.suptitle("Genes lost, broken or possibly gained (gains are read-depth candidates, not confirmed)",
                 x=0.01, ha="left", fontsize=11)
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    fig.savefig(a.out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("wrote", a.out, "and", a.table)


if __name__ == "__main__":
    main()
