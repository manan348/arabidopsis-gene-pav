#!/usr/bin/env python3
"""Figure 17: are the 4 assemblies complete and contiguous enough to judge genes? BUSCO + QUAST, env `map`.

Run (after scripts/run_busco.sh for all 5 genomes and scripts/run_quast.sh):
    python scripts/plot_assembly_quality.py

Inputs (read only): results/assembly/busco/<genome>/short_summary*.txt, results/assembly/quast/report.tsv
Output: results/figures/17_assembly_quality__busco_stacked_bar_cleveland_dot_plot.png (300 dpi, with caption band)

(a) BUSCO 6.1.0, lineage brassicales_odb12 (4,311 conserved single-copy genes), genome mode (miniprot).
    Col-0 = the TAIR10.1 reference itself, the baseline.
(b) QUAST 5.3.0 against Col-0: total length, N50, share of Col-0 covered by alignments, mismatches per 100 kb.
    Much of the "misassembly" and mismatch count is real difference between accessions, not assembly error.
Descriptive only: no statistical test.
"""
import argparse
import glob
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont

GENOMES = ["Col-0", "An-1", "Cvi-0", "Kyoto", "Ler-1"]
ASM = GENOMES[1:]
CATS = [("S", "complete, single copy", "#1c4f8f"), ("D", "complete, duplicated", "#7fb3ea"),
        ("F", "fragmented", "#e3b23c"), ("M", "missing", "#eb6834")]
COL = {"An-1": "#2a78d6", "Ler-1": "#1c4f8f", "Cvi-0": "#eb6834", "Kyoto": "#2a9d8f"}
LIGHT, INK, INK2, GRID, SURFACE = "#c3c2b7", "#0b0b0b", "#52514e", "#e1e0d9", "#fcfcfb"
plt.rcParams.update({
    "font.size": 9, "axes.edgecolor": LIGHT, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "axes.titlesize": 10, "axes.titleweight": "bold", "axes.titlelocation": "left",
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
})


def busco_counts(path):
    txt = open(path).read()
    want = {"S": r"(\d+)\s+Complete and single-copy BUSCOs", "D": r"(\d+)\s+Complete and duplicated BUSCOs",
            "F": r"(\d+)\s+Fragmented BUSCOs", "M": r"(\d+)\s+Missing BUSCOs", "n": r"(\d+)\s+Total BUSCO groups"}
    return {k: int(re.search(p, txt).group(1)) for k, p in want.items()}


def add_caption(path, number, stage, types):
    im = Image.open(path).convert("RGB")
    w, h = im.size
    size = max(28, w // 75)
    fdir = os.path.join(os.path.dirname(matplotlib.__file__), "mpl-data", "fonts", "ttf")
    fr = ImageFont.truetype(os.path.join(fdir, "DejaVuSans.ttf"), size)
    fb = ImageFont.truetype(os.path.join(fdir, "DejaVuSans-Bold.ttf"), size)
    band = int(size * 2.4)
    out = Image.new("RGB", (w, h + band), SURFACE)
    out.paste(im, (0, 0))
    d = ImageDraw.Draw(out)
    d.rectangle([0, h, w, h + band], fill="#eef0f6")
    d.line([0, h, w, h], fill=LIGHT, width=2)
    left = f"Figure {number}  |  Stage {stage}  |  "
    y = h + (band - size) // 2 - 2
    d.text((int(size * 0.8), y), left, font=fr, fill=INK2)
    d.text((int(size * 0.8) + d.textlength(left, font=fr), y), "Plot type: " + types, font=fb, fill="#1c4f8f")
    out.save(path, dpi=(300, 300))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--busco", default="results/assembly/busco")
    ap.add_argument("--quast", default="results/assembly/quast/report.tsv")
    ap.add_argument("--out", default="results/figures/17_assembly_quality__busco_stacked_bar_cleveland_dot_plot.png")
    a = ap.parse_args()

    b = {}
    for g in GENOMES:
        f = glob.glob(f"{a.busco}/{g}/short_summary*.txt")
        assert len(f) == 1, f"BUSCO summary for {g}: {f}"
        b[g] = busco_counts(f[0])
    n = b["Col-0"]["n"]
    assert all(v["n"] == n for v in b.values())
    q = pd.read_csv(a.quast, sep="\t", index_col=0)[ASM]
    ref_len = float(q.loc["Reference length", ASM[0]])
    print(f"BUSCO (n = {n}):")
    for g in GENOMES:
        v = b[g]
        print(f"  {g:<6} C {100 * (v['S'] + v['D']) / n:5.1f}%  S {v['S']}  D {v['D']}  F {v['F']}  M {v['M']}")

    fig = plt.figure(figsize=(13, 5.6))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.5, 1, 1], hspace=0.7, wspace=0.42)
    ax = fig.add_subplot(gs[:, 0])
    ys = np.arange(len(GENOMES))[::-1]
    for y, g in zip(ys, GENOMES):
        left = 0.0
        for k, lab, col in CATS:
            pct = 100 * b[g][k] / n
            ax.barh(y, pct, left=left, color=col, height=0.62, label=lab if g == "Col-0" else None,
                    edgecolor=SURFACE, linewidth=0.6)
            left += pct
        c = 100 * (b[g]["S"] + b[g]["D"]) / n
        ax.text(101, y, f"C {c:.1f}%  (M {b[g]['M']})", va="center", fontsize=8.5, color=INK)
    ax.set_yticks(ys)
    ax.set_yticklabels([g + (" (reference)" if g == "Col-0" else "") for g in GENOMES])
    ax.set_xlim(0, 125)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xlabel(f"% of {n:,} conserved Brassicales genes (BUSCO)")
    ax.set_title("(a) Gene completeness (BUSCO)")
    ax.legend(frameon=False, fontsize=8, loc="upper center", bbox_to_anchor=(0.45, -0.13), ncol=2)
    ax.tick_params(axis="y", length=0)

    metrics = [("Total length", 1e6, "total length (Mb)", ref_len / 1e6, "Col-0 reference"),
               ("N50", 1e6, "N50 (Mb)", None, None),
               ("Genome fraction (%)", 1, "% of Col-0 covered", None, None),
               ("# mismatches per 100 kbp", 1, "mismatches per 100 kb", None, None)]
    titles = ["(b) Size", "(c) Contiguity", "(d) Col-0 covered", "(e) Difference from Col-0"]
    for i, ((row, scale, xl, ref, rlab), t) in enumerate(zip(metrics, titles)):
        axm = fig.add_subplot(gs[i // 2, 1 + i % 2])
        vals = q.loc[row].astype(float) / scale
        yy = np.arange(len(ASM))[::-1]
        axm.hlines(yy, vals.min() - (vals.max() - vals.min()) * 0.6 - 1e-9, vals, color=GRID, lw=1.2)
        axm.scatter(vals, yy, s=55, color=[COL[x] for x in ASM], zorder=3, edgecolor=SURFACE, lw=1)
        for v, y in zip(vals, yy):
            axm.text(v, y + 0.3, f"{v:,.1f}" if v < 1000 else f"{v:,.0f}", ha="center", fontsize=7.5, color=INK2)
        if ref:
            axm.axvline(ref, color=INK2, ls=":", lw=1)
            axm.text(ref, -0.5, " " + rlab, fontsize=7, color=INK2, ha="left", va="bottom")
        axm.set_yticks(yy)
        axm.set_yticklabels(ASM, fontsize=8)
        axm.set_ylim(-0.6, len(ASM) - 0.1)
        span = vals.max() - vals.min()
        lo = min(vals.min(), ref if ref else vals.min()) - span * 0.6 - 0.05
        hi = max(vals.max(), ref if ref else vals.max()) + span * 0.6 + 0.05
        axm.set_xlim(lo, hi)
        axm.set_xlabel(xl, fontsize=8.5)
        axm.set_title(t, fontsize=9.5)
        axm.tick_params(axis="y", length=0)
        axm.grid(axis="x", color=GRID, lw=0.5)
        axm.set_axisbelow(True)
    fig.suptitle("Assembly quality: the 4 assemblies are complete and nearly chromosome-scale, so a gene missing "
                 "from them is really missing", x=0.01, ha="left", fontsize=11.5, fontweight="bold")
    fig.text(0.01, 0.925, f"BUSCO 6.1.0 (brassicales_odb12, {n:,} genes) and QUAST 5.3.0 against Col-0. Dot-plot axes "
             "do not start at 0. Mismatches include real differences between accessions.", fontsize=8.5,
             color=INK2)
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    fig.savefig(a.out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    add_caption(a.out, "17", "6. Check with assemblies", "100% stacked bar chart (BUSCO) + Cleveland dot plots (QUAST)")
    print("wrote", a.out)


if __name__ == "__main__":
    main()
