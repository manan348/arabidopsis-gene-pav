#!/usr/bin/env python3
"""Report-style figures for the work up to PAV calling (Phases 3 and 4).

Run from the project root after scripts/build_pav.py:
    python scripts/plot_project_figures.py

Inputs (read only, never modified):
    results/pav/pav_long.tsv.gz                                 from scripts/build_pav.py
    results/multiqc/multiqc_mapping_data/multiqc_general_stats.txt   MultiQC of results/mapping
    results/mapping/<S>/<S>.flagstat.txt                        read totals and duplicates
Outputs (PNG, 300 dpi) in results/figures/:
    03_mapping_quality__bar_chart_cleveland_dot_plot.png   how good is each accession's mapping (6 panels)
    05_pav_cutoff_check__hexbin_sensitivity_plot.png   Col-0 control vs the other accessions, and how stable the absent calls are
    07_pav_results__stacked_bar_ecdf_scatter.png        what the PAV calls say (4 panels)

The threshold options below are only used to DRAW the reference lines. They must equal the values
used in scripts/build_pav.py (defaults are the same). Colours follow the validated dataviz palette:
blue = accessions, dark grey = Col-0 control, orange = absent, grey = not-clear calls.
"""
import argparse
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, LogNorm

ORDER = ["Col-0", "An-1", "Kyoto", "Kondara", "Cvi-0", "Ler-1", "Se-0", "Ale-Stenar-44-4"]
N_GENES = 27444  # lines of data/reference/genes.bed

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
def read_flagstat(mapping_dir, s):
    """Return (total reads in the BAM, duplicate-flagged reads) from samtools flagstat."""
    total = dup = None
    with open(os.path.join(mapping_dir, s, f"{s}.flagstat.txt")) as fh:
        for line in fh:
            m_total = re.match(r"^(\d+) \+ \d+ in total", line)
            m_dup = re.match(r"^(\d+) \+ \d+ duplicates$", line.strip())
            if m_total:
                total = int(m_total.group(1))
            if m_dup:
                dup = int(m_dup.group(1))
    assert total is not None and dup is not None, f"could not parse flagstat of {s}"
    return total, dup


def load_qc(mqc_path, mapping_dir):
    g = pd.read_csv(mqc_path, sep="\t", index_col=0)
    need = {"samtools_stats-reads_mapped_percent": "mapped_pct",
            "samtools_stats-reads_properly_paired_percent": "proper_pct",
            "samtools_stats-reads_MQ0_percent": "mq0_pct",
            "samtools_stats-error_rate": "error_rate_pct"}
    missing = [c for c in need if c not in g.columns]
    assert not missing, f"MultiQC general stats lacks columns {missing}"
    assert set(ORDER) <= set(g.index), f"samples missing from MultiQC: {set(ORDER) - set(g.index)}"
    qc = g.loc[ORDER, list(need)].rename(columns=need).astype(float)
    tot_dup = [read_flagstat(mapping_dir, s) for s in ORDER]
    qc["pairs_m"] = [t / 2 / 1e6 for t, _ in tot_dup]  # reads / 2 = read pairs (paired-end data)
    qc["dup_pct"] = [100 * d / t for t, d in tot_dup]
    return qc


def load_pav(path):
    long = pd.read_csv(path, sep="\t")
    assert len(long) == N_GENES * len(ORDER), f"expected {N_GENES * len(ORDER)} rows, got {len(long)}"
    assert set(long["accession"]) == set(ORDER)
    assert not long[["norm_depth", "breadth_q10", "breadth_q0", "length"]].isna().any().any(), "NaN in PAV table"
    return long


def label_positions(values, min_gap):
    """Push direct labels apart (upwards) so that they do not overlap."""
    values = np.asarray(values, dtype=float)
    order = np.argsort(values)
    pos = values[order].copy()
    for i in range(1, len(pos)):
        if pos[i] - pos[i - 1] < min_gap:
            pos[i] = pos[i - 1] + min_gap
    out = np.empty_like(pos)
    out[order] = pos
    return out


def save(fig, path):
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


# ----------------------------------------------------------------------------- figure 1
def fig1(qc, medians, out):
    fig, axes = plt.subplots(2, 3, figsize=(11.5, 6.4), sharey=True)
    y = np.arange(len(ORDER))
    colours = [CTRL if a == "Col-0" else BLUE for a in ORDER]
    panels = [  # title, values, style, format
        ("(a) Read pairs in the final data (millions)", qc["pairs_m"].values, "bar", "{:.1f}"),
        ("(b) Median depth of chr1-5 (x)", medians.reindex(ORDER).values, "bar", "{:.0f}"),
        ("(c) Reads mapped (%)", qc["mapped_pct"].values, "dot", "{:.1f}"),
        ("(d) Reads properly paired (%)", qc["proper_pct"].values, "dot", "{:.1f}"),
        ("(e) Reads with MAPQ 0 (%), repeat-like", qc["mq0_pct"].values, "dot", "{:.1f}"),
        ("(f) Duplicate reads (%)", qc["dup_pct"].values, "dot", "{:.1f}"),
    ]
    for ax, (title, v, style, fmt) in zip(axes.ravel(), panels):
        if style == "bar":
            ax.barh(y, v, color=colours, height=0.62)
            ax.set_xlim(0, v.max() * 1.18)
            pad = v.max() * 0.015
        else:
            lo, hi = np.floor(v.min() - 1), np.ceil(v.max() + 1)
            ax.hlines(y, lo, v, color=LIGHT, lw=1.4)
            ax.scatter(v, y, s=46, color=colours, zorder=3)
            ax.set_xlim(lo, hi + (hi - lo) * 0.14)
            pad = (hi - lo) * 0.025
        for yi, vi in zip(y, v):
            ax.text(vi + pad, yi, fmt.format(vi), va="center", ha="left", fontsize=8, color=INK2)
        ax.set_title(title)
        ax.set_yticks(y)
        ax.set_yticklabels(ORDER)
        ax.grid(axis="x", color=GRID, lw=0.6)
        ax.set_axisbelow(True)
    axes[0, 0].invert_yaxis()
    fig.suptitle("Mapping of the 8 accessions to the Col-0 reference (dark grey = Col-0 control; dot axes do not start at 0)",
                 x=0.01, ha="left", fontsize=11, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    save(fig, out)


# ----------------------------------------------------------------------------- figure 2
def hex_panel(ax, d, weight, norm):
    x = d["breadth_q10"].to_numpy()
    yv = d["norm_depth"].clip(upper=3).to_numpy()
    return ax.hexbin(x, yv, C=np.full(len(x), weight), reduce_C_function=np.sum, gridsize=(36, 24),
                     extent=(0, 1, 0, 3), mincnt=1, cmap=SEQ, norm=norm, linewidths=0.2, edgecolors=SURFACE)


def fig2(long, a, out):
    fig = plt.figure(figsize=(14, 4.6))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.25], wspace=0.28)
    ax1, ax2, ax3 = (fig.add_subplot(gs[0, i]) for i in range(3))
    norm = LogNorm(vmin=0.1, vmax=1e4)
    col = long[long["accession"] == "Col-0"]
    oth = long[long["accession"] != "Col-0"]
    n_oth = oth["accession"].nunique()
    hex_panel(ax1, col, 1.0, norm)
    hb = hex_panel(ax2, oth, 1.0 / n_oth, norm)
    for ax, title in ((ax1, "(a) Col-0 (control)"), (ax2, f"(b) {n_oth} other accessions (average per accession)")):
        ax.axvline(a.absent_breadth, color=INK2, ls="--", lw=0.9)
        ax.axvline(a.present_breadth, color=INK2, ls="--", lw=0.9)
        ax.axhline(a.present_depth, color=INK2, ls="--", lw=0.9)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 3)
        ax.set_title(title)
        ax.set_xlabel("breadth: fraction of the gene covered (MAPQ >= 10)")
    ax1.set_ylabel("gene depth / sample median (capped at 3)")
    ax2.text(a.absent_breadth - 0.01, 2.9, "absent", ha="right", va="top", fontsize=8, color=INK2)
    ax2.text(a.present_breadth + 0.01, 2.9, "present", ha="left", va="top", fontsize=8, color=INK2)
    ax2.text(0.5, a.present_depth + 0.04, f"{a.present_depth:g} x median", ha="center", va="bottom", fontsize=8, color=INK2)
    cb = fig.colorbar(hb, ax=[ax1, ax2], fraction=0.035, pad=0.02)
    cb.set_label("genes per cell (log scale)")

    cuts = np.linspace(0.02, 0.30, 15)
    ymax = 0
    ends = {}
    for s in ORDER:
        d = long[long["accession"] == s]
        n_abs = [int(((d["breadth_q10"] < t) & (d["breadth_q0"] < a.lowmap_breadth)).sum()) for t in cuts]
        colour = CTRL if s == "Col-0" else BLUE
        ax3.plot(cuts, n_abs, color=colour, lw=1.6, alpha=1.0 if s == "Col-0" else 0.75)
        ends[s] = n_abs[-1]
        ymax = max(ymax, max(n_abs))
    pos = label_positions([ends[s] for s in ORDER], min_gap=ymax * 0.055)
    for s, p in zip(ORDER, pos):
        ax3.text(cuts[-1] + 0.006, p, s, va="center", fontsize=8, color=INK2)
    ax3.axvline(a.absent_breadth, color=INK2, ls="--", lw=0.9)
    ax3.text(a.absent_breadth + 0.004, ymax * 0.02, "cut-off used", fontsize=8, color=INK2, rotation=90, va="bottom")
    ax3.set_xlim(cuts[0], cuts[-1] + 0.085)
    ax3.set_xlabel("absent if the covered fraction is below this value")
    ax3.set_ylabel("genes called absent")
    ax3.set_title("(c) Absent calls barely change with the cut-off")
    ax3.grid(axis="y", color=GRID, lw=0.6)
    ax3.set_axisbelow(True)
    save(fig, out)


# ----------------------------------------------------------------------------- figure 3
def ecdf(ax, values, colour, label):
    v = np.sort(values)
    ax.step(v, np.arange(1, len(v) + 1) / len(v), where="post", color=colour, lw=1.8, label=label)


def fig3(long, qc, a, out):
    others = [s for s in ORDER if s != "Col-0"]
    calls = long.pivot(index="gene", columns="accession", values="call")[ORDER]
    assert len(calls) == N_GENES
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.2))
    ax1, ax2, ax3, ax4 = axes.ravel()

    # (a) all genes that are not a clear "present", per accession
    cnt = long.groupby(["accession", "call"]).size().unstack(fill_value=0).reindex(ORDER)
    y = np.arange(len(ORDER))
    left = np.zeros(len(ORDER))
    for name, colour in (("absent", ORANGE), ("uncertain", BLUE), ("ambiguous", GREY)):
        ax1.barh(y, cnt[name].values, left=left, color=colour, height=0.62, label=name, edgecolor=SURFACE, linewidth=1)
        left += cnt[name].values
    for yi, tot in zip(y, left):
        ax1.text(tot + 12, yi, f"{int(tot)}", va="center", fontsize=8, color=INK2)
    ax1.set_yticks(y)
    ax1.set_yticklabels(ORDER)
    ax1.invert_yaxis()
    ax1.set_xlim(0, left.max() * 1.12)
    ax1.set_xlabel(f"genes (of {N_GENES:,}) not called present")
    ax1.set_title("(a) Calls other than present")
    ax1.legend(frameon=False, loc="lower right")

    # (b) in how many of the 7 other accessions is each gene absent
    k = (calls[others] == "absent").sum(axis=1)
    freq = k.value_counts().reindex(range(1, len(others) + 1), fill_value=0)
    core = int((calls == "present").all(axis=1).sum())
    bar_colours = [BLUE] * (len(others) - 1) + [ORANGE]
    ax2.bar(freq.index, freq.values, color=bar_colours, width=0.68)
    for xi, vi in zip(freq.index, freq.values):
        ax2.text(xi, vi + freq.max() * 0.015, f"{vi}", ha="center", fontsize=8, color=INK2)
    ax2.set_xticks(freq.index)
    ax2.set_xlabel(f"number of the {len(others)} non-Col-0 accessions in which the gene is absent")
    ax2.set_ylabel("genes")
    ax2.set_ylim(0, freq.max() * 1.15)
    ax2.set_title("(b) How widely each gene is lost")
    ax2.text(0.97, 0.93, f"clean core (present in all 8): {core:,}\norange = absent in every other accession\n(present in Col-0 only)",
             transform=ax2.transAxes, ha="right", va="top", fontsize=8, color=INK2)

    # (c) length of absent genes vs all genes
    length = long[long["accession"] == "Col-0"].set_index("gene")["length"].reindex(calls.index)
    any_abs = (calls == "absent").any(axis=1)
    ecdf(ax3, length.values, GREY, f"all genes (n={len(length):,}), median {int(length.median())} bp")
    ecdf(ax3, length[any_abs].values, ORANGE, f"absent in >= 1 accession (n={int(any_abs.sum())}), median {int(length[any_abs].median())} bp")
    ax3.set_xscale("log")
    ax3.set_xlabel("gene length (bp, log scale)")
    ax3.set_ylabel("fraction of genes (cumulative)")
    ax3.set_title("(c) Absent genes are shorter or longer?")
    ax3.legend(frameon=False, loc="lower right", fontsize=8)
    ax3.grid(color=GRID, lw=0.6)
    ax3.set_axisbelow(True)

    # (d) number absent vs divergence from Col-0
    xs = qc["error_rate_pct"].reindex(ORDER).values
    ys = cnt["absent"].reindex(ORDER).values
    cols = [CTRL if s == "Col-0" else BLUE for s in ORDER]
    ax4.scatter(xs, ys, s=60, color=cols, zorder=3, edgecolor=SURFACE, linewidth=1.2)
    pos = label_positions(ys, min_gap=max(ys) * 0.06)
    for s, xv, yv in zip(ORDER, xs, pos):
        ax4.text(xv + 0.02, yv, s, va="center", fontsize=8, color=INK2)
    rho = pd.Series(xs[1:]).corr(pd.Series(ys[1:]), method="spearman")  # 7 non-Col-0 accessions
    ax4.set_xlim(xs.min() - 0.05, xs.max() + 0.42)
    ax4.set_xlabel("mismatch rate of mapped reads vs Col-0 (samtools error rate, %)")
    ax4.set_ylabel("genes called absent")
    ax4.set_title("(d) Absent genes vs divergence from Col-0")
    ax4.text(0.03, 0.95, f"Spearman rho = {rho:.2f} (n = {len(others)}, descriptive only,\nnot a test; Col-0 excluded)",
             transform=ax4.transAxes, ha="left", va="top", fontsize=8, color=INK2)
    ax4.grid(color=GRID, lw=0.6)
    ax4.set_axisbelow(True)

    fig.suptitle("Gene presence/absence relative to the Col-0 reference (27,444 genes, 8 accessions)", x=0.01, ha="left", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    save(fig, out)


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pav", default="results/pav/pav_long.tsv.gz")
    ap.add_argument("--mqc", default="results/multiqc/multiqc_mapping_data/multiqc_general_stats.txt")
    ap.add_argument("--mapdir", default="results/mapping")
    ap.add_argument("--outdir", default="results/figures")
    # reference lines; must equal the thresholds of scripts/build_pav.py
    ap.add_argument("--present-breadth", type=float, default=0.8)
    ap.add_argument("--present-depth", type=float, default=0.2)
    ap.add_argument("--absent-breadth", type=float, default=0.1)
    ap.add_argument("--lowmap-breadth", type=float, default=0.5)
    a = ap.parse_args()

    long = load_pav(a.pav)
    qc = load_qc(a.mqc, a.mapdir)
    medians = long.groupby("accession")["median_depth"].first()
    os.makedirs(a.outdir, exist_ok=True)
    fig1(qc, medians, os.path.join(a.outdir, "03_mapping_quality__bar_chart_cleveland_dot_plot.png"))
    fig2(long, a, os.path.join(a.outdir, "05_pav_cutoff_check__hexbin_sensitivity_plot.png"))
    fig3(long, qc, a, os.path.join(a.outdir, "07_pav_results__stacked_bar_ecdf_scatter.png"))
    print("wrote 3 figures to", a.outdir)


if __name__ == "__main__":
    main()
