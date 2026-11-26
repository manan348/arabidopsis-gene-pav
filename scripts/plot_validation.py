#!/usr/bin/env python3
"""Figures 12 and 13: read coverage vs assembly, and agreement of three methods (Phase 6), run in env `map`.

Run after scripts/compare_liftoff_pav.py:
    python scripts/plot_validation.py

Inputs (read only): results/pav/pav_long.tsv.gz (read breadth and call), results/assembly/liftoff_vs_pav_genes.tsv.gz
                    (Liftoff class and coverage), results/sv/sv_vs_pav_genes.tsv.gz (share of the gene inside a Delly deletion)
Outputs: results/figures/19_assembly_reads_vs_liftoff__scatter_plots.png, results/figures/20_validation_three_methods__stacked_grouped_bar_chart.png (300 dpi)

fig 19: for every gene in the 4 assembled accessions, read breadth (share of the gene covered by reads with
       MAPQ >= 10) against Liftoff coverage in the assembly (genes Liftoff could not place are drawn in a band
       at the bottom). Colour = read-based call.
fig 20: (a) the genes the reads call absent, split by whether the assembly (not found intact) and Delly
       (>= 50% of the gene inside a PASS deletion) agree; (b) Delly against the assembly: share of genes inside
       a Delly deletion among genes the assembly lacks (not placed) and among genes it has (found intact).
Descriptive only: no statistical test.
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ACCS = ["An-1", "Cvi-0", "Kyoto", "Ler-1"]
BLUE, ORANGE, GREY, LIGHT = "#2a78d6", "#eb6834", "#898781", "#c3c2b7"
INK, INK2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e1e0d9", "#fcfcfb"
BOX = dict(facecolor="#fcfcfb", alpha=0.9, edgecolor="none", pad=1.5)  # keeps corner numbers readable
CALL_COL = {"present": LIGHT, "ambiguous": "#8fb8ea", "uncertain": BLUE, "absent": ORANGE}

plt.rcParams.update({
    "font.size": 9, "axes.edgecolor": LIGHT, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "axes.titlesize": 10, "axes.titleweight": "bold", "axes.titlelocation": "left",
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": False,
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
})


def fig12(d, out):
    rng = np.random.default_rng(1)
    fig, axes = plt.subplots(1, 4, figsize=(15, 4.4), sharey=True, gridspec_kw={"wspace": 0.08})
    for ax, acc in zip(axes, ACCS):
        x = d[d["accession"] == acc].copy()
        placed = x["coverage"].notna()
        y = np.array(x["coverage"], dtype=float)  # a writable copy (pandas 3 can return read-only arrays)
        y[~placed.to_numpy()] = -0.12 + rng.uniform(-0.03, 0.03, (~placed).sum())  # "not placed" band
        xs = x["breadth_q10"].to_numpy() + rng.uniform(-0.006, 0.006, len(x))
        for call in ["present", "ambiguous", "uncertain", "absent"]:  # absent drawn last, on top
            m = (x["call"] == call).to_numpy()
            ax.scatter(xs[m], y[m], s=4 if call == "present" else 7, color=CALL_COL[call], alpha=0.6, lw=0,
                       rasterized=True, label=call)
        lowr = x["breadth_q10"] < 0.1
        lowa = (~placed) | (x["coverage"] < 0.5)
        ax.text(0.02, 0.97, f"{int((lowr & ~lowa).sum())}", transform=ax.transAxes, va="top", fontsize=8, color=INK2, bbox=BOX)
        ax.text(0.98, 0.97, f"{int((~lowr & ~lowa).sum()):,}", transform=ax.transAxes, va="top", ha="right",
                fontsize=8, color=INK2, bbox=BOX)
        ax.text(0.02, 0.02, f"{int((lowr & lowa).sum())}", transform=ax.transAxes, fontsize=8, color=ORANGE,
                fontweight="bold", bbox=BOX)
        ax.text(0.98, 0.02, f"{int((~lowr & lowa).sum())}", transform=ax.transAxes, ha="right", fontsize=8,
                color=INK2, bbox=BOX)
        ax.axvline(0.1, color=GREY, lw=0.8, ls="--")
        ax.axhline(0.5, color=GREY, lw=0.8, ls="--")
        ax.set_xlim(-0.03, 1.03)
        ax.set_ylim(-0.2, 1.06)
        ax.set_title(acc)
        ax.set_xlabel("read breadth (MAPQ >= 10)")
        ax.grid(color=GRID, lw=0.5)
        ax.set_axisbelow(True)
    axes[0].set_yticks([-0.12, 0, 0.25, 0.5, 0.75, 1.0])
    axes[0].set_yticklabels(["not placed", "0", "0.25", "0.5", "0.75", "1"])
    axes[0].set_ylabel("Liftoff coverage in the assembly")
    h, l = axes[0].get_legend_handles_labels()
    order = [l.index(c) for c in ["absent", "uncertain", "ambiguous", "present"]]
    fig.legend([h[i] for i in order], [f"read call: {l[i]}" for i in order], frameon=False, ncol=4,
               loc="lower center", bbox_to_anchor=(0.5, -0.06), markerscale=3, fontsize=8)
    fig.suptitle("Every gene, reads vs assembly (dashed lines: read absent cut-off 0.1, assembly cut-off 0.5; "
                 "corner numbers = genes in each quarter)", x=0.01, ha="left", fontsize=10.5)
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out)


def fig13(d, out):
    ra = d[d["call"] == "absent"].copy()
    ra["asm"] = ra["asm_class"] != "found"
    ra["delly"] = ra["deleted_frac"] >= 0.5
    cats = [("both", "assembly and Delly", "#1c4f8f"), ("asm", "assembly only", BLUE),
            ("delly", "Delly only", "#8fb8ea"), ("none", "neither (likely read error)", ORANGE)]
    ra["cat"] = np.select([ra["asm"] & ra["delly"], ra["asm"], ra["delly"]], ["both", "asm", "delly"], "none")
    t = ra.groupby(["accession", "cat"]).size().unstack(fill_value=0).reindex(index=ACCS, columns=[c[0] for c in cats],
                                                                              fill_value=0)
    print("read-absent genes by support (assembly not found intact / Delly >= 50%):")
    print(t.to_string())

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 4.6), gridspec_kw={"width_ratios": [1.15, 1], "wspace": 0.3})
    ys = np.arange(len(ACCS))
    tot = t.sum(axis=1).to_numpy()
    left = np.zeros(len(ACCS))
    for key, lab, col in cats:
        pct = 100 * t[key].to_numpy() / tot
        ax1.barh(ys, pct, left=left, color=col, height=0.62, label=lab)
        for yi, l0, p, n in zip(ys, left, pct, t[key]):
            if p >= 5:
                ax1.text(l0 + p / 2, yi, f"{n}", ha="center", va="center", fontsize=8,
                         color=INK if col == "#8fb8ea" else "white")
        left += pct
    ax1.set_yticks(ys)
    ax1.set_yticklabels([f"{a}\n(n={n})" for a, n in zip(ACCS, tot)])
    ax1.invert_yaxis()
    ax1.set_xlim(0, 100)
    ax1.set_xlabel("% of genes the reads call absent")
    ax1.set_title("(a) Which other method supports the absent calls?")
    ax1.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2, fontsize=8)

    rows = []
    for acc in ACCS:
        x = d[d["accession"] == acc]
        miss, found = x[x["asm_class"] == "missing"], x[x["asm_class"] == "found"]
        rows.append((acc, 100 * (miss["deleted_frac"] >= 0.5).mean(), 100 * (found["deleted_frac"] >= 0.5).mean()))
    r = pd.DataFrame(rows, columns=["accession", "missing", "found"]).set_index("accession")
    print("\n% of genes inside a Delly deletion, by assembly class:")
    print(r.round(2).to_string())
    w = 0.38
    ax2.bar(ys - w / 2, r["missing"], w, color="#1c4f8f", label="genes the assembly lacks (not placed)")
    ax2.bar(ys + w / 2, r["found"], w, color=LIGHT, label="genes the assembly has (found intact)")
    for xi, v in zip(ys - w / 2, r["missing"]):
        ax2.text(xi, v + 0.6, f"{v:.1f}%", ha="center", fontsize=8, color=INK2)
    for xi, v in zip(ys + w / 2, r["found"]):
        ax2.text(xi, v + 0.6, f"{v:.2f}%", ha="center", fontsize=8, color=INK2)
    ax2.set_xticks(ys)
    ax2.set_xticklabels(ACCS)
    ax2.set_ylim(0, r["missing"].max() * 1.3)
    ax2.set_ylabel("% of genes inside a Delly deletion")
    ax2.set_title("(b) Delly against the assembly")
    ax2.legend(frameon=False, loc="upper left", fontsize=8)
    ax2.grid(axis="y", color=GRID, lw=0.6)
    ax2.set_axisbelow(True)
    fig.suptitle("Three independent views of gene loss: read depth, Delly deletions, PacBio assemblies",
                 x=0.01, ha="left", fontsize=11)
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pav", default="results/pav/pav_long.tsv.gz")
    ap.add_argument("--liftoff", default="results/assembly/liftoff_vs_pav_genes.tsv.gz")
    ap.add_argument("--sv", default="results/sv/sv_vs_pav_genes.tsv.gz")
    ap.add_argument("--out12", default="results/figures/19_assembly_reads_vs_liftoff__scatter_plots.png")
    ap.add_argument("--out13", default="results/figures/20_validation_three_methods__stacked_grouped_bar_chart.png")
    a = ap.parse_args()

    pav = pd.read_csv(a.pav, sep="\t", usecols=["gene", "accession", "breadth_q10"])
    lo = pd.read_csv(a.liftoff, sep="\t", usecols=["gene", "accession", "call", "asm_class", "coverage"])
    sv = pd.read_csv(a.sv, sep="\t", usecols=["gene", "accession", "deleted_frac"])
    d = lo.merge(pav, on=["gene", "accession"], how="left").merge(sv, on=["gene", "accession"], how="left")
    assert d["breadth_q10"].notna().all() and d["deleted_frac"].notna().all(), "a gene is missing from a table"
    assert len(d) == len(lo)
    os.makedirs(os.path.dirname(a.out12) or ".", exist_ok=True)
    fig12(d, a.out12)
    fig13(d, a.out13)


if __name__ == "__main__":
    main()
