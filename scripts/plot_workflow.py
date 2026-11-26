#!/usr/bin/env python3
"""Figure 01: workflow of the project (methods overview), env `map`.

Run:
    python scripts/plot_workflow.py
Output: results/figures/01_overview_workflow__flowchart.png (300 dpi). Numbers are typed in from the research log, not computed here.
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

INK, INK2, SURFACE = "#0b0b0b", "#52514e", "#fcfcfb"
LANE = {"reads": ("#2a78d6", "#e9f1fb"), "var": ("#eb6834", "#fdeee7"), "asm": ("#2a9d8f", "#e6f4f2"),
        "out": ("#1c4f8f", "#eef0f6")}


def box(ax, x, y, w, h, lane, title, body):
    edge, fill = LANE[lane]
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle="round,pad=0.01,rounding_size=0.015",
                                fc=fill, ec=edge, lw=1.2))
    ax.text(x, y + h / 2 - 0.018, title, ha="center", va="top", fontsize=8.6, fontweight="bold", color=INK)
    ax.text(x, y + h / 2 - 0.048, body, ha="center", va="top", fontsize=7.2, color=INK2, linespacing=1.35)


def arrow(ax, p, q, col=INK2):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=10, lw=1, color=col,
                                 shrinkA=2, shrinkB=2, connectionstyle="arc3,rad=0"))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="results/figures/01_overview_workflow__flowchart.png")
    a = ap.parse_args()

    fig = plt.figure(figsize=(12, 10.5))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    W, H = 0.27, 0.085
    X = [0.17, 0.5, 0.83]

    ax.text(0.03, 0.975, "Workflow: gene presence/absence in 8 Arabidopsis thaliana accessions", fontsize=13,
            fontweight="bold", va="top")
    ax.text(0.03, 0.945, "Three independent lines of evidence (reads, variants, assemblies) meet in a joint "
            "validation, then pan-genome analyses.", fontsize=9, color=INK2, va="top")

    # inputs
    box(ax, X[0], 0.86, W, H, "reads", "Short reads (8 accessions)",
        "Illumina WGS, 1001 Genomes (PRJNA273563)\nCol-0 (control) + 7 accessions")
    box(ax, X[1], 0.86, W, H, "out", "Reference genome",
        "TAIR10.1 (Col-0), RefSeq\n27,444 protein-coding genes")
    box(ax, X[2], 0.86, W, H, "asm", "Genome assemblies (4)",
        "An-1, Cvi-0, Kyoto, Ler\nchromosome-level, NCBI")

    # lane 1: reads
    box(ax, X[0], 0.735, W, H, "reads", "Trim, subsample, QC",
        "fastp; ~15x of trimmed bases\nFastQC / MultiQC")
    box(ax, X[0], 0.61, W, H, "reads", "Map to TAIR10.1",
        "bowtie2, mark duplicates\n94-97% mapped, 7-12x median depth")
    box(ax, X[0], 0.485, W, H, "reads", "Gene presence/absence",
        "mosdepth depth + breadth (MAPQ >= 10)\n173-452 absent genes per accession")
    # lane 2: variants
    box(ax, X[1], 0.735, W, H, "var", "Call SNPs and indels",
        "bcftools; 1.85 M sites after filter\nTs/Tv 1.21")
    box(ax, X[1], 0.61, W, H, "var", "Predict effects",
        "SnpEff; 10,215 HIGH-impact variants\n2,070-2,731 'broken' genes per accession")
    box(ax, X[1], 0.485, W, H, "var", "Call deletions",
        "Delly 2.8.2 (PASS)\n5 (Col-0) to 1,370 (Ler-1) deletions")
    # lane 3: assemblies
    box(ax, X[2], 0.735, W, H, "asm", "Assembly quality",
        "QUAST: N50 23-24 Mb, ~90% of Col-0\nBUSCO (brassicales_odb12)")
    box(ax, X[2], 0.61, W, H, "asm", "Lift genes onto assemblies",
        "Liftoff 1.6.3 (gene spans)\nfound / partial / elsewhere / missing")

    # validation + outputs
    box(ax, 0.5, 0.345, 0.6, 0.095, "out", "Three-method validation",
        "89-94% of read-absent genes also missing in the assembly; 10-47x enrichment inside Delly deletions\n"
        "only 21-41 read-absent genes per accession unsupported by both (likely errors)")
    outs = [("Pan-genome", "core 26,533 (96.7%), shell 463,\nprivate 18; accumulation curves"),
            ("Relatedness", "SNP distance, PCoA\nCvi-0 most distant"),
            ("Function", "GO enrichment (clusterProfiler)\nNLR resistance genes, FDR 3e-36"),
            ("Gene change", "lost / broken / gained per accession\nexample gene FAR4 (Cvi-0)")]
    xs = [0.14, 0.38, 0.62, 0.86]
    for x, (t, b) in zip(xs, outs):
        box(ax, x, 0.17, 0.205, 0.085, "out", t, b)
        arrow(ax, (0.5 + (x - 0.5) * 0.55, 0.297), (x, 0.214), LANE["out"][0])

    for x in X:
        arrow(ax, (x, 0.817), (x, 0.778), INK2)
    arrow(ax, (X[0], 0.692), (X[0], 0.653), LANE["reads"][0])
    arrow(ax, (X[0], 0.567), (X[0], 0.528), LANE["reads"][0])
    arrow(ax, (X[1], 0.692), (X[1], 0.653), LANE["var"][0])
    arrow(ax, (X[2], 0.692), (X[2], 0.653), LANE["asm"][0])
    # mapped reads feed variant calling and Delly
    arrow(ax, (X[0] + W / 2, 0.61), (X[1] - W / 2, 0.735), LANE["reads"][0])
    arrow(ax, (X[0] + W / 2, 0.60), (X[1] - W / 2, 0.485), LANE["reads"][0])
    # reference feeds mapping and Liftoff
    arrow(ax, (X[1] - W / 2 + 0.02, 0.817), (X[0] + W / 2 - 0.02, 0.778), INK2)
    arrow(ax, (X[1] + W / 2 - 0.02, 0.817), (X[2] - W / 2 + 0.02, 0.653), INK2)
    # into validation
    arrow(ax, (X[0], 0.442), (0.3, 0.393), LANE["reads"][0])
    arrow(ax, (X[1], 0.442), (0.5, 0.393), LANE["var"][0])
    arrow(ax, (X[2], 0.567), (0.7, 0.393), LANE["asm"][0])

    # broken genes (SnpEff) go straight to the gene-change result, not through the validation
    ax.add_patch(FancyArrowPatch((X[1] + W / 2, 0.60), (xs[3], 0.214), arrowstyle="-|>", mutation_scale=10, lw=1,
                                 color=LANE["var"][0], ls="--", connectionstyle="arc3,rad=-0.25"))
    ax.text(0.03, 0.095, "Lanes: blue = reads, orange = variants, teal = assemblies, dark blue = shared inputs and "
            "results.\nTools: fastp, bowtie2, samtools, mosdepth, bcftools, SnpEff, Delly, Liftoff, minimap2, QUAST, "
            "BUSCO, R clusterProfiler, Python.", fontsize=7.5, color=INK2, va="top", linespacing=1.5)
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    fig.savefig(a.out, dpi=300, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
