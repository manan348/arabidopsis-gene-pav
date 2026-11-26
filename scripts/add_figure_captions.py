"""Add the caption band (figure number, stage, plot type) to the bottom of each figure in results/figures/.

Run after re-making a figure, env `map`:  python scripts/add_figure_captions.py
Files that already have a band are skipped, so it is safe to run more than once.
"""
import os

import matplotlib
import numpy as np
from PIL import Image, ImageDraw, ImageFont

D = "results/figures"
FD = os.path.join(os.path.dirname(matplotlib.__file__), "mpl-data", "fonts", "ttf")
BAND = (0xee, 0xf0, 0xf6)
# (file name without the number, stage, plot type); the number is the position in this list
FIGURES = [
    ("overview_workflow__flowchart", "1. Overview", "Workflow diagram (flowchart)"),
    ("overview_genome__circos_plot", "1. Overview", "Circos plot (circular genome plot)"),
    ("mapping_quality__bar_chart_cleveland_dot_plot", "2. Read mapping", "Horizontal bar charts + Cleveland dot plots"),
    ("pav_gene_depth__histograms", "3. PAV calls", "Histograms (small multiples)"),
    ("pav_cutoff_check__hexbin_sensitivity_plot", "3. PAV calls", "Hexbin plots (2-D density) + sensitivity-analysis line plot"),
    ("pav_calls_per_accession__grouped_bar_chart", "3. PAV calls", "Grouped bar chart"),
    ("pav_results__stacked_bar_ecdf_scatter", "3. PAV calls", "Stacked bar chart + frequency bar chart + ECDF + scatter plot"),
    ("pav_variable_genes__binary_heatmap", "3. PAV calls", "Binary presence/absence heatmap"),
    ("pav_clustered__clustered_heatmap_dendrogram", "3. PAV calls", "Clustered heatmap with dendrograms"),
    ("pav_relative_depth__diverging_heatmap", "3. PAV calls", "Heatmap (diverging colour scale)"),
    ("pav_shared_absences__jaccard_similarity_matrix", "3. PAV calls", "Similarity matrix heatmap (Jaccard index)"),
    ("pav_along_chromosomes__density_heatmap", "3. PAV calls", "Density heatmap along chromosomes"),
    ("variants__bar_chart_pca_distance_heatmap", "4. SNPs and diversity", "Grouped bar chart + PCA plot + distance heatmap + SNP density heatmap"),
    ("variants_snp_density__ecdf_dumbbell_plot", "4. SNPs and diversity", "ECDF + dumbbell plot"),
    ("variants_diversity__sliding_window_plot_violin_plot", "4. SNPs and diversity", "Sliding-window plot + violin plot + Mann-Whitney test"),
    ("deletions_vs_absent__dumbbell_scatter_ecdf", "5. Check with deletions", "Dumbbell plot + scatter plot + ECDF of deletion lengths"),
    ("assembly_quality__busco_stacked_bar_cleveland_dot_plot", "6. Check with assemblies", "100% stacked bar chart (BUSCO) + Cleveland dot plots (QUAST)"),
    ("assembly_liftoff__100pct_stacked_bar_chart", "6. Check with assemblies", "100% stacked bar chart"),
    ("assembly_reads_vs_liftoff__scatter_plots", "6. Check with assemblies", "Scatter plots (small multiples)"),
    ("validation_three_methods__stacked_grouped_bar_chart", "7. Validation", "Stacked bar chart + grouped bar chart"),
    ("validation_length_bias__line_plot_confidence_intervals", "7. Validation", "Line plot + dot plot with 95% confidence intervals (Wilson)"),
    ("validation_absent_genes__chromosome_ideogram", "7. Validation", "Chromosome ideogram (karyotype plot)"),
    ("validation_FAR4_example__coverage_track", "7. Validation", "Coverage track plot (genome-browser style)"),
    ("validation_rearrangements__synteny_riparian_plot", "7. Validation", "Synteny plot (plotsr, riparian plot)"),
    ("geneloss_lost_vs_broken__grouped_stacked_bar_chart", "8. Gene-loss types", "Grouped bar chart + 100% stacked bar chart"),
    ("geneloss_lost_broken_gained__bar_chart_error_bars_scatter", "8. Gene-loss types", "Grouped bar chart with error bars + scatter plot"),
    ("pangenome_categories__donut_chart_histogram", "9. Pan-genome", "Donut chart + gene frequency histogram"),
    ("pangenome_shared_absences__upset_plot", "9. Pan-genome", "UpSet plot"),
    ("pangenome_core_genes__accumulation_rarefaction_curve", "9. Pan-genome", "Gene accumulation (rarefaction) curve + gene frequency spectrum"),
    ("pangenome_gene_patterns__circular_dendrogram_heatmap", "9. Pan-genome", "Circular dendrogram with heatmap rings (circos-style)"),
    ("relatedness_snps__pcoa_plot", "10. Relatedness", "PCoA plot (principal coordinates analysis)"),
    ("relatedness_origins__geographic_map_upgma_dendrogram", "10. Relatedness", "Geographic map + UPGMA dendrogram"),
    ("relatedness_pav_vs_snp__tanglegram_mantel_test", "10. Relatedness", "Tanglegram + Mantel test"),
    ("function_go_enrichment__dot_plot", "11. Function", "GO enrichment dot plot (bubble plot)"),
    ("function_go_genes__cnetplot_network", "11. Function", "Cnetplot (gene-concept network)"),
    ("function_go_genes__chord_diagram", "11. Function", "GO chord diagram (GOChord)"),
    ("function_nlr_genes__chromosome_ideogram_fisher_test", "11. Function", "Chromosome ideogram + Fisher's exact test"),
    ("function_nlr_selection__pnps_dot_plot_ci_scatter_plot", "11. Function", "Dot plot with 95% CI + scatter plot (pN/pS, permutation test)"),
    ("function_chr5_cluster__coverage_track", "11. Function", "Coverage track plot (genome-browser style)"),
]


def has_band(im):
    a = np.asarray(im.convert("RGB"))[-5:, :15].astype(int)
    return bool(np.all(np.abs(a - BAND).max(axis=2) <= 2))


def caption(im, num, stage, types):
    w, h = im.size
    size = max(22, min(w // 75, 60))
    fr = ImageFont.truetype(os.path.join(FD, "DejaVuSans.ttf"), size)
    fb = ImageFont.truetype(os.path.join(FD, "DejaVuSans-Bold.ttf"), size)
    d0 = ImageDraw.Draw(im)
    left = f"Figure {num}  |  Stage {stage}  |  "
    full = d0.textlength(left, font=fr) + d0.textlength("Plot type: " + types, font=fb) + size * 1.6
    two = full > w
    band = int(size * (3.6 if two else 2.4))
    out = Image.new("RGB", (w, h + band), "#fcfcfb")
    out.paste(im, (0, 0))
    d = ImageDraw.Draw(out)
    d.rectangle([0, h, w, h + band], fill="#eef0f6")
    d.line([0, h, w, h], fill="#c3c2b7", width=2)
    pad = int(size * 0.8)
    if two:
        d.text((pad, h + int(size * 0.5)), left.rstrip(" |"), font=fr, fill="#52514e")
        d.text((pad, h + int(size * 1.9)), "Plot type: " + types, font=fb, fill="#1c4f8f")
    else:
        y = h + (band - size) // 2 - 2
        d.text((pad, y), left, font=fr, fill="#52514e")
        d.text((pad + d.textlength(left, font=fr), y), "Plot type: " + types, font=fb, fill="#1c4f8f")
    return out


for i, (stem, stage, types) in enumerate(FIGURES, 1):
    path = f"{D}/{i:02d}_{stem}.png"
    if not os.path.exists(path):
        print("missing", path)
        continue
    im = Image.open(path).convert("RGB")
    if has_band(im):
        continue
    caption(im, f"{i:02d}", stage, types).save(path, dpi=(300, 300))
    print("captioned", path)
