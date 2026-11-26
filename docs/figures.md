# Figure guide

All figures are in `results/figures/`. They are numbered in **pipeline order**: each stage uses the output of the one
before, and each figure answers a question raised by the previous one. In each file name, the word after the number
is the stage and the words after the double underscore (`__`) are the plot types. Each image also has a caption band
at the bottom with its number, stage and plot type.

Only 8 key figures are on GitHub, in `docs/figures/` (smaller copies; see README.md). The full set of 39 is kept
locally in `results/figures/` and is rebuilt with the `scripts/plot_*.py` scripts plus `scripts/add_figure_captions.py`.

```
reads -> mapping -> PAV calls -> SNPs and diversity -> check with deletions -> check with assemblies -> validation
      -> gene-loss types -> pan-genome -> relatedness -> function
```

## Stage 1. Overview

| # | File | Plot type | What it shows | Leads to |
|---|---|---|---|---|
| 01 | 01_overview_workflow__flowchart.png | **Workflow diagram** (flowchart) | The whole pipeline: reads, variants and assemblies feed one validation, then the analyses | What the genome looks like with all data layers |
| 02 | 02_overview_genome__circos_plot.png | **Circos plot** (circular genome plot) | The 5 chromosomes with rings for gene density, SNPs, diversity, absent genes, deletions and NLR genes; gene density drops and diversity rises at the centromeres | Start of the pipeline: mapping the reads |

## Stage 2. Read mapping (input: short reads of 8 accessions)

| # | File | Plot type | What it shows | Leads to |
|---|---|---|---|---|
| 03 | 03_mapping_quality__bar_chart_cleveland_dot_plot.png | **Bar charts** + **Cleveland dot plots** | 94-97% of reads map to Col-0; depth 7-12x; no sample failed | Mapped reads are good enough to measure depth per gene |

## Stage 3. Gene presence/absence (PAV) calls (input: mapped reads)

| # | File | Plot type | What it shows | Leads to |
|---|---|---|---|---|
| 04 | 04_pav_gene_depth__histograms.png | **Histograms** (small multiples) | Per-gene depth relative to each sample's median; absent genes form a peak at 0 | A cut-off near 0 can separate absent genes |
| 05 | 05_pav_cutoff_check__hexbin_sensitivity_plot.png | **Hexbin plots** + **sensitivity-analysis line plot** | Calls barely change when the cut-off moves | The calls are stable, so count them |
| 06 | 06_pav_calls_per_accession__grouped_bar_chart.png | **Grouped bar chart** | Absent / uncertain / ambiguous calls per accession; Col-0 (control) has 0 absent | What the absent genes are like |
| 07 | 07_pav_results__stacked_bar_ecdf_scatter.png | **Stacked bars**, **frequency bars**, **ECDF**, **scatter** | 173-452 absent genes per accession; most lost in one accession | See the 911 variable genes as a matrix |
| 08 | 08_pav_variable_genes__binary_heatmap.png | **Binary PAV heatmap** | The 911 variable genes x 8 accessions | Group them by pattern |
| 09 | 09_pav_clustered__clustered_heatmap_dendrogram.png | **Clustered heatmap** (with dendrograms) | Genes and accessions clustered by absence pattern | Check depth behind the calls |
| 10 | 10_pav_relative_depth__diverging_heatmap.png | **Heatmap** (diverging colours) | Relative depth of the variable genes | Which accessions share absences |
| 11 | 11_pav_shared_absences__jaccard_similarity_matrix.png | **Similarity matrix heatmap** (Jaccard) | Absent genes shared between pairs of accessions | Where in the genome they sit |
| 12 | 12_pav_along_chromosomes__density_heatmap.png | **Density heatmap** along chromosomes | Absent genes along the 5 chromosomes | Are these real deletions or mapping failures? Check with SNPs |

## Stage 4. SNPs and diversity (input: mapped reads)

| # | File | Plot type | What it shows | Leads to |
|---|---|---|---|---|
| 13 | 13_variants__bar_chart_pca_distance_heatmap.png | **Bar chart**, **PCA**, **distance heatmap**, **SNP density heatmap** | 472k-634k SNPs per accession; Cvi-0 most different | Is SNP density higher around absent genes? |
| 14 | 14_variants_snp_density__ecdf_dumbbell_plot.png | **ECDF** + **dumbbell plot** | Absent genes are not in unusually divergent regions | So absences are not just mapping failures; how diverse is the genome overall? |
| 15 | 15_variants_diversity__sliding_window_plot_violin_plot.png | **Sliding-window plot** + **violin plot** + **Mann-Whitney test** | Diversity pi 0.005 per bp (lower bound), highest around centromeres; NLR genes 2.1x more diverse (p = 2e-12) | Check the absent calls with deletions |

## Stage 5. Check with structural deletions (input: mapped reads, Delly)

| # | File | Plot type | What it shows | Leads to |
|---|---|---|---|---|
| 16 | 16_deletions_vs_absent__dumbbell_scatter_ecdf.png | **Dumbbell plot**, **scatter**, **ECDF of deletion lengths** | Absent genes lie inside Delly deletions 8-76 times more often | An independent check with assemblies |

## Stage 6. Check with genome assemblies (input: 4 assemblies)

| # | File | Plot type | What it shows | Leads to |
|---|---|---|---|---|
| 17 | 17_assembly_quality__busco_stacked_bar_cleveland_dot_plot.png | **100% stacked bar chart** (BUSCO) + **Cleveland dot plots** (QUAST) | All 4 assemblies 99.9-100% complete (4,311 genes), 120-122 Mb, N50 23-24 Mb | Complete assemblies, so a gene Liftoff cannot find is really missing |
| 18 | 18_assembly_liftoff__100pct_stacked_bar_chart.png | **100% stacked bar chart** | 89-94% of read-absent genes are also missing in the assembly | Compare every gene, not just absent ones |
| 19 | 19_assembly_reads_vs_liftoff__scatter_plots.png | **Scatter plots** (small multiples) | Read coverage vs assembly coverage for all genes | Combine all three methods |

## Stage 7. Validation summary (input: reads + Delly + assemblies)

| # | File | Plot type | What it shows | Leads to |
|---|---|---|---|---|
| 20 | 20_validation_three_methods__stacked_grouped_bar_chart.png | **Stacked + grouped bar charts** | Only 21-41 read-absent genes per accession lack support from both other methods | Are the errors biased by gene length? |
| 21 | 21_validation_length_bias__line_plot_confidence_intervals.png | **Line plot** + **dot plot with 95% CI** (Wilson) | Short genes called absent more often, but confirmed just as often | Where do the confirmed absences sit? |
| 22 | 22_validation_absent_genes__chromosome_ideogram.png | **Chromosome ideogram** (karyotype plot) | Confirmed absences cluster near centromeres | Look at one gene up close |
| 23 | 23_validation_FAR4_example__coverage_track.png | **Coverage track** (genome-browser style) | FAR4 in Cvi-0: no reads, a Delly deletion, not in the assembly | Beyond single genes: rearrangements |
| 24 | 24_validation_rearrangements__synteny_riparian_plot.png | **Synteny plot** (plotsr, riparian plot) | 92-95% of genomes in the same order; Chr4 inversion in Col-0 | Calls are validated; now: how are genes lost? |

## Stage 8. Types of gene change (input: validated PAV + SnpEff)

| # | File | Plot type | What it shows | Leads to |
|---|---|---|---|---|
| 25 | 25_geneloss_lost_vs_broken__grouped_stacked_bar_chart.png | **Grouped bar chart** + **100% stacked bars** | Broken genes are 6-13 times more common than deleted genes | Add possible gains |
| 26 | 26_geneloss_lost_broken_gained__bar_chart_error_bars_scatter.png | **Grouped bar chart with error bars** + **scatter** | Lost and broken are solid; read-depth gains follow data noise | Summarise across all 8: the pan-genome |

## Stage 9. Pan-genome (input: PAV calls of all 8)

| # | File | Plot type | What it shows | Leads to |
|---|---|---|---|---|
| 27 | 27_pangenome_categories__donut_chart_histogram.png | **Donut chart** + **gene frequency histogram** | Core 96.7%, soft-core 1.4%, shell 1.7%, private 0.07% | Which accessions share the variable genes? |
| 28 | 28_pangenome_shared_absences__upset_plot.png | **UpSet plot** | 389 genes lost in only one accession (145 in Cvi-0); 18 in all 7 | How the core changes with more accessions |
| 29 | 29_pangenome_core_genes__accumulation_rarefaction_curve.png | **Gene accumulation (rarefaction) curve** | Core still shrinking at 8 accessions | See all variable genes by pattern |
| 30 | 30_pangenome_gene_patterns__circular_dendrogram_heatmap.png | **Circular dendrogram with heatmap rings** | Variable genes clustered by which accessions lack them | Do the patterns follow family relationships? |

## Stage 10. Relatedness (input: SNPs + PAV)

| # | File | Plot type | What it shows | Leads to |
|---|---|---|---|---|
| 31 | 31_relatedness_snps__pcoa_plot.png | **PCoA plot** | Cvi-0 stands apart | Where the accessions come from |
| 32 | 32_relatedness_origins__geographic_map_upgma_dendrogram.png | **Geographic map** + **UPGMA dendrogram** | Cape Verde (Cvi-0) and Tajikistan (Kondara) split first | Do gene losses give the same tree? |
| 33 | 33_relatedness_pav_vs_snp__tanglegram_mantel_test.png | **Tanglegram** + **Mantel test** | Yes: r = 0.79 (p = 0.02); without Col-0 r = 0.90 (p = 0.007) | What kinds of genes are lost? |

## Stage 11. Function (input: variable genes + GO annotation + SnpEff)

| # | File | Plot type | What it shows | Leads to |
|---|---|---|---|---|
| 34 | 34_function_go_enrichment__dot_plot.png | **GO enrichment dot plot** (bubble plot) | Disease-resistance (NLR) genes most enriched (FDR 3e-36) | Which genes drive each term? |
| 35 | 35_function_go_genes__cnetplot_network.png | **Cnetplot** (gene-concept network) | Four separate gene groups: NLRs, RNase H, terpene synthases, Chr5 cluster | Same links as ribbons |
| 36 | 36_function_go_genes__chord_diagram.png | **GO chord diagram** (GOChord) | NLRs carry several defence terms at once | Where the NLRs sit |
| 37 | 37_function_nlr_genes__chromosome_ideogram_fisher_test.png | **Chromosome ideogram** + **Fisher's exact test** | 33% of NLRs are variable vs 3.2% of other genes; mostly in clusters | Are the NLRs that remain also changing faster? |
| 38 | 38_function_nlr_selection__pnps_dot_plot_ci_scatter_plot.png | **Dot plot with 95% CI** + **scatter plot** (pN/pS, permutation test) | pN/pS: NLR 0.49 vs core genes 0.25 (p = 1e-4); variable genes 0.60. Weaker constraint or diversifying selection, not proof of positive selection | Check the odd "fertilization" GO term |
| 39 | 39_function_chr5_cluster__coverage_track.png | **Coverage track** | The Chr5 "fertilization" cluster is a multi-copy region; treat that term with caution | End of the pipeline |

## Interactive

| File | Type | What it does |
|---|---|---|
| results/interactive/pav_explorer.html | **Interactive heatmap** (Plotly) | Hover and click 911 genes x 8 accessions; filters by chromosome, GO group, number of accessions |

## Statistical methods used in the figures

- **Mantel test** (33): correlation between two distance matrices, exact permutation p-value.
- **Fisher's exact test** (37): 2 x 2 table, NLR vs other genes, variable vs not.
- **Mann-Whitney U test** (15): gene diversity, NLR vs other genes.
- **Permutation test and bootstrap 95% CI** (38): pN/pS, NLR vs core genes.
- **Hypergeometric test with Benjamini-Hochberg FDR** (34-36): GO over-representation (clusterProfiler).
- **Wilson 95% confidence interval** (21): uncertainty of a proportion.
- **UPGMA** (average-linkage clustering; 30, 32, 33) and **PCoA / classical MDS** (31).
- **Rarefaction** (29): random orders of accessions, median and 5-95% band.
- **Nucleotide diversity pi** (15, 02) and **pN/pS with Nei-Gojobori site counts** (38).
