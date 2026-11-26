#!/usr/bin/env Rscript
# GO enrichment of variable genes (Phase 7), run in conda env `go` (clusterProfiler 4.18, org.At.tair.db 3.22).
#
# Run:
#     Rscript scripts/go_enrichment.R
#
# Inputs (read only): results/pav/pav_long.tsv.gz, results/variants/lof_genes.tsv.gz
# Outputs: results/go/go_enrichment.tsv (all significant terms), results/figures/34_function_go_enrichment__dot_plot.png,
#          results/figures/35_function_go_genes__cnetplot_network.png
#
# Two gene sets, each tested against the same background (all 27,444 genes in the PAV list):
#   absent = called absent in >= 1 accession (911 genes)
#   broken = called present but with a HIGH-impact variant in >= MIN_BROKEN accessions other than Col-0
# Test: one-sided hypergeometric (enrichGO), BH FDR < 0.05, terms with 10-500 annotated background genes.
# Only genes with a GO annotation count; BP terms are made less redundant with simplify() (similarity 0.7).
# GO annotations include electronic ones and are uneven across genes; enrichment shows over-representation, not cause.

suppressMessages({
  library(clusterProfiler)
  library(org.At.tair.db)
  library(enrichplot)
  library(ggplot2)
})

MIN_BROKEN <- 2
FDR <- 0.05
dir.create("results/go", showWarnings = FALSE, recursive = TRUE)
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

pav <- read.delim(gzfile("results/pav/pav_long.tsv.gz"))[, c("gene", "accession", "call")]
universe <- unique(pav$gene)
absent <- unique(pav$gene[pav$call == "absent"])
lof <- read.delim(gzfile("results/variants/lof_genes.tsv.gz"))
lof <- lof[lof$pav_call == "present" & lof$accession != "Col-0" & lof$gene %in% universe, ]
n_acc <- tapply(lof$accession, lof$gene, function(x) length(unique(x)))
broken <- names(n_acc)[n_acc >= MIN_BROKEN]
sets <- list(absent = absent, broken = broken)
stopifnot(length(universe) == 27444, length(absent) == 911)

go_map <- suppressMessages(AnnotationDbi::select(org.At.tair.db, keys = universe, keytype = "TAIR", columns = "GO"))
annotated <- unique(go_map$TAIR[!is.na(go_map$GO)])
cat(sprintf("background: %d genes, %d with a GO annotation (%.0f%%)\n", length(universe), length(annotated),
            100 * length(annotated) / length(universe)))
for (s in names(sets)) {
  cat(sprintf("%-7s %5d genes, %5d with GO (%.0f%%)\n", s, length(sets[[s]]), sum(sets[[s]] %in% annotated),
              100 * mean(sets[[s]] %in% annotated)))
}
cat(sprintf("overlap absent/broken: %d genes\n", length(intersect(absent, broken))))

# gene label = TAIR symbol, or the AT ID when the gene has no symbol (setReadable would give NA)
sym <- suppressMessages(AnnotationDbi::mapIds(org.At.tair.db, keys = universe, keytype = "TAIR", column = "SYMBOL",
                                              multiVals = "first"))
sym[is.na(sym)] <- names(sym)[is.na(sym)]
dup <- duplicated(sym) | duplicated(sym, fromLast = TRUE)
sym[dup] <- paste0(sym[dup], " (", names(sym)[dup], ")")
to_sym <- function(ids) sapply(strsplit(ids, "/"), function(g) paste(sym[g], collapse = "/"))

res <- list()
cc_bp <- NULL
for (ont in c("BP", "MF", "CC")) {
  cc <- compareCluster(geneClusters = sets, fun = "enrichGO", OrgDb = org.At.tair.db, keyType = "TAIR",
                       ont = ont, universe = universe, pAdjustMethod = "BH", pvalueCutoff = FDR,
                       qvalueCutoff = 1, minGSSize = 10, maxGSSize = 500)
  if (is.null(cc) || nrow(as.data.frame(cc)) == 0) {
    cat(ont, ": no significant terms\n")
    next
  }
  if (ont == "BP") {
    cc <- simplify(cc, cutoff = 0.7, by = "p.adjust", select_fun = min)
    cc_bp <- cc
  }
  df <- as.data.frame(cc)
  df$ont <- ont
  df$symbols <- to_sym(df$geneID)
  res[[ont]] <- df
}
tab <- do.call(rbind, res)
rownames(tab) <- NULL
tab <- tab[, c("Cluster", "ont", "ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "Count", "geneID", "symbols")]
names(tab)[1] <- "set"
write.table(tab, "results/go/go_enrichment.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nsignificant terms (FDR < 0.05):\n")
print(table(tab$set, tab$ont))
cat("\ntop 5 per set and ontology:\n")
for (s in names(sets)) for (ont in c("BP", "MF", "CC")) {
  x <- tab[tab$set == s & tab$ont == ont, ]
  if (nrow(x) == 0) next
  x <- head(x[order(x$p.adjust), ], 5)
  for (i in seq_len(nrow(x))) {
    cat(sprintf("  %-6s %s  %-55s %7s  FDR %.1e\n", s, ont, substr(x$Description[i], 1, 55), x$GeneRatio[i],
                x$p.adjust[i]))
  }
}

# fig 34: top terms per set and ontology as a dot plot
ratio <- function(r) sapply(strsplit(r, "/"), function(v) as.numeric(v[1]) / as.numeric(v[2]))
top <- do.call(rbind, lapply(split(tab, list(tab$set, tab$ont), drop = TRUE),
                             function(x) head(x[order(x$p.adjust), ], 8)))
keep <- unique(top[, c("ont", "Description")])
plot_df <- merge(tab, keep, by = c("ont", "Description"))
plot_df$ratio <- ratio(plot_df$GeneRatio)
plot_df$ont <- factor(plot_df$ont, levels = c("BP", "MF", "CC"),
                      labels = c("Biological process", "Molecular function", "Cellular component"))
ord <- aggregate(ratio ~ Description, plot_df, max)
plot_df$Description <- factor(plot_df$Description, levels = ord$Description[order(ord$ratio)])
plot_df$set <- factor(plot_df$set, levels = c("absent", "broken"),
                      labels = c(sprintf("absent\n(%d genes)", length(absent)),
                                 sprintf("broken in >= %d\n(%d genes)", MIN_BROKEN, length(broken))))
p <- ggplot(plot_df, aes(x = set, y = Description, size = ratio, colour = -log10(p.adjust))) +
  geom_point() +
  facet_grid(ont ~ ., scales = "free_y", space = "free_y") +
  scale_colour_gradient(low = "#9ec5f0", high = "#1c4f8f", name = "-log10 FDR") +
  scale_size_area(max_size = 6, name = "share of\nset genes") +
  scale_y_discrete(labels = label_wrap_gen(55)) +
  labs(x = NULL, y = NULL,
       title = "GO terms over-represented among variable genes",
       subtitle = "Background: all 27,444 genes; BH FDR < 0.05; top 8 terms per set and ontology") +
  theme_minimal(base_size = 9) +
  theme(strip.text.y = element_text(angle = 0, hjust = 0, face = "bold"), panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"), plot.background = element_rect(fill = "#fcfcfb", colour = NA))
h <- max(4, 1.2 + 0.22 * length(unique(plot_df$Description)))
ggsave("results/figures/34_function_go_enrichment__dot_plot.png", p, width = 8, height = h, dpi = 300)

# fig 35: which absent genes drive the top terms (BP and MF together). Terms whose gene list overlaps an
# already chosen term by more than half (Jaccard > 0.5) are skipped, so each of the 5 shows a different gene group.
ego <- enrichGO(absent, OrgDb = org.At.tair.db, keyType = "TAIR", ont = "ALL", universe = universe,
                pAdjustMethod = "BH", pvalueCutoff = FDR, qvalueCutoff = 1, minGSSize = 10, maxGSSize = 500)
if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
  r <- ego@result[ego@result$p.adjust < FDR, ]
  r <- r[order(r$p.adjust), ]
  chosen <- character(0)
  for (i in seq_len(nrow(r))) {
    g <- strsplit(r$geneID[i], "/")[[1]]
    ok <- all(vapply(chosen, function(d) {
      h <- strsplit(r$geneID[r$Description == d], "/")[[1]]
      length(intersect(g, h)) / length(union(g, h)) <= 0.5
    }, logical(1)))
    if (ok) chosen <- c(chosen, r$Description[i])
    if (length(chosen) == 5) break
  }
  cat("\nfig19 terms:", paste(chosen, collapse = " | "), "\n")
  ego@result$geneID <- to_sym(ego@result$geneID)
  p2 <- cnetplot(ego, showCategory = chosen) +
    ggtitle("Absent genes behind 5 distinct enriched GO terms (gene symbol, or AT ID if none)") +
    theme(plot.background = element_rect(fill = "#fcfcfb", colour = NA), plot.title.position = "plot",
          plot.margin = margin(10, 10, 10, 20))
  ggsave("results/figures/35_function_go_genes__cnetplot_network.png", p2, width = 11, height = 9, dpi = 300)
} else {
  cat("absent set: no significant terms, fig 35 not made\n")
}
cat("wrote results/go/go_enrichment.tsv and results/figures/34_function_go_enrichment__dot_plot.png (+ fig 35 if terms)\n")
