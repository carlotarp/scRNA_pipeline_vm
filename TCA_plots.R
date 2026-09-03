##
## Plotting functions for scRNA-seq Tumor Cell Annotation pipeline.
## Source this file from 3g_seurat_tumor_analysis.R
##

library(ggplot2)
library(corrplot)
library(dplyr)
library(tidyr)
library(gt)

# ---------------------------------------------------------------
# HTML marker tables per cluster (gt format)
# ---------------------------------------------------------------
plot_marker_tables <- function(results_path) {
  marker_files <- list.files(results_path, pattern = "^cluster_.*_FindMarkers\\.csv$", full.names = FALSE)

  for (f in marker_files) {
    df <- read.csv(paste0(results_path, f))
    colnames(df)[1] <- "gene"
    cl_name <- gsub("cluster_(.*)_FindMarkers\\.csv", "\\1", basename(f))

    table_cl <- df %>%
      select(gene, avg_log2FC, pct.1, pct.2, p_val_adj) %>%
      arrange(desc(avg_log2FC)) %>%
      gt() %>%
      fmt_number(columns = c(avg_log2FC, pct.1, pct.2), decimals = 2) %>%
      fmt_scientific(columns = p_val_adj, decimals = 2) %>%
      tab_header(title = paste0("Top 20 markers - Cluster ", cl_name))

    gtsave(table_cl, paste0(results_path, "Table_FindMarkers_cluster", cl_name, ".html"))
  }
}

# ---------------------------------------------------------------
# PAM50 predicted vs clinical subtype concordance heatmap
# ---------------------------------------------------------------
plot_pam50_concordance <- function(concordance_table, results_path) {
  concordance_df <- as.data.frame(concordance_table) %>%
    group_by(Predicted) %>%
    mutate(pct = Freq / sum(Freq) * 100)

  p <- ggplot(concordance_df, aes(x = Clinical, y = Predicted, fill = pct)) +
    geom_tile(color = "white") +
    geom_text(aes(label = paste0(Freq, "\n(", round(pct, 1), "%)")), size = 3.5) +
    scale_fill_gradient(low = "white", high = "steelblue", name = "% of\npredicted") +
    labs(title = "PAM50 predicted vs clinical Subtype concordance",
         x = "Clinical Subtype", y = "PAM50 predicted") +
    theme_minimal() +
    theme(panel.grid = element_blank())

  ggsave(paste0(results_path, "Heatmap_PAM50_vs_ClinicalSubtype.png"), p,
         width = 7, height = 5, dpi = 300, bg = "white")
}

# ---------------------------------------------------------------
# PAM50 composition stacked bars, faceted by clinical subtype
# ---------------------------------------------------------------
plot_pam50_composition_facet <- function(object, results_path) {
  composition_df <- object@meta.data %>%
    group_by(orig.ident, PAM50_predicted, Subtype) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(orig.ident, Subtype)

  p <- ggplot(composition_df, aes(x = orig.ident, y = n, fill = PAM50_predicted)) +
    geom_col(position = "stack") +
    facet_wrap(~Subtype, scales = "free_x") +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "PAM50 predicted composition per sample, by clinical subtype",
         x = "Sample", y = "Number of cells", fill = "PAM50 predicted")

  ggsave(paste0(results_path, "Composition_PAM50_bySample_facetSubtype.png"), p,
         width = 12, height = 6, dpi = 300, bg = "white")
}

# ---------------------------------------------------------------
# Cell cycle score distribution by cluster, faceted by subtype
# ---------------------------------------------------------------
plot_cellcycle_boxplot <- function(object, results_path) {
  cellcycle_long <- object@meta.data %>%
    select(seurat_clusters, Subtype, S.Score, G2M.Score) %>%
    pivot_longer(cols = c(S.Score, G2M.Score), names_to = "score_type", values_to = "score")

  p <- ggplot(cellcycle_long, aes(x = seurat_clusters, y = score, fill = score_type)) +
    geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.3, position = position_dodge(width = 0.8)) +
    facet_wrap(~Subtype, scales = "free_x") +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Cell cycle score distribution by cluster, by clinical subtype",
         x = "Cluster", y = "Score", fill = "Score type")

  ggsave(paste0(results_path, "CellCycle_scores_boxplot_byCluster_facetSubtype.png"), p,
         width = 12, height = 6, dpi = 300, bg = "white")
}

# ---------------------------------------------------------------
# Differential pathway activity heatmap (PROGENy or DecoupleR-PROGENy)
# ---------------------------------------------------------------
plot_progeny_diff_heatmap <- function(diff_df, results_path, filename, title) {
  p <- ggplot(diff_df, aes(x = cluster, y = pathway, fill = avg_diff)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                         name = "Diff. activity\n(cluster vs rest)") +
    theme_minimal() +
    theme(panel.grid = element_blank(), plot.background = element_rect(fill = "white", color = NA)) +
    labs(title = title, x = "Cluster", y = NULL)

  ggsave(paste0(results_path, filename), p, width = 8, height = 6, dpi = 300, bg = "white")
}

# ---------------------------------------------------------------
# Differential score barplot, one facet per cluster
# ---------------------------------------------------------------
plot_differential_barplot <- function(csv_path, feature_col, title, filename, results_path) {
  df <- read.csv(csv_path)

  p <- ggplot(df, aes(x = .data[[feature_col]], y = avg_diff, fill = avg_diff)) +
    geom_col() +
    facet_wrap(~cluster) +
    scale_fill_gradient2(low = "firebrick", mid = "white", high = "forestgreen", midpoint = 0,
                         name = "Diff. vs\nrest") +
    coord_flip() +
    theme_bw() + theme(panel.grid = element_blank()) +
    labs(title = title, x = NULL, y = "avg_diff (cluster vs rest)")

  ggsave(paste0(results_path, filename), p, width = 12, height = 8, dpi = 300, bg = "white")
}

# ---------------------------------------------------------------
# Spearman correlation matrix heatmap across all continuous scores
# ---------------------------------------------------------------
plot_correlation_matrix <- function(cor_matrix, results_path) {
  png(paste0(results_path, "Correlation_matrix_allScores.png"),
      width = 3000, height = 3000, res = 300)
  corrplot(cor_matrix, method = "color", type = "upper", order = "hclust",
           tl.col = "black", tl.srt = 45, addCoef.col = "black", number.cex = 0.5,
           title = "Correlation across all continuous scores (Spearman)", mar = c(0, 0, 2, 0))
  dev.off()
}

# ---------------------------------------------------------------
# Scatter plot with linear fit for a pair of correlated variables
# ---------------------------------------------------------------
plot_correlation_scatter <- function(data, var1, var2, results_path) {
  r_val <- cor(data[[var1]], data[[var2]], method = "spearman")
  p_val <- cor.test(data[[var1]], data[[var2]], method = "spearman")$p.value

  p <- ggplot(data, aes(x = .data[[var1]], y = .data[[var2]])) +
    geom_point(size = 0.4, alpha = 0.2, color = "steelblue") +
    geom_smooth(method = "lm", color = "firebrick", se = TRUE) +
    theme_bw() + theme(panel.grid = element_blank()) +
    labs(title = paste0(var1, " vs ", var2),
         subtitle = paste0("Spearman r = ", round(r_val, 3), ", p = ", signif(p_val, 3)),
         x = var1, y = var2)

  ggsave(paste0(results_path, "Scatter_", var1, "_vs_", var2, ".png"), p,
         width = 6, height = 5, dpi = 300, bg = "white")
}

# ---------------------------------------------------------------
# Celltype transition tile plot (before vs after DecontX correction)
# ---------------------------------------------------------------
plot_celltype_transition <- function(object, results_path) {
  transition_df <- data.frame(
    cell = colnames(object),
    celltype_before = as.character(object$celltype_cont),
    celltype_after = as.character(object$celltype)
  )

  trans_counts <- transition_df %>%
    dplyr::count(celltype_before, celltype_after, name = "n") %>%
    dplyr::group_by(celltype_before) %>%
    dplyr::mutate(pct = round(n / sum(n) * 100, 1)) %>%
    dplyr::ungroup()

  p <- ggplot(trans_counts, aes(x = celltype_after, y = celltype_before, fill = pct)) +
    geom_tile(color = "white") +
    geom_text(aes(label = pct), size = 2.8) +
    scale_fill_gradient(low = "white", high = "firebrick", name = "% of\nbefore") +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.background = element_rect(fill = "white", color = NA)) +
    labs(title = "Celltype transition (%): before (celltype_cont) vs after (celltype, DecontX)",
         x = "Celltype after DecontX", y = "Celltype before DecontX")

  ggsave(paste0(results_path, "Transition_celltype_beforeAfter_decontX.png"), p,
         width = 10, height = 8, dpi = 300, bg = "white")
}
