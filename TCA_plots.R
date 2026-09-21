##
## Plotting functions for the tumor cluster analysis pipeline (3g, 4a).
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
# Generic categorical x categorical transition/concordance tile heatmap
# ---------------------------------------------------------------
plot_transition_heatmap <- function(counts_df, x_col, y_col, count_col = "n",
                                     title, xlab, ylab, filename, results_path,
                                     mode = c("gradient", "manual"),
                                     fill_mode = c("count", "pct"),
                                     label_mode = c("n", "n_pct"),
                                     pct_group_col = y_col,
                                     low_color = "white", high_color = "firebrick",
                                     group_col = NULL, group_colors = NULL,
                                     legend_name = "n cells",
                                     text_size = 2.8, tile_border = "white",
                                     width = 10, height = 8) {
  mode <- match.arg(mode)
  fill_mode <- match.arg(fill_mode)
  label_mode <- match.arg(label_mode)

  counts_df <- counts_df %>%
    dplyr::group_by(.data[[pct_group_col]]) %>%
    dplyr::mutate(.pct = round(.data[[count_col]] / sum(.data[[count_col]]) * 100, 1)) %>%
    dplyr::ungroup()

  counts_df$.label <- if (label_mode == "n_pct") {
    paste0(counts_df[[count_col]], "\n(", counts_df$.pct, "%)")
  } else {
    counts_df[[count_col]]
  }

  if (mode == "manual") {
    counts_df <- counts_df %>%
      dplyr::group_by(.data[[group_col]]) %>%
      dplyr::mutate(.scaled = (.data[[count_col]] - min(.data[[count_col]])) /
                       (max(.data[[count_col]]) - min(.data[[count_col]]) + 0.001)) %>%
      dplyr::ungroup() %>%
      dplyr::rowwise() %>%
      dplyr::mutate(.tile_color = colorRampPalette(c("white", group_colors[[.data[[group_col]]]]))(100)[round(.scaled * 99) + 1]) %>%
      dplyr::ungroup()

    p <- ggplot(counts_df, aes(x = .data[[x_col]], y = .data[[y_col]], fill = .tile_color)) +
      geom_tile(color = tile_border) +
      geom_text(aes(label = .label), size = text_size) +
      scale_fill_identity()
  } else {
    fill_col <- if (fill_mode == "pct") ".pct" else count_col
    p <- ggplot(counts_df, aes(x = .data[[x_col]], y = .data[[y_col]], fill = .data[[fill_col]])) +
      geom_tile(color = tile_border) +
      geom_text(aes(label = .label), size = text_size) +
      scale_fill_gradient(low = low_color, high = high_color, name = legend_name)
  }

  p <- p +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.background = element_rect(fill = "white", color = NA)) +
    labs(title = title, x = xlab, y = ylab)

  ggsave(paste0(results_path, filename), p, width = width, height = height, dpi = 300, bg = "white")
}

# ---------------------------------------------------------------
# Predicted-subtype composition stacked bars per sample, faceted by Subtype
# ---------------------------------------------------------------
plot_subtype_composition_facet <- function(object, results_path, call_col = "PAM50_predicted",
                                           facet_col = "Subtype", label = call_col, filename = NULL,
                                           colors = NULL) {
  if (is.null(filename)) {
    filename <- paste0("Composition_", call_col, "_bySample_facet", facet_col, ".png")
  }

  composition_df <- object@meta.data %>%
    group_by(orig.ident, .data[[call_col]], .data[[facet_col]]) %>%
    summarise(n = n(), .groups = "drop")

  fill_scale <- if (!is.null(colors)) scale_fill_manual(values = colors) else scale_fill_discrete()

  p <- ggplot(composition_df, aes(x = orig.ident, y = n, fill = .data[[call_col]])) +
    geom_col(position = "stack") +
    fill_scale +
    facet_wrap(vars(.data[[facet_col]]), scales = "free_x") +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste0(label, " composition per sample, by clinical ", facet_col),
         x = "Sample", y = "Number of cells", fill = label)

  ggsave(paste0(results_path, filename), p, width = 12, height = 6, dpi = 300, bg = "white")
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
# Generic score-distribution boxplot: x = cluster_col, y = score
# ---------------------------------------------------------------
plot_score_boxplot <- function(object, score_cols, results_path, filename,
                               cluster_col = "seurat_clusters", fill_col = NULL,
                               fill_labels = NULL, facet_col = NULL,
                               title, y_lab = "Score", colors = NULL) {
  select_cols <- unique(c(cluster_col, score_cols, fill_col, facet_col))

  score_long <- object@meta.data %>%
    select(all_of(select_cols)) %>%
    pivot_longer(cols = all_of(score_cols), names_to = "score_type", values_to = "score")

  if (is.null(fill_col) && !is.null(fill_labels)) {
    score_long$score_type <- unname(fill_labels[score_long$score_type])
  }

  fill_aes <- if (!is.null(fill_col)) fill_col else "score_type"

  p <- ggplot(score_long, aes(x = .data[[cluster_col]], y = score, fill = .data[[fill_aes]])) +
    geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.3, position = position_dodge(width = 0.8)) +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = title, x = "Cluster", y = y_lab, fill = if (!is.null(fill_col)) fill_col else "Score type")

  if (!is.null(facet_col)) p <- p + facet_wrap(vars(.data[[facet_col]]), scales = "free_x")
  if (!is.null(colors)) p <- p + scale_fill_manual(values = colors)

  ggsave(paste0(results_path, filename), p, width = 12, height = 6, dpi = 300, bg = "white")
}

# ---------------------------------------------------------------
# Cell type/call proportion barplot (mean +/- SEM, sample points overlaid)
# ---------------------------------------------------------------
plot_celltype_proportions_bar <- function(prop_df, x_var, results_path, filename,
                                          celltype_col = "celltype", colors = NULL,
                                          ncol_facets = 5) {
  summary_df <- prop_df %>%
    group_by(.data[[x_var]], .data[[celltype_col]]) %>%
    summarise(mean_prop = mean(proportion), sem = sd(proportion) / sqrt(n()), .groups = "drop")

  fill_scale <- if (!is.null(colors)) scale_fill_manual(values = colors) else scale_fill_discrete()

  p <- ggplot(summary_df, aes(x = .data[[x_var]], y = mean_prop, fill = .data[[x_var]])) +
    geom_col(alpha = 0.8) +
    geom_errorbar(aes(ymin = mean_prop - sem, ymax = mean_prop + sem), width = 0.2) +
    geom_jitter(data = prop_df, aes(x = .data[[x_var]], y = proportion), inherit.aes = FALSE,
               width = 0.15, height = 0, size = 1, alpha = 0.5, color = "black") +
    facet_wrap(as.formula(paste("~", celltype_col)), scales = "free_y", ncol = ncol_facets) +
    fill_scale +
    labs(x = x_var, y = "Cell proportion (mean ± SEM)", fill = x_var) +
    theme_bw(base_size = 11) +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1),
         strip.text = element_text(face = "bold"))

  ggsave(paste0(results_path, filename), p, width = 14, height = 10, dpi = 300, bg = "white")
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
                         limits = c(-1, 1), na.value = "white", name = "Diff. vs\nrest") +
    coord_flip(ylim = c(-1, 1)) +
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
# PROGENy pathway score distribution boxplots
# ---------------------------------------------------------------
plot_progeny_score_boxplots <- function(object, results_path) {
  # Long table: one row per cell x pathway, carrying the cell's cluster
  progeny_scores_df <- as.data.frame(t(as.matrix(GetAssayData(object, assay = "progeny", layer = "data"))))
  progeny_scores_df$seurat_clusters <- object$seurat_clusters

  progeny_long <- progeny_scores_df %>%
    pivot_longer(cols = -seurat_clusters, names_to = "pathway", values_to = "score")

  # One panel per pathway, clusters on the x axis
  p_by_cluster <- ggplot(progeny_long, aes(x = seurat_clusters, y = score, fill = seurat_clusters)) +
    geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.3) +
    facet_wrap(~pathway, scales = "free_y") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
          legend.position = "none") +
    labs(title = "PROGENy pathway score distribution by cluster",
         x = "Cluster", y = "Score")

  ggsave(paste0(results_path, "Boxplot_PROGENy_byCluster.png"), p_by_cluster,
         width = 16, height = 12, dpi = 300, bg = "white")

  # One panel per cluster, pathways on the x axis
  p_by_pathway <- ggplot(progeny_long, aes(x = pathway, y = score, fill = pathway)) +
    geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.3) +
    facet_wrap(~seurat_clusters, scales = "free_y") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
          legend.position = "none") +
    labs(title = "PROGENy pathway score distribution by pathway",
         x = "Pathway", y = "Score")

  ggsave(paste0(results_path, "Boxplot_PROGENy_byPathway.png"), p_by_pathway,
         width = 16, height = 12, dpi = 300, bg = "white")
}

# ---------------------------------------------------------------
# CopyKAT prediction composition bars (by sample and by cluster)
# ---------------------------------------------------------------
plot_copykat_composition <- function(object, results_path) {
  meta <- object@meta.data

  # Stacked proportion of CopyKAT calls (aneuploid / diploid / not defined) per sample
  comp_by_sample <- meta %>%
    group_by(orig.ident, copykat_prediction) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(orig.ident) %>%
    mutate(pct = n / sum(n) * 100) %>%
    ungroup()

  p_sample <- ggplot(comp_by_sample, aes(x = orig.ident, y = pct, fill = copykat_prediction)) +
    geom_col(position = "stack") +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "CopyKAT prediction composition per sample",
         x = "Sample", y = "% of cells", fill = "CopyKAT")

  ggsave(paste0(results_path, "BarPlot_CopyKATPrediction_bySample.png"), p_sample,
         width = 10, height = 6, dpi = 300, bg = "white")

  # Same breakdown, one bar per Seurat cluster
  comp_by_cluster <- meta %>%
    group_by(seurat_clusters, copykat_prediction) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(seurat_clusters) %>%
    mutate(pct = n / sum(n) * 100) %>%
    ungroup()

  p_cluster <- ggplot(comp_by_cluster, aes(x = seurat_clusters, y = pct, fill = copykat_prediction)) +
    geom_col(position = "stack") +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "CopyKAT prediction composition per cluster",
         x = "Cluster", y = "% of cells", fill = "CopyKAT")

  ggsave(paste0(results_path, "BarPlot_CopyKATPrediction_byCluster.png"), p_cluster,
         width = 10, height = 6, dpi = 300, bg = "white")
}

# ---------------------------------------------------------------
# Chromosomal instability (CIN_score/PGA_score) distribution by cluster
# ---------------------------------------------------------------
plot_copykat_cna_boxplot <- function(object, results_path, score_col = "CIN_score",
                                     y_lab = "CIN score (SD of CNA profile)", filename = NULL) {
  meta <- object@meta.data[!is.na(object@meta.data[[score_col]]), ]

  if (is.null(filename)) filename <- paste0("Boxplot_", score_col, "_byCluster.png")

  p <- ggplot(meta, aes(x = seurat_clusters, y = .data[[score_col]], fill = seurat_clusters)) +
    geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.3) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    labs(title = paste0(score_col, " distribution by cluster"), x = "Cluster", y = y_lab)

  ggsave(paste0(results_path, filename), p,
         width = 10, height = 6, dpi = 300, bg = "white")
}
