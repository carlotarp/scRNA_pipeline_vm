##
## Plotting functions for the clinical covariate association analysis (4a).
##

library(ggplot2)
library(dplyr)
library(ggpubr)
library(rstatix)

# ---------------------------------------------------------------
# Sample x group_by composition heatmap + stacked proportion barplot
# ---------------------------------------------------------------
generate_heatmap <- function(object, filename, group_by = "celltype", n_clusters = 4) {

    # Build counts matrix
    counts_df <- object@meta.data %>%
    group_by(across(all_of(c("orig.ident", group_by)))) %>%
    summarise(n = n(), .groups = "drop")

    counts_matrix <- counts_df %>%
    pivot_wider(names_from = orig.ident, values_from = n, values_fill = 0) %>%
    column_to_rownames(group_by) %>%
    as.matrix()

    # Normalize counts: proportion by sample, then scaled 0-100 by sample for color
    # (each column/sample independently rescaled 0-100 across its group_by rows,
    # so color intensity is relative to that sample's own min/max, not the
    # group_by row's min/max across samples)
    counts_matrix_prop <- apply(counts_matrix, 2, function(x) x / sum(x) * 100)
    counts_matrix_color <- apply(counts_matrix_prop, 2, function(x) (x - min(x)) / (max(x) - min(x) + 1e-9) * 100)

    # Cluster samples by composition
    col_dist <- dist(t(counts_matrix_prop), method = "euclidean")
    col_clust <- hclust(col_dist, method = "ward.D2")

# Build sample annotation
sample_annotation <- dwAnnotated@meta.data %>%
  distinct(orig.ident, Subtype, pCR, IGG, TLS) %>%
  remove_rownames() %>%
  column_to_rownames("orig.ident")

sample_annotation <- sample_annotation[colnames(counts_matrix), , drop = FALSE]

# Clinical covariate palettes (utils.R)
top_annotation <- HeatmapAnnotation(
  Subtype = sample_annotation$Subtype,
  pCR = sample_annotation$pCR,
  IGG = sample_annotation$IGG,
  TLS = sample_annotation$TLS,
  col = list(Subtype = pam50_subtype_colors, pCR = pcr_colors, IGG = igg_level_colors, TLS = tls_colors),
  annotation_name_side = "left"
)

# Build main heatmap
ht <- Heatmap(counts_matrix_color,
               name = paste0("Relative\n(by ", group_by, ")"),
               col = colorRamp2(c(0, 100), c("white", "#f2a154")),
               top_annotation = top_annotation,
               cluster_columns = col_clust,
               cluster_rows = TRUE,
               column_split = n_clusters,
               cell_fun = function(j, i, x, y, width, height, fill) {
                 grid::grid.text(counts_matrix[i, j], x, y, gp = grid::gpar(fontsize = 7))
               },
               row_names_gp = grid::gpar(fontsize = 9),
               column_names_gp = grid::gpar(fontsize = 9),
               show_heatmap_legend = TRUE)

    # Build stacked proportion barplot
    prop_matrix <- t(counts_matrix_prop / 100)

    celltype_colors <- setNames(pal_igv("default")(ncol(prop_matrix)), colnames(prop_matrix))

    barplot_annotation <- HeatmapAnnotation(
    Proportion = anno_barplot(prop_matrix, gp = grid::gpar(fill = celltype_colors),
                                bar_width = 1, height = unit(4, "cm")),
    which = "column", show_annotation_name = TRUE
    )

    # Empty heatmap, used to stack the barplot below
    ht_barplot <- Heatmap(matrix(nrow = 0, ncol = ncol(counts_matrix),
                                dimnames = list(NULL, colnames(counts_matrix))),
                            top_annotation = barplot_annotation, height = unit(0, "cm"))

    # Draw combined figure
    celltype_legend <- Legend(labels = colnames(prop_matrix),
                                legend_gp = grid::gpar(fill = celltype_colors),
                                title = group_by)

    png(paste0(results_GEMX_CA_path, filename),
        width = 3800, height = 3400, res = 300)
    draw(ht %v% ht_barplot,
        annotation_legend_list = list(celltype_legend),
        merge_legend = TRUE,
        heatmap_legend_side = "right",
        annotation_legend_side = "right")
    dev.off()
}

# ---------------------------------------------------------------
# Stacked proportion barplot: composition of clinical_col within each group_col level
# ---------------------------------------------------------------
plot_association_barplot <- function(meta, group_col, clinical_col, results_path, filename,
                                     colors = NULL, group_label = group_col, clinical_label = clinical_col) {
  df <- meta %>%
    dplyr::filter(!is.na(.data[[group_col]]), !is.na(.data[[clinical_col]])) %>%
    dplyr::count(.data[[group_col]], .data[[clinical_col]], name = "n") %>%
    dplyr::group_by(.data[[group_col]]) %>%
    dplyr::mutate(pct = n / sum(n) * 100) %>%
    dplyr::ungroup()

  fill_scale <- if (!is.null(colors)) scale_fill_manual(values = colors) else scale_fill_discrete()

  p <- ggplot(df, aes(x = .data[[group_col]], y = pct, fill = .data[[clinical_col]])) +
    geom_col(position = "stack") +
    fill_scale +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste0(clinical_label, " composition by ", group_label),
         x = group_label, y = "% of cells", fill = clinical_label)

  ggsave(paste0(results_path, filename), p, width = 10, height = 6, dpi = 300, bg = "white")
}

# ---------------------------------------------------------------
# Barplot of per-group Spearman correlations vs an ordinal clinical covariate
# ---------------------------------------------------------------
plot_correlation_barplot <- function(cor_df, results_path, filename, title, group_label = "Group") {
  cor_df$sig <- ifelse(!is.na(cor_df$p_val) & cor_df$p_val < 0.05, "p < 0.05", "n.s.")

  p <- ggplot(cor_df, aes(x = reorder(group, rho), y = rho, fill = sig)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = c("p < 0.05" = "firebrick", "n.s." = "grey70")) +
    theme_bw() + theme(panel.grid = element_blank()) +
    labs(title = title, x = group_label, y = "Spearman rho (proportion vs clinical score)", fill = NULL)

  ggsave(paste0(results_path, filename), p, width = 8, height = 6, dpi = 300, bg = "white")
}

# ---------------------------------------------------------------
# Cell type proportion boxplots (per-sample proportion), faceted by cell type
# ---------------------------------------------------------------
plot_celltype_proportions <- function(prop_df, x_var, results_path, filename,
                                      celltype_col = "celltype", colors = NULL,
                                      ncol_facets = 5, show_pvalues = TRUE) {
  fill_scale <- if (!is.null(colors)) scale_fill_manual(values = colors) else scale_fill_discrete()

  prop_df[[x_var]] <- factor(prop_df[[x_var]])

  p <- ggplot(prop_df, aes(x = .data[[x_var]], y = proportion, fill = .data[[x_var]])) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.15, height = 0, size = 1.3, alpha = 0.7, color = "black") +
    facet_wrap(as.formula(paste("~", celltype_col)), scales = "free_y", ncol = ncol_facets) +
    fill_scale +
    labs(x = x_var, y = "Cell proportion", fill = x_var) +
    theme_bw(base_size = 11) +
    theme(panel.grid = element_blank(),
         axis.text.x   = element_text(angle = 45, hjust = 1),
         strip.text    = element_text(face = "bold"))

  if (nlevels(prop_df[[x_var]]) >= 2) {
    stat_df <- prop_df %>%
      dplyr::group_by(.data[[celltype_col]]) %>%
      rstatix::wilcox_test(stats::reformulate(x_var, "proportion")) %>%
      rstatix::add_significance("p") %>%
      dplyr::filter(p.signif != "ns")

    if (nrow(stat_df) > 0 & show_pvalues) {
      stat_df <- stat_df %>% rstatix::add_xy_position(x = x_var, scales = "free_y")
      p <- p + ggpubr::stat_pvalue_manual(stat_df, label = "p.signif", tip.length = 0.01)
    }
  }

  ggsave(paste0(results_path, filename), p, width = 14, height = 10, dpi = 300, bg = "white")
}
