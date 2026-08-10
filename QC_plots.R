##
## Plotting functions for scRNA-seq Quality Control pipeline
## Source this file from 1_seurat_quality_control.R .
##

library(ggplot2)
library(Seurat)
library(dplyr)
library(patchwork)

# ---------------------------------------------------------------
# nCount vs nFeature scatter
# ---------------------------------------------------------------
plot_gene_expression_dist <- function(tumor, name, results_path) {
  png(paste0(results_path, "nCountRNA_vs_nFeatureRNA_", name, ".png"), width = 800, height = 600)
  p <- FeatureScatter(tumor, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", smooth = TRUE)
  print(p)
  dev.off()
}

# ---------------------------------------------------------------
# nFeature histogram with QC cutoff lines
# ---------------------------------------------------------------
plot_nfeature_ncount_corr <- function(tumor, name, results_path) {
  png(paste0(results_path, "nFeatureRNA_", name, ".png"), width = 800, height = 600)
  hist(tumor@meta.data$nFeature_RNA,
       main = paste("nFeature_RNA:", name),
       xlab = "Number of Features (Genes)",
       ylab = "Frequency",
       breaks = 50,
       col = "lightgrey",
       border = "white")

  # QC cutoffs (high)
  abline(v = 3500, col = "orange", lwd = 2, lty = 2)
  abline(v = 5500, col = "darkorange", lwd = 2, lty = 2)
  abline(v = 7500, col = "red", lwd = 2, lty = 2)

  # Percentiles (low)
  abline(v = 100, col = "darkgreen", lwd = 3, lty = 1)

  legend("topright",
         legend = c("3500", "5500", "7500", "100"),
         col = c("orange", "darkorange", "red", "darkgreen"),
         lwd = 2, lty = c(2, 2, 2, 1), cex = 0.75)
  dev.off()
}

# ---------------------------------------------------------------
# Violin plots of QC metrics, PRE-filter, per sample
# ---------------------------------------------------------------
plot_qc_vlnplot <- function(object, name, results_path,
                             features = c("nFeature_RNA", "nCount_RNA",
                                          "percent.mt", "percent.hb")) {
  png(paste0(results_path, "VlnPlot_QCmetrics_", name, ".png"), width = 1200, height = 600)
  p <- VlnPlot(object, features = features, ncol = length(features), pt.size = 0.1)
  print(p)
  dev.off()
}

# ---------------------------------------------------------------
# percent.mt vs nFeature_RNA scatter
# ---------------------------------------------------------------
plot_mito_vs_nfeature <- function(object, name, results_path,
                                   mito_cut, feature_lower_cut, feature_upper_cut) {
  df <- object@meta.data
  df$pass_filter <- df$percent.mt < mito_cut &
    df$nFeature_RNA > feature_lower_cut &
    df$nFeature_RNA < feature_upper_cut

  p <- ggplot(df, aes(x = nFeature_RNA, y = percent.mt, color = pass_filter)) +
    geom_point(size = 0.6, alpha = 0.6) +
    geom_vline(xintercept = c(feature_lower_cut, feature_upper_cut),
               linetype = "dashed", color = "black") +
    geom_hline(yintercept = mito_cut, linetype = "dashed", color = "black") +
    scale_color_manual(values = c("TRUE" = "steelblue", "FALSE" = "grey70"),
                        name = "Pass filter") +
    labs(title = paste("Mito % vs nFeature -", name),
         x = "nFeature_RNA", y = "percent.mt") +
    theme_bw() +
    theme(panel.grid = element_blank())

  ggsave(paste0(results_path, "MitoVsNFeature_", name, ".png"), p,
         width = 7, height = 5, dpi = 300)
}

# ---------------------------------------------------------------
# scDblFinder score distribution by class
# ---------------------------------------------------------------
plot_doublet_scores <- function(object, name, results_path) {
  df <- object@meta.data

  p <- ggplot(df, aes(x = scDblFinder.score, fill = scDblFinder.class)) +
    geom_density(alpha = 0.5) +
    scale_fill_manual(values = c("singlet" = "steelblue", "doublet" = "firebrick")) +
    labs(title = paste("scDblFinder score distribution -", name),
         x = "scDblFinder score", y = "Density") +
    theme_bw() +
    theme(panel.grid = element_blank())

  ggsave(paste0(results_path, "DoubletScore_", name, ".png"), p,
         width = 7, height = 5, dpi = 300)
}

# ---------------------------------------------------------------
# Violin plots comparing samples AFTER individual QC (post-merge)
# ---------------------------------------------------------------
plot_vln_compare_samples <- function(object, results_path,
                                      features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                                      group_by = "orig.ident") {
  png(paste0(results_path, "VlnPlot_compare_samples_postQC.png"), width = 1200, height = 600)
  p <- VlnPlot(object, features = features, group.by = group_by,
               ncol = length(features), pt.size = 0)
  print(p)
  dev.off()
}

# ---------------------------------------------------------------
# DimPlot by sample - generic, works for "pca", "harmony" or any reduction
# ---------------------------------------------------------------
plot_dimplot_by_sample <- function(object, reduction, results_path,
                                    group_by = "orig.ident") {
  png(paste0(results_path, "DimPlot_", reduction, "_by_sample.png"), width = 800, height = 600)
  p <- DimPlot(object, reduction = reduction, group.by = group_by) +
    ggtitle(paste("Reduction:", reduction))
  print(p)
  dev.off()
}

# ---------------------------------------------------------------
# QC funnel plot: cells retained per filtering step, one panel per sample
# ---------------------------------------------------------------
plot_qc_funnel <- function(cell_counts_df, results_path) {

  qc_colors <- c("skyblue", "lightgreen", "salmon", "orange", "purple", "cyan")

  samples <- unique(cell_counts_df$SampleID)
  plot_list <- list()

  for (sample in samples) {
    sample_df <- cell_counts_df[cell_counts_df$SampleID == sample, ]
    loaded_cells <- sample_df$Cell_counts[sample_df$Feature_discard == "Loaded cells"]
    qc_df <- sample_df[sample_df$Feature_discard != "Loaded cells", ]
    remaining_cells <- min(qc_df$Cell_counts)

    counts <- loaded_cells
    labels <- "Loaded cells"
    colors <- "black"
    prev <- loaded_cells

    for (i in seq_len(nrow(qc_df))) {
      discarded <- prev - qc_df$Cell_counts[i]
      counts <- c(counts, discarded)
      labels <- c(labels, qc_df$Feature_discard[i])
      colors <- c(colors, qc_colors[i])
      prev <- qc_df$Cell_counts[i]
    }

    counts <- c(counts, remaining_cells)
    labels <- c(labels, "Remaining cells")
    colors <- c(colors, "grey")

    # % always relative to the initial loaded_cells
    perc <- round((counts / loaded_cells) * 100, 1)

    df_plot <- data.frame(
      label = factor(labels, levels = labels),
      count = counts,
      perc = perc,
      color = colors
    )

    p <- ggplot(df_plot, aes(x = label, y = count, fill = label)) +
      geom_col(show.legend = FALSE) +
      geom_text(aes(label = paste0(count, " (", perc, "%)")),
                vjust = -0.3, size = 3) +
      scale_fill_manual(values = setNames(df_plot$color, df_plot$label)) +
      labs(title = paste("Sample", sample), x = NULL, y = "Cells") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            panel.grid = element_blank())

    plot_list[[sample]] <- p
  }

  combined <- wrap_plots(plot_list, ncol = 4)
  ggsave(paste0(results_path, "QC_cell_counts_subplots.png"), combined, width = 24, height = 18, dpi = 300)
}

# ---------------------------------------------------------------
# Percentage of cells filtered per sample, colored by tumor subtype
# ---------------------------------------------------------------
plot_pct_filtered_by_subtype <- function(cell_counts_df, sample_annot, results_path,
                                          id_col = "scRNAseq_ID", subtype_col = "Subtype") {

  subtype_colors <- c("TNBC" = "#e41a1c", "Her2+" = "#377eb8", "Lum" = "#D0A3D9")

  samples <- unique(cell_counts_df$SampleID)
  pct_list <- c()
  remaining_list <- c()

  for (sample in samples) {
    sample_df <- cell_counts_df[cell_counts_df$SampleID == sample, ]
    loaded_cells <- sample_df$Cell_counts[sample_df$Feature_discard == "Loaded cells"]
    remaining_cells <- min(sample_df$Cell_counts[sample_df$Feature_discard != "Loaded cells"])
    pct_list <- c(pct_list, (loaded_cells - remaining_cells) / loaded_cells * 100)
    remaining_list <- c(remaining_list, remaining_cells)
  }

  df_plot <- data.frame(sample = samples, pct_filtered = pct_list, remaining = remaining_list)
  df_plot <- df_plot[order(df_plot$pct_filtered), ]
  df_plot$sample <- factor(df_plot$sample, levels = df_plot$sample)

  # match subtype - trim whitespace to avoid key mismatches (e.g. "Her2+ " vs "Her2+")
  subtype_map <- setNames(trimws(sample_annot[[subtype_col]]), sample_annot[[id_col]])
  df_plot$subtype <- subtype_map[as.character(df_plot$sample)]
  df_plot$subtype[is.na(df_plot$subtype) | !(df_plot$subtype %in% names(subtype_colors))] <- "Other"

  fill_values <- c(subtype_colors, "Other" = "grey50")

  p <- ggplot(df_plot, aes(x = sample, y = pct_filtered, fill = subtype)) +
    geom_col() +
    geom_text(aes(label = paste0(round(pct_filtered, 1), "%"), y = pct_filtered / 2),
              color = "white", fontface = "bold", size = 3.5) +
    geom_text(aes(label = remaining, y = pct_filtered * 1.02),
              vjust = 0, size = 3) +
    scale_fill_manual(values = fill_values, name = "Subtype") +
    labs(title = "Percentage of cells removed by QC per sample",
         x = NULL, y = "Percentage of cells filtered (%)") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(paste0(results_path, "QC_percentage_filtered_with_labels.png"), p, width = 12, height = 6, dpi = 300)
}

# ---------------------------------------------------------------
# scDblFinder.score vs nCount_RNA scatter & nFeature_RNA, colored by singlet/doublet
# ---------------------------------------------------------------
plot_doublet_scatter <- function(df, name, results_path) {

  c <- ggplot(df, aes(x = nCount_RNA, y = scDblFinder.score, color = scDblFinder.class)) +
    geom_point(size = 0.6, alpha = 0.6) +
    scale_color_manual(values = c("singlet" = "steelblue", "doublet" = "firebrick")) +
    labs(title = paste("scDblFinder score vs nCount -", name),
         x = "nCount_RNA", y = "scDblFinder score", color = "scDblFinder class") +
    theme_bw() +
    theme(panel.grid = element_blank())

  f <- ggplot(df, aes(x = nFeature_RNA, y = scDblFinder.score, color = scDblFinder.class)) +
  geom_point(size = 0.6, alpha = 0.6) +
  scale_color_manual(values = c("singlet" = "steelblue", "doublet" = "firebrick")) +
  labs(title = paste("scDblFinder score vs nFeature_RNA -", name),
       x = "nFeature_RNA", y = "scDblFinder score", color = "scDblFinder class") +
  theme_bw() +
  theme(panel.grid = element_blank())

  combined <- (c | f) + patchwork::plot_annotation(title = paste("scDblFinder score vs QC metrics -", name))
  ggsave(paste0(results_path, "DoubletScatter_", name, ".png"), combined,
         width = 15, height = 5, dpi = 300)
}