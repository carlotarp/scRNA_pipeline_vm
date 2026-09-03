##
## Plotting functions for the DecontX contamination-correction step.
## Source this file from 1b_decontx_quality_control.R
##

library(ggplot2)
library(dplyr)

# ---------------------------------------------------------------
# Genome-wide scatter: average expression before vs after DecontX
# ---------------------------------------------------------------
plot_decontx_scatter <- function(gene_change_df, results_path) {
  p <- ggplot(gene_change_df, aes(x = avg_before, y = avg_after)) +
    geom_point(size = 0.4, alpha = 0.3, color = "steelblue") +
    geom_abline(slope = 1, intercept = 0, color = "firebrick", linetype = "dashed") +
    theme_bw() + theme(panel.grid = element_blank()) +
    labs(title = "Genome-wide expression: before vs after DecontX",
         subtitle = "Points below the red line = genes reduced by correction",
         x = "Average expression (before)", y = "Average expression (after)")

  ggsave(paste0(results_path, "GenomeWide_scatter_beforeAfter.png"), p,
         width = 7, height = 7, dpi = 300, bg = "white")
}

# ---------------------------------------------------------------
# MA-plot: expression level vs DecontX correction magnitude
# ---------------------------------------------------------------
plot_decontx_maplot <- function(gene_change_df, results_path) {
  top_labeled <- gene_change_df %>% slice_max(order_by = avg_drop, n = 15)

  p <- ggplot(gene_change_df, aes(x = avg_before, y = avg_drop)) +
    geom_point(size = 0.4, alpha = 0.3, color = "grey50") +
    geom_point(data = top_labeled, color = "firebrick", size = 1.2) +
    geom_text(data = top_labeled, aes(label = gene), size = 3, vjust = -0.6, check_overlap = TRUE) +
    theme_bw() + theme(panel.grid = element_blank()) +
    labs(title = "MA-plot: expression level vs DecontX correction magnitude",
         subtitle = "Top 15 most affected genes highlighted",
         x = "Average expression (before)", y = "avg_before - avg_after") +
    ylim(c(0, 2)) + xlim(c(0, 3))

  ggsave(paste0(results_path, "GenomeWide_MAplot.png"), p, width = 9, height = 7, dpi = 300, bg = "white")
}

# ---------------------------------------------------------------
# Histogram of percentage expression reduction across all genes
# ---------------------------------------------------------------
plot_decontx_hist <- function(gene_change_df, results_path) {
  p <- ggplot(gene_change_df %>% filter(!is.na(avg_drop_pct)), aes(x = avg_drop_pct)) +
    geom_histogram(bins = 60, fill = "steelblue", alpha = 0.7) +
    theme_bw() + theme(panel.grid = element_blank()) +
    labs(title = "Distribution of % expression reduction across all genes",
         x = "% reduction (avg_before -> avg_after)", y = "N genes")

  ggsave(paste0(results_path, "GenomeWide_pctDrop_histogram.png"), p,
         width = 8, height = 5, dpi = 300, bg = "white")
}

# ---------------------------------------------------------------
# Dumbbell plot: top 30 genes most reduced by DecontX
# ---------------------------------------------------------------
plot_decontx_dumbbell <- function(gene_change_df, results_path) {
  top30 <- gene_change_df %>%
    slice_max(order_by = avg_drop, n = 30) %>%
    mutate(gene = reorder(gene, avg_drop))

  p <- ggplot(top30) +
    geom_segment(aes(x = avg_before, xend = avg_after, y = gene, yend = gene), color = "grey60") +
    geom_point(aes(x = avg_before, y = gene), color = "firebrick", size = 2.5) +
    geom_point(aes(x = avg_after, y = gene), color = "steelblue", size = 2.5) +
    theme_bw() + theme(panel.grid.minor = element_blank()) +
    labs(title = "Top 30 genes most affected by DecontX: before (red) vs after (blue)",
         x = "Average expression", y = NULL)

  ggsave(paste0(results_path, "GenomeWide_top30_dumbbell.png"), p,
         width = 8, height = 9, dpi = 300, bg = "white")
}
