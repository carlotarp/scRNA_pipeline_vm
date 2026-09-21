##
## Plotting functions for the clustering pipeline (2a, 2b, 3g).
##

library(ggplot2)
library(Seurat)
library(dplyr)
library(patchwork)
library(networkD3)
library(htmlwidgets)
library(ggsci)

# ---------------------------------------------------------------
# Elbow plot with the chosen number of PCs (co3) marked
# ---------------------------------------------------------------
plot_elbow <- function(object, co3, results_path) {
  p <- ElbowPlot(object, ndims = length(object[["pca_decontX"]]@stdev)) +
    geom_vline(xintercept = co3, linetype = "dashed", color = "red") +
    ggtitle(paste0("Elbow plot (chosen PCs = ", co3, ")")) +
    theme_bw()

  ggsave(paste0(results_path, "ElbowPlot.png"), p, width = 7, height = 5, dpi = 300)
}

# ---------------------------------------------------------------
# Generic DimPlot (wraps the legend into columns past 10 group_by levels)
# ---------------------------------------------------------------
plot_dimplot <- function(object, reduction, group_by, results_path, filename,
                          label = FALSE, cols = NULL, legend_ncol = NULL,
                          width = 8, height = 8) {
  object[[group_by]] <- factor(object[[group_by]][, 1])
  n_levels <- nlevels(object[[group_by]][, 1])

  p <- DimPlot(object, reduction = reduction, group.by = group_by,
               cols = cols, label = label) +
    ggtitle(paste0(reduction, " - grouped by ", group_by))

  if (n_levels > 10) {
    ncol_legend <- if (!is.null(legend_ncol)) legend_ncol else ceiling(n_levels / 15)
    p <- p +
      guides(color = guide_legend(ncol = ncol_legend, override.aes = list(size = 3))) +
      theme(legend.text = element_text(size = 8),
            legend.title = element_text(size = 9),
            legend.key.size = unit(0.4, "cm"))
    width <- width + 1.5 * ncol_legend
  }

  ggsave(paste0(results_path, filename), p, width = width, height = height, dpi = 300)
}

# ---------------------------------------------------------------
# Generic FeaturePlot
# ---------------------------------------------------------------
plot_featureplot <- function(object, reduction, features, results_path, filename,
                              raster = FALSE) {
  p <- FeaturePlot(object, reduction = reduction, features = features, raster = raster)

  n_features <- length(features)
  ggsave(paste0(results_path, filename), p,
         width = 6 * min(3, n_features), height = 6 * ceiling(n_features / 3),
         dpi = 300, limitsize = FALSE)
}

# ---------------------------------------------------------------
# Resolution grid search
# ---------------------------------------------------------------
plot_resolution_grid <- function(object, results_path,
                                  reduction = "umap",
                                  resolutions = c(0.2, 0.4, 0.6, 0.8, 1.0)) {

  plot_list <- list()
  links <- data.frame()
  prev_res <- NULL
  prev_clusters <- NULL

  for (res in resolutions) {
    object <- FindClusters(object, resolution = res, verbose = FALSE)
    current_clusters <- Idents(object)
    n_clusters <- length(unique(current_clusters))

    p <- DimPlot(object, reduction = reduction, label = TRUE) +
      ggtitle(paste0("resolution = ", res, "  (", n_clusters, " clusters)")) +
      theme(legend.position = "none")
    plot_list[[as.character(res)]] <- p

    # Sankey links: can only be built starting from the second resolution
    if (!is.null(prev_clusters)) {
      df <- data.frame(source = prev_clusters, target = current_clusters) %>%
        dplyr::group_by(source, target) %>%
        dplyr::summarise(value = dplyr::n(), .groups = "drop")
      df$source <- paste0("R", prev_res, "_C", df$source)
      df$target <- paste0("R", res, "_C", df$target)
      links <- rbind(links, df)
    }

    prev_res <- res
    prev_clusters <- current_clusters
  }

  # Save UMAP resolution grid
  ncol_grid <- min(3, length(resolutions))
  combined <- wrap_plots(plot_list, ncol = ncol_grid)
  ggsave(paste0(results_path, "ResolutionGrid.png"), combined,
         width = 6 * ncol_grid, height = 6 * ceiling(length(resolutions) / ncol_grid),
         dpi = 300, limitsize = FALSE)

  # Save sankey (only if at least 2 resolutions were tested)
  if (nrow(links) > 0) {
    nodes <- data.frame(name = unique(c(links$source, links$target)))
    links$IDsource <- match(links$source, nodes$name) - 1
    links$IDtarget <- match(links$target, nodes$name) - 1

    sankey_plot <- sankeyNetwork(
      Links = links, Nodes = nodes,
      Source = "IDsource", Target = "IDtarget", Value = "value",
      NodeID = "name", sinksRight = FALSE, fontSize = 14, nodeWidth = 30
    )
    htmlwidgets::saveWidget(sankey_plot, paste0(results_path, "SankeyPlot.html"), selfcontained = FALSE)
  }
}

# ---------------------------------------------------------------
# Composition per sample / per subtype
# ---------------------------------------------------------------
plot_composition <- function(object, group_by, results_path, filename,
                             cluster_col = "seurat_clusters", colors = NULL) {
  df <- object@meta.data %>%
    dplyr::count(.data[[cluster_col]], .data[[group_by]])

  fill_scale <- if (!is.null(colors)) scale_fill_manual(values = colors) else scale_fill_igv()

  p <- ggplot(df, aes(x = .data[[cluster_col]], y = n, fill = .data[[group_by]])) +
    geom_col() +
    fill_scale +
    labs(title = paste0("Composition by ", group_by),
         x = "Cluster", y = "Number of cells") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(paste0(results_path, filename), p, width = 10, height = 5, dpi = 300)
}


# ---------------------------------------------------------------
# QC metrics violin plot, generic grouping
# ---------------------------------------------------------------
plot_vln_qc_by_group <- function(object, results_path, filename,
                                  qc_vars = c("nFeature_RNA", "percent.mt",
                                              "nCount_RNA", "percent.hb"),
                                  group_by = NULL) {
  ncol_grid <- ceiling(sqrt(length(qc_vars)))
  nrow_grid <- ceiling(length(qc_vars) / ncol_grid)

  p <- VlnPlot(object, features = qc_vars, group.by = group_by, pt.size = 0, ncol = ncol_grid)

  ggsave(paste0(results_path, filename), p,
         width = 5 * ncol_grid, height = 4 * nrow_grid, dpi = 300)
}

# ---------------------------------------------------------------
# Differential expression heatmap
# ---------------------------------------------------------------

plot_heatmap <- function(object, results_path, tumor_markers, n=10){
  tumor_markers %>%
      group_by(cluster) %>%
      filter(avg_log2FC > 1) %>%
      slice_head(n = n) %>%
      ungroup() -> top

  object <- ScaleData(object, features = unique(top$gene))
  object <- subset(object, downsample = 50)

  png(file=paste0(results_path,'Integrated_Top', n, 'gene_Heatmap.png'), width = 2000, height = 1600, res = 150)
  p <- DoHeatmap(object, features = top$gene) + NoLegend()
  print(p)
  dev.off()
}