##
## Plotting functions for scRNA-seq Cell Annotation pipeline
## Source this file from 3a_seurat_cell_annotation.R & 3b_seurat_cell_subtype_annotation.R .
##
library(ggplot2)
library(Seurat)
library(patchwork)

# ---------------------------------------------------------------
# Generic DimPlot
# ---------------------------------------------------------------
plot_dimplot <- function(object, reduction, group_by, results_path, filename,
                          cols = NULL, label = FALSE) {
  p <- DimPlot(object, reduction = reduction, group.by = group_by,
               cols = cols, label = label) +
    ggtitle(paste0(reduction, " - grouped by ", group_by))

  ggsave(paste0(results_path, filename), p, width = 8, height = 8, dpi = 300)
}

# ---------------------------------------------------------------
# Generic FeaturePlot
# ---------------------------------------------------------------
plot_featureplot <- function(object, reduction, features, results_path, filename,
                              raster = FALSE) {
  p <- FeaturePlot(object, reduction = reduction, features = features, raster = raster)

  n_features <- length(features)
  ggsave(paste0(results_path, filename), p,
         dpi = 300, limitsize = FALSE)
}

# ---------------------------------------------------------------
# Marker gene dot plot per cluster/group - size
# ---------------------------------------------------------------
plot_marker_dotplot <- function(object, marker_groups, results_path,
                                 filename = "marker_dotplot.png",
                                 group_by = "seurat_clusters",
                                 title = "Marker expression by group") {
  markers <- unique(unlist(marker_groups))
  markers <- markers[markers %in% rownames(object)]  # evita error si algún gen no está en el objeto

  # Map each gene -> "(lineage)_(gene)" label
  gene_label <- sapply(markers, function(g) {
    grp <- names(marker_groups)[sapply(marker_groups, function(x) g %in% x)][1]
    paste0(grp, "_", g)
  })
  gene_label <- gene_label[markers]

  p <- DotPlot(object, features = markers, group.by = group_by) +
    coord_flip() +
    scale_color_gradient(low = "blue", high = "red") +
    scale_x_discrete(labels = gene_label) +
    labs(title = title) +
    theme(plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA),
          axis.text.x = element_text(angle = 0, hjust = 0.5))

  ggsave(paste0(results_path, filename), p,
         width = 10, height = 10, dpi = 300, limitsize = FALSE,
         bg = "white")
}

# ---------------------------------------------------------------
# Resolution grid search
# ---------------------------------------------------------------
plot_resolution_grid <- function(object, results_path,
                                  resolutions = c(0.4, 0.6, 0.8, 1.0),
                                  filename = "resolution_grid.png") {

  plot_list <- list()
  for (res in resolutions) {
    object <- FindClusters(object, resolution = res, verbose = FALSE)
    n_clusters <- length(unique(Idents(object)))

    p <- DimPlot(object, reduction = "umap", label = TRUE) +
      ggtitle(paste0("resolution = ", res, "  (", n_clusters, " clusters)")) +
      theme(legend.position = "none")

    plot_list[[as.character(res)]] <- p
  }

  ncol_grid <- length(resolutions)%/%2
  combined <- wrap_plots(plot_list, ncol = ncol_grid)
  ggsave(paste0(results_path, filename), combined,
         width = 6 * ncol_grid, height = 6 * ceiling(length(resolutions) / ncol_grid),
         dpi = 300, limitsize = FALSE)

}