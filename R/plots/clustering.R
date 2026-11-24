##
## Clustering plots: elbow, resolution grid, UMAPs, composition, heatmap (steps 2-3).
##

library(networkD3)
library(htmlwidgets)
library(ggsci)

# ---------------------------------------------------------------
# Elbow plot with the chosen number of PCs (co3) marked
# ---------------------------------------------------------------
plot_elbow <- function(object, co3, results_path) {
  p <- ElbowPlot(object, ndims = length(object[["pca"]]@stdev)) +
    geom_vline(xintercept = co3, linetype = "dashed", color = "red") +
    ggtitle(paste0("Elbow plot (chosen PCs = ", co3, ")")) +
    theme_bw()

  ggsave(bg = "white", paste0(results_path, "ElbowPlot.png"), p, width = 7, height = 5, dpi = 300)
}

# ---------------------------------------------------------------
# Generic DimPlot (wraps the legend into columns past 10 group_by levels).
# For group_by = "lineage", uses lineage_colors (utils.R) unless an explicit
# `cols` is passed. For group_by = "celltype", orders the levels by lineage
# (instead of the default alphabetical) and uses the fixed celltype_colors
# (utils.R), so each celltype keeps its color in the full object and in the
# lineage subsets.
# ---------------------------------------------------------------
plot_dimplot <- function(object, reduction, group_by, results_path, filename,
                          label = FALSE, cols = NULL, legend_ncol = NULL,
                          width = 8, height = 8) {
  object[[group_by]] <- factor(object[[group_by]][, 1])
  n_levels <- nlevels(object[[group_by]][, 1])

  if (is.null(cols) && group_by == "lineage") {
    cols <- lineage_colors[levels(object[[group_by]][, 1])]
  } else if (is.null(cols) && group_by == "celltype" && "lineage" %in% colnames(object@meta.data)) {
    ordered_levels <- object@meta.data %>%
      dplyr::distinct(celltype, lineage) %>%
      dplyr::arrange(lineage, celltype) %>%
      dplyr::pull(celltype) %>%
      as.character() %>%
      unique()   # a celltype (e.g. "Noise") can appear under several lineages
    object[[group_by]] <- factor(object[[group_by]][, 1], levels = ordered_levels)
    # Fixed celltype_colors (utils.R); celltypes missing from it fall back to
    # pal_igv (deduplicated: it repeats colors 23-24 at positions 26-27)
    cols <- celltype_colors[ordered_levels]
    missing <- is.na(cols)
    if (any(missing)) {
      cols[missing] <- setdiff(unique(pal_igv("default")(51)), cols)[seq_len(sum(missing))]
    }
    names(cols) <- ordered_levels
  }

  p <- DimPlot(object, reduction = reduction, group.by = group_by,
               cols = cols, label = label, repel = label) +
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

  ggsave(bg = "white", paste0(results_path, filename), p, width = width, height = height, dpi = 300)
}

# ---------------------------------------------------------------
# Generic FeaturePlot
# ---------------------------------------------------------------
plot_featureplot <- function(object, reduction, features, results_path, filename,
                              raster = FALSE) {
  features <- features[features %in% rownames(object)]  # e.g. PTPRC is not in the Xenium panel (XeniumPanel run)
  if (length(features) == 0) { cat("  - FeaturePlot skipped (no requested feature in the object):", filename, "\n"); return(invisible(NULL)) }
  p <- FeaturePlot(object, reduction = reduction, features = features, raster = raster)

  n_features <- length(features)
  ggsave(bg = "white", paste0(results_path, filename), p,
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
  ggsave(bg = "white", paste0(results_path, "ResolutionGrid.png"), combined,
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

  ggsave(bg = "white", paste0(results_path, filename), p, width = 10, height = 5, dpi = 300)
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

  ggsave(bg = "white", paste0(results_path, filename), p,
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

  object <- subset(object, downsample = 50)   # 50 cells per cluster; genes scaled across these cells
  object[["RNA"]]$scale.data <- NULL
  object <- ScaleData(object, features = unique(top$gene), verbose = FALSE)
  p <- DoHeatmap(object, features = top$gene) + NoLegend()
  ggsave(paste0(results_path, "Integrated_Top", n, "gene_Heatmap.png"), p, width = 2000 / 150, height = 1600 / 150, dpi = 150, bg = "white")
}
