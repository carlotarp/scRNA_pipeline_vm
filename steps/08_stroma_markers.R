#!/usr/bin/env Rscript
##
## 08_stroma_markers.R - re-cluster the stromal cells and produce the evidence to annotate them
## Usage : Rscript pipeline/steps/08_stroma_markers.R
##         (run once for the resolution evidence, record the decision, run again for the markers)
## Reads : objects/leuko_annotated.rds
## Writes: phase 1: objects/stroma_grid.rds (cache) ; 08_stroma_markers/grid/ ResolutionGrid.png,
##                  SankeyPlot.html, cluster_sizes_per_resolution.csv
##         phase 2: objects/stroma.rds ; annotation/stroma.tsv (template) ;
##                  08_stroma_markers/ DimPlot_UMAP_Stroma.png, Dotplot_Stroma<block>.png,
##                  cluster_<cl>_FindMarkers[_filtered].csv, all_clusters_FindMarkers[_filtered].csv,
##                  top_expressed_genes.csv
## Checkpoint: decision subclustering.stroma.resolution (between phase 1 and phase 2)
##
file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
PIPELINE_DIR <- normalizePath(file.path(dirname(file_arg), ".."))
source(file.path(PIPELINE_DIR, "R/setup.R"))
out_dir <- step_dir("08_stroma_markers")
grid_dir <- step_dir("08_stroma_markers", "grid")

lineage_name <- "Stromal"
grid <- cfg_get("clustering.resolution_grid")
grid_cols <- paste0("RNA_snn_res.", grid)
cache_file <- run_object("stroma_grid")

# ---------------------------------------------------------------
# Phase 1: subset + Harmony + clustering at every grid resolution (cached)
# ---------------------------------------------------------------
if (file.exists(cache_file)) {
  cat("  Loading cached stroma grid:", cache_file, "\n")
  dwStroma <- readRDS(cache_file)
} else {
  dwAnnotated <- read_run_object("leuko_annotated")
  # Stromal cells are the same in lineage_annotated and leuko_annotated (07 only relabels/drops leukocytes)
  dwStroma <- generate_lineage_subset(dwAnnotated, lineage_name = lineage_name, resolution = grid,
                                      n_pcs = cfg_get("subclustering.stroma.n_pcs"))
  rm(dwAnnotated); invisible(gc())
  saveRDS(dwStroma, cache_file)

  # Evidence for the resolution choice
  plot_resolution_grid(dwStroma, grid_dir, reduction = "umap", resolutions = grid)
  sizes <- do.call(rbind, lapply(seq_along(grid), function(i) {
    tab <- table(dwStroma[[grid_cols[i], drop = TRUE]])
    data.frame(resolution = grid[i], cluster = names(tab), n_cells = as.integer(tab))
  }))
  write.csv(sizes, file.path(grid_dir, "cluster_sizes_per_resolution.csv"), row.names = FALSE)
}

# ---------------------------------------------------------------
# Checkpoint: stroma resolution (stops here until decided)
# ---------------------------------------------------------------
res <- get_decision("subclustering.stroma.resolution")
res_col <- paste0("RNA_snn_res.", res)
if (!res_col %in% colnames(dwStroma@meta.data)) {
  cat("  Resolution", res, "is not in the grid: clustering the cached graph at it\n")
  dwStroma <- FindClusters(dwStroma, resolution = res)
}
dwStroma$seurat_clusters <- dwStroma[[res_col, drop = TRUE]]
Idents(dwStroma) <- "seurat_clusters"
cat("  Stroma resolution:", res, "->", nlevels(dwStroma$seurat_clusters), "clusters\n")

# ---------------------------------------------------------------
# Phase 2: markers at the chosen resolution
# ---------------------------------------------------------------
plot_dimplot(dwStroma, reduction = "umap", group_by = "seurat_clusters",
             results_path = out_dir, filename = "DimPlot_UMAP_Stroma.png", label = TRUE)

for (block in names(load_marker_panel("stroma"))) {
  plot_marker_dotplot(dwStroma, marker_groups = load_marker_panel("stroma", block, object = dwStroma),
                      results_path = out_dir, filename = paste0("Dotplot_Stroma", block, ".png"),
                      group_by = "seurat_clusters")
}

markers_filt <- find_markers_for_clusters(dwStroma, levels(dwStroma$seurat_clusters), out_dir)
top_expressed_genes(dwStroma, group_by = "seurat_clusters", results_path = out_dir)

saveRDS(dwStroma, run_object("stroma"))
write_annotation_template(dwStroma, "stroma", markers_df = markers_filt)

cat("\n ---- FINISHED 08_stroma_markers ----\n",
    "Fill", annotation_file("stroma"), "(label per cluster, DROP for noise/doublets; evidence required),\n",
    "then run 09_stroma_annotate.R\n")
