#!/usr/bin/env Rscript
##
## 06_leukocytes_markers.R - re-cluster the leukocytes and produce the evidence to annotate them
## Usage : Rscript pipeline/steps/06_leukocytes_markers.R
##         (run once for the resolution evidence, record the decision, run again for the markers)
## Reads : objects/lineage_annotated.rds
## Writes: phase 1: objects/leukocytes_grid.rds (cache) ; 06_leukocytes_markers/grid/ ResolutionGrid.png,
##                  SankeyPlot.html, cluster_sizes_per_resolution.csv
##         phase 2: objects/leukocytes.rds ; annotation/leukocytes.tsv (template) ;
##                  06_leukocytes_markers/ DimPlot_UMAP_Leukocytes.png, Dotplot_Leukocytes<block>.png,
##                  cluster_<cl>_FindMarkers[_filtered].csv, all_clusters_FindMarkers[_filtered].csv,
##                  top_expressed_genes.csv
## Checkpoint: decision subclustering.leukocytes.resolution (between phase 1 and phase 2)
##
file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
PIPELINE_DIR <- normalizePath(file.path(dirname(file_arg), ".."))
source(file.path(PIPELINE_DIR, "R/setup.R"))
out_dir <- step_dir("06_leukocytes_markers")
grid_dir <- step_dir("06_leukocytes_markers", "grid")

lineage_name <- "Leukocytes"
grid <- cfg_get("clustering.resolution_grid")
grid_cols <- paste0("RNA_snn_res.", grid)
cache_file <- run_object("leukocytes_grid")

# ---------------------------------------------------------------
# Phase 1: subset + Harmony + clustering at every grid resolution (cached)
# ---------------------------------------------------------------
if (file.exists(cache_file)) {
  cat("  Loading cached leukocyte grid:", cache_file, "\n")
  dwLeukocytes <- readRDS(cache_file)
} else {
  dwAnnotated <- read_run_object("lineage_annotated")
  dwLeukocytes <- generate_lineage_subset(dwAnnotated, lineage_name = lineage_name, resolution = grid,
                                          n_pcs = cfg_get("subclustering.leukocytes.n_pcs"))
  rm(dwAnnotated); invisible(gc())
  saveRDS(dwLeukocytes, cache_file)

  # Evidence for the resolution choice
  plot_resolution_grid(dwLeukocytes, grid_dir, reduction = "umap", resolutions = grid)
  sizes <- do.call(rbind, lapply(seq_along(grid), function(i) {
    tab <- table(dwLeukocytes[[grid_cols[i], drop = TRUE]])
    data.frame(resolution = grid[i], cluster = names(tab), n_cells = as.integer(tab))
  }))
  write.csv(sizes, file.path(grid_dir, "cluster_sizes_per_resolution.csv"), row.names = FALSE)
}

# ---------------------------------------------------------------
# Checkpoint: leukocyte resolution (stops here until decided)
# ---------------------------------------------------------------
res <- get_decision("subclustering.leukocytes.resolution")
res_col <- paste0("RNA_snn_res.", res)
if (!res_col %in% colnames(dwLeukocytes@meta.data)) {
  cat("  Resolution", res, "is not in the grid: clustering the cached graph at it\n")
  dwLeukocytes <- FindClusters(dwLeukocytes, resolution = res)
}
dwLeukocytes$seurat_clusters <- dwLeukocytes[[res_col, drop = TRUE]]
Idents(dwLeukocytes) <- "seurat_clusters"
cat("  Leukocyte resolution:", res, "->", nlevels(dwLeukocytes$seurat_clusters), "clusters\n")

# ---------------------------------------------------------------
# Phase 2: markers at the chosen resolution
# ---------------------------------------------------------------
plot_dimplot(dwLeukocytes, reduction = "umap", group_by = "seurat_clusters",
             results_path = out_dir, filename = "DimPlot_UMAP_Leukocytes.png", label = TRUE)

for (block in names(load_marker_panel("leukocytes"))) {
  plot_marker_dotplot(dwLeukocytes, marker_groups = load_marker_panel("leukocytes", block, object = dwLeukocytes),
                      results_path = out_dir, filename = paste0("Dotplot_Leukocytes", block, ".png"),
                      group_by = "seurat_clusters")
}

markers_filt <- find_markers_for_clusters(dwLeukocytes, levels(dwLeukocytes$seurat_clusters), out_dir)
top_expressed_genes(dwLeukocytes, group_by = "seurat_clusters", results_path = out_dir)

saveRDS(dwLeukocytes, run_object("leukocytes"))
write_annotation_template(dwLeukocytes, "leukocytes", markers_df = markers_filt)

cat("\n ---- FINISHED 06_leukocytes_markers ----\n",
    "Fill", annotation_file("leukocytes"), "(label per cluster, DROP for noise/doublets; evidence required),\n",
    "then run 07_leukocytes_annotate.R\n")
