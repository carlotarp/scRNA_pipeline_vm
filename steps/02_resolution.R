#!/usr/bin/env Rscript
##
## 02_resolution.R - evidence for the global clustering resolution (elbow, resolution grid, cluster sizes)
## Usage : Rscript pipeline/steps/02_resolution.R
## Reads : objects/qc_merged.rds
## Writes: 02_resolution/ElbowPlot.png, ResolutionGrid.png, SankeyPlot.html, cluster_sizes_by_resolution.csv
## Checkpoint: none here; it produces the evidence for clustering.resolution (needed by 03_cluster.R):
##   Rscript pipeline/decide.R clustering.resolution <value> "<reason>"
##
file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
PIPELINE_DIR <- normalizePath(file.path(dirname(file_arg), ".."))
source(file.path(PIPELINE_DIR, "R/setup.R"))
out_dir <- step_dir("02_resolution")

obj <- read_run_object("qc_merged")
resolutions <- unlist(cfg_get("clustering.resolution_grid"))

n_pcs <- select_n_pcs(obj)
cat("\n PCs:", n_pcs, "\n")
plot_elbow(obj, n_pcs, out_dir)

# UMAP (PCA init) and graph on the Harmony space, as in 03_cluster.R
obj <- run_umap_pca(obj, reduction = "harmony", dims = 1:n_pcs)
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:n_pcs, assay = "RNA", graph.name = c("RNA_nn", "RNA_snn"))

# Cluster sizes at every grid resolution (table rewritten after each resolution)
sizes <- data.frame()
for (r in resolutions) {
  obj <- FindClusters(obj, resolution = r, graph.name = "RNA_snn", verbose = FALSE)
  cl <- obj[[paste0("RNA_snn_res.", r), drop = TRUE]]
  tab <- table(cluster = cl, sample = obj$orig.ident)
  n <- rowSums(tab)
  top <- apply(tab, 1, which.max)
  sizes <- rbind(sizes, data.frame(resolution = r, n_clusters = nrow(tab), cluster = rownames(tab),
                                   n_cells = as.integer(n), pct_cells = round(100 * n / sum(n), 2),
                                   top_sample = colnames(tab)[top],
                                   top_sample_pct = round(100 * tab[cbind(seq_along(top), top)] / n, 1)))
  write.csv(sizes, file.path(out_dir, "cluster_sizes_by_resolution.csv"), row.names = FALSE)
  cat(sprintf("  resolution %s: %d clusters\n", r, nrow(tab)))
}

# UMAP per resolution + Sankey of cluster splits
plot_resolution_grid(obj, results_path = out_dir, reduction = "umap", resolutions = resolutions)

cat("\n ---- 02_resolution done ----\n Review ResolutionGrid.png, SankeyPlot.html and cluster_sizes_by_resolution.csv, then:\n",
    "  Rscript pipeline/decide.R clustering.resolution <value> \"<reason>\"\n")
