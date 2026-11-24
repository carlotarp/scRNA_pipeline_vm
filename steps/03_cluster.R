#!/usr/bin/env Rscript
##
## 03_cluster.R - global UMAP + clustering at the chosen resolution, QC by cluster, FindAllMarkers
## Usage : Rscript pipeline/steps/03_cluster.R
## Reads : objects/qc_merged.rds
## Writes: objects/clustered.rds ; 03_cluster/ (UMAP/PCA plots, composition, QC violins,
##         cells_per_cluster.tsv, umap_embeddings.tsv, clustered_metadata.tsv,
##         pairwisett_nfeaturerna_subtype.csv, FindAllMarkers.tsv, Integrated_Top10gene_Heatmap.png)
## Checkpoint: clustering.resolution (evidence from 02_resolution.R)
##
file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
PIPELINE_DIR <- normalizePath(file.path(dirname(file_arg), ".."))
source(file.path(PIPELINE_DIR, "R/setup.R"))
out_dir <- step_dir("03_cluster")

res <- as.numeric(get_decision("clustering.resolution"))   # stops here until decided
if (length(res) != 1 || is.na(res)) stop("clustering.resolution must be one number, got: ", paste(res, collapse = ", "))

obj <- read_run_object("qc_merged")
n_pcs <- select_n_pcs(obj)
cat("\n PCs:", n_pcs, "| resolution:", res, "\n")

obj <- run_umap_pca(obj, reduction = "harmony", dims = 1:n_pcs)

# global_clusters keeps the global clusters once seurat_clusters is overwritten by the lineage sub-clustering
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:n_pcs, assay = "RNA", graph.name = c("RNA_nn", "RNA_snn"))
obj <- FindClusters(obj, resolution = res, graph.name = "RNA_snn", cluster.name = "global_clusters")
obj$seurat_clusters <- obj$global_clusters   # cluster.name skips seurat_clusters
has_subtype <- "Subtype" %in% colnames(obj@meta.data) && any(!is.na(obj$Subtype))   # clinical sheet is optional

# Embeddings
for (red in c("harmony", "umap")) {
  tag <- if (red == "harmony") "PCA" else "UMAP"
  if (has_subtype) plot_dimplot(obj, reduction = red, group_by = "Subtype", results_path = out_dir,
                                filename = paste0("DimPlot_", tag, "_bySubtype.png"))
  plot_dimplot(obj, reduction = red, group_by = "seurat_clusters", results_path = out_dir,
               filename = paste0("DimPlot_", tag, "_byCluster.png"))
  plot_featureplot(obj, reduction = red, features = "PTPRC", results_path = out_dir,
                   filename = paste0("FeaturePlot_", tag, "_CD45.png"))
}

# Tables
cells_clusters <- as.data.frame.matrix(table(obj$seurat_clusters, obj$orig.ident))
write.table(cbind(cluster = rownames(cells_clusters), cells_clusters), file.path(out_dir, "cells_per_cluster.tsv"),
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(Embeddings(obj, "umap"), file.path(out_dir, "umap_embeddings.tsv"), sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(obj@meta.data, file.path(out_dir, "clustered_metadata.tsv"), sep = "\t", row.names = TRUE, col.names = TRUE)

# Composition and QC by group
plot_composition(obj, group_by = "orig.ident", results_path = out_dir, filename = "ClusterComposition_bySample.png")
plot_vln_qc_by_group(obj, results_path = out_dir, filename = "VlnPlot_QCmetrics_bySample.png", group_by = "orig.ident")
plot_vln_qc_by_group(obj, results_path = out_dir, filename = "VlnPlot_QCmetrics_byCluster.png", group_by = "seurat_clusters")
if (has_subtype) {
  plot_composition(obj, group_by = "Subtype", results_path = out_dir, filename = "ClusterComposition_bySubtype.png")
  plot_vln_qc_by_group(obj, results_path = out_dir, filename = "VlnPlot_QCmetrics_bySubtype.png", group_by = "Subtype")

  # Pairwise Welch t-tests of nFeature_RNA between subtypes
  subtypes <- levels(factor(obj$Subtype))
  if (length(subtypes) >= 2) {
    pairs <- utils::combn(subtypes, 2, simplify = FALSE)
    pairwise_df <- do.call(rbind, lapply(pairs, function(p) {
      x <- obj$nFeature_RNA[obj$Subtype %in% p[1]]
      y <- obj$nFeature_RNA[obj$Subtype %in% p[2]]
      tt <- t.test(x, y)
      data.frame(group1 = p[1], group2 = p[2], mean_group1 = mean(x, na.rm = TRUE), mean_group2 = mean(y, na.rm = TRUE),
                 t_statistic = unname(tt$statistic), df = unname(tt$parameter), p_value = tt$p.value,
                 conf_low = tt$conf.int[1], conf_high = tt$conf.int[2])
    }))
    write.csv(pairwise_df, file.path(out_dir, "pairwisett_nfeaturerna_subtype.csv"), row.names = FALSE)
  }
}

saveRDS(obj, run_object("clustered"))
cat("\n Saved", run_object("clustered"), "\n")

# Markers of every cluster (Idents = global_clusters) + top-10 heatmap
markers_all <- FindAllMarkers(obj, only.pos = FALSE, logfc.threshold = -1)
write.table(markers_all, file.path(out_dir, "FindAllMarkers.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE)
plot_heatmap(obj, out_dir, markers_all, n = 10)

cat("\n ---- 03_cluster done. Next: 04_lineage_markers.R ----\n")
