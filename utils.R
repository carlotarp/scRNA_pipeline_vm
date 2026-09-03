##
## Shared utility functions for the scRNA-seq annotation pipeline.
## Source this file from: 3a_lineage_annotation, 3b_leukocytes, 3c_stroma, 3e_tumor
##

library(Seurat)
library(dplyr)

# ---------------------------------------------------------------
# Subset object to one lineage, normalize, recluster with Harmony.
# Returns the clustered subset — caller handles plotting.
#
# Args:
#   object       : Seurat object with a $lineage metadata column
#   lineage_name : value to filter on (e.g. "Leukocytes", "Stromal", "Tumor")
#   resolution   : FindClusters resolution
#   assay        : assay to split/normalize (default "RNA_decontX")
# ---------------------------------------------------------------
generate_lineage_subset <- function(object, lineage_name, resolution, assay = "RNA_decontX") {
  keep_cells <- colnames(object)[object$lineage == lineage_name]
  object_subset <- subset(object, cells = keep_cells)

  object_subset[[assay]] <- split(object_subset[[assay]], f = object_subset$orig.ident)
  object_subset <- NormalizeData(object_subset)
  object_subset <- FindVariableFeatures(object_subset)
  object_subset <- ScaleData(object_subset)
  object_subset <- RunPCA(object_subset)

  object_subset <- IntegrateLayers(
    object = object_subset, method = HarmonyIntegration,
    orig.reduction = "pca", new.reduction = "harmony", verbose = FALSE
  )

  pct   <- object_subset[["pca"]]@stdev / sum(object_subset[["pca"]]@stdev) * 100
  cumu  <- cumsum(pct)
  co1   <- which(cumu > 90 & pct < 5)[1]
  co2   <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = TRUE)[1] + 1
  co3   <- min(co1, co2)
  cat(paste0("  Optimal PCs (", lineage_name, "): ", co3, "\n"))

  object_subset <- FindNeighbors(object_subset, reduction = "harmony", dims = 1:co3)
  object_subset <- RunUMAP(object_subset, reduction = "harmony", dims = 1:co3)
  object_subset <- FindClusters(object_subset, resolution = resolution)
  object_subset <- JoinLayers(object_subset)

  return(object_subset)
}

# ---------------------------------------------------------------
# Run FindMarkers for a set of clusters and save one CSV per cluster.
#
# Args:
#   object       : Seurat object
#   clusters     : character/numeric vector of cluster identities to test
#   results_path : directory where CSVs are written
#   top_n        : number of top markers to keep (ranked by avg_log2FC)
# ---------------------------------------------------------------
find_markers_for_clusters <- function(object, clusters, results_path, top_n = 20) {
  clusters <- as.character(clusters)
  object_joined <- JoinLayers(object)

  for (cl in clusters) {
    markers_out <- FindMarkers(object_joined, ident.1 = cl, max.cells.per.ident = 5000)
    markers_out <- head(markers_out[order(-markers_out$avg_log2FC), ], top_n)
    write.csv(markers_out, file.path(results_path, paste0("cluster_", cl, "_FindMarkers.csv")))
    cat(paste0("  - Cluster ", cl, " -> cluster_", cl, "_FindMarkers.csv\n"))
  }
}
