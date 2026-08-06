##
##  Single Cell Analysis Step 3a: Cell Annotation
##

# Import libraries
library("Seurat")
library(dplyr)

# Set paths
project_path <- "/home/usuario/PROJECTS/260724_victor_scRNA/"
wd <- paste0(project_path, "codes/scRNA_pipeline/")
setwd(wd)
results_path <- paste0(project_path, "results/")
results_GEMX_CA_path <- paste0(results_path, "GEMX/CellAnnotation/")

# Import Plot Functions
source(paste0(wd, "CA_plots.R"))

# Load Clustered Data
dwClustered <- readRDS(paste0(results_path, "GEMX/Clustering/clustered_data.rds"))
cat("\n Clustered Data Loaded \n")

# Set Markers
markers <- list(
  Leukocyte_pan = c("PTPRC"),
  Myeloid       = c("CD68", "CD163"),
  Tcell         = c("CD3E", "CD3D"),
  NK            = c("NKG7"),
  DC            = c("FCER1A", "IRF7", "IL3RA", "SPIB"),
  Neutrophil    = c("S100A9", "S100A8"),
  Bcell         = c("CD79A", "MS4A1"),
  plasma_Blast  = c("MZB1", "TNFRSF17", "XBP1")
)

# --- Dotplots & Featureplots ---
plot_marker_dotplot(dwClustered, marker_groups = markers, results_path = results_GEMX_CA_path,
                     filename = "DotPlot_Lineage.png", group_by = "seurat_clusters", title = "Lineage Markers by Cluster")

cat("\n Dot Plots generated \n")

# Find Markers of Ambiguous Clusters
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

clusters_to_check <- c("8", "12")
find_markers_for_clusters(dwClustered, clusters_to_check, results_GEMX_CA_path)
cat("\n Read Ambiguous Cluster CSV (if needed) to complete 'cluster_to_lineage' \n")


# Manual Cluster Annotation
clusters_annotated <- c(
  "0" = "Tumoral", "1" = "Stromal",     "2" = "Tumoral",
  "3" = "Leukocytes",     "4" = "Leukocytes",      "5" = "Tumoral",
  "6" = "Tumoral", "7" = "Leukocytes",      "8" = "Leukocytes",
  "9" = "Stromal", "10" = "Tumoral",
  "11" = "Leukocytes",     "12" = "Leukocytes"
)

dwClustered$lineage <- unname(clusters_annotated[as.character(dwClustered$seurat_clusters)])
dwClustered$lineage <- factor(dwClustered$lineage)

# ---  lineage Annotated Dimplot ---
plot_dimplot(dwClustered, reduction = "umap", group_by = "lineage",
             results_path = results_GEMX_CA_path, filename = "DimPlot_UMAP_Lineage.png")

# Export Annotated Data
saveRDS(dwClustered, file.path(results_GEMX_CA_path, "lineage_annotated_data.rds"))

cat(paste("\n ---- FINISHED CELL ANNOTATION ----
    Generated files:
      · lineage_annotated_data.rds
      · ambiguous_cluster_(cluster)_FindMarkers.csv
    Generated plots:
      · lineage_dotplot.png
      · (lineage)_featureplot.png
      · lineage_UMAP.png
          "))