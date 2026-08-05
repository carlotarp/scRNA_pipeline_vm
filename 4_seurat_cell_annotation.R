##
##  Single Cell Analysis Step 3: Cell Annotation
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
  Tumoral = c("EPCAM", "KRT8", "KRT18", "KRT19", "CLDN4",      # luminal
                 "KRT5", "KRT14", "TP63"),                        # basal/myoTumoral
  Leukocytes     = c("PTPRC",                                         # pan-Leukocytes
                 "CD68", "CD163",                                 # Myeloid
                 "CD3E", "CD3D",                                  # Tcell
                 "NKG7",                                          # NK
                 "FCER1A",                                        # DC
                 "IRF7", "IL3RA", "SPIB",                         # DC (more immature)
                 "S100A9", "S100A8",                              # Neutrophil
                 "CD79A", "MS4A1",                                # Bcell
                 "MZB1", "TNFRSF17", "XBP1"),                     # plasma Blast
  Stroma     = c("PDGFRA", "DCN", "LUM", "COL1A1",                # fibroblast/CAF
                 "RGS5", "MCAM",                                  # pericyte
                 "VWF", "PECAM1")                                 # endothelial
)

# --- Dotplots & Featureplots ---
plot_marker_dotplot(dwClustered, marker_groups = markers, results_path = results_GEMX_CA_path,
                     filename = "DotPlot_Lineage.png", group_by = "seurat_clusters", title = "Lineage Markers by Cluster")

for (lineage in names(markers)) {
  plot_featureplot(dwClustered, reduction = "umap", features = markers[[lineage]],
                    results_path = results_GEMX_CA_path,
                    filename = paste0("", lineage, "_featureplot.png"))
}
cat("\n lineage plots generated \n")

# Diagnosis of Ambiguous Clusters
diagnose_ambiguous_clusters <- function(object, markers, results_path,
                                          z_threshold = 1, top_n = 20) {
  markers_flat <- unique(unlist(markers))
  markers_present <- markers_flat[markers_flat %in% rownames(object)]
  missing <- setdiff(markers_flat, markers_present)
  if (length(missing) > 0) {
    cat(paste("\n Missing Markers:",
              paste(missing, collapse = ", "), "\n"))
  }

  avg_exp <- AverageExpression(object, features = markers_present,
                                group.by = "seurat_clusters", assay = "RNA")$RNA

  # z-score per gene
  avg_exp_z <- t(scale(t(avg_exp)))
  max_z_per_cluster <- apply(avg_exp_z, 2, max, na.rm = TRUE)

  ambiguous_valid <- names(max_z_per_cluster)[max_z_per_cluster < z_threshold]
  ambiguous <- sub("^g", "", ambiguous_valid)

  cat(paste("\n", length(ambiguous), " ambiguous clusters detected",
            paste(ambiguous, collapse = ", "), "\n"))

  if (length(ambiguous) != 0) {
    cat(" Running FindMarkers... \n")
    object_joined <- JoinLayers(object)
    for (cl in ambiguous) {
      markers_out <- FindMarkers(object_joined, ident.1 = cl)
      markers_out <- head(markers_out[order(-markers_out$avg_log2FC), ], top_n)
      write.csv(markers_out, file.path(results_path, paste0("ambiguous_cluster_", cl, "_FindMarkers.csv")))
      cat(paste0("  · Cluster ", cl, " -> ambiguous_cluster_", cl, "_FindMarkers.csv\n"))
    }
  }
}

diagnose_ambiguous_clusters(dwClustered, markers, results_GEMX_CA_path)
cat("\n Read Ambiguous Cluster CSV (if needed) to complete 'cluster_to_lineage' \n")


# Manual Cluster Annotation
cluster_to_lineage <- c(
  "0" = "Tumoral", "1" = "Stromal",     "2" = "Tumoral",
  "3" = "Leukocytes",     "4" = "Leukocytes",      "5" = "Tumoral",
  "6" = "Tumoral", "7" = "Leukocytes",      "8" = "Leukocytes",
  "9" = "Stromal", "10" = "Tumoral",
  "11" = "Leukocytes",     "12" = "Leukocytes"
)

dwClustered$lineage <- unname(cluster_to_lineage[as.character(dwClustered$seurat_clusters)])
dwClustered$lineage <- factor(dwClustered$lineage)

# ---  lineage Annotated Dimplot ---
plot_dimplot(dwClustered, reduction = "umap", group_by = "lineage",
             results_path = results_GEMX_CA_path, filename = "lineage_UMAP.png")

# Export lineage Annotated Dara
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