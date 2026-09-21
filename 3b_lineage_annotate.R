##
##  Single Cell Analysis Step 3b: Lineage annotation
##  Runs AFTER 3a_lineage_markers.R — takes clustered_data.rds as input.
##  Applies the manual lineage dictionary (filled in from 3a's dotplots and
##  cluster marker CSVs) and exports lineage_annotated_data.rds.
##

# Import libraries
library("Seurat")
library(dplyr)

# Set paths
project_path <- "/home/user/PROJECTS/scRNA_vmendez/"
wd <- paste0(project_path, "codes/")
setwd(wd)
results_path <- paste0(project_path, "results/")
results_GEMX_CA_path <- paste0(results_path, "CellAnnotation/")

# Import plot functions and shared utilities
source(paste0(wd, "CA_plots.R"))
source(paste0(wd, "CL_plots.R"))
source(paste0(wd, "utils.R"))

# Manual cluster annotation, filled in from the 3a_lineage_markers.R outputs
lineage_clusters_annotated <- c( #########  DecontX 7500  #########
  "0" = "Tumor", "1" = "Stromal", "2" = "Tumor",
  "3" = "Leukocytes", "4" = "Leukocytes", "5" = "Tumor",
  "6" = "Tumor", "7" = "Stromal", "8" = "Tumor",
  "9" = "Leukocytes", "10" = "Leukocytes", "11" = "Stromal",
  "12" = "Stromal", "13" = "Tumor", "14" = "Leukocytes",
  "15" = "Leukocytes", "16" = "Leukocytes", "17" = "Stromal",
  "18" = "Tumor", "19" = "Tumor", "20" = "Tumor",
  "21" = "Tumor", "22" = "Tumor"
)

# Load clustered data
dwClustered <- readRDS(paste0(results_path, "Clustering/clustered_data.rds"))
dwClustered[["RNA_decontX"]] <- as(dwClustered[["RNA_decontX"]], "Assay5")
cat("\n Clustered Data Loaded \n")

dwClustered$lineage <- unname(lineage_clusters_annotated[as.character(dwClustered$decontX_clusters)])
dwClustered$lineage <- factor(dwClustered$lineage)

# ---  Lineage Annotated Dimplot ---
plot_dimplot(dwClustered, reduction = "umap_decontX", group_by = "lineage",
             results_path = results_GEMX_CA_path, filename = "DimPlot_UMAP_Lineage.png")

# Add celltype column (starts as lineage label, refined by lineage scripts)
dwClustered$celltype <- unname(lineage_clusters_annotated[as.character(dwClustered$decontX_clusters)])
dwClustered$celltype <- factor(dwClustered$celltype)

# Uncorrected (pre-DecontX) baseline clustering, computed in 2b_cluster.R
dwClustered$clusters_cont <- dwClustered$contaminated_clusters

# --- Dimplots comparing pre/post-decontX cluster distributions ---
plot_dimplot(dwClustered, reduction = "umap_decontX", group_by = "clusters_cont", label = T,
             results_path = results_GEMX_CA_path, filename = "UMAP_decontX_byContClusters.png")
plot_dimplot(dwClustered, reduction = "umap_decontX", group_by = "decontX_clusters", label = T,
             results_path = results_GEMX_CA_path, filename = "UMAP_decontX_byDecontClusters.png")

# Export annotated data
saveRDS(dwClustered, file.path(results_GEMX_CA_path, "lineage_annotated_data.rds"))

cat(paste("\n ---- FINISHED LINEAGE ANNOTATION ----
    Run 3c_leukocytes_markers.R next.
    Generated files:
        · lineage_annotated_data.rds
    Generated plots:
        · DimPlot_UMAP_Lineage.png
        · UMAP_decontX_by(Cont|Decont)(Celltype|Clusters).png
        "))
