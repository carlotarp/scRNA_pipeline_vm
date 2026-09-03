##
##  Single Cell Analysis Step 2b: Cluster Marker Genes
##  Runs AFTER 2_seurat_clustering.R — takes clustered_data.rds as input.
##  Separated because FindAllMarkers is slow; run once and keep the output.
##

# Import libraries
library("Seurat")
library(dplyr)

# Set paths
project_path <- "/home/usuario/PROJECTS/260724_victor_scRNA/"
wd <- paste0(project_path, "codes/scRNA_pipeline/")
setwd(wd)
results_path <- paste0(project_path, "results/")
results_GEMX_CL_path <- paste0(results_path, "GEMX/DecontX/Clustering/")

# Import plot functions
source(paste0(wd, "CL_plots.R"))

# Load clustered data
dwIntegrated <- readRDS(file.path(results_GEMX_CL_path, "clustered_data.rds"))
cat(" Clustered data loaded \n")

# Find marker genes per cluster
tumor_markers_all <- FindAllMarkers(dwIntegrated, only.pos = FALSE, logfc.threshold = -1)
write.table(tumor_markers_all, file = paste0(results_GEMX_CL_path, "all_tumor_markers.tsv"),
            sep = "\t", row.names = FALSE, col.names = TRUE)
cat("\n Differential expression per cluster done \n")

# --- Heatmap of top marker genes per cluster ---
plot_heatmap(dwIntegrated, results_GEMX_CL_path, tumor_markers_all, n = 10)

cat(paste("\n ---- FINISHED CLUSTER MARKER GENES ----
    Generated files:
      · all_tumor_markers.tsv
    Generated plots:
      · Integrated_Top10gene_Heatmap.png
          "))
