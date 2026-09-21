##
##  Single Cell Analysis Step 3a: Lineage marker plots
##  Runs AFTER 2c_markers.R — takes clustered_data.rds as input.
##  Produces the lineage marker dotplots and per-cluster FindMarkers CSVs
##  needed to fill in the lineage dictionary in 3b_lineage_annotate.R.
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

# Load clustered data
dwClustered <- readRDS(paste0(results_path, "Clustering/clustered_data.rds"))
dwClustered[["RNA_decontX"]] <- as(dwClustered[["RNA_decontX"]], "Assay5")
cat("\n Clustered Data Loaded \n")

# Set markers
markers_leukos <- list(
  Leukocyte_pan = c("PTPRC"),
  Myeloid       = c("CD68", "CD163"),
  Tcell         = c("CD3E", "CD3D"),
  NK            = c("NKG7"),
  DC            = c("FCER1A", "IRF7", "IL3RA", "SPIB"),
  Neutrophil    = c("S100A9", "S100A8"),
  Bcell         = c("CD79A", "MS4A1"),
  plasma_Blast  = c("MZB1", "TNFRSF17", "XBP1"),
  Mast          = c("CPA3", "GATA2")
)

markers_stromal <- list(
  Fibroblast_general = c("PDGFRA", "DCN", "LUM", "COL1A1", "COL1A2", "DPT", "CD34", "CXCL14", "FBLN1", "MFAP5", "APOD"),
  matrixCAF           = c("FAP", "COL11A1", "POSTN", "CTHRC1", "ASPN", "COMP", "COL10A1", "INHBA", "TNC", "MMP11", "LOXL2", "LRRC15", "THBS2", "FN1", "COL5A2", "COL8A1"),
  contractilCAF        = c("ACTA2", "TAGLN", "MYL9", "TPM2", "CNN", "MYH11", "DES", "CALD1", "COL1A1"),
  Adipocyte            = c("FASN", "GPAM", "LEP", "EBF1", "PDE3B", "PPARG", "CD36")
)

markers_stromaangio <- list(
  Pericyte    = c("RGS5", "CSPG4", "MCAM", "PDGFRB", "NOTCH3", "KCNJ8", "ABCC9", "DES", "CD248", "ANPEP"),
  Endothelial = c("VWF", "EGFL7", "FLT1", "EMCN", "PTPRB", "ENG", "CALCRL", "EPAS1", "ADGRL4", "CLDN5", "CDH5", "CD34", "DLL4", "ACKR1", "PECAM1")
)

pam50_genes <- list(
  Lum   = c("ESR1", "PGR", "BAG1", "MAPT", "NAT1", "ZIP6"),
  Basal = c("MKI67", "CCNE1", "ANLN", "CDC20", "EGFR", "MYC"),
  Her2  = c("ERBB2", "GRB7", "BLVRA", "TMEM45B")
)

# --- Dotplots & Featureplots ---
plot_marker_dotplot(dwClustered, marker_groups = markers_leukos, results_path = results_GEMX_CA_path,
                     filename = "DotPlot_Leukocytes.png", group_by = "seurat_clusters", title = "Leukocyte Markers by Cluster")

plot_marker_dotplot(dwClustered, marker_groups = markers_stromal, results_path = results_GEMX_CA_path,
                     filename = "DotPlot_Stromal1.png", group_by = "seurat_clusters", title = "Stromal Markers by Cluster")

plot_marker_dotplot(dwClustered, marker_groups = markers_stromaangio, results_path = results_GEMX_CA_path,
                     filename = "DotPlot_Stromal2.png", group_by = "seurat_clusters", title = "Stromal Markers by Cluster")

plot_marker_dotplot(dwClustered, marker_groups = pam50_genes, results_path = results_GEMX_CA_path,
                     filename = "DotPlot_Tumor.png", group_by = "seurat_clusters", title = "Tumor Markers by Cluster")

cat("\n Dot Plots generated \n")

# Markers of the clusters that the dotplots leave ambiguous
clusters_to_check <- c("16", "19", "20") # DecontX 7500
find_markers_for_clusters(dwClustered, clusters_to_check, results_GEMX_CA_path)

cat(paste("\n ---- FINISHED LINEAGE MARKER PLOTS ----
    Review the dotplots and cluster CSVs, fill in `lineage_clusters_annotated`
    at the top of 3b_lineage_annotate.R, then run it.
    Generated files:
        · cluster_(cluster)_FindMarkers.csv
        · all_clusters_FindMarkers.csv
    Generated plots:
        · DotPlot_(lineage).png
        "))
