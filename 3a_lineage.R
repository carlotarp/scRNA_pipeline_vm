##
##  Single Cell Analysis Step 3a: Lineage annotation
##  Runs AFTER 2a_cluster.R — takes clustered_data.rds as input.
##  Assigns broad lineages (Tumor / Leukocytes / Stromal) to clusters using
##  marker dotplots and FindMarkers CSVs, then exports lineage_annotated_data.rds.
##

# Import libraries
library("Seurat")
library(dplyr)

# Set paths
project_path <- "/home/usuario/PROJECTS/260724_victor_scRNA/"
wd <- paste0(project_path, "codes/scRNA_pipeline/")
setwd(wd)
results_path <- paste0(project_path, "results/")
results_GEMX_CA_path <- paste0(results_path, "GEMX/DecontX/CellAnnotation/")

# Import plot functions and shared utilities
source(paste0(wd, "CA_plots.R"))
source(paste0(wd, "utils.R"))

# Load clustered data
dwClustered <- readRDS(paste0(results_path, "GEMX/DecontX/Clustering/clustered_data.rds"))
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

clusters_to_check <- c("16", "19", "20") # DecontX 7500
find_markers_for_clusters(dwClustered, clusters_to_check, results_GEMX_CA_path)
cat("\n Read Ambiguous Cluster CSV (if needed) to complete 'cluster_to_lineage' \n")

# Manual cluster annotation
#lineage_clusters_annotated <- c( #########  3500  #########
#  "0" = "Tumor", "1" = "Tumor", "2" = "Stromal",
#  "3" = "Leukocytes", "4" = "Leukocytes", "5" = "Tumor",
#  "6" = "Tumor", "7" = "Tumor", "8" = "Stromal",
#  "9" = "Stromal", "10" = "Leukocytes", "11" = "Leukocytes",
#  "12" = "Stromal", "13" = "Stromal", "14" = "Tumor",
#  "15" = "Leukocytes", "16" = "Leukocytes", "17" = "Leukocytes",
#  "18" = "Stromal"
#)

#lineage_clusters_annotated <- c( #########  5500  #########
#  "0" = "Tumor", "1" = "Tumor", "2" = "Leukocytes",
#  "3" = "Stromal", "4" = "Leukocytes", "5" = "Tumor",
#  "6" = "Stromal", "7" = "Tumor", "8" = "Tumor",
#  "9" = "Stromal", "10" = "Leukocytes", "11" = "Tumor",
#  "12" = "Tumor", "13" = "Leukocytes", "14" = "Stromal",
#  "15" = "Leukocytes", "16" = "Leukocytes", "17" = "Tumor",
#  "18" = "Tumor"
#)

# lineage_clusters_annotated <- c( #########  7500  #########
#  "0" = "Tumor", "1" = "Tumor", "2" = "Leukocytes",
#  "3" = "Leukocytes", "4" = "Tumor", "5" = "Stromal",
#  "6" = "Stromal", "7" = "Tumor", "8" = "Tumor",
#  "9" = "Tumor", "10" = "Stromal", "11" = "Tumor",
#  "12" = "Stromal", "13" = "Leukocytes", "14" = "Leukocytes",
#  "15" = "Stromal", "16" = "Leukocytes", "17" = "Tumor",
#  "18" = "Leukocytes", "19" = "Stromal"
#)

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

dwClustered$lineage <- unname(lineage_clusters_annotated[as.character(dwClustered$decontX_clusters)])
dwClustered$lineage <- factor(dwClustered$lineage)

# ---  Lineage Annotated Dimplot ---
plot_dimplot(dwClustered, reduction = "umap_decontX", group_by = "lineage",
             results_path = results_GEMX_CA_path, filename = "DimPlot_UMAP_Lineage.png")

# Add celltype column (starts as lineage label, refined by lineage scripts)
dwClustered$celltype <- unname(lineage_clusters_annotated[as.character(dwClustered$decontX_clusters)])
dwClustered$celltype <- factor(dwClustered$celltype)

# Load contaminated baseline data (for comparison)
dwContaminated <- readRDS(file.path(results_path, "GEMX/CellAnnotation/7500/notumor_annotated_data.rds"))
dwClustered$clusters_cont <- dwContaminated$seurat_clusters
dwClustered$celltype_cont <- dwContaminated$celltype

# --- Dimplots comparing pre/post-decontX cluster and celltype distributions ---
plot_dimplot(dwClustered, reduction = "umap_decontX", group_by = "celltype_cont", label = T,
             results_path = results_GEMX_CA_path, filename = "UMAP_decontX_byContCelltype.png")
plot_dimplot(dwClustered, reduction = "umap_decontX", group_by = "clusters_cont", label = T,
             results_path = results_GEMX_CA_path, filename = "UMAP_decontX_byContClusters.png")
plot_dimplot(dwClustered, reduction = "umap_decontX", group_by = "decontX_clusters", label = T,
             results_path = results_GEMX_CA_path, filename = "UMAP_decontX_byDecontClusters.png")


# Export annotated data
saveRDS(dwClustered, file.path(results_GEMX_CA_path, "lineage_annotated_data.rds"))

cat(paste("\n ---- FINISHED LINEAGE ANNOTATION ----
    Run 3b_leukocytes.R, 3c_stroma.R, and 3d_tumor.R next (in any order).
    Generated files:
        · lineage_annotated_data.rds
        · cluster_(cluster)_FindMarkers.csv
    Generated plots:
        · DotPlot_(lineage).png
        · DimPlot_(groupedby).png
        "))


