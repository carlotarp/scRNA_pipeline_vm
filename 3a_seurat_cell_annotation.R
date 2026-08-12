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
results_GEMX_CA_path <- paste0(results_path, "GEMX/CellAnnotation/3500/")

# Import Plot Functions
source(paste0(wd, "CA_plots.R"))

# Load Clustered Data
dwClustered <- readRDS(paste0(results_path, "GEMX/Clustering/3500/clustered_data.rds"))
cat("\n Clustered Data Loaded \n")

# Set Markers
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

clusters_to_check <- c("1", "6", "15", "16") # 3500
clusters_to_check <- c("7", "15") # 5500
clusters_to_check <- c("1", "8", "16") # 7500

find_markers_for_clusters(dwClustered, clusters_to_check, results_GEMX_CA_path)
cat("\n Read Ambiguous Cluster CSV (if needed) to complete 'cluster_to_lineage' \n")


# Manual Cluster Annotation
lineage_clusters_annotated <- c( #########  3500  #########
  "0" = "Tumor", "1" = "Tumor", "2" = "Stromal",
  "3" = "Leukocytes", "4" = "Leukocytes", "5" = "Tumor",
  "6" = "Tumor", "7" = "Tumor", "8" = "Stromal",
  "9" = "Stromal", "10" = "Leukocytes", "11" = "Leukocytes",
  "12" = "Stromal", "13" = "Stromal", "14" = "Tumor",
  "15" = "Leukocytes", "16" = "Leukocytes", "17" = "Leukocytes",
  "18" = "Stromal"
)

lineage_clusters_annotated <- c( #########  5500  #########
  "0" = "Tumor", "1" = "Tumor", "2" = "Leukocytes",
  "3" = "Stromal", "4" = "Leukocytes", "5" = "Tumor",
  "6" = "Stromal", "7" = "Tumor", "8" = "Tumor",
  "9" = "Stromal", "10" = "Leukocytes", "11" = "Tumor",
  "12" = "Tumor", "13" = "Leukocytes", "14" = "Stromal",
  "15" = "Leukocytes", "16" = "Leukocytes", "17" = "Tumor",
  "18" = "Tumor"
)

lineage_clusters_annotated <- c( #########  7500  #########
  "0" = "Tumor", "1" = "Tumor", "2" = "Leukocytes",
  "3" = "Leukocytes", "4" = "Tumor", "5" = "Stromal",
  "6" = "Stromal", "7" = "Tumor", "8" = "Tumor",
  "9" = "Tumor", "10" = "Stromal", "11" = "Tumor",
  "12" = "Stromal", "13" = "Leukocytes", "14" = "Leukocytes",
  "15" = "Stromal", "16" = "Leukocytes", "17" = "Tumor",
  "18" = "Leukocytes", "19" = "Stromal"
)

dwClustered$lineage <- unname(lineage_clusters_annotated[as.character(dwClustered$seurat_clusters)])
dwClustered$lineage <- factor(dwClustered$lineage)

# ---  lineage Annotated Dimplot ---
plot_dimplot(dwClustered, reduction = "umap", group_by = "lineage",
             results_path = results_GEMX_CA_path, filename = "DimPlot_UMAP_Lineage.png")

subtype_clusters_annotated <- c( #########  3500  #########
  "0" = "Tumor", "1" = "Tumor", "2" = "Fibroblast // CAF",
  "3" = "TCell // NK", "4" = "Myeloid", "5" = "Tumor",
  "6" = "Tumor", "7" = "Tumor", "8" = "Fibroblast // CAF",
  "9" = "Endothelial", "10" = "BCell", "11" = "PlasmaBlast",
  "12" = "cCAF", "13" = "Fibroblast // CAF // Pericyte", "14" = "Tumor",
  "15" = "Mast", "16" = "Proliferative", "17" = "DC",
  "18" = "Adipocyte"
)

subtype_clusters_annotated <- c( #########  5500  #########
  "0" = "Tumor", "1" = "Tumor", "2" = "TCell",
  "3" = "Fibroblast // CAF // Pericyte", "4" = "Myeloid", "5" = "Tumor",
  "6" = "Fibroblast // CAF", "7" = "Tumor", "8" = "Tumor",
  "9" = "Endothelial", "10" = "BCell", "11" = "Tumor",
  "12" = "Tumor", "13" = "PlasmaBlast", "14" = "cCAF",
  "15" = "Proliferative", "16" = "Mast", "17" = "Tumor",
  "18" = "Tumor"
)

subtype_clusters_annotated <- c( #########  7500  #########
  "0" = "Tumor", "1" = "Tumor", "2" = "TCell // NK",
  "3" = "Myeloid", "4" = "Tumor", "5" = "Fibroblast // CAF",
  "6" = "Fibroblast // CAF", "7" = "Tumor", "8" = "Tumor",
  "9" = "Tumor", "10" = "Fibroblast // CAF", "11" = "Tumor",
  "12" = "Endothelial", "13" = "BCell", "14" = "PlasmaBlast",
  "15" = "cCAF", "16" = "Proliferative", "17" = "Tumor",
  "18" = "Mast", "19" = "Fibroblast // CAF // Pericyte"
)

dwClustered$celltype <- unname(subtype_clusters_annotated[as.character(dwClustered$seurat_clusters)])
dwClustered$celltype <- factor(dwClustered$celltype)

# ---  Stromal Annotated Dimplot ---
plot_dimplot(dwClustered, reduction = "umap", group_by = "celltype", label = T,
             results_path = results_GEMX_CA_path, filename = "DimPlot_UMAP_NoTumor1.png")

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


