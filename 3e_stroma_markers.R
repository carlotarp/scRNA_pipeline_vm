##
##  Single Cell Analysis Step 3e: Stroma marker plots
##  Runs AFTER 3d_leukocytes_annotate.R — takes leuko_annotated_data.rds as input.
##  Subsets and reclusters stromal cells, then produces the marker dotplots and
##  cluster CSVs needed to fill in the dictionary in 3f_stroma_annotate.R.
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
results_GEMX_SCA_path <- paste0(results_GEMX_CA_path, "Stroma/")

# Import plot functions and shared utilities
source(paste0(wd, "CA_plots.R"))
source(paste0(wd, "CL_plots.R"))
source(paste0(wd, "utils.R"))

# Load leukocyte annotated data
dwAnnotated <- readRDS(paste0(results_GEMX_CA_path, "leuko_annotated_data.rds"))
cat("\n leukocyte-annotated data loaded \n")

# Set stroma markers
markers_Stromal <- list(

  # --- Fibroblast / CAF subtypes ---
  CAF = list(
    cCAF = c("ACTA2", "TAGLN", "MYL9", "MYH11", "POSTN"),          # myofibroblastic, contractile
    iCAF = c("CXCL12", "CXCL14", "IL6", "PDGFRA", "CFD"),           # inflammatory, cytokine-secreting
    apCAF = c("CD74", "HLA-DRA", "HLA-DRB1"),                        # antigen-presenting CAF
    matrixCAF = c("COL1A1", "COL1A2", "COL3A1", "FN1", "LUM"),      # ECM-producing, general fibroblast
    vascularCAF = c("PECAM1", "RGS5", "NOTCH3")                      # perivascular-like CAF
  ),

  # --- Endothelial subtypes ---
  Endothelial = list(
    Vascular_general = c("PECAM1", "VWF", "CDH5"),
    Arterial = c("GJA5", "SEMA3G", "HEY1"),
    Venous = c("ACKR1", "SELP", "NR2F2"),
    Capillary = c("CA4", "RGCC"),
    Lymphatic = c("PROX1", "PDPN", "LYVE1", "CCL21"),
    Tip_cell = c("ESM1", "ANGPT2", "APLN")                           # active angiogenesis, vascular sprout tip
  ),

  # --- Mural cells (peri/vascular support) ---
  Mural = list(
    Pericyte = c("RGS5", "PDGFRB", "NOTCH3", "MCAM"),
    Smooth_muscle = c("MYH11", "ACTA2", "DES", "CNN1")
  ),

  # --- Adipocyte ---
  Adipocyte = list( Adipocyte = c("ADIPOQ", "PLIN1", "FABP4", "LEP"))
)

# Subset and recluster stromal cells
dwStroma <- generate_lineage_subset(dwAnnotated, lineage_name = "Stromal", resolution = 0.2)
plot_dimplot(dwStroma, reduction = "umap", group_by = "seurat_clusters",
             results_path = results_GEMX_SCA_path,
             filename = "DimPlot_UMAP_Stroma.png", label = TRUE)

# --- Dotplots by subtype ---
for (subtype in names(markers_Stromal)){
  plot_marker_dotplot(dwStroma, marker_groups = markers_Stromal[[subtype]],
                       results_path = results_GEMX_SCA_path,
                       filename = paste0("Dotplot_Stroma", subtype, ".png"),
                       group_by = "seurat_clusters")
}

# Markers of the clusters that the dotplots leave ambiguous
clusters_to_check <- c("5", "10") # 7500
find_markers_for_clusters(dwStroma, clusters_to_check, results_GEMX_SCA_path)

# Export the reclustered stromal subset
saveRDS(dwStroma, file.path(results_GEMX_SCA_path, "Stroma.rds"))

cat(paste("\n ---- FINISHED STROMA MARKER PLOTS ----
    Review the dotplots and cluster CSVs, fill in `clusters_Stroma_annotated`
    at the top of 3f_stroma_annotate.R, then run it.
    Generated files:
        · Stroma.rds  (reclustered, not yet annotated)
        · cluster_(cluster)_FindMarkers.csv
    Generated plots:
        · DimPlot_UMAP_Stroma.png
        · Dotplot_Stroma(subtype).png
        "))
