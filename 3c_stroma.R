##
##  Single Cell Analysis Step 3c: Stroma annotation
##  Runs AFTER 3b_leukocytes.R — takes the partially annotated object as input.
##  Subsets and reclusters stromal cells, assigns fine-grained cell types using
##  marker dotplots, then transfers labels back to the full object.
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
results_GEMX_SCA_path <- paste0(results_GEMX_CA_path, "Stroma/")

# Import plot functions and shared utilities
source(paste0(wd, "CA_plots.R"))
source(paste0(wd, "utils.R"))

# Load lineage annotated data
dwAnnotated <- readRDS(paste0(results_GEMX_CA_path, "fully_annotated_data.rds"))
cat("\n lineage-annotated data loaded \n")

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

clusters_to_check <- c("5", "10") # 7500
find_markers_for_clusters(dwStroma, clusters_to_check, results_GEMX_SCA_path)
cat("\n Read Ambiguous Cluster CSV (if needed) to complete the annotation \n")

# Manual cluster annotation
# clusters_Stroma_annotated <- c( #########  3500  #########
#  "0" = "DC // TAM",  "1" = "TCell_naive",  "2" = "BCell",
#  "3" = "TCell_cyto",  "4" = "Fibrocyte",  "5" = "PlasmaBlast",
#  "6" = "TCell_ex",  "7" = "DC",  "8" = "Prolifetarive",
#  "9" = "Mast",  "10" = "pDC",  "11" = "PlasmaBlast",
#  "12" = "actDC",  "13" = "PlasmaBlast"
#)

#clusters_Stroma_annotated <- c( #########  5500  #########
#  "0" = "pDC",  "1" = "Tcell",  "2" = "BCell",
#  "3" = "Tcell_naive",  "4" = "Fibrocyte",  "5" = "PlasmaBlast",
#  "6" = "Tcell_ex",  "7" = "pDC",  "8" = "Proliferative",
#  "9" = "Mast",  "10" = "TAM // Monocyte",  "11" = "pDC",
#  "12" = "PlasmaBlast",  "13" = "actDC" , "14" = "PlasmaBlast"
#)

#clusters_Stroma_annotated <- c( #########  7500  #########
#  "0" = "matrixCAF",  "1" = "apCAF",  "2" = "Endothelial",
#  "3" = "cCAF",  "4" = "iCAF",  "5" = "Pericyte",
#  "6" = "Adipocyte"
#)

clusters_Stroma_annotated <- c( #########  DecontX 7500  #########
  "0" = "mCAF",  "1" = "Endothelial",  "2" = "mCAF",
  "3" = "cCAF",  "4" = "iCAF",  "5" = "Pericyte",
  "6" = "Adipocyte", "7" = "Lymphatic"
)

dwStroma$celltype <- unname(clusters_Stroma_annotated[as.character(dwStroma$seurat_clusters)])
dwStroma$celltype <- factor(dwStroma$celltype)

# ---  Visualize Annotated Dimplot ---
plot_dimplot(dwStroma, reduction = "umap", group_by = "celltype", label = T,
             results_path = results_GEMX_SCA_path, filename = "DimPlot_UMAP_Stroma_Annotated.png")

# --- Dimplot comparing contaminated vs decontaminated annotation ---
plot_dimplot(dwStroma, reduction = "umap", group_by = "celltype_cont", label = T,
             results_path = results_GEMX_SCA_path, filename = "DimPlot_UMAP_StromaContaminated.png")

# Export annotated stroma data
saveRDS(dwStroma, file.path(results_path, paste0("Stroma.rds")))

# Label transfer to full object
annotation_vec <- setNames(as.character(dwStroma[["celltype"]][, 1]),
                            colnames(dwStroma))

if ("celltype" %in% colnames(dwAnnotated@meta.data)) {
  existing <- as.character(dwAnnotated[["celltype"]][, 1])
} else {
  existing <- rep(NA_character_, ncol(dwAnnotated))
}
names(existing) <- colnames(dwAnnotated)
existing[names(annotation_vec)] <- annotation_vec
dwAnnotated[["celltype"]] <- factor(existing)

# ---  Visualize Annotated Dimplot ---
plot_dimplot(dwAnnotated, reduction = "umap_decontX", group_by = "celltype", label = T,
             results_path = results_GEMX_CA_path, filename = "DimPlot_UMAP_Stroma.png")

# Export Annotated Data
saveRDS(dwAnnotated, file.path(results_GEMX_CA_path, "notumor_annotated_data.rds"))

cat(paste("\n ---- FINISHED STROMA ANNOTATION ----
    Generated files:
        · Stroma.rds
        · notumor_annotated_data.rds
    Generated plots:
        · DimPlot_UMAP_(groupedby).png
        · DotPlot_(subtype).png
        "))