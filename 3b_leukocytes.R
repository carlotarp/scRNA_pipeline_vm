##
##  Single Cell Analysis Step 3b: Leukocyte annotation
##  Runs AFTER 3a_lineage.R — takes lineage_annotated_data.rds as input.
##  Subsets and reclusters leukocytes, assigns fine-grained cell types using
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
results_GEMX_LCA_path <- paste0(results_GEMX_CA_path, "Leukocytes/")

# Import plot functions and shared utilities
source(paste0(wd, "CA_plots.R"))
source(paste0(wd, "utils.R"))

# Load lineage annotated data
dwAnnotated <- readRDS(paste0(results_GEMX_CA_path, "lineage_annotated_data.rds"))
cat("\n lineage-annotated data loaded \n")

# Set leukocyte markers
markers_leukocytes <- list(
  TcellNK = list(
    "Tcell_general"          = c("CD3E", "CD3D"),
    "CD4_Tcell"               = c("CD4"),
    "CD8_Tcell"               = c("CD8A", "CD8B"),
    "Tcell_activated"         = c("CCL5"),
    "Tcell_cytotoxic"         = c("GZMA", "GZMB", "PRF1", "GNLY"),
    "Tcell_mem_cytotoxic"     = c("TCF7"),
    "Tcell_antigen_reactive"  = c("ENTPD1"),
    "Exhaustion_markers"      = c("TIGIT", "CTLA4", "ICOS", "PDCD1", "EOMES"),
    "Treg"                    = c("CTLA4", "FOXP3"),
    "Tcell_naive"             = c("SELL", "IL7R"),
    "Tcell_mem"               = c("CD44"),
    "Tcell_em"                = c("CD69"),
    "Tcell_cm"                = c("CD27", "CCR7"),
    "Tcell_fh"                = c("CXCR5", "CXCL13"),
    "NK"                      = c("NKG7")
  ),
  Bcell = list(
    "Bcell_general"     = c("CD79A", "CD19", "MS4A1"),
    "Bcell_naive"       = c("IGHD"),
    "Bcell_activated"   = c("IGHM"),
    "plasma_Blast"      = c("IGHG", "TNFRSF17", "POU2AF1", "XBP1", "MZB1", "PIM2", "CD38", "IRF4", "PRDM1", "SDC1")
  ),
  IGG = list(
    "IGG" = c("CD27", "CD79A", "HLA-C", "JCHAIN", "IGKC", "IGLV3-25", "IL2RG", "CXCL8", "LAX1", "NTN3", "PIM2", "POU2AF1", "TNFRSF17")
  ),
  Myeloid = list(
    "TAM_general"       = c("CD163", "CD68"),
    "TAM_unpol"         = c("CCL7", "CCL18"),
    "TAM_proAngio_M1"   = c("CD86", "CD80", "IL1B", "TNF", "CXCL9", "CXCL10", "NOS2"),
    "TAM_proInflamm_M2" = c("CD163", "MRC1", "CD206", "MSR1", "CCL18", "CCL22", "IL10", "TGFB1", "APOE", "MARCO", "VSIG4"),
    "Monocyte"          = c("CD14", "CCR2", "CD64", "CD16", "LYZ")
  ),
  Dendritic = list(
    "Dendritic_general" = c("CD40", "CD83", "HLA-DRA", "ITGAX", "LYZ"),
    "pDC"               = c("CD53", "CLEC4C", "CLEC7A", "CORO1A", "FCER1G", "HLA-DRB1", "IL3RA", "NRP1", "IRF8", "JCHAIN", "IRF7"),
    "cDC"               = c("CADM1", "CD1C", "CLEC10A", "CLEC9A", "FCER1A", "FCER2B", "FLT3", "HLA-DPB1", "HLA-DQA1", "HLA-DQA2"),
    "actDC"             = c("LAMP3", "FSCN1", "IL12B")
  ),
  Mast = list(
    "Mast" = c("CLC", "CPA3", "GATA2", "HDC", "HPGDS", "IL1RL1", "IL5RA", "KIT", "LMO4", "MS4A2", "MS4A3", "PLIN2", "TNFSF10", "TPSAB1", "TPSB2", "SAMSN1", "CD69")
  ),
  TCell_Fibroblast_like = list(
    "Tcell_general"          = c("CD3E", "CD3D"),
    "CD4_Tcell"               = c("CD4"),
    "CD8_Tcell"               = c("CD8A", "CD8B"),
    "NK"                      = c("NKG7"),
    "Fibroblast" = c("COL1A1", "COL1A2"),
    "Treg" = c("FOXP3", "IL2RA", "CTLA4", "IKZF2", "TNFRSF18", "TNFRSF4"),
    "CAF" = c("THBS2", "FN1", "TAGLN", "MYL9", "TPM2", "CALD1", "CXCL12"),                                                                           # contractilCAF
    "Receptor" = c("CD3D", "CD3E", "CD3G", "CD247", "TRAC", "TRBC1", "TRBC2", "TRDC", "TRGC1", "TRGC2", "CD2", "CD5", "CD6")
  )
)

# Subset and recluster leukocytes
dwLeukocytes <- generate_lineage_subset(dwAnnotated, lineage_name = "Leukocytes", resolution = 0.4)
plot_dimplot(dwLeukocytes, reduction = "umap", group_by = "seurat_clusters",
             results_path = results_GEMX_LCA_path,
             filename = "DimPlot_UMAP_Leukocytes.png", label = TRUE)

# --- Dotplots by subtype ---
for (subtype in names(markers_leukocytes)){
  plot_marker_dotplot(dwLeukocytes, marker_groups = markers_leukocytes[[subtype]],
                       results_path = results_GEMX_LCA_path,
                       filename = paste0("Dotplot_Leukocytes", subtype, ".png"),
                       group_by = "seurat_clusters")
}

clusters_to_check <- c("10", "11") # DecontX 7500
find_markers_for_clusters(dwLeukocytes, clusters_to_check, results_GEMX_LCA_path)
cat("\n Read Ambiguous Cluster CSV (if needed) to complete the annotation \n")

# Manual cluster annotation
#clusters_leuko_annotated <- c( #########  3500  #########
#  "0" = "DC // TAM",  "1" = "TCell_naive",  "2" = "BCell",
#  "3" = "TCell_cyto",  "4" = "Fibrocyte",  "5" = "PlasmaBlast",
#  "6" = "TCell_ex",  "7" = "DC",  "8" = "Prolifetarive",
#  "9" = "Mast",  "10" = "pDC",  "11" = "PlasmaBlast",
#  "12" = "actDC",  "13" = "PlasmaBlast"
#)

#clusters_leuko_annotated <- c( #########  5500  #########
#  "0" = "pDC",  "1" = "Tcell",  "2" = "BCell",
#  "3" = "Tcell_naive",  "4" = "Fibrocyte",  "5" = "PlasmaBlast",
#  "6" = "Tcell_ex",  "7" = "pDC",  "8" = "Proliferative",
#  "9" = "Mast",  "10" = "TAM // Monocyte",  "11" = "pDC",
#  "12" = "PlasmaBlast",  "13" = "actDC" , "14" = "PlasmaBlast"
#)

#clusters_leuko_annotated <- c( #########  7500  #########
#  "0" = "pDC",  "1" = "Tcell_cyto",  "2" = "BCell",
#  "3" = "Tcell_naive",  "4" = "PlasmaBlast",  "5" = "Tcell_ex",
#  "6" = "Fibrocyte",  "7" = "pDC",  "8" = "Proliferative",
#  "9" = "Mast",  "10" = "pDC",  "11" = "actDC",
#  "12" = "PlasmaBlast"
#)

clusters_leuko_annotated <- c( #########  DecontX 7500  #########
  "0" = "TAM",  "1" = "TCell_cyto",  "2" = "BCell",
  "3" = "TCell_naive",  "4" = "PlasmaBlast",  "5" = "TCell_ex",
  "6" = "cDC",  "7" = "TAM_fibro",  "8" = "Mast",
  "9" = "pDC",  "10" = "Proliferative",  "11" = "TAM_proInflamm",
  "12" = "actDC", "13" = "PlasmaBlast"
)

dwLeukocytes$celltype <- unname(clusters_leuko_annotated[as.character(dwLeukocytes$seurat_clusters)])
dwLeukocytes$celltype <- factor(dwLeukocytes$celltype)

# ---  Visualize Annotated Dimplot ---
plot_dimplot(dwLeukocytes, reduction = "umap", group_by = "celltype", label = T,
             results_path = results_GEMX_LCA_path, filename = "DimPlot_UMAP_Leuko_Annotated.png")

# --- Dimplot comparing contaminated vs decontaminated annotation ---
plot_dimplot(dwLeukocytes, reduction = "umap", group_by = "celltype_cont", label = T,
             results_path = results_GEMX_LCA_path, filename = "DimPlot_UMAP_LeukoContaminated.png")

# Export annotated leukocyte data
saveRDS(dwLeukocytes, file.path(results_path, paste0("leukocytes.rds")))

# Label transfer to full object
annotation_vec <- setNames(as.character(dwLeukocytes[["celltype"]][, 1]),
                            colnames(dwLeukocytes))

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
             results_path = results_GEMX_CA_path, filename = "DimPlot_UMAP_Leuko.png")

# Export Annotated Data
saveRDS(dwAnnotated, file.path(results_GEMX_CA_path, "leuko_annotated_data.rds"))

cat(paste("\n ---- FINISHED LEUKOCYTE ANNOTATION ----
    Generated files:
        · leukocytes.rds
        · leuko_annotated_data.rds
    Generated plots:
        · DimPlot_UMAP_(groupedby).png
        · DotPlot_(subtype).png
        "))