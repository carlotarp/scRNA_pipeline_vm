##
##  Single Cell Analysis Step 3c: Leukocyte marker plots
##  Runs AFTER 3b_lineage_annotate.R — takes lineage_annotated_data.rds as input.
##  Subsets and reclusters leukocytes, then produces the marker dotplots and
##  cluster CSVs needed to fill in the dictionary in 3d_leukocytes_annotate.R.
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
results_GEMX_LCA_path <- paste0(results_GEMX_CA_path, "Leukocytes/")

# Import plot functions and shared utilities
source(paste0(wd, "CA_plots.R"))
source(paste0(wd, "CL_plots.R"))
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

# Markers of the clusters that the dotplots leave ambiguous
clusters_to_check <- c("10", "11") # DecontX 7500
find_markers_for_clusters(dwLeukocytes, clusters_to_check, results_GEMX_LCA_path)

# Export the reclustered leukocyte subset
saveRDS(dwLeukocytes, file.path(results_GEMX_LCA_path, "leukocytes.rds"))

cat(paste("\n ---- FINISHED LEUKOCYTE MARKER PLOTS ----
    Review the dotplots and cluster CSVs, fill in `clusters_leuko_annotated`
    at the top of 3d_leukocytes_annotate.R, then run it.
    Generated files:
        · leukocytes.rds  (reclustered, not yet annotated)
        · cluster_(cluster)_FindMarkers.csv
    Generated plots:
        · DimPlot_UMAP_Leukocytes.png
        · Dotplot_Leukocytes(subtype).png
        "))
