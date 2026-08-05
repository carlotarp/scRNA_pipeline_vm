##
##  Single Cell Analysis Step 4: Cell Subtype Annotation
##

# Import libraries
library("Seurat")
library(dplyr)

# Set paths
project_path <- "/home/usuario/PROJECTS/260724_victor_scRNA/"
wd <- paste0(project_path, "codes/scRNA_pipeline/")
setwd(wd)
results_path <- paste0(project_path, "results/")
results_GEMX_CSA_path <- paste0(results_path, "GEMX/CellSubtypeAnnotation/")

# Import Plot Functions
source(paste0(wd, "CA_plots.R"))

# Load lineage Annotated Data
dwAnnotated <- readRDS(paste0(results_GEMX_CA_path, "lineage_annotated_data.rds"))
cat("\n lineage-annotated data loaded \n")


# Set Subtype Markers

markers_subtypes <- list(
  Leukocytes = list(
    "Myeloid_general"      = c("CD163", "CD68"),
    "TAM_unpolarized"      = c("CCL7", "CCL18"),
    "TAM_proAngio_M1"      = c("CD86", "CD80", "IL1B", "TNF", "CXCL9", "CXCL10", "NOS2"),
    "TAM_proInflamm_M2"    = c("CD163", "MRC1", "MSR1", "CCL18", "CCL22", "IL10", "TGFB1", "APOE", "MARCO", "VSIG4"),
    "Monocyte"             = c("CD14", "CCR2", "CD64", "CD16", "LYZ"),
    "Neutrophil"           = c("S100A9", "S100A8"),
    "Tcell_general"        = c("CD3E", "CD3D"),
    "CD4_Tcell"            = c("CD4"),
    "CD8_Tcell"            = c("CD8A", "CD8B"),
    "Tcell_activated"      = c("CCL5"),
    "Tcell_cytotoxic"      = c("GZMA", "GZMB", "PRF1", "GNLY"),
    "Tcell_mem_cytotoxic"  = c("TCF7"),
    "Tcell_antigen_reactive" = c("ENTPD1"),
    "Tcell_exhausted"      = c("TIGIT", "CTLA4", "ICOS", "PDCD1", "EOMES"),
    "Treg"                 = c("CTLA4", "FOXP3"),
    "Tcell_naive"          = c("SELL", "IL7R"),
    "Tcell_mem"            = c("CD44"),
    "Tcell_em"             = c("CD69"),
    "Tcell_cm"             = c("CD27", "CCR7"),
    "Tcell_fh"             = c("CXCR5", "CXCL13"),
    "NK"                   = c("NKG7"),
    "Bcell_general"        = c("CD79A", "CD19", "MS4A1"),
    "Bcell_naive"          = c("IGHD"),
    "Bcell_activated"      = c("IGHM"),
    "plasma_Blast"         = c("IGHG", "TNFRSF17", "POU2AF1", "XBP1", "MZB1", "PIM2", "CD38", "IRF4", "PRDM1", "SDC1"),
    "IGG"                  = c("CD27", "CD79A", "HLA-C", "JCHAIN", "IGKC", "IGLV3-25", "IL2RG", "CXCL8", "LAX1", "NTN3", "PIM2", "POU2AF1", "TNFRSF17"),
    "Dendritic_general"    = c("CD40", "CD83", "HLA-DRA", "ITGAX", "LYZ"),
    "pDC"                  = c("CD53", "CLEC4C", "CLEC7A", "CORO1A", "FCER1G", "HLA-DRB1", "IL3RA", "NRP1", "IRF8", "JCHAIN", "IRF7"),
    "cDC"                  = c("CADM1", "CD1C", "CLEC10A", "CLEC9A", "FCER1A", "FCER2B", "FLT3", "HLA-DPB1", "HLA-DQA1", "HLA-DQA2"),
    "actDC"                = c("LAMP3", "FSCN1", "IL12B"),
    "Mast"                 = c("CLC", "CPA3", "GATA2", "HDC", "HPGDS", "IL1RL1", "IL5RA", "KIT", "LMO4", "MS4A2", "MS4A3", "PLIN2", "TNFSF10", "TPSAB1", "TPSB2", "SAMSN1", "CD69")
  ),
  Stromal = list(
    "Fibroblast_general" = c("PDGFRA", "DCN", "LUM", "COL1A1", "COL1A2", "DPT", "CD34", "CXCL14", "FBLN1", "MFAP5", "APOD"),
    "matrixCAF"          = c("FAP", "COL11A1", "POSTN", "CTHRC1", "ASPN", "COMP", "COL10A1", "INHBA", "TNC", "MMP11", "LOXL2", "LRRC15", "THBS2", "FN1", "COL5A2", "COL8A1"),
    "contractilCAF"      = c("ACTA2", "TAGLN", "MYL9", "TPM2", "CNN", "MYH11", "DES", "CALD1", "COL1A1"),
    "Adipocyte"          = c("FASN", "GPAM", "LEP", "EBF1", "PDE3B", "PPARG", "CD36"),
    "Pericyte"           = c("RGS5", "CSPG4", "MCAM", "PDGFRB", "NOTCH3", "KCNJ8", "ABCC9", "DES", "CD248", "ANPEP"),
    "Endothelial"        = c("VWF", "EGFL7", "FLT1", "EMCN", "PTPRB", "ENG", "CALCRL", "EPAS1", "ADGRL4", "CLDN5", "CDH5", "CD34", "DLL4", "ACKR1", "PECAM1")
  )
)


# Generate Subsets
generate_subset <- function(object, lineage_name,
                            results_path) {
  object_subset <- subset(object, subset = lineage == lineage_name)
  object_subset <- NormalizeData(object_subset)
  object_subset <- FindVariableFeatures(object_subset)
  object_subset <- ScaleData(object_subset)
  object_subset <- RunPCA(object_subset)

  object_subset <- IntegrateLayers(
    object = object_subset, method = HarmonyIntegration,
    orig.reduction = "pca", new.reduction = "harmony", verbose = FALSE
  )

  pct <- object_subset[["pca"]]@stdev / sum(object_subset[["pca"]]@stdev) * 100
  cumu <- cumsum(pct)
  co1 <- which(cumu > 90 & pct < 5)[1]
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = TRUE)[1] + 1
  co3_subset <- min(co1, co2)
  cat(paste0("  co3 Computed for ", lineage_name, ": ", co3_subset, "\n"))

  object_subset <- FindNeighbors(object_subset, reduction = "harmony", dims = 1:co3_subset)
  object_subset <- RunUMAP(object_subset, reduction = "harmony", dims = 1:co3_subset)

  return(object_subset)
}

# Final Clustering w/ Optimal Resultion
annotate_lineage_subtypes <- function(object, lineage_name,
                                      subtype_list, results_path, resolution) {
  subset_obj <- FindClusters(object, resolution = resolution)
  subset_obj <- JoinLayers(subset_obj)

  markers_all <- unique(unlist(subtype_list))

  # --- Dotplot by Subtype ---
  plot_marker_dotplot(subset_obj, marker_groups = subtype_list,
                       results_path = results_path,
                       filename = paste0("Dotplot_", lineage_name, ".png"),
                       group_by = "seurat_clusters")

# --- Featureplot w/ All Lineage Markers ---
  plot_featureplot(subset_obj, reduction = "umap", features = markers_all,
                    results_path = results_path,
                    filename = paste0("FeaturePlot_", lineage_name, ".png"))

# --- Visualize UMAP ---
  plot_dimplot(subset_obj, reduction = "umap", group_by = "seurat_clusters",
               results_path = results_path,
               filename = paste0("DimPlot_UMAP_", lineage_name, ".png"), label = TRUE)

  saveRDS(subset_obj, file.path(results_path, paste0(lineage_name, "_subset.rds")))

  return(subset_obj)
}

# Execute Functions by Lineage
dwLeukocytes <- generate_subset(
    object = dwAnnotated, lineage_name = "Leukocytes",
    results_path = results_GEMX_CSA_path
  )

dwStromal <- generate_subset(
    object = dwAnnotated, lineage_name = "Stromal",
    results_path = results_GEMX_CSA_path
  )

cat("\n Check ResolutionGrid_(lineage).png and choose a resolution \n")

# Set Resolution per Lineage
res_leukocytes <- 0.6
res_stromal <- 0.6

annotate_lineage_subtypes(
    object = dwLeukocytes,
    lineage_name = "Leukocytes",
    subtype_list = markers_subtypes[["Leukocytes"]],
    results_path = results_GEMX_CSA_path,
    resolution = res_leukocytes
  )

annotate_lineage_subtypes(
    object = dwStromal,
    lineage_name = "Stromal",
    subtype_list = markers_subtypes[["Stromal"]],
    results_path = results_GEMX_CSA_path,
    resolution = res_stromal
  )


cat(paste("\n ---- FINISHED CELL SUBTYPE ANNOTATION ----
    Generated files:
      · (lineage)_subset.rds
    Generated plots:
      · (lineage)_resolution_grid.png
      · (lineage)_(subtype)_dotplot.png
      · (lineage)_subtype_featureplot.png
      · (lineage)_subtype_UMAP_cluster.png
          "))