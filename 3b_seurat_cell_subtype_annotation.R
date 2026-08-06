##
##  Single Cell Analysis Step 3b: Cell Subtype Annotation
##

# Import libraries
library("Seurat")
library(dplyr)

# Set paths
project_path <- "/home/usuario/PROJECTS/260724_victor_scRNA/"
wd <- paste0(project_path, "codes/scRNA_pipeline/")
setwd(wd)
results_path <- paste0(project_path, "results/")
results_GEMX_CA_path <- paste0(results_path, "GEMX/CellAnnotation/")
results_GEMX_LCA_path <- paste0(results_GEMX_CA_path, "Leukocytes/")

# Import Plot Functions
source(paste0(wd, "CA_plots.R"))

# Load lineage Annotated Data
dwAnnotated <- readRDS(paste0(results_path, "GEMX/CellAnnotation/lineage_annotated_data.rds"))
cat("\n lineage-annotated data loaded \n")


# Set Subtype Markers
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

plot_marker_dotplot(dwAnnotated, marker_groups = markers_stromal, results_path = results_GEMX_CA_path,
                     filename = "DotPlot_Stromal1.png", group_by = "seurat_clusters", title = "Stromal Markers by Cluster")

plot_marker_dotplot(dwAnnotated, marker_groups = markers_stromaangio, results_path = results_GEMX_CA_path,
                     filename = "DotPlot_Stromal2.png", group_by = "seurat_clusters", title = "Stromal Markers by Cluster")

# Manual Subtype Annotation
clusters_stromal_annotated <- c(
  "0" = "Tumoral", "1" = "CAF/Fibroblast/Pericyte",     "2" = "Tumoral",
  "3" = "Leukocytes",     "4" = "Leukocytes",      "5" = "Tumoral",
  "6" = "Tumoral", "7" = "Leukocytes",      "8" = "Leukocytes",
  "9" = "cCAF/Endothelial", "10" = "Tumoral",
  "11" = "Leukocytes",     "12" = "Leukocytes"
)

dwAnnotated$celltype <- unname(clusters_stromal_annotated[as.character(dwAnnotated$seurat_clusters)])
dwAnnotated$celltype <- factor(dwAnnotated$celltype)

# ---  Stromal Annotated Dimplot ---
plot_dimplot(dwAnnotated, reduction = "umap", group_by = "celltype",
             results_path = results_GEMX_CA_path, filename = "DimPlot_UMAP_Stromal.png")

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
  )
)



# Generate Subsets
generate_subset <- function(object, results_path) {
  object_subset <- subset(object, subset = lineage == "Leukocytes")
  object_subset <- NormalizeData(object_subset)
  object_subset <- FindVariableFeatures(object_subset)
  object_subset <- ScaleData(object_subset)
  object_subset <- RunPCA(object_subset)

#  object_subset <- IntegrateLayers(
#    object = object_subset, method = HarmonyIntegration,
#    orig.reduction = "pca", new.reduction = "harmony", verbose = FALSE
#  )

  pct <- object_subset[["pca"]]@stdev / sum(object_subset[["pca"]]@stdev) * 100
  cumu <- cumsum(pct)
  co1 <- which(cumu > 90 & pct < 5)[1]
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = TRUE)[1] + 1
  co3_subset <- min(co1, co2)
  cat(paste0("  co3 Computed: ", co3_subset, "\n"))

  object_subset <- FindNeighbors(object_subset, reduction = "pca", dims = 1:co3_subset) # reduction = "harmony"
  object_subset <- RunUMAP(object_subset, reduction = "pca", dims = 1:co3_subset) # reduction = "harmony"
  object_subset <- FindClusters(object_subset, resolution = 0.4)
  object_subset <- JoinLayers(object_subset)

  # --- Visualize UMAP ---
  plot_dimplot(object_subset, reduction = "umap", group_by = "seurat_clusters",
               results_path = results_path,
               filename = paste0("DimPlot_UMAP_Leukocytes.png"), label = TRUE)

  saveRDS(object_subset, file.path(results_path, paste0("leukocytes.rds")))

  return(object_subset)
}

# Execute Functions by Lineage
dwLeukocytes <- generate_subset(
    object = dwAnnotated,
    results_path = results_GEMX_LCA_path
  )

  # --- Dotplot by Subtype ---
for (subtype in names(markers_leukocytes)){
  plot_marker_dotplot(dwLeukocytes, marker_groups = markers_leukocytes[[subtype]],
                       results_path = results_GEMX_LCA_path,
                       filename = paste0("Dotplot_Leukocytes", subtype, ".png"),
                       group_by = "seurat_clusters")
}

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

clusters_to_check <- c("5", "7", "8", "13")
find_markers_for_clusters(dwLeukocytes, clusters_to_check, results_GEMX_LCA_path)
cat("\n Read Ambiguous Cluster CSV (if needed) to complete the annotation \n")

# Manual Subtype Annotation
clusters_leuko_annotated <- c(
  "0" = "TCell_naive",             "1" = "pDC // TAM",         "2" = "TCell_ct",
  "3" = "pDC // TAM",         "4" = "BCell_act",        "5" = "",
  "6" = "PlasmaBlast",  "7" = "",             "8" = "",
  "9" = "Mast",         "10" = "",            "11" = "pDC // TCell",
  "12" = "pDC // TAM",            "13" = "",
  "14" = "actDC",       "15" = "PlasmaBlast"
)


dwLeukocytes$celltype <- unname(clusters_leuko_annotated[as.character(dwLeukocytes$seurat_clusters)])
dwLeukocytes$celltype <- factor(dwLeukocytes$celltype)


annotation_vec <- setNames(as.character(dwLeukocytes[["celltype"]][, 1]),
                            colnames(dwLeukocytes))

if ("celltype" %in% colnames(dwAnnotated@meta.data)) {
  existing <- as.character(dwAnnotated[["celltype"]][, 1])
} else {
  existing <- rep(NA_character_, ncol(dwAnnotated))
}
names(existing) <- colnames(dwAnnotated)

# Only overwrite the cells present in the subset - everything else (e.g.
# already-annotated Stromal cells) stays exactly as it was
existing[names(annotation_vec)] <- annotation_vec

dwAnnotated[["celltype"]] <- factor(existing)

# ---  Visualize Annotated Dimplot ---
plot_dimplot(dwAnnotated, reduction = "umap", group_by = "celltype",
             results_path = results_GEMX_CA_path, filename = "DimPlot_UMAP_NoTumor.png")

# Export Annotated Data
saveRDS(dwAnnotated, file.path(results_GEMX_CA_path, "notumor_annotated_data.rds"))

cat(paste("\n ---- FINISHED CELL SUBTYPE ANNOTATION ----
    Generated files:
      · (lineage)_subset.rds
    Generated plots:
      · (lineage)_(subtype)_dotplot.png
      · (lineage)_subtype_featureplot.png
      · (lineage)_subtype_UMAP_cluster.png
          "))