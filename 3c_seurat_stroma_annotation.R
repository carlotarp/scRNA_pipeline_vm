##
##  Single Cell Analysis Step 3b: Stroma Annotation
##

# Import libraries
library("Seurat")
library(dplyr)

# Set paths
project_path <- "/home/usuario/PROJECTS/260724_victor_scRNA/"
wd <- paste0(project_path, "codes/scRNA_pipeline/")
setwd(wd)
results_path <- paste0(project_path, "results/")
results_GEMX_CA_path <- paste0(results_path, "GEMX//DecontX/CellAnnotation/")
results_GEMX_SCA_path <- paste0(results_GEMX_CA_path, "Stroma/")

# Import Plot Functions
source(paste0(wd, "CA_plots.R"))

# Load lineage Annotated Data
dwAnnotated <- readRDS(paste0(results_GEMX_CA_path, "fully_annotated_data.rds"))
cat("\n lineage-annotated data loaded \n")

# Set Stroma Markers
markers_Stromal <- list(

  # --- Fibroblast / CAF subtypes ---
  CAF = list(
    cCAF = c("ACTA2", "TAGLN", "MYL9", "MYH11", "POSTN"),          # myofibroblastic, contractile
    iCAF = c("CXCL12", "CXCL14", "IL6", "PDGFRA", "CFD"),           # inflammatory, cytokine-secreting
    apCAF = c("CD74", "HLA-DRA", "HLA-DRB1"),                        # antigen-presenting CAF
    matrixCAF = c("COL1A1", "COL1A2", "COL3A1", "FN1", "LUM"),      # ECM-producing, general fibroblast
    vascularCAF = c("PECAM1", "RGS5", "NOTCH3")                      # perivascular-like CAF (más reciente en literatura)
  ),

  # --- Endothelial subtypes ---
  Endothelial = list(
    Vascular_general = c("PECAM1", "VWF", "CDH5"),
    Arterial = c("GJA5", "SEMA3G", "HEY1"),
    Venous = c("ACKR1", "SELP", "NR2F2"),
    Capillary = c("CA4", "RGCC"),
    Lymphatic = c("PROX1", "PDPN", "LYVE1", "CCL21"),
    Tip_cell = c("ESM1", "ANGPT2", "APLN")                           # angiogénesis activa, punta de brote vascular
  ),

  # --- Mural cells (peri/vascular support) ---
  Mural = list(
    Pericyte = c("RGS5", "PDGFRB", "NOTCH3", "MCAM"),
    Smooth_muscle = c("MYH11", "ACTA2", "DES", "CNN1")
  ),

  # --- Adipocyte ---
  Adipocyte = list( Adipocyte = c("ADIPOQ", "PLIN1", "FABP4", "LEP"))
)

# Generate Subsets
generate_subset <- function(object, results_path) {
  object_subset <- subset(object, subset = lineage == "Stromal")
  object_subset[["RNA_decontX"]] <- split(object_subset[["RNA_decontX"]], f = object_subset$orig.ident)
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
  cat(paste0("  co3 Computed: ", co3_subset, "\n"))

  object_subset <- FindNeighbors(object_subset, reduction = "harmony", dims = 1:co3_subset) # reduction = "harmony"
  object_subset <- RunUMAP(object_subset, reduction = "harmony", dims = 1:co3_subset) # reduction = "harmony"
  object_subset <- FindClusters(object_subset, resolution = 0.2)
  object_subset <- JoinLayers(object_subset)

  # --- Visualize UMAP ---
  plot_dimplot(object_subset, reduction = "umap", group_by = "seurat_clusters",
               results_path = results_path,
               filename = paste0("DimPlot_UMAP_Stroma.png"), label = TRUE)

  return(object_subset)
}

dwStroma <- generate_subset(
    object = dwAnnotated,
    results_path = results_GEMX_SCA_path
  )

  # --- Dotplot by Subtype ---
for (subtype in names(markers_Stromal)){
  plot_marker_dotplot(dwStroma, marker_groups = markers_Stromal[[subtype]],
                       results_path = results_GEMX_SCA_path,
                       filename = paste0("Dotplot_Stroma", subtype, ".png"),
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

# clusters_to_check <- c("5", "7", "8") # 3500
# clusters_to_check <- c("4", "8") # 5500
clusters_to_check <- c("5", "10") # 7500

find_markers_for_clusters(dwStroma, clusters_to_check, results_GEMX_SCA_path)
cat("\n Read Ambiguous Cluster CSV (if needed) to complete the annotation \n")

# Manual Subtype Annotation
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

# -- DimPlot to Compare Contaminated vs Decontaminated
plot_dimplot(dwStroma, reduction = "umap", group_by = "celltype_cont", label = T,
             results_path = results_GEMX_SCA_path, filename = "DimPlot_UMAP_StromaContaminated.png")

# Export Annotated Stroma Data
saveRDS(dwStroma, file.path(results_path, paste0("Stroma.rds")))

# Perform Label Transfer to the Original Seurat Object
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

cat(paste("\n ---- FINISHED Stroma ANNOTATION ----
    Generated files:
      · Stroma.rds
      · notumor_annotated_data.rds
      · notumor_clean_annotated_data.rds
    Generated plots:
      · DimPlot_UMAP_(groupedby).png
      · DotPlot_(groupedby).png
          "))


dwContamina