##
##  Single Cell Analysis Step 3f: Tumor Subtype Annotation
##

# Import libraries
library("Seurat")
library(dplyr)
library(decoupleR)
library(progeny)
library(tidyr)
library(tibble)
library(ggplot2)

# Set paths
project_path <- "/home/usuario/PROJECTS/260724_victor_scRNA/"
wd <- paste0(project_path, "codes/scRNA_pipeline/")
setwd(wd)
results_path <- paste0(project_path, "results/")
results_GEMX_CA_path <- paste0(results_path, "GEMX/DecontX/CellAnnotation/")
results_GEMX_TUMOR_path <- paste0(results_GEMX_CA_path, "Tumor/")

# Load Plot Functions
source(paste0(wd, "CA_plots.R"))
source(paste0(wd, "CL_plots.R"))

# Load Fully Annotated Data
dwAnnotated <- readRDS(paste0(results_GEMX_CA_path, "notumor_annotated_data.rds"))
cat("\n Fully annotated data loaded \n")


# Generate Subsets
generate_subset <- function(object, results_path) {
  object_subset <- subset(object, subset = lineage == "Tumor")
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
  object_subset <- FindClusters(object_subset, resolution = 0.3)
  object_subset <- JoinLayers(object_subset)

  # --- Visualize UMAP ---
  plot_dimplot(object_subset, reduction = "umap", group_by = "seurat_clusters",
               results_path = results_path,
               filename = paste0("DimPlot_UMAP_Tumor.png"), label = TRUE)

  return(object_subset)
}

dwTumoral <- generate_subset(
    object = dwAnnotated,
    results_path = results_GEMX_TUMOR_path
  )

# Set PAM50-like Marker Genes
pam50_genes <- list(
  "Her2+" = c("ERBB2","GRB7","BLVRA","TMEM45B", "FGFR4"),
  Lum = c("ESR1","PGR","BAG1","MAPT","NAT1","ZIP6"),
  TNBC = c("MKI67","CCNE1","ANLN","CDC20","EGFR","MYC")
)

# Compute PAM50 Module Scores
dwTumoral <- AddModuleScore(dwTumoral, features = pam50_genes, name = "PAM50_")

# AddModuleScore names columns PAM50_1, PAM50_2, PAM50_3 in list order - rename
# to the actual subtype names for clarity
pam50_score_cols <- paste0("PAM50_", seq_along(pam50_genes))
pam50_named_cols <- paste0("PAM50_", names(pam50_genes))
colnames(dwTumoral@meta.data)[match(pam50_score_cols, colnames(dwTumoral@meta.data))] <- pam50_named_cols

# Predicted subtype per cell = highest-scoring PAM50 module
pam50_scores <- as.matrix(dwTumoral@meta.data[, pam50_named_cols])
dwTumoral$PAM50_predicted <- names(pam50_genes)[apply(pam50_scores, 1, which.max)]
dwTumoral$PAM50_predicted <- factor(dwTumoral$PAM50_predicted, levels = names(pam50_genes))

cat("\n PAM50 module scores computed \n")
print(table(dwTumoral$PAM50_predicted))

# Write Table of Interest
concordance_table <- table(Predicted = dwTumoral$PAM50_predicted, Clinical = dwTumoral$Subtype)
write.csv(as.data.frame.matrix(concordance_table),
          paste0(results_GEMX_TUMOR_path, "PAM50_vs_Clinical_concordance.csv"))
cat("\n PAM50 predicted vs clinical Subtype concordance table saved \n")


# Compute PROGENy
expr_mat <- as.matrix(GetAssayData(dwTumoral, assay = "RNA_decontX", layer = "data"))
progeny_scores <- progeny::progeny(expr_mat, scale = TRUE, organism = "Human", top = 500, perm = 1)
dwTumoral[["progeny"]] <- CreateAssayObject(data = t(progeny_scores))
dwTumoral <- ScaleData(dwTumoral, assay = "progeny")
cat("\n PROGENy pathway activities computed and added as 'progeny' assay \n")

# Compute Cell Cycle Scores
dwTumoral <- CellCycleScoring(dwTumoral,
                                s.features = cc.genes.updated.2019$s.genes,
                                g2m.features = cc.genes.updated.2019$g2m.genes)
cat("\n Cell cycle scores computed \n")
print(table(dwTumoral$Phase))

# Write Tables of Interest
phase_by_pam50 <- table(PAM50 = dwTumoral$PAM50_predicted, Phase = dwTumoral$Phase)
write.csv(as.data.frame.matrix(phase_by_pam50),
          paste0(results_GEMX_TUMOR_path, "CellCyclePhase_by_PAM50.csv"))

phase_by_cluster <- table(Cluster = dwTumoral$decontX_clusters, Phase = dwTumoral$Phase)
write.csv(as.data.frame.matrix(phase_by_cluster),
          paste0(results_GEMX_TUMOR_path, "CellCyclePhase_by_cluster.csv"))
cat("\n Cell cycle phase tables (by PAM50 and by cluster) saved \n")

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
clusters_to_check <- sort(unique(dwTumoral$seurat_clusters))
find_markers_for_clusters(dwTumoral, clusters_to_check, results_GEMX_TUMOR_path)

# Manual Cluster Annotation
clusters_tumor_annotated <- c( #########  7500  #########
  "0" = "",  "1" = "",
  "4" = "", "7" = "",
  "8" = "", "9" = "",
  "11" = "", "17" = ""
)


dwTumoral$celltype <- dwTumoral$PAM50_predicted
dwTumoral$celltype <- factor(dwTumoral$PAM50_predicted)

dwTumoral$celltype <- unname(clusters_tumor_annotated[as.character(dwTumoral$decontX_clusters)])
dwTumoral$celltype <- factor(dwTumoral$celltype)

# Label Transfer
annotation_vec <- setNames(as.character(dwTumoral[["celltype"]][, 1]),
                            colnames(dwTumoral))

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
             results_path = results_GEMX_CA_path, filename = "DimPlot_UMAP_Annotated.png")


# Export Annoatated Data
saveRDS(dwTumoral, file.path(results_GEMX_TUMOR_path, "tumoral_annotated.rds"))
saveRDS(dwAnnotated, file.path(results_GEMX_CA_path, "fully_annotated_data.rds"))
cat("\n Tumoral subtypes propagated to global 'celltype' column \n")


cat(paste("\n ---- FINISHED TUMOR SUBTYPE ANNOTATION ----
    Generated files:
      · Tumoral_annotated.rds
      · fully_annotated_data.rds
      · PAM50_vs_ClinicalSubtype_concordance.csv
    Generated plots:
      · BarPlot_PAM50_(groupedby).png
      · DimPlot_UMAP_(groupedby).png
      · Dotplot_PROGENy_(groupedby).png
      · HeatMap_PAM50_vs_ClinicalSubtype.png
          "))


