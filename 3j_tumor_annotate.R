##
##  Single Cell Analysis Step 3j: Tumor cluster annotation
##  Runs AFTER 3i_tumor_markers.R — takes tumoral_scored.rds and
##  copykat_annotated_data.rds as input.
##  Applies the manual tumor dictionary, drops the noise clusters, and
##  propagates the labels to the full object.
##

library(Seurat)
library(dplyr)
library(progeny)
library(tidyr)
library(tibble)
library(ggplot2)

# --- Paths ---
project_path <- "/home/user/PROJECTS/scRNA_vmendez/"
wd <- paste0(project_path, "codes/")
setwd(wd)
results_path <- paste0(project_path, "results/")
results_GEMX_CA_path <- paste0(results_path, "CellAnnotation/")
results_GEMX_TUMOR_path <- paste0(results_GEMX_CA_path, "Tumor/")

# --- Plot functions & helpers ---
source(paste0(wd, "TCA_plots.R"))
source(paste0(wd, "CA_plots.R"))
source(paste0(wd, "CL_plots.R"))
source(paste0(wd, "utils.R"))

# =============================================================================
# Manual choices, filled in from the 3i_tumor_markers.R plots and cluster CSVs
# =============================================================================

# Tumor cluster annotation
clusters_tumor_annotated <- c(
  "0" = "Her2+_Quiescent", "1" = "TNBC_Hypoxic", "2" = "ProliferativeTumor_S",
  "3" = "NK/TCell_contamination", "4" = "Lum_Quiescent", "5" = "ProliferativeTumor_G2M",
  "6" = "CAF/TAM_contamination", "7" = "LumB_SC11", "8" = "TNBC_Endothelial",
  "9" = "Inflammatory", "10" = "Interferon", "11" = "Her2+_SC13",
  "12" = "ProliferativeTumor_SC11"
)

# Contaminating clusters dropped from the final annotated objects
noise_decontX_clusters <- c("18", "19", "20", "21", "22")
noise_tumor_clusters <- c("13")

# --- Load scored tumor data and full annotated object ---
dwTumoral <- readRDS(paste0(results_GEMX_TUMOR_path, "tumoral_scored.rds"))
dwAnnotated <- readRDS(paste0(results_GEMX_CA_path, "CopyKAT/copykat_annotated_data.rds"))
cat("\n Tumor scored data and full annotated data loaded \n")


# =============================================================================
# 1. Drop the contaminating (noise) clusters
# =============================================================================

dwAnnotated$celltype <- as.character(dwAnnotated$celltype)
dwAnnotated$lineage <- as.character(dwAnnotated$lineage)
dwAnnotated$celltype[dwAnnotated$decontX_clusters %in% noise_decontX_clusters] <- "Noise"
dwAnnotated$lineage[dwAnnotated$decontX_clusters %in% noise_decontX_clusters] <- "Noise"
dwTumoral <- subset(dwTumoral, subset = !decontX_clusters %in% noise_decontX_clusters)

noise_cells_tumor <- colnames(dwTumoral)[dwTumoral$seurat_clusters %in% noise_tumor_clusters]
dwAnnotated$celltype[colnames(dwAnnotated) %in% noise_cells_tumor] <- "Noise"
dwAnnotated$lineage[colnames(dwAnnotated) %in% noise_cells_tumor] <- "Noise"
dwAnnotated$celltype <- factor(dwAnnotated$celltype)
dwAnnotated$lineage <- factor(dwAnnotated$lineage)
dwTumoral <- subset(dwTumoral, subset = !seurat_clusters %in% noise_tumor_clusters)

plot_dimplot(dwAnnotated, reduction = "umap_decontX", group_by = "celltype",
             results_path = results_GEMX_CA_path, filename = "DimPlot_UMAP_Cleaned.png")
plot_dimplot(dwTumoral, reduction = "umap", group_by = "seurat_clusters", label = TRUE,
             results_path = results_GEMX_CA_path, filename = "DimPlot_UMAP_TumorCleaned.png")
cat("\n Noise clusters dropped \n")

# =============================================================================
# 2. Manual cluster annotation & label propagation to the full object
# =============================================================================

dwTumoral$celltype <- factor(unname(clusters_tumor_annotated[as.character(dwTumoral$seurat_clusters)]),
                             levels = unique(unname(clusters_tumor_annotated)))

# --- Label transfer to full annotated object ---
annotation_vec <- setNames(as.character(dwTumoral[["celltype"]][, 1]), colnames(dwTumoral))

if ("celltype" %in% colnames(dwAnnotated@meta.data)) {
  existing <- as.character(dwAnnotated[["celltype"]][, 1])
} else {
  existing <- rep(NA_character_, ncol(dwAnnotated))
}
names(existing) <- colnames(dwAnnotated)
existing[names(annotation_vec)] <- annotation_vec
dwAnnotated[["celltype"]] <- factor(existing)


# --- Annotated dimplot (full object) ---
dwAnnotated <- subset(dwAnnotated, subset = celltype != "Noise")  # Remove contaminating cluster for final annotated object
plot_dimplot(dwAnnotated, reduction = "umap_decontX", group_by = "celltype",
             results_path = results_GEMX_CA_path, filename = "DimPlot_UMAP_Annotated.png")
plot_dimplot(dwTumoral, reduction = "umap_cna", group_by = "celltype",
             results_path = results_GEMX_TUMOR_path, filename = "DimPlot_UMAP_CNA_celltype_Tumoral.png")

# --- Export final annotated data ---
saveRDS(dwTumoral, file.path(results_GEMX_TUMOR_path, "tumoral_annotated.rds"))
saveRDS(dwAnnotated, file.path(results_GEMX_CA_path, "fully_annotated_data.rds"))
cat("\n Tumor annotations propagated to global 'celltype' column \n")

# --- PAM50 call proportion within each sample's total tumor cells, faceted by celltype ---
pam50_prop_df <- compute_celltype_proportions(dwTumoral, sample_id_col = "celltype",
                                                  celltype_col = "PAM50_predicted",
                                                  clinical_cols = "celltype")
plot_celltype_proportions_bar(pam50_prop_df, x_var = "PAM50_predicted",
                              results_path = results_GEMX_TUMOR_path,
                              filename = "Barplot_PAM50Proportion_byCelltype.png",
                              celltype_col = "celltype", colors = pam50_subtype_colors)

# --- scSubtype call proportion within each sample's total tumor cells, faceted by celltype ---
scsubtype_prop_df <- compute_celltype_proportions(dwTumoral, sample_id_col = "celltype",
                                                  celltype_col = "scSubtype_call",
                                                  clinical_cols = "celltype")
plot_celltype_proportions_bar(scsubtype_prop_df, x_var = "scSubtype_call",
                              results_path = results_GEMX_TUMOR_path,
                              filename = "Barplot_scSubtypeProportion_byCelltype.png",
                              celltype_col = "celltype", colors = scsubtype_colors)

cat(paste("\n ---- FINISHED TUMOR ANNOTATION ----
    Run 4a_clinical_analysis.R next.
    Generated files:
        · tumoral_annotated.rds  (final celltype labels)
        · fully_annotated_data.rds  (full object, all lineages annotated)
    Generated plots:
        · DimPlot_UMAP_Annotated.png, DimPlot_UMAP_(Tumor)Cleaned.png
        · Barplot_(PAM50|scSubtype)Proportion_byCelltype.png
        "))
