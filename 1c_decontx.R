##
##  Single Cell Analysis Step 1c: DecontX contamination correction
##  Runs AFTER 1b_integrate.R — takes sample_annotated_data.rds as input.
##  Per sample, runs DecontX on already QC'd cells and attaches the corrected
##  counts as a new "RNA_decontX" assay. The original "RNA" assay is untouched.
##

# Import libraries
library(Seurat)
library(dplyr)
library(celda)
library(harmony)
library(ggplot2)

# Set paths
project_path <- "/home/usuario/PROJECTS/260724_victor_scRNA/"
wd <- paste0(project_path, "codes/scRNA_pipeline/")
setwd(wd)
results_path <- paste0(project_path, "results/")
results_GEMX_QC_path <- paste0(results_path, "GEMX/QualityControl/7500/")
results_GEMX_DECONTX_path <- paste0(results_path, "GEMX/DecontX/")

# Import plot functions
source(paste0(wd, "DECONTX_plots.R"))

# Load integrated data
dwQC <- readRDS(file.path(results_GEMX_QC_path, "sample_annotated_data.rds"))
cat("\n QC output loaded: ", ncol(dwQC), " cells, ", length(unique(dwQC$orig.ident)), " samples \n")
dwQC <- JoinLayers(dwQC)

sample_names <- sort(unique(as.character(dwQC$orig.ident)))

# Build empty contamination vector and corrected count list (same dims as full object)
all_genes <- rownames(dwQC)
contamination_full <- setNames(rep(NA_real_, ncol(dwQC)), colnames(dwQC))
corrected_list <- list()

# Run DecontX per sample
for (s in sample_names) {

  cat(paste0("\n\n ==== Sample: ", s, " ==== \n"))

  sample_cells <- colnames(dwQC)[dwQC$orig.ident == s]
  toc <- GetAssayData(dwQC, assay = "RNA", layer = "counts")[, sample_cells]

  # Quick clustering to guide DecontX (low-resolution, per-sample only)
  srat_tmp <- CreateSeuratObject(toc)
  srat_tmp <- NormalizeData(srat_tmp, verbose = FALSE)
  srat_tmp <- FindVariableFeatures(srat_tmp, verbose = FALSE)
  srat_tmp <- ScaleData(srat_tmp, verbose = FALSE)
  srat_tmp <- RunPCA(srat_tmp, npcs = 20, verbose = FALSE)
  srat_tmp <- FindNeighbors(srat_tmp, dims = 1:20, verbose = FALSE)
  srat_tmp <- FindClusters(srat_tmp, resolution = 0.5, verbose = FALSE)
  cluster_labels <- as.integer(factor(srat_tmp$seurat_clusters))

  decontx_res <- decontX(x = toc, z = cluster_labels)
  corrected_list[[s]] <- decontx_res$decontXcounts
  contamination_full[sample_cells] <- decontx_res$contamination

  cat(paste0(" Mean estimated contamination: ", round(mean(decontx_res$contamination), 3), "\n"))
}

corrected_full <- do.call(cbind, corrected_list)
corrected_full <- corrected_full[, colnames(dwQC)]

# Attach corrected counts and contamination to the full object
dwQC[["RNA_decontX"]] <- CreateAssayObject(counts = corrected_full)
dwQC$decontX_contamination <- contamination_full
DefaultAssay(dwQC) <- "RNA_decontX"

dwQC <- NormalizeData(dwQC, assay = "RNA_decontX")
dwQC <- FindVariableFeatures(dwQC, assay = "RNA_decontX")
dwQC <- ScaleData(dwQC, assay = "RNA_decontX")
dwQC <- RunPCA(dwQC, assay = "RNA_decontX", reduction.name = "pca_decontX", reduction.key = "PCdX_")

dwQC <- RunHarmony(
  dwQC,
  group.by.vars = "orig.ident",
  reduction.use = "pca_decontX",
  assay.use = "RNA_decontX",
  reduction.save = "harmony_decontX"
)

saveRDS(dwQC, file.path(results_GEMX_DECONTX_path, "decontx_data.rds"))
cat("\n DecontX correction done — data saved \n")

# Compute genome-wide expression change (before vs after correction)
mean_after <- Matrix::rowMeans(GetAssayData(dwQC, assay = "RNA_decontX", layer = "data"))
DefaultAssay(dwQC) <- "RNA"
mean_before <- Matrix::rowMeans(GetAssayData(dwQC, assay = "RNA", layer = "data"))
DefaultAssay(dwQC) <- "RNA_decontX"

shared_genes <- intersect(names(mean_before), names(mean_after))

gene_change_df <- data.frame(
  gene = shared_genes,
  avg_before = mean_before[shared_genes],
  avg_after = mean_after[shared_genes]
) %>%
  mutate(avg_drop = avg_before - avg_after,
         avg_drop_pct = ifelse(avg_before > 0, (avg_before - avg_after) / avg_before * 100, NA))

# --- Genome-wide expression change plots ---
plot_decontx_scatter(gene_change_df, results_GEMX_DECONTX_path)
plot_decontx_maplot(gene_change_df, results_GEMX_DECONTX_path)
plot_decontx_hist(gene_change_df, results_GEMX_DECONTX_path)
plot_decontx_dumbbell(gene_change_df, results_GEMX_DECONTX_path)

cat(paste("\n ---- FINISHED DECONTX CORRECTION ----
    Run 2a_cluster.R next.
    Generated files:
        · decontx_data.rds  (RNA + RNA_decontX assays, Harmony integration on decontX)
    Generated plots:
        · DecontX_scatter_beforeAfter.png
        · DecontX_MAplot.png
        · DecontX_hist_avgDrop.png
        · DecontX_dumbbell_topGenes.png
        "))
