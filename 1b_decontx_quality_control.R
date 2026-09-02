##
##  Single Cell Analysis Step 1c: DecontX correction
##  Runs AFTER seurat_1_QC.R - takes its output (sample_annotated_data.rds:
##  merged, QC'd, doublet-removed, Harmony-integrated, subtype-annotated)
##  as input. Per sample, runs DecontX on that sample's already-QC'd cells,
##  then attaches the corrected counts back onto the SAME object as a new
##  assay ("RNA_decontX"). The original "RNA" assay is left untouched.
##

library(Seurat)
library(dplyr)
library(celda)
library(harmony)
library(ggplot2)

project_path <- "/home/usuario/PROJECTS/260724_victor_scRNA/"
results_path <- paste0(project_path, "results/")
results_GEMX_QC_path <- paste0(results_path, "GEMX/QualityControl/7500/")
results_GEMX_DECONTX_path <- paste0(results_path, "GEMX/DecontX/")

# --- Load the output of seurat_1_QC.R ---
dwQC <- readRDS(file.path(results_GEMX_QC_path, "sample_annotated_data.rds"))
cat("\n QC output loaded: ", ncol(dwQC), " cells, ", length(unique(dwQC$orig.ident)), " samples \n")
dwQC <- JoinLayers(dwQC)

sample_names <- sort(unique(as.character(dwQC$orig.ident)))

# Empty matrix to fill in, same genes/cells as the full object
all_genes <- rownames(dwQC)
contamination_full <- setNames(rep(NA_real_, ncol(dwQC)), colnames(dwQC))
corrected_list <- list()

# =================================================================
# Run DecontX per sample, on cells already QC'd/doublet-filtered by
# seurat_1_QC.R (no re-loading raw CellRanger matrices needed)
# =================================================================
for (s in sample_names) {

  cat(paste0("\n\n ==== Sample: ", s, " ==== \n"))

  sample_cells <- colnames(dwQC)[dwQC$orig.ident == s]
  toc <- GetAssayData(dwQC, assay = "RNA", layer = "counts")[, sample_cells]

  # Quick, low-resolution clustering just to guide DecontX
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

# =================================================================
# Attach corrected counts + contamination onto the full object
# =================================================================
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

# --- Scatter before vs after, con línea 1:1 - visión global de cuánto se aleja cada gen de "sin cambio" ---
p_scatter <- ggplot(gene_change_df, aes(x = avg_before, y = avg_after)) +
  geom_point(size = 0.4, alpha = 0.3, color = "steelblue") +
  geom_abline(slope = 1, intercept = 0, color = "firebrick", linetype = "dashed") +
  theme_bw() + theme(panel.grid = element_blank()) +
  labs(title = "Genome-wide expression: before vs after DecontX",
       subtitle = "Points below the red line = genes reduced by correction",
       x = "Average expression (before)", y = "Average expression (after)")
ggsave(paste0(results_GEMX_DECONTX_path, "GenomeWide_scatter_beforeAfter.png"), p_scatter,
       width = 7, height = 7, dpi = 300, bg = "white")


# --- MA-plot - magnitud de expresión vs magnitud del cambio, con top genes etiquetados ---
top_labeled <- gene_change_df %>% slice_max(order_by = avg_drop, n = 15)

p_ma <- ggplot(gene_change_df, aes(x = avg_before, y = avg_drop)) +
  geom_point(size = 0.4, alpha = 0.3, color = "grey50") +
  geom_point(data = top_labeled, color = "firebrick", size = 1.2) +
  geom_text(data = top_labeled, aes(label = gene), size = 3, vjust = -0.6, check_overlap = TRUE) +
  theme_bw() + theme(panel.grid = element_blank()) +
  labs(title = "MA-plot: expression level vs DecontX correction magnitude",
       subtitle = "Top 15 most affected genes highlighted",
       x = "Average expression (before)", y = "avg_before - avg_after") +
  ylim(c(0,2)) + xlim(c(0,3))
ggsave(paste0(results_GEMX_DECONTX_path, "GenomeWide_MAplot.png"), p_ma, width = 9, height = 7, dpi = 300, bg = "white")


# --- Distribución del % de reducción, a nivel genoma completo ---
p_hist <- ggplot(gene_change_df %>% filter(!is.na(avg_drop_pct)), aes(x = avg_drop_pct)) +
  geom_histogram(bins = 60, fill = "steelblue", alpha = 0.7) +
  theme_bw() + theme(panel.grid = element_blank()) +
  labs(title = "Distribution of % expression reduction across all genes",
       x = "% reduction (avg_before -> avg_after)", y = "N genes")
ggsave(paste0(results_GEMX_DECONTX_path, "GenomeWide_pctDrop_histogram.png"), p_hist,
       width = 8, height = 5, dpi = 300, bg = "white")


# --- Top 30 genes más afectados, como dumbbell (before -> after) ---
top30_dumbbell <- gene_change_df %>%
  slice_max(order_by = avg_drop, n = 30) %>%
  mutate(gene = reorder(gene, avg_drop))

p_top30_dumbbell <- ggplot(top30_dumbbell) +
  geom_segment(aes(x = avg_before, xend = avg_after, y = gene, yend = gene), color = "grey60") +
  geom_point(aes(x = avg_before, y = gene), color = "firebrick", size = 2.5) +
  geom_point(aes(x = avg_after, y = gene), color = "steelblue", size = 2.5) +
  theme_bw() + theme(panel.grid.minor = element_blank()) +
  labs(title = "Top 30 genes most affected by DecontX: before (red) vs after (blue)",
       x = "Average expression", y = NULL)

ggsave(paste0(results_GEMX_DECONTX_path, "GenomeWide_top30_dumbbell.png"), p_top30_dumbbell,
       width = 8, height = 9, dpi = 300, bg = "white")

cat(paste("\n ---- FINISHED DECONTX CORRECTION ----
    Generated files:
      · sample_annotated_data_decontX.rds (RNA + RNA_decontX assays, both on the full merged/integrated object)
      · (sample)_collagen_before_after.csv
          "))