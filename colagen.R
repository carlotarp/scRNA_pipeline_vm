##
##  Collagen distribution analysis - BEFORE DecontX correction
##  By sample, by cluster (global Level 0), and by celltype
##


# Set paths
project_path <- "/home/usuario/PROJECTS/260724_victor_scRNA/"
wd <- paste0(project_path, "codes/scRNA_pipeline/")
setwd(wd)
results_path <- paste0(project_path, "results/")
results_GEMX_CC_path <- paste0(results_path, "GEMX/ContaminationCorrection/")
results_GEMX_CA_path <- paste0(results_path, "GEMX/CellAnnotation/7500/")

# Load Fully Annotated Data
dwAnnotated <- readRDS(file.path(results_GEMX_CA_path, "fully_annotated_data.rds"))
cat("\n Fully annotated data loaded \n")


library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)

collagen_genes <- c("COL1A1", "COL1A2", "COL3A1")
collagen_genes <- collagen_genes[collagen_genes %in% rownames(dwAnnotated)]

dwAnnotated_joined <- JoinLayers(dwAnnotated)

compute_pct_avg_by_group <- function(object, genes, group_var) {
  expr_df <- FetchData(object, vars = c(genes, group_var))
  colnames(expr_df)[colnames(expr_df) == group_var] <- "group"

  expr_df %>%
    pivot_longer(cols = all_of(genes), names_to = "gene", values_to = "expr") %>%
    group_by(group, gene) %>%
    summarise(pct_expressing = round(mean(expr > 0) * 100, 1),
              avg_expr = round(mean(expr), 3),
              n_cells = n(), .groups = "drop")
}

plot_expr_boxplot <- function(object, genes, group_var, title, filename, results_path) {
  expr_df <- FetchData(object, vars = c(genes, group_var))
  colnames(expr_df)[colnames(expr_df) == group_var] <- "group"

  expr_long <- expr_df %>%
    pivot_longer(cols = all_of(genes), names_to = "gene", values_to = "expr")

  p <- ggplot(expr_long, aes(x = group, y = expr, fill = gene)) +
    geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.3) +
    facet_wrap(~gene, scales = "free_y") +
    theme_bw() + theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1),
                         legend.position = "none") +
    labs(title = title, x = NULL, y = "Expression (normalized)")

  n_groups <- length(unique(expr_long$group))
  ggsave(paste0(results_path, filename), p, width = max(8, n_groups * 0.4), height = 6, dpi = 300, bg = "white")
}

# --- By sample ---
by_sample <- compute_pct_avg_by_group(dwAnnotated_joined, collagen_genes, "orig.ident")
write.csv(by_sample, paste0(results_GEMX_CC_path, "Collagen_bySample.csv"), row.names = FALSE)

p_sample <- ggplot(by_sample, aes(x = group, y = pct_expressing, fill = gene)) +
  geom_col(position = "dodge") +
  theme_bw() + theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "% cells expressing collagen genes, by sample", x = "Sample", y = "% expressing")
ggsave(paste0(results_GEMX_CC_path, "Collagen_bySample.png"), p_sample, width = 10, height = 5, dpi = 300, bg = "white")

plot_expr_boxplot(dwAnnotated_joined, collagen_genes, "orig.ident",
                   "Collagen expression distribution, by sample", "Collagen_bySample_boxplot.png",
                   results_GEMX_CC_path)

# --- By cluster (global Level 0) ---
by_cluster <- compute_pct_avg_by_group(dwAnnotated_joined, collagen_genes, "seurat_clusters")
write.csv(by_cluster, paste0(results_GEMX_CC_path, "Collagen_byCluster.csv"), row.names = FALSE)

p_cluster <- ggplot(by_cluster, aes(x = group, y = pct_expressing, fill = gene)) +
  geom_col(position = "dodge") +
  theme_bw() + theme(panel.grid = element_blank()) +
  labs(title = "% cells expressing collagen genes, by cluster", x = "Cluster", y = "% expressing")
ggsave(paste0(results_GEMX_CC_path, "Collagen_byCluster.png"), p_cluster, width = 10, height = 5, dpi = 300, bg = "white")

plot_expr_boxplot(dwAnnotated_joined, collagen_genes, "seurat_clusters",
                   "Collagen expression distribution, by sample", "Collagen_byCluster_boxplot.png",
                   results_GEMX_CC_path)

# --- By Clinic Subtype ---
by_subtype <- compute_pct_avg_by_group(dwAnnotated_joined, collagen_genes, "Subtype")
write.csv(by_subtype, paste0(results_GEMX_CC_path, "Collagen_bySubtype.csv"), row.names = FALSE)

p_subtype <- ggplot(by_subtype, aes(x = group, y = pct_expressing, fill = gene)) +
  geom_col(position = "dodge") +
  theme_bw() + theme(panel.grid = element_blank()) +
  labs(title = "% cells expressing collagen genes, by subtype", x = "Subtype", y = "% expressing")
ggsave(paste0(results_GEMX_CC_path, "Collagen_bySubtype.png"), p_subtype, width = 10, height = 5, dpi = 300, bg = "white")

plot_expr_boxplot(dwAnnotated_joined, collagen_genes, "Subtype",
                   "Collagen expression distribution, by sample", "Collagen_bySubtype_boxplot.png",
                   results_GEMX_CC_path)

# --- By celltype - the most important one: which non-stromal celltypes are
#     showing collagen they shouldn't (the SoupX-paper-style artifact check) ---
by_celltype <- compute_pct_avg_by_group(dwAnnotated_joined, collagen_genes, "celltype")
write.csv(by_celltype, paste0(results_GEMX_CC_path, "Collagen_byCelltype.csv"), row.names = FALSE)

p_celltype <- ggplot(by_celltype, aes(x = reorder(group, pct_expressing), y = pct_expressing, fill = gene)) +
  geom_col(position = "dodge") +
  coord_flip() +
  theme_bw() + theme(panel.grid = element_blank()) +
  labs(title = "% cells expressing collagen genes, by celltype", x = NULL, y = "% expressing")
ggsave(paste0(results_GEMX_CC_path, "Collagen_byCelltype.png"), p_celltype,
       width = 8, height = max(5, length(unique(by_celltype$group)) * 0.3), dpi = 300, bg = "white")

plot_expr_boxplot(dwAnnotated_joined, collagen_genes, "celltype",
                   "Collagen expression distribution, by sample", "Collagen_byCelltype_boxplot.png",
                   results_GEMX_CC_path)

cat("\n ---- FINISHED collagen distribution analysis (pre-DecontX) ---- \n")
cat(" Files: Collagen_bySample.csv/png, Collagen_byCluster.csv/png, Collagen_byCelltype.csv/png \n")

# Print the celltype table sorted - immune/non-stromal celltypes with high
# collagen pct_expressing here are exactly the artifact pattern from the
# SoupX paper (T/MNP clusters showing COL1A1/2/3A1 that shouldn't be there)
cat("\n Celltype table, sorted by pct_expressing (highest first): \n")
print(by_celltype %>% arrange(desc(pct_expressing)))


# =================================================================
# B. Merge DecontX-corrected counts into the fully annotated object
#    Adds a NEW assay ("RNA_decontX") - the original "RNA" assay is kept
#    untouched, so before/after can always be compared directly.
# =================================================================
sample_list <- c("SC7b", "SC8", "SC9", "SC10", "SC11", "SC12", "SC13", "SC14",
                  "SC15", "SC16", "SC17", "SC18", "SC19", "SC20", "SC21", "SC5_SC22")

all_genes <- rownames(dwAnnotated_joined)
corrected_full <- Matrix::Matrix(0, nrow = length(all_genes), ncol = ncol(dwAnnotated_joined),
                                   sparse = TRUE, dimnames = list(all_genes, colnames(dwAnnotated_joined)))

for (s in sample_list) {

  corrected_s <- readRDS(paste0(results_GEMX_CC_path, s, "_corrected_counts.rds"))

  # Cells belonging to this sample in the big annotated object
  sample_cells <- colnames(dwAnnotated_joined)[dwAnnotated_joined$orig.ident == s]

  # Seurat appends a numeric suffix ("_1", "_4", etc.) to barcodes when
  # merging samples - strip it to recover the raw barcode and match against
  # the DecontX matrix's original column names
  raw_barcodes <- sub("_[0-9]+$", "", sample_cells)

  matched_idx <- match(raw_barcodes, colnames(corrected_s))
  valid <- !is.na(matched_idx)

  if (sum(!valid) > 0) {
    warning(s, ": ", sum(!valid), " cells could not be matched to the DecontX output - check barcode format")
  }

  shared_genes <- intersect(rownames(corrected_s), all_genes)
  corrected_full[shared_genes, sample_cells[valid]] <- corrected_s[shared_genes, matched_idx[valid]]

  cat(s, ": ", sum(valid), " cells merged\n")
}

dwAnnotated_joined[["RNA_decontX"]] <- CreateAssayObject(counts = corrected_full)
dwAnnotated_joined <- NormalizeData(dwAnnotated_joined, assay = "RNA_decontX", verbose = FALSE)
cat("\n [B] DecontX-corrected counts merged as assay 'RNA_decontX' \n")


# =================================================================
# Función genérica de comparación antes/después DecontX, para cualquier
# agrupación (celltype, cluster, sample...)
# =================================================================
compare_before_after_decontx <- function(object, genes, group_var, group_label,
                                            results_path_post) {

  # --- Before (assay RNA) ---
  DefaultAssay(object) <- "RNA"
  by_group_before <- compute_pct_avg_by_group(object, genes, group_var)

  # --- After (assay RNA_decontX) ---
  DefaultAssay(object) <- "RNA_decontX"
  by_group_after <- compute_pct_avg_by_group(object, genes, group_var)
  DefaultAssay(object) <- "RNA"

  write.csv(by_group_after, paste0(results_path_post, "Collagen_by", group_label, "_postDecontX.csv"), row.names = FALSE)

  # --- Comparación - basada en avg_expr, no en pct_expressing (DecontX
  #     redistribuye magnitud de counts, no cambia si el gen se detecta o no) ---
  comparison <- by_group_before %>%
    dplyr::rename(pct_before = pct_expressing, avg_before = avg_expr, n_cells_before = n_cells) %>%
    dplyr::left_join(
      by_group_after %>% dplyr::rename(pct_after = pct_expressing, avg_after = avg_expr, n_cells_after = n_cells),
      by = c("group", "gene")
    ) %>%
    dplyr::mutate(avg_drop = avg_before - avg_after,
                  avg_drop_pct = round((avg_before - avg_after) / avg_before * 100, 1))

  write.csv(comparison, paste0(results_path_post, "Collagen_by", group_label, "_before_vs_after.csv"), row.names = FALSE)

  # --- Dumbbell plot (avg_expr, no pct_expressing) ---
  comparison_long <- comparison %>% dplyr::mutate(group = reorder(group, avg_drop))

  p_dumbbell <- ggplot(comparison_long) +
    geom_segment(aes(x = avg_before, xend = avg_after, y = group, yend = group), color = "grey60") +
    geom_point(aes(x = avg_before, y = group), color = "firebrick", size = 2.5) +
    geom_point(aes(x = avg_after, y = group), color = "steelblue", size = 2.5) +
    facet_wrap(~gene, scales = "free_x") +
    theme_bw() + theme(panel.grid.minor = element_blank()) +
    labs(title = paste0("Collagen expression: before (red) vs after (blue), by ", group_label),
         x = "Average expression (normalized)", y = NULL)
  ggsave(paste0(results_path_post, "Collagen_beforeAfter_dumbbell_", group_label, ".png"), p_dumbbell,
         width = 11, height = max(5, length(unique(comparison_long$group)) * 0.35), dpi = 300, bg = "white")

  # --- Drop plot (% de reducción relativa en avg_expr) ---
  p_drop <- ggplot(comparison %>% dplyr::mutate(group = reorder(group, avg_drop_pct)),
                     aes(x = group, y = avg_drop_pct, fill = gene)) +
    geom_col(position = "dodge") +
    coord_flip() +
    theme_bw() + theme(panel.grid = element_blank()) +
    labs(title = paste0("% reduction in collagen expression after DecontX, by ", group_label),
         x = NULL, y = "% reduction in avg_expr")
  ggsave(paste0(results_path_post, "Collagen_avgDrop_sorted_", group_label, ".png"), p_drop,
         width = 8, height = max(5, length(unique(comparison$group)) * 0.3), dpi = 300, bg = "white")
  # --- Drop plot en absoluto (avg_before - avg_after) ---
  p_drop_abs <- ggplot(comparison %>% dplyr::mutate(group = reorder(group, avg_drop)),
                         aes(x = group, y = avg_drop, fill = gene)) +
    geom_col(position = "dodge") +
    coord_flip() +
    theme_bw() + theme(panel.grid = element_blank()) +
    labs(title = paste0("Reduction in collagen expression after DecontX (absolute), by ", group_label),
         x = NULL, y = "avg_before - avg_after")
  ggsave(paste0(results_path_post, "Collagen_avgDrop_absolute_", group_label, ".png"), p_drop_abs,
         width = 8, height = max(5, length(unique(comparison$group)) * 0.3), dpi = 300, bg = "white")

  return(comparison)
}

# --- Aplicar a las 3 dimensiones ---
comparison_celltype <- compare_before_after_decontx(dwAnnotated_joined, collagen_genes, "celltype", "Celltype",
                                                       results_GEMX_CC_path)

comparison_cluster <- compare_before_after_decontx(dwAnnotated_joined, collagen_genes, "seurat_clusters", "Cluster",
                                                       results_GEMX_CC_path)

comparison_sample <- compare_before_after_decontx(dwAnnotated_joined, collagen_genes, "orig.ident", "Sample",
                                                      results_GEMX_CC_path)


DefaultAssay(dwAnnotated_joined) <- "RNA_decontX"

# --- Standard preprocessing on the corrected assay ---
dwAnnotated_joined <- FindVariableFeatures(dwAnnotated_joined, assay = "RNA_decontX", verbose = FALSE)
dwAnnotated_joined <- ScaleData(dwAnnotated_joined, assay = "RNA_decontX", verbose = FALSE)
dwAnnotated_joined <- RunPCA(dwAnnotated_joined, assay = "RNA_decontX",
                               reduction.name = "pca_decontX", reduction.key = "PCdX_", verbose = FALSE)

# --- co3: mismo criterio de selección de dimensiones que ya usas en el resto del pipeline ---
pct <- dwAnnotated_joined[["pca_decontX"]]@stdev / sum(dwAnnotated_joined[["pca_decontX"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = TRUE)[1] + 1
co3 <- min(co1, co2, na.rm = TRUE)

cat("co3 (RNA_decontX):", co3, "\n")

# --- Harmony + clustering + UMAP usando co3 en vez de un valor fijo ---
library(harmony)
dwAnnotated_joined <- RunHarmony(dwAnnotated_joined,
                                    group.by.vars = "orig.ident",
                                    reduction.use = "pca_decontX",
                                    assay.use = "RNA_decontX",
                                    reduction.save = "harmony_decontX")

dwAnnotated_joined <- FindNeighbors(dwAnnotated_joined, reduction = "harmony_decontX", dims = 1:co3)
dwAnnotated_joined <- FindClusters(dwAnnotated_joined, resolution = 0.4, cluster.name = "decontX_clusters")
dwAnnotated_joined <- RunUMAP(dwAnnotated_joined, reduction = "harmony_decontX", dims = 1:co3,
                                 reduction.name = "umap_decontX")

# ¿Cambia el número de clusters?
table(dwAnnotated_joined$seurat_clusters)
table(dwAnnotated_joined$decontX_clusters)

source(paste0(wd, "CA_plots.R"))
plot_dimplot(dwAnnotated_joined, reduction = "umap_decontX", group_by = "celltype",
             results_path = results_GEMX_CC_path, filename = "UMAP_decontX_byCelltype.png")
plot_dimplot(dwAnnotated_joined, reduction = "umap_decontX", group_by = "seurat_clusters",
             results_path = results_GEMX_CC_path, filename = "UMAP_decontX_bySeuratClusters.png")
plot_dimplot(dwAnnotated_joined, reduction = "umap_decontX", group_by = "decontX_clusters", label = T,
             results_path = results_GEMX_CC_path, filename = "UMAP_decontX_byDecontCLusters.png")
plot_dimplot(dwAnnotated_joined, reduction = "pca_decontX", group_by = "orig.ident", label = T,
             results_path = results_GEMX_CC_path, filename = "PCA_decontX_bySample.png")
DefaultAssay(dwAnnotated_joined) <- "RNA"


# =================================================================
# Merge per-cell contamination fraction into the annotated object,
# matching barcodes (same suffix-stripping pattern as the DecontX assay merge)
# =================================================================
dwAnnotated_joined$contamination_fraction <- NA_real_

for (s in sample_list) {

  contamination_df <- read.csv(paste0(results_DECONTX_path, s, "_contamination_fraction.csv"))

  sample_cells <- colnames(dwAnnotated_joined)[dwAnnotated_joined$orig.ident == s]
  raw_barcodes <- sub("_[0-9]+$", "", sample_cells)

  matched_idx <- match(raw_barcodes, contamination_df$cell)
  valid <- !is.na(matched_idx)

  dwAnnotated_joined$contamination_fraction[sample_cells[valid]] <- contamination_df$contamination[matched_idx[valid]]

  cat(s, ": ", sum(valid), " cells matched\n")
}

cat("\n Cells with contamination value assigned: ", sum(!is.na(dwAnnotated_joined$contamination_fraction)), " of ", ncol(dwAnnotated_joined), "\n")


# =================================================================
# Cross-reference tables: mean/median contamination by sample, cluster
# (pre-DecontX, i.e. seurat_clusters), and celltype
# =================================================================
summarise_contamination <- function(object, group_var, group_label) {
  df <- object@meta.data[, c(group_var, "contamination_fraction")]
  colnames(df)[colnames(df) == group_var] <- "group"

  summary_df <- df %>%
    filter(!is.na(contamination_fraction)) %>%
    group_by(group) %>%
    summarise(mean_contamination = round(mean(contamination_fraction), 3),
              median_contamination = round(median(contamination_fraction), 3),
              n_cells = n(), .groups = "drop") %>%
    arrange(desc(mean_contamination))

  write.csv(summary_df, paste0(results_DECONTX_path, "Contamination_by", group_label, ".csv"), row.names = FALSE)

  p <- ggplot(df %>% filter(!is.na(contamination_fraction)),
               aes(x = reorder(group, contamination_fraction, FUN = median), y = contamination_fraction)) +
    geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.3, fill = "steelblue", alpha = 0.5) +
    coord_flip() +
    theme_bw() + theme(panel.grid = element_blank()) +
    labs(title = paste0("Estimated contamination fraction, by ", group_label),
         x = NULL, y = "Contamination fraction")
  ggsave(paste0(results_DECONTX_path, "Contamination_by", group_label, "_boxplot.png"), p,
         width = 8, height = max(5, length(unique(df$group)) * 0.3), dpi = 300, bg = "white")

  return(summary_df)
}

contam_by_sample <- summarise_contamination(dwAnnotated_joined, "orig.ident", "Sample")
contam_by_cluster <- summarise_contamination(dwAnnotated_joined, "seurat_clusters", "Cluster")
contam_by_celltype <- summarise_contamination(dwAnnotated_joined, "celltype", "Celltype")

cat("\n By sample:\n"); print(contam_by_sample)
cat("\n By cluster (pre-DecontX):\n"); print(contam_by_cluster)
cat("\n By celltype:\n"); print(contam_by_celltype)