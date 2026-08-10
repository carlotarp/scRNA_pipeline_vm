
fibrocyte_markers <- list(
  Pan_immune       = c("PTPRC"),                          # CD45 - el marcador definitorio, decisivo
  Fibrocyte_core   = c("COL1A1", "COL3A1", "SPARC", "EGR1"),
  Fibrocyte_specific = c("S100A4", "TM4SF1", "CD34"),
  Macrophage_residual = c("CD68", "ADGRE1", "LYZ"),        # esperado: presente pero MENOR que en macrófagos reales
  APC_function     = c("CD80", "CD86", "CD274")            # costimulación/PD-L1, relevante si es un fibrocito "tipo FC3"
)

plot_marker_dotplot(dwLeukocytes, marker_groups = fibrocyte_markers,
                       results_path = results_GEMX_LCA_path,
                       filename = paste0("Dotplot_Fibrocytes_Leukos_byClusters.png"),
                       group_by = "seurat_clusters")
plot_marker_dotplot(dwAnnotated, marker_groups = fibrocyte_markers,
                       results_path = results_GEMX_LCA_path,
                       filename = paste0("Dotplot_Fibrocytes_All_byClusters.png"),
                       group_by = "seurat_clusters")

png(paste0(results_GEMX_LCA_path, "VlnPlot_All_byCluster.png"), width = 1200, height = 600)
p <- VlnPlot(dwAnnotated, features = c("PTPRC", "CD34", "TM4SF1", "S100A4"), group.by = "seurat_clusters", ncol = 2)
print(p)
dev.off()

png(paste0(results_GEMX_LCA_path, "VlnPlot_All_byCelltype.png"), width = 1200, height = 600)
p <- VlnPlot(dwAnnotated, features = c("PTPRC", "CD34", "TM4SF1", "S100A4"), group.by = "celltype", ncol = 2)
print(p)
dev.off()

png(paste0(results_GEMX_LCA_path, "VlnPlot_Leukos_byCluster.png"), width = 1200, height = 600)
p <- VlnPlot(dwLeukocytes, features = c("PTPRC", "CD34", "TM4SF1", "S100A4"), group.by = "seurat_clusters", ncol = 2)
print(p)
dev.off()


## ============================================================
## 1. PTPRC absoluto: cluster 4 (Leukocytes) vs. CAFs genuinos (Stromal)
## ============================================================

# pct.exp (no color/z-score) de PTPRC en el cluster 4 sospechoso
ptprc_cluster4 <- FetchData(dwLeukocytes, vars = "PTPRC", cells = WhichCells(dwLeukocytes, idents = "6"))
pct_exp_cluster4 <- mean(ptprc_cluster4$PTPRC > 0) * 100

# pct.exp de PTPRC en cada cluster de Stromal (tus CAFs genuinos, referencia negativa)
dwStromal <- subset(dwAnnotated, subset = lineage == "Stromal")
stromal_dot_data <- DotPlot(dwStromal, features = "PTPRC", group.by = "seurat_clusters")$data

cat(paste0("PTPRC en cluster 4 (Leukocytes): ", round(pct_exp_cluster4, 1), "% de células positivas\n\n"))
cat("PTPRC por cluster en Stromal (referencia negativa esperada):\n")
print(stromal_dot_data[, c("id", "pct.exp", "avg.exp")])

# Visual directo, lado a lado
p1 <- VlnPlot(dwLeukocytes, features = "PTPRC", idents = "6", pt.size = 0.1) +
  ggtitle("PTPRC - cluster 4 (Leukocytes)")
p2 <- VlnPlot(dwStromal, features = "PTPRC", group.by = "seurat_clusters", pt.size = 0.1) +
  ggtitle("PTPRC - todos los clusters Stromal")

combined <- p1 | p2
ggsave(paste0(results_GEMX_LCA_path, "PTPRC_cluster4_vs_Stromal.png"), combined,
       width = 14, height = 5, dpi = 300, bg = "white")



## ============================================================
## 2. Sub-clustering del cluster 4, buscando subestructura FC-like
##    (matrix-like vs immune/APC-like, como FC1/FC3 del paper)
## ============================================================

# Aislar solo el cluster 4 y re-clusterizar a mayor resolución
dwCluster4 <- subset(dwLeukocytes, idents = "6")
dwCluster4 <- JoinLayers(dwCluster4)

dwCluster4 <- NormalizeData(dwCluster4)
dwCluster4 <- FindVariableFeatures(dwCluster4, selection.method = "mean.var.plot")
dwCluster4 <- ScaleData(dwCluster4)
dwCluster4 <- RunPCA(dwCluster4, npcs = min(20, ncol(dwCluster4) - 1))

# co3 propio de este sub-subset (mismo criterio que ya usamos)
pct <- dwCluster4[["pca"]]@stdev / sum(dwCluster4[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = TRUE)[1] + 1
co3_c4 <- min(co1, co2, na.rm = TRUE)

dwCluster4 <- FindNeighbors(dwCluster4, dims = 1:co3_c4)
dwCluster4 <- FindClusters(dwCluster4, resolution = 0.4)
dwCluster4 <- RunUMAP(dwCluster4, dims = 1:co3_c4)

# Dotplot con el mismo panel de fibrocitos, ahora a resolución fina
fc_subtype_markers <- list(
  Matrix_like   = c("COL1A1", "COL3A1", "SPARC", "FAP", "POSTN"),
  Immune_APC    = c("CD80", "CD86", "CD274", "PTPRC"),
  Proliferating = c("MKI67", "TOP2A"),
  Fibrocyte_core = c("EGR1", "TM4SF1", "S100A4", "CD34")
)

plot_marker_dotplot(dwCluster4, marker_groups = fc_subtype_markers,
                      results_path = results_GEMX_LCA_path,
                      filename = "Dotplot_Cluster4_FCsubtypes.png",
                      group_by = "seurat_clusters",
                      title = "Sub-clustering of cluster 4 - FC-like subtypes")

plot_dimplot(dwCluster4, reduction = "umap", group_by = "seurat_clusters",
             results_path = results_GEMX_LCA_path,
             filename = "DimPlot_Cluster4_subclusters.png", label = TRUE)

cat(paste0("Sub-clusters encontrados dentro del cluster 4: ", length(unique(dwCluster4$seurat_clusters)), "\n"))
print(table(dwCluster4$seurat_clusters))


##########################################################################333

##
##  Compare Fibrocyte Annotation Across 3 nFeature Upper-Cutoff Versions
##  (7500 / 5500 / 3500) - are the SAME cells being called fibrocytes,
##  or does the population change when the cutoff changes?
##

library(Seurat)
library(dplyr)
library(ggplot2)

project_path <- "/home/usuario/PROJECTS/260724_victor_scRNA/"
wd <- paste0(project_path, "codes/scRNA_pipeline/")
setwd(wd)
results_path <- paste0(project_path, "results/")
results_COMPARE_path <- paste0(results_path, "GEMX/CellAnnotation/Fibrocyte")
dir.create(results_COMPARE_path, recursive = TRUE, showWarnings = FALSE)

# Set paths to each version's fully_annotated_data.rds (must already have
# "celltype" with "Fibrocyte" assigned)
cutoff_runs <- list(
  "cutoff_7500" = paste0(results_path, "GEMX/CellAnnotation/7500/notumor_annotated_data.rds"),
  "cutoff_5500" = paste0(results_path, "GEMX/CellAnnotation/5500/notumor_annotated_data.rds"),
  "cutoff_3500" = paste0(results_path, "GEMX/CellAnnotation/3500/notumor_annotated_data.rds")
)

objects_list <- lapply(cutoff_runs, readRDS)
names(objects_list) <- names(cutoff_runs)
cat("\n All 3 cutoff versions loaded \n")


# ---------------------------------------------------------------
# 1. How many fibrocytes per version
# ---------------------------------------------------------------
fib_cells <- lapply(objects_list, function(obj) {
  colnames(obj)[!is.na(obj$celltype) & obj$celltype == "Fibrocyte"]
})

n_fib_summary <- data.frame(
  cutoff_run = names(fib_cells),
  n_fibrocytes = sapply(fib_cells, length),
  n_total_cells = sapply(objects_list, ncol)
)
n_fib_summary$pct_of_dataset <- round(n_fib_summary$n_fibrocytes / n_fib_summary$n_total_cells * 100, 3)

write.csv(n_fib_summary, paste0(results_COMPARE_path, "Fibrocyte_counts_by_cutoff.csv"), row.names = FALSE)
print(n_fib_summary)


# ---------------------------------------------------------------
# 2. Barcode overlap between the 3 fibrocyte sets (manual 3-way Venn counts)
# ---------------------------------------------------------------
A <- fib_cells[["cutoff_7500"]]
B <- fib_cells[["cutoff_5500"]]
C <- fib_cells[["cutoff_3500"]]

venn_counts <- data.frame(
  region = c("7500 only", "5500 only", "3500 only",
             "7500 & 5500 only", "7500 & 3500 only", "5500 & 3500 only",
             "all three"),
  n_cells = c(
    length(setdiff(A, union(B, C))),
    length(setdiff(B, union(A, C))),
    length(setdiff(C, union(A, B))),
    length(setdiff(intersect(A, B), C)),
    length(setdiff(intersect(A, C), B)),
    length(setdiff(intersect(B, C), A)),
    length(intersect(intersect(A, B), C))
  )
)
write.csv(venn_counts, paste0(results_COMPARE_path, "Fibrocyte_barcode_overlap.csv"), row.names = FALSE)
print(venn_counts)

cat(paste0("\n Jaccard 7500 vs 5500: ", round(length(intersect(A,B))/length(union(A,B)), 3), "\n"))
cat(paste0(" Jaccard 5500 vs 3500: ", round(length(intersect(B,C))/length(union(B,C)), 3), "\n"))
cat(paste0(" Jaccard 7500 vs 3500: ", round(length(intersect(A,C))/length(union(A,C)), 3), "\n"))


# ---------------------------------------------------------------
# 3. Cross-tabulation: for cells present in ALL 3 versions (regardless of
#    label), how consistently are they labeled Fibrocyte in each?
#    This isolates "does the LABEL change" from "does the CELL survive QC".
# ---------------------------------------------------------------
common_barcodes <- Reduce(intersect, lapply(objects_list, colnames))
cat(paste0("\n ", length(common_barcodes), " cells survived all 3 cutoffs (present in all 3 objects) \n"))

celltype_by_version <- data.frame(
  cell = common_barcodes,
  celltype_7500 = as.character(objects_list[["cutoff_7500"]]$celltype[common_barcodes]),
  celltype_5500 = as.character(objects_list[["cutoff_5500"]]$celltype[common_barcodes]),
  celltype_3500 = as.character(objects_list[["cutoff_3500"]]$celltype[common_barcodes])
)

# Keep only cells that were Fibrocyte in AT LEAST ONE version - the ones we care about
ever_fibrocyte <- celltype_by_version[
  celltype_by_version$celltype_7500 == "Fibrocyte" |
  celltype_by_version$celltype_5500 == "Fibrocyte" |
  celltype_by_version$celltype_3500 == "Fibrocyte", ]

write.csv(ever_fibrocyte, paste0(results_COMPARE_path, "Fibrocyte_label_consistency.csv"), row.names = FALSE)

cat("\n Cross-tabulation of labels for cells ever called Fibrocyte (present in all 3 cutoffs):\n")
print(table(ever_fibrocyte$celltype_7500, ever_fibrocyte$celltype_5500,
            dnn = c("7500", "5500")))
print(table(ever_fibrocyte$celltype_5500, ever_fibrocyte$celltype_3500,
            dnn = c("5500", "3500")))

n_stable <- sum(ever_fibrocyte$celltype_7500 == "Fibrocyte" &
                ever_fibrocyte$celltype_5500 == "Fibrocyte" &
                ever_fibrocyte$celltype_3500 == "Fibrocyte")
cat(paste0("\n ", n_stable, " of ", nrow(ever_fibrocyte),
           " ever-fibrocyte cells are labeled Fibrocyte in ALL 3 versions (",
           round(n_stable/nrow(ever_fibrocyte)*100, 1), "%)\n"))


# ---------------------------------------------------------------
# 4. QC profile of fibrocytes per version - did the identity shift the
#    kind of cell being captured as cutoff loosens?
# ---------------------------------------------------------------
qc_by_version <- do.call(rbind, lapply(names(objects_list), function(v) {
  obj <- objects_list[[v]]
  df <- obj@meta.data[fib_cells[[v]], c("nFeature_RNA", "nCount_RNA", "percent.mt")]
  df$cutoff_run <- v
  df
}))

p_qc <- ggplot(qc_by_version, aes(x = cutoff_run, y = nFeature_RNA)) +
  geom_violin(fill = "steelblue", alpha = 0.5) + geom_jitter(size = 0.3, alpha = 0.3, width = 0.15) +
  theme_bw() + theme(panel.grid = element_blank()) +
  labs(title = "nFeature_RNA of Fibrocyte cells by cutoff version")
ggsave(paste0(results_COMPARE_path, "Fibrocyte_nFeature_by_cutoff.png"), p_qc,
       width = 7, height = 5, dpi = 300, bg = "white")

write.csv(qc_by_version %>% group_by(cutoff_run) %>%
            summarise(mean_nFeature = mean(nFeature_RNA), median_nFeature = median(nFeature_RNA),
                      mean_nCount = mean(nCount_RNA), mean_pctmt = mean(percent.mt)),
          paste0(results_COMPARE_path, "Fibrocyte_QC_summary_by_cutoff.csv"), row.names = FALSE)


# ---------------------------------------------------------------
# 5. Marker profile comparison - is the "Fibrocyte" identity signature
#    (COL1A1/PTPRC/etc.) the same across versions, or diluted/shifted?
# ---------------------------------------------------------------
fib_markers_check <- c("PTPRC", "COL1A1", "COL3A1", "SPARC", "CD68", "LYZ", "TM4SF1", "S100A4")

marker_profile <- do.call(rbind, lapply(names(objects_list), function(v) {
  obj <- objects_list[[v]]
  obj_fib <- subset(obj, cells = fib_cells[[v]])
  obj_fib <- JoinLayers(obj_fib)   # <- fix: une counts.SC7b, counts.SC8, etc. en una sola capa

  markers_present <- fib_markers_check[fib_markers_check %in% rownames(obj_fib)]
  avg <- AverageExpression(obj_fib, features = markers_present, assay = "RNA")$RNA
  pct <- rowMeans(GetAssayData(obj_fib, layer = "counts")[markers_present, , drop = FALSE] > 0) * 100
  data.frame(cutoff_run = v, gene = markers_present, avg_expr = as.numeric(avg), pct_expr = pct)
}))

write.csv(marker_profile, paste0(results_COMPARE_path, "Fibrocyte_marker_profile_by_cutoff.csv"), row.names = FALSE)

p_markers <- ggplot(marker_profile, aes(x = gene, y = pct_expr, fill = cutoff_run)) +
  geom_col(position = "dodge") +
  theme_bw() + theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Marker % expressed in Fibrocyte cells, by cutoff version", y = "% cells expressing")
ggsave(paste0(results_COMPARE_path, "Fibrocyte_marker_comparison.png"), p_markers,
       width = 9, height = 5, dpi = 300, bg = "white")

cat(paste("\n ---- FINISHED FIBROCYTE CUTOFF COMPARISON ----
    Generated files:
      · Fibrocyte_counts_by_cutoff.csv
      · Fibrocyte_barcode_overlap.csv
      · Fibrocyte_label_consistency.csv
      · Fibrocyte_QC_summary_by_cutoff.csv
      · Fibrocyte_marker_profile_by_cutoff.csv
    Generated plots:
      · Fibrocyte_nFeature_by_cutoff.png
      · Fibrocyte_marker_comparison.png
          "))