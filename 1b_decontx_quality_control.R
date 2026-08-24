##
##  Parallel script: DecontX ambient RNA correction, per sample
##  Does NOT modify seurat_1_QC.R - produces corrected count matrices that
##  can later be swapped in as the input to that script, if desired.
##
##  Unlike SoupX, DecontX does NOT need the raw/unfiltered CellRanger matrix
##  (no empty droplets required) - it estimates contamination directly from
##  the filtered matrix using a Bayesian mixture model, guided by cluster
##  labels. This is why we switched to it: raw_feature_bc_matrix is not
##  available for these samples.
##

library(Seurat)
library(celda)
library(dplyr)

project_path <- "/home/usuario/PROJECTS/260724_victor_scRNA/"
data_path <- "/home/usuario/DATASETS/scRNAseq/260106_carlota_GEMX/2026_HN00264849/allPool/"
results_path <- paste0(project_path, "results/")
results_DECONTX_path <- paste0(results_path, "GEMX/ContaminationCorrection/")
dir.create(results_DECONTX_path, recursive = TRUE, showWarnings = FALSE)

sample_list <- c("SC7b", "SC8", "SC9", "SC10", "SC11", "SC12", "SC13", "SC14",
                  "SC15", "SC16", "SC17", "SC18", "SC19", "SC20", "SC21", "SC5_SC22")

dwAnnotated <- readRDS(file.path(results_GEMX_CA_path, "fully_annotated_data.rds"))
dwAnnotated <- JoinLayers(dwAnnotated)
all_genes <- rownames(dwAnnotated)
corrected_full <- Matrix::Matrix(0, nrow = length(all_genes), ncol = ncol(dwAnnotated),
                                   sparse = TRUE, dimnames = list(all_genes, colnames(dwAnnotated)))

for (s in sample_list) {

  cat(paste0("\n\n ==== Sample: ", s, " ==== \n"))

  filtered_path <- paste0(data_path, s, "_filtered_feature_barcode_matrix/")
  toc <- Read10X(filtered_path)  # only the filtered matrix is needed

  # Quick, low-resolution clustering - DecontX uses cluster identity to guide
  # its contamination estimate (recommended, not strictly required). This is
  # NOT your final pipeline clustering, just a fast preliminary pass.
  srat_tmp <- CreateSeuratObject(toc)
  srat_tmp <- NormalizeData(srat_tmp, verbose = FALSE)
  srat_tmp <- FindVariableFeatures(srat_tmp, verbose = FALSE)
  srat_tmp <- ScaleData(srat_tmp, verbose = FALSE)
  srat_tmp <- RunPCA(srat_tmp, npcs = 20, verbose = FALSE)
  srat_tmp <- FindNeighbors(srat_tmp, dims = 1:20, verbose = FALSE)
  srat_tmp <- FindClusters(srat_tmp, resolution = 0.5, verbose = FALSE)

  cluster_labels <- as.integer(factor(srat_tmp$seurat_clusters))

  # Run DecontX
  decontx_res <- decontX(x = toc, z = cluster_labels)

  corrected_counts <- decontx_res$decontXcounts
  contamination_per_cell <- decontx_res$contamination

  saveRDS(corrected_counts, paste0(results_DECONTX_path, s, "_corrected_counts.rds"))

  sample_cells <- colnames(dwAnnotated)[dwAnnotated$orig.ident == s]
  raw_barcodes <- sub("_[0-9]+$", "", sample_cells)
  matched_idx <- match(raw_barcodes, colnames(corrected_s))
  shared_genes <- intersect(rownames(corrected_s), all_genes)
  corrected_full[shared_genes, sample_cells[valid]] <- corrected_s[shared_genes, matched_idx[valid]]

  sce <- SingleCellExperiment(assays = list(counts = toc))
  decontXcounts(sce) <- corrected_counts
  sce$decontX_clusters <- cluster_labels
  sce$decontX_contamination <- contamination_per_cell
  reducedDim(sce, "decontX_UMAP") <- decontx_res$estimates$decontX$UMAP

  # --- 1. DecontX's own clusters on its own UMAP ---
  p1 <- plotDimReduceCluster(x = sce$decontX_clusters,
                               dim1 = reducedDim(sce, "decontX_UMAP")[, 1],
                               dim2 = reducedDim(sce, "decontX_UMAP")[, 2])
  ggsave(paste0(results_path_plots, s, "_decontX_clusters_UMAP.png"), p1, width = 6, height = 5, dpi = 300, bg = "white")

  # --- 2. Per-cell contamination fraction on the UMAP - the most informative one ---
  p2 <- plotDecontXContamination(sce)
  ggsave(paste0(results_path_plots, s, "_decontX_contamination_UMAP.png"), p2, width = 6, height = 5, dpi = 300, bg = "white")

  # --- 3. % of marker genes detected, before vs after correction ---
  markers_list <- list(
    Collagen = collagen_genes,
    Tcell = intersect(c("CD3D", "CD3E"), rownames(toc)),
    Myeloid = intersect(c("LYZ", "CD68"), rownames(toc)),
    Bcell = intersect(c("CD79A", "MS4A1"), rownames(toc))
  )
  markers_list <- markers_list[sapply(markers_list, length) > 0]

  p3 <- plotDecontXMarkerPercentage(sce, markers = markers_list, assayName = c("counts", "decontXcounts"))
  ggsave(paste0(results_path_plots, s, "_decontX_markerPct_beforeAfter.png"), p3, width = 8, height = 6, dpi = 300, bg = "white")

  # --- 4. Raw expression of individual marker genes, before vs after ---
  p4 <- plotDecontXMarkerExpression(sce, unlist(markers_list))
  ggsave(paste0(results_path_plots, s, "_decontX_markerExpr_beforeAfter.png"), p4, width = 8, height = 6, dpi = 300, bg = "white")

  cat(paste0("\n Sample ", s, " done -> ", results_DECONTX_path, s, "_corrected_counts.rds \n"))
}
dwAnnotated[["RNA_decontX"]] <- CreateAssayObject(counts = corrected_full)
dwAnnotated <- NormalizeData(dwAnnotated, assay = "RNA_decontX", verbose = FALSE)
saveRDS(dwAnnotated, paste0(results_DECONTX_path, "decontx_corrected_data.rds"))

cat("\n ---- FINISHED DecontX CORRECTION (all samples) ---- \n")