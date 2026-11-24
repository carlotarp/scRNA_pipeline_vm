#!/usr/bin/env Rscript
##
## 15_tumour_copykat_evidence.R - CopyKAT evidence on the annotated tumour clusters + sub-zooms for the hand re-annotation
## Usage : Rscript pipeline/steps/15_tumour_copykat_evidence.R check
##         Rscript pipeline/steps/15_tumour_copykat_evidence.R zoom lp|lowrna
##   check       : CopyKAT call / "flat" CNA per tumour cluster and sample x cluster, epithelial-state genes by
##                 CopyKAT state (diploid_flat, diploid_CNV, aneuploid, not.defined, no call)
##   zoom lp     : sub-zoom of the tumour clusters that may hide normal epithelium (decision
##                 subclustering.tumour_lp.clusters)
##   zoom lowrna : sub-zoom of the low-RNA tumour clusters (decision subclustering.tumour_lowrna.clusters)
## Reads : objects/tumour_annotated.rds, objects/copykat.rds (check);
##         objects/tumour_annotated.rds + 15_tumour_copykat_evidence/check/tumour_cells_copykat.csv.gz (zoom)
## Writes: 15_tumour_copykat_evidence/check/*.csv, *.png
##         objects/tumour_<g>_grid.rds (all grid resolutions), objects/tumour_<g>.rds (chosen resolution);
##         15_tumour_copykat_evidence/zoom_tumour_<g>/ (resolution grid, composition, UMAPs, dot plot, markers,
##         top_expressed_genes.csv) ; annotation/tumour_<g>.tsv template
## Checkpoint: subclustering.tumour_<g>.clusters (cells to zoom into), subclustering.tumour_<g>.resolution;
##             then fill annotation/tumour_lp.tsv and annotation/tumour_lowrna.tsv (step 16)
##
## "flat" CNA = CIN_score <= copykat.flat_quantile of the same sample's T cells (copykat.ref_celltypes).
## subclustering.tumour_<g>.integrate: false = PCA only (small normal groups spread over patients can be dissolved by Harmony),
## true = Harmony per sample (generate_subset).
##
file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
PIPELINE_DIR <- normalizePath(file.path(dirname(file_arg), ".."))
source(file.path(PIPELINE_DIR, "R/setup.R"))
out_dir <- step_dir("15_tumour_copykat_evidence")

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || !(args[1] == "check" || (args[1] == "zoom" && length(args) >= 2 && args[2] %in% c("lp", "lowrna")))) {
  stop("Usage: 15_tumour_copykat_evidence.R check | zoom lp|lowrna")
}
check_dir <- step_dir("15_tumour_copykat_evidence", "check")
ck_table  <- paste0(check_dir, "tumour_cells_copykat.csv.gz")

# ======================================================================================
# check
# ======================================================================================
if (args[1] == "check") {
  ck <- read_run_object("copykat")
  md <- ck@meta.data[, c("orig.ident", "celltype", "lineage", "copykat_prediction", "CIN_score", "PGA_score", "nFeature_RNA", "nCount_RNA")]
  rm(ck); invisible(gc())
  q <- cfg_get("copykat.flat_quantile")
  thr <- md %>% filter(celltype %in% cfg_get("copykat.ref_celltypes"), !is.na(CIN_score)) %>%
    group_by(orig.ident) %>% summarise(cin_ref_T = quantile(CIN_score, q), .groups = "drop")
  write.csv(thr, paste0(check_dir, "flat_threshold_per_sample.csv"), row.names = FALSE)

  tu <- read_run_object("tumour_annotated")
  tm <- tu@meta.data[, c("orig.ident", "seurat_clusters", "celltype")]
  tm$cluster <- as.character(tm$seurat_clusters)
  tm <- cbind(tm, md[rownames(tm), c("copykat_prediction", "CIN_score", "PGA_score", "nFeature_RNA", "nCount_RNA")])
  tm <- left_join(tibble::rownames_to_column(tm, "cell"), thr, by = "orig.ident")
  tm$copykat_prediction[is.na(tm$copykat_prediction)] <- "no call"
  tm$flat <- !is.na(tm$CIN_score) & tm$CIN_score <= tm$cin_ref_T
  tm$state <- ifelse(tm$copykat_prediction == "diploid" & tm$flat, "diploid_flat",
                     ifelse(tm$copykat_prediction == "diploid", "diploid_CNV", tm$copykat_prediction))
  genes <- unique(unlist(load_marker_panel("epithelial_states", object = tu)))
  ex <- as.matrix(t(as.matrix(GetAssayData(tu, assay = "RNA", layer = "data")[genes, tm$cell]))) > 0
  tm <- cbind(tm, ex)
  write.csv(tm, gzfile(ck_table), row.names = FALSE)

  # 1. CopyKAT prediction and flat share per cluster, and per sample x cluster
  pc <- tm %>% count(cluster, copykat_prediction) %>% group_by(cluster) %>% mutate(frac = n / sum(n)) %>% ungroup()
  write.csv(tidyr::pivot_wider(pc %>% select(-frac), names_from = copykat_prediction, values_from = n, values_fill = 0),
            paste0(check_dir, "prediction_by_cluster.csv"), row.names = FALSE)
  pc$cluster <- factor(pc$cluster, levels = as.character(sort(as.integer(unique(pc$cluster)))))
  ggsave(paste0(check_dir, "Barplot_prediction_by_cluster.png"),
         ggplot(pc, aes(cluster, frac, fill = copykat_prediction)) + geom_col() + theme_bw() + labs(y = "fraction of cells"),
         width = 9, height = 4.5, bg = "white")
  fl <- tm %>% group_by(cluster, celltype) %>%
    summarise(cells = n(), diploid = mean(copykat_prediction == "diploid"), flat = mean(flat),
              diploid_flat = mean(state == "diploid_flat"), median_genes = median(nFeature_RNA), .groups = "drop")
  write.csv(fl, paste0(check_dir, "flat_by_cluster.csv"), row.names = FALSE)
  write.csv(tm %>% count(orig.ident, cluster, copykat_prediction) %>%
              tidyr::pivot_wider(names_from = copykat_prediction, values_from = n, values_fill = 0),
            paste0(check_dir, "prediction_by_sample_cluster.csv"), row.names = FALSE)
  write.csv(tm %>% count(cluster, orig.ident, diploid_flat = state == "diploid_flat") %>%
              tidyr::pivot_wider(names_from = diploid_flat, values_from = n, values_fill = 0, names_prefix = "diploid_flat_"),
            paste0(check_dir, "diploidflat_by_sample_cluster.csv"), row.names = FALSE)

  # 2. % of cells expressing each epithelial-state gene, per cluster x CopyKAT state
  grp <- tm %>% group_by(cluster, state) %>%
    summarise(cells = n(), median_genes = median(nFeature_RNA), across(all_of(genes), mean), .groups = "drop")
  write.csv(grp, paste0(check_dir, "markers_by_cluster_state.csv"), row.names = FALSE)

  # 3. tumour UMAP by CopyKAT call and by diploid + flat
  um <- as.data.frame(Embeddings(tu, "umap")[tm$cell, ]); colnames(um) <- c("u1", "u2")
  um <- cbind(um, tm[, c("cluster", "copykat_prediction", "state")])
  ggsave(paste0(check_dir, "UMAP_tumour_prediction.png"),
         ggplot(um, aes(u1, u2, colour = copykat_prediction)) + geom_point(size = 0.05) + theme_void() +
           guides(colour = guide_legend(override.aes = list(size = 3))), width = 7, height = 6, bg = "white")
  um$dflat <- um$state == "diploid_flat"
  ggsave(paste0(check_dir, "UMAP_tumour_diploidflat.png"),
         ggplot(um[order(um$dflat), ], aes(u1, u2, colour = dflat)) + geom_point(size = 0.05) +
           scale_colour_manual(values = c("grey80", "red"), name = "diploid + flat") + theme_void(), width = 7, height = 6, bg = "white")
  print(as.data.frame(fl))
  cat("\n15 check done. Record the zoom groups, e.g.\n",
      "  Rscript pipeline/decide.R subclustering.tumour_lp.clusters <clusters> \"<reason>\"\n",
      "  Rscript pipeline/decide.R subclustering.tumour_lowrna.clusters <clusters> \"<reason>\"\n",
      "then run 'zoom lp' and 'zoom lowrna'.\n")
  quit(save = "no")
}

# ======================================================================================
# zoom lp | lowrna (two phases: grid cache + evidence, then the chosen resolution)
# ======================================================================================
key  <- paste0("tumour_", args[2])
zdir <- step_dir("15_tumour_copykat_evidence", paste0("zoom_", key))
grid <- cfg_get("clustering.resolution_grid")
n_pcs <- cfg_get(paste0("subclustering.", key, ".n_pcs"))
cache <- run_object(paste0(key, "_grid"))
clusters <- sort(as.character(get_decision(paste0("subclustering.", key, ".clusters"))))

# Re-normalise and cluster without integration (PCA neighbours)
subset_no_integration <- function(object, cells, resolution, n_pcs) {
  z <- JoinLayers(subset(object, cells = cells))
  z <- NormalizeData(z); z <- FindVariableFeatures(z); z <- ScaleData(z); z <- RunPCA(z)
  z <- FindNeighbors(z, reduction = "pca", dims = 1:n_pcs)
  z <- run_umap_pca(z, reduction = "pca", dims = 1:n_pcs)
  FindClusters(z, resolution = resolution)
}

z <- if (file.exists(cache)) readRDS(cache) else NULL
if (!is.null(z) && !identical(sort(as.character(z@misc$zoom_clusters)), clusters)) {
  cat("  Cached zoom was built from other clusters: recomputing\n"); z <- NULL
}
if (is.null(z)) {
  if (!file.exists(ck_table)) stop("Missing ", ck_table, ": run '15_tumour_copykat_evidence.R check' first")
  ck <- read.csv(ck_table)
  tu <- read_run_object("tumour_annotated")
  cells <- ck$cell[as.character(ck$cluster) %in% clusters]
  cat(sprintf("  %s: %d cells from tumour clusters %s\n", key, length(cells), paste(clusters, collapse = ", ")))
  z <- if (isTRUE(cfg_get(paste0("subclustering.", key, ".integrate")))) {
    generate_subset(tu, cells, resolution = grid, label = key, n_pcs = n_pcs)
  } else {
    subset_no_integration(tu, cells, resolution = grid, n_pcs = n_pcs)
  }
  rm(tu); invisible(gc())
  m <- ck[match(colnames(z), ck$cell), ]
  z$state <- m$state; z$cluster_prev <- as.character(m$cluster)
  z@misc$zoom_clusters <- clusters
  saveRDS(z, cache)
  plot_resolution_grid(z, zdir, resolutions = grid)
  sizes <- do.call(rbind, lapply(grid, function(r) {
    tb <- table(z[[paste0("RNA_snn_res.", r)]][, 1])
    data.frame(resolution = r, cluster = names(tb), cells = as.integer(tb))
  }))
  write.csv(sizes, paste0(zdir, "cluster_sizes_per_resolution.csv"), row.names = FALSE)
  cat("  Grid cache ->", cache, "\n")
}

res <- get_decision(paste0("subclustering.", key, ".resolution"))
res_col <- paste0("RNA_snn_res.", res)
if (!res_col %in% colnames(z@meta.data)) stop(res_col, " not in the grid cache (clustering.resolution_grid)")
z$seurat_clusters <- factor(z[[res_col]][, 1], levels = sort(as.integer(unique(as.character(z[[res_col]][, 1])))))
Idents(z) <- "seurat_clusters"
cat(sprintf("  %s at resolution %s: %d clusters\n", key, res, nlevels(z$seurat_clusters)))

# composition per sub-cluster: sample, CopyKAT state, previous tumour cluster, depth
tab <- function(col) z@meta.data %>% count(seurat_clusters, .data[[col]]) %>%
  tidyr::pivot_wider(names_from = all_of(col), values_from = n, values_fill = 0)
write.csv(tab("orig.ident"), paste0(zdir, "sub_by_sample.csv"), row.names = FALSE)
write.csv(tab("state"), paste0(zdir, "sub_by_copykat.csv"), row.names = FALSE)
write.csv(tab("cluster_prev"), paste0(zdir, "sub_by_prev_cluster.csv"), row.names = FALSE)
write.csv(tab("celltype"), paste0(zdir, "sub_by_celltype.csv"), row.names = FALSE)
write.csv(z@meta.data %>% group_by(seurat_clusters) %>%
            summarise(cells = n(), median_genes = median(nFeature_RNA), median_counts = median(nCount_RNA), pct_mt = median(percent.mt)),
          paste0(zdir, "sub_depth.csv"), row.names = FALSE)

um <- data.frame(Embeddings(z, "umap")); colnames(um) <- c("u1", "u2")
um <- cbind(um, z@meta.data[, c("seurat_clusters", "orig.ident", "state")])
for (col in c("seurat_clusters", "orig.ident", "state")) {
  p <- ggplot(um, aes(u1, u2, colour = .data[[col]])) + geom_point(size = 0.1) + theme_void() +
    guides(colour = guide_legend(override.aes = list(size = 3)))
  if (col == "seurat_clusters") p <- p + geom_text(data = um %>% group_by(seurat_clusters) %>% summarise(u1 = median(u1), u2 = median(u2)),
                                                   aes(label = seurat_clusters), colour = "black")
  ggsave(paste0(zdir, "UMAP_", col, ".png"), p, width = 7, height = 6, bg = "white")
}

panel <- unique(unlist(load_marker_panel("epithelial_states", object = z)))
ggsave(paste0(zdir, "DotPlot_panel.png"), dot_plot(z, features = panel, group.by = "seurat_clusters") + RotatedAxis(),
       width = 15, height = 5 + 0.2 * nlevels(z$seurat_clusters), bg = "white")

mf <- find_markers_for_clusters(z, levels(z$seurat_clusters), zdir)
top_expressed_genes(z, "seurat_clusters", zdir)
saveRDS(z, run_object(key))
write_annotation_template(z, key, markers_df = mf)
cat(sprintf("\n15 zoom %s done -> %s\nFill annotation/%s.tsv (label KEEP = keep the cell's current label; DROP = remove), then run 16.\n",
            args[2], run_object(key), key))
