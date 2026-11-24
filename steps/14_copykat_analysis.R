#!/usr/bin/env Rscript
##
## 14_copykat_analysis.R - attach CopyKAT calls, CNA matrix, CIN/PGA scores and a CNA-based UMAP
## Usage : Rscript pipeline/steps/14_copykat_analysis.R
## Reads : objects/notumor_annotated.rds ; 13_copykat_run/copykat_reference_cells.txt, 13_copykat_run/<S>/
## Writes: objects/copykat.rds (metadata copykat_prediction, CIN_score, PGA_score; assay copykat_cna;
##         reductions pca_cna, umap_cna) ; 14_copykat_analysis/*.png (all cells) and tumour/*.png, *.csv
##         (cells of lineage Tumor, by global cluster)
## Checkpoint: none
##
## Each sample's calls keep its observed cells and its own reference cells; reference cells of other samples
## (pooled reference) are left out. Samples without CopyKAT output get NA.
## CIN_score = SD of the cell's CNA profile; PGA_score = % of bins with |CNA| > copykat.pga_threshold.
##
file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
PIPELINE_DIR <- normalizePath(file.path(dirname(file_arg), ".."))
source(file.path(PIPELINE_DIR, "R/setup.R"))
out_dir <- step_dir("14_copykat_analysis")
ck_dir  <- step_dir("13_copykat_run")

obj <- read_run_object("notumor_annotated")
ref_cells <- readLines(paste0(ck_dir, "copykat_reference_cells.txt"))
cat("  Reference cells (all samples):", length(ref_cells), "\n")
samples <- sort(unique(as.character(obj$orig.ident)))
ref_sample <- setNames(as.character(obj$orig.ident[ref_cells]), ref_cells)

# ---- read every sample's prediction + CNA bins ----
pred_list <- list(); cna_list <- list()
for (s in samples) {
  sp <- paste0(ck_dir, s)
  if (!file.exists(file.path(sp, paste0(s, "_copykat_prediction.txt")))) {
    cat("  [WARNING] no CopyKAT output for", s, "- left out (NA)\n"); next
  }
  o <- read_copykat_outputs(sp, s)
  p <- o$prediction; rownames(p) <- p$cell.names
  own_ref <- names(ref_sample)[ref_sample == s]
  pred_list[[s]] <- p[!rownames(p) %in% ref_cells | rownames(p) %in% own_ref, ]
  cna_list[[s]] <- o$cna
}
cat("  CopyKAT outputs read for", length(pred_list), "samples\n")
all_pred <- do.call(rbind, unname(pred_list))

matched <- intersect(colnames(obj), rownames(all_pred))
obj$copykat_prediction <- NA
obj$copykat_prediction[matched] <- all_pred[matched, "copykat.pred"]

# ---- CNA matrix (bins shared by every sample x cells) -> assay copykat_cna ----
common_bins <- Reduce(intersect, lapply(cna_list, function(x) paste0(x$chrom, "_", x$chrompos)))
cna_all <- do.call(cbind, lapply(names(cna_list), function(s) {
  x <- cna_list[[s]]
  x <- x[match(common_bins, paste0(x$chrom, "_", x$chrompos)), ]
  cols <- setdiff(colnames(x), c("chrom", "chrompos", "abspos"))
  own_ref <- names(ref_sample)[ref_sample == s]
  keep <- c(cols[!cols %in% ref_cells], cols[cols %in% own_ref])
  as.matrix(x[, keep, drop = FALSE])
}))
rownames(cna_all) <- paste0("chr", common_bins)
cat(sprintf("  CNA matrix: %d bins x %d cells\n", nrow(cna_all), ncol(cna_all)))

cna_full <- matrix(NA_real_, nrow = nrow(cna_all), ncol = ncol(obj), dimnames = list(rownames(cna_all), colnames(obj)))
m <- intersect(colnames(obj), colnames(cna_all))
cna_full[, m] <- cna_all[, m]
obj[["copykat_cna"]] <- CreateAssayObject(data = cna_full)
rm(cna_full, cna_all); invisible(gc())

# ---- CIN / PGA ----
cna_data <- GetAssayData(obj, assay = "copykat_cna", layer = "data")
pga_thr <- cfg_get("copykat.pga_threshold")
obj$CIN_score <- apply(cna_data, 2, sd, na.rm = TRUE)[colnames(obj)]
obj$PGA_score <- apply(cna_data, 2, function(x) mean(abs(x) > pga_thr, na.rm = TRUE) * 100)[colnames(obj)]
obj$CIN_score[!colnames(obj) %in% m] <- NA   # no CNA profile
obj$PGA_score[!colnames(obj) %in% m] <- NA
print(summary(obj$CIN_score))

ggsave(paste0(out_dir, "Boxplot_CINscore_byCopyKATPrediction.png"),
       ggplot(obj@meta.data[!is.na(obj$CIN_score), ], aes(copykat_prediction, CIN_score, fill = copykat_prediction)) +
         geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.3) + theme_bw() + theme(panel.grid = element_blank()) +
         labs(title = "Chromosomal instability (CIN) score by CopyKAT prediction", x = "CopyKAT prediction",
              y = "CIN score (SD of CNA profile)"),
       width = 7, height = 5, dpi = 300, bg = "white")

# CIN/PGA by celltype (facet lineage) and by sample (facet clinical Subtype)
plot_score_by_group <- function(object, score_col, group_col, filename, title, y_lab, facet_col, colors) {
  meta <- object@meta.data[!is.na(object@meta.data[[score_col]]), ]
  p <- ggplot(meta, aes(x = .data[[group_col]], y = .data[[score_col]], fill = .data[[group_col]])) +
    geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.3) + theme_bw() +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
    labs(title = title, x = group_col, y = y_lab) + scale_fill_manual(values = colors) +
    facet_wrap(vars(.data[[facet_col]]), scales = "free_x")
  ggsave(paste0(out_dir, filename), p, width = 10, height = 6, dpi = 300, bg = "white")
}
md <- obj@meta.data
ct_cols <- md %>% distinct(celltype, lineage) %>% mutate(across(everything(), as.character)) %>%
  { setNames(lineage_colors[.$lineage], .$celltype) }
sa_cols <- md %>% distinct(orig.ident, Subtype) %>% mutate(across(everything(), as.character)) %>%
  { setNames(pam50_subtype_colors[.$Subtype], .$orig.ident) }
for (sc in c("CIN_score", "PGA_score")) {
  nm <- sub("_score", "", sc)
  yl <- if (sc == "CIN_score") "CIN score (SD of CNA profile)" else "PGA score (% bins altered)"
  plot_score_by_group(obj, sc, "celltype", paste0("Boxplot_", nm, "score_byCelltype.png"), paste(nm, "by celltype"), yl, "lineage", ct_cols)
  plot_score_by_group(obj, sc, "orig.ident", paste0("Boxplot_", nm, "score_bySample.png"), paste(nm, "by sample"), yl, "Subtype", sa_cols)
}

# ---- UMAP built from the CNA matrix (cells without a profile: NA) ----
n_cna <- cfg_get("copykat.cna_umap_pcs")
cells_cna <- colnames(cna_data)[colSums(!is.na(cna_data)) > 0]
oc <- subset(obj, cells = cells_cna)
DefaultAssay(oc) <- "copykat_cna"
VariableFeatures(oc) <- rownames(oc[["copykat_cna"]])
oc <- ScaleData(oc, features = VariableFeatures(oc))
oc <- RunPCA(oc, features = VariableFeatures(oc), npcs = n_cna, reduction.name = "pca_cna", reduction.key = "PCACNA_")
oc <- run_umap_pca(oc, reduction = "pca_cna", dims = 1:n_cna, reduction.name = "umap_cna")
pad_reduction <- function(red, all_cells) {
  emb <- red@cell.embeddings
  full <- matrix(NA_real_, nrow = length(all_cells), ncol = ncol(emb), dimnames = list(all_cells, colnames(emb)))
  full[rownames(emb), ] <- emb
  CreateDimReducObject(embeddings = full, key = red@key, assay = "copykat_cna")
}
obj[["pca_cna"]]  <- pad_reduction(oc[["pca_cna"]], colnames(obj))
obj[["umap_cna"]] <- pad_reduction(oc[["umap_cna"]], colnames(obj))
rm(oc); invisible(gc())
for (g in c("copykat_prediction", "celltype", "lineage", "orig.ident", "Subtype")) {
  fn <- c(copykat_prediction = "CopyKATPrediction", celltype = "Celltype", lineage = "Lineage",
          orig.ident = "Sample", Subtype = "ClinicalSubtype")[[g]]
  plot_dimplot(obj, reduction = "umap_cna", group_by = g, results_path = out_dir, filename = paste0("DimPlot_UMAP_CNA_", fn, ".png"))
}

vals <- as.vector(cna_data); vals <- vals[!is.na(vals)]
cat("  CNA values: sd =", round(sd(vals), 4), "\n"); print(summary(vals))
ggsave(paste0(out_dir, "Histogram_CNA_values.png"),
       ggplot(data.frame(value = vals), aes(value)) + geom_histogram(bins = 100, fill = "steelblue", color = NA) +
         geom_vline(xintercept = 0, linetype = "dashed", color = "firebrick") + theme_bw() + theme(panel.grid = element_blank()) +
         labs(title = "Distribution of CopyKAT CNA values (all bins x cells)", x = "CNA value", y = "Count"),
       width = 7, height = 5, dpi = 300, bg = "white")
rm(vals, cna_data); invisible(gc())

plot_dimplot(obj, reduction = "umap", group_by = "copykat_prediction", results_path = out_dir,
             filename = "DimPlot_UMAP_CopyKATPrediction.png")
plot_composition(obj, cluster_col = "orig.ident", group_by = "copykat_prediction", results_path = out_dir,
                 filename = "BarPlot_CopyKATPrediction_bySample.png")
plot_composition(obj, cluster_col = "seurat_clusters", group_by = "copykat_prediction", results_path = out_dir,
                 filename = "BarPlot_CopyKATPrediction_byCluster.png")
plot_composition(obj, cluster_col = "celltype", group_by = "copykat_prediction", results_path = out_dir,
                 filename = "BarPlot_CopyKATPrediction_byCelltype.png")

saveRDS(obj, run_object("copykat"))
cat("  ->", run_object("copykat"), "\n")

# ---- tumour-lineage CopyKAT plots (drawn by the tumour markers step in the original pipeline);
#      clusters = global clusters, the tumour sub-clusters do not exist yet at this point ----
tu_dir <- step_dir("14_copykat_analysis", "tumour")
tu <- subset(obj, cells = colnames(obj)[obj$lineage == "Tumor"])
rm(obj); invisible(gc())
tu$seurat_clusters <- droplevels(factor(tu$seurat_clusters))
plot_copykat_composition(tu, tu_dir)
plot_copykat_cna_boxplot(tu, tu_dir, score_col = "CIN_score", y_lab = "CIN score (SD of CNA profile)")
plot_copykat_cna_boxplot(tu, tu_dir, score_col = "PGA_score", y_lab = "PGA score (% bins altered)")
plot_dimplot(tu, reduction = "umap", group_by = "copykat_prediction", results_path = tu_dir,
             filename = "DimPlot_UMAP_CopyKATPrediction_Tumoral.png")
plot_composition(tu, cluster_col = "seurat_clusters", group_by = "copykat_prediction", results_path = tu_dir,
                 filename = "BarPlot_CopyKATPrediction_Tumor_byCluster.png")
for (g in c("copykat_prediction", "seurat_clusters", "orig.ident", "Subtype")) {
  fn <- c(copykat_prediction = "CopyKATPrediction", seurat_clusters = "SeuratClusters", orig.ident = "Sample",
          Subtype = "ClinicalSubtype")[[g]]
  plot_dimplot(tu, reduction = "umap_cna", group_by = g, results_path = tu_dir, filename = paste0("DimPlot_UMAP_CNA_", fn, "_Tumoral.png"))
}
tu_cna <- subset(tu, cells = colnames(tu)[!is.na(tu$CIN_score)])
ck_diff <- differential_scores_by_cluster(tu_cna, c("CIN_score", "PGA_score"))
write.csv(ck_diff, paste0(tu_dir, "Differential_Copykat_byCluster.csv"), row.names = FALSE)
plot_differential_barplot(paste0(tu_dir, "Differential_Copykat_byCluster.csv"), "score",
                          "Differential CopyKAT scores by cluster", "Barplot_Differential_Copykat.png", tu_dir)

cat("\n14_copykat_analysis done. Next: 15_tumour_copykat_evidence.R check (after 12_tumour_annotate.R)\n")
