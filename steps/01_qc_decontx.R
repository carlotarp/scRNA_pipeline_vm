#!/usr/bin/env Rscript
##
## 01_qc_decontx.R - DecontX (raw droplets as ambient background) BEFORE QC, cell QC, doublets,
##                   merge on disk (BPCells), Harmony per sample, optional clinical annotation.
## Usage : Rscript pipeline/steps/01_qc_decontx.R sample <S>   (one process per sample, run concurrently)
##         Rscript pipeline/steps/01_qc_decontx.R merge        (after every sample finished)
## Reads : samples.csv raw_dir (all droplets) + filtered_dir (CellRanger called cells); clinical_sheet (optional)
## Writes: sample -> 01_qc_decontx/per_sample/<S>.rds + per-sample QC plots
##         merge  -> objects/qc_merged.rds ; <bpcells_dir>/<run_name>/<S>/ (on-disk counts) ;
##                   01_qc_decontx/ (cell_filtering.csv, QC plots) ; 01_qc_decontx/decontx/ (contamination, gene change)
## Checkpoint: none
##
## Per sample: quick clustering -> decontX(background = raw droplets) -> rounded counts ->
##   percent.mt < qc.max_percent_mt, qc.min_features < nFeature < qc.max_features,
##   nCount / nFeature outliers at qc.nmads MADs (log scale, per sample) -> scDblFinder singlets.
##
file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
PIPELINE_DIR <- normalizePath(file.path(dirname(file_arg), ".."))
source(file.path(PIPELINE_DIR, "R/setup.R"))
suppressPackageStartupMessages({
  library(celda)
  library(SingleCellExperiment)
  library(scDblFinder)
})
out_dir <- step_dir("01_qc_decontx")
per_sample_dir <- step_dir("01_qc_decontx", "per_sample")
decontx_dir <- step_dir("01_qc_decontx", "decontx")

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[1] else ""
if (!mode %in% c("sample", "merge") || (mode == "sample" && length(args) < 2))
  stop("Usage: Rscript pipeline/steps/01_qc_decontx.R sample <sample> | merge")
if (cfg_get("integration.method") != "harmony") stop("01_qc_decontx: only integration.method = harmony is implemented")

samples <- load_samples(check_dirs = mode == "sample")   # merge only reads the per-sample results
round_counts <- isTRUE(cfg_get("decontx.round_counts"))
mito_cut <- cfg_get("qc.max_percent_mt")
feature_lower_cut <- cfg_get("qc.min_features")
feature_upper_cut <- cfg_get("qc.max_features")
qc_nmads <- cfg_get("qc.nmads")

# ---------------------------------------------------------------
# DecontX on the CellRanger-filtered cells of one sample, with the raw droplets as background
# (decontX drops the barcodes of x from the background). The filtered matrix lacks the probes
# Flex excludes (~900 genes), so the raw droplets are restricted to the filtered genes.
# Returns corrected counts, uncorrected counts and per-cell contamination.
# ---------------------------------------------------------------
run_decontx_sample <- function(name, raw_dir, filtered_dir) {
  cat(paste0("\n\n ==== Sample: ", name, " ==== \n"))
  tod <- Read10X(data.dir = raw_dir)
  filtered <- Read10X(data.dir = filtered_dir)
  stopifnot(all(rownames(filtered) %in% rownames(tod)), all(colnames(filtered) %in% colnames(tod)))
  tod <- tod[rownames(filtered), ]
  toc <- tod[, colnames(filtered)]
  rm(filtered)
  tod <- tod[, Matrix::colSums(tod) > 0]   # zero-count droplets add nothing to the ambient profile
  cat(paste0(" CellRanger-filtered cells: ", ncol(toc), " | raw droplets with counts: ", ncol(tod), "\n"))

  # Quick per-sample clustering to guide DecontX
  dims <- seq_len(cfg_get("decontx.cluster_dims"))
  srat_tmp <- CreateSeuratObject(toc)
  srat_tmp <- NormalizeData(srat_tmp, verbose = FALSE)
  srat_tmp <- FindVariableFeatures(srat_tmp, verbose = FALSE)
  srat_tmp <- ScaleData(srat_tmp, verbose = FALSE)
  srat_tmp <- RunPCA(srat_tmp, npcs = max(dims), verbose = FALSE)
  srat_tmp <- FindNeighbors(srat_tmp, dims = dims, verbose = FALSE)
  srat_tmp <- FindClusters(srat_tmp, resolution = cfg_get("decontx.cluster_res"), verbose = FALSE)
  cluster_labels <- as.integer(factor(srat_tmp$seurat_clusters))
  rm(srat_tmp)

  decontx_res <- decontX(x = toc, z = cluster_labels, background = tod)
  rm(tod); invisible(gc())
  adjusted <- decontx_res$decontXcounts
  if (round_counts) adjusted <- round(adjusted)   # integer counts for QC metrics and scDblFinder
  adjusted <- Matrix::drop0(as(adjusted, "CsparseMatrix"))

  contamination <- data.frame(sample = name, cellranger_cells = ncol(toc),
                              mean_contamination = mean(decontx_res$contamination),
                              median_contamination = median(decontx_res$contamination))
  cat(paste0(" Mean estimated contamination: ", round(mean(decontx_res$contamination), 3), "\n"))
  list(adjusted = adjusted, raw = toc, contamination = setNames(decontx_res$contamination, colnames(toc)),
       contamination_summary = contamination)
}

# ---------------------------------------------------------------
# QC of one sample on the DecontX-corrected counts. Appends one row per filter to `counts`
# (data.frame Cell_counts, Feature_discard, SampleID). Returns the singlets, the counts table
# and the pre-removal doublet metadata.
# ---------------------------------------------------------------
process_sample <- function(name, tumor, counts, results_path) {
  cat(paste("\n Processing Sample", name, "...\n"))
  log_count <- function(n, step) rbind(counts, data.frame(Cell_counts = n, Feature_discard = step, SampleID = name))

  tumor[["log10_nCount_RNA"]] <- log10(tumor$nCount_RNA)
  tumor[["percent.mt"]] <- PercentageFeatureSet(object = tumor, pattern = "^MT-")
  tumor[["percent.hb"]] <- PercentageFeatureSet(object = tumor, pattern = "^HB[^(P|E|S)]")

  # Pre-filter distributions
  plot_qc_vlnplot(tumor, name = name, results_path = results_path)
  plot_mito_vs_nfeature(tumor, name = name, results_path = results_path, mito_cut = mito_cut,
                        feature_lower_cut = feature_lower_cut, feature_upper_cut = feature_upper_cut)

  sample <- subset(x = tumor, subset = percent.mt < mito_cut)
  counts <- log_count(ncol(sample), "Mitoch")
  sample <- subset(x = sample, subset = nFeature_RNA > feature_lower_cut & nFeature_RNA < feature_upper_cut)
  counts <- log_count(ncol(sample), "QC gene num ABS")

  # Library-size and detected-gene outliers: > qc_nmads MADs from the sample median, log scale
  is_outlier <- function(x) abs(x - median(x)) > qc_nmads * mad(x)
  outlier <- is_outlier(log1p(sample$nCount_RNA)) | is_outlier(log1p(sample$nFeature_RNA))
  sample <- subset(x = sample, cells = colnames(sample)[!outlier])
  counts <- log_count(ncol(sample), "QC MAD outliers")

  # Doublets
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = LayerData(sample, layer = "counts")))
  if (!round_counts) SummarizedExperiment::assay(sce, "counts") <- round(SummarizedExperiment::assay(sce, "counts"))  # scDblFinder expects integers
  set.seed(cfg_get("qc.seed"))   # scDblFinder is stochastic
  sce <- scDblFinder(sce)
  sample$scDblFinder.class <- colData(sce)$scDblFinder.class
  sample$scDblFinder.score <- colData(sce)$scDblFinder.score
  rm(sce)
  doublets <- sample@meta.data
  doublets$SampleID <- name

  plot_doublet_scores(sample, name = name, results_path = results_path)
  plot_doublet_scatter(df = sample@meta.data, name = name, results_path = results_path)

  sample <- subset(x = sample, subset = scDblFinder.class == "singlet")
  counts <- log_count(ncol(sample), "scDblFinder")
  list(sample = sample, counts = counts, doublets = doublets)
}

# ===============================================================
# mode "sample": DecontX + QC of ONE sample, saved for the merge
# ===============================================================
if (mode == "sample") {
  name <- args[2]
  if (!name %in% samples$sample) stop("Sample '", name, "' not in samples.csv")
  s <- samples[samples$sample == name, ]
  dx <- run_decontx_sample(name, s$raw_dir, s$filtered_dir)
  counts <- data.frame(Cell_counts = ncol(dx$adjusted), Feature_discard = "Loaded cells", SampleID = name)

  tumor <- CreateSeuratObject(counts = dx$adjusted, project = name, min.cells = cfg_get("decontx.min_cells"),
                              min.features = cfg_get("decontx.min_features"))
  counts <- rbind(counts, data.frame(Cell_counts = ncol(tumor), Feature_discard = "Min features", SampleID = name))

  # DecontX contamination, fraction of counts removed (after rounding) and uncorrected QC metrics
  raw_cells <- dx$raw[, colnames(tumor)]
  rownames(raw_cells) <- gsub("_", "-", rownames(raw_cells))   # same gene names as the Seurat object
  tumor$decontX_contamination <- dx$contamination[colnames(tumor)]
  tumor$decontX_removed <- 1 - tumor$nCount_RNA / pmax(Matrix::colSums(raw_cells), 1)
  tumor$nCount_RNA_raw <- Matrix::colSums(raw_cells)
  tumor$nFeature_RNA_raw <- Matrix::colSums(raw_cells > 0)
  tumor$percent.mt_raw <- Matrix::colSums(raw_cells[grepl("^MT-", rownames(raw_cells)), , drop = FALSE]) /
    pmax(tumor$nCount_RNA_raw, 1) * 100
  contamination_summary <- dx$contamination_summary
  rm(dx); invisible(gc())

  qc <- process_sample(name = name, tumor = tumor, counts = counts, results_path = per_sample_dir)
  raw_cells <- raw_cells[, colnames(qc$sample)]   # all genes: same features in every sample

  # Post-filter plots
  plot_gene_expression_dist(tumor = qc$sample, name = name, results_path = per_sample_dir)
  plot_nfeature_ncount_corr(tumor = qc$sample, name = name, results_path = per_sample_dir)

  # Unique cell names "<barcode>_<sample index in samples.csv>"
  suffix <- paste0("_", match(name, samples$sample))
  qc$sample <- RenameCells(qc$sample, new.names = paste0(colnames(qc$sample), suffix))
  colnames(raw_cells) <- paste0(colnames(raw_cells), suffix)
  rm(tumor); invisible(gc())

  saveRDS(list(sample = qc$sample, raw = raw_cells, contamination = contamination_summary,
               cell_counts = qc$counts, doublets = qc$doublets),
          file.path(per_sample_dir, paste0(name, ".rds")), compress = FALSE)
  cat("\n Sample", name, "saved\n")
  quit(save = "no", status = 0)
}

# ===============================================================
# mode "merge": collect the per-sample results in samples.csv order (cell suffixes follow it)
# ===============================================================
ps_files <- file.path(per_sample_dir, paste0(samples$sample, ".rds"))
missing_samples <- samples$sample[!file.exists(ps_files)]
if (length(missing_samples)) stop("Per-sample results missing (run mode 'sample'): ", paste(missing_samples, collapse = ", "))

sample_list <- list(); raw_list <- list(); contamination_list <- list(); doublet_list <- list(); count_list <- list()
for (i in seq_along(ps_files)) {
  name <- samples$sample[i]
  ps <- readRDS(ps_files[i])
  sample_list[[name]] <- ps$sample; raw_list[[name]] <- ps$raw
  contamination_list[[name]] <- ps$contamination; doublet_list[[name]] <- ps$doublets
  count_list[[name]] <- ps$cell_counts
  rm(ps)
}
cat("\n Per-sample results loaded\n")

write.csv(do.call(rbind, contamination_list), file.path(decontx_dir, "decontx_contamination_bySample.csv"), row.names = FALSE)
plot_doublet_scatter(do.call(rbind, doublet_list), name = "AllSamples", results_path = out_dir)
rm(doublet_list)

cell_counts_df <- do.call(rbind, unname(count_list))
write.csv(cell_counts_df, file = file.path(out_dir, "cell_filtering.csv"))
plot_qc_funnel(cell_counts_df, results_path = out_dir)

# Merge on disk (BPCells): each sample's counts are written once as a matrix directory, padded to the
# union of genes (as merge() does); the object keeps one on-disk counts layer per sample.
bp_path <- project_path(cfg_get("bpcells_dir"), cfg_get("run_name"))
dir.create(bp_path, recursive = TRUE, showWarnings = FALSE)
all_genes <- unique(unlist(lapply(sample_list, rownames)))
bp_mats <- list(); bp_meta <- list()
for (name in names(sample_list)) {
  m <- LayerData(sample_list[[name]], assay = "RNA", layer = "counts")
  missing_genes <- setdiff(all_genes, rownames(m))
  pad <- Matrix::sparseMatrix(i = integer(0), j = integer(0), dims = c(length(missing_genes), ncol(m)),
                              dimnames = list(missing_genes, colnames(m)))
  m <- rbind(m, pad)[all_genes, ]
  unlink(file.path(bp_path, name), recursive = TRUE)
  write_matrix_dir(mat = convert_matrix_type(as(m, "CsparseMatrix"), "uint32_t"), dir = file.path(bp_path, name))
  bp_mats[[name]] <- open_matrix_dir(file.path(bp_path, name))
  bp_meta[[name]] <- sample_list[[name]]@meta.data
}
rm(sample_list, m, pad); invisible(gc())
meta <- do.call(rbind, unname(bp_meta))
merged_seurat <- CreateSeuratObject(counts = bp_mats, meta.data = meta)
merged_seurat$orig.ident <- meta[colnames(merged_seurat), "orig.ident"]   # keep the sample ids
Idents(merged_seurat) <- "orig.ident"
rm(bp_mats, bp_meta, meta); invisible(gc())
cat("\n Samples merged (BPCells:", bp_path, ")\n")

plot_vln_compare_samples(merged_seurat, results_path = out_dir)

merged_seurat <- NormalizeData(merged_seurat)
merged_seurat <- FindVariableFeatures(merged_seurat)
merged_seurat <- ScaleData(merged_seurat)
merged_seurat <- RunPCA(merged_seurat)
plot_dimplot_by_sample(merged_seurat, reduction = "pca", results_path = out_dir)

# Harmony, one batch per sample (the counts layers are split by sample)
integrated <- IntegrateLayers(object = merged_seurat, method = HarmonyIntegration,
                              orig.reduction = "pca", new.reduction = "harmony", verbose = TRUE)
rm(merged_seurat); invisible(gc())
plot_dimplot_by_sample(integrated, reduction = "harmony", results_path = out_dir)

# Clinical annotation (optional): joined by scRNAseq_ID == sample
clinical_file <- project_path(cfg_get("clinical_sheet"))
if (file.exists(clinical_file)) {
  annot_df <- read.csv(clinical_file, stringsAsFactors = FALSE) %>%
    dplyr::select(scRNAseq_ID, Subtype, IGGscale_annot, pCR, TLS) %>%
    dplyr::distinct() %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ trimws(as.character(.x))),   # e.g. "Her2+ " -> "Her2+"
                  IGG = IGGscale_annot)
  integrated@meta.data <- integrated@meta.data %>%
    tibble::rownames_to_column("cell_id") %>%
    dplyr::left_join(annot_df, by = c("orig.ident" = "scRNAseq_ID")) %>%
    tibble::column_to_rownames("cell_id")
  miss <- unique(integrated$orig.ident[is.na(integrated$Subtype)])
  if (length(miss)) cat("\n Samples not in the clinical sheet (Subtype = Unknown):", paste(miss, collapse = ", "), "\n")
  integrated$Subtype <- factor(ifelse(is.na(integrated$Subtype) | integrated$Subtype == "", "Unknown", integrated$Subtype))
} else {
  # Downstream steps group and colour by Subtype: without a sheet every sample is "Unknown"
  cat("\n Clinical sheet not found (", clinical_file, "): Subtype = Unknown for every sample; no pCR / TLS columns\n")
  annot_df <- data.frame(scRNAseq_ID = character(0), Subtype = character(0))
  integrated$Subtype <- factor("Unknown")
}
plot_pct_filtered_by_subtype(cell_counts_df, annot_df, results_path = out_dir)

integrated <- JoinLayers(integrated)
saveRDS(integrated, run_object("qc_merged"))
cat("\n Saved", run_object("qc_merged"), "\n")

# Genome-wide expression change, before vs after DecontX (cells passing QC)
raw_full <- do.call(cbind, unname(raw_list))[rownames(integrated), colnames(integrated)]
mean_before <- Matrix::rowMeans(NormalizeData(raw_full, verbose = FALSE))
mean_after <- Matrix::rowMeans(GetAssayData(integrated, assay = "RNA", layer = "data"))
shared_genes <- intersect(names(mean_before), names(mean_after))
rm(raw_full, raw_list); invisible(gc())
gene_change_df <- data.frame(gene = shared_genes, avg_before = mean_before[shared_genes], avg_after = mean_after[shared_genes]) %>%
  dplyr::mutate(avg_drop = avg_before - avg_after,
                avg_drop_pct = ifelse(avg_before > 0, (avg_before - avg_after) / avg_before * 100, NA))
write.csv(gene_change_df %>% dplyr::arrange(dplyr::desc(avg_drop)), file.path(decontx_dir, "decontx_gene_change.csv"), row.names = FALSE)

plot_decontx_scatter(gene_change_df, decontx_dir, method = "DecontX")
plot_decontx_maplot(gene_change_df, decontx_dir, method = "DecontX")
plot_decontx_hist(gene_change_df, decontx_dir, method = "DecontX")
plot_decontx_dumbbell(gene_change_df, decontx_dir, method = "DecontX")

p <- ggplot(integrated@meta.data, aes(x = orig.ident, y = decontX_contamination)) +
  geom_violin(scale = "width", fill = "grey80") + geom_boxplot(width = 0.1, outlier.shape = NA) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "DecontX: estimated contamination per cell", x = NULL, y = "decontX_contamination")
ggsave(file.path(decontx_dir, "DecontX_contamination_bySample.png"), p, width = 10, height = 5, dpi = 300, bg = "white")

cat("\n ---- 01_qc_decontx merge done. Next: 02_resolution.R ----\n")
