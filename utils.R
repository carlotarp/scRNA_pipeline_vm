##
## Shared utility functions for the scRNA-seq annotation pipeline.
##

library(Seurat)
library(dplyr)

# ---------------------------------------------------------------
# Fixed color palettes for subtype categories
# ---------------------------------------------------------------
pam50_subtype_colors <- c("Lum" = "#00BA38", "TNBC" = "#619CFF",
                          "Her2+" = "#F8766D", "Her2+ " = "#F8766D")

scsubtype_colors <- c("Basal_SC" = "#619CFF", "Her2E_SC" = "#F8766D",
                      "LumA_SC" = "lightgreen", "LumB_SC" = "darkgreen")

# Clinical IGG level palette
igg_level_colors <- c("High" = "#A11212", "Low" = "#F2E205", "Medium" = "#F5A623")

# pCR and TLS clinical covariate palettes
pcr_colors <- c("No" = "white", "Yes" = "black")
tls_colors <- c("0" = "#D9D9D9", "1a" = "#9013FE", "1b" = "#C77DD6", "2a" = "#F9B8D0", "2b" = "#FF00FF")

# Broad lineage palette
lineage_colors <- c("Tumor" = "#875cb8", "Leukocytes" = "#d6790e", "Stromal" = "#0a6807")


# ---------------------------------------------------------------
# Read a single sample's copyKAT output files back from disk.
#
# Args:
#   sample_path : directory containing the two output files (the per-sample
#                 working directory copykat() was run from)
#   sam_name    : the "sam.name" copykat() was run with (the file prefix)
#
# Returns list(
#   prediction = data.frame(cell.names, copykat.pred),
#   cna        = data.frame(chrom, chrompos, abspos, <one column per cell>)
# )
# ---------------------------------------------------------------
read_copykat_outputs <- function(sample_path, sam_name) {
  pred_file <- file.path(sample_path, paste0(sam_name, "_copykat_prediction.txt"))
  cna_file  <- file.path(sample_path, paste0(sam_name, "_copykat_final_results_bin_by_cell.txt"))

  prediction <- read.table(pred_file, header = TRUE, sep = "\t",
                           stringsAsFactors = FALSE, check.names = FALSE)
  cna <- read.table(cna_file, header = TRUE, sep = "\t",
                    stringsAsFactors = FALSE, check.names = FALSE)

  bin_cols <- setdiff(colnames(cna), c("chrom", "chrompos", "abspos"))
  colnames(cna)[colnames(cna) %in% bin_cols] <- sub("\\.", "-", bin_cols)

  list(prediction = prediction, cna = cna)
}

# ---------------------------------------------------------------
# Subset object to a set of cells, normalize, recluster with Harmony.
#
# Args:
#   object     : Seurat object
#   cells      : character vector of cell names to keep
#   resolution : FindClusters resolution
#   assay      : assay to split/normalize (default "RNA_decontX")
#   label      : optional tag used only in the "Optimal PCs" message
# ---------------------------------------------------------------
generate_subset <- function(object, cells, resolution, assay = "RNA_decontX", label = "subset") {
  object_subset <- subset(object, cells = cells)

  object_subset[[assay]] <- split(object_subset[[assay]], f = object_subset$orig.ident)
  object_subset <- NormalizeData(object_subset)
  object_subset <- FindVariableFeatures(object_subset)
  object_subset <- ScaleData(object_subset)
  object_subset <- RunPCA(object_subset)

  object_subset <- IntegrateLayers(
    object = object_subset, method = HarmonyIntegration,
    orig.reduction = "pca", new.reduction = "harmony", verbose = FALSE
  )

  pct   <- object_subset[["pca"]]@stdev / sum(object_subset[["pca"]]@stdev) * 100
  cumu  <- cumsum(pct)
  co1   <- which(cumu > 90 & pct < 5)[1]
  co2   <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = TRUE)[1] + 1
  co3   <- min(co1, co2)
  cat(paste0("  Optimal PCs (", label, "): ", co3, "\n"))

  object_subset <- FindNeighbors(object_subset, reduction = "harmony", dims = 1:co3)
  object_subset <- RunUMAP(object_subset, reduction = "harmony", dims = 1:co3)
  object_subset <- FindClusters(object_subset, resolution = resolution)
  object_subset <- JoinLayers(object_subset)

  return(object_subset)
}

# ---------------------------------------------------------------
# Subset object to one lineage, then recluster (wraps generate_subset()).
#
# Args:
#   object       : Seurat object with a $lineage metadata column
#   lineage_name : value to filter on (e.g. "Leukocytes", "Stromal", "Tumor")
#   resolution   : FindClusters resolution
#   assay        : assay to split/normalize (default "RNA_decontX")
# ---------------------------------------------------------------
generate_lineage_subset <- function(object, lineage_name, resolution, assay = "RNA_decontX") {
  keep_cells <- colnames(object)[object$lineage == lineage_name]
  generate_subset(object, cells = keep_cells, resolution = resolution,
                  assay = assay, label = lineage_name)
}

# ---------------------------------------------------------------
# Run FindMarkers for a set of clusters, save one CSV per cluster plus a combined CSV.
#
# Args:
#   object       : Seurat object
#   clusters     : character/numeric vector of cluster identities to test
#   results_path : directory where CSVs are written
#   top_n        : number of top markers to keep (ranked by avg_log2FC)
#
# Returns (invisibly) the combined data.frame; also written as
#   all_clusters_FindMarkers.csv
# ---------------------------------------------------------------
find_markers_for_clusters <- function(object, clusters, results_path, top_n = 20) {
  clusters <- as.character(clusters)
  object_joined <- JoinLayers(object)

  all_markers <- list()
  for (cl in clusters) {
    markers_out <- FindMarkers(object_joined, ident.1 = cl, max.cells.per.ident = 5000)
    markers_out <- head(markers_out[order(-markers_out$avg_log2FC), ], top_n)
    write.csv(markers_out, file.path(results_path, paste0("cluster_", cl, "_FindMarkers.csv")))
    cat(paste0("  - Cluster ", cl, " -> cluster_", cl, "_FindMarkers.csv\n"))

    markers_out$gene <- rownames(markers_out)
    markers_out$cluster <- cl
    all_markers[[cl]] <- markers_out
  }

  combined <- do.call(rbind, all_markers)
  rownames(combined) <- NULL
  combined <- combined[, c("cluster", "gene",
                           setdiff(colnames(combined), c("cluster", "gene")))]
  write.csv(combined, file.path(results_path, "all_clusters_FindMarkers.csv"), row.names = FALSE)
  cat(paste0("  - Combined -> all_clusters_FindMarkers.csv (", nrow(combined), " rows)\n"))

  invisible(combined)
}

# ---------------------------------------------------------------
# Wilcoxon test of each score column, cluster vs the rest
#
# Args:
#   object     : Seurat object
#   score_cols : numeric metadata columns to test
#   group_by   : clustering column to test each score against (default
#                "seurat_clusters")
# ---------------------------------------------------------------
differential_scores_by_cluster <- function(object, score_cols, group_by = "seurat_clusters") {
  meta <- object@meta.data
  if (!group_by %in% colnames(meta)) {
    stop(sprintf("differential_scores_by_cluster: '%s' is not a metadata column", group_by))
  }
  clusters <- sort(unique(as.character(meta[[group_by]])))

  results <- do.call(rbind, lapply(clusters, function(cl) {
    do.call(rbind, lapply(score_cols, function(sc) {
      in_group  <- meta[[sc]][meta[[group_by]] == cl]
      out_group <- meta[[sc]][meta[[group_by]] != cl]
      wtest <- wilcox.test(in_group, out_group)
      data.frame(
        cluster  = cl,
        score    = sc,
        mean_in  = mean(in_group),
        mean_out = mean(out_group),
        avg_diff = mean(in_group) - mean(out_group),
        p_val    = wtest$p.value
      )
    }))
  }))

  results$p_val_adj <- p.adjust(results$p_val, method = "BH")
  results[order(results$cluster, -abs(results$avg_diff)), ]
}

# ---------------------------------------------------------------
# Fisher's exact test of each category's representation, cluster vs the rest.
#
# Args:
#   meta      : a @meta.data data.frame
#   group_col : clustering column to test (default "orig.ident")
#   comp_col  : categorical column whose composition is tested
#               (default "celltype")
# ---------------------------------------------------------------
differential_composition <- function(meta, group_col = "orig.ident", comp_col = "celltype") {
  clusters <- sort(unique(as.character(meta[[group_col]])))
  celltypes <- sort(unique(as.character(meta[[comp_col]])))

  results <- do.call(rbind, lapply(clusters, function(cl) {
    do.call(rbind, lapply(celltypes, function(ct) {

      # --- 2x2 contingency table: this celltype vs others, in this cluster vs rest ---
      in_cluster <- meta[[group_col]] == cl
      is_celltype <- meta[[comp_col]] == ct

      a <- sum(in_cluster & is_celltype)      # this celltype, this cluster
      b <- sum(in_cluster & !is_celltype)     # other celltypes, this cluster
      c <- sum(!in_cluster & is_celltype)     # this celltype, other clusters
      d <- sum(!in_cluster & !is_celltype)    # other celltypes, other clusters

      pct_in <- a / (a + b) * 100
      pct_out <- c / (c + d) * 100

      ftest <- fisher.test(matrix(c(a, b, c, d), nrow = 2))

      data.frame(cluster = cl, celltype = ct,
                 n_in = a, pct_in = pct_in, pct_out = pct_out,
                 pct_diff = pct_in - pct_out,
                 odds_ratio = unname(ftest$estimate),
                 p_val = ftest$p.value)
    }))
  }))

  results$p_val_adj <- p.adjust(results$p_val, method = "BH")
  results[order(results$cluster, -abs(results$pct_diff)), ]
}

# ---------------------------------------------------------------
# scSubtype scoring (Wu et al. 2021, Nat Genet)
#
# Args:
#   object     : Seurat object (tumor cells)
#   signatures : named list of character vectors, one gene signature each
#                (e.g. scsubtype_signatures from scsubtype_signatures.R)
#   assay      : assay providing gene expression (default "RNA_decontX")
#
# ---------------------------------------------------------------
compute_scsubtype_scores <- function(object, signatures, assay = "RNA_decontX") {
  # Keep only genes present in the object
  check_signature <- function(genes, sig_name) {
    available <- rownames(object[[assay]])
    hit <- genes %in% available
    missing <- genes[!hit]

    cat(sprintf("[%s] %d/%d genes found (%d missing)\n",
                sig_name, sum(hit), length(genes), length(missing)))
    if (length(missing) > 0) {
      cat(sprintf("  Missing from %s: %s\n", sig_name, paste(missing, collapse = ", ")))
    }
    genes[hit]
  }

  signatures_checked <- lapply(names(signatures), function(nm) check_signature(signatures[[nm]], nm))
  names(signatures_checked) <- names(signatures)
  new_names <- names(signatures_checked)
  temp_allgenes <- unique(unlist(signatures_checked))

  # Per-signature mean of scaled expression
  object <- ScaleData(object, assay = assay, features = temp_allgenes)
  tocalc <- as.data.frame(GetAssayData(object, assay = assay, layer = "scale.data"))

  outdat <- matrix(0, nrow = length(signatures_checked), ncol = ncol(tocalc),
                    dimnames = list(new_names, colnames(tocalc)))
  for (nm in new_names) {
    genes <- which(rownames(tocalc) %in% unique(signatures_checked[[nm]]))
    outdat[nm, ] <- apply(tocalc[genes, , drop = FALSE], 2, function(x) mean(as.numeric(x), na.rm = TRUE))
  }

  final <- outdat[rowSums(outdat, na.rm = TRUE) != 0, , drop = FALSE]
  final <- as.data.frame(final)
  is.num <- sapply(final, is.numeric); final[is.num] <- lapply(final[is.num], round, 4)
  finalm <- as.matrix(final)

  # Center scores per cell and call the subtype with the highest centered score
  center_sweep <- function(x, row.w = rep(1, nrow(x)) / nrow(x)) {
    average <- apply(x, 2, function(v) sum(v * row.w) / sum(row.w))
    sweep(x, 2, average)
  }

  finalmt <- as.data.frame(t(finalm))
  finalm.sweep.t <- center_sweep(finalmt)
  calls <- colnames(finalm.sweep.t)[max.col(finalm.sweep.t, ties.method = "first")]

  object$scSubtype_call <- calls
  object@meta.data[, paste0(new_names, "_centered")] <- finalm.sweep.t[colnames(object), new_names]
  object@meta.data[, new_names] <- t(finalm)[colnames(object), new_names]

  object
}


# ---------------------------------------------------------------
# Per-sample proportion of each group level, correlated (Spearman) against
# a numeric-encoded ordinal clinical covariate.
#
# Args:
#   meta            : a @meta.data data.frame with orig.ident, group_col, clinical_col
#   group_col       : grouping variable (e.g. "celltype")
#   clinical_col    : ordinal clinical covariate
#   clinical_levels : named numeric vector mapping clinical_col's values to an
#                     ordinal score (e.g. c(No = 0, Yes = 1))
# ---------------------------------------------------------------
proportion_correlation <- function(meta, group_col, clinical_col, clinical_levels) {
  prop_df <- meta %>%
    dplyr::filter(!is.na(.data[[group_col]]), !is.na(.data[[clinical_col]])) %>%
    dplyr::count(orig.ident, .data[[group_col]], name = "n") %>%
    dplyr::group_by(orig.ident) %>%
    dplyr::mutate(pct = n / sum(n) * 100) %>%
    dplyr::ungroup()
  colnames(prop_df)[colnames(prop_df) == group_col] <- ".group"

  sample_clinical <- meta %>%
    dplyr::distinct(orig.ident, .data[[clinical_col]]) %>%
    dplyr::mutate(.score = unname(clinical_levels[as.character(.data[[clinical_col]])]))

  prop_df <- dplyr::left_join(prop_df, sample_clinical, by = "orig.ident")

  groups <- sort(unique(as.character(prop_df$.group)))
  do.call(rbind, lapply(groups, function(g) {
    sub <- prop_df[as.character(prop_df$.group) == g & !is.na(prop_df$.score), ]
    if (length(unique(sub$.score)) < 2 || nrow(sub) < 3) {
      return(data.frame(group = g, rho = NA_real_, p_val = NA_real_, n_samples = nrow(sub)))
    }
    ctest <- suppressWarnings(cor.test(sub$pct, sub$.score, method = "spearman"))
    data.frame(group = g, rho = unname(ctest$estimate), p_val = ctest$p.value, n_samples = nrow(sub))
  }))
}

# ---------------------------------------------------------------
# Per-sample proportion of each celltype, with clinical covariates attached.
#
# Args:
#   seurat_obj    : Seurat object
#   sample_id_col : column identifying each sample (default "orig.ident")
#   celltype_col  : column identifying the cell type/cluster (default "celltype")
#   clinical_cols : clinical covariate columns to attach (assumed constant
#                   within each sample)
# ---------------------------------------------------------------
compute_celltype_proportions <- function(seurat_obj, sample_id_col = "orig.ident",
                                         celltype_col = "celltype", clinical_cols) {
  md <- seurat_obj@meta.data

  counts <- md %>%
    dplyr::count(.data[[sample_id_col]], .data[[celltype_col]], name = "n_cells")

  totals <- md %>%
    dplyr::count(.data[[sample_id_col]], name = "n_total")

  clinical <- md %>%
    dplyr::distinct(dplyr::across(dplyr::all_of(c(sample_id_col, clinical_cols))))

  counts %>%
    dplyr::left_join(totals, by = sample_id_col) %>%
    dplyr::mutate(proportion = n_cells / n_total) %>%
    dplyr::left_join(clinical, by = sample_id_col)
}

