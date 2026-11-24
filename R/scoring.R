##
## Score-based summaries used as annotation EVIDENCE (never as automatic labels):
## scSubtype scores, per-cluster score tests, composition tests, proportions.
##

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
#                (load_scsubtype_signatures())
#   assay      : assay providing gene expression (default "RNA")
#
# ---------------------------------------------------------------
compute_scsubtype_scores <- function(object, signatures, assay = "RNA") {
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
  object[[assay]]$scale.data <- NULL   # replaced by the signature genes only (previous layer held other genes)
  object <- ScaleData(object, assay = assay, features = temp_allgenes, verbose = FALSE)
  tocalc <- as.data.frame(as.matrix(GetAssayData(object, assay = assay, layer = "scale.data")))   # as.matrix: scale.data may be a BPCells matrix

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
  md <- if (inherits(seurat_obj, "Seurat")) seurat_obj@meta.data else seurat_obj  # object or its meta.data

  # celltypes absent from a sample count as 0 (not dropped)
  counts <- md %>%
    dplyr::count(.data[[sample_id_col]], .data[[celltype_col]], name = "n_cells") %>%
    tidyr::complete(!!rlang::sym(sample_id_col), !!rlang::sym(celltype_col), fill = list(n_cells = 0L))

  totals <- md %>%
    dplyr::count(.data[[sample_id_col]], name = "n_total")

  clinical <- md %>%
    dplyr::distinct(dplyr::across(dplyr::all_of(c(sample_id_col, clinical_cols))))

  counts %>%
    dplyr::left_join(totals, by = sample_id_col) %>%
    dplyr::mutate(proportion = n_cells / n_total) %>%
    dplyr::left_join(clinical, by = sample_id_col)
}
