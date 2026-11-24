##
## Marker panels (resources/*.yaml) and per-cluster differential expression.
##

# ---------------------------------------------------------------
# Load a marker panel from resources/marker_panels.yaml (or another yaml in resources/).
#   load_marker_panel("lineage")                 -> list of blocks (leukocytes, stromal, ...)
#   load_marker_panel("lineage", "leukocytes")   -> named list group -> genes
#   object : when given, genes absent from it are dropped and reported
# ---------------------------------------------------------------
load_marker_panel <- function(panel, block = NULL, object = NULL, file = "marker_panels.yaml") {
  all <- yaml::read_yaml(file.path(PIPELINE_DIR, "resources", file))
  if (!panel %in% names(all)) stop("Panel '", panel, "' not in ", file, " (have: ", paste(names(all), collapse = ", "), ")")
  p <- all[[panel]]
  if (!is.null(block)) p <- p[[block]]
  if (is.null(object)) return(p)
  genes_in <- rownames(object)
  drop_missing <- function(x) {
    if (is.list(x)) return(lapply(x, drop_missing))
    miss <- setdiff(x, genes_in)
    if (length(miss)) cat("  [panel ", panel, "] not in matrix: ", paste(miss, collapse = ", "), "\n", sep = "")
    intersect(x, genes_in)
  }
  drop_missing(p)
}

# scSubtype signatures (Wu et al. 2021) as a named list of gene vectors
load_scsubtype_signatures <- function() yaml::read_yaml(file.path(PIPELINE_DIR, "resources", "scsubtype_wu2021.yaml"))

# ---------------------------------------------------------------
# FindMarkers for each cluster vs the rest. Writes, per cluster and combined:
#   cluster_<cl>_FindMarkers.csv           top_n by avg_log2FC, unfiltered
#   cluster_<cl>_FindMarkers_filtered.csv  same ranking, only padj < markers.padj and pct.1 >= markers.min_pct1
#   all_clusters_FindMarkers[_filtered].csv
# The filtered list is the one to read when annotating: unfiltered top genes are often
# expressed in a handful of cells (ambient RNA, doublets).
#
# Returns (invisibly) the combined filtered data.frame.
# ---------------------------------------------------------------
find_markers_for_clusters <- function(object, clusters, results_path, top_n = 20,
                                      padj = cfg_get("markers.padj"), min_pct1 = cfg_get("markers.min_pct1")) {
  clusters <- as.character(clusters)
  if (length(clusters) == 0) { cat("  - No clusters to test, FindMarkers skipped\n"); return(invisible(NULL)) }
  object <- JoinLayers(object)

  all_raw <- list(); all_filt <- list()
  for (cl in clusters) {
    m <- FindMarkers(object, ident.1 = cl, max.cells.per.ident = 5000)
    m <- m[order(-m$avg_log2FC), ]
    m$gene <- rownames(m); m$cluster <- cl
    raw  <- head(m, top_n)
    filt <- head(m[m$p_val_adj < padj & m$pct.1 >= min_pct1, ], top_n)
    write.csv(raw,  file.path(results_path, paste0("cluster_", cl, "_FindMarkers.csv")), row.names = FALSE)
    write.csv(filt, file.path(results_path, paste0("cluster_", cl, "_FindMarkers_filtered.csv")), row.names = FALSE)
    all_raw[[cl]] <- raw; all_filt[[cl]] <- filt
    cat("  - Cluster", cl, ":", nrow(filt), "filtered markers\n")
  }
  order_cols <- function(df) { rownames(df) <- NULL; df[, c("cluster", "gene", setdiff(colnames(df), c("cluster", "gene")))] }
  combined_raw  <- order_cols(do.call(rbind, all_raw))
  combined_filt <- order_cols(do.call(rbind, all_filt))
  write.csv(combined_raw,  file.path(results_path, "all_clusters_FindMarkers.csv"), row.names = FALSE)
  write.csv(combined_filt, file.path(results_path, "all_clusters_FindMarkers_filtered.csv"), row.names = FALSE)
  invisible(combined_filt)
}

# ---------------------------------------------------------------
# Top expressed genes per cluster (mean normalised expression), excluding mitochondrial,
# ribosomal and MALAT1/NEAT1. Complements FindMarkers: a cluster whose top expressed genes
# are another lineage's markers is a doublet/ambient candidate even if its DE list looks clean.
# Writes top_expressed_genes.csv (cluster, rank, gene, mean_expr, pct_expr).
# ---------------------------------------------------------------
top_expressed_genes <- function(object, group_by = "seurat_clusters", results_path, top_n = 15) {
  object <- JoinLayers(object)
  genes <- rownames(object)
  genes <- genes[!grepl("^MT-|^RP[SL]\\d|^MALAT1$|^NEAT1$", genes)]
  grp <- as.character(object[[group_by, drop = TRUE]])
  data <- LayerData(object, assay = "RNA", layer = "data")[genes, , drop = FALSE]
  out <- do.call(rbind, lapply(sort(unique(grp)), function(g) {
    sub <- data[, grp == g, drop = FALSE]
    mu <- Matrix::rowMeans(sub); pct <- Matrix::rowMeans(sub > 0)
    o <- head(order(-mu), top_n)
    data.frame(cluster = g, rank = seq_along(o), gene = genes[o], mean_expr = round(mu[o], 3), pct_expr = round(pct[o], 3))
  }))
  write.csv(out, file.path(results_path, "top_expressed_genes.csv"), row.names = FALSE)
  invisible(out)
}
