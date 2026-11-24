##
## Annotation dictionaries: cluster -> label, written by hand (human or agent) from the evidence
## of the markers step, never derived automatically from scores.
##
## One TSV per annotation step in results/<run>/annotation/<step>.tsv with columns
##   cluster   cluster id (as in the clustering column)
##   n_cells   cells in the cluster (pre-filled)
##   top_markers  first filtered FindMarkers genes (pre-filled, for reference)
##   label     cell type; the special label DROP removes the cluster (noise, doublets, contamination)
##   evidence  why: marker genes / scores that support the label (required)
##
## Flow:  markers step  -> write_annotation_template()  (never overwrites an edited file)
##        human/agent   -> fills label + evidence
##        annotate step -> read_annotation_dictionary() + apply_annotation()
##

annotation_file <- function(step) paste0(step_dir("annotation"), step, ".tsv")

# ---------------------------------------------------------------
# Write the empty dictionary for a step, pre-filled with cluster sizes and top markers.
#   markers_df : the filtered FindMarkers table (columns cluster, gene), optional
# ---------------------------------------------------------------
write_annotation_template <- function(object, step, cluster_col = "seurat_clusters", markers_df = NULL, n_top = 8) {
  f <- annotation_file(step)
  if (file.exists(f)) { cat("  Dictionary already exists, kept:", f, "\n"); return(invisible(f)) }
  cl <- table(as.character(object[[cluster_col, drop = TRUE]]))
  ord <- names(cl)[order(suppressWarnings(as.numeric(names(cl))), names(cl))]
  top <- vapply(ord, function(c) {
    if (is.null(markers_df)) return("")
    paste(head(markers_df$gene[as.character(markers_df$cluster) == c], n_top), collapse = ",")
  }, "")
  df <- data.frame(cluster = ord, n_cells = as.integer(cl[ord]), top_markers = top, label = "", evidence = "")
  write.table(df, f, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("  Dictionary template ->", f, "\n")
  invisible(f)
}

# ---------------------------------------------------------------
# Read and validate a filled dictionary. Stops at a CHECKPOINT when the file is missing,
# a cluster of the object is not listed, or a label/evidence is empty.
# ---------------------------------------------------------------
read_annotation_dictionary <- function(step, clusters) {
  f <- annotation_file(step)
  if (!file.exists(f)) checkpoint_stop(sprintf("no dictionary for '%s'. Run the markers step first; it writes %s", step, f))
  d <- read.delim(f, colClasses = "character", na.strings = character(0))
  d$label <- trimws(d$label); d$evidence <- trimws(d$evidence)
  clusters <- unique(as.character(clusters))
  problems <- c(
    if (length(miss <- setdiff(clusters, d$cluster))) paste("clusters not in the dictionary:", paste(miss, collapse = ", ")),
    if (length(extra <- setdiff(d$cluster, clusters))) paste("dictionary clusters absent from the object:", paste(extra, collapse = ", ")),
    if (any(dup <- duplicated(d$cluster))) paste("duplicated clusters:", paste(d$cluster[dup], collapse = ", ")),
    if (any(d$label == "")) paste("empty label for clusters:", paste(d$cluster[d$label == ""], collapse = ", ")),
    if (any(d$evidence == "")) paste("empty evidence for clusters:", paste(d$cluster[d$evidence == ""], collapse = ", "))
  )
  if (length(problems)) checkpoint_stop(paste0("dictionary ", f, " is not complete:\n  - ", paste(problems, collapse = "\n  - ")))
  setNames(d$label, d$cluster)
}

# ---------------------------------------------------------------
# Apply a dictionary: sets label_col from cluster_col and removes DROP clusters.
# Returns list(object, dropped = cell names removed), so the step can log them.
# ---------------------------------------------------------------
apply_annotation <- function(object, dict, cluster_col = "seurat_clusters", label_col = "celltype") {
  lab <- unname(dict[as.character(object[[cluster_col, drop = TRUE]])])
  object[[label_col]] <- lab
  dropped <- colnames(object)[lab == "DROP"]
  if (length(dropped)) {
    cat(sprintf("  DROP: %d cells (%s)\n", length(dropped), paste(names(dict)[dict == "DROP"], collapse = ", ")))
    object <- subset(object, cells = setdiff(colnames(object), dropped))
  }
  print(table(object[[label_col, drop = TRUE]]))
  list(object = object, dropped = dropped)
}

# ---------------------------------------------------------------
# Transfer labels of a re-clustered subset back to the parent object by cell name.
# Cells of the parent not in the subset keep their label; cells dropped in the
# subset are removed from the parent too.
# ---------------------------------------------------------------
transfer_labels <- function(parent, subset_obj, dropped = character(0), label_col = "celltype") {
  parent <- subset(parent, cells = setdiff(colnames(parent), dropped))
  idx <- match(colnames(subset_obj), colnames(parent))
  if (anyNA(idx)) stop("transfer_labels: subset has cells absent from the parent")
  lab <- as.character(parent[[label_col, drop = TRUE]])
  lab[idx] <- as.character(subset_obj[[label_col, drop = TRUE]])
  parent[[label_col]] <- lab
  parent
}

# Append dropped cells to results/<run>/annotation/excluded_cells.csv (step, cell, sample, label reason)
log_excluded_cells <- function(object_before, dropped, step, cluster_col = "seurat_clusters") {
  if (!length(dropped)) return(invisible(NULL))
  f <- paste0(step_dir("annotation"), "excluded_cells.csv")
  df <- data.frame(step = step, cell = dropped,
                   sample = as.character(object_before$orig.ident[dropped]),
                   cluster = as.character(object_before[[cluster_col, drop = TRUE]][match(dropped, colnames(object_before))]))
  write.table(df, f, sep = ",", row.names = FALSE, col.names = !file.exists(f), append = file.exists(f))
  invisible(df)
}
