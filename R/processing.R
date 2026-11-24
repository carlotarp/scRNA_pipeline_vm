##
## Dimensionality reduction, integration and (sub-)clustering.
##

# ---------------------------------------------------------------
# Number of PCs. integration.n_pcs in config.yaml (30 in the reference run) wins;
# with n_pcs: null the HBC heuristic is used instead, on each PC's share of the summed
# standard deviations:
#   co1 : first PC with cumulative share > 90 % and own share < 5 %
#   co2 : last PC whose share drops by > 0.1 points to the next one (elbow)
#   n   : min(co1, co2)
# ---------------------------------------------------------------
select_n_pcs <- function(object, reduction = "pca", fixed = cfg_get("integration.n_pcs")) {
  if (!is.null(fixed)) return(as.integer(fixed))
  stdev <- object[[reduction]]@stdev
  pct   <- stdev / sum(stdev) * 100
  co1   <- which(cumsum(pct) > 90 & pct < 5)[1]
  co2   <- sort(which((pct[-length(pct)] - pct[-1]) > 0.1), decreasing = TRUE)[1] + 1
  n_pcs <- min(co1, co2, na.rm = TRUE)
  if (!is.finite(n_pcs)) stop("select_n_pcs: no cutoff found, inspect the elbow plot")
  n_pcs
}

# ---------------------------------------------------------------
# UMAP with PCA initialisation (uwot init = "pca"). Seurat's default spectral init crashed at
# ~1M cells; PCA init was at least as faithful to the Harmony space. Only the 2-D picture depends
# on it: clusters and markers use the Harmony graph. Other settings = RunUMAP defaults.
# ---------------------------------------------------------------
run_umap_pca <- function(object, reduction, dims, reduction.name = "umap", n_threads = 4) {
  set.seed(42)
  emb <- uwot::umap(Embeddings(object, reduction)[, dims, drop = FALSE], n_neighbors = 30, metric = "cosine",
                    min_dist = 0.3, init = "pca", init_sdev = "range", n_threads = n_threads)
  key <- paste0(gsub("[^A-Za-z0-9]", "", reduction.name), "_")
  colnames(emb) <- paste0(key, 1:2)
  object[[reduction.name]] <- CreateDimReducObject(emb, key = key, assay = DefaultAssay(object))
  object
}

# ---------------------------------------------------------------
# Subset to some cells, re-normalise, re-integrate with Harmony (one batch per
# integration.batch_key level) and re-cluster.
#
# Args:
#   cells      : cell names to keep
#   resolution : FindClusters resolution; a vector clusters at every value
#                (columns RNA_snn_res.<r>), the last one becomes seurat_clusters
#   n_pcs      : PCs for neighbours/UMAP (default: select_n_pcs())
# ---------------------------------------------------------------
generate_subset <- function(object, cells, resolution, assay = "RNA", label = "subset", n_pcs = NULL,
                            split_by = cfg_get("integration.batch_key")) {
  sub <- subset(object, cells = cells)
  sub[[assay]]$scale.data <- NULL   # parent's scaled genes; rebuilt below for this subset
  sub[[assay]] <- split(sub[[assay]], f = sub[[split_by]][, 1])
  sub <- NormalizeData(sub)
  sub <- FindVariableFeatures(sub)
  sub <- ScaleData(sub)
  sub <- RunPCA(sub)
  sub <- IntegrateLayers(sub, method = HarmonyIntegration, orig.reduction = "pca",
                         new.reduction = "harmony", verbose = FALSE)
  n_pcs <- n_pcs %||% select_n_pcs(sub)
  cat(sprintf("  PCs (%s): %d\n", label, n_pcs))
  sub <- FindNeighbors(sub, reduction = "harmony", dims = 1:n_pcs)
  sub <- run_umap_pca(sub, reduction = "harmony", dims = 1:n_pcs)
  sub <- FindClusters(sub, resolution = resolution)
  JoinLayers(sub)
}

# Subset to one value of $lineage, then re-cluster (wraps generate_subset)
generate_lineage_subset <- function(object, lineage_name, resolution, assay = "RNA", n_pcs = NULL) {
  keep <- colnames(object)[object$lineage == lineage_name]
  if (!length(keep)) stop("No cells with lineage == '", lineage_name, "'")
  generate_subset(object, cells = keep, resolution = resolution, assay = assay, label = lineage_name, n_pcs = n_pcs)
}
