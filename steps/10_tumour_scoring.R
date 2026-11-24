#!/usr/bin/env Rscript
##
## 10_tumour_scoring.R - re-cluster the tumour cells and score PAM50-like modules, PROGENy, cell cycle, scSubtype
## Usage : Rscript pipeline/steps/10_tumour_scoring.R
## Reads : objects/notumor_annotated.rds
## Writes: objects/tumour_grid.rds (cache), objects/tumour_scored.rds ;
##         10_tumour_scoring/resolution_grid/<grid evidence>, 10_tumour_scoring/<UMAP, phase/PAM50 tables>
## Checkpoint: subclustering.tumour.resolution (phase 1 writes the grid evidence, then stops until decided)
##
file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
PIPELINE_DIR <- normalizePath(file.path(dirname(file_arg), ".."))
source(file.path(PIPELINE_DIR, "R/setup.R"))
out_dir <- step_dir("10_tumour_scoring")

# =============================================================================
# Phase 1: tumour subset clustered at every grid resolution (cached) + evidence
# =============================================================================
grid <- cfg_get("clustering.resolution_grid")
grid_file <- run_object("tumour_grid")
if (file.exists(grid_file)) {
  tumour <- read_run_object("tumour_grid")
} else {
  annotated <- read_run_object("notumor_annotated")
  tumour <- generate_lineage_subset(annotated, lineage_name = "Tumor", resolution = grid,
                                    n_pcs = cfg_get("subclustering.tumour.n_pcs"))
  rm(annotated); invisible(gc())
  saveRDS(tumour, grid_file)
  cat("  Tumour grid cache ->", grid_file, "\n")
}

grid_dir <- step_dir("10_tumour_scoring", "resolution_grid")
sizes_csv <- file.path(grid_dir, "cluster_sizes_per_resolution.csv")
if (!file.exists(sizes_csv)) {
  plot_resolution_grid(tumour, grid_dir, resolutions = grid)
  sizes <- do.call(rbind, lapply(grid, function(r) {
    tb <- table(tumour[[paste0("RNA_snn_res.", r)]][, 1])
    data.frame(resolution = r, cluster = names(tb), n_cells = as.integer(tb))
  }))
  write.csv(sizes, sizes_csv, row.names = FALSE)
}

# =============================================================================
# Phase 2: chosen resolution
# =============================================================================
res <- get_decision("subclustering.tumour.resolution")
res_col <- paste0("RNA_snn_res.", res)
if (!res_col %in% colnames(tumour@meta.data)) tumour <- FindClusters(tumour, resolution = res)
tumour$seurat_clusters <- tumour[[res_col]][, 1]
Idents(tumour) <- "seurat_clusters"
cat(sprintf("  Tumour resolution %s: %d clusters\n", res, nlevels(tumour$seurat_clusters)))

plot_dimplot(tumour, reduction = "umap", group_by = "seurat_clusters",
             results_path = out_dir, filename = "DimPlot_UMAP_Tumor.png", label = TRUE)

# --- PAM50-like module scores; predicted subtype = highest-scoring module ---
pam50_genes <- load_marker_panel("tumour_pam50", object = tumour)
tumour <- AddModuleScore(tumour, features = pam50_genes, name = "PAM50_")
pam50_score_cols <- paste0("PAM50_", seq_along(pam50_genes))
pam50_named_cols <- paste0("PAM50_", names(pam50_genes))
colnames(tumour@meta.data)[match(pam50_score_cols, colnames(tumour@meta.data))] <- pam50_named_cols

pam50_scores <- as.matrix(tumour@meta.data[, pam50_named_cols])
tumour$PAM50_predicted <- factor(names(pam50_genes)[apply(pam50_scores, 1, which.max)], levels = names(pam50_genes))
print(table(tumour$PAM50_predicted))
write.csv(as.data.frame.matrix(table(Predicted = tumour$PAM50_predicted, Clinical = tumour$Subtype)),
          file.path(out_dir, "PAM50_vs_Clinical_concordance.csv"))

# --- PROGENy pathway activity, stored as assay 'progeny' ---
# Scores are per cell, so they are computed in chunks of 5,000 cells (a dense genes x all-cells matrix
# does not fit in memory on large data) and scaled across all cells at the end (= progeny(scale = TRUE)).
expr <- LayerData(tumour, assay = "RNA", layer = "data")
chunks <- split(colnames(tumour), ceiling(seq_len(ncol(tumour)) / 5000))
progeny_scores <- do.call(rbind, lapply(chunks, function(cells)
  progeny::progeny(as.matrix(expr[, cells]), scale = FALSE, organism = "Human", top = 500, perm = 1)))
progeny_scores <- scale(progeny_scores[colnames(tumour), , drop = FALSE])
rm(expr); invisible(gc())
# Scores have no counts layer: Seurat's nCount/nFeature lookup for the new assay warns "No layers found"; that
# one warning is expected and muffled here
withCallingHandlers(tumour[["progeny"]] <- CreateAssay5Object(data = as(t(progeny_scores), "dgCMatrix")),
                    warning = function(w) if (grepl("No layers found", conditionMessage(w))) invokeRestart("muffleWarning"))
tumour <- ScaleData(tumour, assay = "progeny", verbose = FALSE)
cat("  PROGENy pathway activities added as assay 'progeny'\n")

# --- Cell cycle phase ---
cc <- lapply(cc.genes.updated.2019[c("s.genes", "g2m.genes")], function(g) {
  miss <- setdiff(g, rownames(tumour)); if (length(miss)) cat("  cell-cycle genes not in matrix:", paste(miss, collapse = ", "), "\n")
  intersect(g, rownames(tumour))
})
tumour <- CellCycleScoring(tumour, s.features = cc$s.genes, g2m.features = cc$g2m.genes)
print(table(tumour$Phase))
write.csv(as.data.frame.matrix(table(PAM50 = tumour$PAM50_predicted, Phase = tumour$Phase)),
          file.path(out_dir, "CellCyclePhase_by_PAM50.csv"))
write.csv(as.data.frame.matrix(table(Cluster = tumour$global_clusters, Phase = tumour$Phase)),
          file.path(out_dir, "CellCyclePhase_by_cluster.csv"))

# --- scSubtype (Wu et al. 2021): scores + scSubtype_call ---
tumour <- compute_scsubtype_scores(tumour, load_scsubtype_signatures(), assay = "RNA")
print(table(tumour$scSubtype_call))

saveRDS(tumour, run_object("tumour_scored"))
cat("\n ---- FINISHED 10_tumour_scoring: objects/tumour_scored.rds (PAM50, progeny assay, Phase, scSubtype). Next: 11_tumour_markers.R ----\n")
