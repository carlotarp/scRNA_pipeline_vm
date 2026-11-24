#!/usr/bin/env Rscript
##
## 17_subpopulation_evidence.R - generic sub-population study: zoom into some cell types and collect the
## evidence to re-annotate their hidden sub-groups by hand (real subtypes vs doublets vs ambient RNA).
## Each study is an entry of subpopulations.<study> in config.yaml (e.g. "caf", the fibroblasts).
##
## Usage : Rscript pipeline/steps/17_subpopulation_evidence.R <study> <mode> [sample]
##   check        : lineage genes of the focus label next to the compare types and the real owners of those genes;
##                  same-cell co-expression of own + foreign genes
##   profile      : focus label per sample, quality, markers vs the other zoomed labels, foreign genes per cell
##                  (doublet = many foreign genes in one cell, ambient = 1-2)
##   zoom         : all <celltypes> together, joined counts, Harmony per sample; two phases (grid cache, then
##                  subpopulations.<study>.resolution) -> markers + annotation/<study>.tsv template
##   signatures   : published sub-type signatures (panels.signatures) + gate genes per zoom cluster
##   ambient <S>  : one sample: extra RNA of the suspect zoom clusters over the clean ones vs the ambient profile
##                  (empty droplets, 1-100 counts) and the sample's own contaminant cells
##                  (non-negative fit: extra = a * ambient + b * contaminant); one process per sample
##   joint        : zoom cells + contaminant cells in one Harmony map: continuum (bridge) or separate cloud?
##   summary      : per-sample table and figures from the ambient + joint outputs
##   Order: check, profile, zoom (twice: grid, then chosen resolution), signatures, ambient <S> per sample, joint, summary.
##   check/profile/ambient/joint/summary are optional: run them when a group looks like doublets or contamination.
## Reads : objects/<input_object>.rds ; objects/subpop_<study>.rds (signatures) ; the sample's raw/filtered
##         matrices from samples.csv (ambient)
## Writes: 17_subpopulation_evidence/<study>/<mode>/ ; objects/subpop_<study>_grid.rds, objects/subpop_<study>.rds ;
##         annotation/<study>.tsv template
## Checkpoint: subpopulations.<study>.{celltypes, focus_celltype, compare_celltypes, contaminant_celltype,
##             resolution, suspect_groups, clean_groups} (only the ones a mode needs); then annotation/<study>.tsv (step 18)
##
file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
PIPELINE_DIR <- normalizePath(file.path(dirname(file_arg), ".."))
source(file.path(PIPELINE_DIR, "R/setup.R"))
suppressPackageStartupMessages(library(tidyr))

# ---- Marker lists ------------------------------------------------------------------------------------------
# PREFERRED: keep the lists in resources/marker_panels.yaml and point to them from config.yaml
# (subpopulations.<study>.panels: "panel" or "panel.block"). They are then shared, versioned and reused.
# QUICK ALTERNATIVE: write a list here; a non-NULL entry overrides the yaml panel for that study.
# Shapes: lineage_genes / foreign / coexpression / zoom_panel / signatures = named list of gene vectors;
#         own / gate_genes = gene vector.
MANUAL_PANELS <- list(
  # caf = list(
  #   own     = c("COL1A1", "COL1A2", "DCN", "LUM", "PDGFRA", "FAP"),
  #   foreign = list(macrophage = c("C1QA", "C1QB", "CD163", "CD68"), epithelial = c("EPCAM", "KRT8", "KRT18"))
  # )
)
# -------------------------------------------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE); study <- args[1]; mode <- args[2]
modes <- c("check", "profile", "zoom", "signatures", "ambient", "joint", "summary")
if (is.na(study) || is.na(mode) || !mode %in% modes || (mode == "ambient" && length(args) < 3)) {
  stop("Usage: 17_subpopulation_evidence.R <study> check | profile | zoom | signatures | ambient <S> | joint | summary")
}
key <- function(k) paste0("subpopulations.", study, ".", k)
if (is.null(tryCatch(cfg_get(paste0("subpopulations.", study)), error = function(e) NULL))) {
  stop("No study '", study, "' in config.yaml (subpopulations.", study, ")")
}
step <- "17_subpopulation_evidence"
out_dir <- step_dir(step, study, if (mode %in% c("ambient", "joint", "summary")) "contamination" else mode)
zoom_meta_file <- paste0(step_dir(step, study, "zoom"), "zoom_meta.rds")
labels_file    <- paste0(step_dir(step, study, "zoom"), "input_labels.csv.gz")
n_pcs <- cfg_get(key("n_pcs"))

# Decisions are read lazily: each mode only stops for the ones it uses
zoom_types  <- function() as.character(get_decision(key("celltypes")))
focus       <- function() as.character(get_decision(key("focus_celltype")))
compare     <- function() as.character(get_decision(key("compare_celltypes")))
contaminant <- function() as.character(get_decision(key("contaminant_celltype")))

# Marker list: MANUAL_PANELS first, then the yaml panel named in config; NULL if neither exists
study_panel <- function(name, object = NULL, required = TRUE) {
  manual <- MANUAL_PANELS[[study]][[name]]
  if (!is.null(manual)) {
    cat("  [", name, "] using the manual list in 17_subpopulation_evidence.R\n", sep = "")
    p <- manual
  } else {
    ref <- tryCatch(cfg_get(key(paste0("panels.", name))), error = function(e) NULL)
    if (is.null(ref)) {
      if (required) stop("Study '", study, "': no marker list '", name, "'. Add it to resources/marker_panels.yaml and ",
                         key(paste0("panels.", name)), " in config.yaml (preferred), or to MANUAL_PANELS in this script")
      return(NULL)
    }
    parts <- strsplit(ref, ".", fixed = TRUE)[[1]]
    p <- load_marker_panel(parts[1], if (length(parts) > 1) parts[2] else NULL)
  }
  if (is.null(object)) return(p)
  keep <- function(x) if (is.list(x)) lapply(x, keep) else intersect(x, rownames(object))
  keep(p)
}
own_genes     <- function() study_panel("own")
foreign_sets  <- function() study_panel("foreign")
contam_set    <- function() { f <- foreign_sets(); s <- cfg_get(key("contaminant_set"))
                              if (!s %in% names(f)) stop(key("contaminant_set"), " = '", s, "' is not a set of panels.foreign"); f[[s]] }

# Subset, joined counts in memory, Harmony per sample. generate_subset is not used: it selects variable genes
# per sample and fails on samples with very few cells of the zoomed types (e.g. a sample with ~10 fibroblasts).
joined_harmony <- function(object, cells, n_pcs) {
  z <- JoinLayers(subset(object, cells = cells))
  z[["RNA"]]$counts <- as(z[["RNA"]]$counts, "dgCMatrix")
  z[["RNA"]]$scale.data <- NULL   # parent's layer, other genes; rebuilt below
  z <- NormalizeData(z, verbose = FALSE); z <- FindVariableFeatures(z, verbose = FALSE)
  z <- ScaleData(z, verbose = FALSE); z <- RunPCA(z, verbose = FALSE)
  z <- harmony::RunHarmony(z, group.by.vars = cfg_get("integration.batch_key"), verbose = FALSE)
  z <- FindNeighbors(z, reduction = "harmony", dims = 1:n_pcs, verbose = FALSE)
  run_umap_pca(z, reduction = "harmony", dims = 1:n_pcs)
}

# Zoom cluster -> contamination-check group (suspect / clean / other)
group_labels <- function() {
  sus <- as.character(get_decision(key("suspect_groups"))); cln <- as.character(get_decision(key("clean_groups")))
  list(suspect = sus, clean = cln,
       lab_suspect = sprintf("suspect (groups %s)", paste(sus, collapse = ", ")),
       lab_clean = sprintf("clean (groups %s)", paste(cln, collapse = ", ")))
}
zoom_group <- function(sub, g) ifelse(sub %in% g$suspect, g$lab_suspect, ifelse(sub %in% g$clean, g$lab_clean, "other groups"))

# Per-cell number of genes detected (counts > 0) of each gene set
genes_detected <- function(cnt, sets) sapply(sets, function(gs) Matrix::colSums(cnt[intersect(gs, rownames(cnt)), , drop = FALSE] > 0))

# ======================================================================================
if (mode == "check") {
  a <- read_run_object(cfg_get(key("input_object")))
  types <- intersect(unique(c(focus(), zoom_types(), compare())), unique(as.character(a$celltype)))
  set.seed(1)
  cells <- unlist(lapply(types, function(t) { x <- colnames(a)[a$celltype == t]; if (length(x) > 3000) sample(x, 3000) else x }))
  s <- subset(a, cells = cells); rm(a); invisible(gc())
  s <- JoinLayers(s); s <- NormalizeData(s, verbose = FALSE)
  s$celltype <- factor(as.character(s$celltype), levels = rev(types))
  p <- dot_plot(s, features = study_panel("lineage_genes", s), group.by = "celltype") + RotatedAxis() +
    labs(x = NULL, y = NULL, title = paste(focus(), "vs the compare types and the real owners of each gene (<= 3,000 cells per type)"))
  ggsave(paste0(out_dir, "DotPlot_", focus(), "_lineage_genes.png"), p, width = 15, height = 5.5, dpi = 150, bg = "white")
  # same-cell co-expression (counts > 0): an own-lineage gene AND a foreign gene
  co_sets <- study_panel("coexpression", required = FALSE) %||% c(list(own = own_genes()), foreign_sets())
  det <- genes_detected(GetAssayData(s, layer = "counts"), co_sets) > 0
  co <- data.frame(celltype = s$celltype, nFeature_RNA = s$nFeature_RNA)
  for (f in names(co_sets)[-1]) co[[paste0(names(co_sets)[1], "_and_", f, "_pct")]] <- det[, 1] & det[, f]
  co <- co %>% group_by(celltype) %>% summarise(cells = n(), across(ends_with("_pct"), ~ round(100 * mean(.x), 1)),
                                                median_genes = median(nFeature_RNA))
  write.csv(co, paste0(out_dir, "coexpression_by_celltype.csv"), row.names = FALSE); print(co)
}

# ======================================================================================
if (mode == "profile") {
  a <- read_run_object(cfg_get(key("input_object")))
  meta <- a@meta.data; meta$celltype <- as.character(meta$celltype)
  zt <- zoom_types(); fc <- focus()
  # 1. origin: focus label per sample, and as a share of each sample's zoomed cells
  bys <- meta %>% filter(celltype %in% zt) %>% count(orig.ident, celltype) %>%
    pivot_wider(names_from = celltype, values_from = n, values_fill = 0)
  tot <- rowSums(bys[, intersect(zt, colnames(bys)), drop = FALSE])
  bys$pct_focus_of_zoomed <- round(100 * bys[[fc]] / tot, 1)
  bys$pct_of_all_focus <- round(100 * bys[[fc]] / sum(bys[[fc]]), 1)
  bys <- bys[order(-bys[[fc]]), ]
  write.csv(bys, paste0(out_dir, "by_sample.csv"), row.names = FALSE); print(bys, n = 30)
  # 2. quality per type (QC columns that exist)
  qc <- intersect(c("nFeature_RNA", "nCount_RNA", "percent.mt", "decontX_contamination", "scDblFinder.score", "doublet_score"), colnames(meta))
  keep <- intersect(unique(c(fc, zt, compare())), unique(meta$celltype))
  q <- meta %>% filter(celltype %in% keep) %>% select(celltype, all_of(qc)) %>% pivot_longer(-celltype)
  q$celltype <- factor(q$celltype, levels = keep)
  ggsave(paste0(out_dir, "Violin_quality.png"),
         ggplot(q, aes(celltype, value, fill = celltype)) + geom_violin(scale = "width", linewidth = 0.2) +
           geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white") + facet_wrap(~name, scales = "free_y", nrow = 1) +
           theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") + labs(x = NULL, y = NULL),
         width = 4 + 3 * length(qc), height = 4.5, dpi = 150, bg = "white")
  write.csv(meta %>% filter(celltype %in% keep) %>% group_by(celltype) %>% summarise(across(all_of(qc), median), cells = n()),
            paste0(out_dir, "quality_median.csv"), row.names = FALSE)
  # 3. markers focus vs each other zoomed label (both directions), pipeline filter
  s <- subset(a, cells = colnames(a)[as.character(a$celltype) %in% c(zt, contaminant())]); rm(a); invisible(gc())
  s <- JoinLayers(s); s <- NormalizeData(s, verbose = FALSE); Idents(s) <- "celltype"
  for (ref in setdiff(zt, fc)) {
    m <- FindMarkers(s, ident.1 = fc, ident.2 = ref, logfc.threshold = 0.25, min.pct = 0.1, verbose = FALSE)
    m$gene <- rownames(m); write.csv(m, paste0(out_dir, "markers_", fc, "_vs_", ref, ".csv"), row.names = FALSE)
    up <- m %>% filter(p_val_adj < cfg_get("markers.padj"), avg_log2FC > 0, pct.1 >= cfg_get("markers.min_pct1")) %>% arrange(desc(avg_log2FC)) %>% head(40)
    dn <- m %>% filter(p_val_adj < cfg_get("markers.padj"), avg_log2FC < 0, pct.2 >= cfg_get("markers.min_pct1")) %>% arrange(avg_log2FC) %>% head(25)
    cat("\n==", fc, "UP vs", ref, ":", paste(up$gene, collapse = " "), "\n== DOWN vs", ref, ":", paste(dn$gene, collapse = " "), "\n")
  }
  # 4. foreign genes per cell: doublet (many genes) vs ambient (1-2 genes)
  fs <- foreign_sets(); det <- genes_detected(GetAssayData(s, layer = "counts"), fs)
  d <- do.call(rbind, lapply(names(fs), function(f) data.frame(set = sprintf("%s genes (%d)", f, length(fs[[f]])), genes = det[, f],
                                                               celltype = as.character(s$celltype))))
  d$celltype <- factor(d$celltype, levels = unique(c(fc, zt, contaminant())))
  d$class <- cut(d$genes, c(-Inf, 0, 1, 2, Inf), labels = c("none", "1 gene", "2 genes", "3+ genes"))
  tab <- d %>% count(set, celltype, class) %>% group_by(set, celltype) %>% mutate(pct = 100 * n / sum(n))
  write.csv(tab, paste0(out_dir, "foreign_genes_per_cell.csv"), row.names = FALSE); print(as.data.frame(tab))
  ggsave(paste0(out_dir, "Bar_foreign_genes_per_cell.png"),
         ggplot(tab, aes(celltype, pct, fill = class)) + geom_col(width = 0.7) + facet_wrap(~set) +
           scale_fill_manual(values = c("none" = "#e6e6e3", "1 gene" = "#c6b8f0", "2 genes" = "#8f78e0", "3+ genes" = "#4b2fa8")) +
           theme_classic() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
           labs(x = NULL, y = "% of cells", fill = "foreign genes\ndetected in the cell"),
         width = 4 + 3 * length(fs), height = 4.5, dpi = 150, bg = "white")
}

# ======================================================================================
if (mode == "zoom") {
  grid <- cfg_get("clustering.resolution_grid")
  zt <- zoom_types()
  cache <- run_object(paste0("subpop_", study, "_grid"))
  z <- if (file.exists(cache)) readRDS(cache) else NULL
  if (!is.null(z) && !identical(sort(as.character(z@misc$zoom_celltypes)), sort(zt))) {
    cat("  Cached zoom was built from other cell types: recomputing\n"); z <- NULL
  }
  if (is.null(z)) {
    a <- read_run_object(cfg_get(key("input_object")))
    # light per-cell label table for the per-sample modes (they must not load the full object 16 times)
    write.csv(data.frame(cell = colnames(a), sample = a$orig.ident, celltype = as.character(a$celltype)),
              gzfile(labels_file), row.names = FALSE)
    cells <- colnames(a)[as.character(a$celltype) %in% zt]
    cat(sprintf("  Zoom '%s': %d cells (%s)\n", study, length(cells), paste(zt, collapse = ", ")))
    z <- joined_harmony(a, cells, n_pcs); rm(a); invisible(gc())
    z <- FindClusters(z, resolution = grid, verbose = FALSE)
    z@misc$zoom_celltypes <- zt
    saveRDS(z, cache)
    plot_resolution_grid(z, out_dir, resolutions = grid)
    sizes <- do.call(rbind, lapply(grid, function(r) {
      tb <- table(z[[paste0("RNA_snn_res.", r)]][, 1]); data.frame(resolution = r, cluster = names(tb), cells = as.integer(tb))
    }))
    write.csv(sizes, paste0(out_dir, "cluster_sizes_per_resolution.csv"), row.names = FALSE)
    cat("  Grid cache ->", cache, "\n")
  }
  res <- get_decision(key("resolution"))
  res_col <- paste0("RNA_snn_res.", res)
  if (!res_col %in% colnames(z@meta.data)) z <- FindClusters(z, resolution = res, verbose = FALSE)
  z$seurat_clusters <- factor(z[[res_col]][, 1], levels = sort(as.integer(unique(as.character(z[[res_col]][, 1])))))
  Idents(z) <- "seurat_clusters"
  # foreign genes carried by each cell (doublet-like = 3+)
  fs <- foreign_sets(); det <- genes_detected(z[["RNA"]]$counts, fs)
  for (f in names(fs)) z[[paste0(f, "_genes")]] <- det[, f]
  md <- z@meta.data; md$celltype <- as.character(md$celltype)
  write.csv(md %>% count(seurat_clusters, celltype) %>% pivot_wider(names_from = celltype, values_from = n, values_fill = 0),
            paste0(out_dir, "sub_by_celltype.csv"), row.names = FALSE)
  write.csv(md %>% count(seurat_clusters, orig.ident) %>% pivot_wider(names_from = orig.ident, values_from = n, values_fill = 0),
            paste0(out_dir, "sub_by_sample.csv"), row.names = FALSE)
  fcols <- paste0(names(fs), "_genes")
  dep <- md %>% group_by(seurat_clusters) %>%
    summarise(cells = n(), n_samples = n_distinct(orig.ident), top_sample = names(which.max(table(orig.ident))),
              top_sample_pct = round(100 * max(table(orig.ident)) / n(), 0), median_genes = median(nFeature_RNA),
              median_counts = median(nCount_RNA), pct_mt = median(percent.mt), decontX = round(median(decontX_contamination), 3),
              across(all_of(fcols), ~ round(100 * mean(.x >= 3), 1), .names = "pct_{.col}_3plus"),
              across(celltype, ~ names(which.max(table(.x))), .names = "main_celltype"))
  write.csv(dep, paste0(out_dir, "sub_depth.csv"), row.names = FALSE); print(dep)
  um <- data.frame(Embeddings(z, "umap")); colnames(um) <- c("u1", "u2")
  um$sub <- z$seurat_clusters; um$sample <- z$orig.ident; um$celltype <- md$celltype
  ggsave(paste0(out_dir, "UMAP_previous_celltype.png"), ggplot(um, aes(u1, u2, colour = celltype)) + geom_point(size = 0.2) + theme_void() +
           guides(colour = guide_legend(override.aes = list(size = 3))), width = 7.5, height = 6, bg = "white")
  ggsave(paste0(out_dir, "UMAP_sub.png"), ggplot(um, aes(u1, u2, colour = sub)) + geom_point(size = 0.2) + theme_void() +
           geom_text(data = um %>% group_by(sub) %>% summarise(u1 = median(u1), u2 = median(u2)), aes(label = sub), colour = "black") +
           guides(colour = guide_legend(override.aes = list(size = 3))), width = 7, height = 6, bg = "white")
  ggsave(paste0(out_dir, "UMAP_sample.png"), ggplot(um, aes(u1, u2, colour = sample)) + geom_point(size = 0.2) + theme_void() +
           guides(colour = guide_legend(override.aes = list(size = 3))), width = 7.5, height = 6, bg = "white")
  zp <- study_panel("zoom_panel", z, required = FALSE)
  if (!is.null(zp)) ggsave(paste0(out_dir, "DotPlot_panel.png"),
                           dot_plot(z, features = unique(unlist(zp)), group.by = "seurat_clusters") + RotatedAxis(),
                           width = 17, height = 3 + 0.35 * nlevels(z$seurat_clusters), bg = "white")
  mf <- find_markers_for_clusters(z, levels(z$seurat_clusters), out_dir)
  top_expressed_genes(z, "seurat_clusters", out_dir)
  saveRDS(cbind(z@meta.data[, c("orig.ident", "celltype", "seurat_clusters", fcols, "nFeature_RNA", "nCount_RNA", "decontX_contamination")],
                Embeddings(z, "umap")), zoom_meta_file)   # per-cell groups for the follow-up modes
  saveRDS(z, run_object(paste0("subpop_", study)))
  write_annotation_template(z, study, markers_df = mf)
}

# ======================================================================================
if (mode == "signatures") {
  z <- read_run_object(paste0("subpop_", study))
  sig <- study_panel("signatures", z); print(sig)
  z$grp <- factor(as.character(z$seurat_clusters), levels = rev(levels(z$seurat_clusters)))
  ggsave(paste0(out_dir, "DotPlot_signature_genes.png"),
         dot_plot(z, features = sig, group.by = "grp") + RotatedAxis() + labs(x = NULL, y = "zoom cluster", title = "Zoom clusters vs published sub-type markers"),
         width = 19, height = 5.5, dpi = 150, bg = "white")
  z <- AddModuleScore(z, features = sig, name = "Sig", seed = 1)
  sc_cols <- paste0("Sig", seq_along(sig))
  m <- z@meta.data %>% group_by(grp) %>% summarise(across(all_of(sc_cols), mean)) %>% as.data.frame()
  rownames(m) <- m$grp; m$grp <- NULL; colnames(m) <- names(sig)
  mz <- scale(m)   # z-score of each signature across clusters
  d <- data.frame(group = rep(rownames(mz), ncol(mz)), signature = rep(colnames(mz), each = nrow(mz)), z = as.vector(mz))
  d$group <- factor(d$group, levels = levels(z$grp)); d$signature <- factor(d$signature, levels = names(sig))
  ggsave(paste0(out_dir, "Heatmap_signature_scores.png"),
         ggplot(d, aes(signature, group, fill = z)) + geom_tile(colour = "white") + geom_text(aes(label = sprintf("%.1f", z)), size = 3) +
           scale_fill_gradient2(low = "#2a78d6", mid = "white", high = "#e34948") + theme_minimal() +
           theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
           labs(x = NULL, y = "zoom cluster", fill = "score\n(z across\nclusters)", title = "Which published signature fits each zoom cluster best"),
         width = 9, height = 5, dpi = 150, bg = "white")
  write.csv(cbind(group = rownames(m), round(m, 3)), paste0(out_dir, "signature_scores_by_group.csv"), row.names = FALSE)
  print(round(mz, 2))
  gg <- study_panel("gate_genes", z, required = FALSE)
  if (length(gg)) {   # protein-gate analogue at RNA level: % of cells with counts > 0 per cluster
    cnt <- z[["RNA"]]$counts
    pp <- sapply(gg, function(x) tapply(cnt[x, ] > 0, z$grp, mean) * 100)
    write.csv(round(pp, 1), paste0(out_dir, "pct_cells_gate_genes.csv")); print(round(pp, 1))
  }
}

# ======================================================================================
if (mode == "ambient") {
  s <- args[3]
  smp <- load_samples(); smp <- smp[smp$sample == s, ]
  if (!nrow(smp)) stop("Sample ", s, " not in samples.csv")
  dir.create(paste0(out_dir, "per_sample"), showWarnings = FALSE)
  g <- group_labels(); ct <- contaminant(); fs <- foreign_sets(); cs <- contam_set()
  if (!file.exists(labels_file)) stop("Missing ", labels_file, ": run '", study, " zoom' first")
  lab <- read.csv(gzfile(labels_file)); lab <- lab[lab$sample == s, ]; lab$bc <- sub("_[0-9]+$", "", lab$cell)
  tod <- Read10X(smp$raw_dir); flt <- Read10X(smp$filtered_dir)
  tod <- tod[intersect(rownames(flt), rownames(tod)), ]
  flt <- colnames(flt)
  tot <- Matrix::colSums(tod)
  amb_cols <- !(colnames(tod) %in% flt) & tot >= 1 & tot <= 100
  amb <- Matrix::rowSums(tod[, amb_cols])
  zm <- readRDS(zoom_meta_file); zs <- zm[zm$orig.ident == s, ]; zs$bc <- sub("_[0-9]+$", "", rownames(zs))
  zs$sub <- as.character(zs$seurat_clusters)
  prof <- function(bcs) { bcs <- intersect(bcs, colnames(tod)); if (length(bcs) < 10) return(NULL); x <- Matrix::rowSums(tod[, bcs, drop = FALSE]); x / sum(x) }
  ps <- prof(zs$bc[zs$sub %in% g$suspect]); pc <- prof(zs$bc[zs$sub %in% g$clean]); pt <- prof(lab$bc[lab$celltype == ct])
  pa <- amb / sum(amb)
  frac <- function(p, gs) if (is.null(p)) NA else 100 * sum(p[intersect(gs, names(p))])
  res <- data.frame(sample = s, cells = nrow(lab), contaminant_cells = sum(lab$celltype == ct), zoomed = nrow(zs),
                    suspect = sum(zs$sub %in% g$suspect), ambient_droplets = sum(amb_cols), ambient_counts = sum(amb),
                    suspect_median_decontX = median(zs$decontX_contamination[zs$sub %in% g$suspect]))
  for (f in names(fs)) {   # % of each RNA profile that comes from each foreign gene set
    res[[paste0("ambient_", f, "_pct")]] <- frac(pa, fs[[f]]); res[[paste0("suspect_", f, "_pct")]] <- frac(ps, fs[[f]])
    res[[paste0("clean_", f, "_pct")]] <- frac(pc, fs[[f]]);   res[[paste0("contaminant_", f, "_pct")]] <- frac(pt, fs[[f]])
  }
  res[c("a_ambient", "b_contaminant", "fit_r2", "cor_extra_ambient", "cor_extra_contaminant")] <- NA
  if (!is.null(ps) && !is.null(pc) && !is.null(pt)) {
    extra <- pmax(ps - pc, 0)                                   # RNA the suspect groups carry on top of the clean ones
    gx <- names(sort(extra, decreasing = TRUE))[1:300]          # genes carrying the extra RNA
    f <- nnls::nnls(cbind(pa[gx], pt[gx]), extra[gx])
    res$a_ambient <- f$x[1]; res$b_contaminant <- f$x[2]; res$fit_r2 <- 1 - sum(f$residuals^2) / sum((extra[gx] - mean(extra[gx]))^2)
    res$cor_extra_ambient <- cor(extra[gx], pa[gx], method = "spearman"); res$cor_extra_contaminant <- cor(extra[gx], pt[gx], method = "spearman")
    write.csv(data.frame(gene = gx, extra = extra[gx], ambient = pa[gx], contaminant = pt[gx]), paste0(out_dir, "per_sample/", s, "_extra_genes.csv"), row.names = FALSE)
  }
  write.csv(res, paste0(out_dir, "per_sample/", s, "_ambient.csv"), row.names = FALSE); print(t(res))
}

# ======================================================================================
if (mode == "joint") {
  g <- group_labels(); ct <- contaminant(); zt <- zoom_types()
  zm <- readRDS(zoom_meta_file)
  a <- read_run_object(cfg_get(key("input_object")))
  write.csv(data.frame(sample = a$orig.ident, decontX = a$decontX_contamination) %>% group_by(sample) %>%
              summarise(sample_median_decontX = median(decontX), sample_cells = n()),
            paste0(out_dir, "sample_decontX.csv"), row.names = FALSE)
  cells <- colnames(a)[as.character(a$celltype) %in% c(zt, ct)]
  z <- joined_harmony(a, cells, n_pcs); rm(a); invisible(gc())
  grp <- ifelse(as.character(z$celltype) == ct, ct, zoom_group(as.character(zm$seurat_clusters[match(colnames(z), rownames(zm))]), g))
  dat <- GetAssayData(z, layer = "data")
  sc <- function(gs) Matrix::colMeans(dat[intersect(gs, rownames(dat)), , drop = FALSE])
  OWN <- own_genes(); CON <- contam_set()
  d <- data.frame(cell = colnames(z), group = grp, own = sc(OWN), contam = sc(CON), Embeddings(z, "umap"))
  colnames(d)[5:6] <- c("u1", "u2")
  # neighbours in the Harmony space: share of contaminant cells among the 20 nearest cells
  nn <- RANN::nn2(Embeddings(z, "harmony")[, 1:n_pcs], k = 21)$nn.idx[, -1]
  d$contaminant_neighbours <- rowMeans(matrix(d$group[nn] == ct, nrow = nrow(nn)))
  saveRDS(d, paste0(out_dir, "joint_cells.rds"))
  lv <- c(g$lab_clean, "other groups", g$lab_suspect, ct)
  cols <- setNames(c("#1baf7a", "#a9a9a4", "#e34948", "#eda100"), lv); d$group <- factor(d$group, levels = lv)
  set.seed(1); d <- d[sample(nrow(d)), ]
  ggsave(paste0(out_dir, "Joint_UMAP_groups.png"), ggplot(d, aes(u1, u2, colour = group)) + geom_point(size = 0.15) + scale_colour_manual(values = cols) +
           theme_void() + guides(colour = guide_legend(override.aes = list(size = 3))) +
           labs(title = paste("All", paste(zt, collapse = "/"), "+ all", ct, "(Harmony per sample)")),
         width = 8.5, height = 6.5, bg = "white")
  ggsave(paste0(out_dir, "Joint_scores_own_vs_contaminant.png"), ggplot(d, aes(own, contam)) + geom_point(size = 0.15, alpha = 0.4, aes(colour = group)) +
           scale_colour_manual(values = cols) + facet_wrap(~group, nrow = 1) + theme_classic() + theme(legend.position = "none") +
           labs(x = sprintf("own-lineage score (mean of %d genes)", length(OWN)), y = sprintf("%s score (mean of %d genes)", ct, length(CON)),
                title = paste("Each cell's own vs", ct, "programme: a bridge would fill the space between the two clouds")),
         width = 17, height = 4.2, bg = "white")
  ggsave(paste0(out_dir, "Joint_contaminant_neighbours.png"), ggplot(d, aes(group, contaminant_neighbours, fill = group)) + geom_violin(scale = "width") +
           scale_fill_manual(values = cols) + theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1)) +
           labs(x = NULL, y = paste("share of", ct, "among the 20 nearest cells")), width = 8, height = 4.5, bg = "white")
  print(d %>% group_by(group) %>% summarise(cells = n(), own = median(own), contam = median(contam),
                                            nb = round(median(contaminant_neighbours), 2), nb_gt50 = round(100 * mean(contaminant_neighbours > 0.5), 1)))
}

# ======================================================================================
if (mode == "summary") {
  ct <- contaminant(); cs_name <- cfg_get(key("contaminant_set"))
  r <- do.call(rbind, lapply(list.files(paste0(out_dir, "per_sample"), pattern = "_ambient.csv$", full.names = TRUE), read.csv))
  if (is.null(r)) stop("No per-sample ambient tables: run '", study, " ambient <S>' for every sample first")
  sdx <- read.csv(paste0(out_dir, "sample_decontX.csv")); r$sample_median_decontX <- sdx$sample_median_decontX[match(r$sample, sdx$sample)]
  r <- r %>% mutate(suspect_pct_of_zoomed = 100 * suspect / zoomed, contaminant_pct_of_cells = 100 * contaminant_cells / cells)
  write.csv(r, paste0(out_dir, "per_sample_summary.csv"), row.names = FALSE)
  amb_col <- paste0("ambient_", cs_name, "_pct")
  long <- r %>% select(sample, suspect_pct_of_zoomed, contaminant_pct_of_cells, sample_median_decontX, all_of(amb_col)) %>%
    pivot_longer(c(contaminant_pct_of_cells, sample_median_decontX, all_of(amb_col)))
  long$name <- factor(long$name, levels = c("contaminant_pct_of_cells", "sample_median_decontX", amb_col),
                      labels = c(paste0(ct, ", % of the sample's cells"), "median DecontX contamination of the sample", paste(cs_name, "genes, % of the ambient RNA")))
  rho <- long %>% group_by(name) %>% summarise(rho = cor(value, suspect_pct_of_zoomed, method = "spearman"))
  ggsave(paste0(out_dir, "Per_sample_suspect_vs_sample.png"),
         ggplot(long, aes(value, suspect_pct_of_zoomed, label = sample)) + geom_point(colour = "#e34948", size = 2) + geom_text(size = 2.6, vjust = -0.8) +
           facet_wrap(~name, scales = "free_x", labeller = labeller(name = function(x) paste0(x, "\nSpearman rho = ", round(rho$rho[match(x, rho$name)], 2)))) +
           theme_classic() + labs(x = NULL, y = "suspect groups, % of the sample's zoomed cells"), width = 13, height = 4.6, bg = "white")
  print(r %>% select(sample, suspect, zoomed, contaminant_cells, matches("_pct$"), a_ambient, b_contaminant, fit_r2, cor_extra_ambient, cor_extra_contaminant) %>%
          mutate(across(where(is.numeric), ~ round(.x, 3))))
  print(rho)
  comp <- r %>% filter(!is.na(.data[[paste0("suspect_", cs_name, "_pct")]]), suspect >= 20) %>%
    select(sample, matches("^(ambient|clean|suspect|contaminant)_.*_pct$"))
  cl <- comp %>% pivot_longer(-sample) %>%
    mutate(source = sub("_.*", "", name), genes = sub("^[a-z]+_(.*)_pct$", "\\1", name),
           source = factor(source, levels = c("ambient", "clean", "suspect", "contaminant"),
                           labels = c("ambient RNA (empty droplets)", "clean groups", "suspect groups", ct)),
           genes = paste(genes, "genes"))
  ggsave(paste0(out_dir, "Composition_ambient_vs_suspect.png"),
         ggplot(cl, aes(source, value, colour = source)) + geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.15, size = 1.4) +
           facet_wrap(~genes, scales = "free_y") + theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1)) +
           labs(x = NULL, y = "% of the RNA (one dot = one sample with >= 20 cells of the group)"), width = 11, height = 4.8, bg = "white")
}

cat(sprintf("\n17_subpopulation_evidence %s %s done -> %s\n", study, mode, out_dir))
