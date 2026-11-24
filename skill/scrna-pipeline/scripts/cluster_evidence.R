#!/usr/bin/env Rscript
##
## cluster_evidence.R - one evidence table per annotation dictionary, so the annotator (human or agent) reads
## every cluster's facts side by side before writing labels. It does NOT label anything.
##
## Usage: Rscript pipeline/skill/scrna-pipeline/scripts/cluster_evidence.R <dictionary> [object] [cluster_col] [markers_dir]
##   dictionary : lineage | leukocytes | stroma | tumour | tumour_lp | tumour_lowrna | <sub-population study, e.g. caf>
##   object, cluster_col, markers_dir : only needed for a dictionary not in the table below
## Writes: results/<run>/evidence/<dictionary>.md (and prints it); never touches annotation/
##
## Columns per cluster:
##   cells, samples (n samples; top sample and its %), median genes, median DecontX contamination,
##   lin3_<block>: % of cells with >= 3 genes detected of each block of the "lineage_specific" panel
##                 (a cluster with high values for TWO lineages = doublets/multiplets; 1-2 genes = ambient),
##   markers: top filtered FindMarkers genes (pct.1 in brackets), expressed: top mean-expressed genes,
##   tumour only: majority PAM50 / scSubtype / cell-cycle phase (with %) and CopyKAT call when present.
##
args <- commandArgs(trailingOnly = TRUE)
file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
PIPELINE_DIR <- normalizePath(file.path(dirname(file_arg), "..", "..", ".."))
source(file.path(PIPELINE_DIR, "R/setup.R"))
dict <- args[1]
if (is.na(dict)) stop("Usage: cluster_evidence.R <dictionary> [object] [cluster_col] [markers_dir]")

known <- list(
  lineage       = list(obj = "clustered",     col = "global_clusters", mk = "04_lineage_markers"),
  leukocytes    = list(obj = "leukocytes",    col = "seurat_clusters", mk = "06_leukocytes_markers"),
  stroma        = list(obj = "stroma",        col = "seurat_clusters", mk = "08_stroma_markers"),
  tumour        = list(obj = "tumour_scored", col = "seurat_clusters", mk = "11_tumour_markers"),
  tumour_lp     = list(obj = "tumour_lp",     col = "seurat_clusters", mk = "15_tumour_copykat_evidence/zoom_tumour_lp"),
  tumour_lowrna = list(obj = "tumour_lowrna", col = "seurat_clusters", mk = "15_tumour_copykat_evidence/zoom_tumour_lowrna"))
spec <- known[[dict]] %||% list(obj = paste0("subpop_", dict), col = "seurat_clusters",
                                mk = file.path("17_subpopulation_evidence", dict, "zoom"))
if (length(args) >= 2) spec$obj <- args[2]
if (length(args) >= 3) spec$col <- args[3]
if (length(args) >= 4) spec$mk  <- args[4]

o <- read_run_object(spec$obj)
cl <- as.character(o[[spec$col, drop = TRUE]])
md <- o@meta.data
ord <- unique(cl); ord <- ord[order(suppressWarnings(as.numeric(ord)), ord)]

# % of cells with >= 3 genes of each lineage-specific block
blocks <- load_marker_panel("lineage_specific", object = o)
cnt <- LayerData(o, assay = "RNA", layer = "counts")
lin3 <- sapply(blocks, function(b) {
  g <- unique(unlist(b)); m <- as(cnt[g, , drop = FALSE], "dgCMatrix")
  Matrix::colSums(m > 0) >= 3
})
mk_file <- file.path(run_dir(), spec$mk, "all_clusters_FindMarkers_filtered.csv")
mk <- if (file.exists(mk_file)) read.csv(mk_file) else NULL
te_file <- file.path(run_dir(), spec$mk, "top_expressed_genes.csv")
te <- if (file.exists(te_file)) read.csv(te_file) else NULL
maj <- function(x) { t <- sort(table(x), decreasing = TRUE); if (!length(t)) return("-"); sprintf("%s %.0f%%", names(t)[1], 100 * t[1] / sum(t)) }

rows <- lapply(ord, function(k) {
  i <- cl == k; s <- table(md$orig.ident[i])
  r <- list(cluster = k, cells = sum(i),
            samples = sprintf("%d (top %s %.0f%%)", length(s[s > 0]), names(which.max(s)), 100 * max(s) / sum(i)),
            genes = round(median(md$nFeature_RNA[i])),
            decontX = if ("decontX_contamination" %in% colnames(md)) sprintf("%.2f", median(md$decontX_contamination[i])) else "-")
  for (b in colnames(lin3)) r[[paste0("lin3_", b)]] <- sprintf("%.0f", 100 * mean(lin3[i, b]))
  if (!is.null(mk)) { m <- mk[as.character(mk$cluster) == k, ]; m <- head(m[order(-m$avg_log2FC), ], 10)
    r$markers <- paste(sprintf("%s(%.0f)", m$gene, 100 * m$pct.1), collapse = " ") }
  if (!is.null(te)) r$expressed <- paste(head(te$gene[as.character(te$cluster) == k], 6), collapse = " ")
  for (col in intersect(c("PAM50_predicted", "scSubtype_call", "Phase", "copykat_prediction"), colnames(md))) r[[col]] <- maj(md[[col]][i])
  r
})
tab <- do.call(rbind, lapply(rows, as.data.frame))

out <- c(sprintf("# Evidence for dictionary '%s' (object %s, column %s, run %s)", dict, spec$obj, spec$col, cfg_get("run_name")), "",
         "lin3_<block> = % of cells with >= 3 detected genes of that lineage block (two high blocks = doublets/multiplets).",
         "markers = filtered FindMarkers (padj, pct.1 thresholds of config), pct.1 in brackets.", "",
         paste("|", paste(colnames(tab), collapse = " | "), "|"),
         paste("|", paste(rep("---", ncol(tab)), collapse = " | "), "|"),
         apply(tab, 1, function(x) paste("|", paste(x, collapse = " | "), "|")))
f <- paste0(step_dir("evidence"), dict, ".md")
writeLines(out, f)
cat(out, sep = "\n"); cat("\n->", f, "\n")
