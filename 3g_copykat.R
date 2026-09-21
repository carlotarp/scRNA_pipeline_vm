##
##  Single Cell Analysis Step 3g: CopyKAT CNV annotation
##  Runs AFTER 3f_stroma_annotate.R — takes notumor_annotated_data.rds as input.
##  Runs copyKAT once per sample (fixed reference cells pooled across all
##  samples), attaches ploidy predictions, the CNA matrix, the chromosomal
##  instability scores (CIN/PGA) and a CNA-based UMAP.
##

# Import libraries
library(Seurat)
library(copykat)
library(dplyr)
library(ggplot2)
library(ggsci)

# Set paths
project_path <- "/home/user/PROJECTS/scRNA_vmendez/"
wd <- paste0(project_path, "codes/")
setwd(wd)
results_path <- paste0(project_path, "results/")
results_GEMX_COPYKAT_path <- paste0(results_path, "CellAnnotation/CopyKAT/")

# Import plot functions & shared utilities (read_copykat_outputs() lives in utils.R)
source(paste0(wd, "CL_plots.R"))
source(paste0(wd, "utils.R"))
setwd(results_GEMX_COPYKAT_path)

# Load annotated data (output of 3f_stroma_annotate.R)
dwAnnotated <- readRDS(paste0(results_path, "CellAnnotation/notumor_annotated_data.rds"))
cat("\n Data loaded \n")

# Set reference cell types (immune cells — copy-number normal)
ref_groups <- c(
  "TCell_naive", "TCell_cyto", "TCell_ex"
)

dwAnnotated_joined <- JoinLayers(dwAnnotated)
ref_cells <- colnames(dwAnnotated_joined)[dwAnnotated_joined$celltype %in% ref_groups]
cat(paste0("\n Reference cells (pooled across all samples): ", length(ref_cells), " \n"))

samples <- sort(unique(as.character(dwAnnotated_joined$orig.ident)))
cat(paste0("\n Running copyKAT per sample (", length(samples), "): ", paste(samples, collapse = ", "), "\n"))

# Run copyKAT per sample
for (s in samples) {

  cat(paste0("\n\n ==== Sample: ", s, " ==== \n"))
  results_sample_path <- paste0(results_GEMX_COPYKAT_path, s, "/")
  dir.create(results_sample_path, recursive = TRUE, showWarnings = FALSE)
  setwd(results_sample_path)  # copyKAT writes output files to the working directory

  obs_cells_s <- colnames(dwAnnotated_joined)[!(dwAnnotated_joined$celltype %in% ref_groups) &
                                                 dwAnnotated_joined$orig.ident == s]
  cells_this_run <- c(obs_cells_s, ref_cells)

  raw_counts_s <- as.matrix(GetAssayData(dwAnnotated_joined, assay = "RNA_decontX", layer = "counts")[, cells_this_run])

  copykat_result <- tryCatch({
    copykat(
      rawmat = raw_counts_s,
      id.type = "S",
      ngene.chr = 2,
      win.size = 25,
      KS.cut = 0.1,
      sam.name = s,
      distance = "euclidean",
      norm.cell.names = ref_cells,
      output.seg = FALSE,
      plot.genes = FALSE,
      genome = "hg20",
      n.cores = 10
    )
  }, error = function(e) {
    cat(paste0("  [ERROR] copyKAT failed for sample ", s, ": ", conditionMessage(e), "\n"))
    return(NULL)
  })

  cat(paste0(" Sample ", s, " done -> ", results_sample_path, "\n"))
}

setwd(results_GEMX_COPYKAT_path)

# --- Read every successful sample's prediction + CNA files back from disk ---
predictions_list <- list()
cna_list <- list()

for (s in samples) {
  results_sample_path <- paste0(results_GEMX_COPYKAT_path, s)
  out_s <- read_copykat_outputs(results_sample_path, s)

  pred_s <- out_s$prediction
  rownames(pred_s) <- pred_s$cell.names
  predictions_list[[s]] <- pred_s

  cna_list[[s]] <- out_s$cna
}
cat(paste0("\n Read copyKAT outputs back from disk for ", length(predictions_list), " samples \n"))

# Combine predictions across samples
barcodes_by_sample <- lapply(predictions_list, function(pred_s) rownames(pred_s))

ref_cells_by_sample <- setNames(
  as.character(dwAnnotated_joined$orig.ident[ref_cells]),
  ref_cells
)

obs_and_ref_predictions_list <- lapply(names(predictions_list), function(s) {
  pred_s <- predictions_list[[s]]

  obs_s <- pred_s[!(rownames(pred_s) %in% ref_cells), ]

  # Reference cells belonging to this sample only
  ref_cells_this_sample <- names(ref_cells_by_sample)[ref_cells_by_sample == s]
  ref_s <- pred_s[rownames(pred_s) %in% ref_cells_this_sample, ]

  rbind(obs_s, ref_s)
})

all_predictions <- do.call(rbind, obs_and_ref_predictions_list)

rownames(all_predictions) <- sub("^.*?\\.", "", rownames(all_predictions))

cat(paste0("\n Combined predictions: ", nrow(all_predictions), " cells across ", length(obs_and_ref_predictions_list), " samples \n"))

# Attach copyKAT results to dwAnnotated
matched_cells <- intersect(colnames(dwAnnotated), rownames(all_predictions))
dwAnnotated$copykat_prediction <- NA
dwAnnotated$copykat_prediction[matched_cells] <- all_predictions[matched_cells, "copykat.pred"]

# =============================================================================
# Combine the per-sample CNA matrices (genomic bin x cell)
# =============================================================================
common_bins <- Reduce(intersect, lapply(cna_list, function(x) paste0(x$chrom, "_", x$chrompos)))

cna_matrix_list <- lapply(names(cna_list), function(s) {
  cna_s <- cna_list[[s]]
  bin_id <- paste0(cna_s$chrom, "_", cna_s$chrompos)
  cna_s <- cna_s[match(common_bins, bin_id), ]

  cell_cols <- setdiff(colnames(cna_s), c("chrom", "chrompos", "abspos"))
  obs_cols_s <- cell_cols[!(cell_cols %in% ref_cells)]

  # Reference cells belonging to this sample only
  ref_cells_this_sample <- names(ref_cells_by_sample)[ref_cells_by_sample == s]
  ref_cols_s <- cell_cols[cell_cols %in% ref_cells_this_sample]

  as.matrix(cna_s[, c(obs_cols_s, ref_cols_s), drop = FALSE])
})

cna_matrix_all <- do.call(cbind, cna_matrix_list)
rownames(cna_matrix_all) <- paste0("chr", common_bins)
cat(paste0("\n Combined CNA matrix: ", nrow(cna_matrix_all), " bins x ", ncol(cna_matrix_all), " cells \n"))

# Pad out to every cell in each object (NA where copyKAT wasn't run)
attach_copykat_cna_assay <- function(object, cna_matrix) {
  cna_full <- matrix(NA_real_, nrow = nrow(cna_matrix), ncol = ncol(object),
                     dimnames = list(rownames(cna_matrix), colnames(object)))
  matched <- intersect(colnames(object), colnames(cna_matrix))
  cna_full[, matched] <- cna_matrix[, matched]
  object[["copykat_cna"]] <- CreateAssayObject(data = cna_full)
  object
}

dwAnnotated <- attach_copykat_cna_assay(dwAnnotated, cna_matrix_all)
cat("\n copyKAT CNA matrix attached as assay 'copykat_cna' \n")

# =============================================================================
# Chromosomal instability: CIN_score (SD of CNA profile), PGA_score (% bins altered)
# =============================================================================
cna_assay_data <- GetAssayData(dwAnnotated, assay = "copykat_cna", layer = "data")
cin_scores <- apply(cna_assay_data, 2, sd, na.rm = TRUE)
pga_scores <- apply(cna_assay_data, 2, function(x) mean(abs(x) > 0.1, na.rm = TRUE) * 100)

attach_cin_scores <- function(object) {
  matched <- intersect(colnames(object), names(cin_scores))
  object$CIN_score <- NA_real_
  object$CIN_score[matched] <- cin_scores[matched]
  object$PGA_score <- NA_real_
  object$PGA_score[matched] <- pga_scores[matched]
  object
}
dwAnnotated <- attach_cin_scores(dwAnnotated)

cat("\n Chromosomal instability score computed (CIN_score = SD of CNA profile, PGA_score = % bins altered) \n")
print(summary(cin_scores))

p_cin <- ggplot(dwAnnotated@meta.data[!is.na(dwAnnotated$CIN_score), ],
                aes(x = copykat_prediction, y = CIN_score, fill = copykat_prediction)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.3) +
  theme_bw() + theme(panel.grid = element_blank()) +
  labs(title = "Chromosomal instability (CIN) score by CopyKAT prediction",
       x = "CopyKAT prediction", y = "CIN score (SD of CNA profile)")
ggsave(paste0(results_GEMX_COPYKAT_path, "Boxplot_CINscore_byCopyKATPrediction.png"), p_cin,
       width = 7, height = 5, dpi = 300, bg = "white")
cat("\n CIN score boxplot saved -> Boxplot_CINscore_byCopyKATPrediction.png \n")

# --- CIN/PGA score by celltype and by sample ---
celltype_lineage_colors <- dwAnnotated@meta.data %>%
  dplyr::distinct(celltype, lineage) %>%
  dplyr::mutate(celltype = as.character(celltype), lineage = as.character(lineage)) %>%
  { setNames(lineage_colors[.$lineage], .$celltype) }

sample_subtype_colors <- dwAnnotated@meta.data %>%
  dplyr::distinct(orig.ident, Subtype) %>%
  dplyr::mutate(orig.ident = as.character(orig.ident), Subtype = as.character(Subtype)) %>%
  { setNames(pam50_subtype_colors[.$Subtype], .$orig.ident) }

plot_score_by_group <- function(object, score_col, group_col, results_path, filename, title, y_lab,
                                facet_col = NULL, colors = NULL) {
  meta <- object@meta.data[!is.na(object@meta.data[[score_col]]), ]

  p <- ggplot(meta, aes(x = .data[[group_col]], y = .data[[score_col]], fill = .data[[group_col]])) +
    geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.3) +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1),
         legend.position = "none") +
    labs(title = title, x = group_col, y = y_lab)

  if (!is.null(colors)) p <- p + scale_fill_manual(values = colors)
  if (!is.null(facet_col)) p <- p + facet_wrap(vars(.data[[facet_col]]), scales = "free_x")

  ggsave(paste0(results_path, filename), p, width = 10, height = 6, dpi = 300, bg = "white")
}

plot_score_by_group(dwAnnotated, "CIN_score", "celltype", results_GEMX_COPYKAT_path,
                    "Boxplot_CINscore_byCelltype.png",
                    "Chromosomal instability (CIN) score by celltype", "CIN score (SD of CNA profile)",
                    facet_col = "lineage", colors = celltype_lineage_colors)
plot_score_by_group(dwAnnotated, "CIN_score", "orig.ident", results_GEMX_COPYKAT_path,
                    "Boxplot_CINscore_bySample.png",
                    "Chromosomal instability (CIN) score by sample", "CIN score (SD of CNA profile)",
                    facet_col = "Subtype", colors = sample_subtype_colors)
plot_score_by_group(dwAnnotated, "PGA_score", "celltype", results_GEMX_COPYKAT_path,
                    "Boxplot_PGAscore_byCelltype.png",
                    "Percent genome altered (PGA) by celltype", "PGA score (% bins altered)",
                    facet_col = "lineage", colors = celltype_lineage_colors)
plot_score_by_group(dwAnnotated, "PGA_score", "orig.ident", results_GEMX_COPYKAT_path,
                    "Boxplot_PGAscore_bySample.png",
                    "Percent genome altered (PGA) by sample", "PGA score (% bins altered)",
                    facet_col = "Subtype", colors = sample_subtype_colors)
cat("\n CIN/PGA score boxplots by celltype (facet lineage) and by sample (facet Subtype) saved \n")

# =============================================================================
# UMAP built directly from the CNA score matrix
# =============================================================================
run_cna_umap <- function(object, npcs = 20, dims = 1:20) {
  cna_data <- GetAssayData(object, assay = "copykat_cna", layer = "data")
  cells_with_cna <- colnames(cna_data)[colSums(!is.na(cna_data)) > 0]

  object_cna <- subset(object, cells = cells_with_cna)
  DefaultAssay(object_cna) <- "copykat_cna"
  VariableFeatures(object_cna) <- rownames(object_cna[["copykat_cna"]])
  object_cna <- ScaleData(object_cna, features = VariableFeatures(object_cna))
  object_cna <- RunPCA(object_cna, features = VariableFeatures(object_cna), npcs = npcs,
                       reduction.name = "pca_cna", reduction.key = "PCACNA_")
  object_cna <- RunUMAP(object_cna, reduction = "pca_cna", dims = dims,
                        reduction.name = "umap_cna", reduction.key = "UMAPCNA_")

  pad_reduction <- function(reduc_obj, all_cells) {
    emb <- reduc_obj@cell.embeddings
    emb_full <- matrix(NA_real_, nrow = length(all_cells), ncol = ncol(emb),
                       dimnames = list(all_cells, colnames(emb)))
    emb_full[rownames(emb), ] <- emb
    CreateDimReducObject(embeddings = emb_full, key = reduc_obj@key, assay = "copykat_cna")
  }

  object[["pca_cna"]] <- pad_reduction(object_cna[["pca_cna"]], colnames(object))
  object[["umap_cna"]] <- pad_reduction(object_cna[["umap_cna"]], colnames(object))
  return(object)
}

dwAnnotated <- run_cna_umap(dwAnnotated)
cat("\n CNA-based UMAP computed -> reduction 'umap_cna' (cells without a CopyKAT profile are NA) \n")

plot_dimplot(dwAnnotated, reduction = "umap_cna", group_by = "copykat_prediction",
            results_path = results_GEMX_COPYKAT_path, filename = "DimPlot_UMAP_CNA_CopyKATPrediction.png")
plot_dimplot(dwAnnotated, reduction = "umap_cna", group_by = "celltype",
            results_path = results_GEMX_COPYKAT_path, filename = "DimPlot_UMAP_CNA_Celltype.png")
plot_dimplot(dwAnnotated, reduction = "umap_cna", group_by = "lineage",
            results_path = results_GEMX_COPYKAT_path, filename = "DimPlot_UMAP_CNA_Lineage.png")
plot_dimplot(dwAnnotated, reduction = "umap_cna", group_by = "orig.ident",
            results_path = results_GEMX_COPYKAT_path, filename = "DimPlot_UMAP_CNA_Sample.png")
plot_dimplot(dwAnnotated, reduction = "umap_cna", group_by = "Subtype",
            results_path = results_GEMX_COPYKAT_path, filename = "DimPlot_UMAP_CNA_ClinicalSubtype.png")

cna_values <- as.vector(cna_assay_data)
cna_values <- cna_values[!is.na(cna_values)]

cat("\n CNA value summary (all bins x cells):\n")
print(summary(cna_values))
cat(paste0("  sd = ", round(sd(cna_values), 4), "\n"))

p_cna_dist <- ggplot(data.frame(value = cna_values), aes(x = value)) +
  geom_histogram(bins = 100, fill = "steelblue", color = NA) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "firebrick") +
  theme_bw() + theme(panel.grid = element_blank()) +
  labs(title = "Distribution of copyKAT CNA values (all bins x cells)",
       x = "CNA value", y = "Count")
ggsave(paste0(results_GEMX_COPYKAT_path, "Histogram_CNA_values.png"), p_cna_dist,
       width = 7, height = 5, dpi = 300, bg = "white")
cat("\n CNA value distribution histogram saved -> Histogram_CNA_values.png \n")

# --- DimPlots of copyKAT prediction ---
plot_dimplot(dwAnnotated, reduction = "umap_decontX", group_by = "copykat_prediction",
             results_path = results_GEMX_COPYKAT_path, filename = "DimPlot_UMAP_CopyKATPrediction.png")

# --- Plot Ploidy Distribution ---
plot_composition(dwAnnotated, cluster_col = "orig.ident", group_by = "copykat_prediction",
                          filename = "BarPlot_CopyKATPrediction_bySample.png", results_path = results_GEMX_COPYKAT_path
)

plot_composition(dwAnnotated, cluster_col = "decontX_clusters", group_by = "copykat_prediction",
                         filename = "BarPlot_CopyKATPrediction_byCluster.png", results_path = results_GEMX_COPYKAT_path
)

plot_composition(dwAnnotated, cluster_col = "celltype", group_by = "copykat_prediction",
                         filename = "BarPlot_CopyKATPrediction_byCelltype.png", results_path = results_GEMX_COPYKAT_path)

# Export annotated data with copyKAT predictions
saveRDS(dwAnnotated, paste0(results_GEMX_COPYKAT_path, "copykat_annotated_data.rds"))



cat(paste("\n ---- FINISHED COPYKAT ANNOTATION ----
    Run 3h_tumor_scoring.R next.
    Generated files:
        · copykat_annotated_data.rds
        · (sample)/ subfolders with copyKAT per-sample output
    Generated plots:
        · (sample)/(sample)_copykat_heatmap.jpeg
        · BarPlot_CopyKATPrediction_(groupedby).png
        · Boxplot_CopyKAT_CNA_byCluster.png
        · Boxplot_CINscore_byCopyKATPrediction.png (CIN_score/PGA_score in metadata)
        · Boxplot_(CIN|PGA)score_by(Celltype|Sample).png
        · DimPlot_UMAP_CNA_(CopyKATPrediction|Celltype|Lineage|Sample|ClinicalSubtype).png (reduction 'umap_cna')
        "))
