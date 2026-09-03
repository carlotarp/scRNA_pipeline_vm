##
##  Single Cell Analysis Step 3e: CopyKAT CNV annotation
##  Runs AFTER 3a_lineage.R — takes fully_annotated_data.rds as input.
##  Runs copyKAT once per sample (fixed reference cells pooled across all
##  samples) and attaches ploidy predictions + CNA event counts to dwAnnotated.
##

# Import libraries
library(Seurat)
library(copykat)
library(dplyr)
library(ggplot2)

# Set paths
project_path <- "/home/usuario/PROJECTS/260724_victor_scRNA/"
wd <- paste0(project_path, "codes/scRNA_pipeline/")
setwd(wd)
results_path <- paste0(project_path, "results/")
results_GEMX_COPYKAT_path <- paste0(results_path, "GEMX/DecontX/CellAnnotation/CopyKAT/")
setwd(results_GEMX_COPYKAT_path)

# Load annotated data
dwAnnotated <- readRDS(paste0(results_path, "GEMX/DecontX/CellAnnotation/fully_annotated_data.rds"))
cat("\n Data loaded \n")

dwAnnotated <- readRDS(paste0(results_GEMX_COPYKAT_path, "copykat_annotated_data.rds"))

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
predictions_list <- list()
cna_list <- list()

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
      win.size00000000000000000 = 25,
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

  if (is.null(copykat_result)) next

  pred_s <- copykat_result$prediction
  rownames(pred_s) <- pred_s$cell.names
  predictions_list[[s]] <- pred_s

  cna_list[[s]] <- copykat_result$CNAmat

  generated_files <- list.files(results_sample_path, full.names = TRUE)
  files_to_delete <- generated_files[!grepl("\\.jpeg$", generated_files, ignore.case = TRUE)]

  cat(paste0(" Sample ", s, " done -> ", results_sample_path, "\n"))
}

setwd(results_GEMX_COPYKAT_path)

# Combine predictions across samples
barcodes_by_sample <- lapply(predictions_list, function(pred_s) rownames(pred_s))

obs_predictions_list <- lapply(predictions_list, function(pred_s) {
  pred_s[!(rownames(pred_s) %in% ref_cells), ]
})

all_predictions <- do.call(rbind, obs_predictions_list)
rownames(all_predictions) <- sub("^.*?\\.", "", rownames(all_predictions))

cat(paste0("\n Combined predictions: ", nrow(all_predictions), " cells across ", length(obs_predictions_list), " samples \n"))

# Attach copyKAT results to dwAnnotated
matched_cells <- intersect(colnames(dwAnnotated), rownames(all_predictions))
dwAnnotated$copykat_prediction[matched_cells] <- all_predictions[matched_cells, "copykat.pred"]

cat("\n copyKAT classification: \n")
print(table(dwAnnotated$copykat_prediction))
print(table(dwAnnotated$copykat_prediction, dwAnnotated$orig.ident))
print(table(dwAnnotated$copykat_prediction, dwAnnotated$celltype, dwAnnotated$orig.ident))

# Merge per-sample CNA matrices into a single genes × cells matrix
cna_combined <- Reduce(function(x, y) {
  merge(x, y[, !colnames(y) %in% intersect(colnames(x), colnames(y))[-c(1:3)]],
        by = colnames(x)[1:3])
}, cna_list)
cna_values <- cna_combined[, -c(1, 2, 3)]
chrom_vector <- cna_combined[, 1]

count_cna_events <- function(cell_values, chrom_vector, threshold = 0.1) {
  states <- ifelse(cell_values > threshold, 1,
             ifelse(cell_values < -threshold, -1, 0))
  n_events <- 0
  for (chr in unique(chrom_vector)) {
    states_chr <- states[chrom_vector == chr]
    runs <- rle(states_chr)
    # Each contiguous non-neutral run = one CNA event
    n_events <- n_events + sum(runs$lengths[runs$values != 0] > 0)
  }
  return(n_events)
}

n_cna_per_cell <- apply(cna_values, 2, count_cna_events, chrom_vector = chrom_vector)
names(n_cna_per_cell) <- gsub("\\.", "-", names(n_cna_per_cell))
matched_cells <- intersect(colnames(dwAnnotated), names(n_cna_per_cell))
dwAnnotated$copykat_cnas[matched_cells] <- n_cna_per_cell[matched_cells]

saveRDS(dwAnnotated, paste0(results_GEMX_COPYKAT_path, "copykat_annotated_data.rds"))
cat("\n copyKAT results attached to dwAnnotated and saved \n")

cat(paste("\n ---- FINISHED COPYKAT ANNOTATION ----
    Generated files:
        · copykat_annotated_data.rds
        · (sample)/ subfolders with copyKAT per-sample output
    Generated plots:
        · (sample)/(sample)_copykat_heatmap.jpeg
        "))
