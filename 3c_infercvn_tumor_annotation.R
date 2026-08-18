##
##  Single Cell Analysis Step 3c: inferCNV
##

library(dplyr)
library(Seurat)
library(SeuratObject)
library(infercnv)

# Set Paths
project_path <- "/home/usuario/PROJECTS/260724_victor_scRNA/"
wd <- paste0(project_path, "codes/scRNA_pipeline/")
setwd(wd)
results_path <- paste0(project_path, "results/")
results_GEMX_CNV_path <- paste0(results_path, "GEMX/CellAnnotation/7500/InferCNV/")

source(paste0(wd, "CA_plots.R"))

# Load Fully Annotated Data
dwAnnotated <- readRDS(paste0(results_path, "GEMX/CellAnnotation/7500/notumor_clean_annotated_data.rds"))
cat("\n Fully annotated data loaded \n")

cat("\n Available celltype values (use these to fill in ref_groups below): \n")
print(table(dwAnnotated$celltype))


# =================================================================
# MANUAL REFERENCE - edit this list to match your $celltype values exactly
# =================================================================
ref_groups <- c(
  "BCell",
  "Tcell_naive",
  "Tcell_cyto",
  "Tcell_ex",
  "PlasmaBlast",
  "pDC",
  "actDC",
  "Mast"
)


# inferCNV needs raw counts
dwAnnotated <- JoinLayers(dwAnnotated)
counts_matrix <- GetAssayData(dwAnnotated, assay = "RNA", layer = "counts")

# Reference cells - fixed, same set used in every sample's run
ref_cells <- colnames(dwAnnotated)[dwAnnotated$celltype %in% ref_groups]

samples <- sort(unique(as.character(dwAnnotated$orig.ident)))
cat(paste0("\n Running inferCNV per sample (", length(samples), "): ", paste(samples, collapse = ", "), "\n"))

gene_order_file <- paste0(project_path, "data/gencode_v19_gene_pos.txt")

# =================================================================
# Run inferCNV once per sample: ALL reference cells (pooled across every
# sample) vs every NON-reference cell belonging to THIS sample specifically
# (Tumoral, but also Fibrocyte/Proliferative/anything else not in ref_groups
# if that sample has them).
# =================================================================
run_infercnv <- function(s){

  cat(paste0("\n\n ==== Sample: ", s, " ==== \n"))
  results_sample_path <- paste0(results_GEMX_CNV_path, s, "/")
  dir.create(results_sample_path, recursive = TRUE, showWarnings = FALSE)

#  ref_cells_s <- colnames(dwAnnotated)[dwAnnotated$celltype %in% ref_groups & dwAnnotated$orig.ident == s]
  obs_cells_s <- colnames(dwAnnotated)[!(dwAnnotated$celltype %in% ref_groups) & dwAnnotated$orig.ident == s]
  cells_this_run <- c(obs_cells_s, ref_cells)

  cell_annotations_s <- data.frame(
    cell = cells_this_run,
    group = as.character(dwAnnotated$celltype[cells_this_run])
  )

  annotations_file_s <- paste0(results_sample_path, "inferCNV_cell_annotations.txt")
  write.table(cell_annotations_s, annotations_file_s,
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

  infercnv_obj_s <- CreateInfercnvObject(
    raw_counts_matrix = counts_matrix[, cells_this_run],
    annotations_file = annotations_file_s,
    delim = "\t",
    gene_order_file = gene_order_file,
    ref_group_names = ref_groups,
  )

  infercnv_obj_s_run <- infercnv::run(
    infercnv_obj_s,
    cutoff = 0.1,
    out_dir = results_sample_path,
    cluster_by_groups = TRUE,
    denoise = TRUE,
    HMM = TRUE,
    analysis_mode = "samples",
    save_rds = F,
    no_plot = T,
    no_prelim_plot = T,
    num_threads = 10
  )

  cat(paste0("\n Sample ", s, " done -> ", results_sample_path, "\n"))
}
process_infercnv <- function(s) {
  cnv_file <- paste0(results_GEMX_CNV_path, s, "/HMM_CNV_predictions.HMMi6.hmm_mode-samples.Pnorm_0.5.pred_cnv_regions.dat")
  group_file <- paste0(results_GEMX_CNV_path, s, "/17_HMM_predHMMi6.hmm_mode-samples.cell_groupings")

  cnv <- read.table(cnv_file, header = TRUE, sep = "\t")
  groups <- read.table(group_file, header = TRUE, sep = "\t")

  group_ploidy <- cnv %>%
    group_by(cell_group_name) %>%
    summarise(any_altered = any(state != 3), .groups = "drop") %>%
    mutate(ploidy_status = ifelse(any_altered, "aneuploid", "diploid")) %>%
    select(cell_group_name, ploidy_status)

  cell_ploidy <- groups %>%
    left_join(group_ploidy, by = "cell_group_name") %>%
    mutate(ploidy_status = ifelse(is.na(ploidy_status), "diploid", ploidy_status)) %>%
    select(cell, ploidy_status) %>%
    rename(cell_id = cell)

  output_file <- paste0(results_GEMX_CNV_path, s, "_cell_aneuploid_diploid_annotations.txt")
  write.table(cell_ploidy, output_file, sep = "\t", quote = FALSE, row.names = FALSE)

  cat("Cell ploidy file gen", s, "\n")
}


annotate_ploidy <- function (s, object) {

  ploidy_file <- paste0(results_GEMX_CNV_path, s, "_cell_aneuploid_diploid_annotations.txt")
  ploidy_df <- read.table(ploidy_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

  # Match por barcode - solo actualiza las células de ESTA muestra
  matched_idx <- match(ploidy_df$cell_id, colnames(object))
  valid <- !is.na(matched_idx)

  object$ploidy[matched_idx[valid]] <- ploidy_df$ploidy[valid]

  cat(s, ":", sum(valid), "células actualizadas (de", nrow(ploidy_df), "en el archivo)\n")

  return(object)
}

run_infercnv("SC13")
process_infercnv("SC13")
if (!"ploidy" %in% colnames(dwAnnotated@meta.data)) {
  dwAnnotated$ploidy <- NA
}
dwAnnotated <- annotate_ploidy("SC13", dwAnnotated)
dwAnnotated$ploidy <- factor(dwAnnotated$ploidy, levels = c("diploid", "aneuploid"))

# UMAP global coloreado por ploidy
plot_dimplot(dwAnnotated, reduction = "umap", group_by = "ploidy",
             results_path = results_path, filename = "DimPlot_ploidy.png")

# Cross-tab: ¿ploidy coincide con lo esperado por celltype?
# (referencia = diploid, Tumoral = aneuploid; Fibrocyte/Proliferative = lo interesante)
ploidy_by_celltype <- table(dwAnnotated$celltype[dwAnnotated$orig.ident == "SC13"], dwAnnotated$ploidy[dwAnnotated$orig.ident == "SC13"], useNA = "ifany")
print(ploidy_by_celltype)
write.csv(as.data.frame.matrix(ploidy_by_celltype), paste0(results_GEMX_CNV_path, "ploidy_by_celltype.csv"))

# % aneuploide dentro de cada celltype - la tabla más directa para leer
ploidy_pct <- prop.table(ploidy_by_celltype, margin = 1) * 100
print(round(ploidy_pct, 1))