#!/usr/bin/env Rscript
##
## 12_tumour_annotate.R - apply the tumour dictionary (DROP removes noise clusters) and transfer labels to the full object
## Usage : Rscript pipeline/steps/12_tumour_annotate.R
## Reads : objects/tumour_scored.rds, objects/notumor_annotated.rds, annotation/tumour.tsv
## Writes: objects/tumour_annotated.rds, objects/annotated.rds ; annotation/excluded_cells.csv (DROP cells) ;
##         12_tumour_annotate/<UMAPs, PAM50/scSubtype proportion per celltype>
## Checkpoint: annotation/tumour.tsv (written by 11_tumour_markers.R, filled by hand)
##
file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
PIPELINE_DIR <- normalizePath(file.path(dirname(file_arg), ".."))
source(file.path(PIPELINE_DIR, "R/setup.R"))
out_dir <- step_dir("12_tumour_annotate")

tumour <- read_run_object("tumour_scored")
dict <- read_annotation_dictionary("tumour", tumour$seurat_clusters)
annotated <- read_run_object("notumor_annotated")
if (any(annotated$lineage == "Noise")) stop("notumor_annotated still has lineage 'Noise' cells: earlier annotate steps must DROP them")

# --- Tumour labels; DROP clusters removed from both objects ---
res <- apply_annotation(tumour, dict)
log_excluded_cells(tumour, res$dropped, step = "tumour")
tumour <- res$object
tumour$celltype <- factor(tumour$celltype)

annotated <- transfer_labels(annotated, tumour, dropped = res$dropped)
annotated$celltype <- factor(annotated$celltype)
annotated$lineage <- factor(as.character(annotated$lineage))

plot_dimplot(tumour, reduction = "umap", group_by = "seurat_clusters", label = TRUE,
             results_path = out_dir, filename = "DimPlot_UMAP_TumorCleaned.png")
plot_dimplot(tumour, reduction = "umap", group_by = "celltype", label = TRUE,
             results_path = out_dir, filename = "DimPlot_UMAP_TumorCleaned_byCelltype_Tumoral.png")
plot_dimplot(annotated, reduction = "umap", group_by = "celltype",
             results_path = out_dir, filename = "DimPlot_UMAP_Cleaned_byCellType.png")

saveRDS(tumour, run_object("tumour_annotated"))
saveRDS(annotated, run_object("annotated"))
rm(annotated); invisible(gc())

# --- PAM50 / scSubtype call proportions within each tumour celltype ---
pam50_prop_df <- compute_celltype_proportions(tumour, sample_id_col = "celltype",
                                              celltype_col = "PAM50_predicted", clinical_cols = "celltype")
plot_celltype_proportions_bar(pam50_prop_df, x_var = "PAM50_predicted", results_path = out_dir,
                              filename = "Barplot_PAM50Proportion_byCelltype.png",
                              celltype_col = "celltype", colors = pam50_subtype_colors)

scsubtype_prop_df <- compute_celltype_proportions(tumour, sample_id_col = "celltype",
                                                  celltype_col = "scSubtype_call", clinical_cols = "celltype")
plot_celltype_proportions_bar(scsubtype_prop_df, x_var = "scSubtype_call", results_path = out_dir,
                              filename = "Barplot_scSubtypeProportion_byCelltype.png",
                              celltype_col = "celltype", colors = scsubtype_colors)

cat("\n ---- FINISHED 12_tumour_annotate: objects/tumour_annotated.rds, objects/annotated.rds ----\n")
