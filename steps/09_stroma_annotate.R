#!/usr/bin/env Rscript
##
## 09_stroma_annotate.R - apply the stroma dictionary and transfer the labels to the full object
## Usage : Rscript pipeline/steps/09_stroma_annotate.R
## Reads : objects/leuko_annotated.rds ; objects/stroma.rds ; annotation/stroma.tsv
## Writes: objects/stroma_annotated.rds (subset with celltype) ; objects/notumor_annotated.rds (full object) ;
##         09_stroma_annotate/ DimPlot_UMAP_Stroma_Annotated.png, DimPlot_UMAP_Stroma.png ;
##         annotation/excluded_cells.csv (DROP clusters)
## Checkpoint: annotation dictionary "stroma"
##
file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
PIPELINE_DIR <- normalizePath(file.path(dirname(file_arg), ".."))
source(file.path(PIPELINE_DIR, "R/setup.R"))
out_dir <- step_dir("09_stroma_annotate")

dwStroma <- read_run_object("stroma")

# Dictionary: stromal sub-cluster -> cell type (DROP removes the cluster)
stroma_dict <- read_annotation_dictionary("stroma", dwStroma$seurat_clusters)
res <- apply_annotation(dwStroma, stroma_dict, cluster_col = "seurat_clusters", label_col = "celltype")
log_excluded_cells(dwStroma, res$dropped, step = "stroma", cluster_col = "seurat_clusters")
dropped <- res$dropped
dwStroma <- res$object
rm(res); invisible(gc())
dwStroma$celltype <- factor(dwStroma$celltype)

plot_dimplot(dwStroma, reduction = "umap", group_by = "celltype", label = TRUE,
             results_path = out_dir, filename = "DimPlot_UMAP_Stroma_Annotated.png")
saveRDS(dwStroma, run_object("stroma_annotated"))

# Transfer the stroma labels to the full object (other cells keep their celltype; DROP cells are removed)
dwAnnotated <- read_run_object("leuko_annotated")
dwAnnotated <- transfer_labels(dwAnnotated, dwStroma, dropped = dropped, label_col = "celltype")
dwAnnotated$celltype <- factor(dwAnnotated$celltype)

plot_dimplot(dwAnnotated, reduction = "umap", group_by = "celltype", label = TRUE,
             results_path = out_dir, filename = "DimPlot_UMAP_Stroma.png")

saveRDS(dwAnnotated, run_object("notumor_annotated"))

cat("\n ---- FINISHED 09_stroma_annotate ----\n",
    "Wrote", run_object("stroma_annotated"), "and", run_object("notumor_annotated"), "\n",
    "Run 10_tumour_scoring.R next\n")
