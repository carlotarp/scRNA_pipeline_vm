#!/usr/bin/env Rscript
##
## 07_leukocytes_annotate.R - apply the leukocyte dictionary and transfer the labels to the full object
## Usage : Rscript pipeline/steps/07_leukocytes_annotate.R
## Reads : objects/lineage_annotated.rds ; objects/leukocytes.rds ; annotation/leukocytes.tsv
## Writes: objects/leukocytes_annotated.rds (subset with celltype) ; objects/leuko_annotated.rds (full object) ;
##         07_leukocytes_annotate/ DimPlot_UMAP_Leuko_Annotated.png, DimPlot_UMAP_Leuko.png ;
##         annotation/excluded_cells.csv (DROP clusters)
## Checkpoint: annotation dictionary "leukocytes"
##
file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
PIPELINE_DIR <- normalizePath(file.path(dirname(file_arg), ".."))
source(file.path(PIPELINE_DIR, "R/setup.R"))
out_dir <- step_dir("07_leukocytes_annotate")

dwLeukocytes <- read_run_object("leukocytes")

# Dictionary: leukocyte sub-cluster -> cell type (DROP removes the cluster)
leuko_dict <- read_annotation_dictionary("leukocytes", dwLeukocytes$seurat_clusters)
res <- apply_annotation(dwLeukocytes, leuko_dict, cluster_col = "seurat_clusters", label_col = "celltype")
log_excluded_cells(dwLeukocytes, res$dropped, step = "leukocytes", cluster_col = "seurat_clusters")
dropped <- res$dropped
dwLeukocytes <- res$object
rm(res); invisible(gc())
dwLeukocytes$celltype <- factor(dwLeukocytes$celltype)

plot_dimplot(dwLeukocytes, reduction = "umap", group_by = "celltype",
             results_path = out_dir, filename = "DimPlot_UMAP_Leuko_Annotated.png")
saveRDS(dwLeukocytes, run_object("leukocytes_annotated"))

# Transfer the leukocyte labels to the full object (other cells keep their celltype; DROP cells are removed)
dwAnnotated <- read_run_object("lineage_annotated")
dwAnnotated <- transfer_labels(dwAnnotated, dwLeukocytes, dropped = dropped, label_col = "celltype")
dwAnnotated$celltype <- factor(dwAnnotated$celltype)

plot_dimplot(dwAnnotated, reduction = "umap", group_by = "celltype", label = TRUE,
             results_path = out_dir, filename = "DimPlot_UMAP_Leuko.png")

saveRDS(dwAnnotated, run_object("leuko_annotated"))

cat("\n ---- FINISHED 07_leukocytes_annotate ----\n",
    "Wrote", run_object("leukocytes_annotated"), "and", run_object("leuko_annotated"), "\n",
    "Run 08_stroma_markers.R next\n")
