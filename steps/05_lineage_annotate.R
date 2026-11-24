#!/usr/bin/env Rscript
##
## 05_lineage_annotate.R - apply the lineage dictionary (global cluster -> Tumor / Leukocytes / Stromal)
## Usage : Rscript pipeline/steps/05_lineage_annotate.R
## Reads : objects/clustered.rds ; annotation/lineage.tsv
## Writes: objects/lineage_annotated.rds (columns lineage, celltype = lineage) ;
##         05_lineage_annotate/DimPlot_UMAP_Lineage.png ; annotation/excluded_cells.csv (DROP clusters)
## Checkpoint: annotation dictionary "lineage"
##
file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
PIPELINE_DIR <- normalizePath(file.path(dirname(file_arg), ".."))
source(file.path(PIPELINE_DIR, "R/setup.R"))
out_dir <- step_dir("05_lineage_annotate")

cluster_col <- "global_clusters"
allowed_lineages <- c("Tumor", "Leukocytes", "Stromal")   # the lineage sub-clustering steps select on these names

dwClustered <- read_run_object("clustered")

# Dictionary: global cluster -> lineage (DROP removes the cluster)
lineage_dict <- read_annotation_dictionary("lineage", dwClustered[[cluster_col, drop = TRUE]])
bad <- setdiff(unique(lineage_dict), c(allowed_lineages, "DROP"))
if (length(bad)) checkpoint_stop(sprintf("annotation/lineage.tsv: unknown lineage labels: %s (allowed: %s, DROP)",
                                         paste(bad, collapse = ", "), paste(allowed_lineages, collapse = ", ")))

res <- apply_annotation(dwClustered, lineage_dict, cluster_col = cluster_col, label_col = "lineage")
log_excluded_cells(dwClustered, res$dropped, step = "lineage", cluster_col = cluster_col)
dwAnnotated <- res$object
rm(dwClustered, res); invisible(gc())

dwAnnotated$lineage <- factor(dwAnnotated$lineage)

plot_dimplot(dwAnnotated, reduction = "umap", group_by = "lineage",
             results_path = out_dir, filename = "DimPlot_UMAP_Lineage.png")

# celltype starts as the lineage; the leukocyte / stroma / tumour steps refine it
dwAnnotated$celltype <- factor(as.character(dwAnnotated$lineage))

saveRDS(dwAnnotated, run_object("lineage_annotated"))

cat("\n ---- FINISHED 05_lineage_annotate ----\n",
    "Wrote", run_object("lineage_annotated"), "\n",
    "Run 06_leukocytes_markers.R next\n")
