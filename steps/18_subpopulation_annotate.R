#!/usr/bin/env Rscript
##
## 18_subpopulation_annotate.R - apply the hand dictionary of a sub-population study (step 17) to the full object
## Usage : Rscript pipeline/steps/18_subpopulation_annotate.R <study>
## Reads : objects/<input_object>.rds, objects/subpop_<study>.rds (zoom of step 17)
## Writes: objects/<output_object>.rds ; 18_subpopulation_annotate/<study>/labels.csv.gz (cell, sample, zoom cluster,
##         label before/after) + UMAP of the cell types
## Checkpoint: annotation/<study>.tsv (zoom cluster -> label; DROP removes the cluster, e.g. doublets)
##
## Every cell of the zoom gets a dictionary label, so no cell keeps a zoomed label (subpopulations.<study>.celltypes)
## the dictionary does not use. Lineages are unchanged.
##
file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
PIPELINE_DIR <- normalizePath(file.path(dirname(file_arg), ".."))
source(file.path(PIPELINE_DIR, "R/setup.R"))
study <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(study)) stop("Usage: 18_subpopulation_annotate.R <study>")
key <- function(k) paste0("subpopulations.", study, ".", k)
out_dir <- step_dir("18_subpopulation_annotate", study)

z <- read_run_object(paste0("subpop_", study))
dict <- read_annotation_dictionary(study, z$seurat_clusters)
new <- setNames(unname(dict[as.character(z$seurat_clusters)]), colnames(z))
zm <- z@meta.data[, c("orig.ident", "seurat_clusters")]
dropped <- names(new)[new == "DROP"]
log_excluded_cells(z, dropped, study)
rm(z); invisible(gc())

a <- read_run_object(cfg_get(key("input_object")))
if (!all(names(new) %in% colnames(a))) stop("Zoom '", study, "' has cells absent from ", cfg_get(key("input_object")))
before <- setNames(as.character(a$celltype), colnames(a))
write.csv(data.frame(cell = names(new), sample = zm$orig.ident, zoom_group = zm$seurat_clusters,
                     celltype_before = before[names(new)], celltype = new),
          gzfile(paste0(out_dir, "labels.csv.gz")), row.names = FALSE)
print(table(before = before[names(new)], after = new))

ct <- before; ct[names(new)] <- new
a$celltype <- ct
a <- subset(a, cells = setdiff(colnames(a), dropped))
a$celltype <- factor(as.character(a$celltype)); a$lineage <- factor(as.character(a$lineage))
left <- setdiff(as.character(get_decision(key("celltypes"))), dict)
if (any(as.character(a$celltype) %in% left)) stop("Cells outside the zoom still carry: ", paste(intersect(left, a$celltype), collapse = ", "))
out_name <- cfg_get(key("output_object"))
saveRDS(a, run_object(out_name))
cat("  ->", run_object(out_name), "\n")
print(table(a$lineage)); print(sort(table(as.character(a$celltype))))
plot_dimplot(a, reduction = "umap", group_by = "celltype", results_path = out_dir, filename = "DimPlot_UMAP_celltype.png")
cat("\n18_subpopulation_annotate", study, "done\n")
