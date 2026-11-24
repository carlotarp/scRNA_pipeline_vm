#!/usr/bin/env Rscript
##
## 16_tumour_copykat_annotate.R - tumour re-annotation by hand with the CopyKAT sub-zooms (patches the saved objects)
## Usage : Rscript pipeline/steps/16_tumour_copykat_annotate.R
## Reads : objects/annotated.rds, objects/tumour_annotated.rds, objects/tumour_lp.rds, objects/tumour_lowrna.rds
## Writes: objects/annotated_copykat.rds, objects/tumour_annotated_copykat.rds ;
##         16_tumour_copykat_annotate/final_labels.csv.gz (cell, sample, celltype before/after, lineage after) + UMAPs
## Checkpoint: annotation/tumour_lp.tsv, annotation/tumour_lowrna.tsv (cluster -> label, from step 15 evidence);
##             tumour_copykat.normal_epithelium_celltypes (labels moved to lineage Normal_Epithelium)
##
## Dictionary labels: a cell type; KEEP = the cells keep their current label; DROP = cells removed from both
## objects (logged in annotation/excluded_cells.csv). The pre-CopyKAT objects of step 12 are not modified.
##
file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
PIPELINE_DIR <- normalizePath(file.path(dirname(file_arg), ".."))
source(file.path(PIPELINE_DIR, "R/setup.R"))
out_dir <- step_dir("16_tumour_copykat_annotate")

# ---- new label per zoomed cell (KEEP -> no entry) ----
new <- character(0); dropped <- character(0)
for (key in c("tumour_lp", "tumour_lowrna")) {
  z <- read_run_object(key)
  dict <- read_annotation_dictionary(key, z$seurat_clusters)
  lab <- setNames(unname(dict[as.character(z$seurat_clusters)]), colnames(z))
  d <- names(lab)[lab == "DROP"]
  log_excluded_cells(z, d, key)
  dropped <- c(dropped, d)
  lab <- lab[!lab %in% c("KEEP", "DROP")]
  if (any(names(lab) %in% names(new))) stop("Cells relabelled by more than one zoom: ", key)
  new <- c(new, lab)
  cat(sprintf("  %s: %d cells relabelled, %d dropped\n", key, length(lab), length(d)))
  rm(z); invisible(gc())
}
normal_types <- as.character(unlist(get_decision("tumour_copykat.normal_epithelium_celltypes")))

# ---- full object ----
a <- read_run_object("annotated")
a <- subset(a, cells = setdiff(colnames(a), dropped))
before <- setNames(as.character(a$celltype), colnames(a))
ct <- before; hit <- intersect(names(new), colnames(a)); ct[hit] <- new[hit]
lin <- as.character(a$lineage); lin[ct %in% normal_types] <- "Normal_Epithelium"
a$celltype <- factor(ct); a$lineage <- factor(lin)
saveRDS(a, run_object("annotated_copykat"))
write.csv(data.frame(cell = colnames(a), sample = a$orig.ident, celltype_before = before, celltype = ct, lineage = lin),
          gzfile(paste0(out_dir, "final_labels.csv.gz")), row.names = FALSE)
print(table(before = before[before != ct], after = ct[before != ct]))
print(table(a$lineage)); print(sort(table(a$celltype[a$lineage %in% c("Tumor", "Normal_Epithelium")])))
plot_dimplot(a, reduction = "umap", group_by = "celltype", results_path = out_dir, filename = "DimPlot_UMAP_final_celltype.png")
rm(a); invisible(gc())

# ---- tumour object (celltype only, as before) ----
t <- read_run_object("tumour_annotated")
t <- subset(t, cells = setdiff(colnames(t), dropped))
tc <- setNames(as.character(t$celltype), colnames(t)); hit <- intersect(names(new), colnames(t)); tc[hit] <- new[hit]
t$celltype <- factor(tc)
saveRDS(t, run_object("tumour_annotated_copykat"))
plot_dimplot(t, reduction = "umap", group_by = "celltype", label = TRUE, results_path = out_dir,
             filename = "DimPlot_UMAP_tumour_final_celltype.png")
cat("\n16_tumour_copykat_annotate done. Next: 17_subpopulation_evidence.R caf <mode>\n")
