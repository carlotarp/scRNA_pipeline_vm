#!/usr/bin/env Rscript
##
## 04_lineage_markers.R - evidence to label each global cluster with its lineage (Tumor / Leukocytes / Stromal)
## Usage : Rscript pipeline/steps/04_lineage_markers.R
## Reads : objects/clustered.rds
## Writes: annotation/lineage.tsv (template) ;
##         04_lineage_markers/ DotPlot_<block>.png, cluster_<cl>_FindMarkers[_filtered].csv,
##         all_clusters_FindMarkers[_filtered].csv, top_expressed_genes.csv
## Checkpoint: none (fill annotation/lineage.tsv before 05_lineage_annotate.R;
##             labels Tumor, Leukocytes, Stromal or DROP)
##
file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
PIPELINE_DIR <- normalizePath(file.path(dirname(file_arg), ".."))
source(file.path(PIPELINE_DIR, "R/setup.R"))
out_dir <- step_dir("04_lineage_markers")

cluster_col <- "global_clusters"   # set by 03_cluster.R (= seurat_clusters at the chosen resolution)

dwClustered <- read_run_object("clustered")
Idents(dwClustered) <- cluster_col

# Lineage dot plots, one per panel block (leukocytes, stromal, vascular, epithelial)
for (block in names(load_marker_panel("lineage"))) {
  plot_marker_dotplot(dwClustered, marker_groups = load_marker_panel("lineage", block, object = dwClustered),
                      results_path = out_dir, filename = paste0("DotPlot_", block, ".png"),
                      group_by = cluster_col, title = paste0("Lineage markers (", block, ") by cluster"))
}

# Markers of every cluster (the filtered list is the one to annotate from)
markers_filt <- find_markers_for_clusters(dwClustered, levels(Idents(dwClustered)), out_dir)

# Top expressed genes per cluster (ambient RNA / doublet evidence)
top_expressed_genes(dwClustered, group_by = cluster_col, results_path = out_dir)

# Dictionary template: cluster -> lineage
write_annotation_template(dwClustered, "lineage", cluster_col = cluster_col, markers_df = markers_filt)

cat("\n ---- FINISHED 04_lineage_markers ----\n",
    "Fill", annotation_file("lineage"), "(label: Tumor, Leukocytes, Stromal or DROP; evidence required),\n",
    "then run 05_lineage_annotate.R\n")
