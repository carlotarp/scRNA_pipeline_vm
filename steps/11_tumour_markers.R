#!/usr/bin/env Rscript
##
## 11_tumour_markers.R - evidence to annotate the tumour clusters (markers, PAM50/scSubtype, PROGENy, cell cycle, composition)
## Usage : Rscript pipeline/steps/11_tumour_markers.R
## Reads : objects/tumour_scored.rds, objects/notumor_annotated.rds (celltype composition per sample)
## Writes: 11_tumour_markers/<FindMarkers CSVs, top_expressed_genes.csv, figures, differential score tables> ;
##         annotation/tumour.tsv (template)
## Checkpoint: none here; fill annotation/tumour.tsv (label DROP for noise clusters) before 12_tumour_annotate.R
##
file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
PIPELINE_DIR <- normalizePath(file.path(dirname(file_arg), ".."))
source(file.path(PIPELINE_DIR, "R/setup.R"))
out_dir <- step_dir("11_tumour_markers")

tumour <- read_run_object("tumour_scored")
Idents(tumour) <- "seurat_clusters"

pam50_genes <- load_marker_panel("tumour_pam50")
pam50_score_cols <- paste0("PAM50_", names(pam50_genes))
pam50_fill_labels <- setNames(names(pam50_genes), pam50_score_cols)
scsubtype_score_cols <- names(load_scsubtype_signatures())

# =============================================================================
# 0. Cluster markers + annotation template
# =============================================================================
markers_filt <- find_markers_for_clusters(tumour, levels(tumour$seurat_clusters), out_dir, top_n = 50)
top_expressed_genes(tumour, group_by = "seurat_clusters", results_path = out_dir)
write_annotation_template(tumour, "tumour", markers_df = markers_filt)

# =============================================================================
# 1. PAM50 validation & composition
# =============================================================================
plot_marker_dotplot(tumour, marker_groups = pam50_genes, results_path = out_dir,
                    filename = "Dotplot_PAM50.png", group_by = "seurat_clusters")

plot_dimplot(tumour, reduction = "umap", group_by = "Subtype", cols = pam50_subtype_colors,
             results_path = out_dir, filename = "DimPlot_UMAP_Clinical.png")
plot_dimplot(tumour, reduction = "umap", group_by = "PAM50_predicted", cols = pam50_subtype_colors,
             results_path = out_dir, filename = "DimPlot_UMAP_PAM50.png")
plot_dimplot(tumour, reduction = "umap", group_by = "seurat_clusters",
             results_path = out_dir, filename = "DimPlot_UMAP_Clusters.png")

concordance_df <- as.data.frame(table(Predicted = tumour$PAM50_predicted, Clinical = tumour$Subtype))
plot_transition_heatmap(concordance_df, x_col = "Clinical", y_col = "Predicted", count_col = "Freq",
                        title = "PAM50 predicted vs clinical Subtype concordance",
                        xlab = "Clinical Subtype", ylab = "PAM50 predicted",
                        filename = "Heatmap_PAM50_vs_ClinicalSubtype.png", results_path = out_dir,
                        mode = "gradient", fill_mode = "pct", label_mode = "n_pct",
                        pct_group_col = "Predicted", low_color = "white", high_color = "steelblue",
                        legend_name = "% of\npredicted", text_size = 3.5, width = 7, height = 5)

plot_composition(tumour, cluster_col = "orig.ident", group_by = "PAM50_predicted", colors = pam50_subtype_colors,
                 filename = "BarPlot_PAM50Composition_bySample.png", results_path = out_dir)
plot_composition(tumour, cluster_col = "orig.ident", group_by = "Subtype", colors = pam50_subtype_colors,
                 filename = "BarPlot_SubtypeComposition_bySample.png", results_path = out_dir)
plot_composition(tumour, cluster_col = "seurat_clusters", group_by = "PAM50_predicted", colors = pam50_subtype_colors,
                 filename = "BarPlot_PAM50Composition_byCluster.png", results_path = out_dir)
plot_composition(tumour, cluster_col = "seurat_clusters", group_by = "Subtype", colors = pam50_subtype_colors,
                 filename = "BarPlot_SubtypeComposition_byCluster.png", results_path = out_dir)
plot_composition(tumour, cluster_col = "orig.ident", group_by = "seurat_clusters",
                 filename = "BarPlot_ClusterComposition_bySample.png", results_path = out_dir)
plot_composition(tumour, cluster_col = "seurat_clusters", group_by = "orig.ident",
                 filename = "BarPlot_SampleComposition_byCluster.png", results_path = out_dir)

plot_subtype_composition_facet(tumour, out_dir, call_col = "PAM50_predicted", label = "PAM50 predicted",
                               colors = pam50_subtype_colors)

plot_score_boxplot(tumour, score_cols = pam50_score_cols, results_path = out_dir,
                   filename = "PAM50_scores_boxplot_byCluster.png", fill_labels = pam50_fill_labels,
                   colors = pam50_subtype_colors, title = "PAM50 module score distribution by cluster")
plot_score_boxplot(tumour, score_cols = pam50_score_cols, results_path = out_dir,
                   filename = "PAM50_scores_boxplot_bySample.png", cluster_col = "orig.ident", facet_col = "Subtype",
                   fill_labels = pam50_fill_labels, colors = pam50_subtype_colors,
                   title = "PAM50 module score distribution by sample")

pam50_prop_df <- compute_celltype_proportions(tumour, sample_id_col = "seurat_clusters",
                                              celltype_col = "PAM50_predicted", clinical_cols = "seurat_clusters")
plot_celltype_proportions_bar(pam50_prop_df, x_var = "PAM50_predicted", results_path = out_dir,
                              filename = "Barplot_PAM50Proportion_byCluster.png",
                              celltype_col = "seurat_clusters", colors = pam50_subtype_colors)

# =============================================================================
# 2. PROGENy pathway activity
# =============================================================================
DefaultAssay(tumour) <- "progeny"
progeny_groups <- as.list(setNames(rownames(tumour[["progeny"]]), rownames(tumour[["progeny"]])))
progeny_by <- c(bySubtype = "Subtype", byPAM50 = "PAM50_predicted", byCluster = "global_clusters")
for (nm in names(progeny_by)) {
  plot_marker_dotplot(tumour, group_by = progeny_by[[nm]], marker_groups = progeny_groups, results_path = out_dir,
                      filename = paste0("DotPlot_Progeny_", nm, ".png"))
}
plot_progeny_score_boxplots(tumour, out_dir)

progeny_diff <- FindAllMarkers(tumour, assay = "progeny", layer = "scale.data", group.by = "seurat_clusters",
                               test.use = "wilcox", logfc.threshold = 0, min.pct = 0,
                               mean.fxn = rowMeans, fc.name = "avg_diff", verbose = FALSE) %>%
  dplyr::rename(pathway = gene) %>%
  dplyr::select(cluster, pathway, avg_diff, p_val, p_val_adj, pct.1, pct.2)
write.csv(progeny_diff, file.path(out_dir, "Differential_PROGENy_byCluster.csv"), row.names = FALSE)
DefaultAssay(tumour) <- "RNA"

# =============================================================================
# 3. Cell cycle
# =============================================================================
plot_dimplot(tumour, reduction = "umap", group_by = "Phase", results_path = out_dir,
             filename = "DimPlot_CellCyclePhase.png")
plot_cellcycle_boxplot(tumour, out_dir)

# =============================================================================
# 4. scSubtype (scores computed in 10_tumour_scoring.R)
# =============================================================================
concordance_df <- as.data.frame(table(Predicted = tumour$scSubtype_call, Clinical = tumour$Subtype))
plot_transition_heatmap(concordance_df, x_col = "Clinical", y_col = "Predicted", count_col = "Freq",
                        title = "scSubtype predicted vs clinical Subtype concordance",
                        xlab = "Clinical Subtype", ylab = "scSubtype predicted",
                        filename = "Heatmap_scSubtype_vs_ClinicalSubtype.png", results_path = out_dir,
                        mode = "gradient", fill_mode = "pct", label_mode = "n_pct",
                        pct_group_col = "Predicted", low_color = "white", high_color = "steelblue",
                        legend_name = "% of\npredicted", text_size = 3.5, width = 7, height = 5)

trans_counts <- tumour@meta.data %>% dplyr::count(PAM50_predicted, scSubtype_call, name = "n")
plot_transition_heatmap(trans_counts, x_col = "scSubtype_call", y_col = "PAM50_predicted", count_col = "n",
                        title = "PAM50_predicted vs scSubtype_call",
                        xlab = "scSubtype_call", ylab = "PAM50_predicted",
                        filename = "PAM50_vs_scSubtype_call.png", results_path = out_dir,
                        mode = "gradient", fill_mode = "count", label_mode = "n",
                        low_color = "white", high_color = "darkorchid", legend_name = "n cells",
                        tile_border = "black", text_size = 3, width = 10, height = 10)

plot_dimplot(tumour, reduction = "umap", group_by = "scSubtype_call", cols = scsubtype_colors,
             results_path = out_dir, filename = "DimPlot_UMAP_scSubtype_call.png", label = TRUE)

for (fp in list(list(cols = scsubtype_score_cols, file = "FeaturePlot_scSubtype_scores.png"),
                list(cols = pam50_score_cols, file = "FeaturePlot_PAM50_scores.png"))) {
  p <- FeaturePlot(tumour, features = fp$cols, reduction = "umap", ncol = 2, order = TRUE,
                   cols = c("white", "lightgray", "darkblue"))
  ggsave(file.path(out_dir, fp$file), plot = p, width = 10, height = 8, dpi = 300, bg = "white")
}

plot_subtype_composition_facet(tumour, out_dir, call_col = "scSubtype_call", label = "scSubtype predicted",
                               colors = scsubtype_colors)
plot_score_boxplot(tumour, score_cols = scsubtype_score_cols, results_path = out_dir,
                   filename = "scSubtype_scores_boxplot_byCluster.png", colors = scsubtype_colors,
                   title = "scSubtype score distribution by cluster")
plot_score_boxplot(tumour, score_cols = scsubtype_score_cols, results_path = out_dir,
                   filename = "scSubtype_scores_boxplot_bySample.png", cluster_col = "orig.ident", facet_col = "Subtype",
                   colors = scsubtype_colors, title = "scSubtype score distribution by sample")

scsubtype_prop_df <- compute_celltype_proportions(tumour, sample_id_col = "seurat_clusters",
                                                  celltype_col = "scSubtype_call", clinical_cols = "seurat_clusters")
plot_celltype_proportions_bar(scsubtype_prop_df, x_var = "scSubtype_call", results_path = out_dir,
                              filename = "Barplot_scSubtypeProportion_byCluster.png",
                              celltype_col = "seurat_clusters", colors = scsubtype_colors)

# =============================================================================
# 5. Differential scores per cluster (cluster vs rest, Wilcoxon)
# =============================================================================
diff_sets <- list(scSubtype = scsubtype_score_cols,
                  PAM50     = c("PAM50_Lum", "PAM50_TNBC", "PAM50_Her2+"),
                  CellCycle = c("S.Score", "G2M.Score"))
for (nm in names(diff_sets)) {
  csv <- file.path(out_dir, paste0("Differential_", nm, "_byCluster.csv"))
  write.csv(differential_scores_by_cluster(tumour, diff_sets[[nm]]), csv, row.names = FALSE)
  plot_differential_barplot(csv, "score", paste("Differential", nm, "scores by cluster"),
                            paste0("Barplot_Differential_", nm, ".png"), out_dir)
}

plot_progeny_diff_heatmap(progeny_diff, out_dir, filename = "Heatmap_Differential_PROGENy_byCluster.png",
                          title = "Differential PROGENy pathway activity per cluster (vs all other clusters)")
plot_differential_barplot(file.path(out_dir, "Differential_PROGENy_byCluster.csv"),
                          "pathway", "Differential PROGENy pathway activity by cluster",
                          "Barplot_Differential_PROGENy.png", out_dir)

# =============================================================================
# 6. Origin of the tumour clusters (global cluster, sample)
# =============================================================================
trans_counts <- data.frame(cluster = tumour$global_clusters, tumor_cluster = tumour$seurat_clusters) %>%
  dplyr::count(tumor_cluster, cluster, name = "n")
plot_transition_heatmap(trans_counts, x_col = "tumor_cluster", y_col = "cluster", count_col = "n",
                        title = "Origin of tumor cluster cells: cluster vs tumor cluster",
                        xlab = "Tumor cluster", ylab = "Cluster",
                        filename = "Transition_Cluster_x_TumorCluster.png", results_path = out_dir,
                        mode = "gradient", fill_mode = "count", label_mode = "n",
                        low_color = "white", high_color = "firebrick", width = 10, height = 8)

trans_counts <- data.frame(sample = tumour$orig.ident, cluster = tumour$seurat_clusters,
                           Subtype = as.character(tumour$Subtype)) %>%
  dplyr::count(sample, cluster, Subtype, name = "n")
plot_transition_heatmap(trans_counts, x_col = "cluster", y_col = "sample", count_col = "n",
                        title = "Origin of tumor cluster cells: sample vs cluster, by subtype",
                        xlab = "Tumor cluster", ylab = "Sample",
                        filename = "Transition_Sample_x_TumorCluster.png", results_path = out_dir,
                        mode = "manual", group_col = "Subtype", group_colors = pam50_subtype_colors,
                        label_mode = "n", text_size = 2.5, width = 16, height = 6)
rm(tumour); invisible(gc())

# =============================================================================
# 7. Celltype composition per sample, sample vs rest (full object)
# =============================================================================
annotated_meta <- read_run_object("notumor_annotated")@meta.data
comp_diff <- differential_composition(annotated_meta)
write.csv(comp_diff, file.path(out_dir, "Differential_celltypeComposition_bySample.csv"), row.names = FALSE)

p_comp_diff <- ggplot(comp_diff, aes(x = celltype, y = pct_diff, fill = pct_diff)) +
  geom_col() +
  facet_wrap(~cluster, scales = "free_x") +
  scale_fill_gradient2(low = "firebrick", mid = "white", high = "forestgreen", midpoint = 0, name = "Diff. vs\nrest (%)") +
  coord_flip() +
  theme_bw() + theme(panel.grid = element_blank()) +
  labs(title = "Differential celltype composition by sample", x = NULL, y = "pct_in - pct_out")
ggsave(file.path(out_dir, "Barplot_Differential_celltypeComposition.png"), p_comp_diff,
       width = 14, height = 10, dpi = 300, bg = "white")

cat("\n ---- FINISHED 11_tumour_markers: fill annotation/tumour.tsv (label + evidence; DROP for noise), then run 12_tumour_annotate.R ----\n")
