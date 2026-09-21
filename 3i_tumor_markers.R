##
##  Single Cell Analysis Step 3i: Tumor cluster marker & validation plots
##  Runs AFTER 3h_tumor_scoring.R — takes tumoral_scored.rds as input.
##  Scores scSubtype, validates PAM50/scSubtype against the clinical subtype,
##  and produces the differential scores, transition heatmaps and composition
##  plots needed to annotate the tumor clusters in 3j_tumor_annotate.R.
##

library(Seurat)
library(dplyr)
library(progeny)
library(tidyr)
library(tibble)
library(ggplot2)
library(corrplot)
library(gt)

# --- Paths ---
project_path <- "/home/user/PROJECTS/scRNA_vmendez/"
wd <- paste0(project_path, "codes/")
setwd(wd)
results_path <- paste0(project_path, "results/")
results_GEMX_CA_path <- paste0(results_path, "CellAnnotation/")
results_GEMX_TUMOR_path <- paste0(results_GEMX_CA_path, "Tumor/")
results_FINAL_path <- paste0(results_path, "Final_Annotated/")

# --- Plot functions, helpers & scSubtype gene signatures ---
source(paste0(wd, "TCA_plots.R"))
source(paste0(wd, "CA_plots.R"))
source(paste0(wd, "CL_plots.R"))
source(paste0(wd, "utils.R"))
source(paste0(wd, "scsubtype_signatures.R"))

# --- Load scored tumor data and full annotated object ---
dwTumoral <- readRDS(paste0(results_GEMX_TUMOR_path, "tumoral_scored.rds"))
cat("\n Tumor scored data loaded \n")
dwAnnotated <- readRDS(paste0(results_GEMX_CA_path, "CopyKAT/copykat_annotated_data.rds"))
cat("\n Full annotated data loaded \n")

# --- PAM50-like marker genes (for the marker dotplot below) ---
pam50_genes <- list(
  "Her2+" = c("ERBB2", "GRB7", "BLVRA", "TMEM45B", "FGFR4"),
  "Lum" = c("ESR1", "PGR", "BAG1", "MAPT", "NAT1", "ZIP6"),
  "TNBC" = c("MKI67", "CCNE1", "ANLN", "CDC20", "EGFR", "MYC")
)

# =============================================================================
# 1. PAM50 validation & composition plots
# =============================================================================

# --- HTML marker tables per cluster ---
# plot_marker_tables(results_GEMX_TUMOR_path)
cat("\n Cluster marker tables saved \n")

# --- PAM50 dotplot ---
plot_marker_dotplot(dwTumoral, marker_groups = pam50_genes,
                    results_path = results_GEMX_TUMOR_path,
                    filename = "Dotplot_PAM50.png",
                    group_by = "seurat_clusters")

# --- Validation: predicted subtype vs clinical subtype ---
plot_dimplot(dwTumoral, reduction = "umap", group_by = "Subtype", cols = pam50_subtype_colors,
             results_path = results_GEMX_TUMOR_path, filename = "DimPlot_UMAP_Clinical.png")
plot_dimplot(dwTumoral, reduction = "umap", group_by = "PAM50_predicted", cols = pam50_subtype_colors,
             results_path = results_GEMX_TUMOR_path, filename = "DimPlot_UMAP_PAM50.png")
plot_dimplot(dwTumoral, reduction = "umap", group_by = "seurat_clusters",
             results_path = results_GEMX_TUMOR_path, filename = "DimPlot_UMAP_Clusters.png")

# --- PAM50 vs clinical Subtype concordance heatmap ---
concordance_df <- as.data.frame(table(Predicted = dwTumoral$PAM50_predicted, Clinical = dwTumoral$Subtype))
plot_transition_heatmap(concordance_df, x_col = "Clinical", y_col = "Predicted", count_col = "Freq",
                        title = "PAM50 predicted vs clinical Subtype concordance",
                        xlab = "Clinical Subtype", ylab = "PAM50 predicted",
                        filename = "Heatmap_PAM50_vs_ClinicalSubtype.png", results_path = results_GEMX_TUMOR_path,
                        mode = "gradient", fill_mode = "pct", label_mode = "n_pct",
                        pct_group_col = "Predicted",
                        low_color = "white", high_color = "steelblue",
                        legend_name = "% of\npredicted", text_size = 3.5,
                        width = 7, height = 5)

# --- Cluster composition by sample and subtype ---
plot_composition(dwTumoral, cluster_col = "orig.ident", group_by = "PAM50_predicted", colors = pam50_subtype_colors,
                         filename = "BarPlot_PAM50Composition_bySample.png", results_path = results_GEMX_TUMOR_path)
plot_composition(dwTumoral, cluster_col = "orig.ident", group_by = "Subtype", colors = pam50_subtype_colors,
                         filename = "BarPlot_SubtypeComposition_bySample.png", results_path = results_GEMX_TUMOR_path)
plot_composition(dwTumoral, cluster_col = "seurat_clusters", group_by = "PAM50_predicted", colors = pam50_subtype_colors,
                         filename = "BarPlot_PAM50Composition_byCluster.png", results_path = results_GEMX_TUMOR_path)
plot_composition(dwTumoral, cluster_col = "seurat_clusters", group_by = "Subtype", colors = pam50_subtype_colors,
                         filename = "BarPlot_SubtypeComposition_byCluster.png", results_path = results_GEMX_TUMOR_path)
plot_composition(dwTumoral, cluster_col = "orig.ident", group_by = "seurat_clusters",
                         filename = "BarPlot_ClusterComposition_bySample.png", results_path = results_GEMX_TUMOR_path)
plot_composition(dwTumoral, cluster_col = "seurat_clusters", group_by = "orig.ident",
                         filename = "BarPlot_SampleComposition_byCluster.png", results_path = results_GEMX_TUMOR_path)

# --- PAM50 composition stacked bars, faceted by subtype ---
plot_subtype_composition_facet(dwTumoral, results_GEMX_TUMOR_path,
                               call_col = "PAM50_predicted", label = "PAM50 predicted",
                               colors = pam50_subtype_colors)

# --- PAM50 score distribution by cluster, fill by PAM50 module ---
plot_score_boxplot(dwTumoral, score_cols = c("PAM50_Her2+", "PAM50_Lum", "PAM50_TNBC"),
                   results_path = results_GEMX_TUMOR_path,
                   filename = "PAM50_scores_boxplot_byCluster.png",
                   fill_labels = c("PAM50_Her2+" = "Her2+", "PAM50_Lum" = "Lum", "PAM50_TNBC" = "TNBC"),
                   colors = pam50_subtype_colors,
                   title = "PAM50 module score distribution by cluster")

# --- PAM50 score distribution by sample, fill by PAM50 module ---
plot_score_boxplot(dwTumoral, score_cols = c("PAM50_Her2+", "PAM50_Lum", "PAM50_TNBC"),
                   results_path = results_GEMX_TUMOR_path,
                   filename = "PAM50_scores_boxplot_bySample.png",
                   cluster_col = "orig.ident", facet_col = "Subtype",
                   fill_labels = c("PAM50_Her2+" = "Her2+", "PAM50_Lum" = "Lum", "PAM50_TNBC" = "TNBC"),
                   colors = pam50_subtype_colors,
                   title = "PAM50 module score distribution by sample")

# --- PAM50 call proportion within each sample's total tumor cells, faceted by cluster ---
pam50_prop_df <- compute_celltype_proportions(dwTumoral, sample_id_col = "seurat_clusters",
                                              celltype_col = "PAM50_predicted",
                                              clinical_cols = "seurat_clusters")
plot_celltype_proportions_bar(pam50_prop_df, x_var = "PAM50_predicted",
                              results_path = results_GEMX_TUMOR_path,
                              filename = "Barplot_PAM50Proportion_byCluster.png",
                              celltype_col = "seurat_clusters", colors = pam50_subtype_colors)


scsubtype_prop_df <- compute_celltype_proportions(dwTumoral, sample_id_col = "seurat_clusters",
                                                  celltype_col = "scSubtype_call",
                                                  clinical_cols = "seurat_clusters")
plot_celltype_proportions_bar(scsubtype_prop_df, x_var = "scSubtype_call",
                              results_path = results_GEMX_TUMOR_path,
                              filename = "Barplot_scSubtypeProportion_byCluster.png",
                              celltype_col = "seurat_clusters", colors = scsubtype_colors)

# =============================================================================
# 2. PROGENy pathway activity
# =============================================================================

DefaultAssay(dwTumoral) <- "progeny"

# --- Pathway activity dotplots ---
plot_marker_dotplot(dwTumoral, group_by = "Subtype",
                    marker_groups = as.list(setNames(rownames(dwTumoral[["progeny"]]),
                                                     rownames(dwTumoral[["progeny"]]))),
                    filename = "DotPlot_Progeny_bySubtype.png", results_path = results_GEMX_TUMOR_path)
plot_marker_dotplot(dwTumoral, group_by = "PAM50_predicted",
                    marker_groups = as.list(setNames(rownames(dwTumoral[["progeny"]]),
                                                     rownames(dwTumoral[["progeny"]]))),
                    filename = "DotPlot_Progeny_byPAM50.png", results_path = results_GEMX_TUMOR_path)
plot_marker_dotplot(dwTumoral, group_by = "decontX_clusters",
                    marker_groups = as.list(setNames(rownames(dwTumoral[["progeny"]]),
                                                     rownames(dwTumoral[["progeny"]]))),
                    filename = "DotPlot_Progeny_byCluster.png", results_path = results_GEMX_TUMOR_path)

# --- Pathway score distribution boxplots (per cluster and per pathway) ---
plot_progeny_score_boxplots(dwTumoral, results_GEMX_TUMOR_path)

# --- Differential pathway activity per cluster ---
progeny_diff <- FindAllMarkers(dwTumoral,
                               assay = "progeny",
                               slot = "scale.data",
                               group.by = "seurat_clusters",
                               test.use = "wilcox",
                               logfc.threshold = 0,
                               min.pct = 0)
progeny_diff <- progeny_diff %>%
  dplyr::rename(pathway = gene, avg_diff = avg_log2FC) %>%
  dplyr::select(cluster, pathway, avg_diff, p_val, p_val_adj, pct.1, pct.2)
write.csv(progeny_diff, paste0(results_GEMX_TUMOR_path, "Differential_PROGENy_byCluster.csv"), row.names = FALSE)
cat("\n Differential PROGENy pathways per cluster done \n")

DefaultAssay(dwTumoral) <- "RNA_decontX"


# =============================================================================
# 3. Cell cycle & CopyKAT QC plots
# =============================================================================

# --- Cell cycle phase UMAP and boxplot ---
plot_dimplot(dwTumoral, reduction = "umap", group_by = "Phase",
             results_path = results_GEMX_TUMOR_path, filename = "DimPlot_CellCyclePhase.png")
plot_cellcycle_boxplot(dwTumoral, results_GEMX_TUMOR_path)

# --- CopyKAT CNV: prediction composition and CNA burden per cluster ---
plot_copykat_composition(dwTumoral, results_GEMX_TUMOR_path)
plot_copykat_cna_boxplot(dwTumoral, results_GEMX_TUMOR_path, score_col = "CIN_score",
                         y_lab = "CIN score (SD of CNA profile)")
plot_copykat_cna_boxplot(dwTumoral, results_GEMX_TUMOR_path, score_col = "PGA_score",
                         y_lab = "PGA score (% bins altered)")

plot_dimplot(dwTumoral, reduction = "umap", group_by = "copykat_prediction",
             results_path = results_GEMX_TUMOR_path, filename = "DimPlot_UMAP_CopyKATPrediction_Tumoral.png")
plot_composition(dwTumoral, cluster_col = "seurat_clusters", group_by = "copykat_prediction",
                 filename = "BarPlot_CopyKATPrediction_Tumor_byCluster.png", results_path = results_GEMX_TUMOR_path)

# --- CNA-based UMAP (reduction 'umap_cna', inherited from 3g_copykat.R) ---
plot_dimplot(dwTumoral, reduction = "umap_cna", group_by = "copykat_prediction",
             results_path = results_GEMX_TUMOR_path, filename = "DimPlot_UMAP_CNA_CopyKATPrediction_Tumoral.png")
plot_dimplot(dwTumoral, reduction = "umap_cna", group_by = "seurat_clusters",
             results_path = results_GEMX_TUMOR_path, filename = "DimPlot_UMAP_CNA_SeuratClusters_Tumoral.png")
plot_dimplot(dwTumoral, reduction = "umap_cna", group_by = "orig.ident",
             results_path = results_GEMX_TUMOR_path, filename = "DimPlot_UMAP_CNA_Sample_Tumoral.png")
plot_dimplot(dwTumoral, reduction = "umap_cna", group_by = "Subtype",
             results_path = results_GEMX_TUMOR_path, filename = "DimPlot_UMAP_CNA_ClinicalSubtype_Tumoral.png")
plot_dimplot(dwTumoral, reduction = "umap_cna", group_by = "PAM50_predicted", cols = pam50_subtype_colors,
             results_path = results_GEMX_TUMOR_path, filename = "DimPlot_UMAP_CNA_PAM50_Tumoral.png")


# =============================================================================
# 4. scSubtype scoring
# =============================================================================

dwTumoral <- compute_scsubtype_scores(dwTumoral, scsubtype_signatures, assay = "RNA_decontX")
scSubtype_score_cols <- names(scsubtype_signatures)

# --- scSubtype vs clinical Subtype concordance heatmap (mirrors the PAM50 one above) ---
concordance_df <- as.data.frame(table(Predicted = dwTumoral$scSubtype_call, Clinical = dwTumoral$Subtype))
plot_transition_heatmap(concordance_df, x_col = "Clinical", y_col = "Predicted", count_col = "Freq",
                        title = "scSubtype predicted vs clinical Subtype concordance",
                        xlab = "Clinical Subtype", ylab = "scSubtype predicted",
                        filename = "Heatmap_scSubtype_vs_ClinicalSubtype.png", results_path = results_GEMX_TUMOR_path,
                        mode = "gradient", fill_mode = "pct", label_mode = "n_pct",
                        pct_group_col = "Predicted",
                        low_color = "white", high_color = "steelblue",
                        legend_name = "% of\npredicted", text_size = 3.5,
                        width = 7, height = 5)

# --- PAM50_predicted vs scSubtype_call transition heatmap ---
trans_counts <- dwTumoral@meta.data %>%
  dplyr::count(PAM50_predicted, scSubtype_call, name = "n")

plot_transition_heatmap(trans_counts, x_col = "scSubtype_call", y_col = "PAM50_predicted", count_col = "n",
                        title = "PAM50_predicted vs scSubtype_call",
                        xlab = "scSubtype_call (new)", ylab = "PAM50_predicted (old)",
                        filename = "PAM50_vs_scSubtype_call.png", results_path = results_GEMX_TUMOR_path,
                        mode = "gradient", fill_mode = "count", label_mode = "n",
                        low_color = "white", high_color = "darkorchid",
                        legend_name = "n cells", tile_border = "black", text_size = 3,
                        width = 10, height = 10)

# --- UMAP and score feature plots ---
plot_dimplot(dwTumoral, reduction = "umap", group_by = "scSubtype_call", cols = scsubtype_colors,
             results_path = results_GEMX_TUMOR_path,
             filename = "DimPlot_UMAP_scSubtype_call.png", label = TRUE)
plot_dimplot(dwTumoral, reduction = "umap_cna", group_by = "scSubtype_call", cols = scsubtype_colors,
             results_path = results_GEMX_TUMOR_path,
             filename = "DimPlot_UMAP_CNA_scSubtype_call_Tumoral.png")

pam50_score_names <- c("PAM50_Her2+", "PAM50_Lum", "PAM50_TNBC")

p <- FeaturePlot(dwTumoral, features = scSubtype_score_cols, reduction = "umap",
                  ncol = 2, order = TRUE) &
  scale_color_gradientn(colors = c("white", "lightgray", "darkblue"))
ggsave(paste0(results_GEMX_TUMOR_path, "FeaturePlot_scSubtype_scores.png"), plot = p, width = 10, height = 8, dpi = 300)

p <- FeaturePlot(dwTumoral, features = pam50_score_names, reduction = "umap",
                  ncol = 2, order = TRUE) &
  scale_color_gradientn(colors = c("white", "lightgray", "darkblue"))
ggsave(paste0(results_GEMX_TUMOR_path, "FeaturePlot_PAM50_scores.png"), plot = p, width = 10, height = 8, dpi = 300)

# --- scSubtype composition stacked bars, faceted by subtype (mirrors the PAM50 one above) ---
plot_subtype_composition_facet(dwTumoral, results_GEMX_TUMOR_path,
                               call_col = "scSubtype_call", label = "scSubtype predicted",
                               colors = scsubtype_colors)

# --- scSubtype score distribution by cluster, fill by scSubtype signature ---
plot_score_boxplot(dwTumoral, score_cols = scSubtype_score_cols,
                   results_path = results_GEMX_TUMOR_path,
                   filename = "scSubtype_scores_boxplot_byCluster.png",
                   colors = scsubtype_colors,
                   title = "scSubtype score distribution by cluster")

# --- scSubtype score distribution by sample, fill by scSubtype signature ---
plot_score_boxplot(dwTumoral, score_cols = scSubtype_score_cols,
                   results_path = results_GEMX_TUMOR_path,
                   filename = "scSubtype_scores_boxplot_bySample.png",
                   cluster_col = "orig.ident", facet_col = "Subtype",
                   colors = scsubtype_colors,
                   title = "scSubtype score distribution by sample")



# =============================================================================
# 5. Differential module scores per cluster (PAM50, scSubtype, cell cycle, copykat)
# =============================================================================

scSubtype_diff <- differential_scores_by_cluster(dwTumoral, scSubtype_score_cols)
write.csv(scSubtype_diff, paste0(results_GEMX_TUMOR_path, "Differential_scSubtype_byCluster.csv"), row.names = FALSE)
cat("\n Differential scSubtype module scores per cluster done \n")

pam50_score_cols <- c("PAM50_Lum", "PAM50_TNBC", "PAM50_Her2+")
pam50_diff <- differential_scores_by_cluster(dwTumoral, pam50_score_cols)
write.csv(pam50_diff, paste0(results_GEMX_TUMOR_path, "Differential_PAM50_byCluster.csv"), row.names = FALSE)
cat("\n Differential PAM50 module scores per cluster done \n")

cellcycle_score_cols <- c("S.Score", "G2M.Score")
cellcycle_diff <- differential_scores_by_cluster(dwTumoral, cellcycle_score_cols)
write.csv(cellcycle_diff, paste0(results_GEMX_TUMOR_path, "Differential_CellCycle_byCluster.csv"), row.names = FALSE)
cat("\n Differential cell cycle scores per cluster done \n")

# Cells without a CopyKAT profile (CIN_score/PGA_score NA) are excluded
copykat_score_cols <- c("CIN_score", "PGA_score")
dwTumoral_cna <- subset(dwTumoral, cells = colnames(dwTumoral)[!is.na(dwTumoral$CIN_score)])
copykat_diff <- differential_scores_by_cluster(dwTumoral_cna, copykat_score_cols)
write.csv(copykat_diff, paste0(results_GEMX_TUMOR_path, "Differential_Copykat_byCluster.csv"), row.names = FALSE)
cat("\n Differential CIN/PGA scores per cluster done \n")

# --- Heatmap of differential pathway activity ---
plot_progeny_diff_heatmap(progeny_diff, results_GEMX_TUMOR_path,
                          filename = "Heatmap_Differential_PROGENy_byCluster.png",
                          title = "Differential PROGENy pathway activity per cluster (vs all other clusters)")

# --- Differential score barplots per cluster ---
plot_differential_barplot(paste0(results_GEMX_TUMOR_path, "Differential_PAM50_byCluster.csv"),
                          "score", "Differential PAM50 module scores by cluster",
                          "Barplot_Differential_PAM50.png", results_GEMX_TUMOR_path)
plot_differential_barplot(paste0(results_GEMX_TUMOR_path, "Differential_CellCycle_byCluster.csv"),
                          "score", "Differential cell cycle scores by cluster",
                          "Barplot_Differential_CellCycle.png", results_GEMX_TUMOR_path)
plot_differential_barplot(paste0(results_GEMX_TUMOR_path, "Differential_PROGENy_byCluster.csv"),
                          "pathway", "Differential PROGENy pathway activity by cluster",
                          "Barplot_Differential_PROGENy.png", results_GEMX_TUMOR_path)
plot_differential_barplot(paste0(results_GEMX_TUMOR_path, "Differential_Copykat_byCluster.csv"),
                          "score", "Differential Copykat scores by cluster",
                          "Barplot_Differential_Copykat.png", results_GEMX_TUMOR_path)
plot_differential_barplot(paste0(results_GEMX_TUMOR_path, "Differential_scSubtype_byCluster.csv"),
                          "score", "Differential scSubtype scores by cluster",
                          "Barplot_Differential_scSubtype.png", results_GEMX_TUMOR_path)


# =============================================================================
# 6. Cluster / sample transition heatmaps (tumor cluster level)
# =============================================================================

trans_counts <- data.frame(
  cluster = dwTumoral$decontX_clusters,
  tumor_cluster = dwTumoral$seurat_clusters
) %>% dplyr::count(tumor_cluster, cluster, name = "n")

plot_transition_heatmap(trans_counts, x_col = "tumor_cluster", y_col = "cluster", count_col = "n",
                        title = "Origin of tumor cluster cells: cluster vs tumor cluster",
                        xlab = "Tumor cluster", ylab = "Cluster",
                        filename = "Transition_Cluster_x_TumorCluster.png", results_path = results_GEMX_TUMOR_path,
                        mode = "gradient", fill_mode = "count", label_mode = "n",
                        low_color = "white", high_color = "firebrick",
                        width = 10, height = 8)

trans_counts <- data.frame(
  sample = dwTumoral$orig.ident,
  cluster = dwTumoral$seurat_clusters,
  Subtype = as.character(dwTumoral$Subtype)
) %>% dplyr::count(sample, cluster, Subtype, name = "n")

plot_transition_heatmap(trans_counts, x_col = "cluster", y_col = "sample", count_col = "n",
                        title = "Origin of tumor cluster cells: sample vs cluster, by subtype",
                        xlab = "Tumor cluster", ylab = "Sample",
                        filename = "Transition_Sample_x_TumorCluster.png", results_path = results_GEMX_TUMOR_path,
                        mode = "manual", group_col = "Subtype", group_colors = pam50_subtype_colors,
                        label_mode = "n", text_size = 2.5,
                        width = 16, height = 6)


# =============================================================================
# 7. Celltype composition & differential composition (full annotated object)
# =============================================================================

composition_by_cluster <- dwAnnotated@meta.data %>%
  group_by(orig.ident, celltype) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

comp_diff <- differential_composition(dwAnnotated@meta.data)
write.csv(comp_diff, paste0(results_path, "Differential_celltypeComposition_bySample.csv"), row.names = FALSE)

p_comp_diff <- ggplot(comp_diff, aes(x = celltype, y = pct_diff, fill = pct_diff)) +
  geom_col() +
  facet_wrap(~cluster, scales = "free_x") +
  scale_fill_gradient2(low = "firebrick", mid = "white", high = "forestgreen", midpoint = 0,
                        name = "Diff. vs\nrest (%)") +
  coord_flip() +
  theme_bw() + theme(panel.grid = element_blank()) +
  labs(title = "Differential celltype composition by sample",
       x = NULL, y = "pct_in - pct_out")

ggsave(paste0(results_GEMX_CA_path, "Barplot_Differential_celltypeComposition.png"), p_comp_diff,
       width = 14, height = 10, dpi = 300, bg = "white")


# --- Export the scored tumor object (scSubtype scores added) ---
saveRDS(dwTumoral, file.path(results_GEMX_TUMOR_path, "tumoral_scored.rds"))
cat("\n Tumor object re-saved with scSubtype scores \n")

cat(paste("\n ---- FINISHED TUMOR MARKER & VALIDATION PLOTS ----
    Review the plots and cluster CSVs, fill in `clusters_tumor_annotated` and
    the noise cluster lists at the top of 3j_tumor_annotate.R, then run it.
    Generated files:
        · tumoral_scored.rds  (+ scSubtype_call and scSubtype scores)
        · Differential_(PROGENy|PAM50|scSubtype|CellCycle|Copykat)_byCluster.csv
        · Differential_celltypeComposition_bySample.csv
    Generated plots:
        · Heatmap_(PAM50|scSubtype)_vs_ClinicalSubtype.png, PAM50_vs_scSubtype_call.png
        · DimPlot_UMAP_(PAM50|scSubtype_call|Clusters|Clinical).png
        · FeaturePlot_(PAM50|scSubtype)_scores.png
        · BarPlot_(PAM50|Subtype|Cluster|Sample)Composition_(byCluster|bySample).png
        · Composition_(PAM50_predicted|scSubtype_call)_bySample_facetSubtype.png
        · Barplot_(PAM50|scSubtype)Proportion_byCluster.png
        · (PAM50|scSubtype)_scores_boxplot_by(Cluster|Sample).png
        · Barplot_Differential_(PROGENy|PAM50|CellCycle|Copykat|scSubtype).png
        · Transition_(Cluster|Sample)_x_TumorCluster.png
        "))
