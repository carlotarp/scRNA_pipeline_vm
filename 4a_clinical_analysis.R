##
##  Single Cell Analysis Step 4a: Clinical covariate analysis
##  Runs AFTER 3j_tumor_annotate.R — takes fully_annotated_data.rds and
##  tumoral_annotated.rds as input.
##

library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tidyr)
library(tibble)
library(ggsci)
library(ggplot2)
library(Seurat)


# Set paths
project_path <- "/home/user/PROJECTS/scRNA_vmendez/"
wd <- paste0(project_path, "codes/")
setwd(wd)
results_path <- paste0(project_path, "results/")
results_GEMX_CA_path <- paste0(results_path, "CellAnnotation/")
results_GEMX_TUMOR_path <- paste0(results_GEMX_CA_path, "Tumor/")
results_GEMX_CLINICAL_path <- paste0(results_GEMX_CA_path, "Clinical/")

source(paste0(wd, "TCA_plots.R"))
source(paste0(wd, "CA_plots.R"))
source(paste0(wd, "CL_plots.R"))
source(paste0(wd, "Clinical_plots.R"))
source(paste0(wd, "utils.R"))


dwAnnotated <- readRDS(paste0(results_GEMX_CA_path, "fully_annotated_data.rds"))
dwAnnotated <- subset(dwAnnotated, subset = celltype != "Noise")
dwTumoral <- readRDS(paste0(results_GEMX_CA_path, "Tumor/tumoral_annotated.rds"))
cat("\n Full annotated data loaded \n")

dwNoTumoral <- subset(dwAnnotated, subset = lineage != "Tumor")
dwLeukos <- subset(dwAnnotated, subset = lineage == "Leukocytes")
dwStroma <- subset(dwAnnotated, subset = lineage == "Stromal")

generate_heatmap(dwAnnotated, filename = "HeatMapComposition_All.png", n_clusters = 6)
generate_heatmap(dwNoTumoral, filename = "HeatMapComposition_NoTumor.png", n_clusters = 4)
generate_heatmap(dwLeukos, filename = "HeatMapComposition_Leukocytes.png", n_clusters = 4)
generate_heatmap(dwStroma, filename = "HeatMapComposition_Stromal.png", n_clusters = 4)
generate_heatmap(dwTumoral, filename = "HeatMapComposition_Tumor.png", n_clusters = 6)



# =============================================================================
# IGG scoring (plasma cell / immunoglobulin signature) & validation
# =============================================================================

igg_genes <- list(
  IGG = c("CD27", "CD79A", "HLA-C", "JCHAIN", "IGKC", "IGLV3-25", "IL2RG",
          "CXCL8", "LAX1", "NTN3", "PIM2", "POU2AF1", "TNFRSF17")
)

# Compute IGG module score (single signature -> one continuous score column)
dwAnnotated <- AddModuleScore(dwAnnotated, features = igg_genes, name = "IGGscore_")
colnames(dwAnnotated@meta.data)[colnames(dwAnnotated@meta.data) == "IGGscore_1"] <- "IGG_score"

# Predicted IGG level = tertile of the score
igg_tertiles <- quantile(dwAnnotated$IGG_score, probs = c(1 / 3, 2 / 3))
dwAnnotated$IGG_predicted <- cut(dwAnnotated$IGG_score,
                               breaks = c(-Inf, igg_tertiles, Inf),
                               labels = c("Low", "Medium", "High"))
cat("\n IGG module score computed \n")
print(table(dwAnnotated$IGG_predicted))

igg_score_cols <- c("IGG_score")

# --- Validation: predicted IGG level vs clinical IGG level ---
plot_dimplot(dwAnnotated, reduction = "umap_decontX", group_by = "IGG", cols = igg_level_colors,
             results_path = results_GEMX_CLINICAL_path, filename = "DimPlot_UMAP_IGG_Clinical.png")
plot_dimplot(dwAnnotated, reduction = "umap_decontX", group_by = "IGG_predicted", cols = igg_level_colors,
             results_path = results_GEMX_CLINICAL_path, filename = "DimPlot_UMAP_IGG_predicted.png")

# --- IGG predicted vs clinical IGG concordance heatmap ---
igg_concordance_df <- as.data.frame(table(Predicted = dwAnnotated$IGG_predicted, Clinical = dwAnnotated$IGG))
plot_transition_heatmap(igg_concordance_df, x_col = "Clinical", y_col = "Predicted", count_col = "Freq",
                        title = "IGG predicted vs clinical IGG concordance",
                        xlab = "Clinical IGG", ylab = "IGG predicted",
                        filename = "Heatmap_IGG_vs_ClinicalIGG.png", results_path = results_GEMX_CLINICAL_path,
                        mode = "gradient", fill_mode = "pct", label_mode = "n_pct",
                        pct_group_col = "Predicted",
                        low_color = "white", high_color = "steelblue",
                        legend_name = "% of\npredicted", text_size = 3.5,
                        width = 7, height = 5)

# --- IGG score feature plot ---
p_igg_score <- FeaturePlot(dwAnnotated, features = "IGG_score", reduction = "umap_decontX", order = TRUE) &
  scale_color_gradientn(colors = c("lightgray", "blue", "darkblue"))
ggsave(paste0(results_GEMX_CLINICAL_path, "FeaturePlot_IGG_score.png"), plot = p_igg_score, width = 6, height = 5, dpi = 300)

# --- Cluster composition by sample and IGG level ---
plot_composition(dwAnnotated, cluster_col = "orig.ident", group_by = "IGG", colors = igg_level_colors,
                 filename = "BarPlot_IGGComposition_bySample.png", results_path = results_GEMX_CLINICAL_path)
plot_composition(dwAnnotated, cluster_col = "celltype", group_by = "IGG", colors = igg_level_colors,
                 filename = "BarPlot_IGGComposition_byCelltype.png", results_path = results_GEMX_CLINICAL_path)

# --- IGG predicted composition stacked bars, faceted by clinical IGG level ---
plot_subtype_composition_facet(dwAnnotated, results_GEMX_CLINICAL_path,
                               call_col = "IGG_predicted", facet_col = "IGG",
                               label = "IGG predicted", colors = igg_level_colors)

# --- IGG score distribution by cluster, fill by clinical Subtype, faceted by subtype ---
plot_score_boxplot(dwAnnotated, score_cols = "IGG_score",
                   results_path = results_GEMX_CLINICAL_path,
                   filename = "IGG_score_boxplot_byCelltypeSubtype.png",
                   cluster_col = "celltype", fill_col = "Subtype", 
                   colors = pam50_subtype_colors,
                   title = "IGG score distribution by cluster, by clinical IGG level, by subtype")

# --- IGG score distribution by cluster, fill by clinical IGG level, faceted by subtype ---
plot_score_boxplot(dwAnnotated, score_cols = "IGG_score",
                   results_path = results_GEMX_CLINICAL_path,
                   filename = "IGG_score_boxplot_byCelltypeIGG.png",
                   cluster_col = "celltype", fill_col = "IGG", 
                   colors = igg_level_colors,
                   title = "IGG score distribution by cluster, by clinical IGG level, by subtype")

# --- IGG score distribution by sample, fill by clinical IGG level, faceted by subtype ---
plot_score_boxplot(dwAnnotated, score_cols = "IGG_score",
                   results_path = results_GEMX_CLINICAL_path,
                   filename = "IGG_score_boxplot_bySample_facetSubtype.png",
                   cluster_col = "orig.ident",
                   fill_col = "IGG", facet_col = "Subtype",
                   colors = igg_level_colors,
                   title = "IGG score distribution by sample, by clinical IGG level, by subtype")

# --- Differential IGG module score by cluster ---
igg_diff <- differential_scores_by_cluster(dwAnnotated, igg_score_cols, group_by = "decontX_clusters")
write.csv(igg_diff, paste0(results_GEMX_CLINICAL_path, "Differential_IGG_byCluster.csv"), row.names = FALSE)
cat("\n Differential IGG module score per cluster done \n")

plot_differential_barplot(paste0(results_GEMX_CLINICAL_path, "Differential_IGG_byCluster.csv"),
                          "score", "Differential IGG module score by cluster",
                          "Barplot_Differential_IGG.png", results_GEMX_CLINICAL_path)

# =============================================================================
# Clinical covariate association: composition, proportion & correlation
# =============================================================================

dwAnnotated$TLS_binary <- ifelse(dwAnnotated$TLS == "0", "No", "Yes")
dwTumoral$TLS_binary <- ifelse(dwTumoral$TLS == "0", "No", "Yes")

# --- Clinical covariates: palette + ordinal encoding (NULL = no natural order) ---
clinical_vars <- list(
  Subtype    = list(colors = pam50_subtype_colors, ordinal = NULL),
  TLS        = list(colors = tls_colors,           ordinal = c("0" = 0, "1a" = 1, "1b" = 2, "2a" = 3, "2b" = 4)),
  TLS_binary = list(colors = c("No" = unname(tls_colors["0"]), "Yes" = unname(tls_colors["1a"])),
                    ordinal = c("No" = 0, "Yes" = 1)),
  pCR        = list(colors = pcr_colors,           ordinal = c("No" = 0, "Yes" = 1)),
  IGG        = list(colors = igg_level_colors,     ordinal = c("Low" = 1, "Medium" = 2, "High" = 3))
)

# --- Grouping variables: which metadata + which column holds each cluster/call ---
grouping_vars <- list(
  CellType      = list(meta = dwAnnotated@meta.data,      col = "celltype"),
  Cluster       = list(meta = dwAnnotated@meta.data,      col = "decontX_clusters"),
  PAM50         = list(meta = dwTumoral@meta.data,        col = "PAM50_predicted"),
  scSubtype     = list(meta = dwTumoral@meta.data,        col = "scSubtype_call")
)

association_summary <- list()
correlation_summary <- list()

for (grp_name in names(grouping_vars)) {
  g <- grouping_vars[[grp_name]]

  for (clin_name in names(clinical_vars)) {
    c_info <- clinical_vars[[clin_name]]

    # --- Composition + association test (contingency table, chi-square/Fisher, Cramer's V) ---
    assoc <- association_stats(g$meta, group_col = g$col, clinical_col = clin_name,
                               group_label = grp_name, clinical_label = clin_name,
                               results_path = results_GEMX_CLINICAL_path)
    association_summary[[paste0(grp_name, "_", clin_name)]] <- assoc$summary

    plot_transition_heatmap(assoc$table, x_col = "Clinical", y_col = "Group", count_col = "Freq",
                            title = paste0(grp_name, " composition by clinical ", clin_name),
                            xlab = clin_name, ylab = grp_name,
                            filename = paste0("Heatmap_", grp_name, "_x_", clin_name, ".png"),
                            results_path = results_GEMX_CLINICAL_path,
                            mode = "gradient", fill_mode = "pct", label_mode = "n_pct",
                            pct_group_col = "Group",
                            low_color = "white", high_color = "steelblue",
                            legend_name = paste0("% of\n", grp_name), text_size = 2.8,
                            width = 8, height = 7)

    # --- Proportion: stacked composition barplot (% of clin_name within each grouping level) ---
    plot_association_barplot(g$meta, group_col = g$col, clinical_col = clin_name,
                             results_path = results_GEMX_CLINICAL_path,
                             filename = paste0("Barplot_", grp_name, "_x_", clin_name, ".png"),
                             colors = c_info$colors, group_label = grp_name, clinical_label = clin_name)

    # --- Correlation: per-sample proportion vs ordinal clinical score (skipped for Subtype) ---
    if (!is.null(c_info$ordinal)) {
      cor_df <- proportion_correlation(g$meta, group_col = g$col, clinical_col = clin_name,
                                       clinical_levels = c_info$ordinal)
      cor_df$group_var <- grp_name
      cor_df$clinical_var <- clin_name
      correlation_summary[[paste0(grp_name, "_", clin_name)]] <- cor_df

      plot_correlation_barplot(cor_df, results_path = results_GEMX_CLINICAL_path,
                               filename = paste0("Barplot_Correlation_", grp_name, "_x_", clin_name, ".png"),
                               title = paste0(grp_name, " proportion vs ", clin_name, " (Spearman, per sample)"),
                               group_label = grp_name)
    }
  }
}

# --- Per-sample celltype proportion boxplots, faceted by celltype, one per clinical covariate ---
celltype_prop_df <- compute_celltype_proportions(dwAnnotated, sample_id_col = "orig.ident",
                                                 celltype_col = "celltype",
                                                 clinical_cols = names(clinical_vars))
# --- Per-sample tumor cluster proportion boxplots, faceted by cluster, one per clinical covariate ---
cluster_prop_df <- compute_celltype_proportions(dwTumoral, sample_id_col = "orig.ident",
                                                celltype_col = "celltype",
                                                clinical_cols = names(clinical_vars))

for (clin_name in names(clinical_vars)) {
  plot_celltype_proportions(celltype_prop_df, x_var = clin_name,
                            results_path = results_GEMX_CLINICAL_path,
                            filename = paste0("Boxplot_CelltypeProportion_by_", clin_name, ".png"),
                            celltype_col = "celltype",
                            colors = clinical_vars[[clin_name]]$colors)

  plot_celltype_proportions(cluster_prop_df, x_var = clin_name,
                            results_path = results_GEMX_CLINICAL_path,
                            filename = paste0("Boxplot_TumorCelltypeProportion_by_", clin_name, ".png"),
                            celltype_col = "celltype",
                            colors = clinical_vars[[clin_name]]$colors)
}
cat("\n Celltype and tumor cluster proportion boxplots by clinical covariate done \n")

association_summary_df <- do.call(rbind, association_summary)
association_summary_df$p_val_adj <- p.adjust(association_summary_df$p_val, method = "BH")
association_summary_df <- association_summary_df[order(association_summary_df$p_val_adj), ]
write.csv(association_summary_df, paste0(results_GEMX_CLINICAL_path, "Summary_ClinicalAssociation.csv"), row.names = FALSE)
cat("\n Clinical association (composition/chi-square/Cramer's V) summary done \n")

correlation_summary_df <- do.call(rbind, correlation_summary)
correlation_summary_df$p_val_adj <- p.adjust(correlation_summary_df$p_val, method = "BH")
correlation_summary_df <- correlation_summary_df[order(correlation_summary_df$p_val_adj), ]
write.csv(correlation_summary_df, paste0(results_GEMX_CLINICAL_path, "Summary_ProportionCorrelation.csv"), row.names = FALSE)
cat("\n Clinical proportion-correlation (Spearman) summary done \n")

cat(paste("\n ---- FINISHED CLINICAL COVARIATE ASSOCIATION ----
    Generated files (per grouping x clinical combination, in CellAnnotation/Clinical/):
        · CrossTab_<Group>_x_<Clinical>.csv
        · Heatmap_<Group>_x_<Clinical>.png
        · Barplot_<Group>_x_<Clinical>.png
        · Barplot_Correlation_<Group>_x_<Clinical>.png   (TLS/pCR/IGG only)
        · Boxplot_CelltypeProportion_by_<Clinical>.png   (per-sample proportion, faceted by celltype)
        · Boxplot_TumorClusterProportion_by_<Clinical>.png   (per-sample proportion, faceted by tumor cluster)
    Generated summaries:
        · Summary_ClinicalAssociation.csv    (chi-square/Fisher + Cramer's V, all combos)
        · Summary_ProportionCorrelation.csv  (Spearman rho, TLS/pCR/IGG combos)
    "))

