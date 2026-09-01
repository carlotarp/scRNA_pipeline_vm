
# Import libraries
library("Seurat")
library(dplyr)
library(decoupleR)
library(progeny)
library(tidyr)
library(tibble)
library(ggplot2)
library(corrplot)
library(gt)

# Set paths
project_path <- "/home/usuario/PROJECTS/260724_victor_scRNA/"
wd <- paste0(project_path, "codes/scRNA_pipeline/")
setwd(wd)
results_path <- paste0(project_path, "results/")
results_GEMX_CA_path <- paste0(results_path, "GEMX/DecontX/CellAnnotation/")
results_GEMX_TUMOR_path <- paste0(results_GEMX_CA_path, "Tumor/")

# Load Plot Functions
source(paste0(wd, "CA_plots.R"))
source(paste0(wd, "CL_plots.R"))

# Load Tumor Data
dwTumoral <- readRDS(paste0(results_GEMX_TUMOR_path, "tumoral_annotated.rds"))
cat("\n Fully tumor data loaded \n")

# Set PAM50-like Marker Genes
pam50_genes <- list(
  "Her2+" = c("ERBB2","GRB7","BLVRA","TMEM45B", "FGFR4"),
  Lum = c("ESR1","PGR","BAG1","MAPT","NAT1","ZIP6"),
  TNBC = c("MKI67","CCNE1","ANLN","CDC20","EGFR","MYC")
)

marker_files <- list.files(results_GEMX_TUMOR_path, pattern = "^cluster_.*_FindMarkers\\.csv$", full.names = F)
for (f in marker_files) {
  df <- read.csv(paste0(results_GEMX_TUMOR_path, f))
  colnames(df)[1] <- "gene"
  cl_name <- gsub("cluster_(.*)_FindMarkers\\.csv", "\\1", basename(f))

  table_cl <- df %>%
    select(gene, avg_log2FC, pct.1, pct.2, p_val_adj) %>%
    arrange(desc(avg_log2FC)) %>%
    gt() %>%
    fmt_number(columns = c(avg_log2FC, pct.1, pct.2), decimals = 2) %>%
    fmt_scientific(columns = p_val_adj, decimals = 2) %>%
    tab_header(title = paste0("Top 20 markers - Cluster ", cl_name))

  gtsave(table_cl, paste0(results_GEMX_TUMOR_path, "Table_FindMarkers_cluster", cl_name, ".html"))
}

  # --- Dotplot by Subtype ---
plot_marker_dotplot(dwTumoral, marker_groups = pam50_genes,
                     results_path = results_GEMX_TUMOR_path,
                     filename = paste0("Dotplot_PAM50.png"),
                     group_by = "sueurat_clusters")

# --- Validation: compare predicted subtype vs clinical Subtype annotation ---
plot_dimplot(dwTumoral, reduction = "umap", group_by = "Subtype",
             results_path = results_GEMX_TUMOR_path, filename = "DimPlot_UMAP_Clinical.png")
plot_dimplot(dwTumoral, reduction = "umap", group_by = "PAM50_predicted",
             results_path = results_GEMX_TUMOR_path, filename = "DimPlot_UMAP_PAM50.png")
plot_dimplot(dwTumoral, reduction = "umap", group_by = "seurat_clusters",
             results_path = results_GEMX_TUMOR_path, filename = "DimPlot_UMAP_CLusters.png")
plot_dimplot(subset(dwTumoral, subset = dwTumoral$celltype_cont %in% c("Lum", "TNBC", "Her2+")), reduction = "umap", group_by = "celltype_cont",
             results_path = results_GEMX_TUMOR_path, filename = "DimPlot_UMAP_ContPAM50.png")

# --- Visualize Concordance HeatMap ---
concordance_table <- table(Predicted = dwTumoral$PAM50_predicted, Clinical = dwTumoral$Subtype)
concordance_df <- as.data.frame(concordance_table)
concordance_df <- concordance_df %>%
  group_by(Predicted) %>%
  mutate(pct = Freq / sum(Freq) * 100)

p_concordance <- ggplot(concordance_df, aes(x = Clinical, y = Predicted, fill = pct)) +
  geom_tile(color = "white") +
  geom_text(aes(label = paste0(Freq, "\n(", round(pct, 1), "%)")), size = 3.5) +
  scale_fill_gradient(low = "white", high = "steelblue", name = "% of\npredicted") +
  labs(title = "PAM50 predicted vs clinical Subtype concordance",
       x = "Clinical Subtype", y = "PAM50 predicted") +
  theme_minimal() +
  theme(panel.grid = element_blank())

ggsave(paste0(results_GEMX_TUMOR_path, "Heatmap_PAM50_vs_ClinicalSubtype.png"), p_concordance,
       width = 7, height = 5, dpi = 300, bg = "white")

# --- Visualize Cluster Composition ---
plot_cluster_composition(dwTumoral, cluster_col = "orig.ident", group_by = "PAM50_predicted",
                         filename = "BarPlot_PAM50Composition_bySample.png", results_path = results_GEMX_TUMOR_path)
plot_cluster_composition(dwTumoral, cluster_col = "orig.ident", group_by = "Subtype",
                         filename = "BarPlot_SubtypeComposition_bySample.png", results_path = results_GEMX_TUMOR_path)
plot_cluster_composition(dwTumoral, cluster_col = "decontX_clusters", group_by = "PAM50_predicted",
                         filename = "BarPlot_PAM50Composition_byCluster.png", results_path = results_GEMX_TUMOR_path)
plot_cluster_composition(dwTumoral, cluster_col = "decontX_clusters", group_by = "Subtype",
                         filename = "BarPlot_SubtypeComposition_byCluster.png", results_path = results_GEMX_TUMOR_path)

composition_df <- dwTumoral@meta.data %>%
  group_by(orig.ident, PAM50_predicted, Subtype) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(orig.ident, Subtype)

p_composition_facet <- ggplot(composition_df, aes(x = orig.ident, y = n, fill = PAM50_predicted)) +
  geom_col(position = "stack") +
  facet_wrap(~Subtype, scales = "free_x") +
  theme_bw() + theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "PAM50 predicted composition per cluster, by clinical subtype",
       x = "Sample", y = "nº of cells", fill = "PAM50 predicted")

ggsave(paste0(results_GEMX_TUMOR_path, "Composition_PAM50_bySample_facetSubtype.png"), p_composition_facet,
       width = 12, height = 6, dpi = 300, bg = "white")


# --- Pathway activity ---
DefaultAssay(dwTumoral) <- "progeny"
plot_marker_dotplot(dwTumoral, group_by = "Subtype", marker_groups = rownames(dwTumoral[["progeny"]]),
                    filename = "DotPlot_Progeny_bySubtype.png", results_path = results_GEMX_TUMOR_path)
plot_marker_dotplot(dwTumoral, group_by = "PAM50_predicted", marker_groups = rownames(dwTumoral[["progeny"]]),
                    filename = "DotPlot_Progeny_byPAM50.png", results_path = results_GEMX_TUMOR_path)
plot_marker_dotplot(dwTumoral, group_by = "decontX_clusters", marker_groups = rownames(dwTumoral[["progeny"]]),
                    filename = "DotPlot_Progeny_byCluster.png", results_path = results_GEMX_TUMOR_path)

DefaultAssay(dwTumoral) <- "progeny_decoupler"
plot_marker_dotplot(dwTumoral, group_by = "Subtype", marker_groups = rownames(dwTumoral[["progeny_decoupler"]]),
                    filename = "DotPlot_ProgenyDecoupleR_bySubtype.png", results_path = results_GEMX_TUMOR_path)
plot_marker_dotplot(dwTumoral, group_by = "PAM50_predicted", marker_groups = rownames(dwTumoral[["progeny_decoupler"]]),
                    filename = "DotPlot_ProgenyDecoupleR_byPAM50.png", results_path = results_GEMX_TUMOR_path)
plot_marker_dotplot(dwTumoral, group_by = "decontX_clusters", marker_groups = rownames(dwTumoral[["progeny_decoupler"]]),
                    filename = "DotPlot_ProgenyDecoupleR_byCluster.png", results_path = results_GEMX_TUMOR_path)
DefaultAssay(dwTumoral) <- "RNA_decontX"


# --- Cell cycle phase by predicted PAM50 subtype and by cluster ---
plot_dimplot(dwTumoral, reduction = "umap", group_by = "Phase",
             results_path = results_GEMX_TUMOR_path, filename = "DimPlot_CellCyclePhase.png")

cellcycle_long <- dwTumoral@meta.data %>%
  select(seurat_clusters, Subtype, S.Score, G2M.Score) %>%
  pivot_longer(cols = c(S.Score, G2M.Score), names_to = "score_type", values_to = "score")

p_cellcycle_box <- ggplot(cellcycle_long, aes(x = seurat_clusters, y = score, fill = score_type)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.3, position = position_dodge(width = 0.8)) +
  facet_wrap(~Subtype, scales = "free_x") +
  theme_bw() + theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Cell cycle score distribution by cluster, by clinical subtype",
       x = "Cluster", y = "Score", fill = "Score type")

ggsave(paste0(results_GEMX_TUMOR_path, "CellCycle_scores_boxplot_byCluster_facetSubtype.png"), p_cellcycle_box,
       width = 12, height = 6, dpi = 300, bg = "white")



# Diferential Cluster Analisys by PROGENy Pathways
DefaultAssay(dwTumoral) <- "progeny"
progeny_diff <- FindAllMarkers(dwTumoral,
                                 assay = "progeny",
                                 slot = "scale.data",
                                 group.by = "seurat_clusters",
                                 test.use = "wilcox",
                                 logfc.threshold = 0,   # keep all pathways, not just "big" ones
                                 min.pct = 0)

# Rename for clarity
progeny_diff <- progeny_diff %>%
  dplyr::rename(pathway = gene, avg_diff = avg_log2FC) %>%
  dplyr::select(cluster, pathway, avg_diff, p_val, p_val_adj, pct.1, pct.2)

write.csv(progeny_diff, paste0(results_GEMX_TUMOR_path, "Differential_PROGENy_byCluster.csv"), row.names = FALSE)
cat("\n Differential PROGENy pathways per cluster done \n")

# Diferential Cluster Analisys by Progeny-DecoupleR Pathways
DefaultAssay(dwTumoral) <- "progeny_decoupler"
decoupler_progeny_diff <- FindAllMarkers(dwTumoral,
                                 assay = "progeny_decoupler",
                                 slot = "scale.data",
                                 group.by = "seurat_clusters",
                                 test.use = "wilcox",
                                 logfc.threshold = 0,   # keep all pathways, not just "big" ones
                                 min.pct = 0)

# Rename for clarity
decoupler_progeny_diff <- decoupler_progeny_diff %>%
  dplyr::rename(pathway = gene, avg_diff = avg_log2FC) %>%
  dplyr::select(cluster, pathway, avg_diff, p_val, p_val_adj, pct.1, pct.2)

write.csv(decoupler_progeny_diff, paste0(results_GEMX_TUMOR_path, "Differential_PROGENyDecoupleR_byCluster.csv"), row.names = FALSE)
DefaultAssay(dwTumoral) <- "RNA_decontX"
cat("\n Differential PROGENy pathways per cluster done \n")

# Differential Score Function
differential_scores_by_cluster <- function(object, score_cols, group_by = "seurat_clusters") {
  meta <- object@meta.data
  clusters <- sort(unique(as.character(meta[[group_by]])))

  results <- do.call(rbind, lapply(clusters, function(cl) {
    do.call(rbind, lapply(score_cols, function(sc) {
      in_group <- meta[[sc]][meta[[group_by]] == cl]
      out_group <- meta[[sc]][meta[[group_by]] != cl]

      wtest <- wilcox.test(in_group, out_group)

      data.frame(
        cluster = cl,
        score = sc,
        mean_in = mean(in_group),
        mean_out = mean(out_group),
        avg_diff = mean(in_group) - mean(out_group),
        p_val = wtest$p.value
      )
    }))
  }))

  results$p_val_adj <- p.adjust(results$p_val, method = "BH")
  results <- results[order(results$cluster, -abs(results$avg_diff)), ]
  return(results)
}



# Diferential Cluster Analisys by PAM50 Score
DefaultAssay(dwTumoral) <- "RNA_decontX"
pam50_score_cols <- c("PAM50_Lum", "PAM50_TNBC", "PAM50_Her2+")
pam50_diff <- differential_scores_by_cluster(dwTumoral, pam50_score_cols)
write.csv(pam50_diff, paste0(results_GEMX_TUMOR_path, "Differential_PAM50_byCluster.csv"), row.names = FALSE)
cat("\n Differential PAM50 module scores per cluster done \n")

# Diferential Cluster Analisys by Cell Cycle Score
cellcycle_score_cols <- c("S.Score", "G2M.Score")
cellcycle_diff <- differential_scores_by_cluster(dwTumoral, cellcycle_score_cols)
write.csv(cellcycle_diff, paste0(results_GEMX_TUMOR_path, "Differential_CellCycle_byCluster.csv"), row.names = FALSE)
cat("\n Differential cell cycle scores per cluster done \n")

# Diferential Cluster Analisys by Cell Cycle Score
copykat_score_cols <- c("copykat_prediction", "copykat_cnas")
copykat_diff <- differential_scores_by_cluster(dwTumoral, copykat_score_cols)
write.csv(copykat_diff, paste0(results_GEMX_TUMOR_path, "Differential_CellCycle_byCNA.csv"), row.names = FALSE)
cat("\n Differential copykat scores per cluster done \n")

# --- Heatmap of PROGENy Differential scores by Cluster ---
p_progeny_diff <- ggplot(progeny_diff, aes(x = cluster, y = pathway, fill = avg_diff)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                        name = "Diff. activity\n(cluster vs rest)") +
  theme_minimal() +
  theme(panel.grid = element_blank(), plot.background = element_rect(fill = "white", color = NA)) +
  labs(title = "Differential PROGENy pathway activity per cluster (vs all other clusters)",
       x = "Cluster", y = NULL)

ggsave(paste0(results_GEMX_TUMOR_path, "Heatmap_Differential_PROGENy_byCluster.png"), p_progeny_diff,
       width = 8, height = 6, dpi = 300, bg = "white")

# --- Heatmap of DecoupleR-PROGENy Differential scores by Cluster ---
p_progeny_diff <- ggplot(decoupler_progeny_diff, aes(x = cluster, y = pathway, fill = avg_diff)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                        name = "Diff. activity\n(cluster vs rest)") +
  theme_minimal() +
  theme(panel.grid = element_blank(), plot.background = element_rect(fill = "white", color = NA)) +
  labs(title = "Differential DecoupleR-PROGENy pathway activity per cluster (vs all other clusters)",
       x = "Cluster", y = NULL)

ggsave(paste0(results_GEMX_TUMOR_path, "Heatmap_Differential_PROGENyDecoupleR_byCluster.png"), p_progeny_diff,
       width = 8, height = 6, dpi = 300, bg = "white")

plot_differential_barplot <- function(csv_path, feature_col, title, filename, results_path) {
  df <- read.csv(csv_path)

  p <- ggplot(df, aes(x = .data[[feature_col]], y = avg_diff, fill = avg_diff)) +
    geom_col() +
    facet_wrap(~cluster) +
    scale_fill_gradient2(low = "firebrick", mid = "white", high = "forestgreen", midpoint = 0,
                          name = "Diff. vs\nrest") +
    coord_flip() +
    theme_bw() + theme(panel.grid = element_blank()) +
    labs(title = title, x = NULL, y = "avg_diff (cluster vs rest)")

  ggsave(paste0(results_path, filename), p, width = 12, height = 8, dpi = 300, bg = "white")
}

plot_differential_barplot(paste0(results_GEMX_TUMOR_path, "Differential_PAM50_byCluster.csv"),
                            "score", "Differential PAM50 module scores by cluster",
                            "Barplot_Differential_PAM50.png", results_GEMX_TUMOR_path)

plot_differential_barplot(paste0(results_GEMX_TUMOR_path, "Differential_CellCycle_byCluster.csv"),
                            "score", "Differential cell cycle scores by cluster",
                            "Barplot_Differential_CellCycle.png", results_GEMX_TUMOR_path)

plot_differential_barplot(paste0(results_GEMX_TUMOR_path, "Differential_PROGENy_byCluster.csv"),
                            "pathway", "Differential PROGENy pathway activity by cluster",
                            "Barplot_Differential_PROGENy.png", results_GEMX_TUMOR_path)

plot_differential_barplot(paste0(results_GEMX_TUMOR_path, "Differential_PROGENyDecoupleR_byCluster.csv"),
                            "pathway", "Differential DecoupleR-PROGENy pathway activity by cluster",
                            "Barplot_Differential_PROGENyDecoupleR.png", results_GEMX_TUMOR_path)

plot_differential_barplot(paste0(results_GEMX_TUMOR_path, "Differential_Copykat_byCluster.csv"),
                            "score", "Differential Copykat pathway activity by cluster",
                            "Barplot_Differential_Copykat.png", results_GEMX_TUMOR_path)

##
##  Linear relationships & correlations across all continuous scores:
##

# =================================================================
# 1. Build a single data.frame with ALL continuous scores per cell
# =================================================================
pam50_scores <- dwTumoral@meta.data[, c("PAM50_Lum", "PAM50_TNBC", "PAM50_Her2+")]
colnames(pam50_scores) <- c("PAM50_Lum", "PAM50_TNBC", "PAM50_Her2")

cellcycle_scores <- dwTumoral@meta.data[, c("S.Score", "G2M.Score")]


progeny_scores <- t(as.matrix(GetAssayData(dwTumoral, assay = "progeny", layer = "data")))

all_scores <- cbind(pam50_scores, cellcycle_scores, progeny_scores) # infercna_scores
all_scores <- all_scores[complete.cases(all_scores), ]  # drop cells missing any layer

cat(paste0("\n Correlation analysis on ", nrow(all_scores), " cells, ", ncol(all_scores), " variables \n"))


# =================================================================
# 2. Correlation matrix - heatmap
# =================================================================
cor_matrix <- cor(all_scores, method = "spearman")

png(paste0(results_GEMX_TUMOR_path, "Correlation_matrix_allScores.png"),
    width = 3000, height = 3000, res = 300)
corrplot(cor_matrix, method = "color", type = "upper", order = "hclust",
          tl.col = "black", tl.srt = 45, addCoef.col = "black", number.cex = 0.5,
          title = "Correlation across all continuous scores (Spearman)", mar = c(0,0,2,0))
dev.off()


# =================================================================
# 3. Extract and rank the strongest correlations (as a data.frame)
# =================================================================
cor_long <- as.data.frame(as.table(cor_matrix)) %>%
  rename(var1 = Var1, var2 = Var2, r = Freq) %>%
  filter(as.character(var1) < as.character(var2)) %>%  # drop diagonal + duplicates
  arrange(desc(abs(r)))

top_correlations <- head(cor_long, 20)
cat("\n Top 20 strongest correlations: \n")
print(top_correlations)


# =================================================================
# 4. Scatter plots with linear fit, for the top correlated pairs
# =================================================================
plot_correlation_scatter <- function(data, var1, var2, results_path) {
  r_val <- cor(data[[var1]], data[[var2]], method = "spearman")
  p_val <- cor.test(data[[var1]], data[[var2]], method = "spearman")$p.value

  p <- ggplot(data, aes(x = .data[[var1]], y = .data[[var2]])) +
    geom_point(size = 0.4, alpha = 0.2, color = "steelblue") +
    geom_smooth(method = "lm", color = "firebrick", se = TRUE) +
    theme_bw() + theme(panel.grid = element_blank()) +
    labs(title = paste0(var1, " vs ", var2),
         subtitle = paste0("Spearman r = ", round(r_val, 3), ", p = ", signif(p_val, 3)),
         x = var1, y = var2)

  ggsave(paste0(results_path, "Scatter_", var1, "_vs_", var2, ".png"), p,
         width = 6, height = 5, dpi = 300, bg = "white")
}

for (i in 1:min(10, nrow(top_correlations))) {
  plot_correlation_scatter(all_scores, as.character(top_correlations$var1[i]),
                             as.character(top_correlations$var2[i]),
                             results_GEMX_TUMOR_path)
}


write.csv(cor_by_cluster, paste0(results_GEMX_TUMOR_path, "Correlation_byCluster.csv"), row.names = FALSE)

cat("\n ---- FINISHED correlation analysis ----\n")
cat(" Files: Correlation_matrix_allScores.png, Correlation_byCluster.csv\n")
cat(" Plots: top 10 Scatter_*.png (strongest global correlations)\n")



transition_df <- data.frame(
  cell = colnames(dwAnnotated),
  celltype_before = as.character(dwAnnotated$celltype_cont),
  celltype_after = as.character(dwAnnotated$celltype)
)

cat("Total células:", nrow(transition_df), "\n")
cat("NA en 'before':", sum(is.na(transition_df$celltype_before)), "\n")
cat("NA en 'after':", sum(is.na(transition_df$celltype_after)), "\n")



# =================================================================
# Tileplot de transición en % (normalizado por fila - "before")
# =================================================================

transition_df <- data.frame(
  cell = colnames(object),
  celltype_before = as.character(object$celltype_cont),
  celltype_after = as.character(object$celltype)
)

trans_counts <- transition_df %>%
  dplyr::count(celltype_before, celltype_after, name = "n") %>%
  dplyr::group_by(celltype_before) %>%
  dplyr::mutate(pct = round(n / sum(n) * 100, 1)) %>%
  dplyr::ungroup()

p_transition_pct <- ggplot(trans_counts, aes(x = celltype_after, y = celltype_before, fill = pct)) +
  geom_tile(color = "white") +
  geom_text(aes(label = pct), size = 2.8) +
  scale_fill_gradient(low = "white", high = "firebrick", name = "% of\nbefore") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.background = element_rect(fill = "white", color = NA)) +
  labs(title = "Celltype transition (%): before (celltype_cont) vs after (celltype, DecontX)",
       x = "Celltype - after DecontX", y = "Celltype - before DecontX")

ggsave(paste0(results_path, "/GEMX/DecontX/Transition_celltype_beforeAfter_decontX.png"), p_transition_pct,
       width = 10, height = 8, dpi = 300, bg = "white")