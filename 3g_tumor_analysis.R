##
##  Single Cell Analysis Step 4: Tumor cell analysis & annotation
##  Runs AFTER 3d_tumor.R — takes tumoral_scored.rds as input.
##  Computes differential pathway/module scores and correlation analysis,
##  produces all tumor-level summary plots, then performs manual cluster
##  annotation and propagates labels to the full object.
##  Workflow: run once to review plots → fill in clusters_tumor_annotated → re-run.
##

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

# Import plot functions
source(paste0(wd, "TCA_plots.R"))
source(paste0(wd, "CA_plots.R"))
source(paste0(wd, "CL_plots.R"))

# Load scored tumor data and full annotated object
dwTumoral <- readRDS(paste0(results_GEMX_TUMOR_path, "tumoral_scored.rds"))
cat("\n Tumor scored data loaded \n")
dwAnnotated <- readRDS(paste0(results_GEMX_CA_path, "Copykat/copykat_annotated_data.rds"))
cat("\n Full annotated data loaded \n")

# Set PAM50-like marker genes
pam50_genes <- list(
  "Her2+" = c("ERBB2", "GRB7", "BLVRA", "TMEM45B", "FGFR4"),
  Lum = c("ESR1", "PGR", "BAG1", "MAPT", "NAT1", "ZIP6"),
  TNBC = c("MKI67", "CCNE1", "ANLN", "CDC20", "EGFR", "MYC")
)

# --- HTML marker tables per cluster ---
plot_marker_tables(results_GEMX_TUMOR_path)
cat("\n Cluster marker tables saved \n")

# --- PAM50 dotplot ---
plot_marker_dotplot(dwTumoral, marker_groups = pam50_genes,
                    results_path = results_GEMX_TUMOR_path,
                    filename = "Dotplot_PAM50.png",
                    group_by = "seurat_clusters")

# --- Validation: predicted subtype vs clinical subtype ---
plot_dimplot(dwTumoral, reduction = "umap", group_by = "Subtype",
             results_path = results_GEMX_TUMOR_path, filename = "DimPlot_UMAP_Clinical.png")
plot_dimplot(dwTumoral, reduction = "umap", group_by = "PAM50_predicted",
             results_path = results_GEMX_TUMOR_path, filename = "DimPlot_UMAP_PAM50.png")
plot_dimplot(dwTumoral, reduction = "umap", group_by = "seurat_clusters",
             results_path = results_GEMX_TUMOR_path, filename = "DimPlot_UMAP_Clusters.png")
plot_dimplot(subset(dwTumoral, subset = dwTumoral$celltype_cont %in% c("Lum", "TNBC", "Her2+")),
             reduction = "umap", group_by = "celltype_cont",
             results_path = results_GEMX_TUMOR_path, filename = "DimPlot_UMAP_ContPAM50.png")

# --- PAM50 concordance heatmap ---
concordance_table <- table(Predicted = dwTumoral$PAM50_predicted, Clinical = dwTumoral$Subtype)
plot_pam50_concordance(concordance_table, results_GEMX_TUMOR_path)

# --- Cluster composition by sample and subtype ---
plot_cluster_composition(dwTumoral, cluster_col = "orig.ident", group_by = "PAM50_predicted",
                         filename = "BarPlot_PAM50Composition_bySample.png", results_path = results_GEMX_TUMOR_path)
plot_cluster_composition(dwTumoral, cluster_col = "orig.ident", group_by = "Subtype",
                         filename = "BarPlot_SubtypeComposition_bySample.png", results_path = results_GEMX_TUMOR_path)
plot_cluster_composition(dwTumoral, cluster_col = "decontX_clusters", group_by = "PAM50_predicted",
                         filename = "BarPlot_PAM50Composition_byCluster.png", results_path = results_GEMX_TUMOR_path)
plot_cluster_composition(dwTumoral, cluster_col = "decontX_clusters", group_by = "Subtype",
                         filename = "BarPlot_SubtypeComposition_byCluster.png", results_path = results_GEMX_TUMOR_path)

# --- PAM50 composition stacked bars, faceted by subtype ---
plot_pam50_composition_facet(dwTumoral, results_GEMX_TUMOR_path)

# --- PROGENy pathway activity dotplots ---
DefaultAssay(dwTumoral) <- "progeny"
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

# --- DecoupleR pathway activity dotplots ---
DefaultAssay(dwTumoral) <- "progeny_decoupler"
plot_marker_dotplot(dwTumoral, group_by = "Subtype",
                    marker_groups = as.list(setNames(rownames(dwTumoral[["progeny_decoupler"]]),
                                                     rownames(dwTumoral[["progeny_decoupler"]]))),
                    filename = "DotPlot_ProgenyDecoupleR_bySubtype.png", results_path = results_GEMX_TUMOR_path)
plot_marker_dotplot(dwTumoral, group_by = "PAM50_predicted",
                    marker_groups = as.list(setNames(rownames(dwTumoral[["progeny_decoupler"]]),
                                                     rownames(dwTumoral[["progeny_decoupler"]]))),
                    filename = "DotPlot_ProgenyDecoupleR_byPAM50.png", results_path = results_GEMX_TUMOR_path)
plot_marker_dotplot(dwTumoral, group_by = "decontX_clusters",
                    marker_groups = as.list(setNames(rownames(dwTumoral[["progeny_decoupler"]]),
                                                     rownames(dwTumoral[["progeny_decoupler"]]))),
                    filename = "DotPlot_ProgenyDecoupleR_byCluster.png", results_path = results_GEMX_TUMOR_path)
DefaultAssay(dwTumoral) <- "RNA_decontX"

# --- Cell cycle phase UMAP and boxplot ---
plot_dimplot(dwTumoral, reduction = "umap", group_by = "Phase",
             results_path = results_GEMX_TUMOR_path, filename = "DimPlot_CellCyclePhase.png")
plot_cellcycle_boxplot(dwTumoral, results_GEMX_TUMOR_path)

# Compute differential pathway activity per cluster (PROGENy)
DefaultAssay(dwTumoral) <- "progeny"
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

# Compute differential pathway activity per cluster (DecoupleR)
DefaultAssay(dwTumoral) <- "progeny_decoupler"
decoupler_progeny_diff <- FindAllMarkers(dwTumoral,
                                         assay = "progeny_decoupler",
                                         slot = "scale.data",
                                         group.by = "seurat_clusters",
                                         test.use = "wilcox",
                                         logfc.threshold = 0,
                                         min.pct = 0)
decoupler_progeny_diff <- decoupler_progeny_diff %>%
  dplyr::rename(pathway = gene, avg_diff = avg_log2FC) %>%
  dplyr::select(cluster, pathway, avg_diff, p_val, p_val_adj, pct.1, pct.2)
write.csv(decoupler_progeny_diff, paste0(results_GEMX_TUMOR_path, "Differential_PROGENyDecoupleR_byCluster.csv"), row.names = FALSE)
DefaultAssay(dwTumoral) <- "RNA_decontX"
cat("\n Differential DecoupleR-PROGENy pathways per cluster done \n")

# Compute differential module scores per cluster (PAM50, cell cycle, copykat)
differential_scores_by_cluster <- function(object, score_cols, group_by = "seurat_clusters") {
  meta <- object@meta.data
  clusters <- sort(unique(as.character(meta[[group_by]])))

  results <- do.call(rbind, lapply(clusters, function(cl) {
    do.call(rbind, lapply(score_cols, function(sc) {
      in_group  <- meta[[sc]][meta[[group_by]] == cl]
      out_group <- meta[[sc]][meta[[group_by]] != cl]
      wtest <- wilcox.test(in_group, out_group)
      data.frame(
        cluster  = cl,
        score    = sc,
        mean_in  = mean(in_group),
        mean_out = mean(out_group),
        avg_diff = mean(in_group) - mean(out_group),
        p_val    = wtest$p.value
      )
    }))
  }))

  results$p_val_adj <- p.adjust(results$p_val, method = "BH")
  results <- results[order(results$cluster, -abs(results$avg_diff)), ]
  return(results)
}

pam50_score_cols <- c("PAM50_Lum", "PAM50_TNBC", "PAM50_Her2+")
pam50_diff <- differential_scores_by_cluster(dwTumoral, pam50_score_cols)
write.csv(pam50_diff, paste0(results_GEMX_TUMOR_path, "Differential_PAM50_byCluster.csv"), row.names = FALSE)
cat("\n Differential PAM50 module scores per cluster done \n")

cellcycle_score_cols <- c("S.Score", "G2M.Score")
cellcycle_diff <- differential_scores_by_cluster(dwTumoral, cellcycle_score_cols)
write.csv(cellcycle_diff, paste0(results_GEMX_TUMOR_path, "Differential_CellCycle_byCluster.csv"), row.names = FALSE)
cat("\n Differential cell cycle scores per cluster done \n")

copykat_score_cols <- c("copykat_prediction", "copykat_cnas")
copykat_diff <- differential_scores_by_cluster(dwTumoral, copykat_score_cols)
write.csv(copykat_diff, paste0(results_GEMX_TUMOR_path, "Differential_CellCycle_byCNA.csv"), row.names = FALSE)
cat("\n Differential copykat scores per cluster done \n")

# --- Heatmaps of differential pathway activity ---
plot_progeny_diff_heatmap(progeny_diff, results_GEMX_TUMOR_path,
                          filename = "Heatmap_Differential_PROGENy_byCluster.png",
                          title = "Differential PROGENy pathway activity per cluster (vs all other clusters)")
plot_progeny_diff_heatmap(decoupler_progeny_diff, results_GEMX_TUMOR_path,
                          filename = "Heatmap_Differential_PROGENyDecoupleR_byCluster.png",
                          title = "Differential DecoupleR-PROGENy pathway activity per cluster (vs all other clusters)")

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
plot_differential_barplot(paste0(results_GEMX_TUMOR_path, "Differential_PROGENyDecoupleR_byCluster.csv"),
                          "pathway", "Differential DecoupleR-PROGENy pathway activity by cluster",
                          "Barplot_Differential_PROGENyDecoupleR.png", results_GEMX_TUMOR_path)
plot_differential_barplot(paste0(results_GEMX_TUMOR_path, "Differential_Copykat_byCluster.csv"),
                          "score", "Differential Copykat scores by cluster",
                          "Barplot_Differential_Copykat.png", results_GEMX_TUMOR_path)

# Correlation analysis across all continuous scores (Spearman)
pam50_scores <- dwTumoral@meta.data[, c("PAM50_Lum", "PAM50_TNBC", "PAM50_Her2+")]
colnames(pam50_scores) <- c("PAM50_Lum", "PAM50_TNBC", "PAM50_Her2")
cellcycle_scores <- dwTumoral@meta.data[, c("S.Score", "G2M.Score")]
progeny_scores <- t(as.matrix(GetAssayData(dwTumoral, assay = "progeny", layer = "data")))

all_scores <- cbind(pam50_scores, cellcycle_scores, progeny_scores)
all_scores <- all_scores[complete.cases(all_scores), ]
cat(paste0("\n Correlation analysis on ", nrow(all_scores), " cells, ", ncol(all_scores), " variables \n"))

cor_matrix <- cor(all_scores, method = "spearman")
plot_correlation_matrix(cor_matrix, results_GEMX_TUMOR_path)

# Extract unique pairwise correlations, sorted by absolute strength
cor_long <- as.data.frame(as.table(cor_matrix)) %>%
  rename(var1 = Var1, var2 = Var2, r = Freq) %>%
  filter(as.character(var1) < as.character(var2)) %>%
  arrange(desc(abs(r)))

top_correlations <- head(cor_long, 20)
cat("\n Top 20 strongest correlations: \n")
print(top_correlations)

# --- Scatter plots for top 10 pairwise correlations ---
for (i in 1:min(10, nrow(top_correlations))) {
  plot_correlation_scatter(all_scores,
                           as.character(top_correlations$var1[i]),
                           as.character(top_correlations$var2[i]),
                           results_GEMX_TUMOR_path)
}
cat("\n Files: Correlation_matrix_allScores.png\n")
cat(" Plots: top 10 Scatter_*.png (strongest global correlations)\n")

# --- Celltype transition plot (contaminated vs decontaminated) ---
plot_celltype_transition(dwTumoral, paste0(results_path, "GEMX/DecontX/"))

# Manual cluster annotation
# Review FindMarkers CSVs (from 3d_tumor.R) and the plots above before filling in.
# Re-run this script after filling in to generate the final annotated objects.
clusters_tumor_annotated <- c(
  "0" = "",  "1" = "",
  "4" = "",  "7" = "",
  "8" = "",  "9" = "",
  "11" = "", "17" = ""
)

dwTumoral$celltype <- factor(unname(clusters_tumor_annotated[as.character(dwTumoral$decontX_clusters)]))

# Label transfer to full annotated object
annotation_vec <- setNames(as.character(dwTumoral[["celltype"]][, 1]), colnames(dwTumoral))

if ("celltype" %in% colnames(dwAnnotated@meta.data)) {
  existing <- as.character(dwAnnotated[["celltype"]][, 1])
} else {
  existing <- rep(NA_character_, ncol(dwAnnotated))
}
names(existing) <- colnames(dwAnnotated)
existing[names(annotation_vec)] <- annotation_vec
dwAnnotated[["celltype"]] <- factor(existing)

# --- Annotated dimplot (full object) ---
plot_dimplot(dwAnnotated, reduction = "umap_decontX", group_by = "celltype", label = TRUE,
             results_path = results_GEMX_CA_path, filename = "DimPlot_UMAP_Annotated.png")

# Export final annotated data
saveRDS(dwTumoral, file.path(results_GEMX_TUMOR_path, "tumoral_annotated.rds"))
saveRDS(dwAnnotated, file.path(results_GEMX_CA_path, "fully_annotated_data.rds"))
cat("\n Tumor annotations propagated to global 'celltype' column \n")

cat(paste("\n ---- FINISHED TUMOR CELL ANALYSIS ----
    Generated files:
        · tumoral_annotated.rds
        · fully_annotated_data.rds
        · Differential_PROGENy_byCluster.csv
        · Differential_PROGENyDecoupleR_byCluster.csv
        · Differential_PAM50_byCluster.csv
        · Differential_CellCycle_byCluster.csv
        · Differential_CellCycle_byCNA.csv
        · Correlation_matrix_allScores.csv
    Generated plots:
        · DimPlot_UMAP_Annotated.png
        · Heatmap_Differential_PROGENy_byCluster.png
        · Heatmap_Differential_PROGENyDecoupleR_byCluster.png
        · Barplot_Differential_*.png
        · Scatter_*.png
        · Transition_celltype_beforeAfter_decontX.png
        "))
