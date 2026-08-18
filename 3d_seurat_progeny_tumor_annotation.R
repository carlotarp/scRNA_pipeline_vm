##
##  Single Cell Analysis Step 3d: Tumor Subtype Annotation
##

# Import libraries
library("Seurat")
library(dplyr)
library(decoupleR)
library(progeny)
library(tidyr)
library(tibble)
library(ggplot2)

# Set paths
project_path <- "/home/usuario/PROJECTS/260724_victor_scRNA/"
wd <- paste0(project_path, "codes/scRNA_pipeline/")
setwd(wd)
results_path <- paste0(project_path, "results/")
results_GEMX_CA_path <- paste0(results_path, "GEMX/CellAnnotation/7500/")
results_GEMX_TUMOR_path <- paste0(results_GEMX_CA_path, "Tumor/")


# Load Fully Annotated Data
dwAnnotated <- readRDS(paste0(results_GEMX_CA_path, "notumor_clean_annotated_data.rds"))
cat("\n Fully annotated data loaded \n")

dwAnnotated <- readRDS(file.path(results_GEMX_CA_path, "fully_annotated_data.rds"))

# Subset Tumors
dwTumoral <- subset(dwAnnotated, subset = lineage == "Tumor")
dwTumoral <- JoinLayers(dwTumoral)

# Set PAM50-like Marker Genes
pam50_genes <- list(
  "Her2+" = c("ERBB2","GRB7","BLVRA","TMEM45B", "FGFR4"),
  Lum = c("ESR1","PGR","BAG1","MAPT","NAT1","ZIP6"),
  TNBC = c("MKI67","CCNE1","ANLN","CDC20","EGFR","MYC")
)

  # --- Dotplot by Subtype ---
plot_marker_dotplot(dwTumoral, marker_groups = pam50_genes,
                     results_path = results_GEMX_TUMOR_path,
                     filename = paste0("Dotplot_PAM50.png"),
                     group_by = "seurat_clusters")

# Compute PAM50 Module Scores
dwTumoral <- AddModuleScore(dwTumoral, features = pam50_genes, name = "PAM50_")

# AddModuleScore names columns PAM50_1, PAM50_2, PAM50_3 in list order - rename
# to the actual subtype names for clarity
pam50_score_cols <- paste0("PAM50_", seq_along(pam50_genes))
pam50_named_cols <- paste0("PAM50_", names(pam50_genes))
colnames(dwTumoral@meta.data)[match(pam50_score_cols, colnames(dwTumoral@meta.data))] <- pam50_named_cols

# Predicted subtype per cell = highest-scoring PAM50 module
pam50_scores <- as.matrix(dwTumoral@meta.data[, pam50_named_cols])
dwTumoral$PAM50_predicted <- names(pam50_genes)[apply(pam50_scores, 1, which.max)]
dwTumoral$PAM50_predicted <- factor(dwTumoral$PAM50_predicted, levels = names(pam50_genes))

cat("\n PAM50 module scores computed \n")
print(table(dwTumoral$PAM50_predicted))

# --- Validation: compare predicted subtype vs clinical Subtype annotation ---
plot_dimplot(dwTumoral, reduction = "umap", group_by = "Subtype",
             results_path = results_GEMX_TUMOR_path, filename = "DimPlot_UMAP_Clinical.png")
plot_dimplot(dwTumoral, reduction = "umap", group_by = "PAM50_predicted",
             results_path = results_GEMX_TUMOR_path, filename = "DimPlot_UMAP_PAM50.png")
plot_dimplot(dwTumoral, reduction = "umap", group_by = "seurat_clusters",
             results_path = results_GEMX_TUMOR_path, filename = "DimPlot_UMAP_CLusters.png")

concordance_table <- table(Predicted = dwTumoral$PAM50_predicted, Clinical = dwTumoral$Subtype)
write.csv(as.data.frame.matrix(concordance_table),
          paste0(results_GEMX_TUMOR_path, "PAM50_vs_Clinical_concordance.csv"))
cat("\n PAM50 predicted vs clinical Subtype concordance table saved \n")

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

plot_cluster_composition(dwTumoral, cluster_col = "orig.ident", group_by = "PAM50_predicted",
                         filename = "BarPlot_PAM50Composition_bySample.png", results_path = results_GEMX_TUMOR_path)
plot_cluster_composition(dwTumoral, cluster_col = "orig.ident", group_by = "Subtype",
                         filename = "BarPlot_SubtypeComposition_bySample.png", results_path = results_GEMX_TUMOR_path)
plot_cluster_composition(dwTumoral, group_by = "PAM50_predicted",
                         filename = "BarPlot_PAM50Composition_byCluster.png", results_path = results_GEMX_TUMOR_path)
plot_cluster_composition(dwTumoral, group_by = "Subtype",
                         filename = "BarPlot_SubtypeComposition_byCluster.png", results_path = results_GEMX_TUMOR_path)

# Compute PROGENy scores manually - pass a plain matrix (genes x cells) using
# the Seurat v5 `layer` syntax, bypassing progeny's broken .Seurat method
# (which still uses the removed `slot` argument internally)
expr_mat <- as.matrix(GetAssayData(dwTumoral, assay = "RNA", layer = "data"))
 
progeny_scores <- progeny::progeny(expr_mat, scale = TRUE, organism = "Human", top = 500, perm = 1)
 
# progeny()'s matrix method can return either (cells x pathways) or
# (pathways x cells) depending on version - detect orientation automatically
# instead of assuming, and transpose so we end up with pathways x cells
# (what CreateAssayObject expects: features in rows, cells in columns)
if (nrow(progeny_scores) == ncol(dwTumoral)) {
  progeny_mat <- t(progeny_scores)   # was cells x pathways -> transpose
} else {
  progeny_mat <- progeny_scores      # already pathways x cells
}
 
cat(paste0("\n PROGENy matrix dimensions: ", nrow(progeny_mat), " pathways x ", ncol(progeny_mat), " cells \n"))
 
dwTumoral[["progeny"]] <- CreateAssayObject(data = progeny_mat)
 
# Scale the progeny assay (standard step before plotting - puts pathways on
# a comparable, centered scale like z-scores)
dwTumoral <- ScaleData(dwTumoral, assay = "progeny")
cat("\n PROGENy pathway activities computed and added as 'progeny' assay \n")
DefaultAssay(dwTumoral) <- "progeny"

# --- Pathway activity ---
plot_marker_dotplot(dwTumoral, group_by = "Subtype", marker_groups = rownames(dwTumoral[["progeny"]]),
                    filename = "DotPlot_Progeny_bySubtype.png", results_path = results_GEMX_TUMOR_path)
plot_marker_dotplot(dwTumoral, group_by = "PAM50_predicted", marker_groups = rownames(dwTumoral[["progeny"]]),
                    filename = "DotPlot_Progeny_byPAM50.png", results_path = results_GEMX_TUMOR_path)
plot_marker_dotplot(dwTumoral, group_by = "seurat_clusters", marker_groups = rownames(dwTumoral[["progeny"]]),
                    filename = "DotPlot_Progeny_byCluster.png", results_path = results_GEMX_TUMOR_path)
DefaultAssay(dwTumoral) <- "RNA"
 

# Compute Cell Cycle Scores (S.Score, G2M.Score, Phase per cell)
dwTumoral <- CellCycleScoring(dwTumoral,
                                s.features = cc.genes.updated.2019$s.genes,
                                g2m.features = cc.genes.updated.2019$g2m.genes)
cat("\n Cell cycle scores computed \n")
print(table(dwTumoral$Phase))

# --- Cell cycle phase by predicted PAM50 subtype and by cluster ---
plot_dimplot(dwTumoral, reduction = "umap", group_by = "Phase",
             results_path = results_GEMX_TUMOR_path, filename = "DimPlot_CellCyclePhase.png")

phase_by_pam50 <- table(PAM50 = dwTumoral$PAM50_predicted, Phase = dwTumoral$Phase)
write.csv(as.data.frame.matrix(phase_by_pam50),
          paste0(results_GEMX_TUMOR_path, "CellCyclePhase_by_PAM50.csv"))

phase_by_cluster <- table(Cluster = dwTumoral$seurat_clusters, Phase = dwTumoral$Phase)
write.csv(as.data.frame.matrix(phase_by_cluster),
          paste0(results_GEMX_TUMOR_path, "CellCyclePhase_by_cluster.csv"))
cat("\n Cell cycle phase tables (by PAM50 and by cluster) saved \n")



##
##  Single Cell Analysis Step 6: Differential Scores per Cluster
##  Spatial-transcriptomics-style "domain marker" analysis, applied to
##  continuous scores (PROGENy pathways, PAM50 module scores, cell cycle
##  scores) instead of genes: for each cluster, how much MORE is each score
##  representing that cluster vs all other clusters combined?
##


DefaultAssay(dwTumoral) <- "progeny"

progeny_diff <- FindAllMarkers(dwTumoral,
                                 assay = "progeny",
                                 slot = "scale.data",
                                 group.by = "seurat_clusters",
                                 test.use = "wilcox",
                                 logfc.threshold = 0,   # keep all pathways, not just "big" ones
                                 min.pct = 0)

# Rename for clarity - "avg_log2FC" here is really "avg difference in pathway
# activity" (scores aren't logFC-like counts, but Seurat's column name sticks)
progeny_diff <- progeny_diff %>%
  rename(pathway = gene, avg_diff = avg_log2FC) %>%
  select(cluster, pathway, avg_diff, p_val, p_val_adj, pct.1, pct.2)

write.csv(progeny_diff, paste0(results_GEMX_TUMOR_path, "Differential_PROGENy_byCluster.csv"), row.names = FALSE)

DefaultAssay(dwTumoral) <- "RNA"
cat("\n [A] Differential PROGENy pathways per cluster done \n")


# =================================================================
# B. Generic differential score function - for continuous meta.data columns
#    that are NOT stored as a Seurat assay (PAM50 module scores, cell cycle
#    scores). Same logic as FindMarkers (one cluster vs the rest), applied
#    manually since these live in meta.data, not an assay.
# =================================================================
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


# =================================================================
# C. PAM50 module scores - differential per cluster
# =================================================================
DefaultAssay(dwTumoral) <- "RNA"
pam50_score_cols <- c("PAM50_Lum", "PAM50_TNBC", "PAM50_Her2+")
pam50_diff <- differential_scores_by_cluster(dwTumoral, pam50_score_cols)
write.csv(pam50_diff, paste0(results_GEMX_TUMOR_path, "Differential_PAM50_byCluster.csv"), row.names = FALSE)
cat("\n [C] Differential PAM50 module scores per cluster done \n")


# =================================================================
# D. Cell cycle scores - differential per cluster
# =================================================================
cellcycle_score_cols <- c("S.Score", "G2M.Score")
cellcycle_diff <- differential_scores_by_cluster(dwTumoral, cellcycle_score_cols)
write.csv(cellcycle_diff, paste0(results_GEMX_TUMOR_path, "Differential_CellCycle_byCluster.csv"), row.names = FALSE)
cat("\n [D] Differential cell cycle scores per cluster done \n")


# =================================================================
# E. Combined "dominant layer" summary per cluster - the top-scoring,
#    most significant feature from EACH layer, one row per cluster
# =================================================================
top_per_cluster <- function(diff_df, feature_col, layer_name) {
  diff_df %>%
    filter(p_val_adj < 0.05) %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_diff, n = 1) %>%
    ungroup() %>%
    transmute(cluster, layer = layer_name, dominant_feature = .data[[feature_col]],
              avg_diff = round(avg_diff, 3), p_val_adj = signif(p_val_adj, 3))
}

cluster_summary <- bind_rows(
  top_per_cluster(progeny_diff, "pathway", "PROGENy"),
  top_per_cluster(pam50_diff, "score", "PAM50"),
  top_per_cluster(cellcycle_diff, "score", "CellCycle")
) %>%
  arrange(cluster, layer)

write.csv(cluster_summary, paste0(results_GEMX_TUMOR_path, "ClusterAnnotation_DominantFeaturePerLayer.csv"), row.names = FALSE)
cat("\n [E] Combined dominant-feature-per-layer summary saved \n")
print(cluster_summary)


# =================================================================
# F. Heatmap of PROGENy differential scores by cluster (the main plot for
#    this kind of "domain annotation" analysis)
# =================================================================
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


#############

dwTumoral$celltype <- factor(dwTumoral$PAM50_predicted)
annotation_vec <- setNames(as.character(dwTumoral[["celltype"]][, 1]),
                            colnames(dwTumoral))

if ("celltype" %in% colnames(dwAnnotated@meta.data)) {
  existing <- as.character(dwAnnotated[["celltype"]][, 1])
} else {
  existing <- rep(NA_character_, ncol(dwAnnotated))
}
names(existing) <- colnames(dwAnnotated)

# Only overwrite the cells present in the subset - everything else (e.g.
# already-annotated Stromal cells) stays exactly as it was
existing[names(annotation_vec)] <- annotation_vec

dwAnnotated[["celltype"]] <- factor(existing)

# ---  Visualize Annotated Dimplot ---
plot_dimplot(dwAnnotated, reduction = "umap", group_by = "celltype", label = T,
             results_path = results_GEMX_CA_path, filename = "DimPlot_UMAP_Annotated.png")

 
saveRDS(dwTumoral, file.path(results_GEMX_TUMOR_path, "tumoral_annotated.rds"))
saveRDS(dwAnnotated, file.path(results_GEMX_CA_path, "fully_annotated_data.rds"))
cat("\n Tumoral subtypes propagated to global 'celltype' column \n")


cat(paste("\n ---- FINISHED TUMOR SUBTYPE ANNOTATION ----
    Generated files:
      · Tumoral_annotated.rds
      · fully_annotated_data.rds
      · PAM50_vs_ClinicalSubtype_concordance.csv
    Generated plots:
      · BarPlot_PAM50_(groupedby).png
      · DimPlot_UMAP_(groupedby).png
      · Dotplot_PROGENy_(groupedby).png
      · HeatMap_PAM50_vs_ClinicalSubtype.png
          "))


