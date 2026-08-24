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

# Load Plot Functions
source(paste0(wd, "CA_plots.R"))
source(paste0(wd, "CL_plots.R"))

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

# Write Table of Interest
concordance_table <- table(Predicted = dwTumoral$PAM50_predicted, Clinical = dwTumoral$Subtype)
write.csv(as.data.frame.matrix(concordance_table),
          paste0(results_GEMX_TUMOR_path, "PAM50_vs_Clinical_concordance.csv"))
cat("\n PAM50 predicted vs clinical Subtype concordance table saved \n")

# --- Visualize Concordance HeatMap ---
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
plot_cluster_composition(dwTumoral, group_by = "PAM50_predicted",
                         filename = "BarPlot_PAM50Composition_byCluster.png", results_path = results_GEMX_TUMOR_path)
plot_cluster_composition(dwTumoral, group_by = "Subtype",
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

# Compute PROGENy
expr_mat <- as.matrix(GetAssayData(dwTumoral, assay = "RNA", layer = "data"))
progeny_scores <- progeny::progeny(expr_mat, scale = TRUE, organism = "Human", top = 500, perm = 1)
dwTumoral[["progeny"]] <- CreateAssayObject(data = progeny_mat)
dwTumoral <- ScaleData(dwTumoral, assay = "progeny")
cat("\n PROGENy pathway activities computed and added as 'progeny' assay \n")

# --- Pathway activity ---
DefaultAssay(dwTumoral) <- "progeny"
plot_marker_dotplot(dwTumoral, group_by = "Subtype", marker_groups = rownames(dwTumoral[["progeny"]]),
                    filename = "DotPlot_Progeny_bySubtype.png", results_path = results_GEMX_TUMOR_path)
plot_marker_dotplot(dwTumoral, group_by = "PAM50_predicted", marker_groups = rownames(dwTumoral[["progeny"]]),
                    filename = "DotPlot_Progeny_byPAM50.png", results_path = results_GEMX_TUMOR_path)
plot_marker_dotplot(dwTumoral, group_by = "seurat_clusters", marker_groups = rownames(dwTumoral[["progeny"]]),
                    filename = "DotPlot_Progeny_byCluster.png", results_path = results_GEMX_TUMOR_path)
DefaultAssay(dwTumoral) <- "RNA"
 

# Compute Cell Cycle Scores
dwTumoral <- CellCycleScoring(dwTumoral,
                                s.features = cc.genes.updated.2019$s.genes,
                                g2m.features = cc.genes.updated.2019$g2m.genes)
cat("\n Cell cycle scores computed \n")
print(table(dwTumoral$Phase))

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

# Write Tables of Interest
phase_by_pam50 <- table(PAM50 = dwTumoral$PAM50_predicted, Phase = dwTumoral$Phase)
write.csv(as.data.frame.matrix(phase_by_pam50),
          paste0(results_GEMX_TUMOR_path, "CellCyclePhase_by_PAM50.csv"))

phase_by_cluster <- table(Cluster = dwTumoral$seurat_clusters, Phase = dwTumoral$Phase)
write.csv(as.data.frame.matrix(phase_by_cluster),
          paste0(results_GEMX_TUMOR_path, "CellCyclePhase_by_cluster.csv"))
cat("\n Cell cycle phase tables (by PAM50 and by cluster) saved \n")


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
  rename(pathway = gene, avg_diff = avg_log2FC) %>%
  select(cluster, pathway, avg_diff, p_val, p_val_adj, pct.1, pct.2)

write.csv(progeny_diff, paste0(results_GEMX_TUMOR_path, "Differential_PROGENy_byCluster.csv"), row.names = FALSE)
DefaultAssay(dwTumoral) <- "RNA"
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
DefaultAssay(dwTumoral) <- "RNA"
pam50_score_cols <- c("PAM50_Lum", "PAM50_TNBC", "PAM50_Her2+")
pam50_diff <- differential_scores_by_cluster(dwTumoral, pam50_score_cols)
write.csv(pam50_diff, paste0(results_GEMX_TUMOR_path, "Differential_PAM50_byCluster.csv"), row.names = FALSE)
cat("\n Differential PAM50 module scores per cluster done \n")



# Diferential Cluster Analisys by Cell Cycle Score
cellcycle_score_cols <- c("S.Score", "G2M.Score")
cellcycle_diff <- differential_scores_by_cluster(dwTumoral, cellcycle_score_cols)
write.csv(cellcycle_diff, paste0(results_GEMX_TUMOR_path, "Differential_CellCycle_byCluster.csv"), row.names = FALSE)
cat("\n Differential cell cycle scores per cluster done \n")

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


# Manual Cluster Annotation
clusters_tumor_annotated <- c( #########  7500  #########
  "0" = "",  "1" = "",
  "4" = "", "7" = "",
  "8" = "", "9" = "",
  "11" = "", "17" = ""
)

dwTumoral$celltype <- unname(clusters_tumor_annotated[as.character(dwTumoral$seurat_clusters)])
dwTumoral$celltype <- factor(dwTumoral$celltype)
dwTumoral$celltype <- factor(dwTumoral$PAM50_predicted)

# Label Transfer
annotation_vec <- setNames(as.character(dwTumoral[["celltype"]][, 1]),
                            colnames(dwTumoral))

if ("celltype" %in% colnames(dwAnnotated@meta.data)) {
  existing <- as.character(dwAnnotated[["celltype"]][, 1])
} else {
  existing <- rep(NA_character_, ncol(dwAnnotated))
}
names(existing) <- colnames(dwAnnotated)
existing[names(annotation_vec)] <- annotation_vec
dwAnnotated[["celltype"]] <- factor(existing)

# ---  Visualize Annotated Dimplot ---
plot_dimplot(dwAnnotated, reduction = "umap", group_by = "celltype", label = T,
             results_path = results_GEMX_CA_path, filename = "DimPlot_UMAP_Annotated.png")

# Export Annoatated Data
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

