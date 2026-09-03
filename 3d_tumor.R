##
##  Single Cell Analysis Step 3d: Tumor scoring
##  Runs AFTER 3a_lineage.R and 3e_copykat.R.
##  Subsets and reclusters tumor cells, computes PAM50 / PROGENy / DecoupleR
##  pathway scores and cell cycle phase, and finds cluster markers.
##  Manual annotation is done in 4_tumor_analysis.R after reviewing the plots.
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
results_GEMX_CA_path <- paste0(results_path, "GEMX/DecontX/CellAnnotation/")
results_GEMX_TUMOR_path <- paste0(results_GEMX_CA_path, "Tumor/")

# Load plot functions and shared utilities
source(paste0(wd, "CA_plots.R"))
source(paste0(wd, "CL_plots.R"))
source(paste0(wd, "utils.R"))

# Load fully annotated data (output of 3e_copykat.R)
dwAnnotated <- readRDS(paste0(results_GEMX_CA_path, "Copykat/copykat_annotated_data.rds"))
cat("\n Fully annotated data loaded \n")

# Subset and recluster tumor cells
dwTumoral <- generate_lineage_subset(dwAnnotated, lineage_name = "Tumor", resolution = 0.3)
plot_dimplot(dwTumoral, reduction = "umap", group_by = "seurat_clusters",
             results_path = results_GEMX_TUMOR_path,
             filename = "DimPlot_UMAP_Tumor.png", label = TRUE)

# Set PAM50-like marker genes
pam50_genes <- list(
  "Her2+" = c("ERBB2", "GRB7", "BLVRA", "TMEM45B", "FGFR4"),
  Lum     = c("ESR1", "PGR", "BAG1", "MAPT", "NAT1", "ZIP6"),
  TNBC    = c("MKI67", "CCNE1", "ANLN", "CDC20", "EGFR", "MYC")
)

# Compute PAM50 module scores
dwTumoral <- AddModuleScore(dwTumoral, features = pam50_genes, name = "PAM50_")
pam50_score_cols <- paste0("PAM50_", seq_along(pam50_genes))
pam50_named_cols <- paste0("PAM50_", names(pam50_genes))
colnames(dwTumoral@meta.data)[match(pam50_score_cols, colnames(dwTumoral@meta.data))] <- pam50_named_cols

# Predicted subtype = highest-scoring PAM50 module
pam50_scores <- as.matrix(dwTumoral@meta.data[, pam50_named_cols])
dwTumoral$PAM50_predicted <- names(pam50_genes)[apply(pam50_scores, 1, which.max)]
dwTumoral$PAM50_predicted <- factor(dwTumoral$PAM50_predicted, levels = names(pam50_genes))
cat("\n PAM50 module scores computed \n")
print(table(dwTumoral$PAM50_predicted))

concordance_table <- table(Predicted = dwTumoral$PAM50_predicted, Clinical = dwTumoral$Subtype)
write.csv(as.data.frame.matrix(concordance_table),
          paste0(results_GEMX_TUMOR_path, "PAM50_vs_Clinical_concordance.csv"))
cat("\n PAM50 vs clinical Subtype concordance table saved \n")

# Compute PROGENy pathway activity (progeny package)
expr_mat <- as.matrix(GetAssayData(dwTumoral, assay = "RNA_decontX", layer = "data"))
progeny_scores <- progeny::progeny(expr_mat, scale = TRUE, organism = "Human", top = 500, perm = 1)
dwTumoral[["progeny"]] <- CreateAssayObject(data = t(progeny_scores))
dwTumoral <- ScaleData(dwTumoral, assay = "progeny")
cat("\n PROGENy pathway activities added as 'progeny' assay \n")

# Compute PROGENy pathway activity (DecoupleR)
net_progeny <- get_progeny(organism = "human", top = 500)
progeny_acts <- run_mlm(
  mat = expr_mat, network = net_progeny,
  .source = "source", .target = "target", .mor = "weight",
  minsize = 5
)
progeny_mat_decoupler <- progeny_acts %>%
  pivot_wider(id_cols = source, names_from = condition, values_from = score) %>%
  column_to_rownames("source") %>%
  as.matrix()
dwTumoral[["progeny_decoupler"]] <- CreateAssayObject(data = progeny_mat_decoupler)
dwTumoral <- ScaleData(dwTumoral, assay = "progeny_decoupler")
DefaultAssay(dwTumoral) <- "RNA_decontX"
cat("\n DecoupleR-PROGENy pathway activities added as 'progeny_decoupler' assay \n")

# Compute cell cycle scores
dwTumoral <- CellCycleScoring(dwTumoral,
                               s.features = cc.genes.updated.2019$s.genes,
                               g2m.features = cc.genes.updated.2019$g2m.genes)
cat("\n Cell cycle scores computed \n")
print(table(dwTumoral$Phase))

phase_by_pam50 <- table(PAM50 = dwTumoral$PAM50_predicted, Phase = dwTumoral$Phase)
write.csv(as.data.frame.matrix(phase_by_pam50),
          paste0(results_GEMX_TUMOR_path, "CellCyclePhase_by_PAM50.csv"))
phase_by_cluster <- table(Cluster = dwTumoral$decontX_clusters, Phase = dwTumoral$Phase)
write.csv(as.data.frame.matrix(phase_by_cluster),
          paste0(results_GEMX_TUMOR_path, "CellCyclePhase_by_cluster.csv"))
cat("\n Cell cycle phase tables saved \n")

# Find markers of ambiguous clusters
clusters_to_check <- sort(unique(dwTumoral$seurat_clusters))
find_markers_for_clusters(dwTumoral, clusters_to_check, results_GEMX_TUMOR_path)
cat("\n Read cluster CSVs and review 4_tumor_analysis.R plots before filling in annotation \n")

# Export scored data (annotation moved to 4_tumor_analysis.R)
saveRDS(dwTumoral, file.path(results_GEMX_TUMOR_path, "tumoral_scored.rds"))
cat("\n Scored tumor data saved \n")

cat(paste("\n ---- FINISHED TUMOR SCORING ----
    Run 4_tumor_analysis.R next to review cluster plots, then fill in the annotation block.
    Generated files:
        · tumoral_scored.rds  (includes 'progeny' and 'progeny_decoupler' assays)
        · PAM50_vs_Clinical_concordance.csv
        · CellCyclePhase_by_PAM50.csv
        · CellCyclePhase_by_cluster.csv
        · cluster_(cluster)_FindMarkers.csv
    Generated plots:
        · DimPlot_UMAP_Tumor.png
        "))
