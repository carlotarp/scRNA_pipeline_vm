##
##  Single Cell Analysis Step 5: Tumor Subtype Annotation
##  PAM50-like module scores (AddModuleScore) + PROGENy pathway activity (decoupleR)
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
dir.create(results_GEMX_TUMOR_path, recursive = TRUE, showWarnings = FALSE)

# Import Plot Functions
source(paste0(wd, "CA_plots.R"))

# Load Fully Annotated Data
dwAnnotated <- readRDS(paste0(results_GEMX_CA_path, "notumor_annotated_data.rds"))
cat("\n Fully annotated data loaded \n")

# Subset Tumoral Cells
dwTumoral <- subset(dwAnnotated, subset = celltype == "Tumor")
dwTumoral <- JoinLayers(dwTumoral)
cat(paste0("\n Tumoral subset: ", ncol(dwTumoral), " cells \n"))


# Set PAM50-like Marker Genes
pam50_genes <- list(
      "Her2+"  = c("ERBB2", "GRB7", "BLVRA", "TMEM45B"),
      Lum   = c("ESR1", "PGR", "BAG1", "MAPT", "NAT1", "ZIP6"),
      TNBC = c("MKI67", "CCNE1", "ANLN", "CDC20", "EGFR", "MYC")
)


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
plot_dimplot(dwTumoral, reduction = "umap", group_by = "PAM50_predicted",
             results_path = results_GEMX_TUMOR_path, filename = "DimPlot_PAM50_predicted.png")

if ("Subtype" %in% colnames(dwTumoral@meta.data)) {
  plot_dimplot(dwTumoral, reduction = "umap", group_by = "Subtype",
               results_path = results_GEMX_TUMOR_path, filename = "DimPlot_ClinicalSubtype.png")

  concordance_table <- table(Predicted = dwTumoral$PAM50_predicted, Clinical = dwTumoral$Subtype)
  write.csv(as.data.frame.matrix(concordance_table),
            paste0(results_GEMX_TUMOR_path, "PAM50_vs_ClinicalSubtype_concordance.csv"))
  cat("\n PAM50 predicted vs clinical Subtype concordance table saved \n")
}

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
 
# --- Pathway activity by predicted PAM50 subtype ---
DefaultAssay(dwTumoral) <- "progeny"
progeny_dotplot <- DotPlot(dwTumoral, features = rownames(dwTumoral[["progeny"]]),
                            group.by = "PAM50_predicted") +
  RotatedAxis() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(title = "PROGENy pathway activity by predicted PAM50 subtype")
 
ggsave(paste0(results_GEMX_TUMOR_path, "Dotplot_PROGENy_byPAM50.png"), progeny_dotplot,
       width = 10, height = 6, dpi = 300, bg = "white")
DefaultAssay(dwTumoral) <- "RNA"
 

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
      · DimPlot_PAM50_predicted.png
      · DimPlot_ClinicalSubtype.png
      · Dotplot_PROGENy_byPAM50.png
          "))