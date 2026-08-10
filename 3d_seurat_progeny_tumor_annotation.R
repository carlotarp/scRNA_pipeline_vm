##
##  Single Cell Analysis Step 5: Tumor Subtype Annotation
##  PAM50-like module scores (AddModuleScore) + PROGENy pathway activity (decoupleR)
##

# Import libraries
library("Seurat")
library(dplyr)
library(decoupleR)
library(tidyr)
library(tibble)

# Set paths
project_path <- "/home/usuario/PROJECTS/260724_victor_scRNA/"
wd <- paste0(project_path, "codes/scRNA_pipeline/")
setwd(wd)
results_path <- paste0(project_path, "results/")
results_GEMX_CA_path <- paste0(results_path, "GEMX/CellAnnotation/")
results_GEMX_TUMOR_path <- paste0(results_path, "GEMX/TumorAnnotation/")
dir.create(results_GEMX_TUMOR_path, recursive = TRUE, showWarnings = FALSE)

# Import Plot Functions
source(paste0(wd, "CA_plots.R"))

# Load Fully Annotated Data (compartment + Leukocyte/Stromal celltype already in)
dwAnnotated <- readRDS(paste0(results_GEMX_CA_path, "fully_annotated_data.rds"))
cat("\n Fully annotated data loaded \n")

# Subset Tumoral Cells
dwTumoral <- subset(dwAnnotated, subset = compartment == "Tumoral")
dwTumoral <- JoinLayers(dwTumoral)
cat(paste0("\n Tumoral subset: ", ncol(dwTumoral), " cells \n"))


# Set PAM50-like Marker Genes
pam50_genes <- list(
  Lum   = c("ESR1", "PGR", "BAG1", "MAPT", "NAT1", "ZIP6"),
  Basal = c("MKI67", "CCNE1", "ANLN", "CDC20", "EGFR", "MYC"),
  Her2  = c("ERBB2", "GRB7", "BLVRA", "TMEM45B")
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


# Compute PROGENy Pathway Activity (via decoupleR)
net_progeny <- get_progeny(organism = "human", top = 500)

expr_mat <- GetAssayData(dwTumoral, assay = "RNA", layer = "data")

progeny_acts <- run_mlm(mat = expr_mat, network = net_progeny,
                         .source = "source", .target = "target", .mor = "weight",
                         minsize = 5)

# Reshape long (source/condition/score) -> wide matrix (pathway x cell), add as assay
progeny_mat <- progeny_acts %>%
  pivot_wider(id_cols = source, names_from = condition, values_from = score) %>%
  column_to_rownames("source") %>%
  as.matrix()

dwTumoral[["progeny"]] <- CreateAssayObject(data = progeny_mat)
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


# Propagate PAM50_predicted back to the global object, same pattern as Leukocytes
dwAnnotated <- add_subset_annotation_to_global(dwAnnotated, dwTumoral,
                                                 annotation_col = "PAM50_predicted",
                                                 new_col_name = "celltype")

saveRDS(dwTumoral, file.path(results_GEMX_TUMOR_path, "Tumoral_annotated.rds"))
saveRDS(dwAnnotated, file.path(results_GEMX_TUMOR_path, "fully_annotated_data.rds"))
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