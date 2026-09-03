##
##  Single Cell Analysis Step 1b: Integration & Sample Annotation
##  Runs AFTER 1a_qc.R — takes merged_data.rds as input.
##  Normalizes, runs PCA, integrates with Harmony, joins clinical metadata.
##

# Import libraries
library('Seurat')
library(dplyr)
library(tibble)

# Set paths
project_path <- "/home/usuario/PROJECTS/260724_victor_scRNA/"
wd <- paste0(project_path, "codes/scRNA_pipeline/")
setwd(wd)
results_path <- paste0(project_path, "results/")
cellranger_path <- "/home/usuario/DATASETS/scRNAseq/"
results_GEMX_QC_path <- paste0(results_path, "GEMX/QualityControl/")

# Import plot functions
source(paste0(wd, "QC_plots.R"))

# Load merged QC'd data
merged_seurat <- readRDS(file.path(results_GEMX_QC_path, "merged_data.rds"))
cat("\n Merged data loaded \n")

# Load cell filtering table (needed for the subtype QC plot)
cell_counts_df <- read.csv(file.path(results_GEMX_QC_path, "cell_filtering.csv"))

# --- Compare samples after individual QC, before normalizing ---
plot_vln_compare_samples(merged_seurat, results_path = results_GEMX_QC_path)

# Normalize, Scale and Run PCA
merged_seurat <- NormalizeData(merged_seurat)
merged_seurat <- FindVariableFeatures(merged_seurat)
merged_seurat <- ScaleData(merged_seurat)
merged_seurat <- RunPCA(merged_seurat)
cat("\n Normalized, Scaled and ran PCA \n")

# --- Check batch effect in raw PCA, before integration ---
plot_dimplot_by_sample(merged_seurat, reduction = "pca", results_path = results_GEMX_QC_path)

# Integrate w/ Harmony
dwIntegrated <- IntegrateLayers(
    object = merged_seurat, method = HarmonyIntegration,
    orig.reduction = "pca", new.reduction = "harmony",
    verbose = TRUE
)
saveRDS(dwIntegrated, file.path(results_GEMX_QC_path, "integrated_harmony_data.rds"))
cat("\n Data Integrated w/ Harmony \n")

# --- Check batch effect correction after Harmony integration ---
plot_dimplot_by_sample(dwIntegrated, reduction = "harmony", results_path = results_GEMX_QC_path)

# Load sample subtype annotation CSV
annot_file <- paste0(cellranger_path, "260106_carlota_GEMX/Samples_scRNAseq_GEMX-Flex.csv")
annot_df <- read.csv(annot_file, stringsAsFactors = FALSE)
annot_df <- annot_df %>%
    select(scRNAseq_ID, Subtype) %>%
    distinct() %>%
    mutate(
        scRNAseq_ID = as.character(scRNAseq_ID),
        Subtype = as.character(Subtype)
    )

# --- Percentage of cells filtered per sample, colored by tumor subtype ---
plot_pct_filtered_by_subtype(cell_counts_df, annot_df, results_path = results_GEMX_QC_path)
cat("\n Sample annotation data loaded \n")

# Join subtype annotation to object metadata
dwIntegrated@meta.data <- dwIntegrated@meta.data %>%
    rownames_to_column("cell_id") %>%
    left_join(annot_df, by = c("orig.ident" = "scRNAseq_ID")) %>%
    column_to_rownames("cell_id")

dwIntegrated$Subtype <- factor(dwIntegrated$Subtype)
saveRDS(dwIntegrated, file.path(results_GEMX_QC_path, "sample_annotated_data.rds"))
cat("\n Samples Annotated \n")

cat(paste("\n ---- FINISHED INTEGRATION & SAMPLE ANNOTATION ----
    Generated files:
        · integrated_harmony_data.rds
        · sample_annotated_data.rds
    Generated plots:
        · VlnPlot_compare_samples_postQC.png
        · QC_percentage_filtered_with_labels.png
        · DimPlot_(reduction)_bySample.png
        "))
