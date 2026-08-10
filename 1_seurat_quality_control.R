##
##  Single Cell Analysis Step 1: Quality Control, Integration & Sample Annotation
##

# Import libreries
library('Seurat')
library(dplyr)
library(tibble)
library(SingleCellExperiment)
library("scDblFinder")

# Set paths
project_path <- "/home/usuario/PROJECTS/260724_victor_scRNA/"
wd <- paste0(project_path, "codes/scRNA_pipeline/")
setwd(wd)
results_path <- paste0(project_path, "results/")
cellranger_path <- "/home/usuario/DATASETS/scRNAseq/"
GEMX_path <- paste0(cellranger_path, "260106_carlota_GEMX/2026_HN00264849/allPool/")
data_path <- "_filtered_feature_barcode_matrix/"
results_GEMX_QC_path <- paste0(results_path, "GEMX/QualityControl/")

# Import Plot Functions
source(paste0(wd, "QC_plots.R"))

# Set lists
sample_list <- list()
cell_counts <- list()
cell_counts_names <- list()
sample_cell_counts <- list()
doublet_list <- list()

# Load sample Function
load_read10X_sample <- function(cellranger_path, name, data_path, min_cells=3, min_features=200){
    cat(paste("\n Loading Sample", name, "...\n", sep = " "))
    sample_path <- paste0(cellranger_path, name, data_path)

    tumor <- Read10X(data.dir = sample_path)
    tumor <- CreateSeuratObject(counts = tumor, project = name, min.cells = min_cells, min.features = min_features)
    loaded.cells <- as.numeric(nrow(tumor@meta.data))
    cell_counts <<- append(cell_counts, loaded.cells)
    cell_counts_names <<- append(cell_counts_names,'Loaded cells')

    return(tumor)
}

# Process Sample Function
process_sample <- function(name, tumor, mito_cut, feature_lower_cut, feature_upper_cut, results_path){
    cat(paste("\n Processing Sample", name, "...\n", sep = " "))

    tumor[["log10_nCount_RNA"]] <- log10(tumor$nCount_RNA)
    tumor[["percent.mt"]] <- PercentageFeatureSet(object = tumor, pattern = "^MT-")
    tumor[["percent.hb"]] <- PercentageFeatureSet(object = tumor, pattern = "^HB[^(P|E|S)]")

    # --- Plots PRE-filter (raw distributions, before any subset) ---
    plot_qc_vlnplot(tumor, name = name, results_path = results_path)
    plot_mito_vs_nfeature(tumor, name = name, results_path = results_path,
         mito_cut = mito_cut,
         feature_lower_cut = feature_lower_cut,
         feature_upper_cut = feature_upper_cut)

    # Filter by Mitochondrial Genes
    sample <- subset(x = tumor, subset = percent.mt < mito_cut)
    cell_counts <<- append(cell_counts, as.numeric(nrow(sample@meta.data)))
    cell_counts_names <<- append(cell_counts_names, "Mitoch")

    # Filter by Gene Count
    sample <- subset(x = sample, subset = nFeature_RNA > feature_lower_cut & nFeature_RNA < feature_upper_cut)
    cell_counts <<- append(cell_counts, as.numeric(nrow(sample@meta.data)))
    cell_counts_names <<- append(cell_counts_names, "QC gene num ABS")

    # Filter by SD
    upper_bound <- mean(sample$log10_nCount_RNA) + 2*sd(sample$log10_nCount_RNA)
    lower_bound <- mean(sample$log10_nCount_RNA) - 2*sd(sample$log10_nCount_RNA)
    sample <- subset(x = sample, subset = log10_nCount_RNA > lower_bound & log10_nCount_RNA < upper_bound)
    cell_counts <<- append(cell_counts, as.numeric(nrow(sample@meta.data)))
    cell_counts_names <<- append(cell_counts_names, "QC gene num SD")

    # Filter repeats
    sce <- as.SingleCellExperiment(sample)
    sce <- scDblFinder(sce)
    sample$scDblFinder.class <- colData(sce)$scDblFinder.class
    sample$scDblFinder.score <- colData(sce)$scDblFinder.score

    df_snapshot <- sample@meta.data
    df_snapshot$SampleID <- name
    doublet_list[[name]] <<- df_snapshot

    # --- Doublet score distribution, BEFORE removing doublets ---
    plot_doublet_scores(sample, name = name, results_path = results_path)
    plot_doublet_scatter(df = sample@meta.data, name = name, results_path = results_path)

    sample <- subset(x = sample, subset = scDblFinder.class == "singlet")
    cell_counts <<- append(cell_counts, nrow(sample@meta.data))
    cell_counts_names <<- append(cell_counts_names, "scDblFinder")

    # Append data
    sample_cell_counts <<- append(sample_cell_counts, c(name, name, name, name, name))
    sample_list[[name]] <<- sample

    return(tumor)
}

# Set inputs for Load Data Function for GEMX
min_cells <- 3
min_features <- 200

# Set inputs for Preprocess Data Function
mito_cut <- 10
feature_lower_cut <- 100
feature_upper_cut <- 7500

# Load + Process Loop
sample_names <- c('SC7b', 'SC8','SC9','SC10','SC11','SC12','SC13','SC14','SC15','SC16','SC17','SC18','SC19','SC20','SC21','SC5_SC22') #GEMX

for (name in sample_names){
    tumor <- load_read10X_sample(cellranger_path = GEMX_path, name = name, data_path = data_path, min_cells = min_cells, min_features = min_features)
    tumor <- process_sample(name = name, tumor = tumor, mito_cut = mito_cut, feature_lower_cut = feature_lower_cut, feature_upper_cut = feature_upper_cut, results_path = results_GEMX_QC_path)

    # --- Plots POST-filter ---
    plot_gene_expression_dist(tumor = tumor, name = name, results_path = results_GEMX_QC_path)
    plot_nfeature_ncount_corr(tumor = tumor, name = name, results_path = results_GEMX_QC_path)
}
cat("\n Samples Loaded and Processed \n")

# --- Doublet score distribution (all samples) ---
doublet_df <- do.call(rbind, doublet_list)
plot_doublet_scatter(doublet_df, name = "AllSamples", results_path = results_GEMX_QC_path)

# Write Cell filtering Table
cell_counts_df <- data.frame('Cell_counts'=c(unlist(cell_counts)),
                             'Feature_discard'=c(unlist(cell_counts_names)),
                              'SampleID'=c(unlist(sample_cell_counts)))
write.csv(cell_counts_df, file = file.path(results_GEMX_QC_path, "cell_filtering.csv"))
cat("\n Cell Filtering Table Generated \n")

# --- Cells retained per filtering step ---
plot_qc_funnel(cell_counts_df, results_path = results_GEMX_QC_path)

# Clean Environment
rm(tumor)
rm(cell_counts)
rm(sample_cell_counts)
rm(cell_counts_names)

# Merge Samples
merged_seurat <- merge(x = sample_list[[1]],
                       y = sample_list[2:length(sample_list)])
saveRDS(merged_seurat, file.path(results_GEMX_QC_path, "merged_data.rds"))
cat("\n Seurat Objects Merged \n")

# --- Compare samples after individual QC, before normalizing ---
plot_vln_compare_samples(merged_seurat, results_path = results_GEMX_QC_path)

# Normalize, Scale and Run PCA
merged_seurat <- NormalizeData(merged_seurat)
merged_seurat <- FindVariableFeatures(merged_seurat)
merged_seurat <- ScaleData(merged_seurat)
merged_seurat <- RunPCA(merged_seurat)
cat("\n Normalized, Scaled and ran PCA")

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

# Add Sample Subtype Annotations
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

cat("\n Sample Annotation Data Loaded \n")

dwIntegrated@meta.data <- dwIntegrated@meta.data %>%
    rownames_to_column("cell_id") %>%
    left_join(annot_df, by = c("orig.ident" = "scRNAseq_ID")) %>%
    column_to_rownames("cell_id")

dwIntegrated$Subtype <- factor(dwIntegrated$Subtype)
saveRDS(dwIntegrated, file.path(results_GEMX_QC_path, "sample_annotated_data.rds"))
cat("\n Samples Annotated \n")

cat(paste("\n ---- FINISHED QUALITY CONTROL, INTEGRATION & SAMPLE ANNOTATION ----
    Generated files:
        · cell_filtering.csv
        · merged_data.rds
        · integrated_harmony_data.rds
        · annotated_sample_data.rds
    Generated plots:
        · nCountRNA_vs_nFeatureRNA_(sample).png
        · nFeatureRNA_(sample).png
        · VlnPlot_QCmetrics_(sample).png
        · MitoVsNFeature_(sample).png
        · DoubletScore_(sample).png
        · DoubletScatter_(sample).png
        · VlnPlot_compare_samples_postQC.png
        · QC_cell_counts_subplots.png
        · QC_percentage_filtered_with_labels.png
        · DimPlot_(reduction)_by_sample.png
        "))