##
##  Single Cell Analysis Step 3c: Cell Annotation
##

# Import libraries
library(dplyr)
library(Seurat)
library(SeuratObject)
library(infercnv)

# Set Paths
project_path <- "/home/usuario/PROJECTS/260724_victor_scRNA/"
wd <- paste0(project_path, "codes/scRNA_pipeline/")
setwd(wd)
results_path <- paste0(project_path, "results/")
results_GEMX_CNV_path <- paste0(results_path, "GEMX/CellAnnotation/InferCNV/")

# Load Annotated Data
dwAnnotated <- readRDS(paste0(results_path, "GEMX/CellAnnotation/lineage_annotated_data.rds"))
cat("\n Annotated Data Loaded \n")

# Classify cells
clusters_normal <- c(3, 4, 11, 12) # Clusters we know are normal

dwAnnotated$dataset <- ifelse(dwAnnotated$lineage == "Tumoral",
                               paste0("TUMOR_", dwAnnotated$orig.ident),
                               "NORMAL")
dwAnnotated$dataset <- factor(dwAnnotated$dataset, levels = c("NORMAL", "TUMOR"))


# inferCNV needs raw counts
dwAnnotated <- JoinLayers(dwAnnotated)
counts_matrix <- GetAssayData(dwAnnotated, assay = "RNA", layer = "counts")

# Annotation file for inferCNV
cell_annotations <- data.frame(
  cell = colnames(dwAnnotated),
  group = dwAnnotated$dataset
)

write.table(
  cell_annotations,
  file = paste0(results_GEMX_CNV_path, "inferCNV_cell_annotations.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

# Create inferCNV object
infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = counts_matrix,
  annotations_file = paste0(results_GEMX_CNV_path, "inferCNV_cell_annotations.txt"),
  delim = "\t",
  gene_order_file = "/home/usuario/PROJECTS/260724_victor_scRNA/data/gencode_v19_gene_pos.txt",
  ref_group_names = c("NORMAL")
)

# Run inferCNV
infercnv_obj_full_run <- infercnv::run(
  infercnv_obj,
  cutoff = 0.1,
  out_dir = results_GEMX_CNV_path,
  cluster_by_groups = TRUE,
  denoise = TRUE,
  HMM = TRUE,
  analysis_mode = "samples"
)

# Generate Plot
plot_cnv(
  infercnv_obj_full_run,
  out_dir = results_GEMX_CNV_path,
  obs_title = "Input cells",
  ref_title = "References (NORMAL clusters)",
  cluster_by_groups = TRUE,
  cluster_references = TRUE,
  dynamic_resize = 0.2,
  output_format = "png",
  png_res = 300
)