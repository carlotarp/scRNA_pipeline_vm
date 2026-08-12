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
results_GEMX_CNV_path <- paste0(results_path, "GEMX/CellAnnotation/3500/InferCNV/")

# Load Annotated Data
dwAnnotated <- readRDS(paste0(results_path, "GEMX/CellAnnotation/3500/notumor_annotated_data.rds"))
cat("\n Annotated Data Loaded \n")

dwAnnotated <- CellCycleScoring(dwAnnotated,
                                s.features = cc.genes.updated.2019$s.genes,
                                g2m.features = cc.genes.updated.2019$g2m.genes)

# Classify cells
dwAnnotated$dataset <- ifelse(dwAnnotated$celltype == "Tumor",
                               paste0("TUMOR_", dwAnnotated$orig.ident),
                               as.character(dwAnnotated$celltype))
dwAnnotated$dataset <- factor(dwAnnotated$dataset)
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

# Set Reference Cluster
ref_groups <- setdiff(levels(dwAnnotated$dataset), grep("^TUMOR_", levels(dwAnnotated$dataset), value = TRUE))
ref_groups <- setdiff(ref_groups, c("Prolifetarive", "Fibrocyte"))

# Create inferCNV object
infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = counts_matrix,
  annotations_file = paste0(results_GEMX_CNV_path, "inferCNV_cell_annotations.txt"),
  delim = "\t",
  gene_order_file = "/home/usuario/PROJECTS/260724_victor_scRNA/data/gencode_v19_gene_pos.txt",
  ref_group_names = ref_groups
  )

# Run inferCNV
infercnv_obj_full_run <- infercnv::run(
  infercnv_obj,
  cutoff = 0.1,
  out_dir = results_GEMX_CNV_path,
  cluster_by_groups = TRUE,
  denoise = TRUE,
  HMM = TRUE,
  analysis_mode = "samples",
  no_plot = T,
  no_prelim_plot = T,
  output_format = "pdf",
  num_threads = 15

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