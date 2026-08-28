##
##  Single Cell Analysis Step 3f: Tumor Subtype Annotation
##

# Import libraries
library("Seurat")
library(dplyr)
library(decoupleR)
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

# Load Plot Functions
source(paste0(wd, "CA_plots.R"))
source(paste0(wd, "CL_plots.R"))

dwTumoral <- readRDS(paste0(results_GEMX_TUMOR_path, "tumoral_annotated.rds"))

net_progeny <- get_progeny(organism = "human", top = 500)

expr_mat <- as.matrix(GetAssayData(dwTumoral, assay = "RNA_decontX", layer = "data"))

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

saveRDS(dwTumoral, paste0(results_GEMX_TUMOR_path, "tumoral_annotated.rds"))
