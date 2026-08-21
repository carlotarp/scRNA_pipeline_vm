##
##  Single Cell Analysis: Convert .rds Data Format to .h5ad
##

# Set environment
Sys.setenv(RETICULATE_PYTHON = "/home/usuario/miniconda3/envs/seurat5/bin/python")

# Import Libraries
library(reticulate)
library(Seurat)
library(rliger)
library(scCustomize)
library(qs)

# Set Paths
project_path <- "/home/usuario/PROJECTS/260724_victor_scRNA/"
wd <- paste0(project_path, "codes/scRNA_pipeline/")
setwd(wd)
results_path <- paste0(project_path, "results/")

rds_to_h5ad <- function(object_path, file_path, file_name){
    #Load Clustered Sample Data
    object <- readRDS(object_path)
    object <- JoinLayers(object)
    print("\n Samples Loaded \n")

    #Export Normalized Data in h5ad Format
  as.anndata(
    x = object,
    file_path = file_path,
    file_name = file_name,
    assay = "RNA",
    main_layer = "counts",
    other_layers = c("data"),
    verbose = TRUE
  )}

object_path <- paste0(results_path, "GEMX/CellAnnotation/7500/fully_annotated_data.rds")
file_path <- paste0(results_path, "GEMX/CellCommunication/7500/")
file_name <- "fully_annotated_data.h5ad"

rds_to_h5ad(object_path = object_path, file_path = file_path, file_name = file_name)