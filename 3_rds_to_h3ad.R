##
##  Single Cell Analysis Step 2.5: Convert Data Format to .h5ad
##

# Set environment
Sys.setenv(RETICULATE_PYTHON = "/home/usuario/miniconda3/envs/seurat5/bin/python")

# Import Libraries
library(reticulate)
library(Seurat)
library(rliger)
library(scCustomize) #installed with remotes and the rest of packages with conda-forge
library(qs)

# Set Paths
project_path <- "/home/usuario/PROJECTS/260724_victor_scRNA/"
wd <- paste0(project_path, "codes/scRNA_pipeline/")
setwd(wd)
results_path <- paste0(project_path, "results/")
results_GEMX_CL_path <- paste0(results_path, "GEMX/Clustering/")

#Load Clustered Sample Data
dwClustered <- paste0(results_path,"GEMX/Clustering/clustered_data.rds")
dwClustered <- readRDS(dwClustered)
dwClustered <- JoinLayers(dwClustered)
print("\n Samples Loaded \n")

#Export Normalized Data in h5ad Format
as.anndata(x = dwClustered, file_path = results_GEMX_CL_path, file_name = "clustered_normalized_data.h5ad")