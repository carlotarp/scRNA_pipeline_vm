##
##  Single Cell Analysis Step 2a: PC and resolution selection
##  Runs AFTER 1c_decontx.R — takes decontx_data.rds as input.
##  Produces the elbow plot and the clustering resolution grid so the
##  resolution can be chosen before running 2b_cluster.R.
##

# Import libraries
library("Seurat")
library(dplyr)
library(tibble)
library(Matrix)

# Set paths
project_path <- "/home/user/PROJECTS/scRNA_vmendez/"
wd <- paste0(project_path, "codes/")
setwd(wd)
results_path <- paste0(project_path, "results/")
results_GEMX_CL_path <- paste0(results_path, "Clustering/")

# Import plot functions
source(paste0(wd, "CL_plots.R"))

# Load DecontX-corrected data (output of 1c_decontx.R)
dwIntegrated <- readRDS(paste0(results_path, "decontx_data.rds"))
cat(" DecontX data loaded \n")

# Find optimal number of PCs
# co1: first PC where cumulative variance > 90% AND individual contribution < 5%
# co2: last PC before a drop > 0.1% (elbow of the variance curve)
# co3: conservative cutoff — the lesser of co1 and co2
pct <- dwIntegrated[["pca_decontX"]]@stdev / sum(dwIntegrated[["pca_decontX"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = TRUE)[1] + 1
co3 <- min(co1, co2)
cat(paste("\n The Optimal number of PCs is", co3, "\n")) # 12

# --- Elbow plot to visualize co3 ---
plot_elbow(dwIntegrated, co3, results_GEMX_CL_path)

# --- Clustering resolution grid + sankey of cluster splits across resolutions ---
# plot_resolution_grid() draws each resolution on umap_decontX, so the UMAP has
# to exist before the grid (2b_cluster.R recomputes it for the exported object).
dwIntegrated <- RunUMAP(dwIntegrated, reduction = "harmony_decontX", dims = 1:co3,
                        assay = "RNA_decontX", reduction.name = "umap_decontX")
dwIntegrated <- FindNeighbors(dwIntegrated, reduction = "harmony_decontX", dims = 1:co3,
                              assay = "RNA_decontX", graph.name = c("RNA_decontX_nn", "RNA_decontX_snn"))
plot_resolution_grid(dwIntegrated, results_path = results_GEMX_CL_path,
                          reduction = "umap_decontX",
                          resolutions = c(0.2, 0.3, 0.4, 0.5))

cat(paste("\n ---- FINISHED PC & RESOLUTION SELECTION ----
    Review ElbowPlot.png and ResolutionGrid.png, set `res` at the top of
    2b_cluster.R, then run 2b_cluster.R.
    Generated plots:
      · ElbowPlot.png
      · ResolutionGrid.png
      · SankeyPlot.html
          "))
