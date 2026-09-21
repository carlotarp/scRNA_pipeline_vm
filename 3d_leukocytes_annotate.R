##
##  Single Cell Analysis Step 3d: Leukocyte annotation
##  Runs AFTER 3c_leukocytes_markers.R — takes lineage_annotated_data.rds and
##  leukocytes.rds as input.
##  Applies the manual leukocyte dictionary and transfers the labels back to
##  the full object.
##

# Import libraries
library("Seurat")
library(dplyr)

# Set paths
project_path <- "/home/user/PROJECTS/scRNA_vmendez/"
wd <- paste0(project_path, "codes/")
setwd(wd)
results_path <- paste0(project_path, "results/")
results_GEMX_CA_path <- paste0(results_path, "CellAnnotation/")
results_GEMX_LCA_path <- paste0(results_GEMX_CA_path, "Leukocytes/")

# Import plot functions and shared utilities
source(paste0(wd, "CA_plots.R"))
source(paste0(wd, "CL_plots.R"))
source(paste0(wd, "utils.R"))

# Manual cluster annotation, filled in from the 3c_leukocytes_markers.R outputs
clusters_leuko_annotated <- c( #########  DecontX 7500  #########
  "0" = "TAM",  "1" = "TCell_cyto",  "2" = "BCell",
  "3" = "TCell_naive",  "4" = "PlasmaBlast",  "5" = "TCell_ex",
  "6" = "cDC",  "7" = "TAM_fibro",  "8" = "Mast",
  "9" = "pDC",  "10" = "Proliferative",  "11" = "TAM_proInflamm",
  "12" = "actDC", "13" = "PlasmaBlast"
)

# Load lineage annotated data and the reclustered leukocyte subset
dwAnnotated <- readRDS(paste0(results_GEMX_CA_path, "lineage_annotated_data.rds"))
dwLeukocytes <- readRDS(file.path(results_GEMX_LCA_path, "leukocytes.rds"))
cat("\n lineage-annotated data and leukocyte subset loaded \n")

dwLeukocytes$celltype <- unname(clusters_leuko_annotated[as.character(dwLeukocytes$seurat_clusters)])
dwLeukocytes$celltype <- factor(dwLeukocytes$celltype)

# ---  Visualize Annotated Dimplot ---
plot_dimplot(dwLeukocytes, reduction = "umap", group_by = "celltype", label = T,
             results_path = results_GEMX_LCA_path, filename = "DimPlot_UMAP_Leuko_Annotated.png")


# Export annotated leukocyte data
saveRDS(dwLeukocytes, file.path(results_GEMX_LCA_path, "leukocytes.rds"))

# Label transfer to full object
annotation_vec <- setNames(as.character(dwLeukocytes[["celltype"]][, 1]),
                            colnames(dwLeukocytes))

if ("celltype" %in% colnames(dwAnnotated@meta.data)) {
  existing <- as.character(dwAnnotated[["celltype"]][, 1])
} else {
  existing <- rep(NA_character_, ncol(dwAnnotated))
}
names(existing) <- colnames(dwAnnotated)
existing[names(annotation_vec)] <- annotation_vec
dwAnnotated[["celltype"]] <- factor(existing)

# ---  Visualize Annotated Dimplot ---
plot_dimplot(dwAnnotated, reduction = "umap_decontX", group_by = "celltype", label = T,
             results_path = results_GEMX_CA_path, filename = "DimPlot_UMAP_Leuko.png")

# Export Annotated Data
saveRDS(dwAnnotated, file.path(results_GEMX_CA_path, "leuko_annotated_data.rds"))

cat(paste("\n ---- FINISHED LEUKOCYTE ANNOTATION ----
    Run 3e_stroma_markers.R next.
    Generated files:
        · leukocytes.rds
        · leuko_annotated_data.rds
    Generated plots:
        · DimPlot_UMAP_(groupedby).png
        "))
