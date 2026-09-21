##
##  Single Cell Analysis Step 3f: Stroma annotation
##  Runs AFTER 3e_stroma_markers.R — takes leuko_annotated_data.rds and
##  Stroma.rds as input.
##  Applies the manual stroma dictionary and transfers the labels back to the
##  full object, producing notumor_annotated_data.rds.
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
results_GEMX_SCA_path <- paste0(results_GEMX_CA_path, "Stroma/")

# Import plot functions and shared utilities
source(paste0(wd, "CA_plots.R"))
source(paste0(wd, "CL_plots.R"))
source(paste0(wd, "utils.R"))

# Manual cluster annotation, filled in from the 3e_stroma_markers.R outputs
clusters_Stroma_annotated <- c( #########  DecontX 7500  #########
  "0" = "mCAF",  "1" = "Endothelial",  "2" = "mCAF",
  "3" = "cCAF",  "4" = "iCAF",  "5" = "Pericyte",
  "6" = "Adipocyte", "7" = "Lymphatic"
)

# Load leukocyte annotated data and the reclustered stromal subset
dwAnnotated <- readRDS(paste0(results_GEMX_CA_path, "leuko_annotated_data.rds"))
dwStroma <- readRDS(file.path(results_GEMX_SCA_path, "Stroma.rds"))
cat("\n leukocyte-annotated data and stromal subset loaded \n")

dwStroma$celltype <- unname(clusters_Stroma_annotated[as.character(dwStroma$seurat_clusters)])
dwStroma$celltype <- factor(dwStroma$celltype)

# ---  Visualize Annotated Dimplot ---
plot_dimplot(dwStroma, reduction = "umap", group_by = "celltype", label = T,
             results_path = results_GEMX_SCA_path, filename = "DimPlot_UMAP_Stroma_Annotated.png")


# Export annotated stroma data
saveRDS(dwStroma, file.path(results_GEMX_SCA_path, "Stroma.rds"))

# Label transfer to full object
annotation_vec <- setNames(as.character(dwStroma[["celltype"]][, 1]),
                            colnames(dwStroma))

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
             results_path = results_GEMX_CA_path, filename = "DimPlot_UMAP_Stroma.png")

# Export Annotated Data
saveRDS(dwAnnotated, file.path(results_GEMX_CA_path, "notumor_annotated_data.rds"))

cat(paste("\n ---- FINISHED STROMA ANNOTATION ----
    Run 3g_copykat.R next.
    Generated files:
        · Stroma.rds
        · notumor_annotated_data.rds
    Generated plots:
        · DimPlot_UMAP_(groupedby).png
        "))
