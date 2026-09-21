##
##  Single Cell Analysis Step 2b: Clustering
##  Runs AFTER 2a_resolution.R — takes decontx_data.rds as input.
##  Runs UMAP and clustering at the resolution chosen from the grid in 2a,
##  and exports clustered_data.rds.
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

# Clustering resolution chosen from ResolutionGrid.png (2a_resolution.R)
res <- 0.5

# Load DecontX-corrected data (output of 1c_decontx.R)
dwIntegrated <- readRDS(paste0(results_path, "decontx_data.rds"))
cat(" DecontX data loaded \n")

# Optimal number of PCs (same criterion as 2a_resolution.R), per reduction
optimal_pcs <- function(object, reduction) {
  pct <- object[[reduction]]@stdev / sum(object[[reduction]]@stdev) * 100
  cumu <- cumsum(pct)
  co1 <- which(cumu > 90 & pct < 5)[1]
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = TRUE)[1] + 1
  min(co1, co2)
}

co3 <- optimal_pcs(dwIntegrated, "pca_decontX")
co3_rna <- optimal_pcs(dwIntegrated, "pca")
cat(paste("\n Optimal PCs: ", co3, " (RNA_decontX), ", co3_rna, " (RNA) \n"))

# Uncorrected (pre-DecontX) baseline on the RNA assay: same criterion and same
# resolution as the DecontX branch below, but its own assay, PCs, reduction and
# graph, so neither run overwrites the other's neighbour graph.
dwIntegrated <- RunUMAP(dwIntegrated, reduction = "harmony", dims = 1:co3_rna,
                        assay = "RNA", reduction.name = "umap")
dwIntegrated <- FindNeighbors(dwIntegrated, reduction = "harmony", dims = 1:co3_rna,
                              assay = "RNA", graph.name = c("RNA_nn", "RNA_snn"))
dwIntegrated <- FindClusters(dwIntegrated, resolution = res, graph.name = "RNA_snn",
                             cluster.name = "contaminated_clusters")
cat("\n Uncorrected (RNA) baseline UMAP and clusters generated \n")

# Generate UMAP
dwIntegrated <- RunUMAP(dwIntegrated, reduction = "harmony_decontX", dims = 1:co3,
                        assay = "RNA_decontX", reduction.name = "umap_decontX")
cat("\n UMAP Generated \n")

# Generate clusters
dwIntegrated <- FindNeighbors(dwIntegrated, reduction = "harmony_decontX", dims = 1:co3,
                              assay = "RNA_decontX", graph.name = c("RNA_decontX_nn", "RNA_decontX_snn"))
dwIntegrated <- FindClusters(dwIntegrated, resolution = res, graph.name = "RNA_decontX_snn",
                             cluster.name = "decontX_clusters")
cat(paste("\n Clusters Generated w/ resolution", res, "\n"))

# --- Visualize integrated PCA ---
plot_dimplot(dwIntegrated, reduction = "harmony_decontX", group_by = "Subtype",
             results_path = results_GEMX_CL_path, filename = "DimPlot_PCA_bySubtype.png")
plot_dimplot(dwIntegrated, reduction = "harmony_decontX", group_by = "seurat_clusters", label = T,
             results_path = results_GEMX_CL_path, filename = "DimPlot_PCA_byCluster.png")
plot_featureplot(dwIntegrated, reduction = "harmony_decontX", features = "PTPRC",
                  results_path = results_GEMX_CL_path, filename = "FeaturePlot_PCA_CD45.png")

# --- Visualize UMAP ---
plot_dimplot(dwIntegrated, reduction = "umap_decontX", group_by = "Subtype",
             results_path = results_GEMX_CL_path, filename = "DimPlot_UMAP_bySubtype.png")
plot_dimplot(dwIntegrated, reduction = "umap_decontX", group_by = "seurat_clusters", label = T,
             results_path = results_GEMX_CL_path, filename = "DimPlot_UMAP_byCluster.png")
plot_featureplot(dwIntegrated, reduction = "umap_decontX", features = "PTPRC",
                  results_path = results_GEMX_CL_path, filename = "FeaturePlot_UMAP_CD45.png")

# Generate tables of interest
cells_clusters <- table(dwIntegrated@meta.data$seurat_clusters, dwIntegrated@meta.data$orig.ident)
cells_clusters <- cbind(cells_clusters,row.names(cells_clusters))
write.table(cells_clusters,file=paste0(results_GEMX_CL_path,'cells_per_cluster.tsv'),sep="\t",row.names = FALSE, col.names = TRUE)
write.table(dwIntegrated@reductions[["umap_decontX"]]@cell.embeddings,file=paste0(results_GEMX_CL_path,'umap_pvalues.tsv'),sep="\t",row.names = TRUE, col.names = TRUE)
write.table(dwIntegrated@meta.data,file=paste0(results_GEMX_CL_path,'clustered_metadata.tsv'),sep="\t",row.names = TRUE, col.names = TRUE)
cat("\n Tables of Interest Writen \n")

# --- Visualize cluster composition ---
plot_composition(dwIntegrated, group_by = "orig.ident",
                          results_path = results_GEMX_CL_path, filename = "ClusterComposition_bySample.png")
plot_composition(dwIntegrated, group_by = "Subtype",
                          results_path = results_GEMX_CL_path, filename = "ClusterComposition_bySubtype.png")

# --- Quality control plots ---
plot_vln_qc_by_group(dwIntegrated, results_path = results_GEMX_CL_path,
                      filename = "VlnPlot_QCmetrics_bySample.png", group_by = "orig.ident")
plot_vln_qc_by_group(dwIntegrated, results_path = results_GEMX_CL_path,
                      filename = "VlnPlot_QCmetrics_byCluster.png", group_by = "seurat_clusters")
plot_vln_qc_by_group(dwIntegrated, results_path = results_GEMX_CL_path,
                      filename = "VlnPlot_QCmetrics_bySubtype.png", group_by = "Subtype")

# Pairwise t-test for nFeature_RNA across Subtypes
subtypes <- levels(factor(dwIntegrated$Subtype))
pairwise_results <- list()

for (i in 1:(length(subtypes) - 1)) {
  for (j in (i + 1):length(subtypes)) {
    s1 <- subtypes[i]
    s2 <- subtypes[j]

    x <- dwIntegrated$nFeature_RNA[dwIntegrated$Subtype == s1]
    y <- dwIntegrated$nFeature_RNA[dwIntegrated$Subtype == s2]

    tt <- t.test(x, y)

    pairwise_results[[paste(s1, "vs", s2)]] <- data.frame(
      group1 = s1,
      group2 = s2,
      mean_group1 = mean(x, na.rm = TRUE),
      mean_group2 = mean(y, na.rm = TRUE),
      t_statistic = unname(tt$statistic),
      df = unname(tt$parameter),
      p_value = tt$p.value,
      conf_low = tt$conf.int[1],
      conf_high = tt$conf.int[2],
      stringsAsFactors = FALSE
    )
  }
}

pairwise_results_df <- do.call(rbind, pairwise_results)
write.csv(pairwise_results_df, paste0(results_GEMX_CL_path, "pairwisett_nfeaturerna_subtype.csv"), row.names = FALSE)

# Export clustered data
saveRDS(dwIntegrated, file.path(results_GEMX_CL_path, "clustered_data.rds"))

# Generate RNA expression matrix
dwIntegrated <- JoinLayers(dwIntegrated)
count_matrix <- GetAssayData(dwIntegrated, layer = "counts")
umi_counts_df <- as.data.frame(as.matrix(count_matrix))
write.table(umi_counts_df,file=paste0(results_GEMX_CL_path,'umi.tsv'),sep="\t",row.names = TRUE, col.names = TRUE)
rm(umi_counts_df)
rm(count_matrix)
cat("\n RNA expression table generated \n")

cat(paste("\n ---- FINISHED CLUSTERING ----
    Run 2c_markers.R next to compute FindAllMarkers (slow step, separate script).
    Generated files:
      · cells_per_cluster.tsv
      · umap_pvalues.tsv
      · clustered_metadata.tsv
      · pairwisett_nfeaturerna_subtype.csv
      · clustered_data.rds
      · umi.tsv
    Generated plots:
      · DimPlot_(reduction)_(groupedby).png
      · FeaturePlot_(reduction)_(features).png
      · ClusterComposition_(groupedby).png
      · VlnPlot_QCmetrics_(groupedby).png
          "))
