##
##  Single Cell Analysis Step 2a: Sample Annotation & Clustering
##

# Import libraries
library("Seurat")
library(dplyr)
library(tibble)
library(Matrix)

# Set paths
project_path <- "/home/usuario/PROJECTS/260724_victor_scRNA/"
wd <- paste0(project_path, "codes/scRNA_pipeline/")
setwd(wd)
results_path <- paste0(project_path, "results/")
results_GEMX_CL_path <- paste0(results_path, "GEMX/DecontX/Clustering/")

# Import Plot Functions
source(paste0(wd, "CL_plots.R"))

# Load Sample Annotated Data
dwIntegrated <- readRDS(paste0(results_path, "GEMX/QualityControl/7500/sample_annotated_data.rds"))
cat(" Annotated Sample Data Loaded \n")

# Find Optimal Number of PCs
pct <- dwIntegrated[["pca_decontX"]]@stdev / sum(dwIntegrated[["pca_decontX"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = TRUE)[1] + 1
co3 <- min(co1, co2)
cat(paste("\n The Optimal number of PCs is", co3, "\n")) # 12

# --- Elbow plot to visualize co3 ---
plot_elbow(dwIntegrated, co3, results_GEMX_CL_path)

# Generate UMAP
dwIntegrated <- RunUMAP(dwIntegrated, reduction = "harmony_decontX", dims = 1:co3, reduction.name = "umap_decontX")
cat("\n UMAP Generated \n")

# Find Optimal Resolution for Clustering (visually)
dwIntegrated <- FindNeighbors(dwIntegrated, reduction = "harmony_decontX", dims = 1:co3)
plot_resolution_grid(dwIntegrated, results_path = results_GEMX_CL_path,
                          reduction = "harmony_decontX",
                          resolutions = c(0.2, 0.3, 0.4, 0.5))
res <- 0.5
cat(paste("\n The Optimal Resolution is", res, "\n"))

# Generate Clusters
dwIntegrated <- FindClusters(dwIntegrated, resolution = res, cluster.name = "decontX_clusters")
cat("\n Clusters Generated w/ Optimal Resolution \n")


# ---  Visualize Integrated PCA ---
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

# Generate Tables of Interest
cells_clusters <- table(dwIntegrated@meta.data$seurat_clusters, dwIntegrated@meta.data$orig.ident)
cells_clusters <- cbind(cells_clusters,row.names(cells_clusters))
write.table(cells_clusters,file=paste0(results_GEMX_CL_path,'cells_per_cluster.tsv'),sep="\t",row.names = FALSE, col.names = TRUE)
write.table(dwIntegrated@reductions[["umap_decontX"]]@cell.embeddings,file=paste0(results_GEMX_CL_path,'umap_pvalues.tsv'),sep="\t",row.names = TRUE, col.names = TRUE)
write.table(dwIntegrated@meta.data,file=paste0(results_GEMX_CL_path,'clustered_metadata.tsv'),sep="\t",row.names = TRUE, col.names = TRUE)
cat("\n Tables of Interest Writen \n")

# --- Visualize Cluster Composition ---
plot_cluster_composition(dwIntegrated, group_by = "orig.ident",
                          results_path = results_GEMX_CL_path, filename = "ClusterComposition_bySample.png")
plot_cluster_composition(dwIntegrated, group_by = "Subtype",
                          results_path = results_GEMX_CL_path, filename = "ClusterComposition_bySubtype.png")

# --- Quality Control Plots ---
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

# Export Clustered Data
saveRDS(dwIntegrated, file.path(results_GEMX_CL_path, "clustered_data.rds"))

# Generate RNA Expression Matrix
dwIntegrated <- JoinLayers(dwIntegrated)
count_matrix <- GetAssayData(dwIntegrated, layer = "counts")
umi_counts_df <- as.data.frame(as.matrix(count_matrix))
write.table(umi_counts_df,file=paste0(results_GEMX_CL_path,'umi.tsv'),sep="\t",row.names = TRUE, col.names = TRUE)
rm(umi_counts_df)
rm(count_matrix)
print('\n RNA expression table generated \n')

# Identify Cluster de Marker Genes
tumor_markers_all <- FindAllMarkers(dwIntegrated, only.pos = FALSE, logfc.threshold = -1)
write.table(tumor_markers_all,file=paste0(results_GEMX_CL_path,'all_tumor_markers.tsv'),sep="\t",row.names = FALSE, col.names = TRUE)
print('\n Diferential Expression of tumor.markers Performed \n')

# --- Heatmap to Visualize Cluster  ---
plot_heatmap(dwIntegrated, results_GEMX_CL_path, tumor_markers_all, n = 10)

cat(paste("\n ---- FINISHED CLUSTERING ----
    Generated files:
      · cells_per_cluster.tsv
      · umap_pvalues.tsv
      · clustered_metadata.tsv
      · pairwisett_nfeaturerna_subtype.csv
      · clustered_data.rds
      · umi.tsv
      · all_tumor_markers.tsv
    Generated plots:
      · ElbowPlot.png
      · ResolutionGrid.png
      · SankeyPlot.html
      · DimPlot_(reduction)_(groupedby).png
      · FeaturePlot_(reduction)_(features).png
      · ClusterComposition_(groupedby).png
      · VlnPlot_QCmetrics_(groupedby).png
      · Integrated_Top(n)gene_Heatmap.png
          "))