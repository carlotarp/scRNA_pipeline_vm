##
##  Comparative analysis of CEACAM6, FGFR4, ERBB2 across multiple dimensions:
##  cluster (Tumoral AND global), sample, PAM50 subtype, QC quality, and
##  co-expression between the three genes.
##

library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)

# Set Paths
project_path <- "/home/usuario/PROJECTS/260724_victor_scRNA/"
wd <- paste0(project_path, "codes/scRNA_pipeline/")
setwd(wd)
results_path <- paste0(project_path, "results/")
results_GEMX_GC_path <- paste0(results_path, "GEMX/GeneComparison/5500/")

source(paste0(wd, "CA_plots.R"))

# Load  Data
dwAnnotated <- readRDS(paste0(results_path, "GEMX/CellAnnotation/5500/fully_annotated_data.rds"))
dwTumoral <- readRDS(paste0(results_path, "GEMX/CellAnnotation/5500/Tumor/tumoral_annotated.rds"))
cat("\n Data Loaded \n")


genes_check <- c("CEACAM6", "FGFR4", "ERBB2")

compute_pct_avg_by_group <- function(object, genes, group_var) {
  expr_df <- FetchData(object, vars = c(genes, group_var))
  colnames(expr_df)[colnames(expr_df) == group_var] <- "group"

  expr_df %>%
    pivot_longer(cols = all_of(genes), names_to = "gene", values_to = "expr") %>%
    group_by(group, gene) %>%
    summarise(pct_expressing = round(mean(expr > 0) * 100, 1),
              avg_expr = round(mean(expr), 3),
              n_cells = n(), .groups = "drop")
}

# =================================================================
# CORE: % of cells expressing (vs not) each gene, PER CLUSTER - the main
# question. Run on any object with a group_by column (seurat_clusters).
# =================================================================
analyze_pct_expressing_by_cluster <- function(object, genes, results_path, object_label,
                                                group_by = "seurat_clusters") {
  object <- JoinLayers(object)
  genes_present <- genes[genes %in% rownames(object)]

  by_cluster <- compute_pct_avg_by_group(object, genes_present, group_by)
  write.csv(by_cluster, paste0(results_path, "HER2markers_byCluster_", object_label, ".csv"), row.names = FALSE)

  plot_marker_dotplot(object, marker_groups = genes_present, results_path = results_path,
                       filename = paste0("Dotplot_HER2markers_byCluster_", object_label, ".png"),
                       group_by = group_by)

  # Explicit % expressing vs not, as a bar plot per cluster
  p <- ggplot(by_cluster, aes(x = group, y = pct_expressing, fill = gene)) +
    geom_col(position = "dodge") +
    theme_bw() + theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste0("% cells expressing, by cluster (", object_label, ")"),
         x = "Cluster", y = "% expressing")
  ggsave(paste0(results_path, "PctExpressing_byCluster_", object_label, ".png"), p,
         width = max(8, length(unique(by_cluster$group)) * 0.5), height = 5, dpi = 300, bg = "white")

  cat(paste0("\n [CORE] % expressing by cluster done for ", object_label, " \n"))
  return(by_cluster)
}

# --- Run on Tumoral (fine subtype clusters) ---
pct_tumoral <- analyze_pct_expressing_by_cluster(dwTumoral, genes_check, results_GEMX_GC_path,
                                                   object_label = "Tumoral")

# --- Run on the GLOBAL object (all compartments, Level 0 clusters) ---
pct_global <- analyze_pct_expressing_by_cluster(dwAnnotated, genes_check, results_GEMX_GC_path,
                                                  object_label = "Global")

dwTumoral_joined <- JoinLayers(dwTumoral)
dwAnnotated_joined <- JoinLayers(dwAnnotated)
genes_check <- genes_check[genes_check %in% rownames(dwTumoral_joined)]


# =================================================================
# A2. Same question (% expressing), but by sample and by PAM50 subtype -
#     only makes sense within Tumoral (PAM50 doesn't exist globally)
# =================================================================
by_sample <- compute_pct_avg_by_group(dwTumoral_joined, genes_check, "orig.ident")
write.csv(by_sample, paste0(results_GEMX_GC_path, "HER2markers_bySample.csv"), row.names = FALSE)
plot_marker_dotplot(dwTumoral_joined, marker_groups = genes_check, results_path = results_GEMX_GC_path,
                     filename = "Dotplot_HER2markers_bySample.png", group_by = "orig.ident")

p_sample <- ggplot(by_sample, aes(x = group, y = pct_expressing, fill = gene)) +
  geom_col(position = "dodge") +
  theme_bw() + theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "% cells expressing, by sample", x = "Sample", y = "% expressing")
ggsave(paste0(results_GEMX_GC_path, "PctExpressing_bySample.png"), p_sample,
       width = 10, height = 5, dpi = 300, bg = "white")

by_pam50 <- compute_pct_avg_by_group(dwTumoral_joined, genes_check, "PAM50_predicted")
write.csv(by_pam50, paste0(results_GEMX_GC_path, "HER2markers_byPAM50.csv"), row.names = FALSE)
plot_marker_dotplot(dwTumoral_joined, marker_groups = genes_check, results_path = results_GEMX_GC_path,
                   filename = "Dotplot_HER2markers_byPAM50.png", group_by = "PAM50_predicted")

p_sample <- ggplot(by_pam50, aes(x = group, y = pct_expressing, fill = gene)) +
  geom_col(position = "dodge") +
  theme_bw() + theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "% cells expressing, by PAM50 Predicted Subtype", x = "PAM50 Predicted Subtype", y = "% expressing")
ggsave(paste0(results_GEMX_GC_path, "PctExpressing_byPAM50.png"), p_sample,
       width = 10, height = 5, dpi = 300, bg = "white")

by_subtype <- compute_pct_avg_by_group(dwTumoral_joined, genes_check, "Subtype")
write.csv(by_subtype, paste0(results_GEMX_GC_path, "HER2markers_bySubtype.csv"), row.names = FALSE)
plot_marker_dotplot(dwTumoral_joined, marker_groups = genes_check, results_path = results_GEMX_GC_path,
                   filename = "Dotplot_HER2markers_bySubtype.png", group_by = "Subtype")

p_sample <- ggplot(by_subtype, aes(x = group, y = pct_expressing, fill = gene)) +
  geom_col(position = "dodge") +
  theme_bw() + theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "% cells expressing, by Subtype", x = "Subtype", y = "% expressing")
ggsave(paste0(results_GEMX_GC_path, "PctExpressing_bySubtype.png"), p_sample,
       width = 10, height = 5, dpi = 300, bg = "white")
cat("\n [A2] Sample / PAM50 tables and plots done (Tumoral only) \n")

# =================================================================
# A3. Cross-tabulation: sample x subtype (PAM50) and sample x cluster,
#     in ABSOLUTE (n cells expressing) and PROPORTIONAL (% of that
#     sample-group combination) terms
# =================================================================
cross_tab_gene_expression <- function(object, genes, row_var, col_var, results_path, label) {
  expr_df <- FetchData(object, vars = c(genes, row_var, col_var))
  colnames(expr_df)[colnames(expr_df) == row_var] <- "row_group"
  colnames(expr_df)[colnames(expr_df) == col_var] <- "col_group"

  cross_df <- expr_df %>%
    pivot_longer(cols = all_of(genes), names_to = "gene", values_to = "expr") %>%
    group_by(row_group, col_group, gene) %>%
    summarise(n_expressing = sum(expr > 0), n_total = n(),
              pct_expressing = round(n_expressing / n_total * 100, 1), .groups = "drop")

  write.csv(cross_df, paste0(results_path, "CrossTab_", label, ".csv"), row.names = FALSE)

  # Absolute (n cells expressing)
  p_abs <- ggplot(cross_df, aes(x = col_group, y = row_group, fill = n_expressing)) +
    geom_tile(color = "white") +
    geom_text(aes(label = n_expressing), size = 2.5) +
    facet_wrap(~gene) +
    scale_fill_gradient(low = "white", high = "firebrick", name = "n cells") +
    theme_minimal() +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1),
          plot.background = element_rect(fill = "white", color = NA)) +
    labs(title = paste0(label, " - absolute n cells expressing"), x = col_var, y = row_var)
  ggsave(paste0(results_path, "CrossTab_", label, "_absolute.png"), p_abs,
         width = 5 * length(genes), height = max(5, length(unique(cross_df$row_group)) * 0.3),
         dpi = 300, bg = "white", limitsize = FALSE)

  # Proportional (% of that sample-group cell population)
  p_pct <- ggplot(cross_df, aes(x = col_group, y = row_group, fill = pct_expressing)) +
    geom_tile(color = "white") +
    geom_text(aes(label = pct_expressing), size = 2.5) +
    facet_wrap(~gene) +
    scale_fill_gradient(low = "white", high = "steelblue", name = "% expressing") +
    theme_minimal() +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1),
          plot.background = element_rect(fill = "white", color = NA)) +
    labs(title = paste0(label, " - % of cells expressing"), x = col_var, y = row_var)
  ggsave(paste0(results_path, "CrossTab_", label, "_proportional.png"), p_pct,
         width = 5 * length(genes), height = max(5, length(unique(cross_df$row_group)) * 0.3),
         dpi = 300, bg = "white", limitsize = FALSE)

  return(cross_df)
}

# --- sample x cluster ---
crosstab_sample_cluster <- cross_tab_gene_expression(dwAnnotated_joined, genes_check,
                                                        row_var = "orig.ident", col_var = "seurat_clusters",
                                                        results_path = results_GEMX_GC_path,
                                                        label = "Sample_x_Cluster")
# --- sample x tumor cluster ---
crosstab_sample_cluster <- cross_tab_gene_expression(dwTumoral_joined, genes_check,
                                                        row_var = "orig.ident", col_var = "seurat_clusters",
                                                        results_path = results_GEMX_GC_path,
                                                        label = "Sample_x_TumorCluster")
# --- sample x celltype ---
crosstab_sample_subtype <- cross_tab_gene_expression(dwAnnotated_joined, genes_check,
                                                        row_var = "orig.ident", col_var = "celltype",
                                                        results_path = results_GEMX_GC_path,
                                                        label = "Sample_x_Celltype")
# --- sample x PAM50 subtype ---
crosstab_sample_subtype <- cross_tab_gene_expression(dwTumoral_joined, genes_check,
                                                        row_var = "orig.ident", col_var = "PAM50_predicted",
                                                        results_path = results_GEMX_GC_path,
                                                        label = "Sample_x_PAM50")

# --- sample x Subtype ---
crosstab_sample_subtype <- cross_tab_gene_expression(dwAnnotated_joined, genes_check,
                                                        row_var = "orig.ident", col_var = "Subtype",
                                                        results_path = results_GEMX_GC_path,
                                                        label = "Sample_x_Subtype")

# --- tumor cluster x Subtype ---
crosstab_sample_subtype <- cross_tab_gene_expression(dwTumoral_joined, genes_check,
                                                        row_var = "seurat_clusters", col_var = "Subtype",
                                                        results_path = results_GEMX_GC_path,
                                                        label = "TumorCluster_x_Subtype")

# =================================================================
# B. Relationship with QC / technical quality - is expression driven by
#    sequencing depth (more genes detected in general = more likely to
#    detect ANY gene, including these), or is it a real biological signal?
# =================================================================
qc_expr_df <- FetchData(dwTumoral_joined, vars = c(genes_check, "nFeature_RNA", "nCount_RNA", "percent.mt"))

qc_cor <- do.call(rbind, lapply(genes_check, function(g) {
  data.frame(
    gene = g,
    cor_nFeature = round(cor(qc_expr_df[[g]], qc_expr_df$nFeature_RNA, method = "spearman"), 3),
    cor_nCount   = round(cor(qc_expr_df[[g]], qc_expr_df$nCount_RNA, method = "spearman"), 3),
    cor_pctmt    = round(cor(qc_expr_df[[g]], qc_expr_df$percent.mt, method = "spearman"), 3)
  )
}))
write.csv(qc_cor, paste0(results_GEMX_GC_path, "HER2markers_QCcorrelation.csv"), row.names = FALSE)
cat("\n [B] QC correlation done - a high cor_nFeature suggests expression tracks\n",
    "     general library complexity rather than being a specific/independent signal \n")
print(qc_cor)

# Violin: nFeature_RNA split by expressing vs non-expressing, per gene
qc_long <- qc_expr_df %>%
  pivot_longer(cols = all_of(genes_check), names_to = "gene", values_to = "expr") %>%
  mutate(status = ifelse(expr > 0, "expressing", "non_expressing"))

p_qc <- ggplot(qc_long, aes(x = status, y = nFeature_RNA, fill = status)) +
  geom_violin(alpha = 0.6) +
  facet_wrap(~gene) +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = "none") +
  labs(title = "nFeature_RNA: cells expressing vs not, per gene")
ggsave(paste0(results_GEMX_GC_path, "HER2markers_nFeature_byExpressionStatus.png"), p_qc,
       width = 9, height = 5, dpi = 300, bg = "white")


# =================================================================
# C. Co-expression between the 3 genes (e.g. ERBB2/GRB7/FGFR4 sit close on
#    17q12 and are sometimes co-amplified - check if that shows up here)
# =================================================================
if (length(genes_check) >= 2) {
  cor_matrix <- cor(qc_expr_df[, genes_check], method = "spearman")
  write.csv(as.data.frame(cor_matrix), paste0(results_GEMX_GC_path, "HER2markers_coexpression_matrix.csv"))

  pairs_df <- combn(genes_check, 2, simplify = FALSE)
  plots <- lapply(pairs_df, function(pair) {
    ggplot(qc_expr_df, aes(x = .data[[pair[1]]], y = .data[[pair[2]]])) +
      geom_point(size = 0.4, alpha = 0.3) +
      geom_smooth(method = "lm", se = FALSE, color = "firebrick") +
      theme_bw() + theme(panel.grid = element_blank()) +
      labs(title = paste0(pair[1], " vs ", pair[2]))
  })
  combined_scatter <- patchwork::wrap_plots(plots, ncol = length(plots))
  ggsave(paste0(results_GEMX_GC_path, "HER2markers_coexpression_scatter.png"), combined_scatter,
         width = 5 * length(plots), height = 5, dpi = 300, bg = "white")
  cat("\n [C] Co-expression correlation matrix and scatterplots done \n")
  print(cor_matrix)
}

