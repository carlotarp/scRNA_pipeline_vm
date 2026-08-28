##
##  Single Cell Analysis Step 4c: copyKAT (per sample)
##  Same principle as the per-sample inferCNV script: ALL reference cells
##  (pooled across every sample) vs every NON-reference cell belonging to
##  THIS sample specifically. One copyKAT run per sample, same fixed
##  reference reused every time.
##

library(Seurat)
library(copykat)
library(dplyr)
library(ggplot2)

project_path <- "/home/usuario/PROJECTS/260724_victor_scRNA/"
results_path <- paste0(project_path, "results/")
results_GEMX_COPYKAT_path <- paste0(results_path, "GEMX/DecontX/CellAnnotation/CopyKAT/")

# Load data
dwAnnotated <- readRDS(paste0(results_path, "GEMX/DecontX/CellAnnotation/fully_annotated_data.rds"))
dwTumoral <- readRDS(paste0(results_path, "GEMX/DecontX/CellAnnotation/Tumor/tumoral_annotated.rds"))
cat("\n Data loaded \n")


# =================================================================
# MANUAL REFERENCE - same list you already use for inferCNV/infercna.
# Edit to match your $celltype values exactly.
# =================================================================
ref_groups <- c(
  "TCell_naive", "TCell_cyto", "TCell_ex"
)

dwAnnotated_joined <- JoinLayers(dwAnnotated)
ref_cells <- colnames(dwAnnotated_joined)[dwAnnotated_joined$celltype %in% ref_groups]
cat(paste0("\n Reference cells (pooled across all samples): ", length(ref_cells), " \n"))

samples <- sort(unique(as.character(dwAnnotated_joined$orig.ident[dwAnnotated_joined$lineage == "Tumor"])))
cat(paste0("\n Running copyKAT per sample (", length(samples), "): ", paste(samples, collapse = ", "), "\n"))


# =================================================================
# Run copyKAT once per sample: fixed reference vs this sample's
# non-reference cells.
# =================================================================
predictions_list <- list()
cna_list <- list()

for (s in samples) {

  cat(paste0("\n\n ==== Sample: ", s, " ==== \n"))
  results_sample_path <- paste0(results_GEMX_COPYKAT_path, s, "/")
  dir.create(results_sample_path, recursive = TRUE, showWarnings = FALSE)
  setwd(results_sample_path)  # copyKAT writes its own output files to the working directory

  obs_cells_s <- colnames(dwAnnotated_joined)[!(dwAnnotated_joined$celltype %in% ref_groups) &
                                                 dwAnnotated_joined$orig.ident == s]
  cells_this_run <- c(obs_cells_s, ref_cells)

  raw_counts_s <- as.matrix(GetAssayData(dwAnnotated_joined, assay = "RNA_decontX", layer = "counts")[, cells_this_run])

  copykat_result <- tryCatch({
    copykat(
      rawmat = raw_counts_s,
      id.type = "S",
      ngene.chr = 5,
      win.size = 25,
      KS.cut = 0.1,
      sam.name = s,
      distance = "euclidean",
      norm.cell.names = ref_cells,
      output.seg = FALSE,
      plot.genes = TRUE,
      genome = "hg20",
      n.cores = 10
    )
  }, error = function(e) {
    cat(paste0("  [ERROR] copyKAT failed for sample ", s, ": ", conditionMessage(e), "\n"))
    return(NULL)
  })

  if (is.null(copykat_result)) next

  pred_s <- copykat_result$prediction
  rownames(pred_s) <- pred_s$cell.names
  predictions_list[[s]] <- pred_s

  cna_list[[s]] <- copykat_result$CNAmat

  cat(paste0(" Sample ", s, " done -> ", results_sample_path, "\n"))
}

setwd(results_GEMX_COPYKAT_path)

# =================================================================
# Combine results across samples - ONE bind at the end
# =================================================================
all_predictions <- do.call(rbind, predictions_list)

# CNA matrices share the same chrom/pos columns (first 3) but different
# cells (columns) - merge by chrom/pos, keeping reference cells only once
cna_combined <- Reduce(function(x, y) {
  merge(x, y[, !colnames(y) %in% intersect(colnames(x), colnames(y))[-c(1:3)]],
        by = colnames(x)[1:3])
}, cna_list)

cna_values_combined <- cna_combined[, -c(1, 2, 3)]
copykat_signal <- setNames(apply(cna_values_combined, 2, function(x) mean(x^2, na.rm = TRUE)),
                             colnames(cna_values_combined))

cat(paste0("\n Combined predictions: ", nrow(all_predictions), " cells across ", length(predictions_list), " samples \n"))


# =================================================================
# Attach results to dwTumoral
# =================================================================
dwTumoral$copykat_prediction <- NA_character_
dwTumoral$copykat_signal <- NA_real_

matched_cells <- intersect(colnames(dwTumoral), rownames(all_predictions))
dwTumoral$copykat_prediction[matched_cells] <- all_predictions[matched_cells, "copykat.pred"]
dwTumoral$copykat_signal[matched_cells] <- copykat_signal[matched_cells]

dwTumoral$is_malignant_copykat <- dwTumoral$copykat_prediction == "aneuploid"

cat("\n copyKAT classification: \n")
print(table(dwTumoral$copykat_prediction, useNA = "ifany"))

saveRDS(dwTumoral, paste0(results_GEMX_COPYKAT_path, "Tumoral_with_copykat.rds"))
cat("\n copyKAT results attached to dwTumoral and saved \n")


# =================================================================
# Per-cluster summary
# =================================================================
copykat_by_cluster <- dwTumoral@meta.data %>%
  filter(!is.na(copykat_signal)) %>%
  group_by(seurat_clusters) %>%
  summarise(mean_copykat_signal = mean(copykat_signal),
            pct_aneuploid = mean(copykat_prediction == "aneuploid", na.rm = TRUE) * 100,
            pct_not_defined = mean(copykat_prediction == "not.defined", na.rm = TRUE) * 100,
            n_cells = n(), .groups = "drop") %>%
  arrange(desc(mean_copykat_signal))

write.csv(copykat_by_cluster, paste0(results_GEMX_COPYKAT_path, "CopyKAT_summary_byCluster.csv"), row.names = FALSE)
cat("\n Per-cluster copyKAT summary: \n")
print(copykat_by_cluster)


# =================================================================
# Differential analysis: copykat_signal per cluster (one-vs-rest)
# =================================================================
differential_scores_by_cluster <- function(object, score_cols, group_by = "seurat_clusters") {
  meta <- object@meta.data
  meta <- meta[!is.na(meta[[score_cols[1]]]), ]
  clusters <- sort(unique(as.character(meta[[group_by]])))

  results <- do.call(rbind, lapply(clusters, function(cl) {
    do.call(rbind, lapply(score_cols, function(sc) {
      in_group <- meta[[sc]][meta[[group_by]] == cl]
      out_group <- meta[[sc]][meta[[group_by]] != cl]

      wtest <- wilcox.test(in_group, out_group)

      data.frame(
        cluster = cl,
        score = sc,
        mean_in = mean(in_group),
        mean_out = mean(out_group),
        avg_diff = mean(in_group) - mean(out_group),
        p_val = wtest$p.value
      )
    }))
  }))

  results$p_val_adj <- p.adjust(results$p_val, method = "BH")
  results <- results[order(results$cluster, -abs(results$avg_diff)), ]
  return(results)
}

copykat_diff <- differential_scores_by_cluster(dwTumoral, "copykat_signal")
write.csv(copykat_diff, paste0(results_GEMX_COPYKAT_path, "Differential_copykat_byCluster.csv"), row.names = FALSE)


# --- Green-red barplot ---
p_copykat_diff <- ggplot(copykat_diff, aes(x = cluster, y = avg_diff, fill = avg_diff)) +
  geom_col() +
  scale_fill_gradient2(low = "forestgreen", mid = "white", high = "firebrick", midpoint = 0,
                        name = "Diff. vs\nrest") +
  coord_flip() +
  theme_bw() + theme(panel.grid = element_blank()) +
  labs(title = "Differential copyKAT CNA signal by cluster",
       x = "Cluster", y = "avg_diff (cluster vs rest)")
ggsave(paste0(results_GEMX_COPYKAT_path, "Barplot_Differential_copykat.png"), p_copykat_diff,
       width = 7, height = 8, dpi = 300, bg = "white")

# --- % aneuploid per cluster ---
p_pct_aneuploid <- ggplot(copykat_by_cluster, aes(x = reorder(seurat_clusters, pct_aneuploid), y = pct_aneuploid)) +
  geom_col(fill = "firebrick") +
  coord_flip() +
  theme_bw() + theme(panel.grid = element_blank()) +
  labs(title = "% cells classified as aneuploid (copyKAT), by cluster",
       x = "Cluster", y = "% aneuploid")
ggsave(paste0(results_GEMX_COPYKAT_path, "Barplot_pctAneuploid_byCluster.png"), p_pct_aneuploid,
       width = 7, height = 8, dpi = 300, bg = "white")


cat(paste("\n ---- FINISHED SCRIPT 4c: copyKAT (per sample) ----
    Generated files:
      · Tumoral_with_copykat.rds
      · CopyKAT_summary_byCluster.csv
      · Differential_copykat_byCluster.csv
      · <sample>/ subfolders with copyKAT's own per-sample output (heatmap, prediction, CNA matrix)
    Generated plots:
      · Barplot_Differential_copykat.png
      · Barplot_pctAneuploid_byCluster.png
          "))