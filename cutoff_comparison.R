##
##  Full Pipeline Comparison Across nFeature Cutoffs (7500 / 5500 / 3500)
##  Generalizes the fibrocyte-specific comparison to ALL celltypes, plus
##  PAM50 subtype stability in Tumoral - answers: is a stricter cutoff
##  removing noise, or losing real signal?
##

library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)

project_path <- "/home/usuario/PROJECTS/260724_victor_scRNA/"
wd <- paste0(project_path, "codes/scRNA_pipeline/")
setwd(wd)
results_path <- paste0(project_path, "results/")
results_COMPARE_path <- paste0(results_path, "GEMX/FullCutoffComparison/")

# Set paths to each version's FINAL annotated object (after PAM50/PROGENy step)
cutoff_runs <- list(
  "cutoff_7500" = paste0(results_path, "GEMX/CellAnnotation/7500/fully_annotated_data.rds"),
  "cutoff_5500" = paste0(results_path, "GEMX/CellAnnotation/5500/fully_annotated_data.rds"),
  "cutoff_3500" = paste0(results_path, "GEMX/CellAnnotation/3500/fully_annotated_data.rds")
)

objects_list <- lapply(cutoff_runs, readRDS)
names(objects_list) <- names(cutoff_runs)
cat("\n All 3 final annotated versions loaded \n")


# =================================================================
# 0. Normalize spelling/case + build a coarse "celltype_general" column
#    so comparisons across cutoffs (which have different clustering
#    resolution and inconsistent manual naming) are fair. `celltype`
#    itself is kept (only case/typo-fixed), for the fine-grained checks
#    that still need it (block E).
# =================================================================

# --- Collapse into coarse, cutoff-agnostic categories ---
collapse_celltype <- function(ct) {
  dplyr::case_when(
    grepl("^Tcell", ct, ignore.case = TRUE)                  ~ "Tcell",
    grepl("DC|TAM|Monocyte", ct, ignore.case = TRUE)          ~ "Myeloid_DC",
    grepl("^BCell$", ct, ignore.case = TRUE)                  ~ "Bcell",
    grepl("PlasmaBlast", ct, ignore.case = TRUE)              ~ "PlasmaBlast",
    grepl("Mast", ct, ignore.case = TRUE)                     ~ "Mast",
    grepl("Fibrocyte", ct, ignore.case = TRUE)                ~ "Fibrocyte",
    grepl("^Prolif", ct, ignore.case = TRUE)            ~ "Proliferative",
    grepl("CAF|Fibroblast|Pericyte|Adipocyte", ct, ignore.case = TRUE)  ~ "Stromal_fibroblast",
    grepl("Endothelial", ct, ignore.case = TRUE)              ~ "Endothelial",
    grepl("Lum", ct, ignore.case = TRUE)                      ~ "Lum",
    grepl("HER2+", ct, ignore.case = TRUE)                    ~ "Her2+",
    grepl("TNBC", ct, ignore.case = TRUE)                     ~ "TNBC",
    TRUE ~ ct
  )
}

objects_list <- lapply(objects_list, function(obj) {
  obj$celltype_general <- factor(collapse_celltype(as.character(obj$celltype)))
  obj
})

cat("\n [0] celltype normalized + celltype_general (coarse) created \n")
cat(" Mapping check per version:\n")
for (v in names(objects_list)) {
  cat(paste0("  ", v, ":\n"))
  print(table(objects_list[[v]]$celltype_general, objects_list[[v]]$celltype))
}


# =================================================================
# A. Overall cell retention and compartment/celltype composition
# =================================================================
composition_df <- do.call(rbind, lapply(names(objects_list), function(v) {
  obj <- objects_list[[v]]
  data.frame(cutoff_run = v, n_total = ncol(obj),
             table(obj$celltype_general, useNA = "ifany") %>% as.data.frame() %>%
               rename(celltype = Var1, n = Freq))
}))
composition_df$celltype <- as.character(composition_df$celltype)
composition_df$celltype[is.na(composition_df$celltype) | composition_df$celltype == ""] <- "NA"
composition_df$pct <- composition_df$n / composition_df$n_total * 100

write.csv(composition_df, paste0(results_COMPARE_path, "Composition_by_cutoff.csv"), row.names = FALSE)

p_composition <- ggplot(composition_df, aes(x = cutoff_run, y = pct, fill = celltype)) +
  geom_col() +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = "right") +
  labs(title = "Celltype composition (%) by cutoff", x = NULL, y = "% of cells")
ggsave(paste0(results_COMPARE_path, "Composition_by_cutoff.png"), p_composition,
       width = 9, height = 6, dpi = 300, bg = "white")

# Same info, easier to scan: one row per celltype, one column per cutoff
composition_wide <- composition_df %>%
  select(cutoff_run, celltype, pct) %>%
  pivot_wider(names_from = cutoff_run, values_from = pct, values_fill = 0)
write.csv(composition_wide, paste0(results_COMPARE_path, "Composition_wide.csv"), row.names = FALSE)
cat("\n [A] Composition comparison done \n")

p_composition_abs <- ggplot(composition_df, aes(x = celltype, y = n, fill = cutoff_run)) +
  geom_col(position = "dodge") +
  theme_bw() + theme(panel.grid = element_blank(),
                      axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Celltype composition (absolute counts) by cutoff", x = NULL, y = "Number of cells")
ggsave(paste0(results_COMPARE_path, "Composition_by_cutoff_absolute.png"), p_composition_abs,
       width = 10, height = 6, dpi = 300, bg = "white")

# =================================================================
# B + C. Per-celltype barcode overlap (Jaccard) and label stability
#         (generalizing the fibrocyte-specific check to every celltype)
# =================================================================
all_celltypes <- sort(unique(unlist(lapply(objects_list, function(o) as.character(o$celltype_general)))))
all_celltypes <- all_celltypes[!is.na(all_celltypes) & all_celltypes != ""]

common_barcodes <- Reduce(intersect, lapply(objects_list, colnames))
cat(paste0("\n ", length(common_barcodes), " cells present in all 3 cutoffs \n"))

celltype_by_version <- data.frame(
  cell = common_barcodes,
  v7500 = as.character(objects_list[["cutoff_7500"]]$celltype_general[common_barcodes]),
  v5500 = as.character(objects_list[["cutoff_5500"]]$celltype_general[common_barcodes]),
  v3500 = as.character(objects_list[["cutoff_3500"]]$celltype_general[common_barcodes])
)

stability_summary <- lapply(all_celltypes, function(ct) {
  cells_ct <- lapply(objects_list, function(o) colnames(o)[!is.na(o$celltype_general) & o$celltype_general == ct])

  jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))

  ever_ct <- celltype_by_version[celltype_by_version$v7500 == ct |
                                  celltype_by_version$v5500 == ct |
                                  celltype_by_version$v3500 == ct, ]
  n_stable <- if (nrow(ever_ct) > 0) sum(ever_ct$v7500 == ct & ever_ct$v5500 == ct & ever_ct$v3500 == ct) else 0
  pct_stable <- if (nrow(ever_ct) > 0) round(n_stable / nrow(ever_ct) * 100, 1) else NA

  data.frame(
    celltype = ct,
    n_7500 = length(cells_ct[["cutoff_7500"]]),
    n_5500 = length(cells_ct[["cutoff_5500"]]),
    n_3500 = length(cells_ct[["cutoff_3500"]]),
    jaccard_7500_5500 = round(jaccard(cells_ct[["cutoff_7500"]], cells_ct[["cutoff_5500"]]), 3),
    jaccard_5500_3500 = round(jaccard(cells_ct[["cutoff_5500"]], cells_ct[["cutoff_3500"]]), 3),
    jaccard_7500_3500 = round(jaccard(cells_ct[["cutoff_7500"]], cells_ct[["cutoff_3500"]]), 3),
    pct_label_stable_allthree = pct_stable
  )
})
stability_summary <- do.call(rbind, stability_summary)
write.csv(stability_summary, paste0(results_COMPARE_path, "Celltype_stability_summary.csv"), row.names = FALSE)

p_stability <- ggplot(stability_summary, aes(x = reorder(celltype, pct_label_stable_allthree), y = pct_label_stable_allthree)) +
  geom_col(fill = "steelblue") + coord_flip() +
  theme_bw() + theme(panel.grid = element_blank()) +
  labs(title = "% of cells with a stable label across all 3 cutoffs, per celltype",
       x = NULL, y = "% stable")
ggsave(paste0(results_COMPARE_path, "Celltype_label_stability.png"), p_stability,
       width = 8, height = max(5, length(all_celltypes) * 0.3), dpi = 300, bg = "white")
cat("\n [B+C] Per-celltype stability computed - low % / low Jaccard = candidate noise or resolution-sensitive boundary \n")


# =================================================================
# B+C bis. Transition tileplots - for each pair of cutoffs, where do
#          cells "move to" when their celltype_general label changes?
#          Row-normalized (% of the SOURCE version's celltype that ends
#          up in each celltype of the TARGET version).
# =================================================================
plot_transition_tile <- function(df, from_col, to_col, from_label, to_label, results_path, filename) {
  trans_df <- df %>%
    filter(!is.na(.data[[from_col]]), !is.na(.data[[to_col]])) %>%
    count(.data[[from_col]], .data[[to_col]], name = "n") %>%
    group_by(.data[[from_col]]) %>%
    mutate(pct = n / sum(n) * 100) %>%
    ungroup()

  p <- ggplot(trans_df, aes(x = .data[[to_col]], y = .data[[from_col]], fill = pct)) +
    geom_tile(color = "white") +
    geom_text(aes(label = ifelse(n >= 10, n, "")), size = 2.8) +
    scale_fill_gradient(low = "white", high = "firebrick", name = "% of\nsource") +
    labs(title = paste0(from_label, " \u2192 ", to_label),
         x = to_label, y = from_label) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.background = element_rect(fill = "white", color = NA))

  ggsave(paste0(results_path, filename), p, width = 8, height = 7, dpi = 300, bg = "white")
}

plot_transition_tile(celltype_by_version, "v7500", "v5500", "cutoff 7500", "cutoff 5500",
                      results_COMPARE_path, "Transition_7500_to_5500.png")
plot_transition_tile(celltype_by_version, "v5500", "v3500", "cutoff 5500", "cutoff 3500",
                      results_COMPARE_path, "Transition_5500_to_3500.png")
plot_transition_tile(celltype_by_version, "v7500", "v3500", "cutoff 7500", "cutoff 3500",
                      results_COMPARE_path, "Transition_7500_to_3500.png")
cat("\n [B+C bis] 3 transition tileplots generated - off-diagonal cells show where labels move to \n")



# =================================================================
# D. QC / doublet profile of "version-specific" cells (present in only
#    ONE cutoff's celltype call) vs "core stable" cells (all 3 agree)
# =================================================================
qc_vars <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
if ("scDblFinder.score" %in% colnames(objects_list[["cutoff_7500"]]@meta.data)) qc_vars <- c(qc_vars, "scDblFinder.score")

ever_all <- celltype_by_version %>%
  mutate(status = case_when(
    v7500 == v5500 & v5500 == v3500 & !is.na(v7500) ~ "stable_all3",
    TRUE ~ "unstable"
  ))

qc_by_status <- objects_list[["cutoff_5500"]]@meta.data[ever_all$cell, qc_vars, drop = FALSE]
qc_by_status$status <- ever_all$status
qc_by_status$cell <- ever_all$cell

qc_summary <- qc_by_status %>%
  group_by(status) %>%
  summarise(across(all_of(qc_vars),
                    list(mean = \(x) mean(x, na.rm = TRUE),
                         median = \(x) median(x, na.rm = TRUE))))
write.csv(qc_summary, paste0(results_COMPARE_path, "QC_stable_vs_unstable.csv"), row.names = FALSE)

p_qc_status <- ggplot(qc_by_status, aes(x = status, y = nFeature_RNA, fill = status)) +
  geom_violin(alpha = 0.6) +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = "none") +
  labs(title = "nFeature_RNA: stable-label cells vs cutoff-sensitive cells")
ggsave(paste0(results_COMPARE_path, "QC_stable_vs_unstable_nFeature.png"), p_qc_status,
       width = 6, height = 5, dpi = 300, bg = "white")
cat("\n [D] QC profile of stable vs cutoff-sensitive cells done \n")


# =================================================================
# F. PAM50 subtype composition stability within Tumoral
# =================================================================
if (all(sapply(objects_list, function(o) any(grepl("^Lum$|^TNBC$|^Her2+", unique(o$celltype)))))) {
  pam50_df <- do.call(rbind, lapply(names(objects_list), function(v) {
    obj <- objects_list[[v]]
    obj_tumor <- obj@meta.data[obj$compartment == "Tumoral", ]
    data.frame(cutoff_run = v, table(obj_tumor$celltype) %>% as.data.frame() %>%
                 rename(pam50 = Var1, n = Freq))
  }))
  pam50_df <- pam50_df %>% group_by(cutoff_run) %>% mutate(pct = n / sum(n) * 100)
  write.csv(pam50_df, paste0(results_COMPARE_path, "PAM50_composition_by_cutoff.csv"), row.names = FALSE)

  p_pam50 <- ggplot(pam50_df, aes(x = cutoff_run, y = pct, fill = pam50)) +
    geom_col(position = "stack") +
    theme_bw() + theme(panel.grid = element_blank()) +
    labs(title = "PAM50 predicted subtype composition (%) within Tumoral, by cutoff", x = NULL, y = "% of Tumoral cells")
  ggsave(paste0(results_COMPARE_path, "PAM50_composition_by_cutoff.png"), p_pam50,
         width = 7, height = 5, dpi = 300, bg = "white")
  cat("\n [F] PAM50 composition comparison done \n")
} else {
  cat("\n [F] Skipped - PAM50 labels not found as celltype values in one or more versions \n")
}


# =================================================================
# Summary printout
# =================================================================
cat("\n\n ==== SUMMARY: celltypes with LOWEST label stability across cutoffs (top 5) ==== \n")
print(head(stability_summary[order(stability_summary$pct_label_stable_allthree), ], 5))

cat("\n ==== SUMMARY: celltypes with HIGHEST label stability across cutoffs (top 5) ==== \n")
print(head(stability_summary[order(-stability_summary$pct_label_stable_allthree), ], 5))

cat(paste("\n ---- FINISHED FULL CUTOFF COMPARISON ----
    Generated files:
      · Composition_by_cutoff.csv / Composition_wide.csv
      · Celltype_stability_summary.csv
      · QC_stable_vs_unstable.csv
      · Marker_purity_by_cutoff.csv
      · PAM50_composition_by_cutoff.csv
    Generated plots:
      · Composition_by_cutoff.png
      · Celltype_label_stability.png
      · QC_stable_vs_unstable_nFeature.png
      · Marker_purity_by_cutoff.png
      · PAM50_composition_by_cutoff.png
          "))