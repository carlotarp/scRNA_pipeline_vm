# scRNA-seq Pipeline — Breast Cancer (GEMX-Flex)

Single-cell RNA-seq analysis pipeline for breast cancer samples sequenced with 10X Genomics GEM-X Flex. Covers the full workflow from raw CellRanger output to tumor subtype annotation, pathway activity scoring, and cell-cell communication.

---

## Pipeline overview

Scripts run strictly in order, **each exactly once**. Every step that needs a human decision (clustering resolution, cell type dictionaries, noise clusters) is split into a `_markers` script that produces the plots to look at, and an `_annotate` script whose manual choices sit at the top of the file. No script ever needs to be edited and re-run.

```
CellRanger output
      │
1a_qc.R                     Per-sample QC, doublet detection (scDblFinder), merge
1b_integrate.R              Harmony integration, clinical subtype annotation
1c_decontx.R                DecontX contamination correction
      │
2a_resolution.R             Elbow plot + clustering resolution grid      → choose `res`
2b_cluster.R                UMAP + clustering at the chosen resolution
2c_markers.R                FindAllMarkers per cluster (very slow)
      │
3a_lineage_markers.R        Lineage marker dotplots + cluster CSVs       → fill dictionary
3b_lineage_annotate.R       Broad lineage annotation (Tumor/Leukocytes/Stromal)
      │
3c_leukocytes_markers.R     Leukocyte subset, recluster, dotplots        → fill dictionary
3d_leukocytes_annotate.R    Leukocyte subtype annotation
      │
3e_stroma_markers.R         Stromal subset, recluster, dotplots          → fill dictionary
3f_stroma_annotate.R        Stromal subtype annotation
      │
3g_copykat.R                CopyKAT CNV ploidy + CIN/PGA + CNA-based UMAP
3h_tumor_scoring.R          Tumor subset: PAM50, PROGENy, DecoupleR, cell cycle, FindMarkers
      │
3i_tumor_markers.R          scSubtype scoring, PAM50/scSubtype validation,
      │                     differential scores, transition heatmaps      → fill dictionaries
3j_tumor_annotate.R         Tumor annotation + noise removal
      │                     → tumoral_annotated.rds / fully_annotated_data.rds
      │
4a_clinical_analysis.R      Celltype composition heatmaps + IGG scoring +
                            clinical covariate association (Subtype/pCR/IGG/TLS/TLS_binary):
                            composition heatmaps, chi-square/Cramer's V, proportion boxplots
                            (with Wilcoxon significance brackets), Spearman proportion correlation
```

> **Manual decision points** — all live at the top of the corresponding `_annotate` script (or `2b_cluster.R`), filled in from the previous script's plots:
>
> | Script | Variable to fill in | Based on |
> |---|---|---|
> | `2b_cluster.R` | `res` | `ResolutionGrid.png` (2a) |
> | `3b_lineage_annotate.R` | `lineage_clusters_annotated` | dotplots + cluster CSVs (3a) |
> | `3d_leukocytes_annotate.R` | `clusters_leuko_annotated` | dotplots + cluster CSVs (3c) |
> | `3f_stroma_annotate.R` | `clusters_Stroma_annotated` | dotplots + cluster CSVs (3e) |
> | `3j_tumor_annotate.R` | `clusters_tumor_annotated`, `noise_decontX_clusters`, `noise_tumor_clusters` | validation plots + cluster CSVs (3i) |


---

## File structure

```
scRNA_pipeline_vm/
│
├── Pipeline scripts
│   ├── 1a_qc.R
│   ├── 1b_integrate.R
│   ├── 1c_decontx.R
│   ├── 2a_resolution.R
│   ├── 2b_cluster.R
│   ├── 2c_markers.R
│   ├── 3a_lineage_markers.R
│   ├── 3b_lineage_annotate.R
│   ├── 3c_leukocytes_markers.R
│   ├── 3d_leukocytes_annotate.R
│   ├── 3e_stroma_markers.R
│   ├── 3f_stroma_annotate.R
│   ├── 3g_copykat.R
│   ├── 3h_tumor_scoring.R
│   ├── 3i_tumor_markers.R
│   ├── 3j_tumor_annotate.R
│   └── 4a_clinical_analysis.R
│
├── Plot functions (sourced by pipeline scripts)
│   ├── QC_plots.R          → used by 1a, 1b
│   ├── DECONTX_plots.R     → used by 1c
│   ├── CL_plots.R          → used by 2a–2c, 3g–3j (plot_dimplot, plot_featureplot, plot_composition, ...)
│   ├── CA_plots.R          → used by 3a–3h (plot_marker_dotplot)
│   ├── TCA_plots.R         → used by 3i, 3j, 4a (transition heatmaps, score boxplots, composition facets)
│   └── Clinical_plots.R    → used by 4a (generate_heatmap, clinical covariate association barplots,
│                            celltype proportion boxplots with significance brackets)
│
├── Gene signatures
│   └── scsubtype_signatures.R   Basal_SC / Her2E_SC / LumA_SC / LumB_SC gene lists (Wu et al. 2021)
│
├── Shared utilities
│   └── utils.R             color palettes, generate_lineage_subset(), find_markers_for_clusters(),
│                            compute_scsubtype_scores(), differential_scores_by_cluster(),
│                            proportion_correlation(), compute_celltype_proportions(),
│                            read_copykat_outputs() (3g: reads copyKAT's per-sample prediction +
│                            CNA files back from disk after the loop, instead of holding them in
│                            memory throughout)
│
├── Environment
│   └── seurat5.yml         conda env spec for the R/Seurat side of the pipeline
│
└── Auxiliary scripts
    └── rds_to_h3ad.R              Convert Seurat .rds → .h5ad for Python notebooks

```

---

## Dependencies

### R packages

| Package | Version | Purpose |
|---|---|---|
| Seurat | v5 | Core single-cell framework |
| harmony | — | Batch correction |
| celda | — | DecontX contamination correction |
| scDblFinder | — | Doublet detection |
| progeny | — | Pathway activity scoring |
| copykat | — | CNV-based ploidy prediction |
| dplyr, tibble, tidyr | — | Data wrangling |
| ggplot2, corrplot, gt | — | Visualisation |
| Matrix | — | Sparse matrix operations |
| reticulate | — | R–Python bridge (rds_to_h3ad.R) |

---

## How to run

Set `project_path` and `wd` at the top of each script before running. All scripts use the same path convention:

```r
project_path <- "/path/to/project/"
wd <- paste0(project_path, "codes/scRNA_pipeline/")
```

Run scripts sequentially following the pipeline order. Each script prints a summary of generated files on completion.

---

## Key outputs

| File | Generated by | Description |
|---|---|---|
| `merged_data.rds` | 1a | Merged, QC-filtered Seurat object |
| `sample_annotated_data.rds` | 1b | Harmony-integrated + clinical subtype metadata |
| `decontx_data.rds` | 1c | + `RNA_decontX` assay |
| `clustered_data.rds` | 2b | + UMAP and cluster assignments |
| `lineage_annotated_data.rds` | 3b | + `lineage` and `celltype` columns |
| `leukocytes.rds` | 3c / 3d | Reclustered leukocyte subset (annotated in 3d) |
| `leuko_annotated_data.rds` | 3d | Leukocyte fine-grained annotation |
| `Stroma.rds` | 3e / 3f | Reclustered stromal subset (annotated in 3f) |
| `notumor_annotated_data.rds` | 3f | Stromal fine-grained annotation |
| `copykat_annotated_data.rds` | 3g | + `copykat_prediction`, `CIN_score`, `PGA_score`, `copykat_cna` assay, `umap_cna` |
| `tumoral_scored.rds` | 3h / 3i | Tumor subset + PAM50 / PROGENy / DecoupleR assays, + `scSubtype_call` in 3i |
| `tumoral_annotated.rds` | 3j | Tumor subset with final `celltype` labels |
| `fully_annotated_data.rds` | 3j | Full object with all lineages annotated, noise removed |
