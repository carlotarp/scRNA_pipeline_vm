# Steps: what each one does and what you decide

| Step | Does | Checkpoint |
|---|---|---|
| 01 `sample <S>` | DecontX (quick per-sample clustering, 20 PCs, res 0.5, empty droplets as background, rounded counts) → QC: mt < 10 %, 100 < genes < 7500, 3-MAD outliers of log counts/genes per sample, scDblFinder singlets | — |
| 01 `merge` | merge, counts on disk (BPCells), normalise, Harmony per sample; DecontX genome-wide plots, QC funnel | — |
| 02 | elbow + resolution grid (0.2-0.5) on the whole object | — (evidence) |
| 03 | neighbours/UMAP (30 Harmony PCs, PCA-initialised UMAP), clusters at the chosen resolution, FindAllMarkers, heat map | `clustering.resolution` |
| 04 | lineage dot plots (leukocytes / stromal / vascular / epithelial panels), markers per cluster, template | — |
| 05 | lineage = Tumor / Leukocytes / Stromal (or DROP) | dictionary `lineage` |
| 06 / 08 | leukocyte / stromal subset re-integrated, grid, then markers + template at the chosen resolution | `subclustering.{leukocytes,stroma}.resolution` |
| 07 / 09 | cell types into the subset and back into the full object | dictionaries `leukocytes`, `stroma` |
| 10 | tumour subset, grid; PAM50-like scores, PROGENy, cell cycle, scSubtype | `subclustering.tumour.resolution` |
| 11 | tumour markers, validation plots (PAM50 vs clinical, scSubtype, PROGENy, phase), template | — |
| 12 | tumour cell types into the full object | dictionary `tumour` |
| 13 `prep` / `sample <S>` | CopyKAT per sample with T cells as diploid reference (own ≥ 50, else pooled) | — |
| 14 | CNA scores (CIN, PGA), prediction per cell, CNA UMAP | — |
| 15 `check` / `zoom lp` / `zoom lowrna` | CopyKAT state of every tumour cluster; zooms of chosen clusters to find normal epithelium and re-home low-RNA cells | `subclustering.tumour_*.clusters`, `.resolution`; dictionaries `tumour_lp`, `tumour_lowrna` |
| 16 | re-labels tumour cells from the zooms; normal epithelium gets lineage `Normal_Epithelium` | `tumour_copykat.normal_epithelium_celltypes` |
| 17 `<study> <mode>` | sub-population study (config `subpopulations.<study>`): check, profile, zoom, signatures, ambient per sample, joint, summary | study decisions, zoom resolution |
| 18 `<study>` | applies the study dictionary to the full object → `output_object` | dictionary `<study>` |

Why these parameters (reference run, decided with the user): 3 MAD instead of the old ±2 SD (2 SD discarded
~5 % of every sample, mostly large dividing/Lum_Her2 tumour cells; 3 vs 5 MAD barely differ); 30 PCs at every
level (the automatic elbow rule picked 12, which merged populations); PCA-initialised UMAP (spectral init
crashed at ~1M cells; only the picture changes); markers filtered at padj < 0.05 and pct.1 ≥ 33 % (smaller
fractions were dominated by ambient genes).

Objects (`results/<run>/objects/`): qc_merged → clustered → lineage_annotated → leuko_annotated →
notumor_annotated → tumour_scored → annotated (+ tumour_annotated) → copykat → annotated_copykat →
<study output>, `final` in the reference config.
