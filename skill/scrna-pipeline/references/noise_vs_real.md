# Noise or real population?

Three kinds of artefact look like populations. Each leaves a different fingerprint; check them before naming
or dropping a cluster. None of them can be assumed present or absent in a new cohort.

| Artefact | Fingerprint | Evidence (where) |
|---|---|---|
| **Doublets / multiplets** | full programme of lineage A **plus** ≥ 3 specific genes of lineage B in the same cells; more genes and counts than either parent; can escape scDblFinder when one parent is RNA-poor or the pair is frequent | `cluster_evidence.R` columns `lin3_*` (two high blocks); step 17 `profile` (foreign genes per cell) and `joint` (are there intermediate cells between the parents?) |
| **Ambient RNA** | 1-2 foreign genes per cell, at low level; the foreign genes are the sample's soup (its most abundant cells); high DecontX contamination; the group's share follows the sample, not biology | `decontX` column; step 17 `ambient` (extra RNA vs the empty-droplet profile) and `summary` |
| **Low quality** | low genes and counts, no own markers or only mitochondrial / ribosomal / long nuclear transcripts (MALAT1, NEAT1); often one or a few samples | `genes` column, `top_expressed_genes.csv`, QC violins |

## Decision guide
1. Real cell types have **specific** markers expressed by most of their cells (high pct.1) and usually appear
   in several samples (tumour clusters are the exception: they are patient-specific).
2. Foreign genes in 3+ per cell, adding up on top of a full parent programme, and no intermediate cells between
   the parents in a joint map → doublet/multiplet group → `DROP` (state the % of cells and the foreign genes).
3. Foreign genes at 1-2 per cell with high DecontX → ambient: keep the cells under their real label.
4. Low-RNA groups with a recognisable faint programme are real cells sequenced shallowly: label them by that
   programme or `LowRNA`; drop only when there is no programme at all.
5. When unsure, keep the cluster with a provisional `*_Mixed` label and an evidence note, and run a
   sub-population study (step 17) on it. Dropping is irreversible downstream; a provisional label is not.
6. A name inherited from an earlier pass (or another run) is a hypothesis: re-check it against this cluster's
   own evidence before keeping it.

## What a continuum means
A real intermediate state (e.g. a differentiation path) shows cells spread between the two programmes, with
one rising as the other falls. Doublets show the full programme of one parent with the other added on top,
and no cells in between. Step 17 `joint` plots exactly this (own vs contaminant score per cell, and the share
of contaminant cells among each cell's neighbours).

## Running a sub-population study (step 17)
Config `subpopulations.<study>`: `celltypes` (labels zoomed together), `focus_celltype` (the suspicious label),
`compare_celltypes` (context and the real owners of the foreign genes), `contaminant_celltype` (+ its gene
set name), then after the zoom `suspect_groups` / `clean_groups` (zoom clusters). Marker lists: preferably a
panel in `resources/marker_panels.yaml` referenced from the config; `MANUAL_PANELS` at the top of the step is
for quick tests. Order: check → profile → zoom (grid, decision, zoom) → signatures → ambient per sample → joint
→ summary → fill `annotation/<study>.tsv` → step 18. check/profile/ambient/joint are only needed when a group
looks like doublets or contamination; a study can also just look for sub-types (zoom + signatures).
