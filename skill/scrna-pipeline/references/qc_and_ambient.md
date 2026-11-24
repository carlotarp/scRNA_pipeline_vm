# QC and ambient RNA (step 01)

Read after step 01 finishes, before telling the user it is fine. Judge samples against the rest of their own
cohort: absolute "normal" values depend on tissue, preservation and chemistry.

## Where to look
- `01_qc_decontx/cell_filtering.csv` and `QC_funnel*.png`: cells per sample at each filter.
- `01_qc_decontx/decontx/decontx_contamination_bySample.csv` + violin: estimated ambient fraction per cell.
- `01_qc_decontx/per_sample/`: per-sample violins, mito vs genes, doublet scores.
- `01_qc_decontx/decontx/GenomeWide_*`: which genes DecontX removed most. Expected: the soup of each sample's
  most abundant cells. Unexpected genes there (e.g. markers of a rare population) deserve a note.

## Reading them
- A sample with much higher contamination than the others is not a reason to drop it, but its clusters will
  carry more foreign genes: keep that in mind when a mixed cluster comes mostly from it.
- scDblFinder removes part of the doublets only; heterotypic pairs of cell types with very different RNA
  content can pass. Step 17 is the tool to check a suspicious group later.
- "not in matrix" lines in logs list panel genes absent from the probe set / matrix: informational, but a
  panel that loses most of its genes needs replacements in `marker_panels.yaml`.

## When to stop and tell the user (also in autonomous mode)
- A sample loses a much larger share of its cells than the others, or keeps very few cells (< ~500).
- A sample's median genes per cell is far below the cohort's (e.g. < 40 % of the median): its cells may
  cluster by quality instead of biology.
- A sample's composition contradicts its clinical annotation (possible sample swap).
Report the numbers and the options (keep, exclude, re-run with other thresholds as a separate run); do not
change thresholds in the shared config on your own.
