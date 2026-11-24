# CopyKAT: tumour vs normal epithelium (steps 13-16)

CopyKAT infers copy-number profiles per cell against a diploid reference (the sample's own T cells when it has
enough, else a pooled reference; `copykat.*` in config). It is the tool here to separate normal epithelium from
tumour when markers alone cannot (normal and tumour epithelial programmes can overlap).

## Before launching (ask the user, also in autonomous mode)
Hours per sample and tens of GB each; step 14 holds every profile at once. Check `free -g`, other jobs and
`timings.tsv` of earlier runs.

## Reading `15_tumour_copykat_evidence/check/`
- **diploid ≠ normal.** A diploid call can come from low depth or from a tumour with few CNVs. Compare within
  each sample: counts of diploid vs aneuploid cells, and the correlation of each cell's CNA profile with the
  patient's aneuploid consensus (a diploid cell that still carries the patient's CNVs, attenuated, is tumour).
- **flat CNA** = CIN at or below the sample's T-cell reference quantile (`copykat.flat_quantile`). Only flat
  cells with a normal epithelial programme (luminal progenitor, mature luminal, basal/myoepithelial: see the
  `epithelial_states` panel) are candidates for normal epithelium. Their number may be zero: do not force it.
- **Which clusters to zoom** (`subclustering.tumour_lp.clusters`): tumour clusters that hold flat-CNA cells with
  a normal programme, any `*_candidate_normal` cluster, and clusters that mix patients with a programme close
  to a normal one (e.g. basal-like tumour clusters), where small normal groups tend to be absorbed. `subclustering.tumour_lowrna.clusters`: low-RNA tumour clusters, so their cells can be
  re-homed to the programme they faintly express.
- Small normal groups spread over patients can be dissolved by integration: zoom them without Harmony
  (`integrate: false`) and check whether the group survives with integration before trusting either view.

## Dictionaries `tumour_lp` / `tumour_lowrna`
- `KEEP` keeps the cell's current label (e.g. aneuploid single-patient tumour groups).
- A normal group gets a label of its own, named from its markers, and must be listed in
  `tumour_copykat.normal_epithelium_celltypes` so step 16 sets lineage `Normal_Epithelium`.
- Groups too small to characterise (a few dozen cells) are better folded into the closest label than made a
  class; say so in the evidence.
