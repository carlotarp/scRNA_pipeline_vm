# Choosing a resolution

A resolution is good when every cluster can be named from its own markers, and populations the biology
separates are not merged. It is not a quality score: pick the lowest resolution that separates what matters
at that level; the next level (sub-clustering, zooms) handles finer structure. The right value depends on the
cohort (cell numbers, composition, depth), so judge each run on its own evidence.

## Evidence
- `ResolutionGrid.png`: UMAP coloured by clusters at each grid value.
- `cluster_sizes*_resolution.csv`: number and size of clusters per resolution (and top sample share).
- `SankeyPlot.html`: how clusters split from one resolution to the next. Stable splits (a parent dividing into
  children that stay put) are structure; clusters that reshuffle between resolutions are not.
- After choosing, the markers step shows whether a split was worth it: two clusters with the same markers and
  no distinct genes → too high.

## What each level should separate
- **Global (02/03)**: lineages, not cell types. Every cluster should clearly belong to one lineage (or be an
  obvious mixed / low-quality group to examine). Splitting large populations without new markers → too high.
- **Leukocyte and stroma sub-clustering (06, 08)**: the major cell types of the lineage. Increase only if a
  population with distinct markers is hidden inside another cluster; finer states belong to a sub-population
  study.
- **Tumour (10)**: tumour cells cluster by patient, so single-sample clusters are expected here. Choose the
  value where programmes shared across patients (e.g. cycling, interferon, luminal vs basal) appear as their own
  clusters, if they exist in the data.
- **Zooms (15, 17)**: small subsets; every cluster should have enough cells (tens at least) and own markers.

## Red flags
- Many very small clusters (< 0.5 % of the cells) with no specific markers → too high.
- In immune or stromal levels, a cluster made of one sample → batch or quality rather than a cell type, unless
  its markers define a known population (say so in the notes).
- A cluster defined only by low genes / high mitochondrial reads → low quality; see `annotation.md`.

Record: `Rscript pipeline/decide.R <key> <value> "<n clusters; what it separates; why not the next value>"`.
