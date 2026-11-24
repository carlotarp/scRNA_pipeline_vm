# Annotating clusters by hand

Every cohort is different: do not assume any population exists, or looks the way it did elsewhere. Each label
must be supported by the evidence of this run.

## Principles
- **Own evidence first.** Name each cluster from its own markers, top expressed genes and (tumour) scores.
  Labels of earlier runs or the reference dictionaries are vocabulary only: reuse a name when the markers agree,
  never because cells overlap with an old label.
- **No automatic rules**, not even as a first pass. Scores (PAM50, scSubtype, PROGENy, phase, CopyKAT) and
  marker-detection rates are evidence to weigh, not labels: a single shared gene can make a score point to the
  wrong cell type.
- **Filtered markers over raw.** Read `*_FindMarkers_filtered.csv` (padj and pct.1 thresholds of config): top
  log2FC genes expressed in a small fraction of a cluster usually come from ambient RNA or a few doublets.
- **Expressed vs specific.** Top expressed genes say what the cells *are*; markers say how they *differ* from
  the rest. A cluster whose top expressed genes belong to another lineage is suspect (see below).
- **Only the requested label.** One `label` per cluster; propose extra groupings in a sentence, don't add them.
- **Evidence = what a reviewer can check**: 3-6 genes (with % of cells when telling), the sample spread when
  relevant, and doubts (`ambiguous: ...`).

## Where the marker genes are
`pipeline/resources/marker_panels.yaml` holds the curated panels used by the dot plots of each level
(`lineage`, `leukocytes`, `stroma`, `tumour_pam50`, `epithelial_states`, `lineage_specific`, study panels).
Read the panel of the level you annotate; if a population you see is not covered, add its genes there (and
say so in RUN_NOTES) rather than relying on memory. Genes missing from the matrix are listed in the step log.

## Clusters with a mixed profile
A cluster that expresses the programmes of two lineages (e.g. fibroblast + macrophage genes, T cell +
epithelial genes). Possible causes: doublets/multiplets, ambient RNA, or a real intermediate state. Do not
name it as a cell type before telling them apart:
1. Quantify with `scripts/cluster_evidence.R` (`lin3_*`: % of cells with ≥ 3 specific genes of each lineage)
   and the DecontX column; compare with clean clusters of both parent lineages.
2. Read `noise_vs_real.md` for the fingerprints of each cause.
3. Then:
   - clear doublet fingerprint, small group → `DROP` with the numbers in the evidence;
   - ambient only (1-2 foreign genes per cell, high DecontX) → label by the real programme; the foreign genes
     are noise in its markers;
   - unclear, or a large group → keep it with a provisional label (`<parent>_Mixed`) and an evidence note, and
     resolve it in a sub-population study (step 17) before downstream analyses use it.

## Clusters with little RNA
A cluster with low genes/counts per cell and few or no own markers. Possible causes: shallowly sequenced real
cells, damaged/dying cells, or empty droplets with ambient RNA that passed QC.
1. Check `genes`, `top_expressed_genes.csv` and the sample spread. Compare its expressed genes with the other
   clusters of the level.
2. Then:
   - a recognisable programme at low level (it resembles one neighbouring cluster) → label it by that
     programme, or a provisional `LowRNA` if several programmes are mixed; deeper evidence (e.g. CNVs for
     tumour cells, a zoom) can re-home its cells later;
   - only mitochondrial / ribosomal / long nuclear transcripts (MALAT1, NEAT1) and no programme → low quality
     → `DROP`, stating it;
   - mostly the sample's soup genes and the highest DecontX values → ambient-dominated → `DROP`.
3. Never merge it silently into a neighbour: say in the evidence that the cluster is low-RNA.

## Level-specific notes
- **Lineage**: decide cycling clusters by their non-cycling genes. Normal epithelium falls with the tumour at
  this level; separating it needs CNV evidence (`copykat.md`), not markers alone.
- **Sub-lineages** (leukocytes, stroma): a cluster of one sample is a warning sign (batch, quality) unless its
  markers define a known cell type.
- **Tumour**: tumour clusters are largely patient-specific; name them by the programmes they express
  (luminal, HER2, basal, cycling, interferon, stress, ...) from their own markers plus the scores, and only
  when the genes are there.
- **Candidate normal epithelium**: an epithelial cluster with a normal programme (luminal progenitor, mature
  luminal, basal/myoepithelial; `epithelial_states` panel) that **lacks the tumour programme of its own
  sample** (compare with the tumour cells of the same patient: e.g. a marker present in almost all of that
  patient's tumour cells and absent here) is a candidate for normal epithelium. Label it as a candidate
  (`<programme>_candidate_normal`) with that comparison in the evidence; CNV evidence confirms or rejects it
  later (`copykat.md`). The absence of such a cluster proves nothing: in a large cohort, small normal groups
  are often absorbed by tumour clusters with a similar programme from other patients (typically basal-like
  tumours) and only separate in a zoom, so keep looking there.

## Examples
`pipeline/templates/example_annotation.tsv` shows the level of detail expected in the evidence column.
Its clusters are invented: use it for format, not content.
