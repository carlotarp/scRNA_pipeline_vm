# scRNA-seq pipeline (10x GEM-X Flex, breast tumours)

From CellRanger matrices to an annotated Seurat object: ambient-RNA removal (DecontX), cell QC, Harmony
integration, lineage / leukocyte / stroma / tumour annotation, tumour scoring (PAM50-like, scSubtype,
PROGENy, cell cycle), CopyKAT copy-number calls to separate tumour from normal epithelium, and generic
sub-population studies (e.g. CAF sub-types, doublets, contamination).

The pipeline is **config-driven and stops at checkpoints**. The code computes the evidence (resolution
grids, markers, dot plots, scores); a person, or an LLM agent using the included skill, chooses the
resolutions and writes the cluster-to-label dictionaries by hand. Nothing is labelled automatically, and
every decision is written down with its reason. The tools used are cited in [section 12](#12-tools-and-citations).

---

## Contents
1. [Install](#1-install)
2. [Inputs](#2-inputs)
3. [Configure a run](#3-configure-a-run)
4. [Run](#4-run)
5. [Checkpoints and decisions](#5-checkpoints-and-decisions)
6. [Annotation dictionaries](#6-annotation-dictionaries)
7. [Steps](#7-steps)
8. [Outputs](#8-outputs)
9. [Sub-population studies](#9-sub-population-studies)
10. [Running it with an agent](#10-running-it-with-an-agent)
11. [Repository layout](#11-repository-layout)
12. [Tools and citations](#12-tools-and-citations)

---

## 1. Install

Clone the repo **as a folder named `pipeline`** inside a project folder. Inputs (`data/`) and outputs
(`results/`) live next to it, outside the repo, and every relative path in the config is read from the
project folder.

```bash
mkdir my_project && cd my_project
git clone https://github.com/carlotarp/scRNA_pipeline_vm pipeline
```

```
my_project/
├── pipeline/      <- this repo
├── data/          <- your CellRanger outputs (not tracked)
└── results/       <- one folder per run (created by the pipeline)
```

R environment (R 4.4, Seurat 5, Harmony, celda/DecontX, scDblFinder, progeny, uwot…):

```bash
conda env create -f pipeline/envs/seurat5.yml
conda activate seurat5
Rscript pipeline/envs/install_extra.R     # copykat, BPCells, nnls (not on conda); prints a package check
```

`run.sh` activates the conda env `$SCRNA_CONDA_ENV` (default `seurat5`). Without conda, `Rscript` must
already be on the `PATH` with those packages. If BPCells is installed in a separate library, list it in
`r_extra_libs` in the config. Folders in that list that don't exist are ignored.

---

## 2. Inputs

Each sample needs **two CellRanger matrices**:

| | Used for |
|---|---|
| **raw** (all droplets, `raw_feature_bc_matrix`) | DecontX takes the ambient RNA profile from the empty droplets |
| **filtered** (called cells, `filtered_feature_bc_matrix`) | the cells that are analysed |

With GEM-X Flex, the filtered matrix lacks the probes that Flex excludes, so the raw droplets are
restricted to the genes of the filtered matrix.

**`pipeline/samples.csv`** (paths relative to the project folder, or absolute):

```csv
sample,raw_dir,filtered_dir
sample1,data/raw/sample1,data/filtered/sample1_filtered_feature_barcode_matrix
sample2,data/raw/sample2,data/filtered/sample2_filtered_feature_barcode_matrix
```

The row order sets the cell-name suffix (`<barcode>_<index>`), so keep it stable within a cohort.

**Clinical sheet (optional)**, set by `clinical_sheet` in the config. It is a CSV joined by
`scRNAseq_ID == sample` with the columns `Subtype` (e.g. `Lum`, `TNBC`, `Her2+`), `IGGscale_annot`, `pCR`
and `TLS`. Without it, `Subtype` is set to `Unknown` and the plots and tests by subtype are skipped.

---

## 3. Configure a run

All settings are in **`pipeline/config.yaml`**. The main keys:

| Key | What |
|---|---|
| `run_name` | name of the run; outputs go to `results/<run_name>/`. **Use a new name for each run**: past runs are the record of what was decided |
| `samples`, `clinical_sheet`, `bpcells_dir` | inputs; `bpcells_dir` holds the on-disk count matrices (objects keep pointers to it, so don't move it) |
| `decontx.*` | quick per-sample clustering that guides DecontX (20 PCs, res 0.5), count rounding |
| `qc.*` | mt < 10 %, 100 < genes < 7500, 3-MAD outliers of log counts / genes per sample, scDblFinder |
| `integration.*` | Harmony per sample (`orig.ident`), 30 PCs, PCA-initialised UMAP |
| `clustering.resolution_grid` | resolutions computed for the evidence (`[0.2, 0.3, 0.4, 0.5]`) |
| `markers.*` | marker filter: padj < 0.05 and expressed in ≥ 33 % of the cluster |
| `subclustering.*` | PCs and resolution per lineage sub-clustering and CopyKAT zoom |
| `copykat.*` | CopyKAT parameters and the T-cell labels used as diploid reference |
| `subpopulations.*` | sub-population studies (section 9) |
| `resources.*` | parallelism and memory limits |

**`null` means "decided at a checkpoint"**: resolutions, clusters to zoom into, label lists. Leave these
`null` in a new run; the pipeline stops and asks for them. Once a value is chosen, it is stored in the
run's `decisions.yaml`, which takes precedence over the config.

To run a variant without changing the shared config, point the runner at another file:

```bash
SCRNA_CONFIG=path/to/my_config.yaml pipeline/run.sh
```

Marker genes are not in the code: they are in **`pipeline/resources/marker_panels.yaml`** (lineage,
leukocyte, stroma, PAM50-like, epithelial states, lineage-specific gene sets, study panels). The scSubtype
signatures (Wu et al. 2021) are in `pipeline/resources/scsubtype_wu2021.yaml`.

---

## 4. Run

Run every command from the project folder:

```bash
pipeline/run.sh --status        # status of every step (creates the run folder)
pipeline/run.sh                 # run all pending steps; stops at the first checkpoint or error
pipeline/run.sh --to 05         # stop after step 05
pipeline/run.sh --from 06       # start at step 06 (earlier steps must be done)
pipeline/run.sh --only 08       # run a single step (all its sub-steps)
pipeline/run.sh --only 07 --force   # re-run a step that is already done
```

`--from`, `--to` and `--only` take the step number: `13` covers `13a` and `13b`.

**Exit codes:**

| Code | Meaning | What to do |
|---|---|---|
| `0` | all requested steps are done | — |
| `10` | **checkpoint**: a decision is needed | read the message, decide (section 5), run again |
| `1` | error | read `results/<run>/logs/<tag>.log`. An `OOM` status means it ran out of memory: lower `resources.*` |

Steps that are done are skipped, so after each decision just run `pipeline/run.sh` again.

**A typical first session:**

```bash
pipeline/run.sh --to 02          # QC + resolution evidence
# look at results/<run>/02_resolution/, choose a resolution:
Rscript pipeline/decide.R clustering.resolution 0.3 "clean lineage split; 0.4 splits TAM by sample"
pipeline/run.sh                  # clusters, writes the lineage template, stops: dictionary 'lineage'
# fill results/<run>/annotation/lineage.tsv (section 6)
pipeline/run.sh                  # continues to the next checkpoint
```

Test changes on 2–3 downsampled samples first, in a separate `run_name`.

---

## 5. Checkpoints and decisions

When a step needs a decision, it stops with exit code 10 and a message like:

```
=== CHECKPOINT ===
decision 'clustering.resolution' is not set. Inspect the evidence of the previous step and record it with:
  Rscript pipeline/decide.R clustering.resolution <value> "<reason>"
==================
```

Two kinds of checkpoint:

- **A value** (resolution, cluster list, label list). Record it with `decide.R`:
  ```bash
  Rscript pipeline/decide.R <key> <value> "<reason>"
  Rscript pipeline/decide.R subclustering.tumour_lp.clusters 1,4 "flat-CNA cells with LP programme"
  Rscript pipeline/decide.R --show            # list the decisions of the run
  ```
  Numbers are stored as numbers and comma-separated values as lists. Each decision keeps its value,
  reason, author and date in `results/<run>/decisions.yaml`.
- **A dictionary**: fill in the template the step wrote (section 6).

| Checkpoint | Evidence to look at |
|---|---|
| `clustering.resolution`, `subclustering.*.resolution`, `subpopulations.*.resolution` | `ResolutionGrid.png`, `cluster_sizes*_resolution.csv`, `SankeyPlot.html` of the step |
| dictionaries `lineage`, `leukocytes`, `stroma`, `tumour`, `tumour_lp`, `tumour_lowrna`, `<study>` | dot plots, `cluster_*_FindMarkers_filtered.csv`, `top_expressed_genes.csv` |
| `subclustering.tumour_lp/lowrna.clusters`, `tumour_copykat.normal_epithelium_celltypes` | `15_tumour_copykat_evidence/check/` |
| `subpopulations.<study>.*` | `17_subpopulation_evidence/<study>/` |

**Changing a decision already used** (for example a resolution) invalidates the dictionaries built on it.
Rename the old dictionary (`<dict>.tsv.old`) so a new template is written, and re-run from that step with
`--force`.

---

## 6. Annotation dictionaries

At each annotation level the markers step writes a template, `results/<run>/annotation/<dict>.tsv`:

| cluster | n_cells | top_markers | label | evidence |
|---|---|---|---|---|
| 0 | 5210 | CD3E,CD2,TRAC,IL7R | *Leukocytes* | *T cells (CD3E 92 %, TRAC, CD2); all samples* |
| 3 | 611 | COL1A2,C1QA,CD163,APOE | *Stroma_Mixed* | *ambiguous: fibroblast + macrophage genes …* |
| 4 | 248 | MALAT1,NEAT1,MT-CO1 | *DROP* | *low quality: median 310 genes …* |

Fill in **`label`** and **`evidence`** for every cluster (example: `pipeline/templates/example_annotation.tsv`).
The next run checks that the dictionary is complete and stops again if it isn't. An existing dictionary
is never overwritten.

Special labels:
- **`DROP`**: removes the cluster (doublets, low quality, contamination). Say why in the evidence; the
  removed cells are logged.
- **`KEEP`** (CopyKAT zooms only): the cells keep their current label.

Guidelines:
- Name each cluster from **its own** filtered markers and top expressed genes. Labels from earlier runs
  are vocabulary only.
- Scores (PAM50, scSubtype, PROGENy, phase, CopyKAT) are evidence to weigh, not labels.
- For mixed-profile clusters (two lineages) and low-RNA clusters, check whether they are doublets,
  ambient RNA, low quality or real cells before naming them. `skill/scrna-pipeline/references/annotation.md`
  and `noise_vs_real.md` describe the procedure. `skill/scrna-pipeline/scripts/cluster_evidence.R <dict>`
  summarises the evidence of each cluster in one table.

To change labels later, edit the dictionary and re-run that annotation step and the steps after it
(`run.sh --from NN --force`). Never edit a results object by hand.

---

## 7. Steps

The steps are listed in `pipeline/steps.tsv` and run in order.

| Step | Does | Stops for |
|---|---|---|
| 01 `sample` | per sample, in parallel: DecontX → QC (mt, genes, 3-MAD outliers) → scDblFinder singlets | — |
| 01 `merge` | merge (counts on disk, BPCells), normalise, Harmony; DecontX and QC summary plots | — |
| 02 | elbow plot + resolution grid on the whole object | — |
| 03 | neighbours, UMAP, clusters at the chosen resolution, markers, heat map | `clustering.resolution` |
| 04 / 05 | lineage dot plots + markers → Tumor / Leukocytes / Stromal | dictionary `lineage` |
| 06 / 07 | leukocyte subset: re-integration, grid, markers → cell types | resolution, dictionary `leukocytes` |
| 08 / 09 | stromal subset: same | resolution, dictionary `stroma` |
| 10 | tumour subset: grid; PAM50-like, scSubtype, PROGENy, cell-cycle scores | `subclustering.tumour.resolution` |
| 11 / 12 | tumour markers + validation plots → tumour cell types | dictionary `tumour` |
| 13 | CopyKAT per sample, with T cells as the diploid reference (own if ≥ 50, else pooled) | — |
| 14 | CNA scores (CIN, PGA), aneuploid/diploid per cell, CNA UMAP | — |
| 15 | CopyKAT state per tumour cluster; zooms to find normal epithelium (`lp`) and re-home low-RNA cells (`lowrna`) | clusters, resolutions, dictionaries `tumour_lp`, `tumour_lowrna` |
| 16 | re-labels tumour cells; normal epithelium gets lineage `Normal_Epithelium` | `tumour_copykat.normal_epithelium_celltypes` |
| 17 | sub-population study: check, profile, zoom, signatures, ambient per sample, joint, summary | study settings, zoom resolution |
| 18 | applies the study dictionary to the full object | dictionary `<study>` |

CopyKAT (steps 13–16) is by far the slowest and most memory-hungry part. Set `resources.*` (parallel
jobs, memory limit per job, cores) to fit your machine before launching it. Wall time and peak RAM of
every sub-step are written to `timings.tsv`, so a small test run tells you what to expect.

---

## 8. Outputs

Everything from a run goes to `results/<run_name>/`:

```
results/<run>/
├── decisions.yaml          # every checkpoint decision (value, reason, who, when)
├── annotation/<dict>.tsv   # the dictionaries
├── objects/*.rds           # Seurat objects after each stage
├── 01_qc_decontx/ … 18_*/  # figures (PNG) and tables of each step
├── logs/<tag>.log          # one log per sub-step
├── status/<tag>.OK|FAILED|CHECKPOINT|OOM
└── timings.tsv             # wall time and peak RAM per sub-step
```

Objects, in order: `qc_merged` → `clustered` → `lineage_annotated` → `leuko_annotated` →
`notumor_annotated` → `tumour_scored` → `annotated` → `copykat` → `annotated_copykat` → output of each
study (`final` in the default config). Final labels are in `celltype` and `lineage`.

Count matrices stay on disk in `bpcells_dir`, so the objects only work while that folder stays where it is.

---

## 9. Sub-population studies

Steps 17–18 zoom into a set of cell types to look for hidden sub-groups (doublets, ambient contamination,
published sub-types) and re-annotate them by hand. Each study is an entry under `subpopulations:` in the
config. The `caf` entry is the template:

```yaml
subpopulations:
  caf:
    input_object: annotated_copykat   # objects/<name>.rds to start from
    output_object: final              # objects/<name>.rds written by step 18
    celltypes: null                   # labels zoomed together
    focus_celltype: null              # label under suspicion
    compare_celltypes: null           # context cell types
    contaminant_celltype: null        # cell type whose RNA may contaminate the zoom
    contaminant_set: macrophage       # its gene set in the 'foreign' panel
    suspect_groups: null              # zoom clusters suspected of doublets / contamination
    clean_groups: null                # zoom clusters used as clean reference
    resolution: null
    panels:                           # gene lists: panels of resources/marker_panels.yaml
      own: caf_evidence.own
      foreign: caf_evidence.foreign
      signatures: caf_s1
      ...
```

Gene lists are best kept as panels in `marker_panels.yaml` and referenced here (`panel` or
`panel.block`). For quick tests, a hand-written list in `MANUAL_PANELS` at the top of
`steps/17_subpopulation_evidence.R` overrides them.

To add a study, copy the `caf` block under a new name and add its rows to `steps.tsv`: the same 17a–g
modes and step 18, with the study name in `args`.

---

## 10. Running it with an agent

`skill/scrna-pipeline/` is an agent skill. It teaches an LLM to operate the pipeline: choose
resolutions, write the dictionaries from markers and top expressed genes, tell noise (doublets, ambient
RNA, low-quality cells) from real populations, and keep a decision log (`results/<run>/RUN_NOTES.md`).

- **Claude Code**: link the skill from the project folder.
  ```bash
  mkdir -p .claude/skills
  ln -s ../../pipeline/skill/scrna-pipeline .claude/skills/scrna-pipeline
  ```
  Then ask, for example, *"run the pipeline on the samples in pipeline/samples.csv"* or paste a checkpoint
  message.
- **Other agents**: point them to `pipeline/AGENTS.md`, which sends them to `skill/scrna-pipeline/SKILL.md`.

Modes:
- **Semi-autonomous** (default): the agent proposes each decision and waits for your approval.
- **Autonomous** (when you ask for it): the agent decides on its own. It still stops before removing more
  than 10 % of a step's cells, when a sample fails QC badly, before CopyKAT, and when memory is short.

---

## 11. Repository layout

| Path | What |
|---|---|
| `run.sh` | runner (`--status`, `--from`, `--to`, `--only`, `--force`) |
| `steps.tsv` | ordered list of steps (script, per-sample or single, arguments, memory limit) |
| `steps/NN_*.R` | one script per step |
| `R/` | shared code: `setup.R` (entry point), `config.R` (config, decisions, checkpoints), `annotation.R`, `markers.R`, `processing.R`, `scoring.R`, `copykat.R`, `palettes.R`, `plots/` |
| `config.yaml`, `samples.csv` | run configuration and sample sheet |
| `decide.R` | records a checkpoint decision |
| `resources/` | marker panels and scSubtype signatures |
| `templates/example_annotation.tsv` | format example for the dictionaries (invented clusters) |
| `envs/` | conda environment and extra R packages |
| `skill/scrna-pipeline/` | agent skill (`SKILL.md`, `references/`, `scripts/cluster_evidence.R`, `assets/`) |
| `AGENTS.md` | entry point for non-Claude agents |

---

## 12. Tools and citations

If you use this pipeline, please cite the tools it relies on:

| Tool | Used for | Reference |
|---|---|---|
| Cell Ranger (10x Genomics) | GEM-X Flex matrices (input) | 10x Genomics, https://www.10xgenomics.com/support/software/cell-ranger |
| Seurat v5 | normalisation, PCA, neighbours, clustering, markers, module scores | Hao Y, et al. Dictionary learning for integrative, multimodal and scalable single-cell analysis. *Nat Biotechnol* 42, 293–304 (2024) |
| BPCells | on-disk count matrices | Parks B, Greenleaf WJ. Scalable high-performance single cell data analysis with BPCells. *bioRxiv* (2025); https://github.com/bnprks/BPCells |
| DecontX (celda) | ambient RNA removal | Yang S, et al. Decontamination of ambient RNA in single-cell RNA-seq with DecontX. *Genome Biol* 21, 57 (2020) |
| scDblFinder | doublet detection | Germain P-L, et al. Doublet identification in single-cell sequencing data using scDblFinder. *F1000Research* 10, 979 (2021) |
| Harmony | integration across samples | Korsunsky I, et al. Fast, sensitive and accurate integration of single-cell data with Harmony. *Nat Methods* 16, 1289–1296 (2019) |
| Louvain (Seurat `FindClusters`) | graph clustering | Blondel VD, et al. Fast unfolding of communities in large networks. *J Stat Mech* P10008 (2008) |
| UMAP (uwot) | 2-D embedding | McInnes L, Healy J, Melville J. UMAP: Uniform Manifold Approximation and Projection for dimension reduction. *arXiv* 1802.03426 (2018); uwot, https://github.com/jlmelville/uwot |
| PAM50 genes | PAM50-like module scores | Parker JS, et al. Supervised risk predictor of breast cancer based on intrinsic subtypes. *J Clin Oncol* 27, 1160–1167 (2009) |
| scSubtype | single-cell breast cancer subtype calls | Wu SZ, et al. A single-cell and spatially resolved atlas of human breast cancers. *Nat Genet* 53, 1334–1347 (2021) |
| PROGENy | pathway activity | Schubert M, et al. Perturbation-response genes reveal signaling footprints in cancer gene expression. *Nat Commun* 9, 20 (2018) |
| Cell-cycle genes (Seurat `cc.genes.updated.2019`) | S / G2M scores, phase | Tirosh I, et al. Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq. *Science* 352, 189–196 (2016) |
| CopyKAT | copy-number inference, aneuploid vs diploid cells | Gao R, et al. Delineating copy number and clonal substructure in human tumors from single-cell transcriptomes. *Nat Biotechnol* 39, 599–608 (2021) |
| CAF-S1 signatures | CAF sub-type scores and gates (CAF study) | Kieffer Y, et al. Single-cell analysis reveals fibroblast clusters linked to immunotherapy resistance in cancer. *Cancer Discov* 10, 1330–1351 (2020); Croizer H, et al. Deciphering the spatial landscape and plasticity of immunosuppressive fibroblasts in breast cancer. *Nat Commun* 15, 2806 (2024) |
| ggplot2, patchwork | figures | Wickham H. *ggplot2: Elegant Graphics for Data Analysis*. Springer (2016); patchwork, https://github.com/thomasp85/patchwork |
