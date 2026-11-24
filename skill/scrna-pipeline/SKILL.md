---
name: scrna-pipeline
description: Run, resume or extend the breast-tumour scRNA-seq pipeline in pipeline/ (10x GEM-X Flex; DecontX + QC, Harmony, lineage / leukocyte / stroma / tumour annotation, CopyKAT, sub-population studies) as an autonomous or semi-autonomous operator — choosing resolutions, writing the cluster-to-label dictionaries by hand from markers and top expressed genes, and telling noise (doublets, ambient RNA, low-quality cells) from real populations. Use this skill whenever someone asks to run or continue the pipeline, start a new run on new samples, annotate or re-annotate clusters, choose a clustering resolution, resolve a CHECKPOINT message from pipeline/run.sh, investigate a suspicious / mixed cluster, or check whether a population is real — even if they only say "sigue con el pipeline", "anota los clusters", "qué resolución pongo" or paste a checkpoint.
---

# scRNA-seq pipeline operator

You drive `pipeline/run.sh`. The code computes; **you make the judgement calls** at its checkpoints and
leave a written trail so a human can audit every decision. The pipeline never labels cells by itself, on
purpose: automatic score/majority rules can mislabel whole clusters on a single shared gene. Labels come from
reading each cluster's evidence, like a biologist, and nothing is assumed from previous cohorts.

## Layout you work with

| Path | What |
|---|---|
| `pipeline/config.yaml` | defaults (paths, QC, PCs, CopyKAT, sub-population studies). `null` = decided at a checkpoint |
| `pipeline/samples.csv` | sample, raw_dir (all droplets), filtered_dir (CellRanger called cells) |
| `pipeline/steps.tsv` + `pipeline/run.sh` | ordered steps; runner (`--status`, `--from NN`, `--to NN`, `--only NN`, `--force`) |
| `pipeline/decide.R` | records a decision: `Rscript pipeline/decide.R <key> <value> "<reason>"` (`--show` lists them) |
| `results/<run>/decisions.yaml` | decisions of this run (value, reason, who, when) |
| `results/<run>/annotation/<dict>.tsv` | dictionaries: cluster, n_cells, top_markers, **label**, **evidence** |
| `results/<run>/<NN_step>/` | figures and tables of each step; `logs/`, `status/`, `timings.tsv` |
| `pipeline/resources/marker_panels.yaml` | curated marker panels (edit here, not in code) |
| `pipeline/templates/example_annotation.tsv` | invented dictionary: **format example** for the evidence column |

Steps: 01 DecontX+QC (per sample, then merge) → 02 resolution evidence → 03 clustering → 04/05 lineage →
06/07 leukocytes → 08/09 stroma → 10/11/12 tumour → 13/14 CopyKAT → 15/16 tumour re-annotation with CNVs →
17/18 sub-population studies (e.g. `caf`). Details per step: `references/pipeline_steps.md`.

## Mode

- **Semi-autonomous (default).** At each checkpoint gather the evidence, prepare the proposal (filled
  dictionary or value + one-line reasons), show the user a short summary, and continue only after they agree.
  Write the proposed dictionary to the file anyway so they can edit it directly.
- **Autonomous** (only when the user says so: "auto", "hazlo todo tú", "no me preguntes"). Decide, record and
  continue. Still stop and ask when: a decision would DROP > 10 % of the cells of a step, a sample fails QC
  badly (see `references/qc_and_ambient.md`), CopyKAT would start (long and memory-hungry), or the
  machine lacks RAM. Flag every uncertain call in the run notes so the human can review them afterwards.

Mention the mode you are in at the start, and how to switch.

## The loop

1. **Orient.** `pipeline/run.sh --status`; read `results/<run>/RUN_NOTES.md` and `decisions.yaml` if the run
   exists. New run → `references/setup.md` (new `run_name`; never reuse or overwrite an existing results folder:
   past runs are the record of what was decided).
2. **Run** `pipeline/run.sh` and wait for it: it returns by itself at the first checkpoint, error or the end,
   so run it in the foreground when the pending steps take minutes. Use a background job only for steps that
   take hours (CopyKAT, a full cohort) and wait on its exit; do not schedule periodic check-ins, which add
   idle time without adding information. 01 and 13 are per-sample and parallel. Before heavy steps check free
   RAM (`free -g`) and other jobs; peak RAM per step is in `timings.tsv`.
3. **Exit 0** → everything requested is done. **Exit 1** → read `results/<run>/logs/<tag>.log`; fix the cause
   (config, inputs, resources) rather than editing step code mid-run; an `OOM` status means lower
   parallelism or `resources.*` in config. **Exit 10** → a checkpoint; the message says which key or dictionary.
4. **Resolve the checkpoint** with the matching reference (table below), record it, add an entry to RUN_NOTES.
5. **Report after each finished step** in 2-4 lines (cells kept, clusters, anything odd), then go back to 2.

| Checkpoint | Read | Evidence to look at |
|---|---|---|
| `clustering.resolution`, `subclustering.*.resolution`, `subpopulations.*.resolution` | `references/resolution.md` | `*/ResolutionGrid.png`, `cluster_sizes*_resolution.csv`, `SankeyPlot.html` |
| dictionary `lineage`, `leukocytes`, `stroma`, `tumour`, `tumour_lp`, `tumour_lowrna`, `<study>` | `references/annotation.md` + `references/noise_vs_real.md` | `scripts/cluster_evidence.R <dict>`, dot plots, `cluster_*_FindMarkers_filtered.csv`, `top_expressed_genes.csv` |
| `subclustering.tumour_lp/lowrna.clusters`, `tumour_copykat.normal_epithelium_celltypes` | `references/copykat.md` | `15_tumour_copykat_evidence/check/` |
| `subpopulations.<study>.*` (celltypes, focus, compare, contaminant, suspect/clean groups) | `references/noise_vs_real.md` | `17_subpopulation_evidence/<study>/` |

## Writing a dictionary (the core skill)

1. Run `Rscript pipeline/skill/scrna-pipeline/scripts/cluster_evidence.R <dict>`: one table with, per
   cluster, size, sample spread, genes, DecontX contamination, % of cells carrying ≥3 genes of each lineage,
   filtered markers with pct.1, top expressed genes and (tumour) PAM50 / scSubtype / phase.
2. Look at the dot plots of the step (they are images: open them) and the per-cluster marker CSVs when the
   table is not enough.
3. For each cluster decide a label from **its own** markers. Earlier runs give you vocabulary only — reuse a
   name when the markers support it, never because cells overlapped with an old label.
4. Fill `label` and `evidence` (the genes and facts that justify it, e.g. `CD8A CCL5 GZMK; 3 samples`).
   `DROP` removes a cluster (state why: doublets, low quality, contamination); in CopyKAT zooms `KEEP` keeps
   the cell's current label. Use only the label asked for — no extra hierarchy columns.
5. Re-run the runner; it validates that every cluster has a label and evidence.

When markers are ambiguous, say so in the evidence (`ambiguous: ...`) and, in semi mode, ask. A wrong
confident label propagates into every downstream analysis; an honest "unclear" gets fixed.

## RUN_NOTES.md

Keep `results/<run>/RUN_NOTES.md` (template: `assets/RUN_NOTES_template.md`): one section per checkpoint with
the decision, the evidence in one or two lines, alternatives considered, and who decided (user / agent-auto).
It is what a reviewer reads first; keep it factual and short.

## Ground rules

- Do not edit pipeline code while a step runs (Rscript reads files progressively). Code changes are a
  separate task, tested on a small run, and the user agrees to them first. If a step fails on a code bug,
  stop, describe the bug and the fix in RUN_NOTES and the reply, and leave the patch to the user.
- Work inside the project folder. If an optional input (e.g. the clinical sheet) is missing, run without
  it and say so; do not fetch files from elsewhere on the machine.
- Never hand-edit a results object. To change labels, edit the dictionary and re-run that annotate step
  (`run.sh --only NN --force`) and the steps after it.
- Changing an upstream decision (e.g. a resolution) invalidates the dictionaries built on it: rename the old
  dictionary file (`<dict>.tsv.old`) so a fresh template is written, and re-run from that step.
- Verify claims (yours and the user's) in the data before writing them in notes or reports.
- Figures: PNG with white background; one output per sample / cluster written as the loop goes.
- Reply in the user's language; code, notes, dictionaries and figures in English.
