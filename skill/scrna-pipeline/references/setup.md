# Setting up a run

## Environment (once per machine)
- `conda env create -f pipeline/envs/seurat5.yml`, then `Rscript pipeline/envs/install_extra.R` (copykat,
  BPCells, nnls: not on conda). The runner activates `$SCRNA_CONDA_ENV` (default `seurat5`).
- BPCells may live outside the env: `r_extra_libs` in config (default `envs/R_bpcells_lib`).
- Check: `Rscript -e 'for (p in c("Seurat","BPCells","harmony","celda","scDblFinder","progeny","copykat","ragg","yaml")) cat(p, requireNamespace(p, quietly=TRUE), "\n")'`

## Inputs
- Per sample, two CellRanger outputs: the **raw** matrix (all droplets; DecontX uses the empty droplets as the
  ambient profile) and the **filtered** matrix (called cells). GEM-X Flex: the filtered matrix lacks ~900 probes
  that Flex excludes; the step restricts the raw droplets to the filtered genes.
- `pipeline/samples.csv`: `sample,raw_dir,filtered_dir` (relative to the project root or absolute). Sample order
  fixes the cell-name suffix (`<barcode>_<index>`), so keep the order stable within a cohort.
- Optional clinical sheet (`clinical_sheet`), joined by `scRNAseq_ID == sample`; adds `Subtype` etc. Without it
  the subtype plots/tests are skipped.

## A new run
1. Pick a new `run_name` (one folder per run under `results/`; never point at an existing folder unless you
   are resuming that same run). For variants keep a copy of the config (`SCRNA_CONFIG=path/to/config.yaml`)
   rather than editing the shared one.
2. Leave data-dependent values `null` (resolutions, label lists of studies): the run will stop and ask.
   Copy a decision from a previous run only if the user asks for an identical re-run.
3. Dry check: `pipeline/run.sh --status` (creates the run folder), then `pipeline/run.sh --to 02`.
4. Small first: when testing changes, run on 2-3 downsampled samples (1,500 cells each) in a scratch results
   folder; the full cohort (16 samples) takes ~1 h to step 12 and CopyKAT hours.

## Resources (reference cohort: 16 samples, ~180k droplets called, ~145k cells after QC)
- 01 per sample: ~2-4 GB each, 16 in parallel. Merge and clustering: tens of GB (counts stay on disk with BPCells).
- 13 CopyKAT: 1-10 h per sample, up to ~45 GB each (`MemoryMax` per job), 5 in parallel by default.
- 14 CopyKAT analysis peaked at ~220 GB on the full cohort: run it alone.
- If other heavy jobs share the machine, lower `resources.max_parallel_*` instead of risking OOM kills.
