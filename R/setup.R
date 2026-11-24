##
## Entry point of every pipeline step:
##   source(file.path(PIPELINE_DIR, "R/setup.R"))
## (PIPELINE_DIR is set by the step header; see steps/). Loads packages, the run config
## and every helper module, and fixes BLAS threads so parallel jobs do not oversubscribe.
##

if (!exists("PIPELINE_DIR")) PIPELINE_DIR <- Sys.getenv("SCRNA_PIPELINE", "pipeline")
PIPELINE_DIR <- normalizePath(PIPELINE_DIR, mustWork = TRUE)
.pipeline_env <- new.env()

source(file.path(PIPELINE_DIR, "R/config.R"))
globalCallingHandlers(checkpoint = checkpoint_exit_handler)   # outermost: a tryCatch inside a step still wins
# BPCells notes every on-disk -> dense conversion. The pipeline only densifies gene subsets
# (DotPlot / FetchData of a marker panel) or cell chunks (PROGENy), so the note is silenced.
globalCallingHandlers(rlang_message = function(m) {
  if (grepl("Converting to a dense matrix", conditionMessage(m), fixed = TRUE)) invokeRestart("muffleMessage")
})
cfg <- load_config()

# Extra R libraries (e.g. BPCells installed outside the conda env), before attaching packages
extra_libs <- vapply(cfg$r_extra_libs %||% character(0), project_path, "")
if (length(extra_libs)) .libPaths(c(extra_libs[dir.exists(extra_libs)], .libPaths()))

suppressPackageStartupMessages({
  library(Seurat)
  library(BPCells)       # counts on disk (cfg$bpcells_dir); objects keep pointers to those folders
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

blas <- as.character(cfg$resources$threads_blas %||% 1)
Sys.setenv(OPENBLAS_NUM_THREADS = blas, OMP_NUM_THREADS = blas, MKL_NUM_THREADS = blas)

for (f in c("palettes.R", "markers.R", "processing.R", "annotation.R", "scoring.R", "copykat.R",
            "plots/qc.R", "plots/decontx.R", "plots/clustering.R", "plots/markers.R", "plots/tumour.R")) {
  suppressPackageStartupMessages(source(file.path(PIPELINE_DIR, "R", f)))
}

cat(sprintf("\n[pipeline] run = %s | results = %s\n", cfg$run_name, run_dir()))
