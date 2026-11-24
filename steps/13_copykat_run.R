#!/usr/bin/env Rscript
##
## 13_copykat_run.R - CopyKAT (aneuploid vs diploid by inferred CNV), one process per sample
## Usage : Rscript pipeline/steps/13_copykat_run.R prep
##         Rscript pipeline/steps/13_copykat_run.R sample <S> <cores>
##   prep   : reference plan + one input file (counts + reference cells) per sample
##   sample : CopyKAT for one sample; the runner launches the samples concurrently
## Reads : objects/notumor_annotated.rds (prep); 13_copykat_run/_inputs/<S>.rds (sample)
## Writes: 13_copykat_run/copykat_reference_plan.csv, copykat_reference_cells.txt, _inputs/<S>.rds (prep)
##         13_copykat_run/<S>/<S>_copykat_*.txt + heatmap ; 13_copykat_run/_status/<S>.csv (sample, status, hours)
## Checkpoint: none
##
## Reference per sample: all its own T cells (copykat.ref_celltypes) when it has >= copykat.min_ref_tcells;
## otherwise a pooled reference of <= copykat.pooled_ref_per_group cells per sample x T-cell type from every
## sample (seed copykat.seed). copykat_reference_cells.txt = union of every reference cell, so step 14 keeps
## the reference cells of other samples out of each sample's calls.
## Resources (reference cohort, 16 samples, 4 cores each, 5 at a time, MemoryMax 45G per sample; peak RAM not measured):
## time grows with the sample's cells: 3.5k cells 0.3 h, 7k 1.1-1.6 h, 10-12k 3.8-5.8 h, 14-16k 6-10.6 h.
##
file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
PIPELINE_DIR <- normalizePath(file.path(dirname(file_arg), ".."))
source(file.path(PIPELINE_DIR, "R/setup.R"))
out_dir <- step_dir("13_copykat_run")

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || !args[1] %in% c("prep", "sample") || (args[1] == "sample" && length(args) < 3)) {
  stop("Usage: 13_copykat_run.R prep | sample <S> <cores>")
}
in_path <- step_dir("13_copykat_run", "_inputs")
st_path <- step_dir("13_copykat_run", "_status")

ref_groups  <- cfg_get("copykat.ref_celltypes")
min_own_ref <- cfg_get("copykat.min_ref_tcells")
ref_cap     <- cfg_get("copykat.pooled_ref_per_group")

if (args[1] == "prep") {
  o <- JoinLayers(read_run_object("notumor_annotated"))
  counts_all <- GetAssayData(o, assay = "RNA", layer = "counts")
  meta <- o@meta.data[, c("orig.ident", "celltype")]
  meta$celltype <- as.character(meta$celltype)
  rm(o); invisible(gc())

  tcells <- rownames(meta)[meta$celltype %in% ref_groups]
  set.seed(cfg_get("copykat.seed"))
  strata <- split(tcells, paste(meta[tcells, "orig.ident"], meta[tcells, "celltype"]))
  pooled <- unlist(lapply(strata, function(x) if (length(x) > ref_cap) sample(x, ref_cap) else x), use.names = FALSE)

  plan <- do.call(rbind, lapply(sort(unique(as.character(meta$orig.ident))), function(s) {
    in_s <- rownames(meta)[meta$orig.ident == s]
    own <- in_s[meta[in_s, "celltype"] %in% ref_groups]
    if (length(own) >= min_own_ref) {
      ref <- own; cells <- in_s; kind <- "own"
    } else {
      ref <- pooled; cells <- union(setdiff(in_s, pooled), pooled); kind <- "pooled"
    }
    # in-memory copy of the counts (CopyKAT cannot read BPCells matrices)
    saveRDS(list(counts = as(counts_all[, cells], "dgCMatrix"), ref = ref), file.path(in_path, paste0(s, ".rds")))
    cat("  input ->", s, "\n")
    data.frame(sample = s, t_cells = length(own), reference = kind, ref_cells = length(ref), cells = length(cells))
  }))
  write.csv(plan, paste0(out_dir, "copykat_reference_plan.csv"), row.names = FALSE)
  # pooled cells (used by the pooled samples) + every T cell of the own-reference samples
  ref_union <- unique(c(pooled, tcells[meta[tcells, "orig.ident"] %in% plan$sample[plan$reference == "own"]]))
  writeLines(ref_union, paste0(out_dir, "copykat_reference_cells.txt"))
  print(plan)
  cat("\n13_copykat_run prep done. Next: one 'sample <S> <cores>' job per sample, then 14_copykat_analysis.R\n")
  quit(save = "no")
}

# ---- sample mode ----
suppressPackageStartupMessages(library(copykat))
s <- args[2]; cores <- as.integer(args[3])
f_in <- file.path(in_path, paste0(s, ".rds"))
if (!file.exists(f_in)) stop("Missing ", f_in, ": run '13_copykat_run.R prep' first")
x <- readRDS(f_in)
setwd(step_dir("13_copykat_run", s))   # CopyKAT writes its files to the working directory
t0 <- Sys.time()
status <- tryCatch({
  copykat(rawmat = as.matrix(x$counts), id.type = "S", ngene.chr = cfg_get("copykat.ngene_chr"),
          win.size = cfg_get("copykat.win_size"), KS.cut = cfg_get("copykat.ks_cut"), sam.name = s,
          distance = cfg_get("copykat.distance"), norm.cell.names = x$ref, output.seg = FALSE, plot.genes = FALSE,
          genome = cfg_get("copykat.genome"), n.cores = cores)
  "DONE"
}, error = function(e) { cat("CopyKAT failed:", conditionMessage(e), "\n"); "FAILED" })
write.csv(data.frame(sample = s, status = status, cells = ncol(x$counts), ref_cells = length(x$ref),
                     hours = as.numeric(difftime(Sys.time(), t0, units = "hours"))),
          file.path(st_path, paste0(s, ".csv")), row.names = FALSE)
cat("\n13_copykat_run sample", s, status, "\n")
if (status == "FAILED") quit(save = "no", status = 1)
