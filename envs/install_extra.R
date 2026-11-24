## R packages not in envs/seurat5.yml (installed from GitHub / CRAN into the activated env, or into envs/R_bpcells_lib):
##   conda env create -f pipeline/envs/seurat5.yml && conda activate seurat5 && Rscript pipeline/envs/install_extra.R
## Tested versions: copykat 1.2.5 (navinlabcode@01924c0), BPCells 0.3.1, nnls 1.6, yaml 2.3.12
lib <- Sys.getenv("SCRNA_EXTRA_LIB", .libPaths()[1])
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", lib = lib)
if (!requireNamespace("nnls", quietly = TRUE)) install.packages("nnls", lib = lib)
if (!requireNamespace("copykat", quietly = TRUE)) remotes::install_github("navinlabcode/copykat@01924c0999bb660064da8d3e1561615ac6efb76f", lib = lib)
if (!requireNamespace("BPCells", quietly = TRUE)) remotes::install_github("bnprks/BPCells/r@v0.3.1", lib = lib)
for (p in c("Seurat", "BPCells", "harmony", "celda", "scDblFinder", "progeny", "copykat", "nnls", "RANN", "ragg", "yaml", "uwot"))
  cat(sprintf("%-12s %s\n", p, requireNamespace(p, quietly = TRUE)))
