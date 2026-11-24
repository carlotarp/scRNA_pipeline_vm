##
## CopyKAT (CNV-based aneuploid/diploid calls) helpers.
##

# ---------------------------------------------------------------
# Read a single sample's copyKAT output files back from disk.
#
# Args:
#   sample_path : directory containing the two output files (the per-sample
#                 working directory copykat() was run from)
#   sam_name    : the "sam.name" copykat() was run with (the file prefix)
#
# Returns list(
#   prediction = data.frame(cell.names, copykat.pred),
#   cna        = data.frame(chrom, chrompos, abspos, <one column per cell>)
# )
# ---------------------------------------------------------------
read_copykat_outputs <- function(sample_path, sam_name) {
  pred_file <- file.path(sample_path, paste0(sam_name, "_copykat_prediction.txt"))
  cna_file  <- file.path(sample_path, paste0(sam_name, "_copykat_final_results_bin_by_cell.txt"))

  prediction <- read.table(pred_file, header = TRUE, sep = "\t",
                           stringsAsFactors = FALSE, check.names = FALSE)
  cna <- read.table(cna_file, header = TRUE, sep = "\t",
                    stringsAsFactors = FALSE, check.names = FALSE)

  bin_cols <- setdiff(colnames(cna), c("chrom", "chrompos", "abspos"))
  colnames(cna)[colnames(cna) %in% bin_cols] <- sub("\\.", "-", bin_cols)

  list(prediction = prediction, cna = cna)
}
