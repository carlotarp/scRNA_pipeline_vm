#!/usr/bin/env Rscript
## Record a checkpoint decision for the current run (results/<run>/decisions.yaml).
##   Rscript pipeline/decide.R <key> <value> "<reason>"
##   Rscript pipeline/decide.R --show
## Numbers are stored as numbers; comma-separated values as a list.
args <- commandArgs(trailingOnly = TRUE)
file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
PIPELINE_DIR <- dirname(normalizePath(file_arg))
.pipeline_env <- new.env()
source(file.path(PIPELINE_DIR, "R/config.R"))
if (length(args) == 1 && args[1] == "--show") { cat(yaml::as.yaml(read_decisions())); quit(save = "no") }
if (length(args) < 3) stop('Usage: Rscript pipeline/decide.R <key> <value> "<reason>"  |  --show')
parse_value <- function(v) {
  parts <- trimws(strsplit(v, ",", fixed = TRUE)[[1]])
  num <- suppressWarnings(as.numeric(parts))
  if (!anyNA(num)) num else parts
}
set_decision(args[1], parse_value(args[2]), reason = paste(args[-(1:2)], collapse = " "))
cat("Recorded", args[1], "=", args[2], "in", decisions_file(), "\n")
