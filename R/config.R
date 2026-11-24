##
## Configuration, paths and per-run decisions.
##
##   cfg        : pipeline/config.yaml (defaults shared by every run)
##   decisions  : results/<run>/decisions.yaml (choices made at the checkpoints of THIS run:
##                resolutions, excluded clusters...). A decision overrides the config default.
##   dictionaries: results/<run>/annotation/<step>.tsv (cluster -> label, see annotation.R)
##
## Every path in config.yaml is relative to the project root (the folder that contains pipeline/).
##

# ---------------------------------------------------------------
# Load config.yaml once; later calls return the cached copy.
# SCRNA_CONFIG overrides the config file (e.g. a second run with other defaults).
# ---------------------------------------------------------------
load_config <- function(config_file = Sys.getenv("SCRNA_CONFIG", file.path(PIPELINE_DIR, "config.yaml"))) {
  if (!is.null(.pipeline_env$cfg)) return(.pipeline_env$cfg)
  if (!file.exists(config_file)) stop("Config not found: ", config_file)
  cfg <- yaml::read_yaml(config_file)
  if (is.null(cfg$run_name) || cfg$run_name == "") stop("config.yaml: run_name is empty")
  .pipeline_env$cfg <- cfg
  cfg
}

# Project root = parent of pipeline/ (override with SCRNA_PROJECT)
project_root <- function() normalizePath(Sys.getenv("SCRNA_PROJECT", file.path(PIPELINE_DIR, "..")), mustWork = TRUE)

# Resolve a config path: absolute paths stay, relative ones hang from the project root
project_path <- function(...) {
  p <- file.path(...)
  if (grepl("^(/|~)", p)) path.expand(p) else file.path(project_root(), p)
}

# ---------------------------------------------------------------
# Read a nested config value by dotted key, e.g. cfg_get("qc.nmads").
# Stops when the key is missing (a typo must not silently become NULL).
# ---------------------------------------------------------------
cfg_get <- function(key, cfg = load_config()) {
  val <- cfg
  for (k in strsplit(key, ".", fixed = TRUE)[[1]]) {
    if (!is.list(val) || !k %in% names(val)) stop("config.yaml: missing key '", key, "'")
    val <- val[[k]]
  }
  val
}

# ---------------------------------------------------------------
# Run folders. Everything a run produces lives under results/<run_name>/.
#   run_dir()              -> results/<run>/
#   step_dir("2_clustering") -> results/<run>/2_clustering/ (created)
# ---------------------------------------------------------------
run_dir <- function() {
  d <- project_path(cfg_get("results_dir"), cfg_get("run_name"))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

step_dir <- function(...) {
  d <- file.path(run_dir(), ...)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  paste0(d, "/")          # trailing slash: the plot helpers paste file names onto it
}

# ---------------------------------------------------------------
# Sample sheet (pipeline/samples.csv): sample, raw_dir, filtered_dir.
# Checks that every listed folder exists before any work starts.
# ---------------------------------------------------------------
load_samples <- function(check_dirs = TRUE) {
  s <- read.csv(project_path(cfg_get("samples")), stringsAsFactors = FALSE, strip.white = TRUE)
  need <- c("sample", "raw_dir", "filtered_dir")
  if (!all(need %in% colnames(s))) stop("samples.csv needs columns: ", paste(need, collapse = ", "))
  if (anyDuplicated(s$sample)) stop("samples.csv: duplicated sample names")
  s$raw_dir <- vapply(s$raw_dir, project_path, "")
  s$filtered_dir <- vapply(s$filtered_dir, project_path, "")
  if (check_dirs) {
    miss <- c(s$raw_dir, s$filtered_dir)[!dir.exists(c(s$raw_dir, s$filtered_dir))]
    if (length(miss)) stop("samples.csv: missing folders:\n  ", paste(miss, collapse = "\n  "))
  }
  s
}

# ---------------------------------------------------------------
# Per-run decisions (results/<run>/decisions.yaml).
#
# get_decision("clustering.resolution") returns, in order:
#   1. the value recorded for this run in decisions.yaml
#   2. the config.yaml default, when it is not null
#   3. otherwise stops with a CHECKPOINT message: the step needs a human/agent choice.
# set_decision() records a choice (with who/why) so the run is reproducible.
# ---------------------------------------------------------------
decisions_file <- function() file.path(run_dir(), "decisions.yaml")

read_decisions <- function() {
  f <- decisions_file()
  if (file.exists(f)) yaml::read_yaml(f) else list()
}

get_decision <- function(key) {
  d <- read_decisions()
  if (!is.null(d[[key]]$value)) return(d[[key]]$value)
  default <- tryCatch(cfg_get(key), error = function(e) NULL)
  if (!is.null(default)) return(default)
  checkpoint_stop(sprintf("decision '%s' is not set. Inspect the evidence of the previous step and record it with:\n  Rscript pipeline/decide.R %s <value> \"<reason>\"", key, key))
}

set_decision <- function(key, value, reason = "", by = Sys.getenv("USER")) {
  d <- read_decisions()
  d[[key]] <- list(value = value, reason = reason, by = by, date = format(Sys.time(), "%Y-%m-%d %H:%M"))
  yaml::write_yaml(d, decisions_file())
  invisible(d)
}

# ---------------------------------------------------------------
# Stop at a checkpoint: signals a condition of class "checkpoint". Uncaught, the global
# handler registered in setup.R prints the message and exits with status 10, so the runner
# tells "waiting for a decision" apart from a real error (status 1).
# ---------------------------------------------------------------
checkpoint_stop <- function(msg) {
  stop(structure(class = c("checkpoint", "error", "condition"),
                 list(message = msg, call = NULL)))
}

checkpoint_exit_handler <- function(cnd) {
  message("\n=== CHECKPOINT ===\n", conditionMessage(cnd), "\n==================\n")
  if (!interactive()) quit(save = "no", status = 10)
}

# Path of a run object, e.g. run_object("clustered") -> results/<run>/objects/clustered.rds
run_object <- function(name) file.path(step_dir("objects"), paste0(name, ".rds"))

# Read a run object, with a clear message when the producing step has not run
read_run_object <- function(name) {
  f <- run_object(name)
  if (!file.exists(f)) stop("Missing ", f, ": run the step that produces '", name, "' first")
  cat("  Loading", basename(f), "\n")
  readRDS(f)
}
