#!/usr/bin/env bash
##
## Pipeline runner. Runs the steps listed in pipeline/steps.tsv in order, for the run named in
## config.yaml (or $SCRNA_CONFIG), and stops at the first checkpoint or error.
##
##   pipeline/run.sh                 run every step that is not done yet
##   pipeline/run.sh --from 06       start at step 06 (earlier steps must be done; sub-steps 13a/13b = step 13)
##   pipeline/run.sh --to 05         stop after step 05
##   pipeline/run.sh --only 08       run one step (all its sub-steps)
##   pipeline/run.sh --force ...     re-run steps even if they are marked done
##   pipeline/run.sh --status        show the status of every step and exit
##
## Status per sub-step: results/<run>/status/<tag>.{OK,FAILED,CHECKPOINT,OOM}; log per sub-step:
## results/<run>/logs/<tag>.log; wall time and peak RAM: results/<run>/timings.tsv.
## Exit codes: 0 all requested steps done | 10 stopped at a checkpoint (a decision is needed) | 1 error.
##
## Environment: SCRNA_CONDA_ENV (default seurat5) is activated when conda is available; otherwise
## Rscript must already be on PATH with the packages of envs/seurat5.yml.
##
set -uo pipefail
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$PIPELINE_DIR")"
cd "$PROJECT_DIR"

FROM=""; TO=""; ONLY=""; FORCE=0; STATUS=0
while [ $# -gt 0 ]; do
  case "$1" in
    --from) FROM=$2; shift 2 ;;
    --to) TO=$2; shift 2 ;;
    --only) ONLY=$2; shift 2 ;;
    --force) FORCE=1; shift ;;
    --status) STATUS=1; shift ;;
    -h|--help) sed -n '2,20p' "$0" | sed 's/^## \{0,1\}//'; exit 0 ;;
    *) echo "Unknown option: $1"; exit 1 ;;
  esac
done

if command -v conda >/dev/null 2>&1; then
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate "${SCRNA_CONDA_ENV:-seurat5}" || { echo "conda env ${SCRNA_CONDA_ENV:-seurat5} not found"; exit 1; }
fi
command -v Rscript >/dev/null || { echo "Rscript not on PATH"; exit 1; }
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1

# --- config values (read with R so bash never parses YAML) ---
cfg() {
  Rscript --vanilla -e "c <- yaml::read_yaml(Sys.getenv('SCRNA_CONFIG', '$PIPELINE_DIR/config.yaml')); v <- c; for (k in strsplit('$1', '.', fixed = TRUE)[[1]]) v <- v[[k]]; cat(if (is.null(v)) '$2' else v)"
}
RUN=$(cfg run_name "")
RESULTS=$(cfg results_dir results)
case "$RESULTS" in /*) ;; *) RESULTS="$PROJECT_DIR/$RESULTS" ;; esac
RUN_DIR="$RESULTS/$RUN"
[ -n "$RUN" ] || { echo "config.yaml: run_name is empty"; exit 1; }
mkdir -p "$RUN_DIR/status" "$RUN_DIR/logs"
TIMINGS="$RUN_DIR/timings.tsv"
[ -f "$TIMINGS" ] || echo -e "tag\tstatus\twall_s\tpeak_GB\tend" > "$TIMINGS"
SAMPLES=$(Rscript --vanilla -e "c <- yaml::read_yaml(Sys.getenv('SCRNA_CONFIG', '$PIPELINE_DIR/config.yaml')); f <- c\$samples; if (!grepl('^/', f)) f <- file.path('$PROJECT_DIR', f); cat(read.csv(f, strip.white = TRUE)\$sample)")

# --- one sub-step: run, time, classify the exit (OK / CHECKPOINT / OOM / FAILED) ---
run_tag() {  # run_tag <tag> <mem_max or -> <script> [args...]
  local tag=$1 mem=$2 script=$3; shift 3
  local st_dir="$RUN_DIR/status" log="$RUN_DIR/logs/$tag.log" tf; tf=$(mktemp)
  if [ $FORCE -eq 0 ] && [ -f "$st_dir/$tag.OK" ]; then echo "  = $tag already done"; rm -f "$tf"; return 0; fi
  rm -f "$st_dir/$tag".{OK,FAILED,CHECKPOINT,OOM}
  echo "  > $tag  ($(date '+%F %T'))"
  local cmd=(Rscript "$PIPELINE_DIR/steps/$script" "$@")
  if [ "$mem" != "-" ] && command -v systemd-run >/dev/null 2>&1; then
    cmd=(systemd-run --user --scope --quiet -p MemoryMax="$mem" -p MemorySwapMax=0 "${cmd[@]}")
  fi
  /usr/bin/time -f "%e %M" -o "$tf" "${cmd[@]}" > "$log" 2>&1; local rc=$?
  local wall kb; read -r wall kb < <(tail -1 "$tf"); rm -f "$tf"
  local st
  case $rc in
    0) st=OK ;;
    10) st=CHECKPOINT ;;
    137) st=OOM ;;
    *) grep -qi "cannot allocate\|out of memory\|Killed" "$log" && st=OOM || st=FAILED ;;
  esac
  date '+%F %T' > "$st_dir/$tag.$st"
  echo -e "$tag\t$st\t$wall\t$(awk -v k="${kb:-0}" 'BEGIN{printf "%.1f", k/1048576}')\t$(date '+%F %T')" >> "$TIMINGS"
  case $st in
    OK) echo "    ok (${wall}s, $(awk -v k="${kb:-0}" 'BEGIN{printf "%.1f", k/1048576}') GB)" ;;
    CHECKPOINT) echo "    CHECKPOINT - decision needed:"; sed -n '/=== CHECKPOINT ===/,/=================/p' "$log" | sed 's/^/    /' ;;
    *) echo "    $st - last lines of $log:"; tail -15 "$log" | sed 's/^/    /' ;;
  esac
  [ "$st" = OK ]
}
export -f run_tag; export RUN_DIR TIMINGS PIPELINE_DIR FORCE

# --- status table ---
if [ $STATUS -eq 1 ]; then
  echo "Run: $RUN ($RUN_DIR)"
  ls "$RUN_DIR/status" 2>/dev/null | sed 's/\.\([A-Z]*\)$/\t\1/' | sort | column -t
  exit 0
fi

# --- walk the manifest ---
# --from / --to / --only take the step number (01 ... 18); a sub-step tag like 13b counts as its step
FROM=${FROM:0:2}; TO=${TO:0:2}; ONLY=${ONLY:0:2}
in_range() {
  local id=${1:0:2}
  [ -n "$ONLY" ] && { [ "$id" = "$ONLY" ]; return; }
  { [ -z "$FROM" ] || [ "$id" \> "$FROM" ] || [ "$id" = "$FROM" ]; } && { [ -z "$TO" ] || [ "$id" \< "$TO" ] || [ "$id" = "$TO" ]; }
}

echo "Run: $RUN -> $RUN_DIR"
while IFS=$'\t' read -r id script kind args mem parallel_key; do
  case "$id" in ''|\#*) continue ;; esac
  in_range "$id" || continue
  [ "$args" = "-" ] && args=""
  if [ "$kind" = per_sample ]; then
    par=$(cfg "$parallel_key" 4); cores=$(cfg resources.copykat_cores 4)
    echo "[$id] $script $args  (per sample, $par in parallel)"
    printf '%s\n' $SAMPLES | xargs -P "$par" -I{} bash -c \
      "run_tag '${id}_{}' '$mem' '$script' $(echo "$args" | sed "s/{sample}/{}/g; s/{cores}/$cores/g")"
    n_bad=0; for s in $SAMPLES; do [ -f "$RUN_DIR/status/${id}_$s.OK" ] || n_bad=$((n_bad+1)); done
    if [ $n_bad -gt 0 ]; then
      ls "$RUN_DIR/status/${id}_"*.CHECKPOINT >/dev/null 2>&1 && exit 10
      echo "[$id] $n_bad sample(s) not OK - stopping"; exit 1
    fi
  else
    echo "[$id] $script $args"
    # shellcheck disable=SC2086
    run_tag "$id" "$mem" "$script" $args || { [ -f "$RUN_DIR/status/$id.CHECKPOINT" ] && exit 10; exit 1; }
  fi
done < "$PIPELINE_DIR/steps.tsv"
echo "Done: requested steps finished ($(date '+%F %T'))"
