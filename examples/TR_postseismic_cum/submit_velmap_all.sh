#!/bin/bash
set -euo pipefail

# defaults
CONF_DIR="cum_conf"
JOB_SCRIPT="velmap_job.sh"

# parse args
while [[ $# -gt 0 ]]; do
  case "$1" in
    --sboi)
      CONF_DIR="cum_conf_sboi"
      shift
      ;;
    --job-script)
      JOB_SCRIPT="$2"
      shift 2
      ;;
    -h|--help)
      echo "Usage: $0 [--sboi] [--job-script SCRIPT]"
      echo "  --sboi            Use cum_conf_sboi folder instead of cum_conf"
      echo "  --job-script SCR  Path to job script (default: velmap_job.sh)"
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      exit 2
      ;;
  esac
done

mkdir -p logs

# loop over all .conf files
shopt -s nullglob
confs=("$CONF_DIR"/*.conf)
shopt -u nullglob

if (( ${#confs[@]} == 0 )); then
  echo "No configs in $CONF_DIR"
  exit 1
fi

for cfg in "${confs[@]}"; do
  jobname=$(basename "$cfg" .conf)
  echo "Submitting $cfg as job $jobname"
  sbatch \
    --job-name="$jobname" \
    --output="logs/${jobname}.out" \
    --error="logs/${jobname}.err" \
    "$JOB_SCRIPT" "$cfg"
done

