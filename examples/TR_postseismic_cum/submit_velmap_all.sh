#!/bin/bash
set -euo pipefail

CONF_DIR="cum_conf"
JOB_SCRIPT="velmap_job.sh"   # the script that runs MATLAB for one config

mkdir -p logs

# loop over all .conf files in cum_conf/
for cfg in "$CONF_DIR"/*.conf; do
    # skip if no files
    [[ -e "$cfg" ]] || { echo "No configs in $CONF_DIR"; exit 1; }

    # job name = filename without extension
    jobname=$(basename "$cfg" .conf)

    echo "Submitting $cfg as job $jobname"

    # submit to slurm
    sbatch --job-name="$jobname" --output=logs/${jobname}.out --error=logs/${jobname}.err "$JOB_SCRIPT" "$cfg"
done

