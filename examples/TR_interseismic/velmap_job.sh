#!/bin/bash
#SBATCH --job-name=velmap
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --output=logs/velmap_%j.out
#SBATCH --error=logs/velmap_%j.err

set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <config.conf>"
  exit 2
fi

CFG="$1"
mkdir -p logs

# Resolve to absolute path (portable)
if [[ "$CFG" != /* ]]; then
  CFG="$(pwd)/$CFG"
fi

# (optional) set working directory for job outputs
cd "$(dirname "$CFG")/.." 2>/dev/null || true

module load matlab
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"

# Prefer -batch on modern MATLAB
matlab -batch "velmap_aire('$CFG')"
# If -batch not available, uncomment the fallback:
# matlab -nodisplay -nosplash -r "try; velmap_aire('$CFG'); catch e; disp(getReport(e,'extended')); exit(1); end; exit"

