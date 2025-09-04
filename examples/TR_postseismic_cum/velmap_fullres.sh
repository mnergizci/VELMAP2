#!/bin/bash

# === Get argument passed to sbatch ===
OUTDIRFILE=$1

if [ -z "$OUTDIRFILE" ]; then
  echo "Error: No output folder argument supplied."
  echo "Usage: sbatch vvelmap_fullresjob.sh <outdirfile>"
  exit 1
fi


# Run MATLAB with path setup and function call
matlab -nodisplay -nosplash -r "\
VELMAP_path = '/nfs/a1/eemne/VELMAP/'; \
addpath(genpath(fullfile(VELMAP_path, 'v2.1betap', 'pilib'))); \
addpath(genpath(fullfile(VELMAP_path, 'v2.1betap', 'static'))); \
addpath(genpath(fullfile(VELMAP_path, 'v2.1betap', 'mesh'))); \
addpath(genpath(fullfile(VELMAP_path, 'v2.1betap', 'plot'))); \
disp('VELMAP paths added successfully.'); \
post_process_full_res('$OUTDIRFILE'); \
exit;"
