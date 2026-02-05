#!/bin/bash

# === Get argument passed to sbatch ===
OUTDIRFILE=$1

#example of running for i in `ls -d output*`; do echo $i;sbatch_arie_highres velmap_fullres.sh "$i";done
if [ -z "$OUTDIRFILE" ]; then
  echo "Error: No output folder argument supplied."
  echo "Usage: sbatch vvelmap_fullresjob.sh <outdirfile>"
  exit 1
fi


# Run MATLAB with path setup and function call
matlab -nodisplay -nosplash -r "\
VELMAP_path = '/scratch/eemne/github/VELMAP/'; \
addpath(genpath(fullfile(VELMAP_path, 'v2.1betap', 'pilib'))); \
addpath(genpath(fullfile(VELMAP_path, 'v2.1betap', 'static'))); \
addpath(genpath(fullfile(VELMAP_path, 'v2.1betap', 'mesh'))); \
addpath(genpath(fullfile(VELMAP_path, 'v2.1betap', 'plot'))); \
disp('VELMAP paths added successfully.'); \
post_process_full_res('$OUTDIRFILE'); \
exit;"
