#!/bin/bash

# === Get argument passed to sbatch ===
matlab -nodisplay -nosplash -r "\
VELMAP_path = '/nfs/a1/eemne/VELMAP/'; \
addpath(genpath(fullfile(VELMAP_path, 'v2.1betap', 'pilib'))); \
addpath(genpath(fullfile(VELMAP_path, 'v2.1betap', 'static'))); \
addpath(genpath(fullfile(VELMAP_path, 'v2.1betap', 'mesh'))); \
addpath(genpath(fullfile(VELMAP_path, 'v2.1betap', 'plot'))); \
disp('VELMAP paths added successfully.'); \
load insar.mat; \
plot_insar_fields(insar, 'insars');\
exit;"
