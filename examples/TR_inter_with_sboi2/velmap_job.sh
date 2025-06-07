#!/bin/bash
#SBATCH --job-name=velmap_job          # Job name
#SBATCH --time=09:00:00                # Request runtime (hh:mm:ss)
#SBATCH --nodes=1                      # Number of nodes
#SBATCH --ntasks=1                     # One task (serial job)
#SBATCH --cpus-per-task=16             # Number of cores to use in MATLAB
#SBATCH --mem=128G                     # Adjust memory based on your needs
#SBATCH --output=velmap.out            # Stdout file
#SBATCH --error=velmap.err             # Stderr file

# Load MATLAB module (if required)
module load matlab

# Export number of threads for MATLAB parallel pool
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Run the MATLAB script
matlab -nodisplay -nosplash -r "run('velmap.m');exit;"

