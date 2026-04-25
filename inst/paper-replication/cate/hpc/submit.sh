#!/bin/bash
# submit.sh - SLURM submission script for the CATE simulation array.
#
# Each array task runs one (design x n x kx) cell for all n_reps replications.
# The number of array tasks must equal:
#   length(designs) x length(n_vals) x length(kx_vals)
# Default paper grid: 4 x 4 x 1 = 16 cells -> --array=1-16
#
# Usage (from the directory where you want results written):
#   sbatch /path/to/submit.sh
# or copy this file, edit as needed, then:
#   sbatch submit.sh
#
# The script locates task.R via system.file() so RdistNN must be installed
# in the R library visible to the compute nodes. Set R_LIBS_USER below if
# the package lives in a non-default user library.

# SLURM directives
#SBATCH --job-name=cate_sim
#SBATCH --array=1-16              # one task per cell; update if grid changes
#SBATCH --cpus-per-task=8         # cores for within-task rep parallelism
#SBATCH --mem=8G                  # memory per task (adjust for large n)
#SBATCH --time=04:00:00           # wall time per task (adjust as needed)
#SBATCH --output=logs/cell_%a.out
#SBATCH --error=logs/cell_%a.err

# Environment

# Uncomment and edit if RdistNN is installed in a non-default user library:
# export R_LIBS_USER="/path/to/your/R/library:$R_LIBS_USER"

# Output directory for cell .rds files (relative to submission directory):
export SIM_OUTDIR="${SLURM_SUBMIT_DIR}/results"

# Load R module (name varies by cluster; adjust as needed):
# module load R/4.4.0

# Run

mkdir -p logs
mkdir -p "${SIM_OUTDIR}/cells"

TASK_SCRIPT=$(Rscript --vanilla -e \
  'cat(system.file("paper-replication/cate/hpc/task.R", package = "RdistNN"))')

if [ -z "$TASK_SCRIPT" ]; then
  echo "ERROR: Could not locate task.R via system.file(). Is RdistNN installed?" >&2
  exit 1
fi

echo "Running task $SLURM_ARRAY_TASK_ID via $TASK_SCRIPT"
Rscript --vanilla "$TASK_SCRIPT"
