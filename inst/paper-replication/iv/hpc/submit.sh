#!/bin/bash
# submit.sh - SLURM submission script for the IV-GMM simulation array.
#
# Each array task runs one (design x n x kx x q) cell for all n_reps
# replications. The default paper grid has 36 cells:
#   dgp1: 4 n values x 1 q value
#   other 4 designs: 4 n values x 2 q values
#
# Usage, from the directory where you want results written:
#   sbatch /path/to/submit.sh
#
# RdistNN must be installed in the R library visible to the compute nodes.

#SBATCH --job-name=iv_gmm_sim
#SBATCH --array=1-36
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --output=logs/cell_%a.out
#SBATCH --error=logs/cell_%a.err

# Uncomment and edit if RdistNN is installed in a non-default user library:
# export R_LIBS_USER="/path/to/your/R/library:$R_LIBS_USER"

export SIM_OUTDIR="${SLURM_SUBMIT_DIR}/results"

# Load R module if needed by your cluster:
# module load R/4.4.0

mkdir -p logs
mkdir -p "${SIM_OUTDIR}/cells"

TASK_SCRIPT=$(Rscript --vanilla -e \
  'cat(system.file("paper-replication/iv/hpc/task.R", package = "RdistNN"))')

if [ -z "$TASK_SCRIPT" ]; then
  echo "ERROR: Could not locate task.R via system.file(). Is RdistNN installed?" >&2
  exit 1
fi

echo "Running IV-GMM task $SLURM_ARRAY_TASK_ID via $TASK_SCRIPT"
Rscript --vanilla "$TASK_SCRIPT"

