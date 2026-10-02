#!/bin/bash
#SBATCH --job-name=nTE-simarms
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=08:00:00
#SBATCH --output=logs/simarms-%j.out
#SBATCH --error=logs/simarms-%j.err
#
# Run tools/sim_polarity_arms.py over the dnAging replicates prepared by
# slurm/run_dnaging_sim_prep.sbatch. OUTPUT must be a new directory; any extra
# arguments are passed through (e.g. --replicates 1 --tests 2 for a smoke run).

set -euo pipefail
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:?submit from the repository root}}"
: "${OUTPUT:?set OUTPUT to a new results directory}"
cd "$PROJECT"
. "$PROJECT/slurm/slurm_conda_bootstrap.sh"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
date
python -m tools.sim_polarity_arms --output "$OUTPUT" "$@"
date
