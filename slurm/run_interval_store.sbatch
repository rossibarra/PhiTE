#!/bin/bash -l
#SBATCH --job-name=nTE-store
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=interval-store-%j.out
#SBATCH --error=interval-store-%j.err

# Step 1: build the posterior age-interval store for every SNP and TE from the
# posterior ARG draws (normalize_tes.build_snp_interval_store).
#
# Resources: the production 75-draw store (31.2 million sites, 1.87 billion
# intervals) was most likely built by a 1-CPU job that took 2h05m and peaked at
# 34.1 GB. docs/INTERVAL_STORE_BENCHMARKS.md projected 48 GB and up to 16 h.
# The builder is single-threaded and I/O-bound, so more CPUs do not help. It
# writes about 23 GiB of bucket scratch to node-local $TMPDIR; Farm's per-job
# /tmp is large enough. Use `high`: a preempted build restarts from the top.
#
# Example:
#
#   sbatch --export=ALL,TREES="/path/posterior/*.tsz",\
#   CHROM_OFFSETS=/path/chrom_offsets.txt,OUTPUT=results/age_interval_store \
#     slurm/run_interval_store.sbatch
#
# Variables:
#   TREES          glob of posterior ARG draws (*.tsz); quote it      [required]
#   OUTPUT         interval store directory; must not exist           [required]
#   CHROM_OFFSETS  chromosome offsets, unless every ARG embeds them      (unset)
#   MIN_USABLE_FRACTION  minimum fraction of draws with a usable interval (0.1)

set -euo pipefail
# sbatch copies the job script into /var/spool/slurmd, so BASH_SOURCE does not
# point at the checkout. SLURM_SUBMIT_DIR does; the BASH_SOURCE form is the
# fallback for running this script directly.
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}}"
: "${TREES:?set TREES to a quoted glob of posterior ARG draws}"
: "${OUTPUT:?set OUTPUT to the interval store directory to create}"
CHROM_OFFSETS="${CHROM_OFFSETS:-}"
MIN_USABLE_FRACTION="${MIN_USABLE_FRACTION:-0.1}"
SCRATCH="${TMPDIR:?TMPDIR must point at node-local scratch}"
cd "$PROJECT"
# Batch jobs get a non-login shell, where `module` does not exist.
. "$PROJECT/slurm/slurm_conda_bootstrap.sh"

[[ ! -e "$OUTPUT" ]] || { echo "OUTPUT already exists: $OUTPUT" >&2; exit 1; }
# shellcheck disable=SC2206  # TREES is a glob and must expand here
trees=($TREES)
[[ -e "${trees[0]}" ]] || { echo "no files match TREES=$TREES" >&2; exit 1; }
offsets_option=()
if [[ -n "$CHROM_OFFSETS" ]]; then
    [[ -e "$CHROM_OFFSETS" ]] || { echo "missing input: $CHROM_OFFSETS" >&2; exit 1; }
    offsets_option=(--chrom-offsets "$CHROM_OFFSETS")
fi

echo "draws      ${#trees[@]}"
echo "output     $OUTPUT"
echo "scratch    $SCRATCH"
date
python -m normalize_tes.build_snp_interval_store "${trees[@]}" \
    --interval-store "$OUTPUT" \
    --scratch-dir "$SCRATCH" \
    --min-usable-fraction "$MIN_USABLE_FRACTION" \
    "${offsets_option[@]}"
date
