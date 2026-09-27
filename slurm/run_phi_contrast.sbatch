#!/bin/bash -l
#SBATCH --job-name=nTE-phicontrast
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=phi-contrast-%j.out
#SBATCH --error=phi-contrast-%j.err

# Formal contrast between two or more already-published phi-sfs-wasserstein-v2
# result directories (see normalize_tes/phi_contrast.py and
# PHI_SFS_WASSERSTEIN_CODING_PLAN.md section 8). This only reads each
# category's metadata, summary.csv, and null_z_scores.npy -- no VCF scan, no
# spectra recomputation -- so it is cheap and does not need the memory or time
# budget of run_phi_sfs.sbatch.
#
# Example:
#
#   sbatch --export=ALL,RESULTS="in_gene=results/phi_sfs/in_gene control=results/phi_sfs/control",\
#   OUTPUT=results/phi_contrast/in_gene_vs_control \
#     slurm/run_phi_contrast.sbatch
#
# Variables:
#   RESULTS           space-separated LABEL=DIR pairs, at least two   [required]
#   OUTPUT            results directory to publish; must not exist    [required]
#   SEED              base seed for the deterministic null pairing        (1002)
#   PAIRING_REPEATS   independent pairings for the P sensitivity range      (100)

set -euo pipefail
# sbatch copies the job script into /var/spool/slurmd, so BASH_SOURCE does not
# point at the checkout. SLURM_SUBMIT_DIR does; the BASH_SOURCE form is the
# fallback for running this script directly.
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}}"
: "${RESULTS:?set RESULTS to at least two space-separated LABEL=DIR pairs}"
: "${OUTPUT:?set OUTPUT to the results directory to publish}"
SEED="${SEED:-1002}"
PAIRING_REPEATS="${PAIRING_REPEATS:-100}"
cd "$PROJECT"
# Batch jobs get a non-login shell, where `module` does not exist.
. "$PROJECT/slurm/slurm_conda_bootstrap.sh"

# shellcheck disable=SC2206 # intentional word splitting of a space-separated list
result_entries=($RESULTS)
if [[ ${#result_entries[@]} -lt 2 ]]; then
    echo "RESULTS must list at least two LABEL=DIR entries: $RESULTS" >&2
    exit 1
fi

result_args=()
for entry in "${result_entries[@]}"; do
    label="${entry%%=*}"
    directory="${entry#*=}"
    metadata="$directory/metadata.json"
    [[ -e "$metadata" ]] || { echo "missing input: $metadata" >&2; exit 1; }
    echo "result     $label -> $directory"
    result_args+=(--result "$entry")
done

echo "seed       $SEED"
echo "repeats    $PAIRING_REPEATS"
date
python -m normalize_tes.phi_contrast \
    "${result_args[@]}" \
    --output "$OUTPUT" \
    --seed "$SEED" \
    --pairing-repeats "$PAIRING_REPEATS"
date
