#!/bin/bash -l
#SBATCH --job-name=nTE-elig
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=01:00:00
#SBATCH --output=vcf-eligibility-%j.out
#SBATCH --error=vcf-eligibility-%j.err

# Step 3: record which store sites are callable in the VCF and orientable by
# the ARGs (normalize_tes.vcf_eligibility). The artifact is shared by focal sets
# and controls, so both are drawn from the same eligible universe.
#
# Resources: the production artifact (job 38961639) took 9m23s on one CPU and
# peaked at 12.8 GB.
#
# Example:
#
#   sbatch --export=ALL,VCF=/path/variants.vcf.gz,\
#   STORE=results/age_interval_store,ANCESTRAL=results/ancestral_states,\
#   OUTPUT=results/vcf_eligibility \
#     slurm/run_vcf_eligibility.sbatch
#
# Variables:
#   VCF           filtered, biallelic, genome-wide VCF                [required]
#   STORE         interval store from slurm/run_interval_store.sbatch [required]
#   ANCESTRAL     table directory from slurm/run_ancestral_table.sbatch [required]
#   OUTPUT        eligibility directory; must not exist               [required]
#   MIN_CALLABLE  minimum callable individuals per site                   (20)
#   HETEROZYGOUS  error|missing for diploid calls                     (error)

set -euo pipefail
# sbatch copies the job script into /var/spool/slurmd, so BASH_SOURCE does not
# point at the checkout. SLURM_SUBMIT_DIR does; the BASH_SOURCE form is the
# fallback for running this script directly.
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}}"
: "${VCF:?set VCF to the filtered genome-wide VCF}"
: "${STORE:?set STORE to the interval store directory}"
: "${ANCESTRAL:?set ANCESTRAL to the ancestral-state table directory}"
: "${OUTPUT:?set OUTPUT to the eligibility directory to create}"
MIN_CALLABLE="${MIN_CALLABLE:-20}"
HETEROZYGOUS="${HETEROZYGOUS:-error}"
cd "$PROJECT"
# Batch jobs get a non-login shell, where `module` does not exist.
. "$PROJECT/slurm/slurm_conda_bootstrap.sh"

for required in "$VCF" "$STORE/metadata.json" "$ANCESTRAL/metadata.json"; do
    [[ -e "$required" ]] || { echo "missing input: $required" >&2; exit 1; }
done
[[ ! -e "$OUTPUT" ]] || { echo "OUTPUT already exists: $OUTPUT" >&2; exit 1; }

echo "vcf        $VCF"
echo "store      $STORE"
echo "ancestral  $ANCESTRAL"
echo "output     $OUTPUT"
date
python -m normalize_tes.vcf_eligibility \
    --vcf "$VCF" \
    --store "$STORE" \
    --ancestral-table "$ANCESTRAL" \
    --output "$OUTPUT" \
    --min-callable "$MIN_CALLABLE" \
    --heterozygous "$HETEROZYGOUS"
date
