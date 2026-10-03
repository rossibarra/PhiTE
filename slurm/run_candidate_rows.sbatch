#!/bin/bash -l
#SBATCH --job-name=nTE-cand
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=candidate-rows-%j.out
#SBATCH --error=candidate-rows-%j.err

# Step 4: build the control candidate pool (normalize_tes.build_candidate_rows):
# filtered SNPs, minus every TE (all categories) and the focal set A,
# restricted to eligible sites. Every TE category can share one file, because
# its A is already among the TEs; a SNP negative control needs its own, built
# with its A_POSITIONS.
#
# Resources: on the 75-draw store (31.2 million rows, 29.0 million candidates)
# a test run (job 39399950) took 3m11s on one CPU and peaked at 10.6 GB. Its
# output matched the V6 negative control's candidate file byte for byte.
#
# Example:
#
#   sbatch --export=ALL,STORE=results/age_interval_store,\
#   SNP_POSITIONS=/path/snp/all_snp.pos.txt,\
#   ALL_TE_POSITIONS=/path/te/all_te.pos.txt,\
#   A_POSITIONS=/path/te/in_gene.pos.txt,\
#   VCF_ELIGIBILITY=results/vcf_eligibility,OUTPUT=results/candidate_rows.npy \
#     slurm/run_candidate_rows.sbatch
#
# Variables:
#   STORE             interval store directory                        [required]
#   SNP_POSITIONS     filtered SNP positions allowed as controls      [required]
#   ALL_TE_POSITIONS  every TE position, all categories               [required]
#   VCF_ELIGIBILITY   directory from slurm/run_vcf_eligibility.sbatch [required]
#   OUTPUT            candidate .npy to create; must not exist        [required]
#   A_POSITIONS       focal set positions, also excluded               (unset)
#   MIN_RESOLVED_FRACTION  share of listed positions that must
#                     resolve to store rows; production uses 0.70     (0.70)

set -euo pipefail
# sbatch copies the job script into /var/spool/slurmd, so BASH_SOURCE does not
# point at the checkout. SLURM_SUBMIT_DIR does; the BASH_SOURCE form is the
# fallback for running this script directly.
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}}"
: "${STORE:?set STORE to the interval store directory}"
: "${SNP_POSITIONS:?set SNP_POSITIONS to the filtered SNP position file}"
: "${ALL_TE_POSITIONS:?set ALL_TE_POSITIONS to the all-TE position file}"
: "${VCF_ELIGIBILITY:?set VCF_ELIGIBILITY to the eligibility directory}"
: "${OUTPUT:?set OUTPUT to the candidate .npy to create}"
A_POSITIONS="${A_POSITIONS:-}"
MIN_RESOLVED_FRACTION="${MIN_RESOLVED_FRACTION:-0.70}"
cd "$PROJECT"
# Batch jobs get a non-login shell, where `module` does not exist.
. "$PROJECT/slurm/slurm_conda_bootstrap.sh"

exclude=("$ALL_TE_POSITIONS")
[[ -n "$A_POSITIONS" ]] && exclude+=("$A_POSITIONS")
for required in "$STORE/metadata.json" "$SNP_POSITIONS" \
                "$VCF_ELIGIBILITY/metadata.json" "${exclude[@]}"; do
    [[ -e "$required" ]] || { echo "missing input: $required" >&2; exit 1; }
done
[[ ! -e "$OUTPUT" ]] || { echo "OUTPUT already exists: $OUTPUT" >&2; exit 1; }

echo "store      $STORE"
echo "include    $SNP_POSITIONS"
echo "exclude    ${exclude[*]}"
echo "output     $OUTPUT"
date
python -m normalize_tes.build_candidate_rows \
    --store "$STORE" \
    --include-positions "$SNP_POSITIONS" \
    --exclude-positions "${exclude[@]}" \
    --vcf-eligibility "$VCF_ELIGIBILITY" \
    --output "$OUTPUT" \
    --min-resolved-fraction "$MIN_RESOLVED_FRACTION"
date
