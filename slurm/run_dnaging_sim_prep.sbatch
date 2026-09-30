#!/bin/bash
#SBATCH --job-name=nTE-simprep
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --array=1-10
#SBATCH --output=logs/simprep-%A_%a.out
#SBATCH --error=logs/simprep-%A_%a.err
#
# PhiTE steps 1-2 on one dnAging msprime/SINGER simulation replicate, plus the
# per-draw ancestral table that the TE-like arm of tools/sim_polarity_arms.py
# needs for agreeing-draw ages. One array task per replicate (1-10; only these
# ten have SINGER posteriors).
#
# Draws 50-99 are used, the post-burn-in set every dnAging analysis of these
# replicates uses (e.g. analyze_mutation_age_bias_by_frequency.py).
#
# Environment variables, overridable at submit time:
#   SIM_ROOT  dnAging msprime_variable_ne_error directory (read-only)
#   OUT_ROOT  output root inside this repository (default results/sim_dnaging)
#   FIRST_DRAW, LAST_DRAW  inclusive draw range (default 50, 99)

set -euo pipefail

PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:?submit from the repository root}}"
SIM_ROOT="${SIM_ROOT:-/quobyte/jrigrp/beil/logan_collab/dnAging/msprime_variable_ne_error}"
OUT_ROOT="${OUT_ROOT:-results/sim_dnaging}"
FIRST_DRAW="${FIRST_DRAW:-50}"
LAST_DRAW="${LAST_DRAW:-99}"
: "${TMPDIR:?Farm must provide per-job TMPDIR}"

rep=$(printf 'replicate_%03d' "${SLURM_ARRAY_TASK_ID:?run as an array job}")
draw_dir="$SIM_ROOT/argtest/$rep/out/combined"
out="$OUT_ROOT/$rep"

cd "$PROJECT"
. "$PROJECT/slurm/slurm_conda_bootstrap.sh"

trees=()
for d in $(seq "$FIRST_DRAW" "$LAST_DRAW"); do
    f="$draw_dir/run.combined.$d.tsz"
    [[ -f "$f" ]] || { echo "missing draw $f" >&2; exit 3; }
    trees+=("$f")
done
echo "$rep: ${#trees[@]} draws ($FIRST_DRAW-$LAST_DRAW)"
mkdir -p "$out"

echo "=== step 1: interval store ==="; date
python -m normalize_tes.build_snp_interval_store \
    "${trees[@]}" \
    --interval-store "$out/interval_store" \
    --min-usable-fraction 0.1 \
    --num-buckets 8 \
    --bucket-memory-gb 2 \
    --scratch-dir "$TMPDIR"

echo "=== step 2: ancestral-state table ==="; date
python -m normalize_tes.build_ancestral_states \
    --store "$out/interval_store" \
    --output "$out/ancestral_states" \
    "${trees[@]}"

echo "=== per-draw ancestral table ==="; date
python -m normalize_tes.build_draw_polarity \
    --store "$out/interval_store" \
    --output "$out/draw_polarity" \
    "${trees[@]}"
date
