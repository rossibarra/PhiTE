#!/bin/bash -l
#SBATCH --job-name=nTE-unpolSINGER
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --array=1-3
#SBATCH --output=logs/unpol-singer-%A_%a.out
#SBATCH --error=logs/unpol-singer-%A_%a.err
#
# Rerun SINGER on dnAging msprime_variable_ne_error replicates WITHOUT the true
# ancestral state, then the same argtest post-processing, so the posterior ARGs
# carry realistic polarity error. Everything else copies the original run's
# configs (singer/replicate_NNN/config.yaml, argtest/replicate_NNN/config.yaml):
# only the paths change, `polarised` becomes false (singer-snakemake then treats
# every site as unpolarised and SINGER polarises it during MCMC), and argtest's
# validation plots are skipped because they do not change the trees.
#
# The original dnAging files are only read; inputs are symlinked. Outputs:
#   $OUT_ROOT/singer/replicate_NNN   SINGER run (config, input links, output)
#   $OUT_ROOT/argtest/replicate_NNN  argtest run; draws in out/combined/
# after which slurm/run_dnaging_sim_prep.sbatch can build the PhiTE tables with
# SIM_ROOT=$OUT_ROOT.
#
# Reference-haplotype mode (REFERENCE_SEED set): mimic a reference-genome VCF.
# One modern haplotype is chosen per replicate (seed REFERENCE_SEED + replicate)
# and REF is relabelled to its allele (tools/reference_haplotype_vcf.py).
# singer-snakemake then runs unpolarised with its default settings.

set -euo pipefail
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:?submit from the repository root}}"
SIM_ROOT="${SIM_ROOT:-/quobyte/jrigrp/beil/logan_collab/dnAging/msprime_variable_ne_error}"
OUT_ROOT="${OUT_ROOT:-$PROJECT/results/sim_dnaging_unpolarised}"
SINGER_REPO="${SINGER_REPO:-/home/jri/src/singer-snakemake}"
ARGTEST_REPO="${ARGTEST_REPO:-/home/jri/src/argtest}"
REFERENCE_SEED="${REFERENCE_SEED:-}"
rep=$(printf 'replicate_%03d' "${SLURM_ARRAY_TASK_ID:?run as an array job}")

orig_singer="$SIM_ROOT/singer/$rep"
orig_argtest="$SIM_ROOT/argtest/$rep"
singer="$OUT_ROOT/singer/$rep"
argtest="$OUT_ROOT/argtest/$rep"
# A fresh submission refuses existing output. A job that Slurm requeued after
# preemption (SLURM_RESTART_COUNT > 0, e.g. on the `low` partition) resumes its
# own partial run instead: snakemake --rerun-incomplete picks up where it stopped.
restart="${SLURM_RESTART_COUNT:-0}"
if [[ "$restart" == 0 && ( -e "$singer/output" || -e "$argtest/out" ) ]]; then
    echo "$rep: output already exists under $OUT_ROOT; refusing to overwrite" >&2
    exit 4
fi
[[ "$restart" == 0 ]] || echo "$rep: requeued (restart $restart); resuming"

# --- SINGER, unpolarised --------------------------------------------------
mkdir -p "$singer/input"
for f in "$orig_singer"/input/*; do
    if [[ -n "$REFERENCE_SEED" && "$(basename "$f")" == chr1.vcf.gz ]]; then continue; fi
    ln -sfn "$f" "$singer/input/$(basename "$f")"
done
if [[ -n "$REFERENCE_SEED" ]]; then
    # The relabelled VCF is written in place, so a preempted write may be
    # partial; regenerate it (deterministically, from the same seed) on resume.
    rm -f "$singer/input/chr1.vcf.gz" "$singer/input/chr1.reference.json"
    . "$PROJECT/slurm/slurm_conda_bootstrap.sh"
    python -m tools.reference_haplotype_vcf \
        --vcf "$orig_singer/input/chr1.vcf.gz" \
        --output "$singer/input/chr1.vcf.gz" \
        --seed $(( REFERENCE_SEED + SLURM_ARRAY_TASK_ID ))
    conda deactivate
fi
sed -e "s#^input-dir:.*#input-dir: \"$singer/input\"#" \
    -e "s#^output-dir:.*#output-dir: \"$singer/output\"#" \
    -e "s#^polarised:.*#polarised: false#" \
    -e "1i # PhiTE rerun: copied from $orig_singer/config.yaml with polarised: false." \
    "$orig_singer/config.yaml" > "$singer/config.yaml"
grep -q '^polarised: false$' "$singer/config.yaml"
diff "$orig_singer/config.yaml" "$singer/config.yaml" || true

module load conda
conda activate singer-snakemake-env
echo "=== SINGER $rep ==="; date
snakemake \
    --snakefile "$SINGER_REPO/workflow/Snakefile" \
    --directory "$singer" \
    --configfile "$singer/config.yaml" \
    --cores "${SLURM_CPUS_PER_TASK:-32}" \
    --rerun-incomplete
conda deactivate

# --- argtest post-processing ----------------------------------------------
mkdir -p "$argtest/input/chr1"
for f in "$singer"/output/chr1/trees/*.tsz; do
    ln -sfn "$f" "$argtest/input/chr1/$(basename "$f")"
done
ln -sfn "$singer/output/chr1/chr1.mut_rate.p" "$argtest/input/chr1/chr1.mut_rate.p"
cp "$orig_argtest/input/all.hapmap.tsv" "$orig_argtest/input/sim.fai" "$argtest/input/"
ln -sfn "$ARGTEST_REPO/scripts" "$argtest/scripts"
sed -e "s#^root_dir:.*#root_dir: $argtest/input#" \
    -e "s#^hapmap:.*#hapmap: $argtest/input/all.hapmap.tsv#" \
    -e "s#^fai:.*#fai: $argtest/input/sim.fai#" \
    -e "s#^out_dir:.*#out_dir: $argtest/out#" \
    -e "s#^run_validation:.*#run_validation: false#" \
    "$orig_argtest/config.yaml" > "$argtest/config.yaml"
diff "$orig_argtest/config.yaml" "$argtest/config.yaml" || true

conda activate argtest
echo "=== argtest $rep ==="; date
snakemake \
    --snakefile "$ARGTEST_REPO/Snakefile" \
    --directory "$argtest" \
    --configfile "$argtest/config.yaml" \
    --cores "${SLURM_CPUS_PER_TASK:-32}" \
    --rerun-incomplete
ls "$argtest/out/combined" | head -3
date
