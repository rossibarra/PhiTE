#!/bin/bash -l
#SBATCH --job-name=nTE-V6
#SBATCH --account=jrigrp
#SBATCH --partition=low
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=96G
#SBATCH --time=16:00:00
#SBATCH --output=logs/v6-%A_%a.out
#SBATCH --error=logs/v6-%A_%a.err
#
# V6 real-data negative control (docs/REMAINING_VALIDATION_PROPOSAL.md). One
# array task is one test, identified by SLURM_ARRAY_TASK_ID (0-299):
#
#   1. sample exactly M=4000 SNP focal sites from the candidate universe
#      (tools/v6_sample_focal.py), independently per test;
#   2. rebuild the candidate universe with those sites excluded
#      (normalize_tes.build_candidate_rows), so A never enters its own pool;
#   3. build the SNP target and match 110 disjoint sets, 3 restarts, through the
#      production launcher slurm/run_bootstrap_matching.sbatch;
#   4. Phi-SFS, -A SNP, all-mixture default, --max-null-replicates 99.
#
# Seeds are prespecified from the task ID: focal 610000+T, matching 710000+T,
# reference 810000+T. Every step is skipped when its output is complete, so a
# task requeued after preemption on `low` resumes (the matcher with --resume).
# A test whose bundle has fewer than 100 QC-passing sets makes Phi-SFS fail;
# the plan counts that test as a rejection. It is not rerun with a new seed.
#
# Submit from the repository root, for example the pilot:
#   sbatch --array=0-9 slurm/run_v6_negative_control.sbatch

set -euo pipefail
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:?submit from the repository root}}"
cd "$PROJECT"
T="${SLURM_ARRAY_TASK_ID:?run as an array job}"
COMBINED=/quobyte/jrigrp/beil/te_evo/singer_analysis/argtest/results/combined
STORE="${STORE:-$COMBINED/all_te_snp_age_interval_store}"
BASE_CANDIDATES="${BASE_CANDIDATES:-results/candidate-rows-allsnp-v0.8.npy}"
SNP_POSITIONS="${SNP_POSITIONS:-$COMBINED/snp_files/combined.all.chr.all.snp.pos.txt}"
ALL_TE_POSITIONS="${ALL_TE_POSITIONS:-$COMBINED/ages/b73_nam_max6mis_sv_genotyping_format_nodup.final.pos.txt}"
VCF_ELIGIBILITY="${VCF_ELIGIBILITY:-results/vcf_eligibility_75draw}"
VCF="${VCF:-/quobyte/jrigrp/beil/te_evo/singer_analysis/singer-snakemake/nam_te_vcf/combined.vcf}"
ANCESTRAL="${ANCESTRAL:-results/ancestral-rc2}"
V6_ROOT="${V6_ROOT:-results/v6}"
M="${M:-4000}"
SETS="${SETS:-110}"
NULLS="${NULLS:-99}"

dir="$V6_ROOT/test_$(printf '%04d' "$T")"
mkdir -p "$dir"
echo "V6 test $T in $dir (restart ${SLURM_RESTART_COUNT:-0})"; date
. "$PROJECT/slurm/slurm_conda_bootstrap.sh"

# 1. focal set
if [[ ! -s "$dir/focal.pos.txt" ]]; then
    rm -f "$dir/focal.pos.txt.json"
    python -m tools.v6_sample_focal --store "$STORE" \
        --candidate-rows "$BASE_CANDIDATES" --m "$M" \
        --seed $(( 610000 + T )) --output "$dir/focal.pos.txt"
fi

# 2. per-test candidate universe; the report is installed before the rows, so a
#    report without rows is an interrupted write
if [[ ! -s "$dir/candidates.npy" ]]; then
    rm -f "$dir/candidates.npy.json"
    python -m normalize_tes.build_candidate_rows \
        --store "$STORE" \
        --include-positions "$SNP_POSITIONS" \
        --exclude-positions "$ALL_TE_POSITIONS" "$dir/focal.pos.txt" \
        --vcf-eligibility "$VCF_ELIGIBILITY" \
        --output "$dir/candidates.npy" \
        --min-resolved-fraction 0.5
fi
conda deactivate

# 3. target and 110 matched sets through the production launcher
if ! python3 -c "import json,sys; sys.exit(0 if json.load(open('$dir/matches/metadata.json')).get('complete') else 1)" 2>/dev/null; then
    STORE="$STORE" TARGET="$dir/target" A_POSITIONS="$dir/focal.pos.txt" A_TYPE=SNP \
    OUTPUT="$dir/matches" WORK_DIR="$dir/work" CANDIDATE_ROWS="$dir/candidates.npy" \
    VCF_ELIGIBILITY="$VCF_ELIGIBILITY" MISSING_POSITION_POLICY=error \
    REPLICATES="$SETS" RESTARTS=3 SEED=$(( 710000 + T )) \
        bash "$PROJECT/slurm/run_bootstrap_matching.sbatch"
fi

# 4. Phi-SFS with exactly NULLS nulls
if [[ ! -e "$dir/phi/metadata.json" ]]; then
    . "$PROJECT/slurm/slurm_conda_bootstrap.sh"
    python -m normalize_tes.phi_sfs \
        --target "$dir/target" --matches "$dir/matches" \
        --vcf "$VCF" --ancestral-table "$ANCESTRAL" \
        -A SNP -B SNP \
        --reference-seed $(( 810000 + T )) \
        --min-null-replicates "$NULLS" --max-null-replicates "$NULLS" \
        --output "$dir/phi"
fi
echo "V6 test $T done"; date
