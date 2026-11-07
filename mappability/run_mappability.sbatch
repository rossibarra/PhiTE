#!/bin/bash -l
#SBATCH --job-name=nTE-bsmap
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=logs/bsmap-%j.out
#SBATCH --error=logs/bsmap-%j.err

# Bisulfite mappability of B73 v5 (whole genome, so uniqueness is genome-wide).
# GenMap's suffix-array build needs ~10 bytes per base for >2 GB inputs, i.e.
# ~43 GB for the 2 x 2.2 Gb converted reference, hence 96 GB.  The chr10-only
# pilot used 1.9 GB (index) and 1.2 GB (map); 256 is GenMap's default thread
# count, so -T is set explicitly to the CPUs allocated.
#
#   sbatch mappability/run_mappability.sbatch          # K=100, E=2
#   sbatch --export=ALL,K=150,E=0 mappability/run_mappability.sbatch
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
K=${K:-100}
E=${E:-2}
REF=${REF:-/quobyte/jrigrp/DATA/RESOURCES/NAM/fasta/Zm-B73-REFERENCE-NAM-5.0.fa.gz}
W=results/mappability/B73v5_k${K}_e${E}
OUT=results/mappability/B73v5_bisulfite_k${K}_e${E}
mkdir -p "$W" logs
module load conda >/dev/null 2>&1
conda activate "$PWD/.envs/mappability"

python mappability/bisulfite_convert.py "$REF" "$W/ct.fa"
rm -rf "$W/idx"
genmap index -F "$W/ct.fa" -I "$W/idx"
genmap map -I "$W/idx" -O "$W/map" -K "$K" -E "$E" -nc -fs -r -T "$SLURM_CPUS_PER_TASK"
python mappability/genmap_to_bed.py "$W/ct.fa" "$W/map.freq8" "$OUT.mappable.bed" \
    --K "$K" --summary "$OUT.summary.tsv"
echo "done: $OUT.mappable.bed"
