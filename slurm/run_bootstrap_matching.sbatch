#!/usr/bin/env bash
#SBATCH --job-name=nTE-bootstrap-match
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=96G
#SBATCH --time=1-12:00:00
#SBATCH --output=bootstrap-match-%j.out
#SBATCH --error=bootstrap-match-%j.err

# Scheduler output goes to the submission directory, not logs/. SLURM opens
# these files before any line of this script runs, so a directory the script
# creates later cannot exist yet: on a fresh clone the job would fail before
# reaching mkdir. Move the logs afterwards if you want them filed away.

# The canonical production launcher for bootstrap-target control matching.
# It requires mutually disjoint controls across all 1201 sets: one reference,
# 1000 nulls, and 200 spares that replace sets failing matching QC. The matcher
# fails before creating replicate state when the candidate universe cannot hold
# REPLICATES x target-size unique SNPs, or when any target age stratum holds
# fewer than REPLICATES x its quota; it never silently falls back to reuse.
#
# Submit with sbatch, never srun: an srun started from an interactive shell dies
# with that shell, which is how the first production attempt was lost at 48
# minutes. Completed replicate bundles under WORK_DIR survive job death, so the
# identical resubmission continues from where it stopped (--resume is always
# passed below).
#
# Environment variables, all overridable at submit time:
#
#   PROJECT          repository checkout (default: the directory holding this file)
#   STORE            canonical interval store on durable storage        [required]
#   TARGET           target directory built by normalize_tes.te_age_target         [required]
#   OUTPUT           bundle directory to publish; must not exist        [required]
#   CANDIDATE_ROWS   .npy from normalize_tes.build_candidate_rows       [required]
#   CANDIDATE_REPORT provenance JSON (default: CANDIDATE_ROWS.json)
#   WORK_DIR         durable resume directory       (default: OUTPUT.work)
#   REPLICATES       control sets to publish                   (default: 1201)
#   RESTARTS         stratified restarts per replicate            (default: 3)
#   SEED             global seed                                (default: 1002)
#   A_TYPE           TE|SNP for the focal target                         (TE)
#   A_POSITIONS      build TARGET from this position file if TARGET is absent
#   VCF_ELIGIBILITY  shared authenticated callability/orientation artifact [required]
#   ACCEPTANCE_QUANTILE        bootstrap quantile for the target threshold (0.50)
#   MISSING_POSITION_POLICY    error|drop for target construction     (error)
#   SCRATCH_HEADROOM_GB        scratch reserved for the target CDF      (32)
#   TE_POLARITY_MASK           polarity mask dir, when building a target (unset)
#   MAX_FLIPPED_FRACTION       discard TEs flipped above this fraction   (0.5)
#
# Example:
#
#   sbatch --export=ALL,STORE=/quobyte/project/data/snp_interval_store,\
#   TARGET=results/targets/in_gene,OUTPUT=results/bootstrap_matches/in_gene,\
#   CANDIDATE_ROWS=results/candidate-rows.npy,VCF_ELIGIBILITY=results/vcf-eligibility,\
#   WORK_DIR=results/work-in-gene,SEED=1002 \
#     slurm/run_bootstrap_matching.sbatch
#
# Measured on the 75-draw production store (18.2 GB) at this allocation, for a
# 4,067-site target with 100 replicates x 3 restarts: 2 h 12 m wall clock and
# 37.7 GiB peak RSS (SLURM job 37821602). Scaled linearly from that single
# measurement (unverified at this size), the 1201-set default needs about 26 h
# for the same target; the 36 h request leaves room for variation. Target
# construction is the other memory peak: 16.0 GB for a
# 35,512-site target. Scale --time roughly linearly in target size and replicate
# count, and keep --mem at 96G unless a measured figure says otherwise.

set -euo pipefail

# sbatch copies the job script into /var/spool/slurmd, so BASH_SOURCE does not
# point at the checkout. SLURM_SUBMIT_DIR does; the BASH_SOURCE form is the
# fallback for running this script directly.
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}}"
: "${STORE:?set STORE to the completed interval-store directory}"
: "${TARGET:?set TARGET to the target directory}"
: "${OUTPUT:?set OUTPUT to the bundle directory to publish}"
: "${CANDIDATE_ROWS:?set CANDIDATE_ROWS to the VCF-eligible SNP candidate .npy}"
: "${VCF_ELIGIBILITY:?set VCF_ELIGIBILITY to the shared eligibility artifact}"
: "${TMPDIR:?Farm must provide per-job TMPDIR under /local/scratch}"
WORK_DIR="${WORK_DIR:-${OUTPUT}.work}"
REPLICATES="${REPLICATES:-1201}"
RESTARTS="${RESTARTS:-3}"
SEED="${SEED:-1002}"
A_TYPE="${A_TYPE:-TE}"
A_POSITIONS="${A_POSITIONS:-${TE_POSITIONS:-}}"
ACCEPTANCE_QUANTILE="${ACCEPTANCE_QUANTILE:-0.50}"
MISSING_POSITION_POLICY="${MISSING_POSITION_POLICY:-error}"
TE_POLARITY_MASK="${TE_POLARITY_MASK:-}"
MAX_FLIPPED_FRACTION="${MAX_FLIPPED_FRACTION:-0.5}"
CANDIDATE_REPORT="${CANDIDATE_REPORT:-${CANDIDATE_ROWS}.json}"

if [[ "$A_TYPE" != "TE" && "$A_TYPE" != "SNP" ]]; then
    echo "A_TYPE must be TE or SNP" >&2
    exit 5
fi
if [[ "$CANDIDATE_ROWS" == "all" ]]; then
    echo "CANDIDATE_ROWS=all is incompatible with shared SNP eligibility; build the candidate .npy first" >&2
    exit 5
fi
if [[ "$A_TYPE" == "TE" && -z "$TE_POLARITY_MASK" ]]; then
    echo "A_TYPE=TE requires TE_POLARITY_MASK for the inclusive >=50% derived rule" >&2
    exit 5
fi
if [[ "$A_TYPE" == "SNP" && -n "$TE_POLARITY_MASK" ]]; then
    echo "A_TYPE=SNP must not use TE_POLARITY_MASK" >&2
    exit 5
fi
if [[ "$A_TYPE" == "SNP" ]]; then
    MAX_FLIPPED_FRACTION=""
fi

cd "$PROJECT"
mkdir -p "$TMPDIR/interval_store"

# Batch jobs get a non-login shell, where `module` does not exist.
. "$PROJECT/slurm/slurm_conda_bootstrap.sh"

threads="${SLURM_CPUS_PER_TASK:-1}"
export OPENBLAS_NUM_THREADS="$threads"
export OMP_NUM_THREADS="$threads"
export MKL_NUM_THREADS="$threads"

# Staging is mandatory, not an optimization. Against the store on Quobyte the
# run is I/O-bound and projects to 4-21 h of pure store reads for one target;
# staged to node-local scratch the same run measured 97% CPU with 21 major page
# faults. Every job that touches the store stages it first.
# Scratch must hold the staged store AND, on the unmasked path, the temporary
# TE-by-age CDF that target construction writes beside it. Reserving only the
# store plus 20% lets the preflight pass and the job still run out mid-build.
# SCRATCH_HEADROOM_GB is that second term; raise it for large TE categories.
SCRATCH_HEADROOM_GB="${SCRATCH_HEADROOM_GB:-32}"
required_kib=$(du -sk "$STORE" | awk '{print $1}')
headroom_kib=$(( SCRATCH_HEADROOM_GB * 1024 * 1024 ))
free_kib=$(df -Pk "$TMPDIR" | awk 'NR==2 {print $4}')
if (( free_kib < required_kib * 6 / 5 + headroom_kib )); then
    echo "insufficient node-local scratch: need the store size plus 20% " \
         "(${required_kib} KiB) plus ${SCRATCH_HEADROOM_GB}G for the target CDF" >&2
    exit 2
fi

# VCF eligibility is compared by content, not path: a path names a location,
# not what is at it. The identity (VCF digest, genotype policy, callability
# threshold, store and array digests) is computed once from the requested
# artifact's metadata.json and must equal the one recorded by the candidate
# report here and by the target below. The matcher itself also requires the
# candidate report and target to agree, so this launcher check only adds the
# binding to the artifact named by VCF_ELIGIBILITY.
ELIGIBILITY_IDENTITY="$(python - "$VCF_ELIGIBILITY" <<'PYEOF'
import json, sys
from pathlib import Path

metadata_path = Path(sys.argv[1]) / "metadata.json"
if not metadata_path.is_file():
    raise SystemExit(f"missing VCF eligibility metadata: {metadata_path}")
metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
if metadata.get("complete") is not True:
    raise SystemExit(f"{metadata_path} is not a complete eligibility artifact")
from normalize_tes.vcf_eligibility import eligibility_identity
print(json.dumps(eligibility_identity(metadata), sort_keys=True))
PYEOF
)"
export ELIGIBILITY_IDENTITY

python - "$CANDIDATE_REPORT" "$VCF_ELIGIBILITY" <<'PYEOF'
import json, os, sys
from pathlib import Path

report_path, eligibility = map(Path, sys.argv[1:])
wanted = json.loads(os.environ["ELIGIBILITY_IDENTITY"])
if not report_path.is_file():
    raise SystemExit(f"missing candidate provenance report: {report_path}")
report = json.loads(report_path.read_text(encoding="utf-8"))
used = report.get("vcf_eligibility_identity")
if used != wanted:
    differing = sorted(
        key for key in set(wanted) | set(used or {})
        if (used or {}).get(key) != wanted.get(key)
    )
    raise SystemExit(
        f"{report_path} was built from a different VCF eligibility artifact than "
        f"{eligibility} (identity fields differing: {', '.join(differing)})"
    )
PYEOF
# A masked target does not use the scratch-backed CDF writer: te_age_target
# builds the whole block in memory. For those runs --mem, not scratch, is the
# constraint, and the job prints its projected peak before building.
if [[ "$A_TYPE" == "TE" ]]; then
    echo "note: masked target construction builds CDFs in memory; --mem is the" \
         "binding constraint, not scratch"
fi
echo "staging $STORE to $TMPDIR/interval_store"; date
rsync -a --delete "$STORE/" "$TMPDIR/interval_store/"
date

# A draw that mis-polarized a TE placed its mutation on a different branch, so
# the age it recorded belongs to a different event. With a mask, each TE's age
# CDF is built from its agreeing draws only.
polarity_option=()
if [[ "$A_TYPE" == "TE" ]]; then
    polarity_option+=(--te-polarity-mask "$TE_POLARITY_MASK")
    if [[ -n "$MAX_FLIPPED_FRACTION" ]]; then
        polarity_option+=(--max-flipped-fraction "$MAX_FLIPPED_FRACTION")
    fi
elif [[ -n "$MAX_FLIPPED_FRACTION" ]]; then
    echo "MAX_FLIPPED_FRACTION requires TE_POLARITY_MASK" >&2
    exit 5
fi

if [[ ! -d "$TARGET" ]]; then
    if [[ -z "$A_POSITIONS" ]]; then
        echo "TARGET $TARGET does not exist and A_POSITIONS is unset" >&2
        exit 4
    fi
    echo "=== building target $TARGET ==="; date
    python -m normalize_tes.te_age_target \
        --store "$TMPDIR/interval_store" \
        --te-positions "$A_POSITIONS" \
        -A "$A_TYPE" \
        --vcf-eligibility "$VCF_ELIGIBILITY" \
        --output "$TARGET" \
        --scratch-dir "$TMPDIR" \
        --bootstrap-replicates 10000 \
        --acceptance-quantile "$ACCEPTANCE_QUANTILE" \
        --bin-width 1000 \
        --missing-position-policy "$MISSING_POSITION_POLICY" \
        "${polarity_option[@]}" \
        --seed "$SEED"
    date
else
    # The mask builder needs an existing target to supply te_row_indices.npy, so
    # the real workflow is
    #     preliminary unmasked target -> mask -> final masked target -> matching
    # and the two targets must be different directories. An operator who passes
    # the preliminary target here together with TE_POLARITY_MASK would otherwise
    # get matching against the unmasked target: the block above is skipped, the
    # mask is silently unused, and nothing in the output says so.
    echo "=== target $TARGET already exists; verifying it matches the request ==="
    python - "$TARGET" "$A_TYPE" "$VCF_ELIGIBILITY" "$TE_POLARITY_MASK" "$MAX_FLIPPED_FRACTION" <<'PYEOF'
import json, os, sys
from pathlib import Path

target, a_type, eligibility, requested_mask, requested_fraction = sys.argv[1:]
meta = json.loads((Path(target) / "metadata.json").read_text(encoding="utf-8"))
recorded = meta.get("te_polarity") or None

if meta.get("a_type") != a_type:
    raise SystemExit(f"{target} records a_type={meta.get('a_type')!r}, requested {a_type!r}")
eligibility_record = meta.get("vcf_eligibility") or {}
used_identity = eligibility_record.get("identity")
wanted_identity = json.loads(os.environ["ELIGIBILITY_IDENTITY"])
wanted_eligibility = str(Path(eligibility).resolve())
if used_identity != wanted_identity:
    differing = sorted(
        key for key in set(wanted_identity) | set(used_identity or {})
        if (used_identity or {}).get(key) != wanted_identity.get(key)
    )
    raise SystemExit(
        f"{target} was built from a different VCF eligibility artifact than "
        f"VCF_ELIGIBILITY={wanted_eligibility} (identity fields differing: "
        f"{', '.join(differing)})."
    )

if a_type == "SNP":
    if recorded is not None:
        raise SystemExit(f"SNP target {target} unexpectedly records a TE polarity mask")
    print(f"target matches A_TYPE=SNP and eligibility {wanted_eligibility}")
    sys.exit(0)

if recorded is None:
    raise SystemExit(
        f"TE_POLARITY_MASK was set, but {target} already exists and was built "
        "WITHOUT a polarity mask. Matching would silently use the unmasked "
        "target. Point TARGET at a new directory for the masked target; the "
        "unmasked one is the preliminary target the mask was built against."
    )

used = str(Path(recorded["mask"]).resolve())
want = str(Path(requested_mask).resolve())
if used != want:
    raise SystemExit(
        f"{target} was built with polarity mask {used}, but TE_POLARITY_MASK "
        f"requests {want}."
    )

recorded_fraction = recorded.get("max_flipped_fraction")
wanted_fraction = float(requested_fraction) if requested_fraction else None
if recorded_fraction != wanted_fraction:
    raise SystemExit(
        f"{target} was built with max_flipped_fraction={recorded_fraction}, "
        f"but MAX_FLIPPED_FRACTION requests {wanted_fraction}."
    )
print(f"target matches the requested mask: {recorded['sites_kept']:,} TEs kept, "
      f"max_flipped_fraction={recorded_fraction}")
PYEOF
fi

candidate_option=(--candidate-rows "$CANDIDATE_ROWS")

echo "=== matching: $REPLICATES replicates x $RESTARTS restarts ==="; date
python -m normalize_tes.bootstrap_target_matcher \
    --store "$TMPDIR/interval_store" \
    --target "$TARGET" \
    -A "$A_TYPE" \
    "${candidate_option[@]}" \
    --output "$OUTPUT" \
    --work-dir "$WORK_DIR" \
    --resume \
    --replicates "$REPLICATES" \
    --restarts "$RESTARTS" \
    --disjoint-replicates \
    --seed "$SEED"
date
