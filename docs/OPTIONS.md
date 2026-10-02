# Options and less common paths

The [README](../README.md) shows each step with only the flags production changes
from the defaults. This document lists every flag the steps use, the tuning
options, the launchers' variables, and the paths most runs never need. Each
command's `--help` gives its full contract. Variable names follow the README's
Configure block.

## 1. Interval store

```bash
python -m normalize_tes.build_snp_interval_store "$POSTERIOR_DIR"/*.tsz \
  --interval-store "$STORE" --chrom-offsets "$CHROM_OFFSETS" \
  --min-usable-fraction 0.1 --num-buckets 100 --bucket-memory-gb 2 \
  --scratch-dir "${TMPDIR:?}"
```

| flag | purpose | default |
|---|---|---|
| `trees` | posterior ARG draws, one tree sequence per draw | required |
| `--interval-store` | new store directory | required |
| `--chrom-offsets` | chromosome offsets, when compatible metadata is not embedded in every ARG | none |
| `--min-usable-fraction` | minimum fraction of draws with a usable age interval for a row to be eligible | 0.1 |
| `--num-buckets` | temporary row partitions; more buckets means less memory per bucket | 64 |
| `--bucket-memory-gb` | per-bucket sort-memory ceiling | 4 |
| `--scratch-dir` | temporary bucket location; use node-local scratch | required |

The store records a content digest used by downstream identity checks, and a
content identity for each source draw, so later steps can prove they were given
the same posterior wherever the files now live.

**Stores built before v0.7.0** record only a path per draw, so moving the draws
makes every later step reject them. Record their identities in place; this does
not change the store's content digest:

```bash
python -m normalize_tes.record_draw_identities --store "$STORE" "$POSTERIOR_DIR"/*.tsz
```

Pass `--dry-run` first. Files are matched to the recorded draws by name, one to
one; where the store already carries identities, the content must match and only
the path changes.

## 2. Ancestral-state table

| flag | purpose |
|---|---|
| `--store` | store whose rows and source draws the table must match |
| `--output` | new table directory |
| `trees` | the store's complete posterior draw set |
| `--draws START:STOP` | build one part of an array build |
| `--merge PARTS... --expect-draws N` | merge parts once every one has succeeded |

Each draw is authenticated against the store by content, not path. Tables use
schema `ancestral-state-counts-v2`, which records a SHA-256 digest of each array;
the merge and every reader verify them, so a corrupted or replaced array fails
instead of silently mispolarizing sites. A v1 table is refused and must be rebuilt.
The Farm launcher (`slurm/run_ancestral_table.sbatch`) wraps the array build and
merge; see the README.

## 3. VCF eligibility

| flag | purpose | default |
|---|---|---|
| `--vcf` | filtered biallelic analysis VCF | required |
| `--store` | interval store defining the row universe | required |
| `--ancestral-table` | ancestral-state table for site orientation | required |
| `--output` | new eligibility directory | required |
| `--min-callable` | minimum callable individuals per site | 20 |
| `--heterozygous` | `error` rejects heterozygous inbred calls; `missing` treats them as missing | `error` |

The artifact holds a shared genotype and callability mask, the subset of rows with
at least one usable ARG orientation, and each row's posterior orientation
probability. TE and SNP targets and the control universe all use that orientable
subset. Use the same artifact for A and B. If you build it with
`--heterozygous missing`, pass the same policy to Phi-SFS.

## 4. Candidate controls

| flag | purpose | default |
|---|---|---|
| `--store` | store whose rows are being selected | required |
| `--include-positions` | filtered SNP positions allowed as controls | required |
| `--exclude-positions` | all known TEs and every A position | none |
| `--vcf-eligibility` | restrict controls to callable, SNP-orientable rows | required in production |
| `--output` | new candidate `.npy`; a provenance report is written beside it | required |
| `--min-resolved-fraction` | minimum fraction of listed positions that must resolve to store rows; production uses 0.70 | 0.95 |

Candidate rows are store-specific; rebuild them whenever the store changes.
Excluding `A_POSITIONS` is essential for SNP-versus-SNP runs, and harmlessly
redundant when A is a subset of `ALL_TE_POSITIONS`.

## 5. Target and matcher

Target (`normalize_tes.te_age_target`):

| flag | purpose | default |
|---|---|---|
| `--store` | interval store supplying focal-site ages | required |
| `--te-positions` | focal positions, TE or SNP according to `-A` | required |
| `-A`, `--a-type` | focal type, `TE` or `SNP` | `TE` |
| `--vcf-eligibility` | keep only callable, ARG-orientable rows before fixing $M$ | required in production |
| `--output` | new target directory | required |
| `--scratch-dir` | node-local location for the site-by-age CDF matrix | required |
| `--bootstrap-replicates` | focal-site resamples used to calibrate the matching threshold | 10,000 |
| `--acceptance-quantile` | bootstrap-distance quantile used as the threshold | 0.50 |
| `--seed` | bootstrap seed; production uses 1002 | none |

The target streams its site-by-age CDF through `--scratch-dir`, so scratch space is
the limit. `--te-polarity-mask`, `--max-flipped-fraction` and
`normalize_tes.build_te_polarity_mask` remain only as diagnostics: they condition A
alone on ARG polarity, which no control set shares, and the matcher and Phi-SFS
refuse a target built with them.

Matcher (`normalize_tes.bootstrap_target_matcher`):

| flag | purpose | default |
|---|---|---|
| `--store`, `--target` | as above | required |
| `-A`, `--a-type` | must agree with the target | `TE` |
| `--candidate-rows` | control universe and its provenance sidecar | required |
| `--output` | new matched-control bundle | required |
| `--work-dir` | durable per-replicate state for `--resume` | required |
| `--resume` | continue an interrupted compatible run | off |
| `--replicates` | control sets to publish | 500 |
| `--restarts` | optimization restarts per set | 3 |
| `--disjoint-replicates` | never reuse a control SNP between sets; required in production | off |
| `--init-mode` | `median`, `mass` or `random` starting set per restart | `median` |
| `--seed` | global matching seed; production uses 1002 | 0 |

Disjoint mode first checks that the pool can hold `replicates` × $M$ controls and
that every age stratum has enough candidates, then removes each published control
from later pools. If the pool cannot support all sets it fails; it never falls back
to reuse. A resumed run must match the original's parameters, code and inputs
exactly, including whether the checkout had uncommitted edits.

## 6. Phi-SFS

| flag | purpose | default |
|---|---|---|
| `--target`, `--matches` | focal target and its matched bundle | required |
| `--vcf` | the analysis VCF | required |
| `--ancestral-table` | ancestral-state table; orients every TE and SNP | required |
| `-A`, `--a-type` | focal type | `TE` |
| `-B`, `--b-type` | control type; only `SNP` | `SNP` |
| `--reference-seed` | seed, with the target digest, for drawing $B_0$ | 1002 |
| `--reference-replicate` | use this replicate ID as $B_0$ instead of drawing it | none |
| `--min-null-replicates` | fail if fewer QC-passing nulls | 450 |
| `--max-null-replicates N` | fix $R=N$ from the seeded permutation; fails if fewer than $N+1$ pass | every passing set |
| `--reference-sensitivity N` | repeat calibration with the next $N$ sets as $B_0$ | 0 |
| `--asymmetric-polarity-null` | Bernoulli-$q$ option (see [METHODS.md](METHODS.md)) | off |
| `--polarity-imputation-seed` | seed for that option | 2001 |
| `--heterozygous` | must match the eligibility artifact | `error` |
| `--output` | new result directory | required |

`--max-null-replicates` cannot be combined with `--reference-replicate` or
`--reference-sensitivity`. Phi-SFS asserts that A and every accepted set keep
exactly the same $M$ sites; it never drops or downsamples sites.

`--reference-sensitivity N` writes `sensitivity_reference_ids.npy`,
`sensitivity_observed_phi_sfs.npy`, `sensitivity_z_scores.npy` and
`sensitivity_p_values.npy`, and fills the `reference_sensitivity_*` columns of
`summary.csv`. The primary result is unchanged by requesting it.
`comparisons.csv` reports each set's largest control-reuse count (all 1 in a valid
disjoint bundle).

With `--asymmetric-polarity-null`, Phi-SFS also writes the hard null spectra
(`b_bernoulli_q_*.npy`), and `metadata.json` records the null design, the seed, the
imputation algorithm, and the separate A, reference and null polarity rules.

## Farm launchers

Launchers in `slurm/` activate the conda environment themselves and are configured
by environment variables passed with `sbatch --export=ALL,...`. Each script's header
lists every variable with its default. The main ones:

| launcher | variables (default) |
|---|---|
| `run_bootstrap_matching.sbatch` | `STORE`, `TARGET`, `OUTPUT`, `CANDIDATE_ROWS`, `VCF_ELIGIBILITY` (required); `A_POSITIONS` (builds a missing target); `A_TYPE` (TE); `WORK_DIR` (`OUTPUT.work`); `REPLICATES` (500); `RESTARTS` (3); `SEED` (1002); `INIT_MODE` (median); `ACCEPTANCE_QUANTILE` (0.50); `MISSING_POSITION_POLICY` (error); `SCRATCH_HEADROOM_GB` (32) |
| `run_phi_sfs.sbatch` | `TARGET`, `MATCHES`, `VCF`, `ANCESTRAL`, `OUTPUT` (required); `A_TYPE` (TE); `B_TYPE` (SNP); `REFERENCE_SEED` (1002); `MIN_NULL_REPLICATES` (450); `MAX_NULL_REPLICATES` (unset); `REFERENCE_SENSITIVITY` (0); `ASYMMETRIC_POLARITY_NULL` (false); `POLARITY_IMPUTATION_SEED` (2001); `HETEROZYGOUS` (error) |
| `run_ancestral_table.sbatch` | `STORE`, `TREES`, `OUTPUT`, `PER_TASK` for array parts; `MERGE=1`, `PARTS`, `EXPECT_DRAWS` for the merge |

The matching launcher builds a missing target from `A_POSITIONS`. When `TARGET`
exists, it checks the recorded A type and eligibility artifact. It refuses a target
built with the TE polarity mask, exits if `TE_POLARITY_MASK` or
`MAX_FLIPPED_FRACTION` is set, rejects an unrestricted candidate universe, and
checks that the candidate provenance names the same eligibility artifact. It
always passes `--resume`, so an identical resubmission continues after preemption.
Measured runtimes and memory are in
[BOOTSTRAP_HPC_VALIDATION.md](BOOTSTRAP_HPC_VALIDATION.md).

## Many categories

The store, ancestral table, eligibility artifact and candidate universe are shared.
Each category needs its own target, bundle, durable work directory and seed. Write a
tab-separated manifest:

```text
label	positions	target	matches	work_dir	seed
all_te	/quobyte/project/te/all.pos.txt	/quobyte/project/targets/all_te	/quobyte/project/matches/all_te	/quobyte/project/work/all_te	1001
in_gene	/quobyte/project/te/in_gene.pos.txt	/quobyte/project/targets/in_gene	/quobyte/project/matches/in_gene	/quobyte/project/work/in_gene	1002
young	/quobyte/project/te/young.pos.txt	/quobyte/project/targets/young	/quobyte/project/matches/young	/quobyte/project/work/young	1003
```

- The first line is the header shown.
- Fields are separated by literal tabs; paths must not contain tabs, newlines or
  commas.
- Every category path must be unique; target and match paths must not exist on a
  first submission.
- A `work_dir` may be reused only to resume the identical matching run.
- Fix seeds before submission and keep them on resubmission.

Then submit one matching job per row:

```bash
PROJECT=/quobyte/project/PhiTE
STORE=/quobyte/project/data/age_interval_store
CANDIDATES=/quobyte/project/data/candidate_rows.npy
VCF_ELIGIBILITY=/quobyte/project/data/vcf_eligibility
MANIFEST=/quobyte/project/manifests/te_categories.tsv

while IFS=$'\t' read -r label positions target matches work seed; do
  [[ "$label" == label ]] && continue
  [[ -n "$label" ]] || continue
  match_job=$(sbatch --parsable --job-name="match-${label}" \
    --export=ALL,PROJECT="$PROJECT",STORE="$STORE",TARGET="$target",\
A_POSITIONS="$positions",A_TYPE=TE,OUTPUT="$matches",CANDIDATE_ROWS="$CANDIDATES",\
VCF_ELIGIBILITY="$VCF_ELIGIBILITY",WORK_DIR="$work",SEED="$seed" \
    slurm/run_bootstrap_matching.sbatch)
  printf '%s\tmatch=%s\n' "$label" "${match_job%%;*}"
done < "$MANIFEST"
```

Save the printed job IDs and confirm with `sacct` that every row completed. The
loop does not skip existing outputs; for a partial rerun, submit only the missing
categories, or resubmit an interrupted matcher with its original paths and seed.
Run Phi-SFS and the drift check for each category once its bundle is published.

## Output schemas

Outputs are published atomically and never overwritten. Phi-SFS writes schema
`phi-sfs-wasserstein-v3`; its generic `a_*` and `b_*` arrays serve both TE-SNP and
SNP-SNP runs. Version 3 polarizes TEs from the ARG posterior; version 2 gave TEs
hard biological polarity after a derived-support filter; version 1 used an earlier
symmetric-mixture null. No two versions may be combined. Matched control sets are
Monte Carlo null replicates, not independent biological samples.
