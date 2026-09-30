# PhiTE v0.9.0-rc1

PhiTE builds neutral SNP control sets matched to the posterior ages of a focal
variant category, then compares their unfolded site-frequency spectra. Dataset A may
be either TEs or SNPs; dataset B is currently SNPs. The production workflow therefore
supports the primary TE-versus-SNP analysis and SNP-versus-SNP negative controls.

The alternative pipeline for estimating the PDF of derived allele ages for a
sample is found [here](derived_distribution_readme.md).

This README is the operator guide. Method definitions and the evidence behind the
production settings are linked under [Methods and validation](#methods-and-validation).

## Citation

If you use PhiTE or the Phi-SFS method, please cite:

> Liu, B., Munasinghe, M., Fairbanks, R. A., Hirsch, C. N., and Ross-Ibarra, J. (2025).
> Genome-wide selection on transposable elements in maize. bioRxiv 2025.09.16.676665.
> <https://doi.org/10.1101/2025.09.16.676665>

## Repository layout

- `normalize_tes/` contains the production package and supported command modules.
- `tools/` contains development benchmarks, diagnostics, probes, and simulations.
- `slurm/` contains scheduler launchers and their shared conda bootstrap helper.
- `tests/` contains the complete test suite.
- `docs/` contains methods, validation records, plans, reviews, and the changelog.

Run production commands from the repository root with
`python -m normalize_tes.COMMAND`; examples below use that form throughout.

## Before you start

Create the environment and run the tests on a Linux compute node:

```bash
conda env create -f environment.yml
conda activate normalizeTE
python -m pytest -q tests
```

Use an immutable release tag or commit for production:

```bash
git fetch --tags
git checkout COMMIT_HASH
```

The workflow expects:

- posterior ARG draws; the ancestry stage requires tszip archives;
- SNP and TE position files with two whitespace-separated columns: chromosome and
  1-based VCF position;
- chromosome labels matching the ARG metadata or a compatible chromosome-offset
  file supplied when the store is built;
- a filtered, genome-wide, biallelic VCF in which every TE is an ACGT-coded site
  with the same alleles as in the ARGs (this dataset uses `A` for absence and `G`
  for presence), so the ancestral table can orient it like a SNP;
- node-local `$TMPDIR` for interval-store and target scratch data;
- durable storage for published outputs and matcher resume state.

Blank lines and `#` comments are allowed in position files. Run each script with
`--help` for its complete input contract and advanced options.

## Configure one category

Set these paths once and use them throughout the commands below:

```bash
POSTERIOR_DIR=/path/project-data/posterior
CHROM_OFFSETS=/path/project-data/chrom_offsets.txt
SNP_POSITIONS=/path/project-data/snp/all_snp.pos.txt
ALL_TE_POSITIONS=/path/project-data/te/all_te.pos.txt
TE_POSITIONS=/path/project-data/te/in_gene.pos.txt
A_POSITIONS="$TE_POSITIONS"
A_TYPE=TE
B_TYPE=SNP
VCF=/path/variants.vcf.gz

STORE=results/age_interval_store
CANDIDATES=results/candidate_rows.npy
TARGET=results/targets/in_gene
MATCHES=results/bootstrap_matches/in_gene
WORK_DIR=results/work/in_gene
ANCESTRAL=results/ancestral_states
ELIGIBILITY=results/vcf_eligibility
PHI=results/phi_sfs/in_gene

mkdir -p results results/targets \
  results/bootstrap_matches results/work results/phi_sfs
```

Every published output path must be new. The tools refuse to overwrite existing
artifacts.

## Run the pipeline

The commands below show the production TE-versus-SNP path and note the changes for a
SNP-versus-SNP negative control. Run heavy commands inside a scheduled compute
allocation, not on a login/head node.

### 1. Build the interval store

Build one reusable posterior-age store for all categories:

```bash
python -m normalize_tes.build_snp_interval_store \
  "$POSTERIOR_DIR"/*.tsz \
  --interval-store "$STORE" \
  --chrom-offsets "$CHROM_OFFSETS" \
  --min-usable-fraction 0.1 \
  --num-buckets 100 \
  --bucket-memory-gb 2 \
  --scratch-dir "${TMPDIR:?TMPDIR is not set}"
```

| flag | purpose |
|---|---|
| `trees` | posterior ARG draws, one tree sequence per draw |
| `--interval-store` | new store directory to publish |
| `--chrom-offsets` | chromosome offsets when compatible metadata is not embedded in every ARG |
| `--min-usable-fraction` | minimum fraction of draws with a usable age interval for an eligible row |
| `--num-buckets` | number of temporary row partitions; increase this to reduce per-bucket memory |
| `--bucket-memory-gb` | per-bucket sort-memory ceiling |
| `--scratch-dir` | temporary bucket location; use node-local scratch |

Omit `--chrom-offsets` only when every draw contains compatible chromosome metadata.
The completed store records a content digest used by downstream identity checks, and a
content identity for each source draw so that later steps can prove they were handed the
same posterior, wherever those files now live.

A store built before v0.7.0 records only a path per draw, so moving the draws makes every
later step reject them. Record their identities in place -- this does not change the
store's content digest, and leaves artifacts already stamped with it valid:

```bash
python -m normalize_tes.record_draw_identities \
  --store "$STORE" \
  "$POSTERIOR_DIR"/*.tsz
```

Pass `--dry-run` first to see what it would change. The files are matched to the store's
recorded draws by name, one to one, because an unauthenticated store offers nothing
better; where the store already carries identities, the content must match and only the
path moves.

### 2. Build the ancestral-state table

Build one store-aligned ancestral-state table. It supplies the posterior
orientation probability of every site, TE and SNP, in focal set A and in the B
controls:

```bash
python -m normalize_tes.build_ancestral_states \
  --store "$STORE" \
  --output "$ANCESTRAL" \
  "$POSTERIOR_DIR"/*.tsz
```

| flag | purpose |
|---|---|
| `--store` | store whose rows and source draws the table must match |
| `--output` | new ancestral-state table directory |
| `trees` | the store's complete posterior draw set |

Each draw is authenticated against the store by content rather than path. For an
array build, use `--draws START:STOP` for each part and merge the parts with
`--merge ... --expect-draws N`; the launcher example below shows this pattern.
Tables use schema `ancestral-state-counts-v2`, which records a SHA-256 digest of
each array. The merge and every reader verify the digests, so a corrupted or
replaced array fails instead of silently mispolarizing sites. A v1 table is
refused and must be rebuilt.

### 3. Build the shared VCF eligibility artifact

Scan the analysis VCF once before fixing the focal-site count or matching controls:

```bash
python -m normalize_tes.vcf_eligibility \
  --vcf "$VCF" \
  --store "$STORE" \
  --ancestral-table "$ANCESTRAL" \
  --output "$ELIGIBILITY" \
  --min-callable 20 \
  --heterozygous error
```

| flag | purpose |
|---|---|
| `--vcf` | filtered biallelic analysis VCF |
| `--store` | interval store defining the row universe |
| `--ancestral-table` | authenticated posterior ancestral-state counts for site orientation |
| `--output` | new eligibility directory |
| `--min-callable` | minimum callable individuals per site; production uses 20 |
| `--heterozygous` | reject heterozygous inbred calls, or treat them as missing |

The artifact contains a shared genotype/callability mask and the subset of rows
with at least one usable ARG orientation, plus each row's full posterior orientation
probability. TE and SNP targets and the control universe all use that orientable
subset, because Phi-SFS polarizes every site from the ARG. Use this same artifact
for A and B. Here lowercase $m=20$ is the number
of individuals used by the SFS projection; uppercase $M$ below is the number of sites
in each compared set.

### 4. Build the candidate control universe

Restrict controls to the filtered SNP list and exclude all known TE positions:

```bash
python -m normalize_tes.build_candidate_rows \
  --store "$STORE" \
  --include-positions "$SNP_POSITIONS" \
  --exclude-positions "$ALL_TE_POSITIONS" "$A_POSITIONS" \
  --vcf-eligibility "$ELIGIBILITY" \
  --output "$CANDIDATES" \
  --min-resolved-fraction 0.70
```

| flag | purpose |
|---|---|
| `--store` | store whose canonical rows are being selected |
| `--include-positions` | filtered SNP positions allowed in the control universe |
| `--exclude-positions` | all known TEs and every A position to remove from B |
| `--vcf-eligibility` | restrict B to callable, SNP-orientable rows before matching |
| `--output` | new candidate-row `.npy`; a provenance report is written beside it |
| `--min-resolved-fraction` | minimum fraction of requested positions that must resolve to store rows |

Candidate rows are store-specific. Rebuild this artifact whenever the store changes.
Excluding `A_POSITIONS` is essential for SNP-versus-SNP runs; it is harmlessly
redundant when A is a subset of `ALL_TE_POSITIONS`.
The justification for the production resolution threshold belongs in the validation
record, not in this how-to.

### 5. Build the target and match controls

Build the focal target from all posterior draws and construct the matched control
sets:

```bash
python -m normalize_tes.te_age_target \
  --store "$STORE" \
  --te-positions "$A_POSITIONS" \
  --output "$TARGET" \
  --scratch-dir "${TMPDIR:?TMPDIR is not set}" \
  --a-type "$A_TYPE" \
  --vcf-eligibility "$ELIGIBILITY" \
  --bootstrap-replicates 10000 \
  --acceptance-quantile 0.50 \
  --seed 1002

python -m normalize_tes.bootstrap_target_matcher \
  --store "$STORE" \
  --target "$TARGET" \
  -A "$A_TYPE" \
  --candidate-rows "$CANDIDATES" \
  --output "$MATCHES" \
  --work-dir "$WORK_DIR" \
  --resume \
  --replicates 1001 \
  --restarts 3 \
  --disjoint-replicates \
  --seed 1002
```

Target flags:

| flag | purpose |
|---|---|
| `--store` | interval store supplying focal-site ages |
| `--te-positions` | focal positions to resolve and summarize (TE or SNP, per `--a-type`) |
| `--output` | new target directory |
| `--scratch-dir` | node-local location for the temporary site-by-age CDF matrix |
| `--a-type` | focal dataset type, `TE` or `SNP` |
| `--vcf-eligibility` | keep only callable, ARG-orientable rows before fixing $M$ |
| `--bootstrap-replicates` | focal-site resamples used to calibrate the matching threshold |
| `--acceptance-quantile` | bootstrap-distance quantile used as that threshold |
| `--seed` | bootstrap random seed |

Sizing note: the target streams its site-by-age CDF through `--scratch-dir`, so
scratch is the constraint. Measured resource figures are in
[BOOTSTRAP_HPC_VALIDATION.md](docs/BOOTSTRAP_HPC_VALIDATION.md).

Matcher flags:

| flag | purpose |
|---|---|
| `--store` | store supplying candidate SNP ages |
| `--target` | focal target built above |
| `-A`, `--a-type` | focal type, `TE` or `SNP`; must agree with the target metadata |
| `--candidate-rows` | TE-excluded control universe and its provenance sidecar |
| `--output` | new matched-control bundle |
| `--work-dir` | durable per-replicate state used by `--resume` |
| `--resume` | continue an interrupted compatible run |
| `--replicates` | total matched control sets; production uses 1001, from which Phi-SFS draws $B_0$ and takes every other QC-passing set as a null |
| `--restarts` | optimization restarts per control set |
| `--disjoint-replicates` | prevent reuse of a control SNP between published sets |
| `--seed` | matching random seed |

Keep `WORK_DIR` on durable storage and repeat the identical command after
preemption. Disjoint mode preflights the necessary pool size: at least `replicates`
times $M$ eligible candidates must exist. It then removes every published control
from later candidate pools. If the pool cannot support all sets, the run fails; it
never falls back to reuse.

A TE target is built exactly as a SNP target is: ages from all usable posterior
draws, no filter on ARG polarity, and the same orientable eligibility rows. The
matcher and Phi-SFS refuse a target built with `--te-polarity-mask` or
`--max-flipped-fraction`, for either type. Those options, and
`normalize_tes.build_te_polarity_mask`, remain only as diagnostics: the
at-least-50%-derived filter and agreeing-draw ages conditioned A alone on ARG
polarity, which no SNP control set shared (see
[Polarity](#polarity-a-single-posterior-rule)). In both cases the shared
eligibility filter must have been applied before $M$ and the focal-age CDF were
fixed. The age matcher does not use allele frequency and remains SFS-blind.

For a SNP focal set, set `A_TYPE=SNP`. Despite the historical flag and array
names, `--te-positions` then defines the SNP focal set A, whose rows must already
have been excluded from `CANDIDATES`.

### 6. Calculate Phi-SFS

Calculate the unfolded SFS comparison for focal set A and every matched B set:

```bash
python -m normalize_tes.phi_sfs \
  --target "$TARGET" \
  --matches "$MATCHES" \
  --vcf "$VCF" \
  --ancestral-table "$ANCESTRAL" \
  -A "$A_TYPE" \
  -B "$B_TYPE" \
  --reference-seed 1002 \
  --min-null-replicates 900 \
  --output "$PHI"
```

| flag | purpose |
|---|---|
| `--target` | final focal A target |
| `--matches` | matched-control bundle from step 5 |
| `--vcf` | filtered genome-wide biallelic VCF covering all requested sites |
| `--ancestral-table` | store-aligned posterior ancestral-state table; orients every TE and SNP |
| `-A`, `--a-type` | focal type: `TE` (default) or `SNP` |
| `-B`, `--b-type` | control type; currently `SNP` only |
| `--reference-seed` | seed, combined with the target digest, for drawing $B_0$ uniformly from the QC-passing sets; default 1002 |
| `--reference-replicate` | optional explicit $B_0$ replicate ID, overriding the draw |
| `--min-null-replicates` | floor on $R$: every QC-passing non-reference set is a null, and the run fails if fewer than this pass; default 900 |
| `--reference-sensitivity` | optionally repeat calibration with $N$ alternative references, the next $N$ sets of the same seeded permutation; default 0 |
| `--max-null-replicates` | optionally fix $R=N$: $B_0$ is the first QC-passing set of the seeded permutation and the nulls are the next $N$; fails if fewer than $N+1$ pass; cannot be combined with `--reference-replicate` or `--reference-sensitivity`; default unset, meaning every QC-passing set |
| `--asymmetric-polarity-null`, `--no-asymmetric-polarity-null` | hard-orient A and every null by one Bernoulli-$q$ draw per site, keeping $B_0$ a mixture; off by default, when A, $B_0$ and every null are posterior mixtures |
| `--polarity-imputation-seed` | seed for the coordinate-keyed Bernoulli-$q$ orientations; used only with `--asymmetric-polarity-null`; default 2001 |
| `--output` | new Phi-SFS result directory |

The default rejects heterozygous calls. Use `--heterozygous missing` only when the
eligibility artifact was built with the same policy. Eligibility is fixed upstream;
the calculation asserts that A and every accepted B set retain exactly the same $M$
sites rather than silently dropping or downsampling sites.

`--reference-sensitivity` reruns the primary selection rule once per alternative
reference (each uses every other QC-passing set as its nulls, so the primary $B_0$
becomes a null in that rerun) and publishes
`sensitivity_reference_ids.npy`, `sensitivity_observed_phi_sfs.npy`,
`sensitivity_z_scores.npy`, and `sensitivity_p_values.npy`, plus
`reference_sensitivity_n`/`_z_min`/`_z_max`/`_p_min`/`_p_max` columns in
`summary.csv` (empty when $N=0$). Primary outputs are unaffected by requesting it.
`comparisons.csv` also reports `left_max_control_reuse` and
`right_max_control_reuse`, each set's largest global control-reuse count from the
match bundle's `reuse_row_indices.npy`/`reuse_counts.npy` (empty for A; all 1 in a
valid disjoint bundle).

With `--asymmetric-polarity-null`, PhiTE additionally writes
`b_bernoulli_q_raw_sfs.npy`, `b_bernoulli_q_normalized_sfs.npy`, and
`b_bernoulli_q_cdf.npy`. `metadata.json` records `null_polarity_design`,
`asymmetric_polarity_null`, `polarity_imputation_seed`,
`polarity_imputation_algorithm`, and the separate A, B-reference, and null-left
polarity rules, and `te_sites_polarized`, `te_usable_arg_draws`, and
`te_unusable_arg_draws` for the focal TEs.

#### Wasserstein definition

The revised statistic compares the normalized cumulative unfolded SFS of a focal
set $A$ with that of an age-matched neutral SNP set, $B_0$. Let $F_A(x)$
and $F_{B_0}(x)$ be their CDFs on the derived-allele-frequency (DAF) axis. Define

$$
\Phi_{\mathrm{SFS}}(A,B_0)
= W_1(A,B_0)
= \int_0^1 \left|F_A(x)-F_{B_0}(x)\right|\,dx.
$$

For equally spaced projected DAF bins $x_j=j/m$, where $m=20$ is the projection
sample size, calculate this exactly as

$$
\Phi_{\mathrm{SFS}}(A,B_0)
= \frac{1}{m}\sum_{j=1}^{m-1}
  \left|F_A(x_j)-F_{B_0}(x_j)\right|.
$$

Thus, Φ-SFS is the shaded area between the two CDFs. It is zero only when the
spectra are identical and grows as probability mass must move farther along the DAF
axis. It is unsigned: the CDFs and bin-level residuals show whether the focal set has
an excess of rare or high-frequency derived alleles.

![Schematic definition of Phi-SFS as the area between focal and neutral SFS cumulative distribution functions](figures/phi_sfs_definition_schematic.png)

#### Polarity: a single posterior rule

Every site, TE or SNP, focal or control, is polarized the same way. Let $p$ be the
observed ALT frequency and let

$$
q=P(\mathrm{ALT\ is\ derived}\mid\mathrm{usable\ ARG\ draws})
$$

be the fraction of usable ARG draws in which ALT is derived. The site contributes
the posterior mixture of its two orientations,

$$
q\,h(k,n)+(1-q)\,h(n-k,n).
$$

All $q\in[0,1]$ are retained, and there is no polarity filter. ARG draws that cannot
orient either observed allele are reported as unusable, not counted toward either
direction. A TE is an ACGT-coded site in the ARGs like any other, so the same
ancestral table gives its $q$. Where the ARG calls absence derived, the TE
contributes at the absence frequency, just as a SNP whose REF is derived contributes
at $1-p$.

TE presence is known biologically to be the derived state for a true single
insertion, and earlier designs used that: a TE contributed a hard insertion-derived
spectrum, TEs with less than 50% derived support were discarded, and TE ages came
from agreeing draws only. That made the observed comparison differ in construction
from every null comparison. A had true hard polarity, while the nulls carried the
ARG's posterior polarity, so any miscalibration of $q$ entered
$\Phi_{\mathrm{obs}}$ and no $\Phi_i^0$
([review 12](docs/CODE_REVIEW_ROUND12.md), findings 1 and 3). In the dnAging
simulations, where every site is neutral, that construction rejected far above the
nominal rate once the ARG was inferred without ancestral information. The table
gives rejection fractions at $\alpha=0.05$, with 150 tests per cell, from
`results/sim_dnaging_unpolarised/arms_v1/summary.csv`:

| focal construction (nulls Bernoulli-$q$, $B_0$ mixture) | $M=250$ | $M=1000$ | $M=4000$ |
|---|---|---|---|
| true polarity everywhere (oracle) | 0.09 | 0.05 | 0.09 |
| A true-hard, no filter, true ages | 0.24 | 0.44 | 0.83 |
| A true-hard, filtered, agreeing-draw ages, inferred ages (former production) | 0.27 | 0.75 | 0.97 |
| A Bernoulli-$q$, like the nulls | 0.09 | 0.10 | 0.07 |

The reference-haplotype replicates (`results/sim_dnaging_refhap/arms_v1/`) show the
same pattern. Treating A like the nulls removes the excess, so the production
design now treats TEs as SNPs are treated. The cost is that biological knowledge of
TE polarity is not used. The test asks whether TE sites differ from age-matched SNPs
as both are seen through the ARG, not what the TE's true spectrum is, and signal is
probably attenuated where the ARG mis-polarizes TEs.

**Default: posterior mixture throughout.** A, $B_0$ and every null set $B_i$ are
posterior mixtures:

$$
\Phi_{\mathrm{obs}}=\Phi_{\mathrm{SFS}}(A_{\mathrm{mix}},B_{0,\mathrm{mix}}),
\qquad
\Phi_i^0=\Phi_{\mathrm{SFS}}(B_{i,\mathrm{mix}},B_{0,\mathrm{mix}}).
$$

Because $E[h(\mathrm{DAF})]=q\,h(k,n)+(1-q)\,h(n-k,n)$, a Bernoulli-$q$ hard spectrum
is the mixture plus imputation noise that carries no information about the data.
The mixture avoids that noise and has no dependence on an imputation seed. The
simulation rows above use Bernoulli-$q$ A and nulls, not mixtures. The all-mixture
design itself was run on the original dnAging replicates (10 replicates, the
simulation where the former TE design also passed). There it rejected 0.8–4.0% of
tests at $\alpha=0.05$. That is valid but conservative, and it held under
sequential depletion wherever tests had at least 19 nulls (see
[All-mixture and depletion simulation](#all-mixture-and-depletion-simulation)).
It has **not yet** been run on the unpolarised or reference-haplotype replicates,
where the former design failed. That test, and the production-matcher negative
control, are blocking items in
[REMAINING_VALIDATION_PROPOSAL.md](docs/REMAINING_VALIDATION_PROPOSAL.md).

**Option: Bernoulli-$q$ hard orientation.** With `--asymmetric-polarity-null`, A and
every $B_i$ are hard-oriented by one Bernoulli-$q$ draw per site, and $B_0$ stays a
mixture:
$\Phi_i^0=\Phi_{\mathrm{SFS}}(B_{i,\mathrm{Bernoulli}(q)},B_{0,\mathrm{mix}})$.
ALT is declared derived when a reproducible coordinate-keyed uniform variate is less
than $q$, keyed by
`sha256(phi-sfs-bernoulli-q-v1, seed, chromosome, position)`, so a site receives the
same orientation wherever it is used. Results then condition on one
`--polarity-imputation-seed`.

#### Null calibration, Z-scores, and P-values

Finite site counts make the distance positive even under neutrality. Calibrate that
sampling floor separately for every focal category rather than comparing raw
Φ-SFS values across categories:

1. For a focal set $A$ containing $M$ variants, generate one reference SNP set $B_0$
   and $R$ additional SNP sets $B_1,\ldots,B_R$. The existing matcher generates all
   $R+1$ sets by the same procedure: each independently bootstraps the observed
   A-site ages and matches exactly $M$ SNPs to that bootstrap age CDF. Because
   disjoint matching draws later sets from a depleted pool, the sets are not
   identically distributed. Drawing $B_0$ at random (below) makes the choice of
   reference uniform and SFS-blind among the accepted sets. It does not make the
   depleted sets identically distributed, and it does not show that the
   $A$-versus-$B_0$ comparison is exchangeable with the $B_i$-versus-$B_0$
   comparisons; that remains an open validation item
   ([review 12](docs/CODE_REVIEW_ROUND12.md), finding 7). All sets use the shared samples,
   callability, and data-quality rules, with the single polarity rule above.
   The focal-age bootstrap resamples sites iid. This adopts the Poisson-random-field
   approximation that local LD averages out for genome-wide control pools containing
   millions of SNPs; local linkage is therefore a documented assumption, not a current
   production blocker. Revisit it for spatially restricted or strongly clustered sets.
2. Keep $B_0$ mixture-polarized and calculate the observed distance
   $\Phi_{\mathrm{obs}}=\Phi_{\mathrm{SFS}}(A,B_{0,\mathrm{mixture}})$, with A
   polarized as the nulls are.
3. Calculate $\Phi_i^0=\Phi_{\mathrm{SFS}}(B_i,B_{0,\mathrm{mixture}})$ for
   $i=1,\ldots,R$, with each $B_i$ a posterior mixture by default, or
   Bernoulli-$q$ hard-oriented under `--asymmetric-polarity-null`. Here $\Phi_i^0$
   is a raw Φ-SFS distance between two neutral SNP sets, not a Z-score.
4. Let $μ_0$ and $s_0$ be the mean and sample standard deviation of the
   $\Phi_i^0$. Report the standardized test statistic

   $$
   Z_A=\frac{\Phi_{\mathrm{obs}}-\mu_0}{s_0}.
   $$

5. Report the one-sided Monte Carlo P-value

   $$
   P_A=\frac{1+\sum_{i=1}^{R}
   \mathbf{1}\!\left(\Phi_i^0\ge \Phi_{\mathrm{obs}}\right)}{R+1}.
   $$

The matcher publishes 1001 disjoint sets. Phi-SFS draws the mixture-polarized $B_0$
uniformly from the
sets that pass matching QC, with a seed derived from `--reference-seed` and the
target digest, and uses every other QC-passing set as a null, so $R$ is whatever
passes QC and may differ between categories. A floor fixed before the run
(`--min-null-replicates`, 900 by default) guards against a coarse P-value. The
add-one P-value is valid for any such $R$. QC is computed only from age matching and does not directly inspect the
SFS. However, because allele age and allele frequency are related, selection on
age-matching QC is not guaranteed to be neutral with respect to the resulting SFS.
This is an inherent limitation of the selection scheme and should be considered when
interpreting calibrated results. All published sets must be globally disjoint.
Thus every control SNP has maximum reuse one and every $B_i$ has zero overlap with
$B_0$. $B_0$ is drawn before any SFS is examined, so the choice of reference is
SFS-blind; this does not by itself make the matched sets exchangeable (finding 7 of
[review 12](docs/CODE_REVIEW_ROUND12.md)). Replicate 0 is no longer the default $B_0$, because it is
matched first, from the undepleted pool, and so is not a typical set; it can still
be drawn as $B_0$ or used as a null like any other QC-passing set. Within the
calculation, $B_0$ differs from the other sets in being held fixed as the reference.

With $R$ near 1000 the minimum attainable P-value is about $1/(R+1)\approx10^{-3}$. Plot one equal-size point per focal category at
$Z_A$, color it by $-\log_{10}P_A$, and show its category-specific null Z-score
distribution in gray. Cap the displayed color scale at $\log_{10}(R+1)$, about 3.
Because $Z_A$ grows with $M$, this plot shows test strength, not effect size; do not
compare $Z_A$ across categories of different $M$ as an effect size.

![Illustrative category-specific null distributions, standardized Phi-SFS test statistics, and P-value colors](figures/phi_sfs_null_standardization_example.png)

The P-value tests whether a focal spectrum is farther from its matched neutral
background than expected from two finite neutral samples of the same size. The
Z-score is a standardized test statistic: it expresses that departure in
category-specific null standard deviations, and because $s_0$ is expected to shrink
as $M$ grows, it grows with $M$ for a fixed spectral difference. It is not a
magnitude: report $\Phi_{\mathrm{obs}}-\mu_0$ (in DAF units), the signed CDF and bin
residuals, and $M$ as the magnitude of a departure. Neither Z nor P identifies the direction of the SFS
shift; the CDFs and signed bin residuals do.

#### SNP type-I pilot

A preliminary empirical negative control treated 100 held-out real SNP sets as
Bernoulli-$q$-hard focal sets and compared each with a mixture-polarized SNP
reference using 344 Bernoulli-$q$-hard SNP null sets. Each set contained 4,023
sites. Seven of 100 tests rejected at $\alpha=0.05$, for an observed rejection
fraction of 0.07 and a descriptive Wilson 95% interval of 0.034--0.137.

This pilot is encouraging for the polarity construction but is not final pipeline
validation. All 100 tests shared one reference and one empirical null vector, so the
P-values are dependent. The sets came from an interrupted, unpublished matcher work
prefix, and none passed the production matching-error-ratio threshold: the median
ratio was 2.86 and 0/445 completed sets were below 0.5. The result must therefore be
reported as a preliminary SNP negative control, not evidence that the complete
matching-and-polarity pipeline has a precisely estimated 7% type-I error rate. Its
machine-readable provenance is in
`results/phi_sfs/snp_type1_asymmetric_100/summary.json`.

The pilot used the Bernoulli-$q$ design, not the current all-mixture default, and
SNP focal sets only. TE focal sets now share the SNP polarity and age construction,
so the pilot no longer misses a TE-specific polarity path. It still cannot test
properties specific to TE sites, such as TE genotyping error or how the ARG handles
TE sites ([review 12](docs/CODE_REVIEW_ROUND12.md), finding 2).

**Between-category contrasts are not currently supported.** The repository retains
`normalize_tes.phi_contrast` as experimental code, but its random pairing of category
nulls has not been validated and ignores covariance between nested or overlapping
categories. Do not interpret its P-values or use it in the production workflow.
For now, report each category's Phi-SFS result and uncertainty separately. A visual
difference between two points is not a formal between-category test.

#### All-mixture and depletion simulation

`tools/sim_polarity_arms.py` directly tests the current all-mixture estimand. For
both true and inferred ages it compares independently reusable controls with
sequentially depleted, globally disjoint controls. Every control set receives its
own iid bootstrap of the focal ages; the depleted arm removes selected controls
before matching the next bootstrap target and draws $B_0$ uniformly from the
completed sets. The simulator reports type-I error, Phi-SFS drift with control-set
order, the late-minus-early distance shift, matching error, and maximum reuse.

The depleted sampler is a binned diagnostic, not the production greedy CDF
optimizer. Its purpose is to determine whether the all-mixture statistic is
calibrated in this controlled setting and whether sequential depletion alone creates
an order effect. Simulation replicates are the independent units for uncertainty;
repeated focal draws within one replicate provide Monte Carlo precision but are not
counted as independent biological replicates. Cells with too few nulls to attain the
requested alpha are recorded as diagnostic-only and excluded from rejection-rate
summaries.

Run the full ten-replicate validation on a compute node. The simulator is not
checkpointed, so this forces Farm's non-preemptible `high` partition rather than
risk restarting an eight-hour `low` job from the beginning:

```bash
HPC_LOW=high HPC_HIGH=high HPC_CPUS=1 HPC_MEM=16G HPC_TIME=08:00:00 \
  ~/.claude/bin/hpc_run \
  'python -m tools.sim_polarity_arms \
    --prep-root results/sim_dnaging \
    --output results/sim_dnaging/all_mixture_depletion_v1 \
    --replicates 1 2 3 4 5 6 7 8 9 10 \
    --sizes 250 1000 4000 \
    --tests 50 --nulls 199 --mixture-sets 80'
```

The focal A set is observed once and remains fixed. Small or unusual focal sets can
therefore yield unstable results even after null calibration. Always report $M$,
the raw distance, null mean and standard deviation, Z-score, Monte Carlo P-value,
replicate count, and matching diagnostics. Results use the
`phi-sfs-wasserstein-v3` schema and cannot be silently combined with older Phi-SFS
outputs. The implementation design is recorded in
[PHI_SFS_WASSERSTEIN_CODING_PLAN.md](docs/PHI_SFS_WASSERSTEIN_CODING_PLAN.md).

## Farm/Quobyte launchers

Submit launchers with `sbatch` from the repository checkout. The launchers activate
the conda environment themselves. `$TMPDIR` is node-local scratch; matcher
`WORK_DIR` must remain on Quobyte or other durable storage.

Build the target and match controls:

```bash
sbatch --export=ALL,STORE="$STORE",TARGET="$TARGET",A_POSITIONS="$A_POSITIONS",\
OUTPUT="$MATCHES",CANDIDATE_ROWS="$CANDIDATES",WORK_DIR="$WORK_DIR",\
VCF_ELIGIBILITY="$VCF_ELIGIBILITY",A_TYPE="$A_TYPE",\
REPLICATES=1001,RESTARTS=3,SEED=1002,SCRATCH_HEADROOM_GB=32 \
  slurm/run_bootstrap_matching.sbatch
```

The launcher builds a missing target from `A_POSITIONS` with the A type and shared
eligibility artifact. When `TARGET` already exists, it verifies the recorded A type
and eligibility artifact and refuses a target built with the TE polarity mask. It
exits if `TE_POLARITY_MASK` or `MAX_FLIPPED_FRACTION` is set. It also rejects an
unrestricted candidate universe and verifies that the candidate provenance names the
same eligibility artifact.

Build the ancestral table as an array and merge it after every array task succeeds:

```bash
sbatch --array=0-14 --export=ALL,STORE="$STORE",\
TREES="$POSTERIOR_DIR/*.tsz",OUTPUT=results/ancestral-parts,PER_TASK=5 \
  slurm/run_ancestral_table.sbatch

sbatch --export=ALL,STORE="$STORE",MERGE=1,\
PARTS="results/ancestral-parts/part-*",OUTPUT="$ANCESTRAL",EXPECT_DRAWS=75 \
  slurm/run_ancestral_table.sbatch
```

Calculate Phi-SFS:

```bash
sbatch --export=ALL,TARGET="$TARGET",MATCHES="$MATCHES",VCF="$VCF",\
ANCESTRAL="$ANCESTRAL",OUTPUT="$PHI",A_TYPE="$A_TYPE",B_TYPE=SNP \
  slurm/run_phi_sfs.sbatch
```

This runs either TE-versus-SNP or SNP-versus-SNP according to `A_TYPE`; `B_TYPE`
is currently constrained to `SNP`, matching the command-line interface. The launcher
defaults to `ASYMMETRIC_POLARITY_NULL=false`, the all-mixture design. Set
`ASYMMETRIC_POLARITY_NULL=true` (with `POLARITY_IMPUTATION_SEED`, default 2001) for
the Bernoulli-$q$ option.

Scheduler allocations, measured resource use, scratch sizing, and parameter evidence
are recorded in [BOOTSTRAP_HPC_VALIDATION.md](docs/BOOTSTRAP_HPC_VALIDATION.md).

### Submit many TE categories

The store, candidate universe, and ancestral table are shared across categories.
Give every category its own target, matched bundle, durable work directory, and
seed. First create a tab-separated manifest:

```text
label	positions	target	matches	work_dir	seed
all_te	/quobyte/project/te/all.pos.txt	/quobyte/project/targets/all_te	/quobyte/project/matches/all_te	/quobyte/project/work/all_te	1001
in_gene	/quobyte/project/te/in_gene.pos.txt	/quobyte/project/targets/in_gene	/quobyte/project/matches/in_gene	/quobyte/project/work/in_gene	1002
young	/quobyte/project/te/young.pos.txt	/quobyte/project/targets/young	/quobyte/project/matches/young	/quobyte/project/work/young	1003
```

The manifest rules are:

- the first line is the header shown above;
- fields are separated by literal tabs and paths must not contain tabs, newlines, or
  commas;
- every category-specific path must be unique; target and match paths must not
  exist on a first submission;
- `work_dir` is durable and may be reused only to resume the identical matching run;
- seeds should be fixed before submission and remain unchanged on resubmission.

Set the shared inputs, then submit one target/matching job per manifest row:

```bash
PROJECT=/quobyte/project/PhiTE
STORE=/quobyte/project/data/age_interval_store
CANDIDATES=/quobyte/project/data/candidate_rows.npy
VCF_ELIGIBILITY=/quobyte/project/data/vcf_eligibility
MANIFEST=/quobyte/project/manifests/te_categories.tsv

while IFS=$'\t' read -r label positions target matches work seed; do
  [[ "$label" == label ]] && continue
  [[ -n "$label" ]] || continue

  match_job=$(sbatch --parsable \
    --job-name="match-${label}" \
    --export=ALL,PROJECT="$PROJECT",STORE="$STORE",TARGET="$target",\
A_POSITIONS="$positions",A_TYPE=TE,OUTPUT="$matches",CANDIDATE_ROWS="$CANDIDATES",\
VCF_ELIGIBILITY="$VCF_ELIGIBILITY",WORK_DIR="$work",\
REPLICATES=1001,RESTARTS=3,SEED="$seed",SCRATCH_HEADROOM_GB=32 \
    slurm/run_bootstrap_matching.sbatch)
  match_job=${match_job%%;*}

  printf '%s\tmatch=%s\n' "$label" "$match_job"
done < "$MANIFEST"
```

The categories run concurrently. Save the printed job IDs and use `squeue`,
`sacct`, and the scheduler logs to confirm that every manifest row completed. The
loop does not silently skip existing outputs; for a partial rerun, submit only the
missing categories or resubmit an interrupted matcher with its original target,
output, work directory, and seed.

## Verify a production run

Before accepting the results:

1. Confirm every artifact records the expected release version, Git commit, and
   non-null input identities in `metadata.json`.
2. Confirm the candidate-row report meets the requested resolution threshold and is
   bound to `STORE` and the intended VCF eligibility artifact.
3. Confirm the target records the intended `a_type`, the eligibility artifact, and
   plausible kept/removed counts, and records no `te_polarity` mask.
4. Confirm the matcher published 1001 identically generated sets in disjoint mode,
   maximum control reuse is one, every overlap with $B_0$ is zero, and all sets used
   in Phi-SFS pass matching QC.
5. Confirm the `phi-sfs-wasserstein-v3` result records the intended A/B types,
   `null_polarity_design=posterior-mixture-vs-posterior-mixture` (or, for the
   Bernoulli option, `bernoulli-q-hard-vs-posterior-mixture` with its seed and
   algorithm), identical A and null-left polarity rules, exactly equal site count $M$, reference replicate, null count,
   raw distance, null mean and sample SD, Z-score, exceedances, and add-one P-value.

The exact acceptance criteria and the tests supporting them are in
[BOOTSTRAP_HPC_VALIDATION.md](docs/BOOTSTRAP_HPC_VALIDATION.md).

## Outputs

| artifact | purpose |
|---|---|
| `age_interval_store/` | reusable posterior age intervals and store identity |
| `ancestral_states/` | posterior ancestral-base counts for store rows |
| `vcf_eligibility/` | shared callable rows plus the ARG-orientable subset and its posterior orientation |
| `candidate_rows.npy` plus `.json` | store-bound, eligibility-filtered control universe excluding known TEs and A |
| `targets/CATEGORY/` | focal age target, from all posterior draws, and acceptance threshold |
| `bootstrap_matches/CATEGORY/` | disjoint $B_0,\ldots,B_R$ sets, bootstrap targets, restart traces, reuse checks, and QC |
| `phi_sfs/CATEGORY/` | A/B spectra and CDFs, Wasserstein distances, null Z-scores, summary tables, and provenance |

Outputs are published atomically and are never overwritten. The Phi-SFS output schema
is `phi-sfs-wasserstein-v3`; its generic `a_*` and `b_*` arrays support both TE-SNP
and SNP-SNP analyses, and under `--asymmetric-polarity-null` its Bernoulli-$q$ arrays
record the hard null-left spectra. Version 3 polarizes TEs from the ARG posterior;
version 2 gave TEs hard biological polarity after a derived-support filter, and
version 1 used the earlier symmetric-mixture null. None may be silently combined
with another. Matched-control sets are Monte Carlo null replicates, not independent
biological samples; see the validation report for the correct interpretation of
their spread.

## Methods and validation

- [BOOTSTRAP_HPC_VALIDATION.md](docs/BOOTSTRAP_HPC_VALIDATION.md) — production settings,
  validation tests, measured resources, decision evidence, and acceptance criteria.
- [BOOTSTRAP_TARGET_MATCHING_PLAN.md](docs/BOOTSTRAP_TARGET_MATCHING_PLAN.md) — bootstrap
  target and matching design.
- [PHI_SFS_WASSERSTEIN_CODING_PLAN.md](docs/PHI_SFS_WASSERSTEIN_CODING_PLAN.md) —
  Wasserstein definition, polarity paths, null calibration, and output schema.
- [BOOTSTRAP_DISCARDED_APPROACHES.md](docs/BOOTSTRAP_DISCARDED_APPROACHES.md) — evaluated
  approaches that are not part of the production route.
- [CHANGELOG.md](docs/CHANGELOG.md) — release-level behavior changes.
- [CODE_REVIEW_ROUND12.md](docs/CODE_REVIEW_ROUND12.md) — latest review (Claude and
  Codex): statistical design of the Phi-SFS calibration, the matching-QC blocker, and
  the proposed simulation validation study.

Historical `INTERVAL_STORE_*`, `GLOBAL_QUANTILE_*`, sampler plans, and older code
reviews document development history; they are not operator instructions.
