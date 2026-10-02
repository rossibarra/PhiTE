# Code Review — Round 10

Date: 2026-09-25  
Reviewed revision: `d07c333` ("Implement calibrated Wasserstein Phi-SFS", `bootstrap-target-hpc-validation`)  
Scope: `96a8611..d07c333` — the Wasserstein Phi-SFS overhaul, measured against
[PHI_SFS_WASSERSTEIN_CODING_PLAN.md](PHI_SFS_WASSERSTEIN_CODING_PLAN.md).

Files reviewed: `normalize_tes/phi_sfs.py`, `normalize_tes/vcf_eligibility.py` (new),
`normalize_tes/te_age_target.py`, `normalize_tes/build_candidate_rows.py`,
`normalize_tes/bootstrap_target_matcher.py`, `slurm/run_bootstrap_matching.sbatch`,
`slurm/run_phi_sfs.sbatch`, the four affected test modules, and README sections 3–8.

## Summary

The numerical core is correct and well tested. `phi_sfs()` computes the discrete
Wasserstein-1 distance exactly as the plan specifies, `calibrate_phi()` uses `ddof=1`,
the `>=` tail rule and the add-one P-value, and the vectorized null distances in
`calculate()` agree with the pure function. The polarity mixture
`q·h(k,n) + (1−q)·h(n−k,n)` is kept, and SNP-versus-SNP is supported end to end.
Most of the plan's guards are present: a disjoint bundle is required, maximum reuse
must be 1, rows in A may not appear in B, overlap with B0 must be 0, site counts must
be equal, and the reference is fixed before the SFS scan. The full suite passes
(292 tests, run on a compute node).

The pipeline is not ready for a production run. Three problems come first:

1. `phi_sfs` accepts the A type from its command line and never checks it against the
   target. The launchers' defaults make this mismatch easy to hit, and it produces a
   silently wrong result (Finding 1).
2. The default design (1001 sets published, at least 1000 null sets required) only
   completes if every replicate passes matching QC. The recent disjoint run passed
   95/100, so I expect a default run to fail at the last step (Finding 2).
3. The new disjoint capacity preflight checks the total pool size, not how many
   candidates each age stratum holds. The project's own record puts the in-gene target's
   scarcest decile at about 787 sets' worth, which is fewer than 1001 (Finding 3).

## Decisions after review (2026-09-25)

| finding | decision | status |
|---|---|---|
| 1 | Both `phi_sfs` and the matcher check the A type, the TE `max_flipped_fraction` of 0.5 and the eligibility record against target (and match) metadata before any work is done. | implemented |
| 2 | Revised 2026-09-26: publish 1001 disjoint sets, draw B0 at random (seeded, before the SFS scan) from the QC-passing sets, and use every other QC-passing set as a null, so R varies by category above a floor of 900. The earlier Option A (1201 sets, R fixed at 1000) was dropped because the in-gene target's youngest ages hold only about 1,045 disjoint sets' worth of candidate age mass. | implemented |
| 3 | A per-stratum disjoint capacity preflight runs before any matching work. It now compares candidate and target age **mass** per stratum; the first version compared median counts with mass quotas and reported about 357 sets for the in-gene target instead of about 1,045. The depletion re-measurement is still outstanding. | implemented; depletion check open |
| 4 | The eligibility artifact is identified by content: VCF digest, array digests, heterozygous policy and `min_callable`. The identity is recorded in the target, the candidate report and the match metadata. It is checked by the matcher, the launcher and `phi_sfs`, which also checks its VCF digest against it. | implemented |
| 5 | Skipped by decision. | won't fix |
| 6 | TEs with no usable draws are discarded and counted in `sites_with_no_usable_draw`. The existing tests cover this. | implemented |
| 7 | Measured with `tools/benchmark_phi_sfs_scale.py` (synthetic, 1201 sets): peak 2.5 GiB and about 1.4 min at M = 4,067; 11.2 GiB and about 4 min at M = 19,000. This is within 48G / 6 h. The VCF scan and ancestral-table memory were not measured, so a real end-to-end run is still needed. | partly done |
| 8 | Implemented: `--reference-sensitivity`, per-row reuse columns, the byte-reproducibility and injected-failure tests (the scipy cross-check is written but skipped because scipy is not in the environment), `normalize_tes.phi_contrast`, the calibration simulation study ([PHI_SFS_CALIBRATION_VALIDATION.md](PHI_SFS_CALIBRATION_VALIDATION.md)), and the CHANGELOG entry and v0.8.0 bump. The study finds the test correctly sized under an exchangeable null. In stylised models, however, the TE-versus-SNP polarity difference produces false signals. | implemented |
| 9 | VCF reading and genotype decoding are in `normalize_tes.vcf_io`, shared by the eligibility scan and `phi_sfs`. The scan parses CHROM and POS first. | implemented |
| 10 | Removed the dead schema entry and the duplicate row loads, fixed the docstrings, and split the eligibility exclusion counts. | implemented |

## Findings

### 1. High — `phi_sfs` trusts `-A` and never checks it against the target

`calculate()` chooses how to polarize A from `args.a_type` alone
(`phi_sfs.py:905–906`). It never reads `a_type`, `te_polarity`, or `vcf_eligibility`
from the target metadata. `te_age_target` now records all three
(`te_age_target.py:804–806`). The only `target_meta` fields `phi_sfs` reads are the
store digests and `schema_version`.

Failure scenario: a user builds a SNP target (`A_TYPE=SNP` in
`run_bootstrap_matching.sbatch`) and then runs `run_phi_sfs.sbatch` without also
setting `A_TYPE`. That launcher defaults to `TE` (`run_phi_sfs.sbatch:65`). Every A
SNP is then added to `PolarityResolver`'s TE set and gets `p_alt_derived = 1.0`, so ALT
is treated as derived. The run completes and publishes a complete, plausible result
with the wrong A polarity. The metadata then states that the "at-least-50%-derived
retention" rule was applied (`phi_sfs.py:1138–1142`), although nothing checked it.

The mismatch in the other direction (a TE target run with `-A SNP`) probably fails
loudly, on a non-ACGT TE allele or a missing orientation. I have not verified this.

The test suite runs the same target bundle once with `-A TE` and once with `-A SNP`
(`test_end_to_end_metadata_and_diagnostics`, `test_end_to_end_snp_a_uses_posterior_polarity`),
and both pass. This confirms that `phi_sfs` does not check the target's type.

Recommended correction: make the target the authority. Require
`target_meta["a_type"] == args.a_type` (or drop `-A` and read it from the target).
For `TE`, require `te_polarity` to be present with `max_flipped_fraction <= 0.5`.
For both types, require a non-null `vcf_eligibility` record. Write the
`a_polarity_rule` metadata only after these checks pass. Add a test that a SNP target
run with `-A TE` fails.

### 2. High — the default replicate design requires 100% matching-QC pass

The matcher publishes 1001 sets by default (`bootstrap_target_matcher.py:43`,
`run_bootstrap_matching.sbatch:79`). `phi_sfs` requires
`--min-null-replicates 1000` by default, counted after B0 is reserved
(`phi_sfs.py:762–768`, `:1239`; `run_phi_sfs.sbatch:68`). A default run therefore
succeeds only if all 1001 sets pass QC. Sets that fail QC are not redrawn.

Historical pass rates of existing bundles (read from their `qc_pass.npy`):

| bundle | QC pass |
|---|---|
| `in_gene_75draw_disjoint` | 95/100 |
| `in_gene_75draw` | 96/100 |
| `in_gene_75draw_restarts6` | 97/100 |
| `in_gene_75draw_logage_disjoint`, `_loggrid_disjoint`, `init_stratified` | 100/100 |

Prediction, assuming failures are independent: at a 5% per-set failure rate, the
chance that all 1001 sets pass is about e⁻⁵¹. At 1% it is about e⁻¹⁰. Even the 100/100
bundles are consistent with a failure rate of a few percent. Pool depletion
(Finding 3) would push the failure rate up in later replicates. If these rates carry
over, a default production run would spend about a day matching and then stop in
`phi_sfs`.

Recommended correction: settle the rule before any SFS is examined. One option is to
publish a margin, for example 1001 + K sets, and take B0 plus the first 1000 QC-passing
null sets in replicate order. QC is SFS-blind, so this choice is outcome-independent.
Another is to set R to the number of passing sets, with a prespecified floor. Either
way, record the rule in the metadata and make the matcher default, the `phi_sfs`
default and both launchers agree.

### 3. High — the disjoint preflight checks total capacity, not stratum capacity

The new preflight (`bootstrap_target_matcher.py:909–919`) requires
`replicates × M <= |candidates|`. The per-replicate guard (`:1031–1039`) only checks
that at least M candidates remain. Neither check considers where in the age
distribution the candidates are.

Commit `d2b5774`, which introduced disjoint mode, gave this justification:
"the scarcest decile holds ~787 sets' worth of candidates against the 100 needed".
This commit removed that comment. At 1001 sets, that decile would run out at around
replicate 787, even though the total pool is only 17.7% consumed. The 787 figure was
measured on the pre-eligibility candidate pool. The VCF-eligible pool is a subset of
it, so the true capacity is probably lower. I have not verified this.

When a stratum runs out, `stratified_initial_set` silently fills the gap from other
strata (`bootstrap_target_matcher.py:184–187`). The optimizer then does the best it
can with a pool that is thin in that age range. The expected result is that later
replicates match progressively worse: they either fail QC (feeding Finding 2) or pass
with matching error that trends upward. B0 is replicate 0, the best-supplied set drawn
from the full pool, so a trend in the null sets makes the B_i-to-B0 null distances
depend on replicate index.

The plan made re-measuring depletion at the production replicate count a
precondition (§2.2, §9.5, acceptance criterion 8). That re-measurement has not been
done.

Recommended correction: extend the preflight to check capacity per stratum, using the
target's `interval_boundary_ages`/`interval_quotas` and each candidate's median-age
stratum (the same assignment `stratified_initial_set` uses). Fail if any stratum holds
fewer than `replicates × quota` candidates. Then run the planned depletion check:
`E_r` and `D_i^0` against replicate index at R+1 = 1001, on a VCF-eligible pool.

### 4. Medium — the eligibility artifact is identified by path, and nothing binds it to the VCF used by `phi_sfs`

Candidate rows and the target record the eligibility artifact only as a resolved path
(`build_candidate_rows.py`, `"vcf_eligibility": str(args.vcf_eligibility.resolve())`;
`te_age_target.py`, `eligibility_report["mask"]`). The only cross-check is in the
launcher, which compares these paths (`run_bootstrap_matching.sbatch:147`, `:219–226`).
The matcher CLI, which is what README step 7 invokes directly, performs no check. The
v0.7.0 changelog retired path equality as an identity test because a path identifies a
location, not the content at it. This new artifact reintroduces that pattern.

Downstream, `phi_sfs` computes `vcf_sha256` during its scan but never compares it with
the mask's `vcf_sha256`. It also never checks that its `--heterozygous` policy matches
the mask's. The README asks the operator to keep these consistent by hand ("Use
`--heterozygous missing` only when the eligibility artifact was built with the same
policy").

Most mismatches fail loudly through the equal-M assertion. A different VCF in which
every requested site is still callable does not. The published result then comes from
a VCF other than the one that defined eligibility, and the metadata does not show it.

Recommended correction: record the mask's `array_sha256["row_indices"]` (or a digest
of the whole artifact), its `vcf_sha256`, `heterozygous` and `min_callable` in both the
candidate report and the target metadata. Have the matcher require that they agree.
Have `phi_sfs` require that its own VCF digest and heterozygous policy match.

### 5. Medium — allele validation uses a substring test (skipped by decision)

Both places that validate SNP alleles test membership in a string, not in a set of
single bases:

- `vcf_eligibility.py:253`: `if ref not in "ACGT" or alt not in "ACGT"`
- `phi_sfs.py:504`: `if ref not in self.BASES or alt not in self.BASES` with `BASES = "ACGT"`

`"AC" in "ACGT"`, `"GT" in "ACGT"` and `"CGT" in "ACGT"` are all `True` (checked).
`"ACGT".index("GT")` returns 2, so a multi-base allele is silently read as its first
base. For example, a REF=`GT`, ALT=`G` deletion would be looked up as G versus G. It
would then pass as an orientable SNP with a meaningless q.

I have not measured the impact. The candidate universe is restricted to the filtered
SNP position list, which probably excludes indels. However, `scan_vcf` evaluates every
callable catalog row, and `-A SNP` targets come from a user-supplied list.

Recommended correction: use `frozenset("ACGT")` (or require `len(allele) == 1`) in
both places, and add a test with a multi-base allele.

### 6. Medium — TEs with no usable draws pass the "at least 50% derived" rule

`load_polarity_selection` keeps every TE with `usable == 0`, regardless of the
threshold (`te_age_target.py:435`). The plan's rule is to retain an insertion "only
when it is inferred to be derived in at least 50% of the ARG draws". The launcher's
error text ("the inclusive >=50% derived rule", `run_bootstrap_matching.sbatch:99`),
README step 6 and the `phi_sfs` metadata all describe the rule without this exception.

The exception predates this change. What is new is that it now sits inside the
declared Phi-SFS rule. The report does not count these sites separately:
`sites_with_no_agreeing_draw` is computed before thresholding and also includes sites
that the threshold then discards.

Recommended correction: decide whether zero-evidence TEs belong in A. If they do,
document the exception wherever the rule is stated and report a
`sites_kept_without_usable_draws` count. If they do not, drop them under the rule.

### 7. Low — resource figures at 1001 replicates are extrapolated, not measured

- The matcher's `--time=1-12:00:00` rests on "about 22 h", which is one 100-replicate
  measurement scaled ×10 (`run_bootstrap_matching.sbatch:58–61`). That is a single-point
  linear extrapolation. It leaves out the per-replicate `np.isin` over a growing claimed
  set and any extra optimizer work caused by depletion.
- `run_phi_sfs.sbatch` still requests 48G and 6 h. Memory "scales with the number of
  requested sites" (its own comment), and that number has grown about 10×.
  `calculate()` holds (R+1)·M coordinate tuples twice (`all_snp_coordinates` and
  `requested`), plus a `SiteCount` dict and a row map. My rough estimate, not measured,
  is 10–12 GB at about 23 k sites × 1001 sets. That fits, but the figure should be
  measured.
- Each result also republishes `b_bootstrap_counts.npy`, a (1001, M) int64 array that
  is already present in the match bundle.

Acceptance criterion 10 of the plan requires a representative end-to-end run inside
the resource envelope. That run has not been done.

### 8. Low — plan items not yet implemented

These are consistent with the plan's staged sequence (§10, steps 6–8), but the
acceptance criteria are not yet met:

- `--reference-sensitivity` and the reference-sensitivity range in `summary.csv`
  (metadata hard-codes `reference_sensitivity_run: false`);
- the between-category contrast module (§8);
- the statistical calibration validation (§9.5: P-value uniformity, type-I error, and
  depletion under disjoint sampling);
- tests for: agreement with `scipy.stats.wasserstein_distance`, byte-reproducible
  arrays from fixed inputs, and cleanup of the staging directory after an injected
  failure;
- reuse diagnostics per row of `comparisons.csv` (only overlap with B0 is written);
- a CHANGELOG entry and version bump for the breaking schema change (criterion 9).

### 9. Low — duplicated VCF reading code that must stay in step

`vcf_eligibility.py` contains its own copies of `_HashingStream`, `_open_vcf` and the
genotype decoder from `phi_sfs.py`. The equal-M invariant depends on the two decoders
making identical callability decisions. They agree today, but nothing forces them to
stay that way. Move them into one shared helper.

Separately, `scan_vcf` splits every sample column of every record
(`vcf_eligibility.py:184`), including records that are not in the catalog. `phi_sfs`
deliberately avoids this because it dominates scan time on a sample-rich VCF. Parsing
only CHROM and POS first, as `read_site_counts` does, would bring the same saving.

### 10. Low — smaller cleanups

- The `swap-age-matched-controls-v1` entry in `MATCH_IDENTIFIERS` (`phi_sfs.py:57–60`)
  is now dead code, because `calculate()` rejects every schema except
  `bootstrap-target-matches-v1`.
- `row_indices.npy` and `te_row_indices.npy` are loaded twice (inside
  `_load_coordinates` and again at `phi_sfs.py:774–779`).
- `excluded_by_reason` in the eligibility report mixes rows removed from the callable
  mask with rows removed only from the SNP-orientable subset (`snp_non_acgt_alleles`,
  `snp_no_usable_orientation`). The second group remains in `row_indices.npy`.
  Reporting the two groups separately would make `eligible_rows` reconcile with the
  reasons.
- The `SiteCount` and `PolarityResolver` docstrings still describe SNP polarity as
  applying to "control" sites only. With `-A SNP` it applies to A as well.

## Design observations (not code defects)

These follow from the approved plan and are recorded so that the calibration work in
§9.5 tests them explicitly:

- **Asymmetric error budget.** D_obs = W(A, B0) contains one focal-bootstrap draw
  and one matching error. Each null D_i = W(B_i, B0) contains two of each. Under a true
  null I therefore expect the null distances to be slightly larger than the observed
  one, which would make the test conservative. The simulation study should measure
  how large this effect is.
- **TE and SNP polarity differ by construction.** TEs contribute a point-polarized
  projection (weight 1). SNPs contribute a q-mixture, with ARG orientation about 91%
  accurate against TE ground truth. This gap in method can by itself separate the TE
  and SNP spectra, so the `-A SNP` negative control is the right comparison for
  interpreting a TE-versus-SNP Z-score.

## Verification performed

- Full test suite on a compute node through `hpc_run`: **292 passed**, 4 unrelated
  fork deprecation warnings.
- Read all of `phi_sfs.py` and `vcf_eligibility.py`, and the full diffs of the other
  changed modules and both launchers.
- Checked the substring behaviour in Finding 5 directly in Python.
- Counted QC pass rates from `qc_pass.npy` in every bundle under `results/bootstrap_matches/`.
- Confirmed with grep that `phi_sfs.py` reads no `a_type`, `te_polarity`, or
  `vcf_eligibility` field from the target metadata.
- Not run: any end-to-end production-scale job. The resource figures in Finding 7 and
  the failure-rate projection in Finding 2 are estimates, labelled as such.
