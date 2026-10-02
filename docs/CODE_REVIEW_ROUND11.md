# Code Review — Round 11

Date: 2026-09-25
Reviewed revision: `05f0d41` (v0.8.0 release commit on
`bootstrap-target-hpc-validation`)
Scope: committed changes from `v0.7.0` through `05f0d41`, with emphasis on the
feasibility and statistical validity of disjoint SNP age matching. This review also
reassesses the open items recorded in
[CODE_REVIEW_ROUND10.md](CODE_REVIEW_ROUND10.md).

The current worktree contains an uncommitted extension to
`tools/validate_phi_calibration.py`; it is not part of revision `05f0d41` and is not
treated as released behavior here.

## Summary

The v0.8 numerical implementation is substantially stronger than v0.7. The discrete
Wasserstein calculation, deterministic projection to 20 individuals, add-one
P-value, fixed null count, content-addressed VCF eligibility, A-type checks, equal
site-count checks, atomic publication, and provenance guards are all sensible. Most
of the direct software defects found in Round 10 were corrected.

The production inferential design is not ready. The main blocker is not the total
number of SNPs. For the legacy in-gene target, 1,201 sets of 4,067 sites require
4,884,467 controls, well below the 23,026,051-row legacy candidate pool. The problem
is the combination of age-specific scarcity and global disjointness. A capacity scan
on that legacy target/pool estimated only about 357 target quotas in the scarcest
median-age stratum. However, that scan predates v0.8 VCF eligibility and does not by
itself prove that 1,201 good aggregate-CDF matches are mathematically impossible on
the final inputs. It proves only that the current exact-quota preflight rejects that
legacy combination.

More importantly, even a passing preflight would not establish validity. The matcher
generates sets sequentially from changing candidate pools, so replicate 0 and later
replicates are not identically distributed. The add-one Monte Carlo P-value is
validated only under an exchangeable null, while the finite-pool depletion regime is
explicitly untested. The current preflight also does not enforce the quota assumptions
it uses: initialization can back-fill a deficient sampled stratum, and optimization
can swap freely across strata.

## Findings

### 1. High — sequential global disjointness breaks the claimed exchangeability

In disjoint mode, the matcher publishes replicate 0 from the full candidate pool,
claims its selected rows, removes them, and repeats. Every replicate therefore sees a
different candidate universe. The comments in `bootstrap_target_matcher.run()`
correctly acknowledge this dependence, but the README and coding plan describe the
sets as identically generated and say that B0 differs only by being held fixed.

That claim is false for the implemented algorithm:

- B0, normally replicate 0, receives the largest and least-depleted pool;
- later sets are generated after increasingly many selected rows have been removed;
- the first R QC-passing sets are preferentially early sets, while unused spares are
  preferentially late, depleted sets; and
- any age scarcity or optimizer degradation can therefore make matching error and
  SFS distance depend on replicate ID.

The add-one formula is exact under the exchangeability assumptions used in the
simulation study. The study explicitly did not simulate sequential finite-pool
depletion. Until exchangeability or calibrated type-I error is demonstrated for the
real matching process, the reported P-value should not be described as calibrated.

Recommended correction: settle the null design before inspecting any production SFS.
Reasonable candidates include:

1. make every Bi use the same candidate universe after excluding B0, while allowing
   reuse among the Bi; then measure overlap, dependence, and effective Monte Carlo
   size rather than requiring global maximum reuse one;
2. construct a globally disjoint collection jointly and randomly assign the reference
   and null labels only after the collection is complete, so generation order is not
   confounded with B0; or
3. use a smaller prespecified R where a fully disjoint design is feasible, report its
   coarser attainable P-value, and validate depletion at that R.

Whichever design is chosen needs a neutral simulation or negative-control experiment
that includes the actual finite candidate pool, matcher, QC selection, and reference
selection—not only `calibrate_phi()` applied to i.i.d. spectra.

### 2. High — the stratum preflight is neither necessary nor sufficient for matching feasibility

The new preflight assigns every candidate to the stratum containing its median age and
requires `candidates_in_stratum >= replicates * target_quota`. This is a useful warning
about obvious local scarcity, but the code and documentation give it a stronger
interpretation than it supports.

It is not sufficient because the run does not preserve those quotas:

- `stratified_initial_set()` first samples only `init_oversample * M` candidates from
  the whole remaining pool. If that sample underrepresents a scarce stratum, it fills
  the missing positions from other strata without reporting the deviation.
- `optimize_restart()` proposes replacements from the entire replicate candidate
  universe. It does not require an accepted swap to remain in the removed row's
  median-age stratum.
- The selected optimized set, not its quota-conforming initial state, is claimed from
  the global pool. Early sets can therefore consume more than the assumed quota from a
  scarce stratum, leaving less than the preflight predicted for later initializations.

It is not necessary because the objective matches the aggregate age CDF, not a vector
of exact median-age-stratum counts. Posterior age distributions span boundaries; a set
that violates the median-bin quota can still be a good aggregate-CDF match. Failure of
the quota preflight proves that the current implementation refuses the run, not that
no scientifically acceptable 1,201-set solution exists.

Recommended correction: first decide whether exact stratum quotas are a real design
constraint or only an initialization heuristic.

- If they are a constraint, cache the candidate stratum assignments from the
  preflight, draw directly within strata, restrict or account for cross-stratum swaps,
  and verify the selected set's stratum counts before claiming it.
- If they are only a heuristic, report capacity as a diagnostic rather than a hard
  proof of infeasibility, and use pilot matching plus certified W1/QC trends to decide
  feasibility.

Add a test in which optimization changes a selected set's stratum composition and
show that later initialization either remains valid or fails before partial work is
accepted. The current tests cover static counts but not this consumption path.

### 3. High — 1,201-set feasibility has not been measured on the final v0.8 inputs

The 23,026,051-row `results/candidate-rows-75draw.npy` and the legacy
`in_gene_75draw` target predate the v0.8 VCF-eligibility artifact. The measured
approximately 357-set stratum capacity therefore does not describe the exact final
target/candidate pair required by the README.

VCF filtering can only remove rows from a fixed candidate array, but the final focal
filter can also change M, age boundaries, and quotas. Consequently, one cannot infer
the final stratum capacities merely by applying the old 357 figure. The exact check
requires:

- the final target after callability and TE-polarity filtering;
- the final SNP candidate array intersected with the same eligibility artifact; and
- the same interval store used to build both.

There is also a scale contradiction in the operator guide. It recommends 1,201 sets
for every category, including large categories. Any target with
`1201 * M > candidate_count` fails the total-capacity guard regardless of age
strata. Category-specific feasibility must therefore precede submission, and R or the
null design cannot be a universal constant unless the smallest candidate-to-target
ratio supports it.

Recommended correction: publish a read-only capacity-report command that records the
target and candidate identities, total capacity, per-stratum diagnostic capacity, and
an explicit statement that the latter is a quota-model bound rather than a proof of
aggregate matching impossibility. Run it on every final category before launching the
matcher.

### 4. Medium — the formal category contrast rests on an unvalidated arbitrary pairing

`phi_contrast` combines two unrelated null vectors by permuting one vector once and
pairing it elementwise with the other. Repeat 0 supplies the formal P-value and the
Holm/BH corrections; other random pairings are reported only as a sensitivity range.

This procedure assumes that each category's null entries are exchangeable. That is
already unresolved under sequential depletion. Even if exchangeability held, the
reported primary result depends on an arbitrary pairing seed and uses only R of the
R-squared possible empirical differences. The range across repeated pairings is
evidence that pairing uncertainty exists, not a validation of repeat 0.

Recommended correction: validate the complete contrast procedure under paired and
independent-category simulations that reproduce matching, QC, shared inputs, and
multiple testing. Prespecify whether categories have a genuine shared replicate
identity. If they do not, use a method for independent empirical null distributions
rather than presenting an arbitrary bijection as a formal pairing.

### 5. Low — production resource and end-to-end validation remains incomplete

Round 10 Finding 7 is only partly closed. Synthetic Phi-SFS benchmarks fit within the
48 GiB / 6 h request, but they did not measure the complete VCF scan and ancestral
table path. Matcher time for 1,201 sets is a linear extrapolation from one 100-set run.
No representative final v0.8 run has demonstrated:

- capacity on the exact eligible target and candidate pool;
- 1,201-set matcher completion within the 96 GiB / 36 h allocation;
- acceptable QC pass rate and no replicate-index trend;
- end-to-end Phi-SFS completion on the real genome-wide VCF; or
- calibrated negative-control behavior.

This remains acceptance criterion 10 of the coding plan and should be completed on a
Farm compute node. The matcher is resumable, but a resumable extrapolation is not a
resource measurement.

## Round 10 carry-forward status

| Round 10 item | Round 11 status |
|---|---|
| Finding 1: A type trusted from CLI | Closed. Matcher and Phi-SFS validate target authority. |
| Finding 2: exactly 1,001 sets required 100% QC pass | Software behavior closed by 200 spares and fixed-R selection. Its extra capacity cost and selection behavior remain part of Findings 1 and 3 above. |
| Finding 3: no per-stratum capacity preflight | Static preflight implemented, but the final eligible inputs were not measured and the preflight assumptions are not enforced; Findings 2 and 3 remain open. |
| Finding 4: eligibility identified only by path | Closed by content identities and VCF digest checks. |
| Finding 5: substring allele validation | Closed by upstream prefiltering that removes indels and multiallelic variants. |
| Finding 6: TEs with no usable draws retained | Closed; such TEs are discarded and counted. |
| Finding 7: production resource figures extrapolated | Partly closed by synthetic Phi-SFS benchmarks; real end-to-end validation remains open as Finding 5. |
| Finding 8: missing plan features and calibration | Most software features were implemented. Finite-pool depletion calibration remains absent, and the SciPy cross-check remains skipped in the declared environment. |
| Finding 9: duplicated VCF reading | Closed through `normalize_tes.vcf_io`. |
| Finding 10: cleanups | Closed as recorded in Round 10. |

One unresolved Round 10 design observation also remains material: the observed
comparison and null comparisons have asymmetric bootstrap and matching error budgets.
Simulations suggest conservatism in stylised cases but do not measure the real
pipeline.

## Site-count policy clarification

v0.8 deliberately rejects post-hoc downsampling of every set to the smallest surviving
site count. This does not remove the deterministic projection of each site from n
callable individuals to m=20; it removes only after-the-fact downsampling of the number
of sites.

The intended order is:

1. determine callability and SNP orientability;
2. filter the focal set and candidate pool;
3. fix focal site count M and its age target; and
4. match exactly M eligible controls in every set.

For data in which every SNP has at least 20 callable samples, VCF eligibility should
remove no SNP for callability. It can still remove SNPs for absent records, malformed
genotypes, or zero usable ARG orientation. Uncertain or incorrectly oriented SNPs are
not excluded when q exists; their uncertainty is represented by the q-mixture. A
smaller M remains possible, but it must be chosen and applied before target construction
with a prespecified seed or stratified rule—not after observing matching or SFS results.

## Verification performed

- Read Round 10, the v0.7-to-v0.8 diffs, the final matcher, Phi-SFS, eligibility,
  category-contrast, launchers, tests, README, coding plan, and calibration report.
- Confirmed that the legacy candidate report contains 23,026,051 rows and predates
  the v0.8 eligibility identity.
- Confirmed that the legacy in-gene target has M=4,067, so total demand at 1,201 sets
  is 4,884,467 rather than more than 20 million.
- `python -m compileall -q normalize_tes tools` succeeded in the repository conda
  environment.
- A fresh full test run was attempted through the required Farm compute-node helper,
  but Slurm controller discovery failed with a DNS resolution error before a job was
  allocated. No current test result is claimed. Round 10's 292-passing result was for
  an earlier reviewed revision.
- Not run: the capacity scan on a final v0.8 target/candidate pair, production-scale
  matching, depletion analysis, or real end-to-end Phi-SFS.

## Release recommendation

Do not use v0.8 P-values as confirmatory biological inference yet. The software is
suitable for continued method development, capacity measurement, and negative-control
work. Production sign-off requires, in order:

1. a prespecified null design that resolves global-disjoint feasibility and reference
   exchangeability;
2. capacity and pilot matching on each final VCF-eligible category;
3. finite-pool type-I-error validation including matching and QC selection;
4. a representative end-to-end Farm run within the requested resource envelope.
