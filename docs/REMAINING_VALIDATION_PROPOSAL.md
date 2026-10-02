# Remaining validation before v0.9.0

Status: **frozen, 2026-09-29.** Claude drafted it; Codex reviewed it in two
rounds and approved it with the edits incorporated below. Point E (the V6 test
count) was settled at n = 300 by the user and Codex. No criterion may be changed
after its results have been seen.

Code at `718436f`. Blocking validation runs against the `v0.9.0-rc2` tag (see
"Order").

## Facts checked (2026-09-29)

- **Replicate counts.** `results/sim_dnaging_unpolarised` and
  `results/sim_dnaging_refhap` each contain prepared replicates 1–3 only, and
  `arms_v1` used those three. `results/sim_dnaging` has 10.
- **Depleted-arm capacity.** In
  `results/sim_dnaging/all_mixture_depletion_v1/tests.csv`:
  - M = 4000: every depleted test completed 4–17 sets (3–16 nulls).
  - M = 1000: 100 and 102 tests completed 15–19 sets.
  - M = 250: tests completed 42–80 sets.
  - Tests below 19 nulls are labelled `diagnostic_only_coarse_p`.
- **Ancestral-table integrity.** `results/ancestral-75draw` metadata records the
  store digest but no array digests. Phi-SFS checks only shape and dtype
  (`_checked_table_array` in `normalize_tes/phi_sfs.py`).
- **Stale wording.** README line 464 and
  `results/init_compare/RESULTS_mask_vs_unmasked_qc.md` still describe the
  all-mixture design as unvalidated in simulation.

## Where jobs run

Every job in this plan runs on the `low` partition under account jrigrp (user
decision, 2026-09-29). `low` is preemptible, and a preempted job restarts from
the top. Matcher runs therefore always use a durable `WORK_DIR` with
`--resume`, and a preempted job is resumed with its original seed.

## Blocking status

- **Release blockers:** V1, V2, V4, V5, V6, V9, V10.
- **Non-blocking:** V3, V11.
- **Out of scope for v0.9.0:** V8.

## V1. All-mixture arms where the former design failed (blocking)

**Preparation.** Prepare replicates 4–10 for the unpolarised and refhap
conditions.

- The target is 10 paired replicates; the floor is 6.
- Time one replicate first.
- Choose between 10 and 6 on runtime alone, before any V1 outcome has been seen.
- Whatever the count, use the same underlying msprime replicates in both
  conditions.

**Run.** The four mixture arms (`mixture_{true,inferred}_{reuse,depleted}`) at
M = 250, 1000 and 4000.

**Units.**
- The two conditions are paired, not independent.
- Each cell reports equal-weighted per-replicate rejection rates, with a
  replicate-level (or hierarchical) bootstrap interval.

**Pass criterion.** One-sided non-inferiority, applied per cell:
- mean rejection rate ≤ 0.075;
- one-sided 95% upper bound ≤ 0.10.

**Scope of the pass criterion.**
- It applies to the reuse arms at every M.
- It applies to the depleted arms only for tests with at least 19 nulls.
- The bounds are not relaxed if only 6 replicates are available. If the upper
  bound cannot pass at the available n, V1 is inconclusive and remains blocking.

**Reported, not pass criteria.** The mean Φ_obs − Φ̄⁰ by M, and whether it
grows with M.

These limits guard against serious inflation. They do not establish exact 5%
calibration.

## V2. Depleted-arm capacity report (blocking for interpreting V1)

The cause is established: disjoint sampling exhausts age bins. The report must
give, per replicate and M:

- the distribution of completed set counts;
- the fraction of tests with at least 20 total sets;
- the first limiting age bin and its remaining capacity;
- an explicit denominator for every rejection-rate summary;
- the largest M at which at least 95% of focal tests reach 19 nulls.

## V3. Power (non-blocking; must be reported)

**Designs compared.** All-mixture, snp_pilot-style Bernoulli, and all-Bernoulli,
run on the same focal sets and analysed as paired replicate-level power
differences.

**Effect sizes.** Fixed weighting parameters, set in advance and not tuned to
the observed power. Report the true-DAF displacement each one induces.

**"Materially higher power" requires all three:**
- at least 10 percentage points higher, at two or more of the prespecified
  effect sizes;
- a one-sided replicate-bootstrap lower bound above 5 percentage points;
- the alternative design passes V1.

## V4. Full production matching and depletion (blocking)

> **Amended 2026-10-01 after results were seen; see Amendment A.**

**Run.** `results/targets/in_gene_v0.9_posterior`, 1,001 disjoint sets,
3 restarts, seed 1002.

**Blocking criteria.**
- Maximum reuse is exactly 1, and no set contains a duplicate.
- At least 901 sets pass QC.
- Every complete block of 100 sets has at least 80 QC passes.
- The QC pass rate falls by at most 10 percentage points from the first block
  to the last.
- Reference-free SFS drift. Each set's distance from the pooled control spectrum
  must have |ρ| < 0.1 with replicate index, and an early-versus-late
  standardized mean difference below 0.2.

**Reported diagnostics, not blocking.**
- W1 and error ratio against replicate index (ρ), and the first-to-last block
  shifts in their medians. Some rise is expected by construction under
  depletion.
- Age-residual trends by stratum.
- Signed SFS-bin residual plots, which explain any drift.
- KS statistics.

## V5. Production Phi-SFS on V4 (blocking; integrity only)

> **Amended 2026-10-01 after results were seen; see Amendment A.**

- R ≥ 900.
- Every set used passes QC.
- Maximum reuse is 1, and no set overlaps B0.
- M is exactly equal across A and every set.
- The null mean is finite and the null SD is nonzero.
- The exceedance count and add-one P-value, recomputed, agree with the output.
- Reference sensitivity is reported.
- Every input identity and digest agrees.

The biological P-value is **not** a pass or fail criterion.

## V6. Real-data negative control through the production matcher (blocking)

**Design.**
- Exactly **300 evaluable** held-out SNP focal sets, each of exactly M = 4,000
  sites with no duplicates.
- Each focal set is sampled independently from the restored genome-wide pool.
  Focal sets are **not** forced to be disjoint from one another, because forcing
  it creates finite-population dependence. Overlap between focal sets is
  reported as a diagnostic.
- Each test has its own matched bundle, its own B0, and an independent seed
  fixed in advance.
- Each test's focal sites are excluded from its own control pool.
- The matcher runs 3 restarts.
- R = 99 nulls per test, chosen with `--max-null-replicates 99` (V10b).

- Every test must produce R = 99 valid nulls.
- Jobs run on the `low` partition (account jrigrp), like every run in this
  plan.
  - Preemption and other infrastructure failures are rerun or resumed with the
    same seed.
  - Scientific matching or QC failures are reported, and are never silently
    replaced or excluded.

**Pilot.** Time 5–10 complete bundles first, to confirm wall time and memory.

**Pass criterion.** At most **21 rejections of 300** at α = 0.05. That is the
largest count whose one-sided 95% Clopper–Pearson upper bound (0.0992) is
below 0.10. A calibrated pipeline passes with probability 0.951, and one with a
10% error rate passes with probability 0.046. The observed rate must also be at
most 0.075; that rule is redundant but harmless, because the bound is stricter.

**Reported.** The rank and P-value distribution, rejection rates in blocks of 50
tests, and rejection rate by matcher-order block.

**Scope of the inference.** It is conditional on this genome and this candidate
pool.

**Cost.** The planning estimate is 660–720 job-hours, depending on whether
bundles need 100 or about 110 sets. Each job needs about 40 GB of memory (6
CPUs). These are extrapolations from single runs, to be replaced by the pilot's
measurements.

**Tests that yield fewer than 100 QC-passing sets (decided by the user,
2026-09-29).** Every bundle publishes a fixed number of sets, set in advance at
**110**. Any test with fewer than 100 QC-passing sets counts as a **rejection**,
which is conservative, and is reported. It is not replaced or excluded.

## V7. Exchangeability (review 12, finding 7)

There is no separate study. V6 tests null calibration, and V4's SFS-drift
criteria test deep depletion. Together with V6's rank diagnostics and the
reference sensitivity from V5, they cover the question.

## V8. Nested-category contrasts (out of scope)

`phi_contrast` stays labelled experimental. It is not part of the production
route.

## V9. Release housekeeping (blocking)

- Fix the three tests that assert `0.8.0`.
- The full test suite passes in the declared conda environment.
- `git diff --check` is clean.
- Validation outputs record the exact `v0.9.0-rc2` commit.
- The README and the results note no longer describe the all-mixture design as
  unvalidated in simulation.
- The tagged worktree is clean.
- Every Codex commit has its co-author trailer.
- Done: stale jobs 38983318 and 38983319 were cancelled on 2026-09-29.

## V10. Code required before rc2 (blocking)

### V10a. Ancestral-table integrity

- Record SHA-256 digests for `ancestral_counts.npy` and `present_draw_count.npy`,
  and verify them when Phi-SFS loads the table.
- Test each failure: corruption, swapped arrays, wrong store, incomplete merge,
  a duplicate draw, and a missing draw.
- Confirm that the production table covers exactly the intended posterior draws.

### V10b. `--max-null-replicates N` (optional; absent means use every passing set)

1. Apply QC.
2. Make one deterministic, seeded, SFS-blind permutation of every QC-passing
   replicate ID.
3. The first ID is B0.
4. The next N IDs are the nulls.
5. Fail if fewer than N + 1 sets pass QC.

Metadata records the permutation algorithm, the seed, the selected IDs, the
unused IDs, and the total QC-passing count. With N = 99 the P-value grid is
exactly 0.01, and the nominal 0.05 test has exact attainable size.

## V11. Documented limitations (non-blocking)

None of these validations tests:
- TE genotyping error;
- recurrent insertion or deletion;
- selection in real TE categories;
- ARG failure that differs between TE and SNP sites;
- demography outside the simulated models.

## Settled point E: the V6 test count

Settled at n = 300 by the user and Codex, who checked the exact binomial
figures:

| n tests | largest passing k | P(pass \| error 0.05) | P(pass \| error 0.10) |
|---|---|---|---|
| 200 | 12 | 0.796 | 0.032 |
| 300 | 21 | 0.951 | 0.046 |
| 400 | 29 | 0.981 | 0.036 |

## Amendment A (2026-10-01): V4 and V5 on the first 500 sets

**This breaks the freeze rule above.** The user made the change after seeing the
V4 drift result. It is recorded here so the release does not present the
500-set criteria as prespecified.

**What was seen.** The 1,001-set V4 run (`results/bootstrap_matches/in_gene_rc2`)
passed the five matching criteria: reuse 1, no duplicates, 996 QC passes, at
least 95 passes in every block of 100, and a 5-point fall. It **failed** SFS
drift: rho = 0.164 and early-versus-late SMD = 0.476
(`results/v4/in_gene_rc2`). Mean distance to the pooled spectrum is flat through
replicate 799 and rises in blocks 800-899 and 900-999. Median W1 against the
bootstrap target follows the same pattern, at 50-54 through block 6 and then
85, 185 and 320, and all five QC failures fall in block 900-999. Late sets
appear to be matched from a depleted pool. That has not been checked per
stratum.

The definitions `tools/v4_depletion_report.py` uses were fixed before the drift
was computed: Phi distance to the pooled spectrum, Spearman rho, and SMD between
the first and last quarters. Under the alternatives checked afterwards, the first
500 sets pass with Phi, but with L1 distance the quarter SMD is 0.279 and fails.
The first 400 sets pass under both.

**Change.**
- V4 and V5 use replicates 0-499 of the same run, written by
  `tools/trim_match_bundle.py` to `results/bootstrap_matches/in_gene_rc2_first500`.
  Seeds depend only on the global seed, the target digest, the replicate and the
  restart, and depletion is sequential. The first 500 sets are therefore the
  sets a 500-set run would produce. That follows from the code; no separate
  500-set run was made.
- V4: at least 451 QC passes (901/1,001 scaled to 500). The block, QC-fall and
  drift criteria are unchanged and apply to the five complete blocks.
- V5: R >= 450 replaces R >= 900. Reference sensitivity uses N = 10 alternative
  references (`REFERENCE_SENSITIVITY=10`). The proposal had not fixed N.
- The smallest attainable P rises from 1/996 (R = 995) to 1/500 (R = 499).

The 1,001-set results stay on record as a V4 failure.

## Amendment B (2026-10-01): which commit counts as rc2 (V9)

**This changes a V9 item after results were seen.** The user made the change on
2026-10-01.

**What was seen.** No validation output records the `v0.9.0-rc2` commit itself.
Each output records `git describe` at the time it was written, and these runs span
several later commits:

| output | recorded |
|---|---|
| V4 bundle, V5 (1,001 sets) | `v0.9.0-rc2-4-g9ec39f5` |
| V5 (first 500 sets), V1 arms | `v0.9.0-rc2-4-g9ec39f5-dirty` |
| V6 (300 tests) | 277 at `-4-g9ec39f5`, 2 at `-3-gf780721`, 6 at `-2-g5889409`, 15 dirty (below) |

The commits after rc2 change only `slurm/` and `tools/`:
`git diff --stat v0.9.0-rc2 9ec39f5 -- normalize_tes/` is empty.

**Change.** The V9 item "Validation outputs record the exact `v0.9.0-rc2`
commit" becomes: every validation output records `v0.9.0-rc2` or a later commit
whose `normalize_tes/` is identical to rc2's, shown by an empty
`git diff v0.9.0-rc2 <commit> -- normalize_tes/`.

**Dirty trees.** A `-dirty` label means uncommitted edits were present, and the
outputs do not record what they were. The matcher keeps a hash of its loaded
source modules in `work/identity.json`, but that directory is deleted when a run
completes, so it is not available afterwards.

- V5 (first 500 sets) and the V1 arms ran with only
  `docs/REMAINING_VALIDATION_PROPOSAL.md` and `slurm/run_phi_sfs.sbatch`
  modified. Those edits were checked at submission, and neither file is in
  `normalize_tes/`.
- V6 tests 3 and 5 (`-2-g5889409-dirty`) and tests 162, 163, 168, 175, 181, 189,
  216, 244, 245, 274, 279, 291 and 298 (`-4-g9ec39f5-dirty`) have no record of
  what was modified. They do not yet meet the amended item.

**Resolution (2026-10-02).** The 15 tests were rerun on a clean checkout with the
same seeds. Ten reproduced exactly and one reproduced its P. The other four
differed because their age targets were rebuilt with a different thread count or
CPU type, which changes the bootstrap distances in the last digits and, through
the target digest, the choice of $B_0$; a same-node test confirmed this. The
uncommitted edits did not change the analysis. V6 is judged on the original runs
(11/300); the reruns give 10/300. Details are in `docs/VALIDATION.md`.

**Default values (2026-10-02).** After validation, the package defaults changed
to match production: `bootstrap_target_matcher --replicates` from 1001 to 500 and
`phi_sfs --min-null-replicates` from 900 to 450. Every validation run passed these
values explicitly, so no result depends on the defaults. The release commit's
`normalize_tes/` therefore differs from rc2's in these two default values (and
their help text) and in the version string (`0.9.0-rc2` to `0.9.0`) only;
`git diff v0.9.0-rc2 <release> -- normalize_tes/` shows exactly that.

## Order

1. Commit this frozen file.
2. Implement V10a and V10b, and make the V9 code and wording fixes.
3. Pass the full test suite. Tag that exact commit `v0.9.0-rc2`.
4. Against rc2, run V2 and the V1 replicate timing and preparation, then V1 and
   V4 in parallel.
5. Run V5, then the V6 pilot, then V6.
6. Run V3 when convenient, and complete V11.
7. Tag `v0.9.0` when V1, V2, V4, V5, V6, V9 and V10 all pass.
