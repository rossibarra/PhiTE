# Code Review — Round 12

Date: 2026-09-27
Reviewed revision: `1c95b81` on `bootstrap-target-hpc-validation` (v0.9.0-rc1,
`8ebd681`, plus `--init-mode` from `5268ff4` and the README wording fixes of
`1c95b81`). Every `file:line` citation below refers to `1c95b81`.
Scope: the statistical design of the Wasserstein Phi-SFS calibration as documented in
`README.md` sections "Polarity and the asymmetric null", "Null calibration, Z-scores,
and P-values" and "SNP type-I pilot"; the category contrast in
`normalize_tes/phi_contrast.py`; and the matching-QC evidence the design depends on.

The worktree contains an uncommitted extension to `tools/validate_phi_calibration.py`;
it is not part of `1c95b81` and is not reviewed here.

This review supersedes the first round-12 draft (`7491f71`) and the separate
Claude + Codex consensus note. Codex (`codex-cli 0.155.0-alpha.16.3`, read-only
sandbox) reviewed the draft independently against the code in two rounds. Claude
checked every factual correction Codex made against the source before accepting it;
the corrections that changed a finding are marked in that finding. The prompts and
Codex's replies are in the session scratchpad (`codex/round{1,2}_{prompt,reply}.md`),
not in the repository.

**Evidence status.** Findings 1–5, 7 and 8 are reasoning about the method; none of
them has been tested, and statements in them about what a failure *would* do are
predictions. Finding 6 is measured. The section "Proposed validation study" is the
agreed design for testing findings 1–5, 7 and 8.

## Summary

The Bernoulli-q asymmetric null improves on the legacy mixture-versus-mixture null:
it reproduces the hard-versus-uncertain polarity architecture of the observed TE
comparison. Whether it also reproduces that comparison's null distribution is not
established. A TE's hard orientation is treated as biological truth, while a null
set's hard orientation is a draw from the ARG posterior, so the two constructions
agree only if the ARG's polarity, age and frequency inference are jointly calibrated
and the matched sets are exchangeable with the focal comparison. The TE polarity
filter and agreeing-draw TE ages are applied to A only. The SNP type-I pilot builds
its focal sets and nulls with the same Bernoulli-q treatment, so it cannot detect a
failure specific to the TE path.

Two further inferential gaps remain: drawing B₀ at random does not by itself make
the A-versus-B₀ comparison exchangeable with the null comparisons (carried over from
Round 11), and production conditions on one polarity-imputation seed. The category
contrast pairs nulls by random permutation, which is defensible for independent
categories but ignores covariance between nested ones. Z is a test statistic that
grows with M and is documented as an effect size.

Separately, and measured: no matched set currently passes matching QC. In one
controlled 10-replicate experiment the starting set did not change that.

## Findings

### 1. High — the asymmetric null is valid only if ARG polarity, age and frequency inference are jointly calibrated and matching is exchangeable

The observed statistic is
$\Phi_{\mathrm{obs}}=\Phi(A_{\mathrm{TE,hard}},B_{0,\mathrm{mixture}})$. A retained
TE is oriented with insertion presence derived, which the README treats as known, so
under the null the spectrum of A is a draw around the true spectrum of an age-matched
neutral set. A null distance is
$\Phi_i^0=\Phi(B_{i,\mathrm{Bernoulli}(q)},B_{0,\mathrm{mixture}})$, and by the
README's identity $E[h]=q\,h(k,n)+(1-q)\,h(n-k,n)$ the expected spectrum of $B_i$ is
the posterior mixture, not the truth.

A calibrated marginal q is necessary for the two to agree in expectation, but it is
not sufficient. The same ARG supplies q, the allele ages, and therefore the matching
variables, and errors in polarity, age and frequency can be correlated (for example,
a mis-polarized high-frequency allele also receives a mis-estimated age). A null set
matched on inferred ages then carries those correlated errors, while A carries hard
TE orientation and (finding 3) filtered ages.

If these differ, $\Phi_{\mathrm{obs}}$ plausibly contains a systematic term (true
spectrum versus posterior mixture, conditional on matching) that no $\Phi_i^0$
contains. The first draft said a neutral TE category would then show a positive Z and
excess rejections. That is softened here: a mismatch probably adds distance, but its
sign after normalization and matching is not guaranteed. Untested.

Recommendation: state the joint-calibration and exchangeability assumptions in the
README. Test them where truth is known: simulation arms 1, 2 and 4 of the proposed
validation study compare a true-hard focal set with the production null.

### 2. High — the SNP type-I pilot cannot validate the TE path

For `-A SNP`, focal A receives the same Bernoulli-q hard treatment as every null set
(`normalize_tes/phi_sfs.py:1123-1129`), and the pilot tool splits one matched pool
into held-out A sets, one reference and the nulls
(`tools/empirical_asymmetric_null.py:99-117`). The focal sets and nulls therefore
share the polarity and SFS construction. The pilot (7 of 100 rejected at 0.05;
`results/phi_sfs/snp_type1_asymmetric_100/summary.json`) cannot test TE hard
polarity, the TE polarity filter, or agreeing-draw TE ages.

Codex correction, accepted: the first draft said the pilot "should pass by design".
It does not. It can still detect SNP-path SFS bugs and some matching or
non-exchangeability failures. What it cannot see is anything that affects only a TE
focal set.

The pilot has two further limitations, both already reported in the README and in
`summary.json`: all 100 tests share one reference (replicate 43) and one null vector
of 344 sets, so the P-values are dependent; and none of the 445 completed sets pass
the matching-error-ratio criterion (median ratio 2.86).

Recommendation: describe the pilot as validating the SNP path only (done in the
README, below). The TE-path negative control is simulation arm 2 versus arm 3.

### 3. High — the TE polarity filter and agreeing-draw ages apply to A only

A TE is retained only when its flipped fraction among usable ARG draws is at most
`--max-flipped-fraction` (0.5 in production), that is, when at least half of usable
draws support presence as derived (`normalize_tes/te_age_target.py:424-436`;
44 of 4,067 in-gene TEs discarded, 4,023 kept, per
`results/targets/in_gene_v0.8/metadata.json`). SNPs keep every q. If the ARG
mis-polarizes TEs more often at high derived frequency, the filter would remove
high-frequency TEs selectively and shift A's spectrum towards rare alleles, and
nothing in B mimics that. The README describes the filter as removing homoplasy
(deletion, recurrent movement); it may also remove genuine single insertions that
the ARG mis-polarizes. Untested.

A masked TE's age CDF is built from agreeing draws only, while SNP ages use all
draws, so A and B are age-matched on differently conditioned quantities. Matching is
already the binding constraint (finding 6).

Codex corrections, both accepted after checking the code:

- (a) Agreeing-draw ages have fallbacks. At the mask level, a TE with no agreeing
  draw keeps all present draws (`normalize_tes/te_age_target.py:438-440`); at the
  interval level, a TE whose agreeing draws supply no usable age interval keeps all
  of its intervals (`normalize_tes/te_age_target.py:271-279`). For in-gene v0.8 the
  metadata records 0 TEs with no agreeing draw, so the fallbacks did not change that
  target; the asymmetry applies to every other retained TE.
- (b) The first draft recommended a real-data SNP `q >= 0.5` filter on B as a
  matched analogue. That is **withdrawn**. SNP q is P(ALT is derived), so such a
  filter depends on arbitrary REF/ALT labels, not on orientation reliability. A
  symmetric filter needs known truth, which only simulation supplies.

Recommendation: in simulation, compare true DAF and true age of retained and
discarded A sites (arm 4 and the secondary filter arms). On the real data, report
the observed derived-frequency distribution of discarded versus retained TEs as a
descriptive check only.

### 4. Medium — Z is a test statistic that grows with M, not an effect size

$s_0$ is expected to shrink roughly as $1/\sqrt{M}$, so for a fixed true spectral
difference $Z_A$ grows with category size. Plotting one Z per category will then rank
large categories as more different at equal effect. At `1c95b81` the README called Z
both a "null-standardized effect size" and a description of "the magnitude of that
departure".

The first draft recommended the finite-sample correction
$\Phi_\infty=\sqrt{\Phi_{\mathrm{obs}}^2-\Phi_{\mathrm{floor}}^2}$ (validated to about
3% for $n_{\mathrm{TE}}\ge250$ under the earlier design) as the effect size. That is
**not** recommended until it has been revalidated under the asymmetric null and the
current matching design.

Recommendation: call Z a standardized test statistic. Report
$\Phi_{\mathrm{obs}}-\mu_0$ (DAF units), signed CDF and bin residuals, and M as the
magnitudes. Do not compare Z across categories of different M as an effect size.
Measure the M dependence in simulation (M = 250, 1,000, 4,000).

### 5. Medium — category contrasts ignore covariance between nested or overlapping categories

**Replaced.** The first draft (and Round 11, finding 4) said the contrast paired
nulls by index. That is wrong for the code. `contrast()` in
`normalize_tes/phi_contrast.py:95-141` sets $R=\min(R_1,R_2)$, draws a random
permutation of category 1's null Z-scores truncated to R, pairs it with category 2's
nulls (a random subset of R when $R_2>R$), and computes a two-sided add-one P-value
(`phi_contrast.py:142-143`). The pairing is seeded from `--seed` and the sorted
category labels (`phi_contrast.py:73-92`); repeat 0 is the reported result, and
`--pairing-repeats` (default 100) gives a P-value sensitivity range.

For two independent categories this is a defensible sample from the product null,
though it uses only R of the $R_1R_2$ possible differences. The open problem is
nested or overlapping categories (in-gene TEs versus all TEs): they share focal sites
and so have correlated Z, and independent pairing discards that covariance. I expect
the contrast to be conservative for positively correlated categories, but that is
unverified. At `1c95b81` the README formula $\Delta_i^0=Z_{1i}^0-Z_{2i}^0$
misdescribed the code.

Recommendation: fix the README description (done, below) and warn that nested
categories are not handled. Use nested pseudo-categories in simulation to measure
the calibration of the current contrast before designing a covariance-aware one.

### 6. High — no matched set passes QC (measured)

Production requires at least 900 QC-passing null sets (`--min-null-replicates`,
`normalize_tes/phi_sfs.py:893-897`). A set passes QC only if its matching-error
ratio is below `--qc-max-ratio` (0.5) **and** its best W1 is at most
`--qc-max-absolute-fraction` (0.34) × the target's acceptance threshold
(`normalize_tes/bootstrap_target_matcher.py:940-947`). For in-gene v0.8 the
acceptance threshold is 1,356.06 generations (`logs/bootstrap-match-38964576.out`),
so the absolute cap is about 461 generations.

Codex correction, accepted: the first draft compared best W1 with 1,356 generations,
which is the acceptance threshold, not the QC cap. With the correct cap the gap is
larger: the interrupted in-gene v0.8 run (job 38964576) completed 445 replicates whose
best W1 ranged 2,943–4,805 generations (median about 3,755), so every one fails both
criteria.

To test whether the median-age starting set was the cause, 10 replicates × 3
restarts were run with seed 1002 against the same target, varying only
`--init-mode` (jobs 39016374–39016376, outputs in `results/init_compare/`). W1
columns summarize all 30 restarts per arm (`restarts.csv`); QC is per replicate
(`replicates.csv`):

| init mode | initial W1, median | best W1, median (range) | QC-passing replicates |
|---|---|---|---|
| median (current) | 74,341 | 3,834 (3,414–4,527) | 0/10 |
| mass | 18,480 | 3,859 (3,408–4,533) | 0/10 |
| random | 877,653 | 3,844 (3,420–4,581) | 0/10 |

The starting distance differs almost 50-fold between arms, but the final distances
agree closely, and each replicate's selected matching-error ratio differs by at most
2.8% across the three arms (median-arm ratios 1.71–5.27). In this experiment the
search reached the same plateau from every start.

That is a result for this 10-replicate, 3-restart experiment on one target; it is
not a general result about initialization. The plateau may reflect optimizer
settings (proposal count, epochs, the material-improvement stopping rule), the
support of the candidate pool, or both. The per-stratum capacity preflight passing
(minimum 1,212 sets) does not guarantee that a set can match the full CDF shape.

Recommendation: before any further production run, find out which part of the age
CDF the best sets miss (compare `restart_best_cdfs.npy` with
`bootstrap_target_cdfs.npy` by age) and whether the SNP pool has enough mass there
with the right shape, and separately test longer or differently tuned searches.
Until QC passes, the README's production route cannot run.

### 7. High — exchangeability is still unproven (carried over from Round 11, finding 1)

Drawing B₀ uniformly from the QC-passing sets with a prespecified seed
(`normalize_tes/phi_sfs.py:854-868`) makes the reference role uniform and SFS-blind
within the accepted sets. It does not make sequentially depleted sets identically
distributed, and it does not make the A-versus-B₀ comparison exchangeable with the
Bᵢ-versus-B₀ comparisons, which is what the add-one P-value needs. A is a fixed
observed set with its own construction (findings 1 and 3); the nulls are matched
SNP sets.

The README sentence added in `1c95b81` ("drawing $B_0$ at random (below) is what
makes it exchangeable with the nulls") overclaimed and is reworded below. The code
comment at `normalize_tes/phi_sfs.py:858` ("making B0 exchangeable with the rest")
makes the same claim about the reference role only; it is accurate in that narrow
sense but should not be read as the P-value's exchangeability condition.

Recommendation: validate exchangeability with the depletion-diagnostic grid of the
simulation study, which uses fully disjoint nulls through the production matcher.

### 8. Medium — polarity-imputation uncertainty is not reported

Production hard-orients every null set with one coordinate-keyed Bernoulli-q draw
per SNP from a single seed (`--polarity-imputation-seed`, default 2001,
`normalize_tes/phi_sfs.py:1607-1609`; used at `phi_sfs.py:1113-1116`). Z and P are
therefore conditional on one imputation. Whether that imputation drives them is
unmeasured.

Recommendation: report Z and P across several prespecified seeds (a secondary
simulation arm first, then production).

## Proposed validation study

This is the design Claude and Codex agreed. None of it has been run.

**Source.** `/quobyte/jrigrp/beil/logan_collab/dnAging/msprime_variable_ne_error/`
(read-only; owned by another user). `simulations/` holds 100 msprime replicates;
`singer/` holds SINGER output for 10 of them (`replicate_001`–`replicate_010`), each
with 100 posterior tree draws. The study uses those 10. Each has its own random
piecewise-constant Ne history, 26 modern and 10 ancient haploids, and infinite-sites
mutations. The consensus pass recorded 25,525–80,703 modern-polymorphic sites per
replicate (363,113 in total), ascertained to derived counts 1–25, exact mutation times
and derived counts (`mutation_truth.tsv.gz`), and the true modern ARG
(`known_modern_arg.trees`). The VCF writes the ancestral allele as REF, so true
polarity is known. Codex verified the site counts and polarity encoding; Claude
verified the sample sizes, draw count and file layout. Treat replicates as
independent units and combine estimates across them afterwards, not by pooling
sites.

**Per-site quantities**, computed from all 100 draws: true age and true derived
count; posterior support for the true orientation; an age CDF from all draws; and an
age CDF from truth-agreeing draws only.

**Allele-label conventions** (for q-reliability):

1. native (REF = ancestral);
2. independent 50% REF/ALT swaps;
3. reference-genome analogue: REF = the allele carried by one designated modern
   haplotype, repeated for several haplotypes. This mimics B73-referenced maize SNPs,
   where REF correlates with ancestral state and frequency.

**Grid.** No site appears twice within a set, and A and B₀ sites are excluded from
every null.

| Purpose | M | R nulls | Null reuse |
|---|---:|---:|---|
| Primary polarity / type-I study | 250, 1,000, 4,000 | 199 | reused across nulls |
| Depletion diagnostic | 250 | 79 | fully disjoint |
| | 1,000 | 19 | fully disjoint |
| | 4,000 | 4 | fully disjoint; matching diagnostic only |

With reuse at M = 4,000, two nulls drawn from about 17k remaining sites share
roughly a quarter of their sites, so the nulls are correlated. R = 19 is the smallest
R that can reject at α = 0.05 with the add-one P-value. If a filter leaves too few
sites, reduce M, not R.

**Priority arms.** Arms 1–3 match on true age, to isolate polarity; arm 4 uses
inferred ages.

1. **Oracle:** true-hard A, B₀ and Bᵢ. Gives the attainable type-I rate.
2. **Production asymmetry, unfiltered:** true-hard A, mixture B₀,
   Bernoulli-q-hard Bᵢ. Tests finding 1.
3. **SNP-pilot construction:** Bernoulli-q-hard A and Bᵢ, mixture B₀, on the same
   partitions as arm 2. The arm 2 − arm 3 difference is what the pilot cannot see
   (finding 2).
4. **Full TE-like:** keep A sites whose support for the true orientation is ≥ 0.5,
   give A agreeing-draw ages, and leave B unfiltered with all-draw ages, matched on
   inferred ages. Compared with arm 2, this measures the combined filter and
   age-conditioning effect (findings 3 and 1).

**Secondary arms:** selection-only versus agreeing-age-only; a truth-aware symmetric
filter on both sides; the legacy mixture-versus-mixture null; several Bernoulli
seeds (finding 8); nested pseudo-categories for the contrast (finding 5); and
inferred-age matching with the production matcher and init modes (finding 6,
algorithmic aspects only), which with the depletion grid also addresses finding 7.

**Outcomes per arm:** $\Phi_{\mathrm{obs}}-\overline{\Phi^0}$, Z and P; the type-I
fraction and P-value uniformity across replicates; dependence on M (finding 4);
q-reliability and hard-orientation error by true DAF and true age; true DAF and age
of retained versus discarded A sites; inferred-age and true-age W1.

**Cannot test:** TE-specific homoplasy (deletion, recurrent insertion), TE genotyping
error, selection, real demographic misspecification, and production-scale matching
feasibility on the 23-million-site maize pool.

## Documentation corrections

Applied to `README.md` in `1c95b81` (backup `README.md.pre-round12.bak`):

- The distance symbols $D_{\mathrm{obs}}$ and $D_i^0$ are now $\Phi_{\mathrm{obs}}$
  and $\Phi_i^0$, removing the clash with the age distance $D(\cdot,\cdot)$ in the
  B/E/O/R definitions.
- Reference sensitivity was described as using "its own first-$`R`$-QC-passing" nulls;
  each alternative reference now reads as using every other QC-passing set as its
  nulls (`normalize_tes/phi_sfs.py:906-916`).
- "Replicate 0 is not used by default" read as though replicate 0 were excluded; it
  can be drawn as $B_0$ or used as a null like any other QC-passing set.
- "Generates all $R+1$ sets identically" was replaced by a statement that disjoint
  matching makes the sets not identically distributed. The sentence added with it,
  that drawing $B_0$ at random makes it exchangeable with the nulls, overclaimed and
  is corrected below (finding 7).
- `m` in the Wasserstein sum is now defined (the projection size, 20).
- "A masked target (step 5)" is now step 7; "using step 3" for preliminary targets is
  now step 5.
- The methods list named Round 9 as the latest review; it now names Round 12.

Proposed with this review (draft README):

- "Null calibration, Z-scores, and P-values", step 1: random $B_0$ is described as
  making the reference choice uniform and SFS-blind among accepted sets; the README
  now says this does not make depleted sets identically distributed or show that the
  A-versus-$`B_0`$ comparison is exchangeable with the null comparisons, and that this
  is an open validation item (finding 7).
- Step 4 reports $Z_A$ as a "standardized test statistic", not a "null-standardized
  effect size". The plot guidance is kept, with a warning that Z is not comparable as
  an effect size across categories of different M. The closing paragraph says Z grows
  with M for a fixed spectral difference and names $\Phi_{\mathrm{obs}}-\mu_0$, signed
  CDF and bin residuals, and M as the magnitudes (finding 4).
- "SNP type-I pilot": added that the pilot shares the Bernoulli-q construction with
  the nulls, so it tests the SNP path only and cannot validate TE hard polarity, the
  TE polarity filter, or agreeing-draw ages (finding 2).
- The between-category contrast now describes `phi_contrast.py`: random-permutation
  pairing, truncation to $\min(R_1,R_2)$, two-sided add-one P-value, seeded pairing
  with `--pairing-repeats` sensitivity; and warns that it ignores covariance between
  nested or overlapping categories (finding 5).
- "Polarity and the asymmetric null": added that production results condition on a
  single `--polarity-imputation-seed` and that sensitivity across seeds is not yet
  reported (finding 8).
- The methods-list entry for this review mentions the Codex review and the proposed
  validation study.

Also corrected in the same pass: the README sentence "$`B_0`$ is drawn before any SFS is
examined, so it is exchangeable with the other QC-passing sets" now says only that the
reference choice is SFS-blind (finding 7); the figure alt text "standardized Phi-SFS
effects" now reads "test statistics" (finding 4); and two `normalize_tes/phi_sfs.py`
comments were aligned with the code, one that claimed random $B_0$ makes it
exchangeable and one that described reference-sensitivity nulls as "independently
selected" when every other QC-passing set is used (`phi_sfs.py:906-909`).

## Verification performed

- Read `README.md` at `1c95b81` and the `1c95b81` diff; read the B₀, null and
  reference-sensitivity selection (`normalize_tes/phi_sfs.py:850-920`), the SNP
  focal-set construction (`phi_sfs.py:1110-1135`) and the seed option
  (`phi_sfs.py:1607-1609`).
- Read `contrast()`, `paired_contrast()` and the repeat handling in
  `normalize_tes/phi_contrast.py`.
- Read the QC rule and defaults in `normalize_tes/bootstrap_target_matcher.py`
  (lines 54–56 and 940–947), the TE polarity filter and fallbacks in
  `normalize_tes/te_age_target.py` (lines 271–279, 424–440), and the pilot split in
  `tools/empirical_asymmetric_null.py:99-117`.
- Read `results/phi_sfs/snp_type1_asymmetric_100/summary.json` and
  `results/targets/in_gene_v0.8/metadata.json`.
- From `logs/bootstrap-match-38964576.out`: acceptance threshold 1,356.06; per-replicate
  best W1 2,943–4,805 (the log also contains one partial replicate).
- Recomputed the init-mode table from `results/init_compare/*/restarts.csv` and
  `replicates.csv` with shell tools, and the per-replicate ratio spread across arms
  (at most 2.8%); read `metadata.json` for mode, seed 1002, 3 restarts and QC settings.
- Listed the simulation source: 100 simulated replicates, SINGER output for 10, 100
  posterior draws in `replicate_001`. Site counts were verified by Codex in the
  consensus pass and not re-checked here.
- Codex (`codex-cli 0.155.0-alpha.16.3`) reviewed in two rounds; Claude checked each
  of its factual corrections (findings 2, 3a, 3b, 5, 6) in the code.
- No tests, Python or Slurm jobs were run for this review. Findings 1–5, 7 and 8 were
  not tested.

## Release recommendation

Do not tag v0.9.0 as a production release. Finding 6 blocks the production route
outright: no matched set passes QC. Findings 1, 2, 3 and 7 determine whether a TE
P-value can be interpreted at all and should be tested with the proposed simulation
study (priority arms 1–4 and the depletion grid) before any TE result is reported.
Findings 4, 5 and 8 affect how results are reported and compared; the README changes
above address their documentation, and their measurement belongs in the same study.
