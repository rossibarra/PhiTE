# Code Review — Round 12

Date: 2026-09-27
Reviewed revision: `8ebd681` (v0.9.0-rc1 on `bootstrap-target-hpc-validation`)
Scope: the statistical design of the Wasserstein Phi-SFS calibration as documented in
`README.md` sections "Polarity and the asymmetric null", "Null calibration, Z-scores,
and P-values" and "SNP type-I pilot", and the matching-QC evidence it depends on.
This is a design review. Findings 1–5 are reasoning about the method and have not
been tested unless a measurement is cited. Finding 6 is measured.

## Summary

The Bernoulli-q asymmetric null improves on the legacy mixture-versus-mixture
null: it reproduces the variance difference between a hard-oriented focal spectrum
and a posterior-mixture reference. It does not yet reproduce the observed TE
comparison in expectation. A TE's hard orientation is treated as biological truth,
while a null set's hard orientation is a draw from the ARG posterior, so the two
agree only if the posterior q is calibrated. The SNP type-I pilot cannot detect a
failure of that assumption, because on the `-A SNP` path the focal set and the nulls
are built the same way. The TE polarity filter is also applied to A and not to B.
Separately, no matched set currently passes matching QC, and a controlled
comparison of starting sets shows this is not a problem with how the search starts.

## Findings

### 1. High — the asymmetric null matches the observed comparison only if q is calibrated

The observed statistic is
$\Phi_{\mathrm{obs}}=\Phi(A_{\mathrm{TE,hard}},B_{0,\mathrm{mixture}})$. A retained
TE is oriented with insertion presence derived, which the README treats as known, so
$E[\text{spectrum of }A]$ is the true spectrum of an age-matched neutral set under
the null. A null distance is
$\Phi_i^0=\Phi(B_{i,\mathrm{Bernoulli}(q)},B_{0,\mathrm{mixture}})$, and by the
README's own identity $E[h]=q\,h(k,n)+(1-q)\,h(n-k,n)$ the expected spectrum of
$B_i$ is the posterior mixture, not the truth.

If the ARG's q is biased (for example, if it pulls high-frequency derived alleles
towards ancestral), the two differ in expectation. $\Phi_{\mathrm{obs}}$ then contains a
term (true spectrum versus posterior mixture) that no $\Phi_i^0$ contains, and a neutral
TE category would be expected to show a positive Z and too many rejections. The null
matches the variance structure; it matches the mean only when q is calibrated.

Recommendation: state the calibration assumption explicitly in the README. Test it
on SNPs whose orientation is known independently of the ARG (outgroup-polarized
sites treated as "hard truth", or forward simulations), comparing their hard spectrum
with the Bernoulli-q hard spectrum of the same sites.

### 2. High — the SNP type-I pilot cannot validate the TE path

For `-A SNP`, focal A receives the same Bernoulli-q hard treatment as every null set,
so the focal set and the nulls are the same construction. Apart from matching error,
the pilot should pass by design (7 of 100 rejected at 0.05; `summary.json` for
`results/phi_sfs/snp_type1_asymmetric_100`). It tests the SNP negative-control path;
it says nothing about finding 1 or finding 3, which affect only a TE focal set.

The pilot has two further limitations, and the README already reports both: all 100
tests share one reference (replicate 43) and one null vector of 344 sets, so the
P-values are dependent, and none of the 445 completed sets pass matching QC (median
error ratio 2.86).

Recommendation: describe the pilot as validating the SNP path only. Add a TE-path
negative control in which the focal set's hard orientation comes from a source other
than the ARG posterior (see finding 1).

### 3. High — the TE polarity filter is applied to one side only

A TE is retained only when at least 50% of usable ARG draws support presence as
derived. SNPs keep every q. If the ARG mis-polarizes TEs more often at high derived
frequency, the filter removes high-frequency TEs selectively and shifts A's
spectrum towards rare alleles, and nothing in B mimics that. The README describes the
filter as removing homoplasy (deletion, recurrent movement); it also removes genuine
single insertions that the ARG mis-polarizes.

A related asymmetry: a masked TE's age CDF is built from agreeing draws only, while
SNP ages use all draws. That changes what is being age-matched, and matching is
already the binding constraint (finding 6).

Recommendation: compare the observed derived-frequency distribution of discarded and
retained TEs (44 of 4,067 discarded for in-gene v0.8). If discarded TEs are enriched
at high frequency, either apply a matched filter to B (drop SNPs with $q<0.5$) as a
sensitivity analysis, or model the filter's effect on the null.

### 4. Medium — Z is presented as an effect size but grows with M

$s_0$ is expected to shrink roughly as $1/\sqrt{M}$, so for a fixed true spectral
difference $Z_A$ grows with category size. Plotting one Z per category (README
"one equal-size point per focal category at $Z_A$") will rank large categories as
more different even at equal effect. The README calls Z both a "null-standardized
effect size" and a description of "the magnitude of that departure".

The finite-sample correction developed earlier,
$\Phi_{\infty}=\sqrt{\Phi_{\mathrm{obs}}^2-\Phi_{\mathrm{floor}}^2}$ (validated to
about 3% for $n_{\mathrm{TE}}\ge250$), estimates the effect size and is no longer
mentioned.

Recommendation: call Z a standardized test statistic, report $\Phi_\infty$ (or
$\Phi_{\mathrm{obs}}-\mu_0$ in DAF units) as the effect size, and compare effect sizes
across categories rather than Z.

### 5. Medium — the between-category contrast pairs nulls by index (carried forward from Round 11, finding 4)

$\Delta_i^0=Z_{1i}^0-Z_{2i}^0$ assumes null $i$ in category 1 corresponds to null
$i$ in category 2. The categories have different matched sets, $R$ now differs by
category, and overlapping categories (in-gene versus all TEs) share focal sites, so
the pairing is arbitrary and may not be defined. This is still open.

Recommendation: remove the contrast from the README until it has a design, or
replace it with an unpaired comparison of effect sizes with bootstrap intervals.

### 6. High — no matched set passes QC, and the starting set does not change that (measured)

Production requires at least 900 QC-passing sets (`--min-null-replicates`). The
interrupted in-gene v0.8 run (job 38964576) completed 445 replicates with none below
the acceptance threshold of 1,356 generations: best W1 per replicate ranged
2,943–4,805, median 3,755.

To test whether the median-age starting set was the cause, 10 replicates × 3
restarts were run with seed 1002 against the same target, varying only
`--init-mode` (jobs 39016374–39016376, outputs in `results/init_compare/`):

| init mode | initial W1, median | best W1, median (range) | QC-passing replicates |
|---|---|---|---|
| median (current) | 74,341 | 3,834 (3,414–4,527) | 0/10 |
| mass | 18,480 | 3,859 (3,408–4,533) | 0/10 |
| random | 877,653 | 3,844 (3,420–4,581) | 0/10 |

The starting distance differs almost 50-fold between arms, but the final distance
and each replicate's matching-error ratio agree to within about 3% (median-arm
ratios 1.71–5.27). The search therefore reaches the same plateau from any start.
This suggests (not yet tested) that the limit lies in the candidate pool or in the
swap search itself, not in initialization. The per-stratum mass preflight passing
(minimum 1,212 sets) does not guarantee that a set can match the full CDF shape.

Recommendation: before any further production run, find out which part of the age
CDF the best sets miss (compare `restart_best_cdfs.npy` with
`bootstrap_target_cdfs.npy` by age), and whether the SNP pool contains enough mass
there with the right shape. Until then, the README's production route cannot run.

## Documentation corrections

These were found in the same pass and are wording-only. All were applied to
`README.md` with this review (backup `README.md.pre-round12.bak`); the distance
symbols $D_{\mathrm{obs}}$ and $D_i^0$ are now $\Phi_{\mathrm{obs}}$ and $\Phi_i^0$.

- Reference sensitivity is described as using "its own first-$R$-QC-passing" nulls;
  the code uses every other QC-passing set (`phi_sfs.py`, reference-sensitivity block).
- "Replicate 0 is not used by default" reads as though replicate 0 is excluded; in
  the code it can be drawn as $B_0$ or used as a null like any other set.
- The R+1 sets are described as generated "identically". Under disjoint matching they
  are not identically distributed; the random draw of $B_0$ makes $B_0$ exchangeable
  with the nulls, which is what the P-value requires.
- "A masked target (step 5)" should be step 7; "using step 3" for preliminary targets
  should be step 5.
- The methods list names Round 9 as the latest review.
- `m` in the Wasserstein sum is not defined in the README (it is the projection size,
  20), and `D` for spectral distances collides with the age distance $D(\cdot,\cdot)$
  in the B/E/O/R definitions.

## Verification performed

- Read `README.md` at `8ebd681` and the B₀/null selection code in
  `normalize_tes/phi_sfs.py`.
- Read `results/phi_sfs/snp_type1_asymmetric_100/summary.json`.
- Summarized `logs/bootstrap-match-38964576.out` and `results/init_compare/*/restarts.csv`
  and `replicates.csv`.
- Findings 1–5 were not tested.

## Release recommendation

Do not tag v0.9.0 as a production release. Finding 6 blocks the production route
outright. Findings 1–3 determine whether a TE P-value can be interpreted at all and
should be tested before any TE result is reported.
