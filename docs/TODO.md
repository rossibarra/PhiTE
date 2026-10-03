# To do

Open work after v0.9.0, most important first.

## V3: power of the all-mixture design (tabled 2026-10-02)

Non-blocking item V3 in
[REMAINING_VALIDATION_PROPOSAL.md](REMAINING_VALIDATION_PROPOSAL.md). Compare
all-mixture, snp_pilot-style Bernoulli and all-Bernoulli designs on the same focal
sets, as paired replicate-level power differences.

Before running:

- **Fix the effect sizes first**, before any power result is seen: weighting
  parameters for focal-site selection (for example weight ∝ exp(β·DAF) at three
  values of β), and report the true-DAF shift each induces.
- Choose the conditions: unpolarised, reference haplotype, or both.

Code needed in `tools/sim_polarity_arms.py`: weighted focal sampling, the
true-DAF shift report, and an all-Bernoulli arm (A, $B_0$ and nulls all
Bernoulli-$q$). Then a paired power report like `tools/v1_report.py`: per-replicate
power differences, one-sided bootstrap lower bound, and the "at least 10 points
at two or more effect sizes" rule.

Estimate: 1–2 hours of code, then about 1 hour of compute (the V1 arms took about
40 minutes per condition).

## Check the floor correction for $W_1$

`PHI_SFS_SAMPLE_SIZE_BIAS.md` validated
$\sqrt{\Phi_{\mathrm{obs}}^2-\Phi_{\mathrm{floor}}^2}$ for the earlier
total-variation statistic only. Simulate known spectral differences at several
$M$ and check whether it, or $\Phi_{\mathrm{obs}}-\mu_0$, recovers the true $W_1$.

## Merge the release branch into main

The GitHub front page shows `main`. Open a PR from
`bootstrap-target-hpc-validation`, and optionally publish a GitHub release from
the `v0.9.0` tag.
