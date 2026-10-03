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

## Production-matched floor correction

See [PHI_SFS_FLOOR_CORRECTION_W1.md](PHI_SFS_FLOOR_CORRECTION_W1.md). Quadrature
is now the recommended effect size. The completed check used i.i.d. sites;
repeat it with matched, depleted control sets to test whether production $\mu_0$
behaves like the i.i.d. floor.
