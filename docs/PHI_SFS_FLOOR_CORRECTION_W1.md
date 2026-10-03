# Floor correction for the $W_1$ Phi-SFS: validation by simulation

Date: 2026-10-02. Script: `tools/validate_floor_correction.py` (commit `2208d11`).
Raw outputs: `logs/floor_correction/` (`meta.json`, `results.json`,
`summary.tsv`; not tracked). Run on one Farm node (job 39398879): 54 CPU-minutes, about 4 minutes on 16 CPUs.

## Question

[PHI_SFS_SAMPLE_SIZE_BIAS.md](PHI_SFS_SAMPLE_SIZE_BIAS.md) validated the
quadrature correction $\sqrt{\Phi_{\mathrm{obs}}^2-\Phi_{\mathrm{floor}}^2}$
for the earlier total-variation statistic. Production now uses the Wasserstein
distance $W_1$ ([METHODS.md](METHODS.md)) and reported
$\Phi_{\mathrm{obs}}-\mu_0$ as the effect size until the correction was
checked for $W_1$. This note checks, against an exactly known true distance
$W_{1,\infty}$, how well four estimates recover it as the set size $M$ changes:

| name | estimate |
|---|---|
| raw | $\Phi_{\mathrm{obs}}=W_1(A,B_0)$ |
| subtraction | $\Phi_{\mathrm{obs}}-\mu_0$ |
| $\hat\Phi_{\mathrm{SFS}}$ (quadrature) | $\sqrt{\max(\Phi_{\mathrm{obs}}^2-\mu_0^2,\,0)}$ |
| extrapolation | route 2 of PHI_SFS_SAMPLE_SIZE_BIAS.md: subsample $A$ and $B_0$ to sizes $L=M/2$, $M/4$, and $M/8$ (levels $L\ge250$; 400 draws each), regress the squared median distance on $1/L$, then take the square root of the intercept |

$\mu_0$ is the mean of $W_1(B_i,B_0)$ over $R=200$ nulls, computed as production
computes it, with the same $B_0$ used for $\Phi_{\mathrm{obs}}$.

## Design

**Site models.** Both use $n_{\mathrm{hap}}=26$ haploids and known polarity.

- **dnAging.** The 241,623 sites with true mutation age in $[10^4, 3\times10^6]$
  generations from `msprime_variable_ne_error` replicates 1–10 of the dnAging
  project (`/quobyte/jrigrp/beil/logan_collab/dnAging/.../per_site.tsv`), the
  same pool as PHI_SFS_SAMPLE_SIZE_BIAS.md (median age 45,695 generations). $B$
  sets are drawn from the pool unweighted. $A$ is drawn with an age weight that
  multiplies the rate of sites younger than the 25% or 50% age quantile by 1.25,
  1.5, 2, 3 or 10, or of sites older than the 75% quantile by 1.5 or 3: 12
  effects plus the null.
- **Toy.** The derived-count weight is
  $f_\theta(d)\propto(1/d)\,e^{\theta d/n_{\mathrm{hap}}}$, as in
  [PHI_SFS_CALIBRATION_VALIDATION.md](PHI_SFS_CALIBRATION_VALIDATION.md), with
  $\theta=0$ for $B$ and $\theta\in\{\pm0.1,\pm0.2,\pm0.4,\pm0.8\}$ for $A$: 8
  effects plus the null.

True distances range from 0.006 to 0.095. The grid uses
$M\in\{100,250,500,1000,2500,5000,10000,20000\}$: 2,000 analyses per cell and
176 cells in total.

**Truth.** $W_{1,\infty}$ is $W_1$ between the two models' normalized expected
spectra, computed exactly. It is not read off a large-$M$ plateau, which was the
weakness of the earlier analysis.

**Production code used unchanged:** `hypergeometric_projection`,
`normalized_spectrum`, `phi_sfs`. The loops use a vectorized $W_1$; a self-check at
start-up confirms it matches `phi_sfs` to $10^{-13}$.

**Simplifications.** Sites are drawn i.i.d., so a set is a multinomial vector of
derived counts. Matching, pool depletion, disjoint sets and polarity
uncertainty are not modelled.

## Results

### 1. The floor scales as $1/\sqrt{M}$

Under the null, $M\mu_0^2$ is constant from $M=100$ to $20{,}000$:

| $M$ | $\mu_0$ (dnAging) | $M\mu_0^2$ (dnAging) | $\mu_0$ (toy) | $M\mu_0^2$ (toy) |
|---:|---:|---:|---:|---:|
| 100 | 0.0360 | 0.129 | 0.0350 | 0.122 |
| 250 | 0.0231 | 0.133 | 0.0220 | 0.121 |
| 500 | 0.0161 | 0.129 | 0.0155 | 0.120 |
| 1,000 | 0.0114 | 0.130 | 0.0109 | 0.119 |
| 2,500 | 0.0072 | 0.128 | 0.0069 | 0.121 |
| 5,000 | 0.0050 | 0.127 | 0.0049 | 0.119 |
| 10,000 | 0.0036 | 0.129 | 0.0035 | 0.120 |
| 20,000 | 0.0025 | 0.130 | 0.0025 | 0.121 |

The constant depends on the spectrum (0.129 versus 0.120), so the floor must
still be computed for the data at hand, as production does.

### 2. Under $W_1$ the floor does not add to a real difference

This is the main result, and it differs from total variation. Mean bias in units
of the floor, by the ratio of the true distance to the floor (all 176 cells):

| $W_1/\mu_0$ | cells | raw | subtraction | $\hat\Phi_{\mathrm{SFS}}$ | extrapolation |
|---|---:|---:|---:|---:|---:|
| 0 (null) | 16 | +1.00 | +0.00 | +0.40 | +0.50 |
| 0–0.5 | 16 | +0.70 | -0.30 | +0.11 | +0.23 |
| 0.5–1 | 20 | +0.45 | -0.55 | -0.09 | +0.03 |
| 1–2 | 34 | +0.17 | -0.83 | -0.25 | -0.16 |
| 2–4 | 34 | +0.02 | -0.98 | -0.22 | -0.10 |
| ≥ 4 | 56 | -0.00 | -1.00 | -0.08 | -0.01 |

- **Raw:** $\Phi_{\mathrm{obs}}$ is biased upward by the full floor under the
  null, but the bias disappears once the true distance exceeds about twice the
  floor. Where the focal CDF lies clearly on one side of the neutral CDF, sampling
  noise moves $|F_A-F_{B_0}|$ up as often as down and averages out, instead of
  adding as it does for total variation.
- **Subtraction** is unbiased only under the null. For any resolvable effect it
  removes a floor that is not there, and underestimates the true distance by
  almost the whole of $\mu_0$.
- **Quadrature:** $\hat\Phi_{\mathrm{SFS}}$ stays within $+0.40$ to $-0.25$
  floors of the truth across the whole range: the smallest worst-case bias of
  the three single-number estimates.
- **Extrapolation** behaves like $\hat\Phi_{\mathrm{SFS}}$ with somewhat smaller
  bias for resolvable effects, but needs $M\ge500$ and is biased $+0.5$ floors
  under the null.

### 3. Relative error by set size

Absolute bias as a percentage of the true distance, median and maximum over the
20 non-null effects:

| $M$ | raw median / max | subtraction median / max | $\hat\Phi_{\mathrm{SFS}}$ median / max | extrapolation median / max |
|---:|---:|---:|---:|---:|
| 100 | 69% / 474% | 82% / 98% | 19% / 132% | — |
| 250 | 23% / 270% | 71% / 89% | 13% / 60% | — |
| 500 | 8% / 162% | 58% / 92% | 14% / 21% | 4% / 64% |
| 1,000 | 4% / 99% | 43% / 82% | 9% / 24% | 4% / 19% |
| 2,500 | 1% / 43% | 30% / 71% | 6% / 18% | 3% / 15% |
| 5,000 | 0% / 18% | 21% / 63% | 3% / 19% | 1% / 17% |
| 10,000 | 0% / 5% | 16% / 52% | 2% / 16% | 0% / 14% |
| 20,000 | 0% / 2% | 10% / 39% | 1% / 10% | 0% / 7% |

At $M\ge500$, $\hat\Phi_{\mathrm{SFS}}$ is within 24% of the truth for every effect
tested, and within 14% for the median effect. Subtraction is off by 43% for the
median effect at $M=1000$. At $M\le250$ no estimate is reliable for effects near
the floor.

### 4. Two examples

**dnAging, 3× rate in youngest quarter**, true $W_1=0.0452$. Bias (RMSE):

| $M$ | $\mu_0$ | raw | subtraction | $\hat\Phi_{\mathrm{SFS}}$ | extrapolation |
|---:|---:|---:|---:|---:|---:|
| 100 | 0.0360 | +0.0076 (0.0281) | -0.0284 (0.0390) | -0.0088 (0.0339) | — |
| 250 | 0.0228 | +0.0016 (0.0203) | -0.0212 (0.0295) | -0.0061 (0.0248) | — |
| 500 | 0.0160 | +0.0001 (0.0152) | -0.0158 (0.0222) | -0.0037 (0.0177) | -0.0007 (0.0163) |
| 1,000 | 0.0113 | +0.0000 (0.0109) | -0.0113 (0.0159) | -0.0017 (0.0119) | -0.0002 (0.0112) |
| 2,500 | 0.0072 | -0.0000 (0.0072) | -0.0073 (0.0104) | -0.0007 (0.0074) | -0.0001 (0.0072) |
| 5,000 | 0.0051 | +0.0000 (0.0049) | -0.0050 (0.0071) | -0.0003 (0.0049) | +0.0000 (0.0049) |
| 10,000 | 0.0036 | -0.0001 (0.0034) | -0.0037 (0.0051) | -0.0003 (0.0034) | -0.0001 (0.0034) |
| 20,000 | 0.0026 | -0.0000 (0.0025) | -0.0026 (0.0037) | -0.0001 (0.0025) | -0.0000 (0.0025) |

**dnAging, 1.5× rate in youngest quarter**, true $W_1=0.0150$. Bias (RMSE):

| $M$ | $\mu_0$ | raw | subtraction | $\hat\Phi_{\mathrm{SFS}}$ | extrapolation |
|---:|---:|---:|---:|---:|---:|
| 100 | 0.0358 | +0.0227 (0.0301) | -0.0130 (0.0220) | +0.0020 (0.0229) | — |
| 250 | 0.0229 | +0.0107 (0.0174) | -0.0122 (0.0178) | -0.0020 (0.0164) | — |
| 500 | 0.0161 | +0.0050 (0.0122) | -0.0111 (0.0153) | -0.0032 (0.0136) | -0.0003 (0.0145) |
| 1,000 | 0.0114 | +0.0024 (0.0091) | -0.0090 (0.0126) | -0.0026 (0.0109) | -0.0013 (0.0113) |
| 2,500 | 0.0072 | +0.0004 (0.0064) | -0.0067 (0.0093) | -0.0019 (0.0078) | -0.0014 (0.0077) |
| 5,000 | 0.0051 | +0.0001 (0.0050) | -0.0050 (0.0071) | -0.0010 (0.0057) | -0.0005 (0.0056) |
| 10,000 | 0.0036 | +0.0002 (0.0035) | -0.0034 (0.0050) | -0.0004 (0.0037) | -0.0000 (0.0036) |
| 20,000 | 0.0025 | -0.0001 (0.0025) | -0.0026 (0.0037) | -0.0003 (0.0026) | -0.0001 (0.0025) |

In every cell where the true distance is at least the floor, subtraction has the
highest RMSE of raw, subtraction and $\hat\Phi_{\mathrm{SFS}}$. Where it is at
least twice the floor, raw or extrapolation has the lowest RMSE and
$\hat\Phi_{\mathrm{SFS}}$ is close behind.

## Conclusions

1. **The quadrature correction carries over to** $W_1$.
   $\hat\Phi_{\mathrm{SFS}}$ is the best single-number estimate of the true
   distance: worst-case bias about $-0.25\mu_0$ for resolvable effects and
   $+0.4\mu_0$ under the null, within 24% of the truth at $M\ge500$ for every
   effect here. The earlier note's claim of "a few percent" error for sets of at
   least 250 sites does not hold for $W_1$ near the floor.
2. **Subtraction should not be used as the effect size.** $\Phi_{\mathrm{obs}}-\mu_0$ is
   unbiased only when there is no effect, and otherwise underestimates the true
   distance by nearly $\mu_0$, so it penalizes small categories most.
3. Raw $\Phi_{\mathrm{obs}}$ is nearly unbiased once the effect is more than
   about twice the floor, so for clearly significant categories raw and
   $\hat\Phi_{\mathrm{SFS}}$ agree. They differ for effects near the floor, where
   $\hat\Phi_{\mathrm{SFS}}$ is less biased.
4. Extrapolation is a useful cross-check for $M\ge1000$, not a replacement.

## Not tested

- Matching, disjoint depletion and the non-identical distribution of production
  control sets. In production $\mu_0$ comes from depleted, matched sets, so it may
  differ from the i.i.d. floor modelled here.
- Polarity uncertainty: production SNP sites use the posterior-weighted mixture
  and TE sites are not masked. Those smooth the projected spectra; I expect the
  $1/\sqrt{M}$ scaling to hold but the constant to change (unverified).
- Effects outside the age-weighting and tilt families above.
