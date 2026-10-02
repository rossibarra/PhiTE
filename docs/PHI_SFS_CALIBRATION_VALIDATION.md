# Φ-SFS calibration: validation by simulation (plan §9.5)

Date: 2026-09-25. Script: `tools/validate_phi_calibration.py`. Raw outputs (per-analysis
arrays and summary JSON) were written to the session scratchpad
(`/tmp/claude-9002/.../scratchpad/calibration/`) and are not part of the repository.

## 1. Purpose

The calibrated test compares an observed distance D_obs = W(A, B0) with R null distances
D_i = W(B_i, B0) and reports an add-one Monte Carlo P-value and a Z-score (plan §2,
§3.4). This note checks, by simulation and independently of the real data, whether that
procedure:

1. gives uniform P-values and nominal type-I error when A and every B set really are
   draws from the same distribution;
2. stays calibrated as M (sites per set) changes, although the null mean distance rises
   as M falls;
3. is affected by two design features flagged in CODE_REVIEW_ROUND10 ("Design
   observations"): the unequal bootstrap-plus-matching error in D_obs and D_i, and the
   different polarity rules for TEs and SNPs; and
4. has power that increases with effect size and M, with `mean_daf_difference` giving
   the right direction.

**Not covered:** the check that sequential pool depletion under `--disjoint-replicates`
leaves type-I error calibrated at 1001 sets (the last bullet of §9.5). It needs the real
matcher and candidate pool. Every simulation here draws sites i.i.d. from an infinite
pool, so it cannot say anything about depletion.

## 2. What is real code and what is a model

**Production code used unchanged** (imported read-only from `normalize_tes.phi_sfs`):
`project_sites` (and through it `hypergeometric_projection`) to project each site to 20
with the SNP q-mixture `q·h(k,n) + (1−q)·h(n−k,n)`, `normalized_spectrum`, `phi_sfs`
(distance and `mean_daf_difference`), and `calibrate_phi` (Z, add-one P). Those five
functions are identical in HEAD (`12361c0`) and in the working tree used for the runs.
Other agents were editing other parts of `phi_sfs.py` at the time.

**Simulated:** everything upstream of the per-site counts, meaning the site
distributions, polarity posteriors, ages, bootstrapping and matching.

**Sampling shortcut.** Every site type (k, n, q) has q on a 0.05 grid, so a set of M
i.i.d. sites is a multinomial vector of type counts, and its raw spectrum is that vector
times the per-type projection matrix. This gives exactly the same distribution as
drawing M sites one at a time. The self-check in the script confirms that it reproduces
`project_sites` on individually listed sites to 1e-12.

## 3. Models and assumptions

### 3.1 Neutral exchangeable model (Part 1)

- **True derived count** d among n callable inbred lines: standard neutral, P(d) ∝ 1/d,
  d = 1…n−1.
- **n:** fixed at 26, the maximum callable count in the production VCF (as cited in
  `PHI_SFS_SAMPLE_SIZE_BIAS.md`). A second model draws n uniformly per site from
  20…40.
- **SNP polarity posterior:** the ARG confidence c = max(q, 1−q) is a mixture: 85% of
  sites have c ~ Beta(40, 1) restricted to [0.5, 1] (mean about 0.976), and 15% have
  c ~ U(0.5, 1). With probability 0.7 the ARG favours ALT (q_ALT = c), otherwise
  q_ALT = 1−c. The posterior is assumed **calibrated**: ALT is truly derived with
  probability q_ALT, so k_ALT = d or n−d. Under this model the wrong-call rate
  E[1−c] is 5.5%. The mis-orientation mass the q-mixture carries,
  m_B = 2E[c(1−c)], is 8.27%.
- **A and all B_r** are drawn from this same distribution (a SNP-vs-SNP negative
  control). This is the exchangeable null, so the procedure should be exactly calibrated
  here apart from Monte Carlo error.

### 3.2 Stylised age-covariate model of the bootstrap asymmetry (Part 2)

This model is **stylised**. It isolates one mechanism and is not a model of the maize
data.

- **Ages:** 20 age bins. P(d | age bin g) ∝ (1/d)·exp(s·z_g·d/n), with z_g running from
  −0.5 (youngest) to +0.5 (oldest), n = 26, and SNP polarity as in §3.1. Two strengths
  of the age–frequency relationship were used. At s = 4, expected DAF runs from 0.23 in
  the youngest bin to 0.44 in the oldest. At s = 12 it runs from 0.17 to 0.70.
- **A:** M = 4000 site ages drawn from a young-skewed distribution, P(g) ∝ exp(−g/6).
  Each A frequency is drawn given its age.
- **Each B_r (r = 0…R):** multinomial-bootstrap A's age histogram to get the target
  h*_r, then produce the matched age histogram as h*_r·exp(σ·ξ_r + b·z), rounded back
  to M sites (largest remainder), with ξ_r ~ N(0,1) independently per bin and per
  replicate. B_r frequencies are then drawn given the matched ages.
  - σ is independent, zero-mean matching error.
  - b is a **systematic** matching bias shared by all B_r (b > 0 means the matched sets
    are older than their target). This was added beyond the brief because random error
    and shared bias act in opposite directions.
- **Reference arm, "observed target":** every B_r takes A's exact age histogram, with no
  bootstrap and no matching error. A and the B_r are then exchangeable given the ages.
- **Size of the perturbation:** reported as the median of E_r/V_r, where
  V_r = W1(h*_r, h_A) and E_r = W1(matched_r, h*_r) are measured on age-bin units. This
  mirrors the production QC ratio `R_r` (median 0.181, gate < 0.5, in
  BOOTSTRAP_HPC_VALIDATION T5). It does so only by analogy: the production ratio is
  computed in generations on a 36,746-point grid, and production E_r has a floor
  (223–337 generations) rather than being proportional to V_r.

### 3.3 Stylised TE-vs-SNP polarity model (Part 3)

- **A (TE-like):** weight-1 polarity, so p_alt_derived = 1. A fraction ε of the retained
  insertions are truly ancestral, so the observed count is n−d.
- **B:** SNP q-mixture as in §3.1.
- **True d:** identically neutral for both A and B.
- **Stylised:** ε is a free parameter. ε = 0.0827 is the value at which A's expected
  mis-orientation matches B's (m_B), and the two expected spectra then agree to
  W = 1.6e-5. ε = 0.09 roughly corresponds to the review's "about 91% accurate", but
  that figure measures ARG accuracy on TEs and is not a measured TE mis-orientation rate
  after the ≥50%-derived filter.

### 3.4 Power model (Part 4)

A is drawn from the SNP model of §3.1 with the true-count distribution tilted:
P_A(d) ∝ (1/d)·exp(θ·d/n). θ < 0 gives an excess of rare variants and θ > 0 a shift
toward high frequency. B is drawn from the untilted model, and A and B use the same
q-mixture. The tilts were fixed before the full run from a 2-analysis pilot, to span the
power range. Because the tilt is monotone the two CDFs never cross, so the expected
W equals |expected mean-DAF difference|.

### 3.5 Replication

- **Replicates R:** 200 for the sweeps, with the smallest attainable P of 1/201 and exact
  attainable levels 10/201 = 0.0498 and 2/201 = 0.00995. Production R = 1000 was used
  for three null configurations.
- **Analyses per configuration:** 2000 (null, n = 26, R = 200), 1000 (other null and
  Parts 2–3), or 500 (power).
- **Seeds:** independent per configuration and chunk.
- **Intervals:** all type-I and power intervals are 95% Clopper–Pearson.

## 4. Results

### 4.1 Exchangeable null: P-values, type-I error, Z (measured)

| config | n | M | R | analyses | reject α=0.05 [95% CI] | reject α=0.01 [95% CI] | KS P | decile χ² P | Z mean ± SE | Z sd |
|---|---|---:|---:|---:|---|---|---:|---:|---|---:|
| n26 M500 | 26 | 500 | 200 | 2000 | 0.0465 [0.038, 0.057] | 0.0115 [0.007, 0.017] | 0.51 | 0.36 | −0.006 ± 0.023 | 1.015 |
| n26 M1000 | 26 | 1000 | 200 | 2000 | 0.0485 [0.040, 0.059] | 0.0110 [0.007, 0.017] | 0.16 | 0.52 | 0.034 ± 0.023 | 1.010 |
| n26 M4000 | 26 | 4000 | 200 | 2000 | 0.0475 [0.039, 0.058] | 0.0060 [0.003, 0.011] | 0.41 | 0.25 | −0.004 ± 0.022 | 1.001 |
| n26 M20000 | 26 | 20000 | 200 | 2000 | 0.0540 [0.045, 0.065] | 0.0095 [0.006, 0.015] | 0.32 | 0.58 | 0.025 ± 0.022 | 0.998 |
| n20–40 M500 | 20–40 | 500 | 200 | 1000 | 0.0640 [0.050, 0.081] | 0.0110 [0.006, 0.020] | 0.85 | 0.52 | 0.005 ± 0.032 | 1.027 |
| n20–40 M1000 | 20–40 | 1000 | 200 | 1000 | 0.0550 [0.042, 0.071] | 0.0160 [0.009, 0.026] | 0.12 | 0.63 | −0.019 ± 0.034 | 1.065 |
| n20–40 M4000 | 20–40 | 4000 | 200 | 1000 | 0.0650 [0.051, 0.082] | 0.0100 [0.005, 0.018] | 0.62 | 0.85 | 0.000 ± 0.033 | 1.055 |
| n20–40 M20000 | 20–40 | 20000 | 200 | 1000 | 0.0430 [0.031, 0.058] | 0.0090 [0.004, 0.017] | 0.57 | 0.46 | −0.013 ± 0.031 | 0.979 |
| **n26 M500, R=1000** | 26 | 500 | 1000 | 1000 | 0.0450 [0.033, 0.060] | 0.0100 [0.005, 0.018] | **0.0004** | 0.008 | −0.068 ± 0.032 | 1.008 |
| ↳ replicate seed 1 | 26 | 500 | 1000 | 1000 | 0.0580 [0.044, 0.074] | 0.0160 [0.009, 0.026] | 0.49 | 0.77 | 0.023 ± 0.033 | 1.032 |
| ↳ replicate seed 2 | 26 | 500 | 1000 | 1000 | 0.0660 [0.051, 0.083] | 0.0130 [0.007, 0.022] | 0.75 | 0.63 | 0.014 ± 0.033 | 1.028 |
| **n26 M4000, R=1000** | 26 | 4000 | 1000 | 1000 | 0.0450 [0.033, 0.060] | 0.0120 [0.006, 0.021] | 0.38 | 0.35 | 0.018 ± 0.031 | 0.990 |
| **n26 M20000, R=1000** | 26 | 20000 | 1000 | 1000 | 0.0420 [0.030, 0.056] | 0.0090 [0.004, 0.017] | 0.35 | 0.29 | −0.036 ± 0.031 | 0.970 |

**Pooled** over the 11 primary configurations (15,000 analyses): type-I error was
**0.0501 [0.0467, 0.0537]** at α = 0.05 and **0.0102 [0.0087, 0.0119]** at α = 0.01.

Decile counts of P at production R = 1000 (expected 100 per decile):

- M = 4000: 98, 97, 96, 120, 101, 107, 98, 88, 111, 84.
- M = 20000: 85, 108, 95, 98, 91, 98, 104, 122, 109, 90.

**The M = 500, R = 1000 outlier.** The first run's KS P = 0.0004 is the only departure
among 14 null runs. The deciles show a mild deficit at P < 0.4 and an excess at 0.5–0.8,
which is the conservative direction, and its rejection rates are nominal. Two
independent re-runs of the identical configuration were clean (KS P = 0.49 and 0.75).
Pooling all three runs gives type-I error 169/3000 = 0.056 [0.048, 0.065]. I read the
outlier as probably a chance fluctuation, but three runs cannot rule out a small
seed-specific or M-specific effect.

**Null mean distance against M.** The mean null distance falls as 1/√M. Mean
null_mean·√M was 0.337–0.346 in every configuration:

| M | mean of null_mean | mean of null_sd |
|---:|---:|---:|
| 500 | 0.0153 | 0.0066 |
| 1000 | 0.0108 | 0.0047 |
| 4000 | 0.0054 | 0.0023 |
| 20000 | 0.0024 | 0.0011 |

Type-I error did not trend with M. So the "floor" rises roughly 6-fold from M = 20000 to
M = 500, and calibrating against same-M nulls absorbs it, as the design intends.

**Z-score scale.** Z_obs had mean ≈ 0 and sd ≈ 1 in every null configuration (0.97–1.07),
but the null distance distribution is right-skewed (mean skewness 0.88–0.92). As a
result, reading Z against a standard normal is anticonservative:

- Pooled P(Z ≥ 1.645) was **0.071 [0.067, 0.075]**, against a nominal 0.05.
- P(Z ≥ 2.326) was 0.021–0.037 per configuration, against a nominal 0.01.
- The empirical 95th percentile of null Z at M = 4000, R = 1000 was 1.80.

Z is a useful effect-size scale, but the Monte Carlo P-value is the test statistic.
Normal-theory P-values should not be quoted from Z.

### 4.2 Bootstrap-target asymmetry (stylised model)

All configurations use M = 4000 and R = 200 with 1000 analyses, unless marked otherwise.

| age effect | B-set generation | median E/V | reject α=0.05 [95% CI] | reject α=0.01 [95% CI] | Z mean ± SE | Z sd | D_obs / null mean |
|---|---|---:|---|---|---|---:|---:|
| s=4 | observed target (reference) | — | 0.062 [0.048, 0.079] | 0.013 [0.007, 0.022] | 0.014 ± 0.033 | 1.05 | 1.006 |
| s=4 | bootstrap, σ=0 | 0 | 0.030 [0.020, 0.043] | 0.010 [0.005, 0.018] | −0.053 ± 0.030 | 0.94 | 0.981 |
| s=4 | bootstrap, σ=0.01 | 0.14 | 0.050 [0.037, 0.065] | 0.004 [0.001, 0.010] | −0.012 ± 0.031 | 0.99 | 0.993 |
| s=4 | bootstrap, σ=0.02 | 0.29 | 0.042 [0.030, 0.056] | 0.008 [0.004, 0.016] | −0.033 ± 0.031 | 0.97 | 0.987 |
| s=4 | bootstrap, σ=0.035 | 0.50 | 0.054 [0.041, 0.070] | 0.005 [0.002, 0.012] | 0.026 ± 0.032 | 1.01 | 1.010 |
| s=4 | bootstrap, σ=0.07 | 1.01 | 0.051 [0.038, 0.067] | 0.014 [0.008, 0.023] | −0.012 ± 0.032 | 1.02 | 0.991 |
| s=4 | bootstrap, σ=0.15 | 2.16 | 0.036 [0.025, 0.050] | 0.002 [0.000, 0.007] | −0.068 ± 0.029 | 0.91 | 0.976 |
| s=12 | observed target (reference) | — | 0.052 [0.039, 0.068] | 0.014 [0.008, 0.023] | 0.008 ± 0.032 | 1.03 | 1.002 |
| s=12 | bootstrap, σ=0 | 0 | 0.023 [0.015, 0.034] | 0.003 [0.001, 0.009] | −0.134 ± 0.028 | 0.87 | 0.941 |
| s=12 | bootstrap, σ=0.01 | 0.14 | 0.030 [0.020, 0.043] | 0.007 [0.003, 0.014] | −0.061 ± 0.029 | 0.92 | 0.976 |
| s=12 | bootstrap, σ=0.02 | 0.29 | 0.023 [0.015, 0.034] | 0.003 [0.001, 0.009] | −0.123 ± 0.028 | 0.88 | 0.950 |
| s=12 | bootstrap, σ=0.035 | 0.50 | 0.028 [0.019, 0.040] | 0.003 [0.001, 0.009] | −0.128 ± 0.028 | 0.88 | 0.946 |
| s=12 | bootstrap, σ=0.07 | 1.01 | 0.021 [0.013, 0.032] | 0.000 [0, 0.004] | −0.174 ± 0.027 | 0.84 | 0.929 |
| s=12 | bootstrap, σ=0.15 | 2.16 | 0.009 [0.004, 0.017] | 0.000 [0, 0.004] | −0.281 ± 0.022 | 0.70 | 0.874 |
| s=12, M=1000 | bootstrap, σ=0 | 0 | 0.028 [0.019, 0.040] | 0.004 [0.001, 0.010] | −0.077 ± 0.029 | 0.91 | 0.969 |
| s=12, M=20000 | bootstrap, σ=0 | 0 | 0.032 [0.022, 0.045] | 0.005 [0.002, 0.012] | −0.107 ± 0.029 | 0.92 | 0.956 |
| s=12 | bootstrap, σ=0.01, **shared bias b=0.05** | 0.74 | 0.032 [0.022, 0.045] | 0.005 [0.002, 0.012] | −0.064 ± 0.030 | 0.94 | 0.971 |
| s=12 | bootstrap, σ=0.01, **shared bias b=0.1** | 1.50 | 0.064 [0.050, 0.081] | 0.013 [0.007, 0.022] | 0.087 ± 0.034 | 1.08 | 1.029 |
| s=12 | bootstrap, σ=0.01, **shared bias b=0.2** | 3.05 | 0.227 [0.201, 0.254] | 0.074 [0.059, 0.092] | 0.809 ± 0.047 | 1.49 | 1.329 |
| s=12 | bootstrap, σ=0.01, **shared bias b=0.4** | 6.29 | 0.644 [0.613, 0.674] | 0.423 [0.392, 0.454] | 2.662 ± 0.064 | 2.02 | 2.130 |

**What was measured, in this model:**

- **Random error.** The asymmetry makes the test conservative, as the review predicted.
  The size depends on how strongly age predicts frequency.
  - **Strong age effect (s = 12)**, in the six configurations with E/V ≤ 0.5 (σ ≤ 0.035
    at M = 4000, plus σ = 0 at M = 1000 and 20000): pooled type-I error at α = 0.05 is
    **0.027 [0.023, 0.032]**, about half nominal.
  - **Weaker age effect (s = 4)**, bootstrap configurations with σ ≤ 0.07: pooled type-I
    error is 0.045 [0.040, 0.052]. That is not clearly different from 0.05, although
    individual cells range from 0.030 to 0.054.
  - **Reference arm.** The observed-target configurations are nominal (0.052 and 0.062),
    which confirms that the conservatism comes from the per-replicate bootstrap.
  - **Larger random error** makes the test more conservative still (s = 12, σ = 0.15,
    E/V ≈ 2.2: 0.009). It never makes it anticonservative.
- **Systematic matching bias** works the other way. If every B_r is displaced in the
  same direction from A's age distribution, A differs from all the B's in a way the
  B's do not differ from each other, and type-I error inflates:
  - 0.032 at E/V = 0.74;
  - 0.064 at E/V = 1.5;
  - 0.23 at E/V = 3.1;
  - 0.64 at E/V = 6.3.

  In this model the conservatism from the bootstrap asymmetry absorbs a shared bias up
  to about E/V ≈ 1. Below the production QC gate (E/V < 0.5), I saw no inflation.

**What this does not establish:**

- How strongly real allele age predicts DAF in the maize ARG, relative to s = 4 or 12.
- Whether production matching error is random or shared across replicates. The T5 floor
  of 223–337 generations, present in every replicate, could be either.
- Whether the E/V scale here corresponds to the production ratio.

The cheapest production check is to test whether the signed matching residual
(G_{B_r} − G*_r) has a consistent sign across replicates, and to report the SNP-vs-SNP
negative-control P-value, which carries both effects.

### 4.3 TE-vs-SNP polarity rules (stylised model)

A uses weight-1 polarity with a fraction ε truly flipped. B uses the calibrated
q-mixture (m_B = 0.0827). All configurations use R = 200 and 1000 analyses.

| ε | expected W (A vs B) | M | reject α=0.05 [95% CI] | Z mean ± SE | Z sd | fraction with mean_daf_difference > 0 |
|---:|---:|---:|---|---|---:|---:|
| 0 | 0.0384 | 1000 | 0.936 [0.919, 0.950] | 6.22 ± 0.09 | 2.69 | 0.000 |
| 0 | 0.0384 | 4000 | 1.000 [0.996, 1] | 14.67 ± 0.12 | 3.93 | 0.000 |
| 0 | 0.0384 | 20000 | 1.000 [0.996, 1] | 36.22 ± 0.25 | 7.87 | 0.000 |
| 0.03 | 0.0245 | 1000 | 0.710 [0.681, 0.738] | 3.31 ± 0.07 | 2.35 | 0.017 |
| 0.03 | 0.0245 | 4000 | 0.979 [0.968, 0.987] | 8.54 ± 0.10 | 3.11 | 0.000 |
| 0.03 | 0.0245 | 20000 | 1.000 [0.996, 1] | 22.55 ± 0.17 | 5.40 | 0.000 |
| 0.0827 (matched) | 1.6e-5 | 1000 | 0.079 [0.063, 0.098] | 0.13 ± 0.04 | 1.14 | 0.53 |
| 0.0827 (matched) | 1.6e-5 | 4000 | 0.068 [0.053, 0.085] | 0.11 ± 0.04 | 1.13 | 0.52 |
| 0.0827 (matched) | 1.6e-5 | 20000 | 0.087 [0.070, 0.106] | 0.21 ± 0.04 | 1.17 | 0.50 |
| 0.09 | 0.0034 | 1000 | 0.098 [0.080, 0.118] | 0.23 ± 0.04 | 1.25 | 0.63 |
| 0.09 | 0.0034 | 4000 | 0.143 [0.122, 0.166] | 0.48 ± 0.04 | 1.37 | 0.74 |
| 0.09 | 0.0034 | 20000 | 0.421 [0.390, 0.452] | 1.65 ± 0.06 | 1.94 | 0.89 |

**What was measured, in this model:** two separate effects.

- **Bias.** When A's mis-orientation rate differs from the SNP mixture's, the difference
  looks to the test like a real SFS difference.
  - A perfectly polarized TE set (ε = 0) against a q-mixed SNP set is "detected" in 94%
    of analyses at M = 1000 and in every analysis at M ≥ 4000.
  - The spurious signal is always toward an **apparent excess of rare variants** in A,
    because the SNP mixture is partly folded. It grows with M, like any real effect.
  - A mismatch of only 0.7 percentage points (ε = 0.09 against m_B = 0.083) gives 42%
    rejection at M = 20000.
- **Variance.** This effect is independent of M. Even when the expected spectra are
  matched (ε = m_B), type-I error is **0.068–0.087** (pooled 234/3000 = 0.078), with
  Z sd ≈ 1.15. A weight-1 site is a "hard" 0/1 orientation, so an A set has more
  site-to-site variance than a B set of fractional q-mixtures. D_obs is therefore
  stochastically larger than the B_i-to-B0 distances.

**Interpretation, which is an inference:** a TE-vs-SNP Z-score by itself does not
separate biology from polarity method. The SNP-vs-SNP control (`-A SNP`) removes both
effects by construction, since A and B use the same rule, so the review's
recommendation looks right. The size of both effects in the real data depends on
unmeasured quantities: the TE mis-orientation rate after the ≥50% filter, and how
calibrated the ARG q is. This model shows the mechanism, not the real magnitude.

### 4.4 Power (exchangeable SNP-vs-SNP model plus tilt)

Rejection at α = 0.05 [95% CI] and mean Z. Each cell has 500 analyses with R = 200.
The expected shift in mean DAF (A − B) is given for each θ. The last column, "sign
correct", is the fraction of analyses in which `mean_daf_difference` has the sign of
the true shift.

| θ | expected Δ mean DAF | M=1000 power | Z | M=4000 power | Z | M=20000 power | Z | sign correct (1000 / 4000 / 20000) |
|---:|---:|---|---:|---|---:|---|---:|---|
| −0.8 | −0.0365 | 0.938 [0.913, 0.958] | 5.90 | 1.000 [0.993, 1] | 14.15 | 1.000 [0.993, 1] | 34.04 | 1.00 / 1.00 / 1.00 |
| −0.4 | −0.0194 | 0.542 [0.497, 0.586] | 2.17 | 0.962 [0.941, 0.977] | 6.44 | 1.000 [0.993, 1] | 17.27 | 0.97 / 1.00 / 1.00 |
| −0.2 | −0.0100 | 0.228 [0.192, 0.267] | 0.71 | 0.576 [0.531, 0.620] | 2.35 | 0.970 [0.951, 0.983] | 7.73 | 0.84 / 0.97 / 1.00 |
| −0.1 | −0.0051 | 0.098 [0.073, 0.128] | 0.19 | 0.204 [0.170, 0.242] | 0.73 | 0.686 [0.643, 0.727] | 2.93 | 0.68 / 0.87 / 0.98 |
| +0.1 | +0.0052 | 0.096 [0.072, 0.125] | 0.24 | 0.248 [0.211, 0.288] | 0.87 | 0.706 [0.664, 0.746] | 3.08 | 0.68 / 0.86 / 0.99 |
| +0.2 | +0.0106 | 0.224 [0.188, 0.263] | 0.70 | 0.612 [0.568, 0.655] | 2.45 | 0.982 [0.966, 0.992] | 8.23 | 0.82 / 0.98 / 1.00 |
| +0.4 | +0.0217 | 0.640 [0.596, 0.682] | 2.73 | 0.946 [0.922, 0.964] | 7.40 | 1.000 [0.993, 1] | 19.23 | 0.98 / 1.00 / 1.00 |
| +0.8 | +0.0458 | 0.976 [0.959, 0.988] | 8.00 | 1.000 [0.993, 1] | 18.09 | 1.000 [0.993, 1] | 43.63 | 1.00 / 1.00 / 1.00 |

**What was measured:**

- Power and mean Z increase monotonically with |θ| at every M, and with M at every θ.
- Mean Z grows somewhat faster than √M: at θ = −0.8 it is 5.9, 14.2 and 34.0. That is
  what is expected if Z ≈ (W_true − null mean)/null sd, with both null terms ∝ 1/√M.
- Among the 8,317 analyses that reject at α = 0.05, `mean_daf_difference` has the
  correct sign in 8,310 (99.9%). The worst cell is 0.92, at M = 1000 with |θ| = 0.1,
  where power is about 0.1. Across all analyses, rejecting or not, the sign is correct in
  68–100% of cases, rising with effect size and M.
- A shift of 0.01 in mean DAF (θ = ±0.2) is detected with power of about 0.6 at
  M = 4000, about the in-gene target size, and about 0.98 at M = 20000.
- For the same |θ|, high-frequency shifts are marginally easier to detect than
  rare-variant excess. The expected shift is also slightly larger for θ > 0.

These power figures apply to this tilt family, to n = 26, and to the exchangeable null.
Under the conservatism of §4.2 I expect power to be somewhat lower, but that was not
measured.

## 5. Conclusions

**Measured** (simulation, exchangeable null, production scoring code):

1. P-values from `calibrate_phi` are uniform, and type-I error is nominal across
   M = 500–20000, R = 200 and 1000, and fixed or variable n. Pooled type-I error was
   0.0501 [0.047, 0.054] at α = 0.05 and 0.0102 [0.009, 0.012] at α = 0.01.
2. The null mean distance rises as 1/√M (null_mean·√M ≈ 0.34), and calibration absorbs
   it: type-I error does not trend with M.
3. Z_obs has mean 0 and sd 1 under the null, but the null distribution is right-skewed.
   Normal-theory tail probabilities from Z are anticonservative: P(Z ≥ 1.645) was 0.071.
   Report the Monte Carlo P-value, not a normal P-value from Z.
4. Power and Z increase with effect size and M, and `mean_daf_difference` recovers the
   direction of rare-variant and high-frequency shifts.

**Stylised findings** (mechanisms demonstrated; production magnitude unknown):

5. Bootstrapping each B_r's age target makes the test **conservative** when age predicts
   frequency. Type-I error at α = 0.05 fell to about 0.027 with a strong age–frequency
   relationship and about 0.045 with a weaker one. Independent matching error adds to
   the conservatism.
6. A matching bias **shared** by all B_r pushes the other way and inflates type-I error.
   In this model the inflation appeared only at matching-error ratios above the
   production QC gate (E/V ≳ 1). Whether production matching error is shared across
   replicates has not been checked.
7. Different polarity rules for TEs (weight 1) and SNPs (q-mixture) can by themselves
   produce a significant TE-vs-SNP Φ-SFS. They do so through a bias that grows with M
   and points toward apparent rare-variant excess in the TEs, and through a variance
   mismatch that gave type-I error of about 0.08 even when the mean polarity error was
   matched. This supports interpreting TE-vs-SNP results against the `-A SNP` negative
   control.

**Not tested:** sequential pool depletion under `--disjoint-replicates` at 1001 sets,
the finite-pool part of §9.5, which needs the real matcher; real ARG q calibration; and
real age–frequency structure.

## 6. Reproduction and runtime

```
python -m tools.validate_phi_calibration check --out DIR     # self-checks only
python -m tools.validate_phi_calibration pilot --out DIR --workers 8
python -m tools.validate_phi_calibration run --parts null,asym,polar,power --out DIR --workers 8
python -m tools.validate_phi_calibration report --out DIR
```

Runs used `~/.claude/bin/hpc_run` with HPC_CPUS=8 and HPC_MEM=16G, on the `low`
partition, one job per part. Base seed 20260925; the M = 500 re-runs used seeds 1 and
2.

**Pilot.** Two analyses per configuration predicted about 3.0 CPU-h for the original
grid.

**Measured runtime:**

| job | wall | CPU (sum over workers) | peak RSS |
|---|---:|---:|---:|
| null | 7.7 min | 3660 s | 90 MB |
| asym | 10.5 min | 4974 s | 64 MB |
| polar | 3.8 min | 1792 s | 59 MB |
| power | 3.6 min | 1673 s | 60 MB |
| null re-runs (seeds 1, 2) | 2.0 min | 907 s | – |
| shared-bias arm | 2.8 min | 1309 s | – |
| **total** | – | **≈ 4.0 CPU-h** (≈ 4.05 CPU-h allocated) | – |
