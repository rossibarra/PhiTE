# Phi-SFS is biased upward at small target sizes

## Summary

Phi-SFS is a distance, so it is strictly positive between any two finite sets of
sites even when the underlying spectra are identical. The size of that floor
depends on how many sites went into the comparison. Two Phi-SFS values computed
from different numbers of sites are therefore **not directly comparable**, and
the smaller set will score higher for no biological reason.

This document defines the quantities involved, gives the size of the effect, and
records two validated ways to correct for it.

## Notation

| symbol | meaning |
|---|---|
| `n_TE` | number of variant sites in a target set, and in each matched control set |
| `Phi_obs` | the Phi-SFS actually computed from finite data — what `normalize_tes.phi_sfs` returns |
| `Phi_inf` | the asymptotic Phi-SFS: total variation between the two *true* spectra, the value `Phi_obs` converges to as `n_TE` grows. The estimand. |
| `Phi_floor` | `Phi_obs` between two same-size sets drawn from the *same* spectrum; pure sampling noise |

`n_TE` is a count of sites, not of samples. It is unrelated to the projection
size `m = 20` of `PHI_SFS_IMPLEMENTATION_PLAN.md` section 2.1, which is a number
of individuals and fixes the 19 spectrum bins. `Phi_inf` is likewise unrelated to
the `D(.,.)` age-distribution distance of `BOOTSTRAP_HPC_VALIDATION.md`, which is
measured in generations.

## Why the bias exists

Phi-SFS is the total variation distance between two normalized 19-bin spectra
(`normalize_tes/phi_sfs.py`). Each spectrum is estimated from a finite number of
sites, so each carries multinomial sampling error. Total variation between two
noisy estimates cannot be smaller than that noise, so

```
E[Phi_obs] > Phi_inf   always, with equality only in the limit.
```

The bias is one-sided and does not average out across replicates or across
matched control sets.

## Size of the floor

Measured on 241,623 neutral sites (see Provenance), drawing two independent
same-size sets from one spectrum:

| `n_TE` | `Phi_floor` |
|---:|---:|
| 100 | 0.1248 |
| 250 | 0.0781 |
| 500 | 0.0540 |
| 1,000 | 0.0410 |
| 2,500 | 0.0251 |
| 5,000 | 0.0175 |
| 10,000 | 0.0132 |

`Phi_floor^2` scales as `1/n_TE`: the product `n_TE * Phi_floor^2` stays within
1.46-1.74 across the whole range.

These are indicative magnitudes, not a lookup table for production data. The
floor depends on the spectrum's own shape and on the callable sample size, so
compute it for the data at hand rather than reading it from here.

### What this costs in resolution

Against a known spectral difference of `Phi_inf = 0.1456` — the effect of a
tenfold recent rise in insertion rate:

| `n_TE` | `Phi_floor` | `Phi_obs` | ratio |
|---:|---:|---:|---:|
| 100 | 0.1248 | 0.1759 | 1.41 |
| 250 | 0.0781 | 0.1611 | 2.06 |
| 500 | 0.0540 | 0.1598 | 2.96 |
| 1,000 | 0.0410 | 0.1487 | 3.63 |
| 2,500 | 0.0251 | 0.1507 | 6.00 |
| 5,000 | 0.0175 | 0.1475 | 8.41 |
| 10,000 | 0.0132 | 0.1456 | 10.99 |

At `n_TE = 100` a large perturbation is only 1.4x the noise. At `n_TE >= 500`
it is comfortably resolved. Treat `n_TE = 250` as the lower bound for a group
worth interpreting, and `n_TE = 100` as uninterpretable on its own.

## When this matters

It does **not** matter when every quantity being compared shares one `n_TE`.
Comparing a TE target against its matched controls at fixed `n_TE` is unaffected:
both sides carry the same floor.

It **does** matter whenever `n_TE` differs across the things being compared:

- TE categories of different abundance;
- distance-to-gene or chromatin bins, where TE density varies systematically;
- any rate-change or selection comparison, since both change how many TEs segregate;
- comparisons across replicates with variable site counts.

In all of these the group with fewer sites scores higher, and the direction of
that artifact often coincides with the biological hypothesis.

Note that the spread across the 100 matched control sets does **not** protect
against this. Those error bars vary the control set while holding the TE target
fixed, so they capture control-resampling noise only. The floor is driven by
finite sampling in both sets, and a tight error bar around a value that is
entirely floor is a consistent outcome.

## Route 1: quadrature correction

`Phi_obs` behaves close to a quadrature sum of the true distance and the floor:

```
Phi_obs^2 ~= Phi_inf^2 + Phi_floor^2
  =>  Phi_inf_hat = sqrt(max(Phi_obs^2 - Phi_floor^2, 0))
```

Accuracy against a known `Phi_inf = 0.1456`:

| `n_TE` | `Phi_inf_hat` | error |
|---:|---:|---:|
| 100 | 0.1240 | -14.9% |
| 250 | 0.1409 | -3.2% |
| 500 | 0.1504 | +3.3% |
| 1,000 | 0.1429 | -1.8% |
| 2,500 | 0.1486 | +2.1% |
| 5,000 | 0.1465 | +0.6% |
| 10,000 | 0.1450 | -0.4% |

Good to a few percent for `n_TE >= 250`; unusable at 100.

### Estimating `Phi_floor` correctly

Do **not** estimate `Phi_floor` from pairs of the 100 matched control sets. They
are drawn from one candidate pool and matched to the same target, so they
overlap. Overlapping sets resemble each other more than independent draws do,
which deflates `Phi_floor` and leaves `Phi_inf_hat` too high.

Instead split the control pool into two disjoint halves and draw an age-matched
set of `n_TE` sites from each; Phi between those is the honest floor. If the pool
cannot supply two disjoint matched sets, control reuse is already heavy enough to
matter, which is what the reuse diagnostics of
[BOOTSTRAP_TARGET_MATCHING_PLAN.md](BOOTSTRAP_TARGET_MATCHING_PLAN.md) are for.

## Route 2: subsample and extrapolate

Since `Phi_obs^2` is linear in `1/n_TE`, regressing it on `1/n_TE` gives
`Phi_inf^2` as the intercept. This needs no separate floor estimate.

1. Take the target and its matched control sets at their full `n_TE`.
2. Subsample both to `n_TE/2`, `n_TE/4`, `n_TE/8`, keeping only levels at or
   above 250 sites.
3. Draw ~1,000 subsamples at each level and take the median `Phi_obs`.
4. Regress the squared medians on `1/n_TE`.
5. `Phi_inf_hat = sqrt(intercept)`.

Accuracy against the same known `Phi_inf = 0.1456`:

| points used | intercept | `Phi_inf_hat` | error |
|---|---:|---:|---:|
| all 7 levels | 0.02195 | 0.1481 | +1.8% |
| `n_TE > 100` | 0.02164 | 0.1471 | +1.0% |
| `n_TE > 250` | 0.02120 | 0.1456 | 0.0% |

Route 2 is preferred where the data allow it. It avoids the `Phi_floor`
estimation pitfall entirely, and it is self-diagnosing: if `Phi_obs^2` is not
linear in `1/n_TE`, the fit shows it and the correction should not be applied.

### Cost

Seconds per group on one core. Measured on the reference data:

| `n_TE` | levels | Phi computations | time |
|---:|---|---:|---:|
| 500 | 500, 250 | 2,000 | 0.52 s |
| 2,000 | 2000, 1000, 500, 250 | 4,000 | 1.52 s |
| 8,000 | 8000, 4000, 2000, 1000 | 4,000 | 4.07 s |

Two things keep it cheap. Each site's 19-bin hypergeometric projection depends
only on its `(k, n)` pair, so it is computed once into an `(n_TE, 19)` matrix and
every subsample is a sum of rows. And subsampling operates on the
**already-matched** target and control sets, so nothing upstream of
`normalize_tes.phi_sfs` reruns — no ARG inference, no interval store, no
bootstrap matching. A set-level rather than pairwise match is preserved in
expectation, since both sides are subsampled without bias.

### Limitation

Subsampling only goes downward, and the fit extrapolates to `1/n_TE = 0`, outside
the observed range. A group with `n_TE = 8,000` yields four levels spanning 8x
and a solid estimate; a group with `n_TE = 500` yields two levels, the bare
minimum for a two-parameter fit, with no way to check linearity. **The groups
that most need the correction are the ones where Route 2 has the least room to
work.** For those, fall back to Route 1 with a disjoint-pool floor, and treat
agreement between the two routes as the evidence that either is trustworthy.

Given the cost, run both on every group. Disagreement between them is itself a
diagnostic.

## Recommended reporting

- Report `n_TE` for every group. A Phi-SFS value without its site count cannot be
  interpreted.
- Report `Phi_floor` alongside `Phi_obs`, as a band or a second series. Points on
  the band are uninformative.
- For comparisons across groups, either correct to `Phi_inf` by one of the routes
  above, or subsample every group to a common `n_TE` and say which was used.
- Do not subtract `Phi_floor` from `Phi_obs`. The two do not add; use the
  quadrature form.

## Provenance

Every number here was measured on simulated neutral data, not on maize:

- Source: `msprime_variable_ne_error` replicates 1-10 of the `dnAging` project,
  10 Mb, `mu = r = 1e-8`, piecewise-constant Ne, 26 modern haploid samples,
  SINGER posteriors with 100 MCMC draws.
- Sites restricted to true mutation ages in `[1e4, 3e6]` generations, the range
  of interest for maize TEs: 241,623 sites, median true age 45,695 generations.
- Spectra projected 26 -> 20 individuals and scored with `normalize_tes.phi_sfs`,
  the production statistic.
- The known `Phi_inf = 0.1456` was produced by imposing a tenfold recent rise in
  insertion rate through age-weighted resampling. Thinning a Poisson process by a
  time-dependent weight reproduces the point process of a genuine rate change
  exactly, so the imposed history is not an approximation. `Phi_inf` was taken
  from the `n_TE = 10,000` plateau, so that row anchors the fit rather than
  testing it; the rows at or below 5,000 are what validate both routes.

Because the floor depends on spectrum shape and callable sample size, the
magnitudes above should be re-measured on production data before being quoted for
it. The scaling law, both correction routes, and the recommended reporting are
expected to carry over unchanged.
