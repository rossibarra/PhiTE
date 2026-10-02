# Phi-SFS methods

This document defines the statistic PhiTE reports, how sites are polarized, and
how the null is calibrated. The operator guide is the [README](../README.md); the
evidence that the design is calibrated is in [VALIDATION.md](VALIDATION.md); the
implementation design is in
[PHI_SFS_WASSERSTEIN_CODING_PLAN.md](PHI_SFS_WASSERSTEIN_CODING_PLAN.md).

Notation: uppercase $M$ is the number of sites in each compared set; lowercase
$m=20$ is the number of individuals in the SFS projection. $R$ is the number of
null sets.

## The statistic

Phi-SFS compares the normalized cumulative unfolded SFS of a focal set $A$ with
that of an age-matched neutral SNP set, $B_0$. Let $F_A(x)$ and $F_{B_0}(x)$ be
their CDFs on the derived-allele-frequency (DAF) axis. Then

$$
\Phi_{\mathrm{SFS}}(A,B_0)
= W_1(A,B_0)
= \int_0^1 \left|F_A(x)-F_{B_0}(x)\right|\,dx.
$$

For equally spaced projected DAF bins $x_j=j/m$, this is exactly

$$
\Phi_{\mathrm{SFS}}(A,B_0)
= \frac{1}{m}\sum_{j=1}^{m-1}
  \left|F_A(x_j)-F_{B_0}(x_j)\right|.
$$

Phi-SFS is the area between the two CDFs. It is zero only when the spectra are
identical, and grows as probability mass must move farther along the DAF axis. It
is unsigned: the CDFs and signed bin residuals show whether the focal set has an
excess of rare or high-frequency derived alleles.

![Schematic definition of Phi-SFS as the area between focal and neutral SFS cumulative distribution functions](../figures/phi_sfs_definition_schematic.png)

## Matched controls

Each control set independently bootstraps the focal set's site ages (iid, from all
usable posterior draws) and matches exactly $M$ SNPs to that bootstrap age CDF.
The matcher never looks at allele frequency. Production publishes 500 sets in
disjoint mode: once a SNP is used it is removed from the pool, so no SNP appears
in two sets.

Disjoint matching draws later sets from a depleted pool, so the sets are not
identically distributed. In a 1,001-set run on the in-gene category
($M=4{,}067$), sets after about 800 drifted away from the earlier ones; the first
500 did not. Production therefore uses 500 sets, and every category must pass the
drift check in the README before its result is used (see
[VALIDATION.md](VALIDATION.md#v4-and-v5-production-matching-and-phi-sfs)). A
larger category depletes the pool faster, so its drift may start earlier.

TE targets are built exactly as SNP targets are: ages from all usable posterior
draws, no filter on ARG polarity, and the same eligibility rows. The matcher and
Phi-SFS refuse a target built with a TE polarity mask.

## Polarity: one posterior rule for every site

Every site, TE or SNP, focal or control, is polarized the same way. Let $p$ be the
observed ALT frequency and

$$
q=P(\mathrm{ALT\ is\ derived}\mid\mathrm{usable\ ARG\ draws}),
$$

the fraction of usable ARG draws in which ALT is derived. The site contributes the
posterior mixture of its two orientations,

$$
q\,h(k,n)+(1-q)\,h(n-k,n).
$$

All $q\in[0,1]$ are kept; there is no polarity filter. ARG draws that cannot
orient either allele are reported as unusable and count toward neither direction.
A TE is an ACGT-coded site in the ARGs like any other, so the same ancestral table
gives its $q$. Where the ARG calls absence derived, the TE contributes at the
absence frequency, just as a SNP whose REF is derived contributes at $1-p$.

**Why TEs are not given their known polarity.** TE presence is biologically the
derived state for a true single insertion, and earlier designs used that: TEs got
a hard insertion-derived spectrum, TEs with under 50% derived support were
dropped, and TE ages came from agreeing draws only. That made the observed
comparison differ in construction from every null comparison. A had true hard
polarity while the nulls carried the ARG's posterior polarity, so any
miscalibration of $q$ entered $\Phi_{\mathrm{obs}}$ and no $\Phi_i^0$
([review 12](CODE_REVIEW_ROUND12.md), findings 1 and 3). In neutral simulations
that construction rejected far above the nominal rate
([VALIDATION.md](VALIDATION.md#polarity-construction)). The cost of the current
rule is that biological knowledge of TE polarity is not used. The test asks
whether TE sites differ from age-matched SNPs as both are seen through the ARG,
not what the TE's true spectrum is, and signal is probably attenuated where the
ARG mis-polarizes TEs.

**Default: posterior mixture throughout.** A, $B_0$ and every null set $B_i$ are
posterior mixtures:

$$
\Phi_{\mathrm{obs}}=\Phi_{\mathrm{SFS}}(A_{\mathrm{mix}},B_{0,\mathrm{mix}}),
\qquad
\Phi_i^0=\Phi_{\mathrm{SFS}}(B_{i,\mathrm{mix}},B_{0,\mathrm{mix}}).
$$

Because $E[h(\mathrm{DAF})]=q\,h(k,n)+(1-q)\,h(n-k,n)$, a Bernoulli-$q$ hard
spectrum is the mixture plus imputation noise that carries no information about
the data. The mixture avoids that noise and does not depend on an imputation seed.

**Option: Bernoulli-$q$ hard orientation** (`--asymmetric-polarity-null`). A and
every $B_i$ are hard-oriented by one Bernoulli-$q$ draw per site, and $B_0$ stays a
mixture:
$\Phi_i^0=\Phi_{\mathrm{SFS}}(B_{i,\mathrm{Bernoulli}(q)},B_{0,\mathrm{mix}})$. ALT
is declared derived when a reproducible coordinate-keyed uniform variate is less
than $q$, keyed by `sha256(phi-sfs-bernoulli-q-v1, seed, chromosome, position)`,
so a site gets the same orientation wherever it is used. Results then condition on
one `--polarity-imputation-seed`.

## Null calibration

Finite site counts make the distance positive even under neutrality, so the
sampling floor is calibrated separately for every focal category rather than
comparing raw Phi-SFS values across categories.

1. Apply matching QC. Draw $B_0$ uniformly from the QC-passing sets, with a seed
   derived from `--reference-seed` and the target digest. Every other QC-passing
   set is a null, so $R$ is whatever passes QC; production requires $R\ge450$.
2. $\Phi_{\mathrm{obs}}=\Phi_{\mathrm{SFS}}(A,B_0)$.
3. $\Phi_i^0=\Phi_{\mathrm{SFS}}(B_i,B_0)$ for $i=1,\ldots,R$. Each $\Phi_i^0$ is
   a raw distance between two neutral SNP sets, not a Z-score.
4. With $\mu_0$ and $s_0$ the mean and sample standard deviation of the
   $\Phi_i^0$,

   $$
   Z_A=\frac{\Phi_{\mathrm{obs}}-\mu_0}{s_0}.
   $$

5. The one-sided Monte Carlo P-value is

   $$
   P_A=\frac{1+\sum_{i=1}^{R}
   \mathbf{1}\!\left(\Phi_i^0\ge \Phi_{\mathrm{obs}}\right)}{R+1}.
   $$

The add-one P-value is valid for any $R$. With $R\approx500$ the smallest
attainable value is about $2\times10^{-3}$.

$B_0$ is drawn before any SFS is examined, so the choice of reference is SFS-blind.
Replicate 0 is not the default reference, because it is matched first, from the
undepleted pool, and so is not a typical set. `--reference-sensitivity N` repeats
the calibration with the next $N$ sets of the same seeded permutation as $B_0$,
which shows how much the result depends on that one draw.

`--max-null-replicates N` instead fixes $R=N$: $B_0$ is the first QC-passing set of
the seeded permutation and the nulls are the next $N$. The negative control uses
$N=99$, which gives an exact P-value grid of 0.01.

## Effect size, P and Z

**Effect size: $\Phi_{\mathrm{obs}}-\mu_0$**, in DAF units. Raw
$\Phi_{\mathrm{obs}}$ is not an effect size on its own: two finite sets drawn from
the same spectrum still have $\Phi>0$, and that floor, $\mu_0$, is larger for
smaller $M$. A small category can therefore have the largest raw
$\Phi_{\mathrm{obs}}$ and the smallest departure. Subtracting the category's own
null mean removes the floor. It is an approximate correction: under an
alternative, the observed distance and the floor need not simply add, so the
excess is not necessarily the distance between the true spectra. Report it with $M$
and the signed CDF and bin residuals, which give the direction of the shift.

**P** tests whether the focal spectrum is farther from its matched neutral
background than two finite neutral samples of the same size would be. It is the
significance of the departure, not its size.

**$Z_A$** expresses the departure in null standard deviations. Because $s_0$
shrinks as $M$ grows, $Z_A$ grows with $M$ for a fixed spectral difference, so it
measures test strength, not effect size. It is kept in `summary.csv` but should
not be plotted or compared across categories.

To plot many categories, use the $\Phi$ scale. For each category, draw its null
distribution of $\Phi_i^0$ in grey, mark $\mu_0$, and draw $\Phi_{\mathrm{obs}}$
as a point coloured by $-\log_{10}P_A$, with a segment from $\mu_0$ to
$\Phi_{\mathrm{obs}}$ for the effect size. Cap the colour scale at
$\log_{10}(R+1)$, about 2.7 for $R\approx500$.

![Illustrative Phi-SFS null distributions by category, with observed values, null means, effect sizes and P-value colours](../figures/phi_sfs_null_example.png)

In this synthetic example, the >5 kb category has the highest raw
$\Phi_{\mathrm{obs}}$ but the smallest effect, and is not significant, because
its small $M$ gives it the highest floor.

`docs/PHI_SFS_SAMPLE_SIZE_BIAS.md` gives a floor correction,
$\sqrt{\Phi_{\mathrm{obs}}^2-\Phi_{\mathrm{floor}}^2}$, that recovers the true
distance more closely. It was validated for the earlier total-variation
statistic, not for the current Wasserstein distance, so it is not used here until
it has been checked for $W_1$.

## Assumptions and limits

- **Exchangeability.** Drawing $B_0$ at random makes the choice of reference
  uniform among accepted sets. It does not make depleted sets identically
  distributed, and does not show that the $A$-versus-$B_0$ comparison is
  exchangeable with the $B_i$-versus-$B_0$ comparisons
  ([review 12](CODE_REVIEW_ROUND12.md), finding 7). The negative control, the
  drift check and reference sensitivity together address this
  ([VALIDATION.md](VALIDATION.md)).
- **QC selection.** Matching QC uses only age matching, but age and allele
  frequency are related, so selecting sets on QC is not guaranteed to be neutral
  with respect to their spectra.
- **Linkage.** The focal-age bootstrap resamples sites iid. This adopts the
  Poisson-random-field approximation that local LD averages out for genome-wide
  control pools of millions of SNPs. Revisit it for spatially restricted or
  strongly clustered sets.
- **A is observed once.** The focal set is fixed, so small or unusual focal sets
  can give unstable results even after calibration. Always report $M$, the raw
  distance, null mean and SD, Z, P, $R$ and the matching diagnostics.
- **Between-category contrasts are not supported.** `normalize_tes.phi_contrast`
  is experimental: its random pairing of category nulls is not validated and
  ignores covariance between nested or overlapping categories. Report each
  category's result separately; a visual difference between two points is not a
  formal test.
- **Not tested by any validation:** TE genotyping error, recurrent insertion or
  deletion, selection in real TE categories, ARG failure that differs between TE
  and SNP sites, and demography outside the simulated models.

Results use the `phi-sfs-wasserstein-v3` schema and cannot be combined with older
Phi-SFS outputs (see [OPTIONS.md](OPTIONS.md#output-schemas)).
