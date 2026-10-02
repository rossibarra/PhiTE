# Validation

This document collects the evidence that PhiTE's Phi-SFS test is calibrated. The
statistic is defined in [METHODS.md](METHODS.md). The prespecified criteria, and the
two amendments made after results were seen, are in
[REMAINING_VALIDATION_PROPOSAL.md](REMAINING_VALIDATION_PROPOSAL.md) (items V1–V11).
Production settings and measured resources are in
[BOOTSTRAP_HPC_VALIDATION.md](BOOTSTRAP_HPC_VALIDATION.md).

## Status for v0.9.0

| item | what it tests | result |
|---|---|---|
| V1 | all-mixture design on simulations where the former design failed | **pass**: 20/20 applicable cells |
| V2 | how far disjoint matching can go before age bins run out | done; reproduces the earlier capacity table |
| V3 | power of the all-mixture design against alternatives | not run (non-blocking) |
| V4 | production matching under depletion | **pass on 500 sets** (Amendment A); the 1,001-set run failed drift |
| V5 | integrity of production Phi-SFS | **pass** |
| V6 | real-data negative control through the production matcher | **pass**: 11/300 rejections |
| V7 | exchangeability | covered by V4, V5 and V6 |
| V8 | between-category contrasts | out of scope; `phi_contrast` stays experimental |
| V9 | release housekeeping | **pass**; tests pass, and the 15 dirty-checkout V6 tests were rerun (Amendment B) |
| V10 | ancestral-table digests, fixed-$R$ option | done |
| V11 | documented limits | see [METHODS.md](METHODS.md#assumptions-and-limits) |

## Polarity construction

The dnAging simulations are neutral throughout, so any rejection rate above
$\alpha$ is miscalibration. With nulls hard-oriented by Bernoulli-$q$ and $B_0$ a
posterior mixture, giving the focal set its true polarity rejected far above the
nominal rate once the ARG was inferred without ancestral information. Rejection
fractions at $\alpha=0.05$, 150 tests per cell, from
`results/sim_dnaging_unpolarised/arms_v1/summary.csv`:

| focal construction | $M=250$ | $M=1000$ | $M=4000$ |
|---|---|---|---|
| true polarity everywhere (oracle) | 0.09 | 0.05 | 0.09 |
| A true-hard, no filter, true ages | 0.24 | 0.44 | 0.83 |
| A true-hard, filtered, agreeing-draw ages, inferred ages (former production) | 0.27 | 0.75 | 0.97 |
| A Bernoulli-$q$, like the nulls | 0.09 | 0.10 | 0.07 |

The reference-haplotype replicates (`results/sim_dnaging_refhap/arms_v1/`) show the
same pattern. Treating A like the nulls removes the excess, which is why TEs are
now polarized as SNPs are.

## All-mixture design in simulation (V1, V2)

`tools/sim_polarity_arms.py` tests the all-mixture design directly. For both true
and inferred ages it compares independently reusable controls with sequentially
depleted, globally disjoint ones. Every control set gets its own iid bootstrap of
the focal ages; the depleted arm removes selected controls before matching the
next set and draws $B_0$ uniformly from the completed sets. The depleted sampler is
a binned stand-in, not the production optimizer. Simulation replicates are the
independent units; repeated focal draws within a replicate add Monte Carlo
precision but are not independent. Tests with fewer than 19 nulls cannot reach
$p\le0.05$ and are excluded from rejection rates.

**Original dnAging replicates** (10 replicates, where the former design also
passed): the all-mixture arms rejected 0.8–4.0% of tests at $\alpha=0.05$. That is
valid but conservative, and it held under depletion wherever tests had at least 19
nulls.

**Unpolarised and reference-haplotype replicates (V1)**, where the former design
failed. On the same 10 paired replicates the simulator's `production` arm (A true,
$B_0$ a mixture, Bernoulli-$q$ nulls) rejected 14–83% of tests. The four
all-mixture arms, with equal replicate weights:

| condition | rejection rate per cell | one-sided 95% upper bound |
|---|---|---|
| unpolarised | 0.8–4.6% | ≤ 0.060 |
| reference haplotype | 0.6–4.2% | ≤ 0.058 |

The criterion was a rate of at most 0.075 and an upper bound (replicate bootstrap)
of at most 0.10 in every cell. The depleted arms could not reach 19 nulls at
$M=4000$ in any replicate, as the V2 capacity report predicted, so those four
cells are untested. Report: `results/v1_report`.

To reproduce, on a compute node (each condition took about 40 minutes; the
simulator is not checkpointed, so a preempted job must be restarted with a new
output directory):

```bash
for c in unpolarised refhap; do
  python -m tools.sim_polarity_arms \
    --prep-root results/sim_dnaging_$c/prep \
    --output results/sim_dnaging_$c/arms_v1_rc2 \
    --replicates 1 2 3 4 5 6 7 8 9 10
done
python -m tools.v1_report \
  --condition unpolarised=results/sim_dnaging_unpolarised/arms_v1_rc2 \
  --condition refhap=results/sim_dnaging_refhap/arms_v1_rc2 \
  --output results/v1_report
```

The original-replicate run used
`--prep-root results/sim_dnaging --sizes 250 1000 4000 --tests 50 --nulls 199
--mixture-sets 80`, which are the defaults.

## V4 and V5: production matching and Phi-SFS

**V4** matched 1,001 disjoint sets for the in-gene TE category ($M=4{,}067$, 3
restarts, seed 1002). It passed every matching criterion: maximum reuse 1, no
duplicates, 996/1,001 QC passes, at least 95 in every block of 100, and a 5-point
fall in pass rate. It **failed** SFS drift: each set's distance from the pooled
spectrum correlated with matching order ($\rho=0.164$), and the first and last
quarters differed by $\mathrm{SMD}=0.48$ (limits 0.1 and 0.2). The distance was
flat through about set 800 and rose after it. Median matching error showed the
same jump, and all five QC failures were in the last block.

Under Amendment A, V4 and V5 are judged on the first 500 sets. Seeds depend only
on the global seed, target, replicate and restart, and depletion is sequential, so
these are the sets a 500-set run produces. On those sets drift passes
($\rho=0.027$, $\mathrm{SMD}=0.014$), all 500 pass QC, and maximum reuse is 1. This
change was made after the result was seen; the 1,001-set failure stays on record.
Reports: `results/v4/in_gene_rc2` (1,001 sets) and
`results/v4/in_gene_rc2_first500`.

**V5** ran Phi-SFS on the 500 sets: $R=499$, every set passes QC, reuse 1, no
overlap with $B_0$, $M$ identical across A and every set, finite null mean and
nonzero SD, a recomputed add-one P-value that agrees with the output, and matching
input digests. Reference sensitivity with 10 alternative references gave Z from
9.3 to 15.8; every P-value was at its floor of 0.002. Output:
`results/phi_sfs/in_gene_rc2_first500`.

## Real-data negative control (V6)

300 independent SNP-versus-SNP tests through the production matcher and the
all-mixture design. Each test sampled its own focal set of 4,000 SNPs from the
genome-wide pool, excluded those sites from its own controls, matched 110 disjoint
sets with 3 restarts and its own seed, and ran Phi-SFS with exactly $R=99$ nulls.

- All 300 tests produced 99 valid nulls; every bundle passed QC on all 110 sets.
- **11 of 300 rejected at $\alpha=0.05$ (0.037).** The one-sided 95%
  Clopper–Pearson upper bound is 0.060. The prespecified limit was at most 21
  rejections, the largest count whose bound stays below 0.10.
- P-values were close to uniform: decile counts 23–38 against 30 expected.
- Rejections were 1–3 in each block of 50 tests, and 2–5% by the position of
  $B_0$ in matcher order.
- Focal sets were not forced to be disjoint; two sets shared at most 7 of their
  4,000 sites.

The result is conditional on this genome, this candidate pool, and $M=4{,}000$,
and uses SNP focal sets only. Report: `results/v6_report`, from
`tools/v6_report.py`.

**Earlier pilot (superseded).** 100 held-out real SNP sets, Bernoulli-$q$ focal
sets against a mixture reference with 344 Bernoulli-$q$ nulls, rejected 7 of 100
(Wilson 95% interval 0.034–0.137). All 100 tests shared one reference and one null
vector, so the P-values were dependent, and none of the sets met the production
matching threshold. Provenance: `results/phi_sfs/snp_type1_asymmetric_100/summary.json`.

## Provenance (V9, V10)

Ancestral tables record a SHA-256 digest of each array, and every reader verifies
it. No validation output records the `v0.9.0-rc2` commit itself: the runs span
later commits that change only `slurm/` and `tools/`. Under Amendment B, an output
counts if its commit's `normalize_tes/` is identical to rc2's.

Fifteen V6 tests had run on a checkout with unrecorded uncommitted edits. They
were rerun on a clean checkout with the same seeds (`results/v6_clean_rerun`):

- 10 reproduced bit for bit, and test 163 reproduced its $B_0$, nulls and P, with
  the null mean and SD differing only in the last digits.
- Tests 3, 5, 245 and 279 drew a different $B_0$. Their focal sets, candidate rows
  and matching QC were identical, but the age target's bootstrap distances
  differed by about 1 part in $10^5$. The target digest hashes those values, and
  the $B_0$ draw is seeded from the digest, so the difference changed $B_0$.
- A follow-up built test 279's target repeatedly on one node
  (`results/target_determinism_39392444`). Builds with the same thread count were
  bit-identical on that node and on a second node of the same kind; changing the
  thread count (1, 2 or 6) or running on a GPU node changed the values. Tests 3
  and 5 originally ran with 6 CPUs and tests 245 and 279 on a GPU node, which
  accounts for all four. The uncommitted edits were not the cause.
- Test 279's P moved from 0.04 to 0.61. V6 counts the original runs (11/300
  rejections); with the reruns it would be 10/300. Both pass.

**Known limitation.** Because $B_0$ is seeded from a digest of floating-point
target arrays, rebuilding a target with a different thread count or CPU type can
select a different $B_0$ and change P. Rerunning from the published target is
reproducible. The fix, seeding $B_0$ from the target's inputs instead, is planned
for the next release.
