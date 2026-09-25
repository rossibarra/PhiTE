# Simulation plan for validating PhiTE

## Objective

Use forward simulations to test whether age-conditioned Phi-SFS distinguishes
selection on transposable-element (TE) insertions from changes in the rate at
which new insertions arise.

The key comparison is between:

1. frequency spectra computed against randomly sampled neutral SNPs;
2. frequency spectra computed against SNPs matched to the true TE ages; and
3. the production PhiTE workflow using genealogy-based posterior age
   distributions and bootstrap-matched SNP controls.

The simulations should demonstrate two properties:

- neutral TEs with nonconstant insertion rates can produce misleading
  unmatched Phi-SFS values; and
- conditioning the SNP controls on TE age reduces those false signals while
  retaining power to detect selection.

Phi-SFS is a symmetric distance, so the signed bin residuals must accompany it
to distinguish an excess of low-frequency TE insertions from an excess of
high-frequency insertions.

## TE model

The initial implementation is `simulations/te_forward.slim`. It is based on
the conceptual TE model in section 14.10 of the SLiM manual and the simulations
of Horvath et al. (2022), but makes the transposition mechanism, selection
coefficient, and rate history command-line parameters.

### Insertion identity and state

Every insertion locus receives a unique integer `te_id`. An insertion can be
in either of two allelic states on a haplotype:

- `m2`: functional and able to transpose;
- `m3`: disabled and unable to transpose.

Disabling changes `m2` to `m3` at the same position while preserving the
insertion ID, original insertion tick, parent insertion ID, and selection
coefficient. Both states therefore represent the same biological insertion.

At the end of a run, the two states are collapsed into one insertion-presence
allele:

```text
inserted = functional OR disabled
```

The population output reports functional, disabled, and total inserted copy
counts separately. The sampled presence matrix contains only the merged binary
insertion state. This matrix, rather than the raw `m2`/`m3` mutation records,
is the authoritative TE genotype source for downstream VCF construction.

The raw tree sequence is intentionally not yet a PhiTE input. It can
contain both functional and disabled mutation records at one insertion site,
and cut-and-paste movement can introduce loss events at donor sites. A later
canonicalization step will build the biallelic insertion VCF and either infer
posterior ARGs from it or construct a validated single-mutation representation.

Insertion IDs, not coordinates, are the permanent event identities. A new
destination must differ from every currently extant TE locus and every other
destination allocated in the same tick, but a coordinate may be reused after
the previous insertion there has been lost from the population. The allocator
uses a per-tick occupancy bitmap, draws destinations in batches, and stops with
an explicit error if its configurable attempt cap is reached. It never retains
or linearly scans a lifetime history of insertion positions.

### Copy-and-paste

When a functional source copy transposes:

1. the source remains present on the transmitting haplotype;
2. a globally unused destination coordinate is chosen;
3. a new functional insertion is created at the destination; and
4. the destination receives a new insertion ID whose `parent_te_id` is the
   source ID.

Copy number can increase under this model.

### Cut-and-paste

When a functional source copy transposes:

1. a globally unused destination coordinate is chosen;
2. a new functional insertion with a new insertion ID is created there; and
3. the source insertion is removed from the transmitting haplotype.

The movement event itself conserves copy number. The donor insertion can remain
segregating on other haplotypes. Disabled copies do not move in either model.

This is an insertion-locus model, not a nucleotide-level excision model. The
first analysis should treat source loss and destination gain as the relevant
presence/absence changes and should not attempt to model an excision footprint.

## Initial parameterization

The defaults reproduce the broad scale of the published TE simulations:

| Parameter | Default |
|---|---:|
| Diploid population size | 500 |
| Chromosome length | 3 Mb |
| Burn-in | 5,000 ticks |
| Rate-change interval | 250 ticks |
| Initial fixed functional TEs | 5,000 |
| Neutral SNP overlay rate | `1e-8` per site per generation |
| Recombination rate | `1e-8` per site per generation |
| Baseline TE movement probability | `1e-4` per functional copy per tick |
| TE disabling probability | `5e-5` per functional copy per generation |
| TE dominance coefficient | `0.5` |
| Sample | one haplotype from each of 25 individuals |

SLiM records ancestry and TE dynamics but has a forward neutral mutation rate
of zero. After the run, the ancestry is recapitated and neutral mutations are
overlaid with msprime. The core design uses `r / mu = 1`; the previous
metacentric map with `r / mu` as high as 100 would leave too few mutations per
tree for a fair evaluation of ARG inference.

The core control pool is generated on a separate, unlinked neutral sequence
under the same demographic history. Neutral marker mutations may also be
overlaid on the TE-bearing sequence to provide local information for ARG
inference, but those linked markers are not the primary neutral control pool.
This separates the control SFS from linked selection without making the
TE-bearing genealogy nearly uninformative.

The 5,000-tick burn-in remains a starting value for TE dynamics, not a device
for accumulating neutral SNPs. Pilot equilibrium diagnostics should determine
whether it can be shortened. Recapitation supplies ancestral genealogy before
the forward simulation; it does not by itself replace the time needed for TE
copy-number and age distributions to approach their intended state.

The default dominance coefficient is 0.5 and is configurable with
`TE_DOMINANCE`. The model rejects `TE_SELECTION_COEFF <= -1`, but that bound
does not by itself guarantee a comfortable numerical margin for multiplicative
fitness. The default parameter grid remains within double precision; any
increase in founder count or selection magnitude requires a pilot check of the
fitness distribution before the full factorial run.

## Core factorial design

| Factor | Values |
|---|---|
| Movement mechanism | copy-and-paste, cut-and-paste |
| Scaled selection, `4 Ne s` | 0, -2, -10, -20 |
| `s` when `Ne = 500` | 0, -0.001, -0.005, -0.01 |
| Post-burn-in rate multiplier | 1, 10, 0.1 |
| Interpretation | constant, recent increase, recent decrease |

The pilot should use ten independent runs per cell. Runtime, memory, insertion
counts, neutral-SNP counts, and coordinate-collision rates from the pilot will
determine the main replicate count. The intended main analysis is 100--200
independent populations per cell.

Because copy-and-paste and cut-and-paste can produce very different numbers of
segregating insertions, the primary comparisons should use fixed target sizes
such as 100, 250, 500, and 1,000 TEs, whenever a replicate contains enough
sites. Total TE abundance remains a separate biological outcome.

## Simulation outputs

For an output prefix `RUN`, the initial SLiM model writes:

- `RUN.trees`: raw SLiM tree sequence;
- `RUN.parameters.tsv`: complete model parameters and realized seed;
- `RUN.samples.tsv`: sampled individuals and haplotype pedigree IDs;
- `RUN.te_population.tsv`: population-wide functional, disabled, and merged
  insertion counts;
- `RUN.te_presence.tsv`: merged haploid insertion calls for loci polymorphic in
  the sampled haplotypes; and
- `RUN.complete`: written last, only after every other output succeeds.

`RUN.te_population.tsv` retains all extant population loci, including loci
absent or fixed in the sample. `RUN.te_presence.tsv` excludes sample count zero
and sample count `SAMPLE_SIZE`, and records `sample_inserted_count` explicitly.
The SLiM positions in both tables are zero-based; the later VCF exporter must
convert them to one-based coordinates.

The output writer refuses to overwrite a prefix carrying `RUN.complete`.
If a previous attempt was preempted or failed before writing the marker, the
model removes only its known partial files and reruns the exact prefix. Every
`writeFile()` return value is checked. The tree sequence carries the model
parameters as top-level user metadata so it remains self-describing if moved
away from its tabular companions. Each Slurm array task must use its own output
prefix because the project resides on Quobyte, where concurrent writers must
never share a file.

The output prefix's parent directory must already exist. A representative
neutral copy-and-paste run with a recent tenfold rate increase is:

```bash
slim \
  -d "OUT_PREFIX='results/simulations/pilot/copy_up_neutral_001/run'" \
  -s 1001 \
  -d "MOVEMENT_MODE='copy'" \
  -d TE_SELECTION_COEFF=0.0 \
  -d POST_BURNIN_MULTIPLIER=10.0 \
  simulations/te_forward.slim
```

For a cut-and-paste run, set `MOVEMENT_MODE='cut'`. The rate-history settings
are `POST_BURNIN_MULTIPLIER=1.0`, `10.0`, and `0.1` for constant, increased,
and decreased movement rates. `RANDOM_SEED` defaults to SLiM's realized seed,
so `-s` works normally; an explicit `-d RANDOM_SEED=...` intentionally
overrides it. These commands are examples only; full runs must
be launched on a compute node after SLiM has been pinned and the script has
passed the small functional tests below.

## SLiM model validation gates

Before the factorial pilot, run short models with a small population,
chromosome, initial TE count, burn-in, and change interval. A model is ready
for the pilot only if all of these checks pass for both movement modes:

1. the script parses and runs under the pinned SLiM version;
2. rerunning the same seed reproduces byte-identical tabular outputs;
3. every `te_id` maps to exactly one genomic position;
4. functional plus disabled count equals merged inserted count;
5. no sampled haplotype is called inserted twice at one `te_id`;
6. disabling preserves position, insertion tick, parent ID, and fitness effect;
7. copy-and-paste adds one destination without removing its source;
8. cut-and-paste adds one destination and removes its source from the moving
   haplotype;
9. no two extant insertion IDs share a coordinate, same-tick destinations are
   unique, and coordinates of extinct insertions can be reused;
10. constant, increased, and decreased rate schedules activate on the intended
    ticks; and
11. tree-sequence sample identifiers agree with `RUN.samples.tsv`.

The performance gate should also confirm that each haplotype uses one binomial
draw for disabling and one for movement, and that destination allocation is
batched once per tick. Reintroducing per-copy random vectors or a lifetime
position list would make the default model computationally infeasible.

The later recapitation/overlay stage requires `pyslim` in addition to the
repository's existing msprime and tskit dependencies. The later VCF exporter
must remove recurrent or ambiguous SNP coordinates and any SNP coordinate
overlapping a TE insertion. It must retain only biallelic sites that are
polymorphic in the sampled haplotypes.

## Analysis tiers

### Tier 1: oracle ages

Use the recorded insertion ticks and overlaid SNP mutation times directly.
Recapitate the SLiM ancestry, overlay neutral marker mutations on the
TE-bearing sequence, and generate a separate unlinked neutral-control sequence
under the same demography. Do not accumulate neutral SNPs forward in SLiM.

For each independently simulated population, construct 100 control sets by:

1. sampling neutral SNPs without age matching; and
2. matching neutral SNPs to the exact TE-age distribution.

This tier establishes whether age conditioning works independently of ARG
inference and PhiTE's matching optimizer.

### Tier 2: known-tree age intervals

Derive mutation-age intervals from the known simulated genealogy. This tests
the effect of branch-width uncertainty even when the genealogy is correct.
Raw TE state mutations must first be canonicalized to one insertion-presence
mutation per locus; otherwise the raw tree sequence violates PhiTE's
single-event site assumptions.

### Tier 3: posterior ARGs and production PhiTE

For a stratified subset of simulations:

1. export a biallelic VCF containing neutral SNPs and merged TE-presence calls;
2. infer posterior ARG draws using the empirical analysis workflow;
3. build the PhiTE interval store;
4. build TE targets and bootstrap age-matched controls; and
5. calculate Phi-SFS with the production implementation.

This tier quantifies the additional loss from ARG inference, posterior age
uncertainty, and imperfect control matching.

## Primary performance analysis

The independent simulated population is the replicate. The 100 matched-control
sets within a population are dependent and are not a null distribution or 100
independent tests.

For each simulation and control method:

1. retain all Phi-SFS values and matching diagnostics;
2. use the median Phi-SFS as the population-level score;
3. retain signed SFS residuals to establish the direction of the difference;
4. calibrate a 5% rejection threshold from independent neutral,
   constant-rate simulations; and
5. apply that fixed threshold to all other histories and selection strengths.

Report:

- false-positive rate for neutral constant, increased, and decreased insertion
  rates;
- power for each nonzero selection coefficient and rate history;
- ROC curves and area under the curve within and across rate histories;
- sensitivity to TE target size;
- Phi-SFS versus age-match Wasserstein distance; and
- control reuse, overlap, and effective replicate diagnostics.

The expected result is inflated false-positive rates for unmatched controls
after neutral rate changes, approximately calibrated age-matched results for
those same simulations, and increasing age-matched power as selection becomes
more deleterious.

## Follow-up robustness analyses

After the core experiment succeeds, add:

- samples taken during the rate change and 250, 1,000, and 5,000 ticks after
  restoration of the baseline rate;
- a population bottleneck, alone and combined with each rate change;
- reference-genome ascertainment of TE insertions;
- missing genotype calls;
- different disabling rates;
- an excision-footprint model for cut-and-paste TEs;
- alternative recombination maps and linked versus unlinked control pools;
- nearly neutral selection (`4 Ne s = -0.2`); and
- true ages, known-tree intervals, and inferred posterior ages on the same
  simulated populations.

## Farm execution plan

SLiM simulations are compute-node work. After a one-replicate functional test,
run the pilot as a one-CPU Slurm array with one output directory per task. The
`low` partition is appropriate for short, restartable simulations, with the
tradeoff that it is preemptible and restarts a task from the beginning.

The live Farm association query was unavailable while this plan was written.
Existing repository launchers use account `jrigrp`, but the association and QOS
must be checked again before a job script is submitted. No Farm SLiM module was
found, so the executable must be supplied by a pinned Conda environment or an
Apptainer image before validation runs.

## References

- Haller and Messer. SLiM manual, section 14.10, "Modeling transposable
  elements": <https://benhaller.com/slim/SLiM_Manual.pdf>
- Horvath et al. 2022. "Controlling for Variable Transposition Rate with an
  Age-Adjusted Site Frequency Spectrum":
  <https://doi.org/10.1093/gbe/evac016>
- Published simulation source for Horvath et al.:
  <https://github.com/Roberthorv/age-adjusted_site_frequency_spectrum>
