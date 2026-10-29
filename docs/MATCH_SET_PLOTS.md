# Match-set plots

Diagnostic figures for one matcher run (`normalize_tes.bootstrap_target_matcher`):
how every tested set, the accepted final matches and the observed TEs compare,
optionally against an unmatched random-SNP baseline. They are for inspection,
not part of the pipeline's outputs.

## What the figures show

Both figures use the same groups.

| Group | Age CDF figure | W1 figure |
|---|---|---|
| tested (grey) | best state of every restart of every replicate (`restart_best_cdfs.npy`) | W1 of each of those to its bootstrap target (`restart_match_to_bootstrap_w1.npy`) |
| accepted (blue) | final match of each QC-passing replicate (`cdfs.npy`, `qc_pass.npy`) | its W1 to the bootstrap target (`match_to_bootstrap_w1.npy`) |
| observed (vermillion) | the observed TE age CDF (`target_cdf.npy`) | `bootstrap_to_observed_w1.npy`: how far each bootstrap target lies from the observed TEs |
| random SNPs (green, optional) | age CDF of random sets | W1 of random sets to the observed CDF |

Replicates that fail QC appear only in the grey group.

Matches are made to a bootstrap target, not to the observed CDF, so no match has
a W1 to compare with an "observed W1" directly. In the W1 figure the vermillion
curve is the sampling scale instead: how far resampling alone moves the target.
Matched sets sit well below it when matching worked. `match_to_observed_w1.npy`
holds each match's distance to the observed CDF if you want that.

The age CDF figure has a second panel of each CDF minus the observed CDF,
because the curves are otherwise indistinguishable. With a random baseline it
has a third panel at full scale, since the random sets are far larger than the
matched differences and the second panel is zoomed to the matched ones.

## Environment

`matplotlib` is in `environment.yml`. Update the env once:
`conda env update -f environment.yml`.

Run everything through `~/.claude/bin/hpc_run` (or your own `srun`/`sbatch`),
not on a login node.

## 1. Plot a matcher run

```
~/.claude/bin/hpc_run 'python tools/plot_match_sets.py MATCH_DIR --name NAME'
```

`MATCH_DIR` is a matcher output directory, the product of one data store and one
TE target. The script reads its saved arrays; it does not open the store.

Writes `NAME_age_cdf.{png,pdf}` and `NAME_w1.{png,pdf}` to `--out-dir`
(default `figures/match_sets`). `--name` defaults to the directory name.
`--linear-age` uses a linear age axis; `--age-max G` drops bins above `G`
generations.

## 2. Add a random-SNP baseline for your data

`tools/random_control_sets.py` draws random SNP sets of the target's size,
uniformly and ignoring age, from the same candidate universe the matcher used,
and scores their age CDFs. It needs the interval store, so it is a heavy job.

```
~/.claude/bin/hpc_run 'python -m tools.random_control_sets \
    --store STORE --target TARGET --candidate-rows CANDIDATES.npy \
    --output random_baseline --n-sets 200'
~/.claude/bin/hpc_run 'python tools/plot_match_sets.py MATCH_DIR --random-dir random_baseline'
```

Use the same `--store`, `--target` and `--candidate-rows` as the matcher run.
Without `--candidate-rows`, every eligible store row outside the target is used,
as in the matcher. `-A` sets the A type (default `TE`) and must agree with the
target. `--seed` seeds the draws.

Output in `random_baseline/`: `random_cdfs.npy`, `random_w1.npy`,
`random_seeds.npy`, `age_bins.npy`, `metadata.json`. `plot_match_sets.py`
refuses a baseline whose age grid differs from the match's.

The random sets' W1 is to the observed CDF. Nothing is matched, so there is no
bootstrap target.

If you only have W1 values, `--random-w1 FILE.npy` adds them to the W1 figure
alone (no age CDFs), and `--random-label TEXT` sets the legend text.

## Examples

`figures/match_sets/` has figures for `results/bootstrap_matches/in_gene_75draw_restarts6`.

- `*_age_cdf.*`, `*_w1.*`: matched sets and observed only.
- `*_synthetic_random_*`: the same plus an **invented** random baseline from
  `tools/synthetic_random_baseline.py`. Its age curve shape and spread are made
  up; only its location is tuned so the median W1 is about 8.8e5 generations, the
  median seen in `results/init_compare/random`, which used a different target.
  Treat it as an illustration of the layout, not as a measurement.
- `*_with_random_w1.*`: the W1 figure with the real `init_compare/random` initial
  W1 values, labelled as a different target.

## Not yet tested

`random_control_sets.py` has been checked to import and print `--help`, and its
output layout is what `plot_match_sets.py` reads (tested with the synthetic
data). It has not been run against a real store.
