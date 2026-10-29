"""Draw random SNP sets from the control universe and score their age CDFs.

This is the unmatched baseline for tools/plot_match_sets.py. Each set is a
uniform draw, ignoring age, of the target's size from the same candidate
universe the matcher used, so it shows what a control set looks like with no
age matching at all.

Reads the interval store, a target built by normalize_tes.te_age_target, and
(optionally) the candidate-rows file passed to bootstrap_target_matcher.py;
without --candidate-rows every eligible store row not in the target is used,
as in the matcher. Writes to --output:

  random_cdfs.npy   (n_sets, n_age_bins) age CDF of each random set
  random_w1.npy     (n_sets,) W1 of each set's CDF to the observed target CDF
  random_seeds.npy  (n_sets,) seed of each draw
  age_bins.npy      the target's age grid
  metadata.json     inputs and settings

The W1 is to the observed target, not to a bootstrap target: nothing is
matched here, so there is no bootstrap target to compare with.

Usage (the store scan is heavy; run it through hpc_run):
    python -m tools.random_control_sets --store STORE --target TARGET \\
        --candidate-rows CANDIDATES.npy --output random_baseline [--n-sets 200]
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from normalize_tes.bootstrap_target_matcher import (
    _candidate_rows, _validate_inputs, derive_seed, target_digest_for,
)
from normalize_tes.release_provenance import software_provenance
from normalize_tes.sample_age_matched_controls import target_seed_identity
from normalize_tes.snp_age_store import is_interval_store, open_snp_age_store
from normalize_tes.swap_control_sampler import aggregate_cdf, analysis_points
from normalize_tes.te_age_target import wasserstein_1


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--store", type=Path, required=True)
    p.add_argument("--target", type=Path, required=True)
    p.add_argument("--candidate-rows", type=Path,
                   help="candidate rows file used for matching (default: all "
                        "eligible rows outside the target)")
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--n-sets", type=int, default=200)
    p.add_argument("--seed", type=int, default=1002)
    p.add_argument("-A", "--a-type", choices=("TE", "SNP"), default="TE")
    args = p.parse_args()
    if args.n_sets < 1:
        p.error("--n-sets must be positive")
    if args.output.exists():
        raise FileExistsError(f"output already exists: {args.output}")

    store = open_snp_age_store(args.store)
    if not is_interval_store(store):
        raise ValueError("random baseline requires an interval store")
    _, target_rows, observed, age_bins, _, target_meta = target_digest_for(args.target)
    _validate_inputs(store, target_meta, args.a_type)
    candidates, _, _ = _candidate_rows(args, store, target_rows)
    if candidates.size <= target_rows.size:
        raise ValueError("candidate universe must exceed target set size")
    points = analysis_points(age_bins)
    seed_identity = target_seed_identity(target_rows, target_meta)

    cdfs = np.empty((args.n_sets, age_bins.size), dtype=np.float32)
    w1 = np.empty(args.n_sets, dtype=np.float64)
    seeds = np.empty(args.n_sets, dtype=np.uint64)
    started = time.perf_counter()
    for i in range(args.n_sets):
        seeds[i] = derive_seed(args.seed, f"{seed_identity}\0random-baseline", i)
        rng = np.random.default_rng(int(seeds[i]))
        rows = np.sort(rng.choice(candidates, size=target_rows.size, replace=False))
        cdf = aggregate_cdf(store, rows, points)
        cdfs[i] = cdf
        w1[i] = wasserstein_1(cdf, observed, age_bins)
        if (i + 1) % 10 == 0 or i + 1 == args.n_sets:
            print(f"set={i + 1}/{args.n_sets} elapsed={time.perf_counter() - started:.0f}s "
                  f"median_w1={np.median(w1[:i + 1]):.6g}", flush=True)

    args.output.mkdir(parents=True)
    np.save(args.output / "random_cdfs.npy", cdfs)
    np.save(args.output / "random_w1.npy", w1)
    np.save(args.output / "random_seeds.npy", seeds)
    np.save(args.output / "age_bins.npy", np.asarray(age_bins, dtype=np.float64))
    (args.output / "metadata.json").write_text(json.dumps({
        "n_sets": args.n_sets,
        "set_size": int(target_rows.size),
        "candidate_universe_size": int(candidates.size),
        "store": str(args.store.resolve()),
        "target": str(args.target.resolve()),
        "candidate_rows": (
            str(args.candidate_rows.resolve()) if args.candidate_rows else None),
        "global_seed": args.seed,
        "a_type": args.a_type,
        "draw": "uniform without replacement from the candidate universe, ignoring age",
        "w1_reference": "observed target CDF",
        "software": software_provenance(),
    }, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(args.output)


if __name__ == "__main__":
    main()
