"""Write an INVENTED random-SNP baseline for illustrating plot_match_sets.py.

Not real data. Use tools/random_control_sets.py with a real store for that.

Random SNPs ignore age, so their age distribution is the control universe's,
not the TEs'. The invented universe here is lognormal; its location is tuned so
that the median W1 of the random sets to the observed target equals
--target-w1 (default 8.8e5 generations, the median seen in
results/init_compare/random, a different target). Each set is the empirical
age CDF of --set-size draws, so the spread between sets is real sampling noise
for that set size. Output has the same layout as random_control_sets.py.

Usage:
    python tools/synthetic_random_baseline.py MATCH_DIR --output DIR
        [--n-sets 100] [--set-size N] [--target-w1 8.8e5] [--sigma 1.3]
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def w1(cdf: np.ndarray, target: np.ndarray, ages: np.ndarray) -> float:
    # Same discrete W1 as normalize_tes.te_age_target.wasserstein_1.
    return float(np.sum(np.abs(cdf[:-1] - target[:-1]) * np.diff(ages)))


def draw_cdfs(mu, sigma, n_sets, set_size, ages, seed):
    rng = np.random.default_rng(seed)
    samples = np.sort(rng.lognormal(mu, sigma, size=(n_sets, set_size)), axis=1)
    return np.stack([
        np.searchsorted(row, ages, side="right") / set_size for row in samples
    ])


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("match_dir", type=Path)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--n-sets", type=int, default=100)
    p.add_argument("--set-size", type=int)
    p.add_argument("--target-w1", type=float, default=8.8e5)
    p.add_argument("--sigma", type=float, default=1.3)
    p.add_argument("--seed", type=int, default=1)
    args = p.parse_args()
    if args.output.exists():
        raise FileExistsError(f"output already exists: {args.output}")

    ages = np.load(args.match_dir / "age_bins.npy")
    observed = np.load(args.match_dir / "target_cdf.npy")
    set_size = args.set_size or int(
        json.loads((args.match_dir / "metadata.json").read_text())["set_size"])

    def median_w1(mu: float) -> float:
        cdfs = draw_cdfs(mu, args.sigma, 20, set_size, ages, args.seed)
        return float(np.median([w1(c, observed, ages) for c in cdfs]))

    lo, hi = np.log(1e3), np.log(1e8)          # W1 rises with mu
    for _ in range(40):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if median_w1(mid) < args.target_w1 else (lo, mid)
    mu = (lo + hi) / 2
    cdfs = draw_cdfs(mu, args.sigma, args.n_sets, set_size, ages, args.seed + 1)
    dist = np.array([w1(c, observed, ages) for c in cdfs])

    args.output.mkdir(parents=True)
    np.save(args.output / "random_cdfs.npy", cdfs.astype(np.float32))
    np.save(args.output / "random_w1.npy", dist)
    np.save(args.output / "age_bins.npy", ages)
    (args.output / "metadata.json").write_text(json.dumps({
        "SYNTHETIC": "invented illustrative data, not computed from a store",
        "n_sets": args.n_sets, "set_size": set_size,
        "lognormal_mu": mu, "lognormal_sigma": args.sigma,
        "target_median_w1": args.target_w1,
        "achieved_median_w1": float(np.median(dist)),
        "seed": args.seed,
    }, indent=2) + "\n")
    print(f"mu={mu:.3f} median_w1={np.median(dist):.4g} "
          f"range=({dist.min():.4g}, {dist.max():.4g}) -> {args.output}")


if __name__ == "__main__":
    main()
