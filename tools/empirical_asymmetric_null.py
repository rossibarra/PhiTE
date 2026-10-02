#!/usr/bin/env python3
"""Estimate SNP type-I error for the Bernoulli-q asymmetric Phi-SFS null.

This is a held-out empirical check, not a demographic simulation. It reads
completed, globally disjoint matcher work files, randomly reserves real matched
SNP sets as pseudo-focal A sets, and uses the remaining sets as one q-mixture
reference plus Bernoulli(q)-hard null-left sets.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

from normalize_tes.phi_sfs import PROJECTION_SIZE, hypergeometric_projection


def load_best_rows(work_dir: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    paths = sorted((work_dir / "replicates").glob("replicate-*.npz"))
    if not paths:
        raise ValueError(f"no completed replicate files under {work_dir / 'replicates'}")
    rows, best_distance, bootstrap_distance = [], [], []
    expected = list(range(len(paths)))
    observed = [int(path.stem.split("-")[1]) for path in paths]
    if observed != expected:
        raise ValueError("replicate files are not a complete zero-based prefix")
    for path in paths:
        with np.load(path, allow_pickle=False) as bundle:
            metadata = json.loads(str(bundle["metadata"].item()))
            distances = np.asarray([entry["best_distance"] for entry in metadata])
            best = int(np.argmin(distances))
            rows.append(np.asarray(bundle["rows"][best], dtype=np.int64))
            best_distance.append(float(distances[best]))
            bootstrap_distance.append(float(bundle["bootstrap_distance"]))
    return (
        np.stack(rows),
        np.asarray(best_distance),
        np.asarray(bootstrap_distance),
    )


def site_arrays(rows: np.ndarray, eligibility: Path) -> tuple[np.ndarray, ...]:
    eligible_rows = np.load(eligibility / "row_indices.npy", mmap_mode="r")
    flat = rows.reshape(-1)
    index = np.searchsorted(eligible_rows, flat)
    if np.any(index >= eligible_rows.size) or np.any(eligible_rows[index] != flat):
        raise ValueError("a matched row is absent from the VCF eligibility artifact")
    shape = rows.shape
    return tuple(
        np.asarray(np.load(eligibility / name, mmap_mode="r")[index]).reshape(shape)
        for name in ("alt_counts.npy", "callable_counts.npy", "p_alt_derived.npy")
    )


def spectra(
    alt: np.ndarray,
    callable_count: np.ndarray,
    q: np.ndarray,
    *,
    seed: int,
) -> tuple[np.ndarray, np.ndarray]:
    pairs = np.stack([alt.reshape(-1), callable_count.reshape(-1)], axis=1)
    unique, inverse = np.unique(pairs, axis=0, return_inverse=True)
    projected = np.stack([
        hypergeometric_projection(int(k), int(n))[1:PROJECTION_SIZE]
        for k, n in unique
    ])
    base = projected[inverse].reshape(*alt.shape, PROJECTION_SIZE - 1)
    reverse = base[..., ::-1]
    q3 = q[..., None]
    mixture = (q3 * base + (1.0 - q3) * reverse).sum(axis=1)
    rng = np.random.default_rng(seed)
    alt_derived = rng.random(q.shape) < q
    hard = np.where(alt_derived[..., None], base, reverse).sum(axis=1)
    mixture /= mixture.sum(axis=1, keepdims=True)
    hard /= hard.sum(axis=1, keepdims=True)
    return mixture, hard


def phi_to_reference(sfs: np.ndarray, reference: np.ndarray) -> np.ndarray:
    cdf = np.cumsum(sfs, axis=1)
    reference_cdf = np.cumsum(reference)
    return np.abs(cdf[:, :-1] - reference_cdf[:-1]).sum(axis=1) / PROJECTION_SIZE


def wilson(x: int, n: int, z: float = 1.959963984540054) -> list[float]:
    p = x / n
    denominator = 1.0 + z * z / n
    center = (p + z * z / (2 * n)) / denominator
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denominator
    return [center - half, center + half]


def run(args: argparse.Namespace) -> dict:
    rows, match_error, bootstrap_error = load_best_rows(args.work_dir)
    if rows.shape[0] < args.a_sets + 3:
        raise ValueError("not enough completed sets for held-out A and calibration sets")
    alt, callable_count, q = site_arrays(rows, args.eligibility)
    mixture, hard = spectra(alt, callable_count, q, seed=args.orientation_seed)

    rng = np.random.default_rng(args.split_seed)
    order = rng.permutation(rows.shape[0])
    a_index = order[:args.a_sets]
    calibration = order[args.a_sets:]
    reference_index = int(calibration[0])
    null_index = calibration[1:]
    null_distance = phi_to_reference(hard[null_index], mixture[reference_index])
    observed_distance = phi_to_reference(hard[a_index], mixture[reference_index])
    p_values = np.asarray([
        (1 + np.count_nonzero(null_distance >= value)) / (null_distance.size + 1)
        for value in observed_distance
    ])
    rejected = p_values <= args.alpha
    x = int(rejected.sum())
    ratio = match_error / np.maximum(bootstrap_error, np.finfo(float).tiny)
    result = {
        "design": "held-out-real-SNP Bernoulli(q)-hard vs q-mixture asymmetric null",
        "work_dir": str(args.work_dir.resolve()),
        "eligibility": str(args.eligibility.resolve()),
        "completed_sets": int(rows.shape[0]),
        "sites_per_set": int(rows.shape[1]),
        "held_out_a_sets": int(a_index.size),
        "null_sets": int(null_index.size),
        "reference_replicate": reference_index,
        "split_seed": args.split_seed,
        "orientation_seed": args.orientation_seed,
        "alpha": args.alpha,
        "rejections": x,
        "rejection_rate": x / a_index.size,
        "wilson_95_descriptive": wilson(x, a_index.size),
        "p_value_min": float(p_values.min()),
        "p_value_median": float(np.median(p_values)),
        "p_value_max": float(p_values.max()),
        "matching_error_ratio_min": float(ratio.min()),
        "matching_error_ratio_median": float(np.median(ratio)),
        "matching_error_ratio_max": float(ratio.max()),
        "matching_qc_ratio_below_0.5": int(np.count_nonzero(ratio < 0.5)),
        "limitations": [
            "all 100 tests share one mixture reference and one empirical null vector, so their p-values are dependent",
            "the interrupted work prefix was not published as a completed match bundle",
            "matching QC is reported and must be considered before interpreting the rejection rate",
        ],
    }
    args.output.mkdir(parents=True, exist_ok=False)
    with (args.output / "tests.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=[
            "a_replicate", "observed_phi", "p_value", "reject",
            "matching_error", "bootstrap_error", "matching_error_ratio",
        ])
        writer.writeheader()
        for index, distance, p_value, reject in zip(
            a_index, observed_distance, p_values, rejected
        ):
            writer.writerow({
                "a_replicate": int(index),
                "observed_phi": float(distance),
                "p_value": float(p_value),
                "reject": bool(reject),
                "matching_error": float(match_error[index]),
                "bootstrap_error": float(bootstrap_error[index]),
                "matching_error_ratio": float(ratio[index]),
            })
    np.save(args.output / "null_phi.npy", null_distance, allow_pickle=False)
    np.save(args.output / "a_replicate.npy", a_index, allow_pickle=False)
    np.save(args.output / "null_replicate.npy", null_index, allow_pickle=False)
    with (args.output / "summary.json").open("w", encoding="utf-8") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--eligibility", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--a-sets", type=int, default=100)
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--split-seed", type=int, default=20260926)
    parser.add_argument("--orientation-seed", type=int, default=2001)
    return parser.parse_args()


def main() -> int:
    result = run(parse_args())
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
