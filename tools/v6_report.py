#!/usr/bin/env python3
"""V6 report: real-data negative control through the production matcher.

Validation item V6 (docs/REMAINING_VALIDATION_PROPOSAL.md) is judged on the 300
tests written by slurm/run_v6_negative_control.sbatch under results/v6.

Definitions the proposal leaves open are fixed here:
  rejection     p <= alpha (0.05); a test whose bundle has fewer than 100
                QC-passing sets, or that has no Phi-SFS output, also counts as a
                rejection, as the proposal requires
  upper bound   one-sided 95% Clopper-Pearson bound on the rejection rate
  test blocks   tests 0-49, 50-99, ... in test-ID order
  matcher-order blocks
                the position of each test's B0 in its own matcher order
                (reference_replicate_id 0-109), in blocks of 22
  p distribution
                counts of p in deciles, and of each grid value up to 0.10; with
                R = 99 the null distribution is uniform on {0.01, ..., 1.00}
  focal overlap the number of sites two focal sets share, over all pairs

Outputs (new directory):
  criteria.csv     one row per blocking criterion: value, threshold, pass
  tests.csv        one row per test
  blocks.csv       rejection rates by test block and by matcher-order block
  report.json      everything above plus the p distribution and overlap summary
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

TESTS = 300
NULLS = 99
MIN_QC_SETS = 100
MAX_REJECTIONS = 21
MAX_RATE = 0.075
MAX_UPPER = 0.10
TEST_BLOCK = 50
ORDER_BLOCK = 22


def binomial_cdf(k: int, n: int, p: float) -> float:
    return sum(math.comb(n, i) * p**i * (1 - p) ** (n - i) for i in range(k + 1))


def clopper_pearson_upper(k: int, n: int, level: float = 0.95) -> float:
    """One-sided upper bound: the p at which P(X <= k) = 1 - level."""
    if k >= n:
        return 1.0
    lo, hi = k / n, 1.0
    for _ in range(200):
        mid = (lo + hi) / 2
        if binomial_cdf(k, n, mid) > 1 - level:
            lo = mid
        else:
            hi = mid
    return hi


def read_test(test_dir: Path, alpha: float) -> dict:
    row = {"test": int(test_dir.name.split("_")[1]), "status": "ok"}
    matches = json.loads((test_dir / "matches" / "metadata.json").read_text())
    row["qc_passes"] = int(matches["qc_passes"])
    summary = test_dir / "phi" / "summary.csv"
    if not summary.exists():
        row.update(status="no_phi", null_replicates=0, reject=True)
        return row
    with summary.open() as handle:
        record = next(csv.DictReader(handle))
    row.update(
        null_replicates=int(record["null_replicates_r"]),
        reference_replicate_id=int(record["reference_replicate_id"]),
        observed_phi_sfs=float(record["observed_phi_sfs"]),
        z_score=float(record["z_score"]),
        exceedances=int(record["exceedances"]),
        p_value=float(record["p_value"]),
    )
    # The add-one P-value is recomputed from the exceedance count.
    row["p_recomputed_agrees"] = math.isclose(
        row["p_value"], (1 + row["exceedances"]) / (row["null_replicates"] + 1)
    )
    if row["qc_passes"] < MIN_QC_SETS:
        row["status"] = "too_few_qc_sets"
    row["reject"] = row["status"] != "ok" or row["p_value"] <= alpha
    return row


def focal_overlap(test_dirs: list[Path]) -> dict:
    sets = []
    for test_dir in test_dirs:
        sites = np.loadtxt(test_dir / "focal.pos.txt", dtype=np.int64, ndmin=2)
        keys = sites[:, 0] * 10**10 + sites[:, 1]
        if np.unique(keys).size != keys.size:
            raise ValueError(f"duplicate focal sites in {test_dir}")
        sets.append(keys)
    sizes = sorted({keys.size for keys in sets})
    universe, inverse = np.unique(np.concatenate(sets), return_inverse=True)
    membership = np.zeros((len(sets), universe.size), dtype=np.float32)
    membership[np.repeat(np.arange(len(sets)), [k.size for k in sets]), inverse] = 1
    shared = membership @ membership.T
    pairs = shared[np.triu_indices(len(sets), k=1)].astype(np.int64)
    return {
        "focal_set_sizes": sizes,
        "distinct_sites_across_tests": int(universe.size),
        "pairs": int(pairs.size),
        "shared_sites_mean": float(pairs.mean()),
        "shared_sites_median": float(np.median(pairs)),
        "shared_sites_max": int(pairs.max()),
        "pairs_sharing_any_site": int((pairs > 0).sum()),
        "most_tests_containing_one_site": int(membership.sum(axis=0).max()),
    }


def block_rates(rows: list[dict], key: str, width: int, upper: int) -> list[dict]:
    out = []
    for start in range(0, upper, width):
        members = [r for r in rows if r.get(key) is not None and start <= r[key] < start + width]
        rejected = sum(r["reject"] for r in members)
        out.append({
            "grouping": key, "first": start, "last": min(start + width, upper) - 1,
            "tests": len(members), "rejections": rejected,
            "rate": rejected / len(members) if members else float("nan"),
        })
    return out


def write_csv(path: Path, rows: list[dict]) -> None:
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def run(args: argparse.Namespace) -> dict:
    if args.output.exists():
        raise SystemExit(f"output already exists: {args.output}")
    test_dirs = sorted(args.root.glob("test_[0-9][0-9][0-9][0-9]"))
    rows = [read_test(d, args.alpha) for d in test_dirs]
    n = len(rows)
    rejections = sum(r["reject"] for r in rows)
    upper = clopper_pearson_upper(rejections, n) if n else float("nan")
    p_values = np.array([r["p_value"] for r in rows if "p_value" in r])

    criteria = [
        ("evaluable tests", n, f"== {TESTS}", n == TESTS),
        ("tests with exactly R = 99 nulls",
         sum(r["null_replicates"] == NULLS for r in rows), f"== {TESTS}",
         all(r["null_replicates"] == NULLS for r in rows) and n == TESTS),
        ("recomputed add-one P agrees",
         sum(r.get("p_recomputed_agrees", False) for r in rows), f"== {n}",
         all(r.get("p_recomputed_agrees", False) for r in rows)),
        (f"rejections at p <= {args.alpha}", rejections, f"<= {MAX_REJECTIONS}",
         rejections <= MAX_REJECTIONS),
        ("rejection rate", round(rejections / n, 4), f"<= {MAX_RATE}",
         rejections / n <= MAX_RATE),
        ("one-sided 95% Clopper-Pearson upper bound", round(upper, 4),
         f"< {MAX_UPPER}", upper < MAX_UPPER),
    ]
    blocks = (block_rates(rows, "test", TEST_BLOCK, TESTS)
              + block_rates(rows, "reference_replicate_id", ORDER_BLOCK, 110))
    deciles = np.histogram(p_values, bins=np.linspace(0, 1, 11) + 1e-9)[0]
    grid = {f"{v:.2f}": int(np.isclose(p_values, v).sum())
            for v in np.arange(1, 11) / 100}
    report = {
        "root": str(args.root),
        "alpha": args.alpha,
        "tests": n,
        "status_counts": {s: sum(r["status"] == s for r in rows)
                          for s in sorted({r["status"] for r in rows})},
        "qc_passes_min": min(r["qc_passes"] for r in rows),
        "rejections": rejections,
        "clopper_pearson_upper_95": upper,
        "p_deciles_counts": deciles.tolist(),
        "p_deciles_expected_each": n / 10,
        "p_grid_counts_to_0.10": grid,
        "p_mean": float(p_values.mean()),
        "z_mean": float(np.mean([r["z_score"] for r in rows if "z_score" in r])),
        "focal_overlap": focal_overlap(test_dirs),
        "criteria": [dict(zip(("criterion", "value", "threshold", "pass"), c))
                     for c in criteria],
        "blocks": blocks,
    }
    args.output.mkdir(parents=True)
    write_csv(args.output / "criteria.csv", report["criteria"])
    write_csv(args.output / "tests.csv", rows)
    write_csv(args.output / "blocks.csv", blocks)
    (args.output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    for c in criteria:
        print(f"{str(c[3]):8} {c[0]}: {c[1]} ({c[2]})")
    return report


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--root", type=Path, default=Path("results/v6"))
    p.add_argument("--output", type=Path, required=True, help="new output directory")
    p.add_argument("--alpha", type=float, default=0.05)
    return p.parse_args(argv)


def main(argv=None) -> int:
    report = run(parse_args(argv))
    print(json.dumps({k: report[k] for k in (
        "status_counts", "qc_passes_min", "p_deciles_counts", "p_grid_counts_to_0.10",
        "p_mean", "z_mean", "focal_overlap")}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
