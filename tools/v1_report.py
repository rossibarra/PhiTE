#!/usr/bin/env python3
"""V1 report: all-mixture arms on the unpolarised and refhap replicates.

Validation item V1 (docs/REMAINING_VALIDATION_PROPOSAL.md) is judged on the
tests.csv written by tools/sim_polarity_arms.py for each condition.

Definitions the proposal leaves open are fixed here:
  cell          condition x mixture arm x M
  evaluable     status "ok": at least 19 nulls, so p <= 0.05 is attainable;
                the reuse arms always qualify, the depleted arms may not
  rejection     p <= alpha (0.05)
  rate          the mean over replicates of each replicate's rejection rate
                among its evaluable tests (equal replicate weights); a
                replicate with no evaluable test in a cell is left out of that
                cell and counted in replicates_without_tests
  upper bound   the 95th percentile of 20,000 replicate-bootstrap means
                (one-sided 95%), seeded from the cell name
  pass          rate <= 0.075 and upper bound <= 0.10; "n/a" for a cell with no
                evaluable test, which the criterion does not cover (V2 shows the
                depleted arms cannot reach 19 nulls at M = 4000)

Outputs (new directory):
  criteria.csv     one row per cell: denominators, rate, upper bound, pass
  report.json      the rows plus mean Phi_obs - null mean by M (reported only)
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

ARMS = ("mixture_true_reuse", "mixture_true_depleted",
        "mixture_inferred_reuse", "mixture_inferred_depleted")
MAX_RATE = 0.075
MAX_UPPER = 0.10
DRAWS = 20_000


def one_sided_upper(rates: np.ndarray, seed: int, level: float = 0.95) -> float:
    rng = np.random.default_rng(seed)
    means = rng.choice(rates, size=(DRAWS, rates.size), replace=True).mean(axis=1)
    return float(np.quantile(means, level))


def cells(condition: str, tests_csv: Path, alpha: float) -> list[dict]:
    with tests_csv.open() as handle:
        rows = [r for r in csv.DictReader(handle) if r["arm"] in ARMS]
    replicates = sorted({r["replicate"] for r in rows})
    out = []
    for arm in ARMS:
        for m in sorted({int(r["M"]) for r in rows}):
            group = [r for r in rows if r["arm"] == arm and int(r["M"]) == m]
            ok = [r for r in group if r["status"] == "ok"]
            by_rep: dict[str, list[bool]] = {}
            for r in ok:
                by_rep.setdefault(r["replicate"], []).append(float(r["p"]) <= alpha)
            rates = np.array([np.mean(v) for v in by_rep.values()])
            name = f"{condition}\0{arm}\0{m}\0v1-upper"
            seed = int.from_bytes(hashlib.sha256(name.encode()).digest()[:8], "little")
            upper = one_sided_upper(rates, seed) if rates.size else float("nan")
            rate = float(rates.mean()) if rates.size else float("nan")
            excess = [float(r["phi_obs"]) - float(r["null_mean"]) for r in ok]
            out.append({
                "condition": condition, "arm": arm, "M": m,
                "tests": len(group), "evaluable_tests": len(ok),
                "replicates": len(replicates),
                "replicates_with_tests": int(rates.size),
                "replicates_without_tests": len(replicates) - int(rates.size),
                "rejections": int(sum(float(r["p"]) <= alpha for r in ok)),
                "rate": rate, "upper_95": upper,
                "replicate_rate_min": float(rates.min()) if rates.size else float("nan"),
                "replicate_rate_max": float(rates.max()) if rates.size else float("nan"),
                "mean_phi_excess": float(np.mean(excess)) if excess else float("nan"),
                "pass": (bool(rate <= MAX_RATE and upper <= MAX_UPPER)
                         if rates.size else "n/a"),
            })
    return out


def run(args: argparse.Namespace) -> list[dict]:
    if args.output.exists():
        raise SystemExit(f"output already exists: {args.output}")
    rows = []
    for spec in args.condition:
        condition, path = spec.split("=", 1)
        rows += cells(condition, Path(path) / "tests.csv", args.alpha)
    args.output.mkdir(parents=True)
    with (args.output / "criteria.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    report = {
        "conditions": args.condition, "alpha": args.alpha,
        "max_rate": MAX_RATE, "max_upper": MAX_UPPER, "bootstrap_draws": DRAWS,
        "applicable_cells": sum(r["pass"] != "n/a" for r in rows),
        "all_applicable_cells_pass": all(r["pass"] for r in rows if r["pass"] != "n/a"),
        "cells": rows,
    }
    (args.output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    return rows


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--condition", action="append", required=True,
                   help="NAME=ARMS_DIR, repeatable")
    p.add_argument("--output", type=Path, required=True, help="new output directory")
    p.add_argument("--alpha", type=float, default=0.05)
    return p.parse_args(argv)


def main(argv=None) -> int:
    for r in run(parse_args(argv)):
        print(f"{str(r['pass']):6} {r['condition']:12} {r['arm']:26} M={r['M']:<5}"
              f" eval {r['evaluable_tests']:>3}/{r['tests']:<3}"
              f" reps {r['replicates_with_tests']:>2}/{r['replicates']}"
              f" rate {r['rate']:.3f} up95 {r['upper_95']:.3f}"
              f" excess {r['mean_phi_excess']:+.5f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
