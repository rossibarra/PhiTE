#!/usr/bin/env python3
"""V2 capacity report for the depleted all-mixture simulation arms.

Validation item V2 (docs/REMAINING_VALIDATION_PROPOSAL.md) asks, for every
replicate and M, how many control sets sequential depletion completed, how often
at least 20 were completed, and which age bin blocked the next set and how much
capacity it had left.

`tools.sim_polarity_arms` stops at the first infeasible set but does not record
why. Rather than change that validated code, this reruns it unmodified, with
`bootstrap_matched_sets` wrapped. After each call, the wrapper recomputes the
remaining pool from the sets that call returned and identifies the first bin
whose remaining members fall short of the next target. The wrapper makes no
random draws, so the run is the original run. That is checked: every depleted
row's completed-set count must equal the reference `tests.csv`.

Outputs (new directory): the simulator's usual files, plus
  capacity_by_test.csv   one row per depleted test
  capacity_summary.csv   per replicate x arm x M
  capacity_overall.csv   per arm x M, pooled over replicates
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from tools import sim_polarity_arms as sim

MIN_NULLS = 19  # nulls needed for P <= 0.05 to be attainable


def _limiting_bin(pool, bins, completed, targets, disjoint):
    """First bin that blocks set len(completed), with its remaining capacity."""
    k = completed.shape[0]
    if k >= targets.shape[0]:
        return None
    available = np.asarray(pool, dtype=np.int64)
    if disjoint and k:
        available = available[~np.isin(available, completed.ravel())]
    want = targets[k]
    for bin_index in np.flatnonzero(want):
        have = int(np.count_nonzero(bins[available] == bin_index))
        if have < int(want[bin_index]):
            return {"limiting_bin": int(bin_index), "remaining": have,
                    "needed": int(want[bin_index])}
    raise RuntimeError("stopped early without an infeasible bin; wrapper out of sync")


def run(args: argparse.Namespace) -> None:
    records: list[dict] = []
    original_bootstrap = sim.bootstrap_matched_sets
    original_design = sim.run_mixture_design
    last: dict = {}

    def recording_bootstrap(focal, pool, site_bins, n_sets, rng, n_bins, *,
                            disjoint, target_counts=None):
        completed = original_bootstrap(focal, pool, site_bins, n_sets, rng, n_bins,
                                       disjoint=disjoint, target_counts=target_counts)
        if target_counts is None:
            raise RuntimeError("mixture arms always pass target_counts")
        last.clear()
        last.update(completed=int(completed.shape[0]), limit=_limiting_bin(
            pool, np.asarray(site_bins), completed,
            np.asarray(target_counts), disjoint))
        return completed

    def recording_design(*, arm, **kwargs):
        result = original_design(arm=arm, **kwargs)
        records.append({"arm": arm, "M": int(kwargs["a"].size),
                        "completed": last["completed"], **(last["limit"] or {})})
        return result

    sim.bootstrap_matched_sets = recording_bootstrap
    sim.run_mixture_design = recording_design
    try:
        sim.main([
            "--sim-root", str(args.sim_root), "--prep-root", str(args.prep_root),
            "--output", str(args.output),
            "--replicates", *map(str, args.replicates),
            "--sizes", *map(str, args.sizes),
            "--tests", str(args.tests), "--nulls", str(args.nulls),
            "--mixture-sets", str(args.mixture_sets),
            "--min-mixture-sets", str(args.min_mixture_sets),
            "--n-bins", str(args.n_bins), "--polarity-seed", str(args.polarity_seed),
            "--seed", str(args.seed), "--alpha", str(args.alpha),
        ])
    finally:
        sim.bootstrap_matched_sets = original_bootstrap
        sim.run_mixture_design = original_design

    rows = [r for r in csv.DictReader((args.output / "tests.csv").open())
            if r["arm"].startswith("mixture_")]
    if len(rows) != len(records):
        raise RuntimeError(f"{len(rows)} mixture rows but {len(records)} records")
    by_test = []
    for row, rec in zip(rows, records):
        if row["arm"] != rec["arm"] or int(row["M"]) != rec["M"]:
            raise RuntimeError("record order does not match tests.csv")
        if int(row["control_sets"]) != rec["completed"]:
            raise RuntimeError("completed-set count disagrees with tests.csv")
        if row["arm"].endswith("_depleted"):
            by_test.append({"replicate": row["replicate"], "arm": row["arm"],
                            "M": rec["M"], "test": row["test"], "status": row["status"],
                            "control_sets": rec["completed"],
                            "null_replicates": int(row["null_replicates"]),
                            "limiting_bin": rec.get("limiting_bin", ""),
                            "limiting_bin_remaining": rec.get("remaining", ""),
                            "limiting_bin_needed": rec.get("needed", "")})

    if args.reference_tests:
        reference = {(r["replicate"], r["M"], r["test"], r["arm"]): r["control_sets"]
                     for r in csv.DictReader(args.reference_tests.open())
                     if r["arm"].endswith("_depleted")}
        mismatched = [t for t in by_test
                      if reference.get((t["replicate"], str(t["M"]), t["test"], t["arm"]))
                      != str(t["control_sets"])]
        if len(reference) != len(by_test) or mismatched:
            raise RuntimeError(
                f"rerun does not reproduce {args.reference_tests}: "
                f"{len(mismatched)} of {len(by_test)} depleted rows differ")

    def summarize(group):
        sets = np.array([t["control_sets"] for t in group])
        nulls = np.array([t["null_replicates"] for t in group])
        limits = [t["limiting_bin"] for t in group if t["limiting_bin"] != ""]
        values, counts = np.unique(limits, return_counts=True) if limits else ([], [])
        order = np.argsort(counts)[::-1]
        return {
            "tests": len(group),
            "sets_min": int(sets.min()), "sets_q25": float(np.quantile(sets, 0.25)),
            "sets_median": float(np.median(sets)), "sets_q75": float(np.quantile(sets, 0.75)),
            "sets_max": int(sets.max()),
            "fraction_ge_20_sets": float(np.mean(sets >= 20)),
            "fraction_ge_19_nulls": float(np.mean(nulls >= MIN_NULLS)),
            "tests_evaluable": int(np.sum(nulls >= MIN_NULLS)),
            "tests_stopped_early": len(limits),
            "limiting_bins": ";".join(
                f"{int(values[i])}:{int(counts[i])}" for i in order),
            "limiting_remaining_median": (
                float(np.median([t["limiting_bin_remaining"] for t in group
                                 if t["limiting_bin"] != ""])) if limits else ""),
            "limiting_needed_median": (
                float(np.median([t["limiting_bin_needed"] for t in group
                                 if t["limiting_bin"] != ""])) if limits else ""),
        }

    summary, overall = [], []
    keys = sorted({(t["replicate"], t["arm"], t["M"]) for t in by_test})
    for rep, arm, m in keys:
        group = [t for t in by_test if (t["replicate"], t["arm"], t["M"]) == (rep, arm, m)]
        summary.append({"replicate": rep, "arm": arm, "M": m, **summarize(group)})
    for arm, m in sorted({(t["arm"], t["M"]) for t in by_test}):
        group = [t for t in by_test if (t["arm"], t["M"]) == (arm, m)]
        overall.append({"arm": arm, "M": m, **summarize(group)})
    for name, table in (("capacity_by_test.csv", by_test),
                        ("capacity_summary.csv", summary),
                        ("capacity_overall.csv", overall)):
        with (args.output / name).open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(table[0]))
            writer.writeheader()
            writer.writerows(table)
    supported = {}
    for arm in sorted({r["arm"] for r in overall}):
        ok = [r["M"] for r in overall if r["arm"] == arm and r["fraction_ge_19_nulls"] >= 0.95]
        supported[arm] = max(ok) if ok else None
    (args.output / "capacity_run.json").write_text(json.dumps({
        "reference_tests": str(args.reference_tests) if args.reference_tests else None,
        "reproduced_reference": bool(args.reference_tests),
        "largest_tested_M_with_95pct_tests_ge_19_nulls": supported,
        "tested_sizes": args.sizes,
    }, indent=2) + "\n")


def parse_args(argv=None) -> argparse.Namespace:
    base = sim.parse_args(["--output", "unused"])
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--output", type=Path, required=True, help="new output directory")
    p.add_argument("--reference-tests", type=Path,
                   help="tests.csv of the run being explained; the rerun must match it")
    for name in ("sim_root", "prep_root"):
        p.add_argument(f"--{name.replace('_', '-')}", type=Path, default=getattr(base, name))
    p.add_argument("--replicates", type=int, nargs="+", default=base.replicates)
    p.add_argument("--sizes", type=int, nargs="+", default=base.sizes)
    for name in ("tests", "nulls", "mixture_sets", "min_mixture_sets", "n_bins",
                 "polarity_seed", "seed"):
        p.add_argument(f"--{name.replace('_', '-')}", type=int, default=getattr(base, name))
    p.add_argument("--alpha", type=float, default=base.alpha)
    return p.parse_args(argv)


def main(argv=None) -> int:
    run(parse_args(argv))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
