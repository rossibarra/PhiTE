#!/usr/bin/env python3
"""Simulation test of Phi-SFS polarity and filtering choices (review 12).

Uses the dnAging msprime/SINGER replicates, where each site's true mutation age
and true polarity are known (the simulated VCF writes the ancestral allele as
REF), and PhiTE's own interval store, ancestral-state table and per-draw
ancestral table built from the posterior ARG draws by
`slurm/run_dnaging_sim_prep.sbatch`.

Every simulated site is neutral, so every test here is a null test: a
well-calibrated design rejects at the nominal rate. A pseudo-TE focal set A is a
random set of sites given its TRUE hard polarity, which is what PhiTE assumes of
a TE. Controls are matched to A on age by stratified sampling in log-age
quantile bins. That is a stand-in for the production CDF matcher, not the
matcher itself.

Arms (review 12, "Proposed validation study"):

  oracle              true-age matching; A, B0 and every B_i true-hard
  production          true-age matching; A true-hard, B0 q-mixture, B_i Bernoulli(q)-hard
  snp_pilot           as production but A is Bernoulli(q)-hard (the -A SNP design)
  production_inferred inferred all-draw ages on both sides; polarity as production
  te_filter_only      as production_inferred, A kept only where support for its true
                      orientation is >= 0.5 (the TE rule); A ages still all-draw
  te_full             as te_filter_only, and A ages from agreeing draws only, with
                      PhiTE's fallback to all draws when none agrees

oracle/production/snp_pilot share each test's A, B0 and nulls, so they differ
only in polarity. production_inferred/te_filter_only/te_full share A before
filtering. Bernoulli orientations use PhiTE's coordinate-keyed draw, so a site
has one orientation wherever it is used, as in production.

Outputs (new directory): tests.csv, summary.csv, q_calibration.csv,
age_error.csv, filter_by_count.csv, run.json.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import time
from pathlib import Path

import numpy as np

from normalize_tes.build_draw_polarity import open_draw_polarity
from normalize_tes.phi_sfs import (
    PROJECTION_SIZE,
    RETAINED_BINS,
    _site_uniform,
    calibrate_phi,
    hypergeometric_projection,
    phi_sfs,
)
from normalize_tes.snp_age_store import open_snp_age_store

BASES = "ACGT"
CHROM = "chr1"
TRUE_AGE_ARMS = ("oracle", "production", "snp_pilot")
INFERRED_ARMS = ("production_inferred", "te_filter_only", "te_full")
ARMS = TRUE_AGE_ARMS + INFERRED_ARMS
COUNT_BINS = ((1, 1), (2, 2), (3, 4), (5, 8), (9, 13), (14, 18), (19, 22), (23, 25))


def load_replicate(sim_root: Path, prep_root: Path, replicate: int) -> dict:
    """Join truth, genotypes, store rows, polarity support and inferred ages."""
    label = f"replicate_{replicate:03d}"
    sim = sim_root / "simulations" / label
    truth = {}
    with gzip.open(sim / "mutation_truth.tsv.gz", "rt") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            truth[int(row["site_id"])] = (float(row["mutation_time"]),
                                          int(row["modern_derived_count"]))
    singer_pos = {}
    with (sim_root / "singer" / label / "position_map.tsv").open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            singer_pos[int(row["site_id"])] = int(row["singer_pos"])
    genotypes = {}
    n_modern = None
    with gzip.open(sim / "all_samples.vcf.gz", "rt") as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            fields = line.rstrip("\n").split("\t")
            if line.startswith("#"):
                modern = [i for i, name in enumerate(fields) if name.startswith("modern_")]
                n_modern = len(modern)
                continue
            alt_count = sum(fields[i] == "1" for i in modern)
            genotypes[int(fields[2])] = (fields[3], fields[4], alt_count)

    prep = prep_root / label
    store = open_snp_age_store(prep / "interval_store")
    ancestral = np.load(prep / "ancestral_states" / "ancestral_counts.npy")
    draw_base, _ = open_draw_polarity(prep / "draw_polarity", store)
    positions = np.asarray(store.positions, dtype=np.float64)
    eligible = np.asarray(store.eligible)

    site_ids, rows, ref_index, alt_index = [], [], [], []
    ages, counts = [], []
    dropped = {"no_store_row": 0, "ineligible": 0, "genotype_mismatch": 0,
               "non_acgt": 0, "unoriented": 0}
    for site_id, (age, derived) in truth.items():
        if site_id not in genotypes or site_id not in singer_pos:
            dropped["no_store_row"] += 1
            continue
        ref, alt, alt_count = genotypes[site_id]
        if alt_count != derived:
            dropped["genotype_mismatch"] += 1
            continue
        if ref not in BASES or alt not in BASES:
            dropped["non_acgt"] += 1
            continue
        target = float(singer_pos[site_id])
        index = int(np.searchsorted(positions, target))
        if index >= positions.size or positions[index] != target:
            dropped["no_store_row"] += 1
            continue
        if not eligible[index]:
            dropped["ineligible"] += 1
            continue
        oriented = ancestral[index, BASES.index(ref)] + ancestral[index, BASES.index(alt)]
        if oriented <= 0:
            dropped["unoriented"] += 1
            continue
        site_ids.append(site_id)
        rows.append(index)
        ref_index.append(BASES.index(ref))
        alt_index.append(BASES.index(alt))
        ages.append(age)
        counts.append(derived)
    if not rows:
        raise ValueError(f"{label}: no usable sites")
    rows = np.asarray(rows, dtype=np.int64)
    ref_index = np.asarray(ref_index)
    alt_index = np.asarray(alt_index)
    ref_calls = ancestral[rows, ref_index].astype(np.float64)
    alt_calls = ancestral[rows, alt_index].astype(np.float64)
    # REF is the true ancestral allele, so q = P(ALT derived) is also the
    # posterior support for the true orientation.
    q = ref_calls / (ref_calls + alt_calls)

    # Inferred ages: interval midpoints averaged over all draws, and over the
    # draws whose ancestral call is the true ancestral allele. A site with no
    # agreeing interval keeps all of its intervals, as te_age_target does.
    batch = store.intervals(rows)
    all_age = np.empty(rows.size)
    agree_age = np.empty(rows.size)
    fell_back = 0
    for i in range(rows.size):
        start, stop = int(batch.offsets[i]), int(batch.offsets[i + 1])
        mid = (np.asarray(batch.below[start:stop], dtype=np.float64)
               + np.asarray(batch.above[start:stop], dtype=np.float64)) / 2.0
        draws = np.asarray(batch.draw_id[start:stop], dtype=np.int64)
        all_age[i] = mid.mean()
        agree = np.asarray(draw_base[rows[i], draws]) == ref_index[i]
        if not agree.any():
            agree = np.ones(draws.size, dtype=bool)
            fell_back += 1
        agree_age[i] = mid[agree].mean()

    return {
        "label": label,
        "n_modern": int(n_modern),
        "site_id": np.asarray(site_ids, dtype=np.int64),
        "position": np.asarray([singer_pos[s] for s in site_ids], dtype=np.int64),
        "true_age": np.asarray(ages, dtype=np.float64),
        "derived": np.asarray(counts, dtype=np.int64),
        "q": q,
        "all_age": all_age,
        "agree_age": agree_age,
        "fell_back": fell_back,
        "dropped": dropped,
        "n_draws": int(draw_base.shape[1]),
    }


def site_spectra(data: dict, seed: int) -> dict[str, np.ndarray]:
    """Per-site retained projections (bins 1..19) under each polarity rule."""
    n = data["n_modern"]
    cache = {k: hypergeometric_projection(int(k), n) for k in np.unique(data["derived"])}
    true = np.stack([cache[int(k)] for k in data["derived"]])
    q = data["q"][:, None]
    mixture = q * true + (1.0 - q) * true[:, ::-1]
    uniform = np.asarray([_site_uniform(seed, (CHROM, int(p))) for p in data["position"]])
    bern = np.where((uniform < data["q"])[:, None], true, true[:, ::-1])
    keep = slice(1, PROJECTION_SIZE)
    return {"true": true[:, keep], "mixture": mixture[:, keep], "bernoulli": bern[:, keep]}


def phi_many(a: np.ndarray, spectra: np.ndarray) -> np.ndarray:
    """Vectorized phi_sfs(a, b).value for normalized a and each row of spectra."""
    grid = RETAINED_BINS.astype(np.float64) / PROJECTION_SIZE
    residual = np.cumsum(spectra, axis=1) - np.cumsum(a)[None, :]
    return np.sum(np.abs(residual[:, :-1]) * np.diff(grid)[None, :], axis=1)


def normalize(raw: np.ndarray) -> np.ndarray:
    return raw / raw.sum(axis=-1, keepdims=True)


def age_bins(key: np.ndarray, n_bins: int) -> np.ndarray:
    """Assign each site to a log-age quantile bin computed over all sites."""
    edges = np.quantile(np.log(key), np.linspace(0.0, 1.0, n_bins + 1)[1:-1])
    return np.searchsorted(edges, np.log(key), side="right")


def matched_sets(focal_bins: np.ndarray, pool: np.ndarray, pool_bins: np.ndarray,
                 n_sets: int, rng: np.random.Generator, n_bins: int) -> np.ndarray | None:
    """Draw n_sets sets matching focal_bins' per-bin counts, without replacement
    within a set. Returns (n_sets, M) site indices, or None if a bin is short."""
    want = np.bincount(focal_bins, minlength=n_bins)
    columns = []
    for b in np.flatnonzero(want):
        members = pool[pool_bins == b]
        if members.size < want[b]:
            return None
        pick = np.argpartition(rng.random((n_sets, members.size)), want[b] - 1, axis=1)
        columns.append(members[pick[:, :want[b]]])
    return np.concatenate(columns, axis=1)


def true_age_w1(a: np.ndarray, b: np.ndarray, age: np.ndarray) -> float:
    """W1 (generations) between the true ages of two equal-size sets."""
    return float(np.mean(np.abs(np.sort(age[a]) - np.sort(age[b]))))


def run_test(arm_polarity: dict, spectra: dict, a: np.ndarray, b0: np.ndarray,
             nulls: np.ndarray) -> dict:
    """Phi_obs, null distances and calibration for one test."""
    a_spec = normalize(spectra[arm_polarity["a"]][a].sum(axis=0))
    b0_spec = normalize(spectra[arm_polarity["b0"]][b0].sum(axis=0))
    null_spec = normalize(spectra[arm_polarity["null"]][nulls].sum(axis=1))
    observed = phi_sfs(a_spec, b0_spec)
    null = phi_many(b0_spec, null_spec)
    check = phi_sfs(null_spec[0], b0_spec).value
    if not math.isclose(null[0], check, rel_tol=1e-9, abs_tol=1e-12):
        raise AssertionError("vectorized Phi disagrees with phi_sfs")
    cal = calibrate_phi(observed.value, null)
    return {"phi_obs": cal.observed, "null_mean": cal.null_mean, "null_sd": cal.null_sd,
            "z": cal.z_score, "p": cal.p_value,
            "mean_daf_difference": observed.mean_daf_difference}


POLARITY = {
    "oracle": {"a": "true", "b0": "true", "null": "true"},
    "production": {"a": "true", "b0": "mixture", "null": "bernoulli"},
    "snp_pilot": {"a": "bernoulli", "b0": "mixture", "null": "bernoulli"},
    "production_inferred": {"a": "true", "b0": "mixture", "null": "bernoulli"},
    "te_filter_only": {"a": "true", "b0": "mixture", "null": "bernoulli"},
    "te_full": {"a": "true", "b0": "mixture", "null": "bernoulli"},
}


def replicate_tests(data: dict, args: argparse.Namespace, rng: np.random.Generator,
                    writer: csv.DictWriter) -> None:
    spectra = site_spectra(data, args.polarity_seed)
    n_sites = data["q"].size
    everyone = np.arange(n_sites)
    true_bins = age_bins(data["true_age"], args.n_bins)
    all_bins = age_bins(data["all_age"], args.n_bins)
    # Agreeing-draw ages are binned on the all-draw edges, so A and B share one
    # age scale; te_full differs from te_filter_only only in A's age values.
    edges = np.quantile(np.log(data["all_age"]),
                        np.linspace(0.0, 1.0, args.n_bins + 1)[1:-1])
    agree_bins = np.searchsorted(edges, np.log(data["agree_age"]), side="right")
    kept = data["q"] >= 0.5

    for m in args.sizes:
        if 2 * m + m > n_sites:
            print(f"{data['label']}: skip M={m}, only {n_sites} sites", flush=True)
            continue
        for t in range(args.tests):
            base = {"replicate": data["label"], "M": m, "test": t}
            a = rng.choice(everyone, size=m, replace=False)

            # True-age arms: one A, B0 and null collection shared by all three.
            rest = np.setdiff1d(everyone, a, assume_unique=True)
            b0 = matched_sets(true_bins[a], rest, true_bins[rest], 1, rng, args.n_bins)
            rows = []
            if b0 is not None:
                b0 = b0[0]
                pool = np.setdiff1d(rest, b0, assume_unique=True)
                nulls = matched_sets(true_bins[a], pool, true_bins[pool], args.nulls,
                                     rng, args.n_bins)
                if nulls is not None:
                    w1 = true_age_w1(a, b0, data["true_age"])
                    for arm in TRUE_AGE_ARMS:
                        rows.append({**base, "arm": arm, "M_eff": m, "status": "ok",
                                     "true_age_w1_a_b0": w1,
                                     **run_test(POLARITY[arm], spectra, a, b0, nulls)})
            if not rows:
                rows = [{**base, "arm": arm, "M_eff": m, "status": "bin_short"}
                        for arm in TRUE_AGE_ARMS]

            # Inferred-age arms share the same A before filtering.
            variants = {
                "production_inferred": (a, all_bins),
                "te_filter_only": (a[kept[a]], all_bins),
                "te_full": (a[kept[a]], agree_bins),
            }
            for arm, (focal, focal_bin_source) in variants.items():
                row = {**base, "arm": arm, "M_eff": int(focal.size)}
                if focal.size < 20:
                    rows.append({**row, "status": "too_few"})
                    continue
                rest = np.setdiff1d(everyone, focal, assume_unique=True)
                fb = focal_bin_source[focal]
                b0 = matched_sets(fb, rest, all_bins[rest], 1, rng, args.n_bins)
                if b0 is None:
                    rows.append({**row, "status": "bin_short"})
                    continue
                b0 = b0[0]
                pool = np.setdiff1d(rest, b0, assume_unique=True)
                nulls = matched_sets(fb, pool, all_bins[pool], args.nulls, rng, args.n_bins)
                if nulls is None:
                    rows.append({**row, "status": "bin_short"})
                    continue
                rows.append({**row, "status": "ok",
                             "true_age_w1_a_b0": true_age_w1(focal, b0, data["true_age"]),
                             **run_test(POLARITY[arm], spectra, focal, b0, nulls)})
            for row in rows:
                writer.writerow(row)


def wilson(k: int, n: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if n == 0:
        return (math.nan, math.nan)
    p = k / n
    centre = (p + z * z / (2 * n)) / (1 + z * z / n)
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return (centre - half, centre + half)


def summarize(tests_path: Path, alpha: float) -> list[dict]:
    with tests_path.open() as handle:
        rows = [r for r in csv.DictReader(handle) if r["status"] == "ok"]
    out = []
    for arm in ARMS:
        for m in sorted({int(r["M"]) for r in rows}):
            group = [r for r in rows if r["arm"] == arm and int(r["M"]) == m]
            if not group:
                continue
            p = np.asarray([float(r["p"]) for r in group])
            excess = np.asarray([float(r["phi_obs"]) - float(r["null_mean"]) for r in group])
            z = np.asarray([float(r["z"]) for r in group])
            reject = int(np.count_nonzero(p <= alpha))
            by_rep = {}
            for r in group:
                by_rep.setdefault(r["replicate"], []).append(float(r["p"]) <= alpha)
            rep_rates = [float(np.mean(v)) for v in by_rep.values()]
            lo, hi = wilson(reject, p.size)
            out.append({
                "arm": arm, "M": m, "tests": int(p.size),
                "median_M_eff": float(np.median([int(r["M_eff"]) for r in group])),
                "rejections": reject, "rejection_rate": reject / p.size,
                "wilson_lo": lo, "wilson_hi": hi,
                "replicate_rate_min": min(rep_rates), "replicate_rate_max": max(rep_rates),
                "mean_z": float(z.mean()), "mean_phi_excess": float(excess.mean()),
                "median_true_age_w1": float(np.median(
                    [float(r["true_age_w1_a_b0"]) for r in group])),
            })
    return out


def diagnostics(data: dict) -> tuple[list[dict], list[dict], list[dict]]:
    q_rows, age_rows, filter_rows = [], [], []
    for lo, hi in COUNT_BINS:
        sel = (data["derived"] >= lo) & (data["derived"] <= hi)
        if not sel.any():
            continue
        q = data["q"][sel]
        base = {"replicate": data["label"], "derived_count": f"{lo}-{hi}",
                "sites": int(sel.sum())}
        q_rows.append({**base, "mean_support_for_truth": float(q.mean()),
                       "fraction_support_ge_0.5": float(np.mean(q >= 0.5)),
                       "fraction_support_lt_0.1": float(np.mean(q < 0.1))})
        true = data["true_age"][sel]
        age_rows.append({**base,
                         "median_log2_all_over_true": float(np.median(
                             np.log2(data["all_age"][sel] / true))),
                         "median_log2_agree_over_true": float(np.median(
                             np.log2(data["agree_age"][sel] / true)))})
        filter_rows.append({**base, "kept": int(np.sum(q >= 0.5)),
                            "discarded": int(np.sum(q < 0.5))})
    return q_rows, age_rows, filter_rows


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sim-root", type=Path, default=Path(
        "/quobyte/jrigrp/beil/logan_collab/dnAging/msprime_variable_ne_error"))
    p.add_argument("--prep-root", type=Path, default=Path("results/sim_dnaging"),
                   help="per-replicate store, ancestral and draw-polarity tables")
    p.add_argument("--output", type=Path, required=True, help="new output directory")
    p.add_argument("--replicates", type=int, nargs="+", default=list(range(1, 11)))
    p.add_argument("--sizes", type=int, nargs="+", default=[250, 1000, 4000])
    p.add_argument("--tests", type=int, default=50, help="focal sets per replicate and M")
    p.add_argument("--nulls", type=int, default=199, help="null sets per test (R)")
    p.add_argument("--n-bins", type=int, default=20, help="log-age quantile bins")
    p.add_argument("--polarity-seed", type=int, default=2001)
    p.add_argument("--seed", type=int, default=20260927)
    p.add_argument("--alpha", type=float, default=0.05)
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    if args.output.exists():
        raise SystemExit(f"output already exists: {args.output}")
    args.output.mkdir(parents=True)
    started = time.time()
    rng = np.random.default_rng(args.seed)
    fields = ["replicate", "M", "test", "arm", "M_eff", "status", "true_age_w1_a_b0",
              "phi_obs", "null_mean", "null_sd", "z", "p", "mean_daf_difference"]
    inputs, q_rows, age_rows, filter_rows = {}, [], [], []
    tests_path = args.output / "tests.csv"
    with tests_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for replicate in args.replicates:
            data = load_replicate(args.sim_root, args.prep_root, replicate)
            print(f"{data['label']}: {data['q'].size:,} sites, {data['n_draws']} draws, "
                  f"dropped {data['dropped']}, agreeing-age fallback {data['fell_back']}",
                  flush=True)
            inputs[data["label"]] = {"sites": int(data["q"].size), "dropped": data["dropped"],
                                     "n_draws": data["n_draws"], "n_modern": data["n_modern"],
                                     "agreeing_age_fallback": data["fell_back"]}
            q, age, filt = diagnostics(data)
            q_rows += q
            age_rows += age
            filter_rows += filt
            replicate_tests(data, args, rng, writer)
            handle.flush()
    summary = summarize(tests_path, args.alpha)
    write_csv(args.output / "summary.csv", summary)
    write_csv(args.output / "q_calibration.csv", q_rows)
    write_csv(args.output / "age_error.csv", age_rows)
    write_csv(args.output / "filter_by_count.csv", filter_rows)
    (args.output / "run.json").write_text(json.dumps({
        "arguments": {k: (str(v) if isinstance(v, Path) else v) for k, v in vars(args).items()},
        "polarity": POLARITY, "inputs": inputs, "elapsed_seconds": time.time() - started,
    }, indent=2) + "\n")
    for row in summary:
        print(f"{row['arm']:<20} M={row['M']:<5} tests={row['tests']:<4} "
              f"reject={row['rejection_rate']:.3f} [{row['wilson_lo']:.3f},{row['wilson_hi']:.3f}] "
              f"mean_z={row['mean_z']:+.2f}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
