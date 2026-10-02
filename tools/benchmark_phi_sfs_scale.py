#!/usr/bin/env python3
"""Benchmark the Python-object-heavy parts of ``normalize_tes.phi_sfs.calculate``.

Round-10 code review Finding 7 flags that ``run_phi_sfs.sbatch``'s ``--mem=48G``
``--time=06:00:00`` rest on an unmeasured guess ("10-12 GB, not measured") for
``calculate()`` at production scale: about 1201 matched sets (the null
replicates plus disjoint spares) of M requested control sites each, alongside
the M focal-set (A) sites.

This script does not call ``calculate()`` or the CLI (both are being edited
concurrently elsewhere in this checkout, and the CLI may start requiring extra
metadata mid-run). Instead it reproduces, using the real functions imported
from ``normalize_tes.phi_sfs``, the sequence of memory-relevant steps that
``calculate()`` performs once the VCF scan (out of scope here -- that part is
I/O-bound and already measured separately) has produced a ``counts`` dict:

  1. synthesize coordinates: a chromosome-label pool of ~10 chromosomes and
     globally unique integer positions, split into R matched-control sets of
     M sites each (``snp_chrom``/``snp_pos``, shape (R, M)) and one focal set
     of M sites (``te_chrom``/``te_pos``), exactly as ``_load_coordinates``
     returns them (``labels[codes]`` fancy indexing for the chromosome array).
  2. build ``all_snp_coordinates`` and ``te_coordinates`` the same way
     ``calculate()`` does: ``list(zip(snp_chrom[row].tolist(),
     snp_pos[row].tolist()))`` per row.
  3. build the ``requested`` set the same way: ``set(te_coordinates)`` unioned
     with every matched set's coordinates.
  4. build a ``counts`` dict of real ``phi_sfs.SiteCount`` instances, one per
     requested site, with synthetic but realistic-range ``alt``/``callable``
     draws and ``p_alt_derived`` drawn from a small discrete grid (``k / 75``)
     so that ``project_sites``'s (k, n)-pair caching behaves the way it does
     on real ARG-derived posteriors instead of degenerating to "everything
     distinct" or "everything collapses".
  5. call the real ``project_sites`` and, per set (A and every B), the real
     ``accumulate_spectrum`` and ``normalized_spectrum``.
  6. build the (R, M) int64 ``bootstrap_counts`` array (each row summing to M,
     matching the invariant ``calculate()`` checks).
  7. reproduce the Wasserstein null-vector computation (CDF residuals against
     one reference replicate) using the real ``phi_sfs`` distance function and
     ``calibrate_phi``.

Nothing is freed between stages (no ``del``, no ``gc.collect()``), matching
``calculate()``, which holds every one of these structures simultaneously
until the function returns. Peak RSS is therefore cumulative across stages,
and is sampled via ``resource.getrusage(RUSAGE_SELF).ru_maxrss`` after each
stage (a running maximum on Linux, so the value never decreases and the final
sample is the whole run's peak). Wall time is measured per stage.

What this deliberately does NOT measure:
  - The VCF scan itself (``read_site_counts``): I/O-bound, out of scope here,
    and already measured separately per the task that requested this script.
  - Whether a real ``counts`` dict built from actual VCF genotypes and a real
    ancestral-state table has the same (k, n)-pair duplication structure as
    the synthetic draws here. Real data may collapse to fewer or more
    distinct rows in ``project_sites``, which changes its output size only
    (not the dominant ``counts``-dict and coordinate-list costs measured
    here).

Run only through the project's ``hpc_run`` helper; this allocates real
production-scale arrays (tens of millions of Python tuple/dataclass objects)
and must not run on a login node.

Example:
    HPC_MEM=64G HPC_TIME=01:00:00 ~/.claude/bin/hpc_run \\
        'python tools/benchmark_phi_sfs_scale.py --M 19000 --R 1201 \\
         --output /tmp/.../bench_M19000.json'
"""
from __future__ import annotations

import argparse
import json
import platform
import resource
import sys
import time
from pathlib import Path

import numpy as np

from normalize_tes.phi_sfs import (
    PROJECTION_SIZE,
    RETAINED_BINS,
    SiteCount,
    accumulate_spectrum,
    calibrate_phi,
    normalized_spectrum,
    project_sites,
)
from normalize_tes.phi_sfs import phi_sfs as compute_phi_sfs


def peak_rss_kb() -> int:
    """Return the process's maximum resident set size so far, in KiB.

    ``ru_maxrss`` is cumulative (a running maximum) on Linux, so sampling it
    repeatedly and taking the last value gives the whole run's peak; it never
    needs resetting between stages.
    """
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss


class StageRecorder:
    def __init__(self) -> None:
        self.stages: list[dict[str, object]] = []
        self._t0 = time.perf_counter()

    def record(self, name: str) -> None:
        now = time.perf_counter()
        self.stages.append({
            "stage": name,
            "wall_seconds": now - self._t0,
            "peak_rss_kib": peak_rss_kb(),
        })
        self._t0 = now


def synthesize_coordinates(
    rng: np.random.Generator, *, m_sites: int, n_sets: int, n_chrom: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Build (te_chrom, te_pos, snp_chrom, snp_pos) with globally unique sites.

    Mirrors ``_load_coordinates``'s output shapes and dtypes: ``te_chrom`` is
    a 1-D string array (like ``.astype(str)`` of a chromosome-name array),
    ``snp_chrom`` is ``labels[codes]`` fancy-indexed from a small label pool
    (real chromosome counts are ~10), and both position arrays are int64.

    Positions are assigned as one running counter partitioned round-robin
    across chromosomes, which trivially guarantees every (chrom, pos) pair
    across the whole run -- the M focal sites plus all R * M matched-control
    sites -- is unique. That matches the disjoint-replicate invariant the
    matcher enforces (``maximum_control_reuse == 1``): a real production bundle
    at R ~= 1201, M = 19000 draws almost the entire ~23M-candidate pool, so
    real cross-set duplication is negligible too. Realistic genomic spacing is
    not needed here: only the object count and hashing cost of (str, int)
    tuples matter for this benchmark, not genomic distances.
    """
    total = m_sites * (n_sets + 1)
    labels = np.array([f"chr{i + 1}" for i in range(n_chrom)])
    all_codes = (np.arange(total, dtype=np.int64) % n_chrom)
    all_pos = (np.arange(total, dtype=np.int64) // n_chrom) + 1

    order = rng.permutation(total)
    all_codes = all_codes[order]
    all_pos = all_pos[order]

    te_codes, te_pos = all_codes[:m_sites], all_pos[:m_sites]
    rest_codes = all_codes[m_sites:].reshape(n_sets, m_sites)
    rest_pos = all_pos[m_sites:].reshape(n_sets, m_sites)

    te_chrom = labels[te_codes].astype(str)
    snp_chrom = labels[rest_codes].astype(str)
    return te_chrom, te_pos.astype(np.int64), snp_chrom, rest_pos.astype(np.int64)


def build_counts(
    rng: np.random.Generator,
    requested: set,
    *,
    callable_min: int,
    callable_max: int,
    p_draws: int,
) -> dict:
    """Build a ``{coordinate: SiteCount}`` dict of the real dataclass.

    ``callable`` is drawn uniformly in ``[callable_min, callable_max]`` (a
    stand-in for realistic inbred-panel callable-individual counts; real
    values come from actual genotype calls and were not measured for this
    benchmark). ``alt`` is drawn uniformly in ``[0, callable]``.
    ``p_alt_derived`` is drawn from the discrete grid ``{k / p_draws : k in
    [0, p_draws]}``, matching the task's request that it come from a small
    value set so ``project_sites``'s (k, n)-pair de-duplication behaves
    realistically rather than either fully collapsing or never colliding.
    """
    n = len(requested)
    callable_vals = rng.integers(callable_min, callable_max + 1, size=n)
    alt_vals = np.floor(rng.random(n) * (callable_vals + 1)).astype(np.int64)
    alt_vals = np.minimum(alt_vals, callable_vals)
    k_vals = rng.integers(0, p_draws + 1, size=n)
    p_vals = k_vals.astype(np.float64) / p_draws

    counts: dict = {}
    for index, coordinate in enumerate(requested):
        counts[coordinate] = SiteCount(
            int(alt_vals[index]), int(callable_vals[index]), float(p_vals[index])
        )
    return counts


def run(args: argparse.Namespace) -> dict:
    rng = np.random.default_rng(args.seed)
    recorder = StageRecorder()
    recorder.record("start")

    te_chrom, te_pos, snp_chrom, snp_pos = synthesize_coordinates(
        rng, m_sites=args.m, n_sets=args.r, n_chrom=args.n_chrom,
    )
    recorder.record("synthesize_coordinates(setup, not part of calculate())")

    # Mirrors calculate() lines building all_snp_coordinates / te_coordinates.
    te_coordinates = list(zip(te_chrom.tolist(), te_pos.tolist()))
    all_snp_coordinates = [
        list(zip(snp_chrom[row].tolist(), snp_pos[row].tolist()))
        for row in range(snp_pos.shape[0])
    ]
    recorder.record("build_coordinate_lists")

    requested: set = set(te_coordinates)
    for row in all_snp_coordinates:
        requested.update(row)
    recorder.record("build_requested_set")

    counts = build_counts(
        rng, requested,
        callable_min=args.callable_min, callable_max=args.callable_max,
        p_draws=args.p_draws,
    )
    recorder.record("build_counts_dict")

    site_rows, projections, endpoints = project_sites(counts)
    recorder.record("project_sites")

    a_counts, a_endpoint, a_eligible = accumulate_spectrum(
        te_coordinates, site_rows, projections, endpoints
    )
    a_raw, a_normalized = normalized_spectrum(a_counts)

    n_sets = len(all_snp_coordinates)
    b_raw = np.empty((n_sets, PROJECTION_SIZE - 1), dtype=np.float64)
    b_normalized = np.empty_like(b_raw)
    b_endpoints = np.empty(n_sets, dtype=np.float64)
    for replicate, coordinates in enumerate(all_snp_coordinates):
        counts_vector, endpoint, eligible = accumulate_spectrum(
            coordinates, site_rows, projections, endpoints
        )
        raw, normalized = normalized_spectrum(counts_vector)
        b_raw[replicate] = raw
        b_normalized[replicate] = normalized
        b_endpoints[replicate] = endpoint
    recorder.record("accumulate_spectrum_per_set")

    # Mirrors the (R, M) bootstrap_counts.npy array held alongside everything
    # else; each row sums to M like calculate()'s validated invariant.
    bootstrap_counts = rng.multinomial(
        args.m, np.full(args.m, 1.0 / args.m), size=n_sets
    ).astype(np.int64)
    recorder.record("bootstrap_counts_array")

    daf = RETAINED_BINS.astype(np.float64) / PROJECTION_SIZE
    b_cdf = np.cumsum(b_normalized, axis=1)
    reference_index = args.reference_index
    reference_sfs = b_normalized[reference_index]
    reference_cdf = b_cdf[reference_index]
    observed_result = compute_phi_sfs(a_normalized, reference_sfs, daf=daf)
    null_indices = np.delete(np.arange(n_sets), reference_index)
    null_phi = (
        np.abs(b_cdf[null_indices, :-1] - reference_cdf[:-1]) * np.diff(daf)
    ).sum(axis=1)
    calibration = calibrate_phi(observed_result.value, null_phi)
    recorder.record("wasserstein_null_vector")

    return {
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "args": {
            "m_sites": args.m,
            "n_sets": args.r,
            "n_chrom": args.n_chrom,
            "seed": args.seed,
            "callable_min": args.callable_min,
            "callable_max": args.callable_max,
            "p_draws": args.p_draws,
            "reference_index": args.reference_index,
        },
        "derived": {
            "requested_sites": len(requested),
            "distinct_kn_pairs": int(projections.shape[0]),
            "eligible_sites_in_A": int(a_eligible),
            "observed_phi_sfs": float(observed_result.value),
            "null_mean": float(calibration.null_mean),
            "null_sd": float(calibration.null_sd),
        },
        "stages": recorder.stages,
        "final_peak_rss_kib": recorder.stages[-1]["peak_rss_kib"],
        "final_peak_rss_gib": recorder.stages[-1]["peak_rss_kib"] / (1024 ** 2),
        "total_wall_seconds": sum(s["wall_seconds"] for s in recorder.stages),
    }


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--m", type=int, required=True,
                         help="M: sites per set (e.g. 4067 in-gene, 19000 genome-wide)")
    parser.add_argument("--r", type=int, default=1201,
                         help="number of matched control sets (default: 1201)")
    parser.add_argument("--n-chrom", type=int, default=10)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--callable-min", type=int, default=20,
                         help="minimum callable-individual draw (>= PROJECTION_SIZE)")
    parser.add_argument("--callable-max", type=int, default=150)
    parser.add_argument("--p-draws", type=int, default=75,
                         help="p_alt_derived is drawn from {k/p_draws : k in [0, p_draws]}")
    parser.add_argument("--reference-index", type=int, default=0)
    parser.add_argument("--output", type=Path, default=None,
                         help="write the JSON report here instead of stdout")
    return parser.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    if args.callable_min < PROJECTION_SIZE:
        raise ValueError(
            f"--callable-min must be >= PROJECTION_SIZE ({PROJECTION_SIZE}) or every "
            "site will be filtered out by project_sites and eligibility checks will fail"
        )
    report = run(args)
    text = json.dumps(report, indent=2)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
        print(f"wrote {args.output}", flush=True)
    print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
