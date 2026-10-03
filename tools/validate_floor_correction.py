"""Check Phi-SFS (W1) floor corrections against a known true distance.

Phi_obs = W1(A, B0) between two finite site sets is biased upward by sampling
noise. This script measures how well four estimates recover the true distance
W1_inf between the two underlying spectra, as a function of set size M:

  raw   Phi_obs
  sub   Phi_obs - mu_0                       (effect size recommended before this check)
  quad  sqrt(max(Phi_obs^2 - mu_0^2, 0))     (PHI_SFS_SAMPLE_SIZE_BIAS.md route 1)
  ext   subsample A and B0 to M/2, M/4, M/8 (levels >= 250), regress the squared
        median Phi on 1/L, sqrt of the intercept    (route 2)

mu_0 is the mean of W1(B_i, B0), i = 1..R, exactly as production computes it.

Two site models, both with n = 26 callable haploids and known polarity:

  dnaging  the 241,623 msprime sites with true mutation age in [1e4, 3e6]
           generations (dnAging msprime_variable_ne_error, replicates 1-10), the
           pool used by PHI_SFS_SAMPLE_SIZE_BIAS.md. B sets are drawn from the
           pool unweighted; A is drawn with an age weight that multiplies the
           rate of sites younger (or older) than an age quantile by a fold change.
  toy      P(d) proportional to (1/d) exp(tilt d/n), as in
           tools/validate_phi_calibration.py; B uses tilt 0.

Sites are drawn i.i.d. (with replacement), so a set is a multinomial vector of
derived counts and its raw spectrum is that vector times the per-count
projection matrix. Matching, depletion and polarity error are not modelled.
W1_inf is exact: W1 between the normalized expected spectra of the two models.

Production code used unchanged: hypergeometric_projection, normalized_spectrum,
phi_sfs (self-checked against the vectorized W1 used in the loops).

Usage:
  python -m tools.validate_floor_correction run --out DIR --workers 8
  python -m tools.validate_floor_correction report --out DIR
"""

from __future__ import annotations

import argparse
import csv
import json
import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from normalize_tes.phi_sfs import (
    PROJECTION_SIZE,
    hypergeometric_projection,
    normalized_spectrum,
    phi_sfs,
)

N = 26
KS = np.arange(1, N)  # derived counts 1..25
DNAGING_SITES = Path(
    "/quobyte/jrigrp/beil/logan_collab/dnAging/msprime_variable_ne_error/"
    "mutation_age_frequency/per_site.tsv"
)
AGE_RANGE = (1e4, 3e6)
M_GRID = (100, 250, 500, 1000, 2500, 5000, 10000, 20000)
R = 200
EXT_MIN_LEVEL = 250
EXT_DRAWS = 400

# (k_rows, 19): retained bins 1..19 of each count's projection to 20.
H = np.stack([hypergeometric_projection(int(k), N)[1:PROJECTION_SIZE] for k in KS])
STEP = 1.0 / PROJECTION_SIZE


def w1_rows(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """W1 between rows of normalized spectra (broadcasting), on the 1/20 grid."""
    return np.abs(np.cumsum(x, axis=-1) - np.cumsum(y, axis=-1))[..., :-1].sum(-1) * STEP


def norm_rows(counts: np.ndarray) -> np.ndarray:
    raw = counts @ H
    return raw / raw.sum(axis=-1, keepdims=True)


def expected_spectrum(p: np.ndarray) -> np.ndarray:
    return normalized_spectrum(p @ H)[1]


# ---------------------------------------------------------------------------
# Models
# ---------------------------------------------------------------------------


def load_dnaging() -> tuple[np.ndarray, np.ndarray]:
    ages, ks = [], []
    with DNAGING_SITES.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            t = float(row["known_mutation_age"])
            if AGE_RANGE[0] <= t <= AGE_RANGE[1]:
                ages.append(t)
                ks.append(int(row["modern_derived_count"]))
    ages, ks = np.asarray(ages), np.asarray(ks)
    if ks.min() < 1 or ks.max() > N - 1:
        raise ValueError("expected derived counts 1..25")
    return ages, ks


def dnaging_probs(ages, ks, direction: str, q: float, fold: float) -> np.ndarray:
    """Count distribution after multiplying the rate of young/old sites by fold."""
    cut = np.quantile(ages, q)
    w = np.ones_like(ages)
    if direction == "recent":
        w[ages < cut] = fold
    elif direction == "old":
        w[ages > cut] = fold
    elif direction != "none":
        raise ValueError(direction)
    p = np.bincount(ks - 1, weights=w, minlength=KS.size)
    return p / p.sum()


def toy_probs(tilt: float) -> np.ndarray:
    p = (1.0 / KS) * np.exp(tilt * KS / N)
    return p / p.sum()


@dataclass
class Effect:
    model: str        # dnaging | toy
    name: str
    direction: str = "none"
    q: float = 0.5
    fold: float = 1.0
    tilt: float = 0.0


def effects() -> list[Effect]:
    out = [Effect("dnaging", "null")]
    for q in (0.25, 0.5):
        for fold in (1.25, 1.5, 2.0, 3.0, 10.0):
            out.append(Effect("dnaging", f"recent_q{q:g}_x{fold:g}", "recent", q, fold))
    for fold in (1.5, 3.0):
        out.append(Effect("dnaging", f"old_q0.75_x{fold:g}", "old", 0.75, fold))
    out.append(Effect("toy", "tilt+0"))
    for tilt in (-0.8, -0.4, -0.2, -0.1, 0.1, 0.2, 0.4, 0.8):
        out.append(Effect("toy", f"tilt{tilt:+g}", tilt=tilt))
    return out


def model_probs(e: Effect, pool) -> tuple[np.ndarray, np.ndarray]:
    if e.model == "dnaging":
        ages, ks = pool
        return (dnaging_probs(ages, ks, e.direction, e.q, e.fold),
                dnaging_probs(ages, ks, "none", e.q, 1.0))
    return toy_probs(e.tilt), toy_probs(0.0)


# ---------------------------------------------------------------------------
# Simulation
# ---------------------------------------------------------------------------


def ext_levels(M: int) -> list[int]:
    return [L for L in (M, M // 2, M // 4, M // 8) if L >= EXT_MIN_LEVEL]


def extrapolate(rng, a_counts, b_counts, M) -> float:
    """Route 2 on one (A, B0) pair; nan when fewer than two levels exist."""
    levels = ext_levels(M)
    if len(levels) < 2:
        return float("nan")
    meds = []
    for L in levels:
        if L == M:
            meds.append(float(w1_rows(norm_rows(a_counts), norm_rows(b_counts))))
            continue
        sa = rng.multivariate_hypergeometric(a_counts, L, size=EXT_DRAWS)
        sb = rng.multivariate_hypergeometric(b_counts, L, size=EXT_DRAWS)
        meds.append(float(np.median(w1_rows(norm_rows(sa), norm_rows(sb)))))
    x = 1.0 / np.asarray(levels, dtype=float)
    slope, intercept = np.polyfit(x, np.square(meds), 1)
    return float(np.sqrt(max(intercept, 0.0)))


def run_config(args) -> dict:
    e_dict, M, analyses, seed, pA, pB = args
    rng = np.random.default_rng(seed)
    pA, pB = np.asarray(pA), np.asarray(pB)
    t0 = time.time()
    rec = {k: np.empty(analyses) for k in ("obs", "mu0", "sd0", "ext")}
    for i in range(analyses):
        a_c = rng.multinomial(M, pA)
        b_c = rng.multinomial(M, pB, size=R + 1)
        a, b = norm_rows(a_c), norm_rows(b_c)
        null = w1_rows(b[1:], b[0])
        rec["obs"][i] = w1_rows(a, b[0])
        rec["mu0"][i] = null.mean()
        rec["sd0"][i] = null.std(ddof=1)
        rec["ext"][i] = extrapolate(rng, a_c, b_c[0], M)
    return dict(effect=e_dict, M=M, seconds=time.time() - t0,
                **{k: v.tolist() for k, v in rec.items()})


def self_check(pool) -> None:
    rng = np.random.default_rng(1)
    for e in effects()[:3] + effects()[-2:]:
        pA, pB = model_probs(e, pool)
        a = norm_rows(rng.multinomial(300, pA, size=5))
        b = norm_rows(rng.multinomial(300, pB, size=5))
        fast = w1_rows(a, b)
        slow = np.array([phi_sfs(x, y).value for x, y in zip(a, b)])
        assert np.allclose(fast, slow, rtol=0, atol=1e-13), (fast, slow)
    print("self-check: vectorized W1 matches normalize_tes.phi_sfs")


def run(out: Path, workers: int, analyses: int, seed: int) -> None:
    out.mkdir(parents=True, exist_ok=True)
    pool = load_dnaging()
    print(f"dnaging pool: {pool[0].size} sites")
    self_check(pool)
    effs = effects()
    truth = {}
    jobs = []
    ss = np.random.SeedSequence(seed)
    children = iter(ss.spawn(len(effs) * len(M_GRID)))
    for e in effs:
        pA, pB = model_probs(e, pool)
        truth[f"{e.model}/{e.name}"] = float(phi_sfs(expected_spectrum(pA),
                                                     expected_spectrum(pB)).value)
        for M in M_GRID:
            jobs.append((asdict(e), M, analyses, next(children), pA.tolist(), pB.tolist()))
    results = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for r in ex.map(run_config, jobs):
            print(f"{r['effect']['model']}/{r['effect']['name']} M={r['M']} "
                  f"{r['seconds']:.1f}s", flush=True)
            results.append(r)
    meta = dict(seed=seed, analyses=analyses, R=R, M_grid=M_GRID, n=N,
                age_range=AGE_RANGE, pool_sites=int(pool[0].size),
                age_quantiles={q: float(np.quantile(pool[0], q)) for q in (0.25, 0.5, 0.75)},
                ext_min_level=EXT_MIN_LEVEL, ext_draws=EXT_DRAWS,
                source=str(DNAGING_SITES), truth=truth)
    (out / "meta.json").write_text(json.dumps(meta, indent=1))
    (out / "results.json").write_text(json.dumps(results))
    report(out)


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------


def report(out: Path) -> None:
    meta = json.loads((out / "meta.json").read_text())
    results = json.loads((out / "results.json").read_text())
    rows = []
    for r in results:
        key = f"{r['effect']['model']}/{r['effect']['name']}"
        w = meta["truth"][key]
        obs, mu0 = np.asarray(r["obs"]), np.asarray(r["mu0"])
        est = dict(
            raw=obs,
            sub=obs - mu0,
            quad=np.sqrt(np.maximum(obs**2 - mu0**2, 0.0)),
            ext=np.asarray(r["ext"]),
        )
        row = dict(model=r["effect"]["model"], effect=r["effect"]["name"], M=r["M"],
                   w1_true=w, floor=float(mu0.mean()),
                   M_floor2=float(r["M"] * mu0.mean() ** 2),
                   quad_clipped=float(np.mean(obs < mu0)))
        for name, v in est.items():
            v = v[np.isfinite(v)]
            if v.size == 0:
                row.update({f"{name}_mean": np.nan, f"{name}_bias": np.nan,
                            f"{name}_rmse": np.nan})
                continue
            row[f"{name}_mean"] = float(v.mean())
            row[f"{name}_bias"] = float(v.mean() - w)
            row[f"{name}_rmse"] = float(np.sqrt(np.mean((v - w) ** 2)))
        rows.append(row)
    rows.sort(key=lambda d: (d["model"], d["effect"], d["M"]))
    with (out / "summary.tsv").open("w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        wr.writeheader()
        for row in rows:
            wr.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in row.items()})
    print(f"wrote {out / 'summary.tsv'} ({len(rows)} rows)")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--out", type=Path, required=True)
    r.add_argument("--workers", type=int, default=1)
    r.add_argument("--analyses", type=int, default=400)
    r.add_argument("--seed", type=int, default=20261002)
    p = sub.add_parser("report")
    p.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    if a.cmd == "run":
        run(a.out, a.workers, a.analyses, a.seed)
    else:
        report(a.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
