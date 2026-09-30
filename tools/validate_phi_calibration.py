"""Simulation validation of the calibrated Wasserstein Phi-SFS test (plan §9.5).

This is independent of the real data. Sites are drawn from explicit generative
models, projected with the production `project_sites` (hypergeometric projection
to 20 plus the SNP q-mixture), normalized with `normalized_spectrum`, scored with
`phi_sfs` and calibrated with `calibrate_phi`. Only those pure functions are
used; the matcher, VCF reading and bundle code are not exercised.

Sampling shortcut (exact, not an approximation): every site is one of a finite
set of types (k, n, q) with q on a 0.05 grid. A set of M sites drawn i.i.d. from
a type distribution is therefore a multinomial type-count vector, and its raw
spectrum is that count vector times the per-type projection matrix. This is the
same distribution as drawing M sites one by one from an infinite pool. It does
not model a finite pool depleted without replacement (plan §2.2 depletion
check); that needs the real matcher and is out of scope.

Parts (see docs/PHI_SFS_CALIBRATION_VALIDATION.md):
  null   exchangeable neutral model; P-value uniformity, type-I, Z, null mean vs M
  asym   stylised age-covariate model of the bootstrap-target asymmetry
  polar  stylised TE (weight-1, fraction eps flipped) vs SNP (q-mixture) polarity
  power  exponential-tilt rare-variant / high-frequency shifts in A

Usage:
  python -m tools.validate_phi_calibration pilot --out DIR
  python -m tools.validate_phi_calibration run --parts null,asym,polar,power \
      --out DIR --workers 8
  python -m tools.validate_phi_calibration report --out DIR
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
import zlib
from dataclasses import dataclass, field, asdict
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from statistics import NormalDist

from normalize_tes.phi_sfs import (
    PROJECTION_SIZE,
    SiteCount,
    calibrate_phi,
    normalized_spectrum,
    phi_sfs,
    project_sites,
)

Q_GRID = np.round(np.linspace(0.0, 1.0, 21), 10)  # polarity weights on a 0.05 grid
C_GRID = Q_GRID[Q_GRID >= 0.5]                    # confidence max(q, 1-q)
DAF = np.arange(1, PROJECTION_SIZE) / PROJECTION_SIZE


# ---------------------------------------------------------------------------
# Small statistics helpers (the project environment has no scipy)
# ---------------------------------------------------------------------------


def _gammainc_upper_reg(a: float, x: float) -> float:
    """Regularized upper incomplete gamma Q(a, x) (Numerical Recipes gser/gcf)."""
    if x <= 0:
        return 1.0
    gln = math.lgamma(a)
    if x < a + 1:
        ap, total, delta = a, 1.0 / a, 1.0 / a
        for _ in range(1000):
            ap += 1
            delta *= x / ap
            total += delta
            if abs(delta) < abs(total) * 1e-15:
                break
        return 1.0 - total * math.exp(-x + a * math.log(x) - gln)
    b = x + 1 - a
    c, d = 1e300, 1 / b
    h = d
    for i in range(1, 1000):
        an = -i * (i - a)
        b += 2
        d = an * d + b
        d = 1e-300 if abs(d) < 1e-300 else d
        c = b + an / c
        c = 1e-300 if abs(c) < 1e-300 else c
        d = 1 / d
        delta = d * c
        h *= delta
        if abs(delta - 1) < 1e-15:
            break
    return math.exp(-x + a * math.log(x) - gln) * h


def chi2_sf(x: float, df: int) -> float:
    return _gammainc_upper_reg(df / 2.0, x / 2.0)


def chisquare_uniform_p(counts: np.ndarray) -> float:
    counts = np.asarray(counts, dtype=float)
    e = counts.sum() / counts.size
    return chi2_sf(float(((counts - e) ** 2 / e).sum()), counts.size - 1)


def ks_uniform(x: np.ndarray) -> tuple[float, float]:
    """One-sample KS against U(0,1); asymptotic Kolmogorov P with Stephens' correction."""
    x = np.sort(np.asarray(x, dtype=float))
    n = x.size
    i = np.arange(1, n + 1)
    d = float(max(np.max(i / n - x), np.max(x - (i - 1) / n)))
    t = (math.sqrt(n) + 0.12 + 0.11 / math.sqrt(n)) * d
    if t < 0.27:
        return d, 1.0
    p = 2 * sum((-1) ** (k - 1) * math.exp(-2 * k * k * t * t) for k in range(1, 101))
    return d, float(min(max(p, 0.0), 1.0))


def _binom_cdf(x: int, n: int, p: float) -> float:
    if p <= 0:
        return 1.0
    if p >= 1:
        return 1.0 if x >= n else 0.0
    lp, lq = math.log(p), math.log1p(-p)
    terms = [math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1) + k * lp
             + (n - k) * lq for k in range(0, x + 1)]
    m = max(terms)
    return min(1.0, math.exp(m) * sum(math.exp(t - m) for t in terms))


def _bisect(f, lo=0.0, hi=1.0, it=100):
    for _ in range(it):
        mid = (lo + hi) / 2
        if f(mid):
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def skewness(v: np.ndarray) -> float:
    v = np.asarray(v, dtype=float)
    c = v - v.mean()
    return float(np.mean(c ** 3) / np.mean(c ** 2) ** 1.5)


NORM = NormalDist()

# ---------------------------------------------------------------------------
# Generative model pieces
# ---------------------------------------------------------------------------


def confidence_distribution(high_frac: float = 0.85, beta_a: float = 40.0) -> np.ndarray:
    """Discretized distribution of ARG confidence c = max(q, 1 - q) on C_GRID.

    Mixture: `high_frac` of sites c ~ Beta(beta_a, 1) restricted to [0.5, 1]
    (confident, mean about 0.976), remainder c ~ Uniform(0.5, 1). Each grid point
    takes the mass of its +/-0.025 cell (clipped to [0.5, 1]).
    """
    edges = np.concatenate([[0.5], (C_GRID[:-1] + C_GRID[1:]) / 2, [1.0]])
    # Beta(a, 1) CDF is x**a; restrict to [0.5, 1] by renormalizing.
    b_mass = np.diff(edges ** beta_a)
    b_mass /= b_mass.sum()
    u_mass = np.diff(edges) / 0.5
    probs = high_frac * b_mass + (1.0 - high_frac) * u_mass
    return probs / probs.sum()


def neutral_d(n: int, tilt: float = 0.0) -> np.ndarray:
    """P(true derived count d), d = 1..n-1: standard neutral 1/d, times exp(tilt d/n)."""
    d = np.arange(1, n)
    w = (1.0 / d) * np.exp(tilt * d / n)
    return w / w.sum()


def snp_types(n: int, tilt: float = 0.0, *, alt_called_derived: float = 0.7,
              conf: np.ndarray | None = None) -> dict[tuple[int, int, float], float]:
    """Joint probability of observed SNP types (k_alt, n, q_alt).

    Generative story (q is a calibrated posterior):
      d ~ neutral_d(n, tilt); c ~ confidence_distribution();
      with prob `alt_called_derived` q_alt = c, else q_alt = 1 - c;
      ALT is truly derived with prob q_alt, so k_alt = d or n - d.
    """
    conf = confidence_distribution() if conf is None else conf
    pd = neutral_d(n, tilt)
    out: dict[tuple[int, int, float], float] = {}
    for c, pc in zip(C_GRID, conf):
        for q, ps in ((c, alt_called_derived), (round(1.0 - c, 10), 1.0 - alt_called_derived)):
            if c == 0.5 and q != c:
                q = 0.5
            for d, p in zip(range(1, n), pd):
                base = p * pc * ps
                for k, pt in ((d, q), (n - d, 1.0 - q)):
                    if pt <= 0:
                        continue
                    key = (k, n, float(q))
                    out[key] = out.get(key, 0.0) + base * pt
    return out


def te_types(n: int, eps: float, tilt: float = 0.0) -> dict[tuple[int, int, float], float]:
    """TE-like types: weight-1 polarity, a fraction eps truly flipped.

    Insertion presence count k is the observed derived count; with prob eps the
    insertion is actually ancestral so k = n - d.
    """
    pd = neutral_d(n, tilt)
    out: dict[tuple[int, int, float], float] = {}
    for d, p in zip(range(1, n), pd):
        for k, w in ((d, 1.0 - eps), (n - d, eps)):
            if w > 0:
                key = (k, n, 1.0)
                out[key] = out.get(key, 0.0) + p * w
    return out


def mix_over_n(builder, ns, **kw) -> dict:
    out: dict = {}
    for n in ns:
        for key, p in builder(n, **kw).items():
            out[key] = out.get(key, 0.0) + p / len(ns)
    return out


class TypeSpace:
    """Union of type keys with their production projection rows."""

    def __init__(self, dists: list[dict]):
        keys = sorted({key for dist in dists for key in dist})
        self.keys = keys
        self.index = {key: i for i, key in enumerate(keys)}
        counts = {("t", i): SiteCount(k, n, q) for i, (k, n, q) in enumerate(keys)}
        rows, proj, _ = project_sites(counts)
        if len(rows) != len(keys):
            raise RuntimeError("some types have n < 20")
        self.P = np.stack([proj[rows[("t", i)]] for i in range(len(keys))])

    def probs(self, dist: dict) -> np.ndarray:
        p = np.zeros(len(self.keys))
        for key, v in dist.items():
            p[self.index[key]] += v
        return p / p.sum()

    def expected_spectrum(self, p: np.ndarray) -> np.ndarray:
        return normalized_spectrum(p @ self.P)[1]


def draw_counts(rng, M: int, p: np.ndarray, size: int, cdf: np.ndarray | None = None):
    """(size, T) type counts for `size` sets of M i.i.d. sites."""
    T = p.size
    if M < T // 4:
        cdf = np.cumsum(p) if cdf is None else cdf
        cdf = cdf / cdf[-1]
        idx = np.searchsorted(cdf, rng.random((size, M)), side="right")
        idx = np.minimum(idx, T - 1)
        out = np.zeros((size, T), dtype=np.int64)
        for s in range(size):
            out[s] = np.bincount(idx[s], minlength=T)
        return out
    return rng.multinomial(M, p, size=size)


def spectra(counts: np.ndarray, P: np.ndarray) -> np.ndarray:
    raw = counts @ P
    return np.stack([normalized_spectrum(r)[1] for r in raw])


def calibrate_sets(a: np.ndarray, b: np.ndarray) -> dict:
    """a: (19,), b: (R+1, 19) with b[0] = B0. Returns summary of one analysis."""
    b0 = b[0]
    obs = phi_sfs(a, b0)
    null = np.array([phi_sfs(bi, b0).value for bi in b[1:]])
    cal = calibrate_phi(obs.value, null)
    return dict(
        p=cal.p_value, z=cal.z_score, d_obs=cal.observed, null_mean=cal.null_mean,
        null_sd=cal.null_sd, mean_daf_diff=obs.mean_daf_difference,
        null_skew=skewness(null),
    )


# ---------------------------------------------------------------------------
# Configurations
# ---------------------------------------------------------------------------


@dataclass
class Config:
    name: str
    part: str
    M: int
    R: int
    analyses: int
    n_mode: str = "26"            # "26" or "20-40"
    a_model: str = "snp"          # snp | te
    eps: float = 0.0              # te flip fraction
    tilt: float = 0.0             # A tilt (power)
    age_strength: float = 0.0     # asym: s
    sigma: float = 0.0            # asym: matching perturbation (independent per replicate)
    bias: float = 0.0             # asym: systematic matching bias shared by all B_r (older if > 0)
    bootstrap: bool = True        # asym: bootstrap A ages for each B_r
    chunk: int = 50
    extra: dict = field(default_factory=dict)


def ns_for(cfg: Config):
    return [26] if cfg.n_mode == "26" else list(range(20, 41))


G_AGE = 20


def age_z():
    return np.linspace(-0.5, 0.5, G_AGE)


def focal_age_probs():
    """A's age distribution over G_AGE bins: skewed young (stylised TE-like)."""
    g = np.arange(G_AGE)
    w = np.exp(-g / 6.0)
    return w / w.sum()


def build_model(cfg: Config):
    """Return (space, pA, pB) or for asym (space, PgA (G,T))."""
    ns = ns_for(cfg)
    if cfg.part == "asym":
        dists = [mix_over_n(snp_types, ns, tilt=cfg.age_strength * z) for z in age_z()]
        space = TypeSpace(dists)
        return space, np.stack([space.probs(d) for d in dists])
    b = mix_over_n(snp_types, ns)
    if cfg.a_model == "te":
        a = mix_over_n(te_types, ns, eps=cfg.eps, tilt=cfg.tilt)
    else:
        a = mix_over_n(snp_types, ns, tilt=cfg.tilt)
    space = TypeSpace([a, b])
    return space, space.probs(a), space.probs(b)


def round_to_total(x: np.ndarray, M: int) -> np.ndarray:
    """Largest-remainder rounding of each row of x (nonneg) to integers summing to M."""
    scaled = x / x.sum(axis=1, keepdims=True) * M
    fl = np.floor(scaled).astype(np.int64)
    short = M - fl.sum(axis=1)
    rem = scaled - fl
    order = np.argsort(-rem, axis=1)
    for i in range(x.shape[0]):
        fl[i, order[i, : short[i]]] += 1
    return fl


def w1_hist(h1: np.ndarray, h2: np.ndarray) -> np.ndarray:
    """W1 (in age-bin units) between rows of two histogram arrays."""
    c1 = np.cumsum(h1 / h1.sum(axis=-1, keepdims=True), axis=-1)
    c2 = np.cumsum(h2 / h2.sum(axis=-1, keepdims=True), axis=-1)
    return np.abs(c1 - c2)[..., :-1].sum(axis=-1)


def run_chunk(args):
    cfg_dict, chunk_id, seed = args
    cfg = Config(**cfg_dict)
    rng = np.random.default_rng(np.random.SeedSequence([seed, chunk_id]))
    n_here = min(cfg.chunk, cfg.analyses - chunk_id * cfg.chunk)
    out: list[dict] = []
    t0 = time.time()
    if cfg.part == "asym":
        space, Pg = build_model(cfg)
        fA = focal_age_probs()
        for _ in range(n_here):
            hA = rng.multinomial(cfg.M, fA)
            a_counts = sum(rng.multinomial(hA[g], Pg[g]) for g in range(G_AGE))
            a = spectra(a_counts[None, :], space.P)[0]
            S = cfg.R + 1
            if cfg.bootstrap:
                target = rng.multinomial(cfg.M, hA / cfg.M, size=S)
            else:
                target = np.repeat(hA[None, :], S, axis=0)
            if cfg.sigma > 0 or cfg.bias != 0:
                w = target * np.exp(cfg.sigma * rng.standard_normal(target.shape)
                                    + cfg.bias * age_z()[None, :])
                matched = round_to_total(w, cfg.M)
            else:
                matched = target
            b_counts = np.zeros((S, space.P.shape[0]), dtype=np.int64)
            for g in range(G_AGE):
                b_counts += rng.multinomial(matched[:, g], Pg[g])
            b = spectra(b_counts, space.P)
            rec = calibrate_sets(a, b)
            V = w1_hist(target, np.repeat(hA[None, :], S, axis=0))
            E = w1_hist(matched, target)
            with np.errstate(divide="ignore", invalid="ignore"):
                ratio = np.where(V > 0, E / V, np.nan)
            rec["E_over_V_median"] = (float(np.nanmedian(ratio)) if np.any(np.isfinite(ratio))
                                      else float("nan"))
            out.append(rec)
    else:
        space, pA, pB = build_model(cfg)
        cA, cB = np.cumsum(pA), np.cumsum(pB)
        for _ in range(n_here):
            a = spectra(draw_counts(rng, cfg.M, pA, 1, cA), space.P)[0]
            b = spectra(draw_counts(rng, cfg.M, pB, cfg.R + 1, cB), space.P)
            out.append(calibrate_sets(a, b))
    return cfg.name, chunk_id, out, time.time() - t0


# ---------------------------------------------------------------------------
# Config sets
# ---------------------------------------------------------------------------


def configs_for(part: str, scale: float = 1.0) -> list[Config]:
    def k(n):  # analyses, scaled for pilot
        return max(2, int(round(n * scale)))

    cfgs: list[Config] = []
    if part == "null":
        for M in (500, 1000, 4000, 20000):
            cfgs.append(Config(f"null_n26_M{M}_R200", "null", M, 200, k(2000)))
        for M in (500, 1000, 4000, 20000):
            cfgs.append(Config(f"null_n20-40_M{M}_R200", "null", M, 200, k(1000), n_mode="20-40"))
        for M in (500, 4000, 20000):
            cfgs.append(Config(f"null_n26_M{M}_R1000", "null", M, 1000, k(1000)))
    elif part == "asym":
        for s in (4.0, 12.0):
            cfgs.append(Config(f"asym_s{s:g}_obsTarget_sig0", "asym", 4000, 200, k(1000),
                               age_strength=s, bootstrap=False))
            for sig in (0.0, 0.01, 0.02, 0.035, 0.07, 0.15):
                cfgs.append(Config(f"asym_s{s:g}_boot_sig{sig:g}", "asym", 4000, 200, k(1000),
                                   age_strength=s, sigma=sig))
        for bias in (0.05, 0.1, 0.2, 0.4):
            cfgs.append(Config(f"asym_s12_boot_sig0.01_bias{bias:g}", "asym", 4000, 200, k(1000),
                               age_strength=12.0, sigma=0.01, bias=bias))
        cfgs.append(Config("asym_s12_boot_sig0_M1000", "asym", 1000, 200, k(1000),
                           age_strength=12.0))
        cfgs.append(Config("asym_s12_boot_sig0_M20000", "asym", 20000, 200, k(1000),
                           age_strength=12.0))
    elif part == "polar":
        for M in (1000, 4000, 20000):
            for eps in (0.0, 0.03, 0.0827, 0.09):
                cfgs.append(Config(f"polar_M{M}_eps{eps:g}", "polar", M, 200, k(1000),
                                   a_model="te", eps=eps))
    elif part == "power":
        for M in (1000, 4000, 20000):
            for tilt in (-0.8, -0.4, -0.2, -0.1, 0.1, 0.2, 0.4, 0.8):
                cfgs.append(Config(f"power_M{M}_tilt{tilt:+g}", "power", M, 200, k(500),
                                   tilt=tilt))
    else:
        raise ValueError(part)
    return cfgs


# ---------------------------------------------------------------------------
# Summaries
# ---------------------------------------------------------------------------


def clopper_pearson(x: int, n: int, level: float = 0.95):
    a = (1 - level) / 2
    # lower: P(X >= x | lo) = a ; upper: P(X <= x | hi) = a
    lo = 0.0 if x == 0 else _bisect(lambda q: 1 - _binom_cdf(x - 1, n, q) < a)
    hi = 1.0 if x == n else _bisect(lambda q: _binom_cdf(x, n, q) > a)
    return float(lo), float(hi)


def summarize(cfg: Config, recs: list[dict], space_info: dict) -> dict:
    p = np.array([r["p"] for r in recs])
    z = np.array([r["z"] for r in recs])
    n = p.size
    s: dict = {"config": asdict(cfg), "n": n}
    s.update(space_info)
    for alpha in (0.05, 0.01):
        x = int(np.count_nonzero(p <= alpha))
        s[f"rej_{alpha}"] = x / n
        s[f"rej_{alpha}_ci"] = clopper_pearson(x, n)
        s[f"exact_level_{alpha}"] = math.floor(alpha * (cfg.R + 1) + 1e-9) / (cfg.R + 1)
        zc = NORM.inv_cdf(1 - alpha)
        s[f"rej_normalZ_{alpha}"] = float(np.mean(z >= zc))
    s["ks_D"], s["ks_p"] = ks_uniform(p)
    dec = np.histogram(p, bins=np.linspace(0, 1, 11))[0]
    s["deciles"] = dec.tolist()
    s["deciles_chi2_p"] = chisquare_uniform_p(dec)
    s["mean_p"] = float(p.mean())
    s["z_mean"], s["z_sd"] = float(z.mean()), float(z.std(ddof=1))
    s["z_mean_se"] = float(z.std(ddof=1) / math.sqrt(n))
    s["z_median"] = float(np.median(z))
    for key in ("d_obs", "null_mean", "null_sd", "null_skew", "mean_daf_diff"):
        v = np.array([r[key] for r in recs])
        s[f"{key}_mean"] = float(v.mean())
        s[f"{key}_sd"] = float(v.std(ddof=1))
    s["d_obs_over_null_mean"] = s["d_obs_mean"] / s["null_mean_mean"]
    s["null_mean_x_sqrtM"] = s["null_mean_mean"] * math.sqrt(cfg.M)
    s["frac_daf_diff_positive"] = float(np.mean([r["mean_daf_diff"] > 0 for r in recs]))
    if "E_over_V_median" in recs[0]:
        ev = np.array([r["E_over_V_median"] for r in recs])
        s["E_over_V_median"] = float(np.nanmedian(ev)) if np.any(np.isfinite(ev)) else float("nan")
    return s


def space_info(cfg: Config) -> dict:
    """Population-level quantities of the model (not simulated)."""
    info: dict = {}
    if cfg.part == "asym":
        space, Pg = build_model(cfg)
        means = [float(space.expected_spectrum(Pg[g]) @ DAF) for g in range(G_AGE)]
        info["mean_daf_youngest"], info["mean_daf_oldest"] = means[0], means[-1]
        fA = focal_age_probs()
        info["mean_daf_A_expected"] = float(space.expected_spectrum(fA @ Pg) @ DAF)
        return info
    space, pA, pB = build_model(cfg)
    ea, eb = space.expected_spectrum(pA), space.expected_spectrum(pB)
    info["expected_W"] = phi_sfs(ea, eb).value
    info["expected_mean_daf_diff"] = phi_sfs(ea, eb).mean_daf_difference
    info["expected_mean_daf_B"] = float(eb @ DAF)
    # SNP mis-orientation mass m_B = sum over types of P(type) * weight on the wrong
    # orientation, computed from the generative model directly.
    conf = confidence_distribution()
    info["snp_misorientation_mass"] = float(np.sum(conf * 2 * C_GRID * (1 - C_GRID)))
    info["snp_wrong_call_rate"] = float(np.sum(conf * (1 - C_GRID)))
    return info


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------


def execute(cfgs: list[Config], out: Path, workers: int, seed: int, tag: str) -> dict:
    out.mkdir(parents=True, exist_ok=True)
    tasks = []
    for ci, cfg in enumerate(cfgs):
        nchunks = math.ceil(cfg.analyses / cfg.chunk)
        for j in range(nchunks):
            tasks.append((asdict(cfg), j, seed * 100_000 + zlib.crc32(cfg.name.encode()) % 100_000))
    # longest first for load balance
    def cost(t):
        c = t[0]
        return c["M"] * (c["R"] + 2) * min(c["chunk"], c["analyses"])
    tasks.sort(key=cost, reverse=True)
    results: dict[str, dict[int, list]] = {c.name: {} for c in cfgs}
    cpu: dict[str, float] = {c.name: 0.0 for c in cfgs}
    t0 = time.time()
    with Pool(workers) as pool:
        for name, j, recs, dt in pool.imap_unordered(run_chunk, tasks):
            results[name][j] = recs
            cpu[name] += dt
    wall = time.time() - t0
    summaries = {}
    for cfg in cfgs:
        recs = [r for j in sorted(results[cfg.name]) for r in results[cfg.name][j]]
        np.savez(out / f"{cfg.name}{tag[len(cfg.part):] if tag.startswith(cfg.part) else ''}.npz", **{k: np.array([r[k] for r in recs]) for k in recs[0]})
        s = summarize(cfg, recs, space_info(cfg))
        s["cpu_seconds"] = cpu[cfg.name]
        summaries[cfg.name] = s
    meta = {"wall_seconds": wall, "workers": workers, "seed": seed,
            "cpu_seconds_total": sum(cpu.values()), "summaries": summaries}
    (out / f"summary_{tag}.json").write_text(json.dumps(meta, indent=1))
    print(f"[{tag}] wall {wall:.1f}s cpu {sum(cpu.values()):.1f}s", flush=True)
    return meta


def fmt_ci(ci):
    return f"[{ci[0]:.4f}, {ci[1]:.4f}]"


def report(out: Path) -> None:
    for f in sorted(out.glob("summary_*.json")):
        meta = json.loads(f.read_text())
        print(f"\n## {f.name}  wall {meta['wall_seconds']:.0f}s  cpu {meta['cpu_seconds_total']:.0f}s")
        for name, s in meta["summaries"].items():
            c = s["config"]
            line = (f"{name}: n={s['n']} rej05={s['rej_0.05']:.4f} {fmt_ci(s['rej_0.05_ci'])} "
                    f"rej01={s['rej_0.01']:.4f} {fmt_ci(s['rej_0.01_ci'])} KS_D={s['ks_D']:.4f} "
                    f"KSp={s['ks_p']:.3f} dec_p={s['deciles_chi2_p']:.3f} z={s['z_mean']:.3f}"
                    f"+-{s['z_mean_se']:.3f} sd={s['z_sd']:.3f} normZrej05={s['rej_normalZ_0.05']:.4f}"
                    f" nullmean={s['null_mean_mean']:.3e} x sqrtM={s['null_mean_x_sqrtM']:.4f}"
                    f" dobs/null={s['d_obs_over_null_mean']:.3f} skew={s['null_skew_mean']:.2f}"
                    f" daf+={s['frac_daf_diff_positive']:.3f} dafdiff={s['mean_daf_diff_mean']:.4f}")
            for key in ("expected_W", "expected_mean_daf_diff", "E_over_V_median",
                        "mean_daf_youngest", "mean_daf_oldest", "snp_misorientation_mass",
                        "snp_wrong_call_rate", "expected_mean_daf_B"):
                if key in s:
                    line += f" {key}={s[key]:.4g}"
            line += f" cpu={s['cpu_seconds']:.0f}s dec={s['deciles']}"
            print(line)


def self_check() -> None:
    """Sanity checks of the type-space shortcut against direct per-site projection."""
    space = TypeSpace([snp_types(26)])
    rng = np.random.default_rng(1)
    p = space.probs(snp_types(26))
    counts = draw_counts(rng, 300, p, 1)[0]
    sites = {}
    i = 0
    for t, c in enumerate(counts):
        k, n, q = space.keys[t]
        for _ in range(c):
            sites[("s", i)] = SiteCount(k, n, q)
            i += 1
    rows, proj, _ = project_sites(sites)
    direct = normalized_spectrum(sum(proj[rows[key]] for key in sites))[1]
    via = spectra(counts[None, :], space.P)[0]
    assert np.allclose(direct, via, atol=1e-12), "type shortcut disagrees with project_sites"
    # model sums to one and the searchsorted/multinomial paths agree in mean
    assert abs(p.sum() - 1) < 1e-12
    m1 = draw_counts(rng, 50, p, 4000).mean(0)
    m2 = rng.multinomial(50, p, size=4000).mean(0)
    se = np.sqrt(2 * 50 * p * (1 - p) / 4000) + 1e-12
    assert np.max(np.abs(m1 - m2) / se) < 6, "sampler paths disagree"
    # statistics helpers against known values
    assert abs(chi2_sf(16.919, 9) - 0.05) < 1e-4
    lo, hi = clopper_pearson(50, 1000)
    assert abs(lo - 0.03736) < 2e-4 and abs(hi - 0.06546) < 2e-4, (lo, hi)
    assert abs(ks_uniform(np.linspace(0.0005, 0.9995, 1000))[1] - 1.0) < 1e-6
    t_crit = 1.3581 / (math.sqrt(1000) + 0.12 + 0.11 / math.sqrt(1000))
    shifted = np.clip(np.linspace(0.0005, 0.9995, 1000) + t_crit - 0.0005, 0, 1)
    assert abs(ks_uniform(shifted)[1] - 0.05) < 2e-3, ks_uniform(shifted)
    print("self-check ok; T(n=26) =", len(space.keys),
          "misorientation mass", float(np.sum(confidence_distribution() * 2 * C_GRID * (1 - C_GRID))))


# ---------------------------------------------------------------------------
# Polarity simulations ("polarity2"): mechanistic ARG draws, TE and SNP weights
# ---------------------------------------------------------------------------
#
# Truth (n = 26 inbred lines, standard neutral P(d) proportional to 1/d):
#   SNP: REF is the reference line's allele (one of the 26 lines), so ALT is
#        derived iff the reference line is ancestral. P(k_ALT = k, ALT derived)
#        is proportional to (1/k)(n-k)/n, P(k_ALT = k, ALT ancestral) to 1/n, so
#        P(ALT derived | k) = (n-k)/n and the ALT-count spectrum is 1/k.
#   TE:  insertion count k. Insertion derived with prob 1-delta (d = k); with
#        prob delta the insertion fixed and was later deleted, so it is the
#        ancestral allele and the derived (deletion) count is d = n-k.
# ARG: per-site evidence e (a log-likelihood ratio for "focal allele derived",
#   focal = ALT for SNPs, insertion for TEs) ~ N(s*mu, 2*mu), s = +1/-1 for the
#   true orientation, mu ~ Gamma(kappa, mean m). Each of T = 75 draws calls the
#   focal allele derived with prob pi = sigmoid(lam*e + L(k)); q = X/75.
#   L(k) = log((n-k)/k) ("aware": the coalescent frequency prior, which makes q
#   a calibrated posterior given (e, k) when lam = 1); L = 0 ("blind-flat");
#   L = L0 = logit P(ALT derived) ("blind-marginal", calibrated given e only,
#   the analogue of the Part 3 model).

N_POL = 26
T_DRAWS = 75
K_POL = np.arange(1, N_POL)


def _sigmoid(x):
    return 1.0 / (1.0 + np.exp(-x))


def _logit(p):
    return np.log(p) - np.log1p(-p)


def pol_truth(kind: str, delta: float = 0.0) -> np.ndarray:
    """P[k-1, s] for k = 1..n-1, s = 1 focal allele derived, s = 0 ancestral."""
    k = K_POL.astype(float)
    n = N_POL
    P = np.zeros((k.size, 2))
    if kind == "snp":
        P[:, 1] = (1.0 / k) * (n - k) / n
        P[:, 0] = 1.0 / n
    elif kind == "te":
        P[:, 1] = (1.0 - delta) / k
        P[:, 0] = delta / (n - k)
    else:
        raise ValueError(kind)
    return P / P.sum()


def snp_marginal_logodds() -> float:
    P = pol_truth("snp")
    return float(_logit(P[:, 1].sum()))


def _gamma_nodes(kappa: float, mean: float, Q: int = 24) -> np.ndarray:
    """Deterministic quantile midpoints of Gamma(kappa, mean/kappa)."""
    rng = np.random.default_rng(12345)
    x = rng.gamma(kappa, mean / kappa, size=400_000)
    return np.quantile(x, (np.arange(Q) + 0.5) / Q)


_GH_X, _GH_W = np.polynomial.hermite.hermgauss(40)
_LOGC = np.array([math.lgamma(T_DRAWS + 1) - math.lgamma(x + 1) - math.lgamma(T_DRAWS - x + 1)
                  for x in range(T_DRAWS + 1)])


def arg_prior(model: dict) -> np.ndarray:
    k = K_POL.astype(float)
    if model["prior"] == "aware":
        return np.log((N_POL - k) / k)
    if model["prior"] == "flat":
        return np.zeros(k.size)
    if model["prior"] == "marginal":
        return np.full(k.size, snp_marginal_logodds())
    raise ValueError(model["prior"])


def arg_x_given_ks(model: dict) -> np.ndarray:
    """P(X = x | k, s) for x = 0..75, shape (n-1, 2, 76)."""
    mu = _gamma_nodes(model["kappa"], model["m"])                  # (Q,)
    L = arg_prior(model)                                           # (K,)
    s = np.array([-1.0, 1.0])                                      # index 0: ancestral
    e = (s[None, :, None, None] * mu[None, None, :, None]
         + 2.0 * np.sqrt(mu)[None, None, :, None] * _GH_X[None, None, None, :])  # (1,2,Q,G)
    lin = model["lam"] * e + L[:, None, None, None]               # (K,2,Q,G)
    pi = np.clip(_sigmoid(lin), 1e-15, 1 - 1e-15)
    x = np.arange(T_DRAWS + 1)
    logp = (_LOGC + x * np.log(pi)[..., None] + (T_DRAWS - x) * np.log1p(-pi)[..., None])
    w = (_GH_W / math.sqrt(math.pi))[None, None, None, :, None] / mu.size
    return (np.exp(logp) * w).sum(axis=(2, 3))                    # (K,2,76)


def pol_joint(kind: str, model: dict, delta: float = 0.0) -> np.ndarray:
    """P(k, s, X), shape (n-1, 2, 76)."""
    return pol_truth(kind, delta)[:, :, None] * arg_x_given_ks(model)


def _dist_from(weights_kx: np.ndarray, prob_kx: np.ndarray) -> dict:
    out: dict = {}
    for i, k in enumerate(K_POL):
        for x in range(T_DRAWS + 1):
            p = prob_kx[i, x]
            if p <= 0:
                continue
            key = (int(k), N_POL, float(round(weights_kx[i, x], 12)))
            out[key] = out.get(key, 0.0) + p
    tot = sum(out.values())
    return {k: v / tot for k, v in out.items()}


def snp_dist(model: dict, rule: str = "mean") -> dict:
    """SNP site types (k_ALT, n, weight). rule: mean (production) | impute."""
    pk = pol_joint("snp", model).sum(axis=1)                       # (K,76)
    q = np.arange(T_DRAWS + 1) / T_DRAWS
    if rule == "mean":
        return _dist_from(np.broadcast_to(q, pk.shape), pk)
    if rule == "impute":  # per-set posterior imputation: hard orientation drawn from q
        out: dict = {}
        for i, k in enumerate(K_POL):
            p1 = float((pk[i] * q).sum()); p0 = float((pk[i] * (1 - q)).sum())
            out[(int(k), N_POL, 1.0)] = p1
            out[(int(k), N_POL, 0.0)] = p0
        tot = sum(out.values())
        return {k: v / tot for k, v in out.items()}
    raise ValueError(rule)


def te_weights(rule: str, delta_assumed: float = 0.03) -> np.ndarray:
    """Weight (prob insertion derived) per (k, X) for TE rules w1 | q | bayes."""
    x = np.arange(T_DRAWS + 1)
    if rule == "w1":
        return np.ones((K_POL.size, x.size))
    if rule == "q":
        return np.broadcast_to(x / T_DRAWS, (K_POL.size, x.size)).copy()
    if rule == "bayes":
        # Strip the ARG's own (coalescent) frequency prior from a smoothed q and
        # replace it with the TE prior 1 - delta_assumed.
        qt = (x + 0.5) / (T_DRAWS + 1)
        L_aware = np.log((N_POL - K_POL) / K_POL)
        return _sigmoid(_logit(qt)[None, :] - L_aware[:, None] + _logit(1 - delta_assumed))
    raise ValueError(rule)


def te_dist(model: dict, delta: float, rule: str, filt: bool = True,
            delta_assumed: float = 0.03) -> dict:
    pk = pol_joint("te", model, delta).sum(axis=1)                 # (K,76)
    if filt:  # retain a TE only if >= 50% of draws call the insertion derived
        pk = pk * (np.arange(T_DRAWS + 1) >= math.ceil(T_DRAWS / 2))[None, :]
    return _dist_from(te_weights(rule, delta_assumed), pk)


def te_filter_stats(model: dict, delta: float) -> dict:
    J = pol_joint("te", model, delta)                              # (K,2,76)
    keep = np.arange(T_DRAWS + 1) >= math.ceil(T_DRAWS / 2)
    kept = J[:, :, keep].sum(axis=2)
    tot = J.sum(axis=2)
    return dict(retained=float(kept.sum()),
                deleted_retained=float(kept[:, 0].sum() / max(tot[:, 0].sum(), 1e-300)),
                derived_retained=float(kept[:, 1].sum() / tot[:, 1].sum()),
                retained_by_k=(kept.sum(1) / tot.sum(1)).tolist())


def true_spectrum() -> np.ndarray:
    """Normalized projected spectrum of the true derived counts, d ~ 1/d, n = 26."""
    d = K_POL.astype(float)
    pd = (1 / d) / (1 / d).sum()
    space = TypeSpace([{(int(k), N_POL, 1.0): float(p) for k, p in zip(K_POL, pd)}])
    return space.expected_spectrum(space.probs({(int(k), N_POL, 1.0): float(p) for k, p in zip(K_POL, pd)}))


def arg_diagnostics(model: dict, delta: float = 0.03) -> dict:
    """Predicted observables comparable with the chr10 real-data check."""
    q = np.arange(T_DRAWS + 1) / T_DRAWS
    out = {}
    for kind in ("snp", "te"):
        J = pol_joint(kind, model, delta if kind == "te" else 0.0)
        pkx = J.sum(axis=1)
        pk = pkx.sum(axis=1)
        out[f"{kind}_mean_q_by_k"] = ((pkx * q).sum(1) / pk).tolist()
        unan = pkx[:, [0, T_DRAWS]].sum()
        out[f"{kind}_frac_unanimous"] = float(unan)
        out[f"{kind}_mean_q"] = float((pkx * q).sum())
    Jt = pol_joint("te", model, delta).sum(axis=1)
    out["te_unanimous_correct"] = float(Jt[:, T_DRAWS].sum() / Jt[:, [0, T_DRAWS]].sum())
    out["te_majority_correct"] = float(Jt[:, math.ceil(T_DRAWS / 2):].sum())
    # Calibration of the SNP q given (k, q): P(ALT derived | k, X) vs X/75
    Js = pol_joint("snp", model)
    post = Js[:, 1, :] / np.maximum(Js.sum(axis=1), 1e-300)
    wts = Js.sum(axis=1)
    out["snp_mean_abs_calibration_gap"] = float((np.abs(post - q[None, :]) * wts).sum())
    return out


REAL_TE_Q_BY_K = [0.959, 0.94, 0.915, 0.937, 0.899, 0.826, 0.862, 0.841, 0.74, 0.755, 0.799,
                  0.81, 0.637, 0.703, 0.675, 0.703, 0.65, 0.623, 0.574, 0.595, 0.637, 0.303,
                  0.415, 0.279, 0.217]
REAL_TE_N_BY_K = [1077, 452, 261, 175, 95, 103, 80, 80, 46, 52, 58, 35, 28, 34, 27, 37, 24, 44,
                  25, 27, 24, 28, 27, 34, 35]
REAL_SNP_Q_BY_K = [0.787, 0.825, 0.774, 0.829, 0.756, 0.715, 0.693, 0.683, 0.603, 0.587, 0.566,
                   0.495, 0.52, 0.48, 0.431, 0.437, 0.339, 0.373, 0.193, 0.304, 0.225, 0.19,
                   0.241, 0.098, 0.045]
REAL_TE_UNANIMOUS = 0.754
REAL_TE_UNANIMOUS_CORRECT = 0.908


def fit_arg(prior: str, delta: float = 0.03) -> dict:
    """Grid fit of (kappa, m, lam) to the chr10 TE observables (count-weighted)."""
    wk = np.array(REAL_TE_N_BY_K, float); wk /= wk.sum()
    tq = np.array(REAL_TE_Q_BY_K)
    best = None
    for kappa in (0.3, 0.6, 1.0, 2.0):
        for m in np.geomspace(0.5, 60, 18):
            for lam in (0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0):
                model = dict(prior=prior, kappa=kappa, m=float(m), lam=lam)
                d = arg_diagnostics(model, delta)
                sse = float((wk * (np.array(d["te_mean_q_by_k"]) - tq) ** 2).sum()) \
                    + 0.05 * (d["te_frac_unanimous"] - REAL_TE_UNANIMOUS) ** 2 \
                    + 0.05 * (d["te_unanimous_correct"] - REAL_TE_UNANIMOUS_CORRECT) ** 2
                if best is None or sse < best[0]:
                    best = (sse, model, d)
    sse, model, d = best
    return dict(model=model, sse=sse, rmse_te_q=float(np.sqrt((wk * (np.array(d["te_mean_q_by_k"]) - tq) ** 2).sum())),
                diag=d)



# ---- polarity2 simulation driver ------------------------------------------

# ARG models. "fit" values are filled from `polfit` (chr10 TE observables); the
# calibrated / under / over arms share the fitted evidence distribution.
ARG_MODELS: dict[str, dict] = {}

# Best grid fit (polfit, frequency-aware prior, delta = 0.03) to the chr10 TE
# mean-q-by-insertion-count curve, unanimity and unanimous accuracy.
POL_FIT = dict(kappa=0.3, m=8.35637259521861, lam=4.0)


def tune_over_lam(base: dict, target: float = 0.908, delta: float = 0.03) -> float:
    """Evidence scale lam at which TE unanimous-draw accuracy equals `target`."""
    lo, hi = 1.0, 64.0
    for _ in range(30):
        mid = math.sqrt(lo * hi)
        acc = arg_diagnostics(dict(prior="aware", lam=mid, **base), delta)["te_unanimous_correct"]
        if acc > target:
            lo = mid
        else:
            hi = mid
    return float(math.sqrt(lo * hi))


def _register_models(fit: dict) -> None:
    base = dict(kappa=fit["kappa"], m=fit["m"])
    ARG_MODELS.update({
        "aware_cal": dict(prior="aware", lam=1.0, **base),
        "aware_under": dict(prior="aware", lam=0.5, **base),
        "aware_fit": dict(prior="aware", lam=fit["lam"], **base),
        "aware_over": dict(prior="aware", lam=tune_over_lam(base), **base),
        "blind_flat": dict(prior="flat", lam=1.0, **base),
        "blind_marg": dict(prior="marginal", lam=1.0, **base),
    })


@dataclass
class PolConfig:
    name: str
    kind: str                 # te_snp | snp_snp | contrast_te | contrast_snp
    arg: str
    M: int
    R: int
    analyses: int
    delta: float = 0.03
    te_rule: str = "w1"       # w1 | q | bayes
    te_filter: bool = True
    snp_rule: str = "mean"    # mean | impute
    M2: int = 0               # contrast: second category size
    chunk: int = 50


_POL_CACHE: dict = {}


def pol_dists(cfg: PolConfig) -> tuple[dict, dict]:
    key = (cfg.kind, cfg.arg, cfg.delta, cfg.te_rule, cfg.te_filter, cfg.snp_rule)
    if key not in _POL_CACHE:
        model = ARG_MODELS[cfg.arg]
        b = snp_dist(model, cfg.snp_rule)
        if cfg.kind in ("te_snp", "contrast_te"):
            a = te_dist(model, cfg.delta, cfg.te_rule, cfg.te_filter)
        else:
            a = b
        _POL_CACHE[key] = (a, b)
    return _POL_CACHE[key]


def pol_space(cfg: PolConfig):
    a, b = pol_dists(cfg)
    key = ("space",) + (cfg.kind, cfg.arg, cfg.delta, cfg.te_rule, cfg.te_filter, cfg.snp_rule)
    if key not in _POL_CACHE:
        space = TypeSpace([a, b])
        _POL_CACHE[key] = (space, space.probs(a), space.probs(b))
    return _POL_CACHE[key]


def _one_category(rng, space, pA, pB, M, R):
    a = spectra(draw_counts(rng, M, pA, 1), space.P)[0]
    b = spectra(draw_counts(rng, M, pB, R + 1), space.P)
    b0 = b[0]
    obs = phi_sfs(a, b0)
    null = np.array([phi_sfs(bi, b0).value for bi in b[1:]])
    return obs, calibrate_phi(obs.value, null)


def run_pol_chunk(args):
    cfg_dict, chunk_id, seed, fit = args
    if not ARG_MODELS:
        _register_models(fit)
    cfg = PolConfig(**cfg_dict)
    rng = np.random.default_rng(np.random.SeedSequence([seed, chunk_id]))
    n_here = min(cfg.chunk, cfg.analyses - chunk_id * cfg.chunk)
    space, pA, pB = pol_space(cfg)
    out = []
    t0 = time.time()
    for _ in range(n_here):
        if cfg.kind.startswith("contrast"):
            _, c1 = _one_category(rng, space, pA, pB, cfg.M, cfg.R)
            _, c2 = _one_category(rng, space, pA, pB, cfg.M2, cfg.R)
            d_obs = c1.z_score - c2.z_score
            d_null = c1.null_z_scores - c2.null_z_scores[rng.permutation(cfg.R)]
            p = (1 + np.count_nonzero(np.abs(d_null) >= abs(d_obs))) / (cfg.R + 1)
            out.append(dict(p=float(p), z=float(d_obs), z1=c1.z_score, z2=c2.z_score,
                            d_obs=float(abs(d_obs)), null_mean=float(np.mean(np.abs(d_null))),
                            null_sd=float(np.std(d_null, ddof=1)), mean_daf_diff=0.0,
                            null_skew=0.0))
        else:
            obs, cal = _one_category(rng, space, pA, pB, cfg.M, cfg.R)
            out.append(dict(p=cal.p_value, z=cal.z_score, d_obs=cal.observed,
                            null_mean=cal.null_mean, null_sd=cal.null_sd,
                            mean_daf_diff=obs.mean_daf_difference,
                            null_skew=skewness(cal.null)))
    return cfg.name, chunk_id, out, time.time() - t0


def pol_expected(cfg: PolConfig) -> dict:
    """Population-level (not simulated) spectra: bias of A and B vs the truth."""
    space, pA, pB = pol_space(cfg)
    ea, eb = space.expected_spectrum(pA), space.expected_spectrum(pB)
    t = true_spectrum()
    r = phi_sfs(ea, eb)
    info = dict(expected_W=r.value, expected_mean_daf_diff=r.mean_daf_difference,
                W_A_truth=phi_sfs(ea, t).value, W_B_truth=phi_sfs(eb, t).value,
                daf_A_minus_truth=phi_sfs(ea, t).mean_daf_difference,
                daf_B_minus_truth=phi_sfs(eb, t).mean_daf_difference)
    if cfg.kind in ("te_snp", "contrast_te"):
        info.update({f"filter_{k}": v for k, v in te_filter_stats(ARG_MODELS[cfg.arg], cfg.delta).items()
                     if k != "retained_by_k"})
    return info


def pol_summarize(cfg: PolConfig, recs: list[dict]) -> dict:
    p = np.array([r["p"] for r in recs]); z = np.array([r["z"] for r in recs])
    n = p.size
    x = int(np.count_nonzero(p <= 0.05))
    s = dict(config=asdict(cfg), n=n, rej_05=x / n, rej_05_ci=clopper_pearson(x, n),
             exact_level=math.floor(0.05 * (cfg.R + 1) + 1e-9) / (cfg.R + 1),
             z_mean=float(z.mean()), z_sd=float(z.std(ddof=1)),
             z_mean_se=float(z.std(ddof=1) / math.sqrt(n)),
             frac_daf_pos=float(np.mean([r["mean_daf_diff"] > 0 for r in recs])),
             mean_daf_diff=float(np.mean([r["mean_daf_diff"] for r in recs])),
             ks_p=ks_uniform(p)[1])
    x1 = int(np.count_nonzero(p <= 0.01))
    s.update(rej_01=x1 / n, rej_01_ci=clopper_pearson(x1, n))
    if cfg.kind.startswith("contrast"):
        s["z1_mean"] = float(np.mean([r["z1"] for r in recs]))
        s["z2_mean"] = float(np.mean([r["z2"] for r in recs]))
    return s


def pol_configs(stage: str, scale: float = 1.0) -> list[PolConfig]:
    def k(n):
        return max(2, int(round(n * scale)))
    C: list[PolConfig] = []
    args = ("aware_cal", "aware_fit", "aware_over", "aware_under", "blind_flat", "blind_marg")
    if stage == "snp":            # SNP-vs-SNP negative control
        for arg in args:
            C.append(PolConfig(f"snpsnp_{arg}_M4000", "snp_snp", arg, 4000, 200, k(1000)))
        for arg in ("aware_fit",):
            for M in (1000, 20000):
                C.append(PolConfig(f"snpsnp_{arg}_M{M}", "snp_snp", arg, M, 200, k(1000)))
    elif stage == "te":           # TE-vs-SNP, all ARG models x delta x TE rule, M = 4000
        for arg in args:
            for delta in (0.0, 0.03, 0.09):
                for rule, filt in (("w1", True), ("q", True), ("bayes", True),
                                   ("w1", False), ("bayes", False)):
                    tag = f"{rule}{'' if filt else '_nofilt'}"
                    C.append(PolConfig(f"tesnp_{arg}_d{delta:g}_{tag}_M4000", "te_snp", arg, 4000,
                                       200, k(1000), delta=delta, te_rule=rule, te_filter=filt))
    elif stage == "teM":          # M sweep for the realistic arm
        for M in (1000, 20000):
            for rule, filt in (("w1", True), ("q", True), ("bayes", True), ("bayes", False)):
                tag = f"{rule}{'' if filt else '_nofilt'}"
                C.append(PolConfig(f"tesnp_aware_fit_d0.03_{tag}_M{M}", "te_snp", "aware_fit", M,
                                   200, k(1000), delta=0.03, te_rule=rule, te_filter=filt))
    elif stage == "contrast":     # between-category Z contrast, M 1000 vs 20000
        for rule, filt in (("w1", True), ("bayes", True), ("bayes", False)):
            tag = f"{rule}{'' if filt else '_nofilt'}"
            C.append(PolConfig(f"contrast_te_aware_fit_d0.03_{tag}_M1000vs20000", "contrast_te",
                               "aware_fit", 1000, 200, k(1000), delta=0.03, te_rule=rule,
                               te_filter=filt, M2=20000))
        C.append(PolConfig("contrast_snp_aware_fit_M1000vs20000", "contrast_snp", "aware_fit",
                           1000, 200, k(1000), M2=20000))
    elif stage == "R1000":        # production R for the key configurations
        C.append(PolConfig("snpsnp_aware_fit_M4000_R1000", "snp_snp", "aware_fit", 4000, 1000, k(1000)))
        C.append(PolConfig("snpsnp_aware_cal_M4000_R1000", "snp_snp", "aware_cal", 4000, 1000, k(1000)))
        for rule, filt in (("w1", True), ("bayes", False)):
            tag = f"{rule}{'' if filt else '_nofilt'}"
            C.append(PolConfig(f"tesnp_aware_cal_d0.03_{tag}_M4000_R1000", "te_snp", "aware_cal",
                               4000, 1000, k(1000), delta=0.03, te_rule=rule, te_filter=filt))
    elif stage == "impute":       # low-priority diagnostic: isolates the variance effect
        for arg in ("aware_cal", "aware_fit"):
            C.append(PolConfig(f"snpsnp_{arg}_impute_M4000", "snp_snp", arg, 4000, 200, k(1000),
                               snp_rule="impute"))
            for rule, filt in (("w1", False), ("bayes", False)):
                tag = f"{rule}{'' if filt else '_nofilt'}"
                C.append(PolConfig(f"tesnp_{arg}_d0.03_{tag}_impute_M4000", "te_snp", arg, 4000,
                                   200, k(1000), delta=0.03, te_rule=rule, te_filter=filt,
                                   snp_rule="impute"))
    else:
        raise ValueError(stage)
    return C


def pol_verify() -> dict:
    """Analytic check of the tower-property claim for each ARG model (no sampling)."""
    t = true_spectrum()
    out = {}
    for name, model in ARG_MODELS.items():
        sp = snp_dist(model, "mean")
        space = TypeSpace([sp])
        e = space.expected_spectrum(space.probs(sp))
        r = phi_sfs(e, t)
        d = arg_diagnostics(model, 0.03)
        out[name] = dict(model=model, W_snp_mixture_vs_truth=r.value,
                         daf_snp_mixture_minus_truth=r.mean_daf_difference,
                         snp_calibration_gap=d["snp_mean_abs_calibration_gap"],
                         snp_mean_q_by_k=d["snp_mean_q_by_k"], te_mean_q_by_k=d["te_mean_q_by_k"],
                         te_frac_unanimous=d["te_frac_unanimous"],
                         te_unanimous_correct=d["te_unanimous_correct"],
                         te_majority_correct=d["te_majority_correct"],
                         snp_frac_unanimous=d["snp_frac_unanimous"])
        print(f"{name:12s} W(SNP mixture, truth)={r.value:.2e} dDAF={r.mean_daf_difference:+.4f} "
              f"calib_gap={d['snp_mean_abs_calibration_gap']:.4f} TE unan={d['te_frac_unanimous']:.3f} "
              f"unan_ok={d['te_unanimous_correct']:.3f} maj_ok={d['te_majority_correct']:.3f}", flush=True)
    return out


def pol_report(out: Path) -> None:
    for f in sorted(out.glob("polsummary_*.json")):
        meta = json.loads(f.read_text())
        print(f"\n## {f.name} wall {meta['wall_seconds']:.0f}s cpu {meta['cpu_seconds_total']:.0f}s")
        for name, s in meta["summaries"].items():
            line = (f"{name}: n={s['n']} rej05={s['rej_05']:.3f} [{s['rej_05_ci'][0]:.3f},{s['rej_05_ci'][1]:.3f}]"
                    f" rej01={s['rej_01']:.3f} z={s['z_mean']:.2f}+-{s['z_mean_se']:.2f} sd={s['z_sd']:.2f}"
                    f" daf+={s['frac_daf_pos']:.2f} dafdiff={s['mean_daf_diff']:+.4f}"
                    f" EW={s['expected_W']:.2e} EdDAF={s['expected_mean_daf_diff']:+.4f}"
                    f" WA={s['W_A_truth']:.2e} WB={s['W_B_truth']:.2e}")
            for key in ("filter_retained", "filter_deleted_retained", "filter_derived_retained",
                        "z1_mean", "z2_mean"):
                if key in s:
                    line += f" {key}={s[key]:.3f}"
            print(line)


def pol_execute(cfgs, out: Path, workers: int, seed: int, tag: str, fit: dict) -> dict:
    out.mkdir(parents=True, exist_ok=True)
    tasks = []
    for cfg in cfgs:
        for j in range(math.ceil(cfg.analyses / cfg.chunk)):
            tasks.append((asdict(cfg), j, seed * 100_000 + zlib.crc32(cfg.name.encode()) % 100_000, fit))
    tasks.sort(key=lambda t: (t[0]["M"] + t[0]["M2"]) * t[0]["R"], reverse=True)
    results = {c.name: {} for c in cfgs}
    cpu = {c.name: 0.0 for c in cfgs}
    t0 = time.time()
    with Pool(workers) as pool:
        for name, j, recs, dt in pool.imap_unordered(run_pol_chunk, tasks):
            results[name][j] = recs
            cpu[name] += dt
    if not ARG_MODELS:
        _register_models(fit)
    summaries = {}
    for cfg in cfgs:
        recs = [r for j in sorted(results[cfg.name]) for r in results[cfg.name][j]]
        np.savez(out / f"{cfg.name}.npz", **{k: np.array([r[k] for r in recs]) for k in recs[0]})
        sm = pol_summarize(cfg, recs)
        sm.update(pol_expected(cfg))
        sm["cpu_seconds"] = cpu[cfg.name]
        summaries[cfg.name] = sm
    meta = dict(wall_seconds=time.time() - t0, cpu_seconds_total=sum(cpu.values()),
                workers=workers, seed=seed, fit=fit,
                models=ARG_MODELS, summaries=summaries)
    (out / f"polsummary_{tag}.json").write_text(json.dumps(meta, indent=1))
    print(f"[{tag}] wall {meta['wall_seconds']:.1f}s cpu {meta['cpu_seconds_total']:.1f}s", flush=True)
    return meta


# ---------------------------------------------------------------------------
# Real-data check: does the ARG's polarity posterior depend on allele frequency?
# ---------------------------------------------------------------------------

REALQ_VCF = ("/quobyte/jrigrp/beil/te_evo/singer_analysis/argtest/results/combined/"
             "nam_te_vcf/chr.10.combined.snp.te.sorted.vcf")
REALQ_STORE = ("/quobyte/jrigrp/beil/te_evo/singer_analysis/argtest/results/combined/"
               "all_te_snp_age_interval_store")


def real_q_check(out: Path, vcf: str = REALQ_VCF, store: str = REALQ_STORE,
                 ancestral: str = "results/ancestral-75draw") -> dict:
    """Relate the ARG polarity proportion to allele frequency on chromosome 10.

    Reads (read-only) the chromosome-10 SINGER input VCF, the interval-store
    positions and the 75-draw ancestral table. SNPs: q_ALT = draws calling REF
    ancestral / draws calling REF or ALT ancestral. TEs (REF=A absence,
    ALT=G insertion): q_ins = draws calling A ancestral / draws calling A or G.
    Genotypes are the haploid calls '0'/'1'; other calls are treated as missing,
    as in the earlier chr10 probes. Sites need n >= 20 and 0 < k < n.
    """
    bases = "ACGT"
    pos, k_alt, n_call, ref_i, alt_i, is_te = [], [], [], [], [], []
    b73_alt = 0
    with open(vcf) as h:
        for raw in h:
            if raw[0] == "#":
                continue
            f = raw.rstrip("\n").split("\t")
            if len(f[3]) != 1 or len(f[4]) != 1 or f[3] not in bases or f[4] not in bases:
                continue
            gts = f[9:]
            k = gts.count("1")
            m = k + gts.count("0")
            if m < PROJECTION_SIZE or k == 0 or k == m:
                continue
            te = f[2] != "."
            if not te and gts[0] == "1":
                b73_alt += 1
            pos.append(int(f[1])); k_alt.append(k); n_call.append(m)
            ref_i.append(bases.index(f[3])); alt_i.append(bases.index(f[4])); is_te.append(te)
    pos = np.array(pos, dtype=np.float64)
    k_alt = np.array(k_alt); n_call = np.array(n_call)
    ref_i = np.array(ref_i); alt_i = np.array(alt_i); is_te = np.array(is_te)
    meta = json.loads((Path(store) / "metadata.json").read_text())
    offset = {c["chrom"]: c["offset"] for c in meta["chromosomes"]}["10"]
    cat = np.load(Path(store) / "positions.npy", mmap_mode="r")
    g = pos + offset
    lo = int(np.searchsorted(cat, g.min())); hi = int(np.searchsorted(cat, g.max(), side="right"))
    sub = np.asarray(cat[lo:hi])
    ins = np.searchsorted(sub, g)
    ok = ins < sub.size
    ok[ok] &= sub[ins[ok]] == g[ok]
    rows = ins[ok] + lo
    A = np.load(Path(ancestral) / "ancestral_counts.npy", mmap_mode="r")
    cnt = np.asarray(A[rows], dtype=np.float64)
    k_alt, n_call, ref_i, alt_i, is_te = k_alt[ok], n_call[ok], ref_i[ok], alt_i[ok], is_te[ok]
    c_ref = cnt[np.arange(rows.size), ref_i]
    c_alt = cnt[np.arange(rows.size), alt_i]
    usable = (c_ref + c_alt) > 0
    q_alt = np.where(usable, c_ref / np.maximum(c_ref + c_alt, 1), np.nan)
    draws = c_ref + c_alt
    result: dict = {"vcf": vcf, "records_parsed_polymorphic_n20": int(ok.size),
                    "matched_store_rows": int(ok.sum()), "snp_b73_carries_alt": b73_alt}

    def table(mask, freq_count, qv, nv, label):
        rows_out = []
        # bin by projected-scale frequency: count/n in 13 bins of width 1/13 (~2 lines at n=26)
        edges = np.linspace(0, 1, 14)
        f = freq_count / nv
        for i in range(13):
            sel = mask & usable & (f > edges[i]) & (f <= edges[i + 1])
            if sel.sum() == 0:
                continue
            qq = qv[sel]
            rows_out.append(dict(
                bin=f"({edges[i]:.3f},{edges[i+1]:.3f}]", n=int(sel.sum()),
                mean_freq=float(f[sel].mean()), mean_q=float(qq.mean()),
                frac_q1=float(np.mean(qq == 1.0)), frac_q0=float(np.mean(qq == 0.0)),
                frac_q_ge_half=float(np.mean(qq > 0.5)),
                frac_unanimous=float(np.mean((qq == 1.0) | (qq == 0.0))),
                mean_draws=float(draws[sel].mean()),
                neutral_prior_alt_derived=float(np.mean(1 - f[sel]))))
        result[label] = rows_out

    snp = ~is_te
    table(snp, k_alt, q_alt, n_call, "snp_q_alt_by_alt_freq")
    # TE: insertion = ALT (G), q_ins = q_alt
    table(is_te, k_alt, q_alt, n_call, "te_q_ins_by_ins_freq")
    for label, m in (("snp", snp), ("te", is_te)):
        sel = m & usable
        f = k_alt[sel] / n_call[sel]
        qq = q_alt[sel]
        # Spearman correlation of q with frequency (ranks via argsort)
        def rank(x):
            r = np.empty(x.size); r[np.argsort(x, kind="mergesort")] = np.arange(x.size); return r
        rs = float(np.corrcoef(rank(f), rank(qq))[0, 1])
        result[f"{label}_summary"] = dict(n=int(sel.sum()), mean_q=float(qq.mean()),
                                          spearman_q_vs_freq=rs,
                                          frac_unanimous=float(np.mean((qq == 1) | (qq == 0))))
    te = is_te & usable
    qt = q_alt[te]
    unan = (qt == 1) | (qt == 0)
    result["te_summary"].update(
        majority_correct=float(np.mean(qt > 0.5)),
        unanimous_correct=float(np.mean(qt[unan] == 1)) if unan.any() else None,
        n_unanimous=int(unan.sum()),
        mean_per_draw_accuracy=float(qt.mean()))
    # SNP ALT-count spectrum at n = 26 for reference
    s26 = snp & usable & (n_call == 26)
    result["snp_alt_count_spectrum_n26"] = np.bincount(k_alt[s26], minlength=26).tolist()
    t26 = is_te & usable & (n_call == 26)
    result["te_ins_count_spectrum_n26"] = np.bincount(k_alt[t26], minlength=26).tolist()
    # Mean q by exact count at n = 26
    result["snp_mean_q_by_k_n26"] = [float(q_alt[s26 & (k_alt == k)].mean()) if np.any(s26 & (k_alt == k)) else None for k in range(26)]
    result["te_mean_q_by_k_n26"] = [float(q_alt[t26 & (k_alt == k)].mean()) if np.any(t26 & (k_alt == k)) else None for k in range(26)]
    out.mkdir(parents=True, exist_ok=True)
    (out / "real_q_check.json").write_text(json.dumps(result, indent=1))
    print(json.dumps({k: v for k, v in result.items() if "summary" in k or k.startswith("rec") or k.startswith("matched")}, indent=1))
    return result


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("mode", choices=["pilot", "run", "report", "check", "realq", "polfit",
                                     "polpilot", "polrun", "polreport", "polverify"])
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--parts", default="null,asym,polar,power")
    ap.add_argument("--workers", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 1)))
    ap.add_argument("--seed", type=int, default=20260925)
    ap.add_argument("--scale", type=float, default=1.0, help="multiply analyses per config")
    ap.add_argument("--only", default="", help="comma-separated config names to run (run mode)")
    ap.add_argument("--tag", default="", help="suffix for the summary file name")
    a = ap.parse_args(argv)
    if a.mode == "check":
        self_check()
        return 0
    if a.mode == "realq":
        real_q_check(a.out)
        return 0
    if a.mode in ("polpilot", "polrun", "polreport", "polverify"):
        fit = dict(POL_FIT)
        _register_models(fit)
        a.out.mkdir(parents=True, exist_ok=True)
        if a.mode == "polverify":
            (a.out / "polverify.json").write_text(json.dumps(pol_verify(), indent=1))
            return 0
        if a.mode == "polreport":
            pol_report(a.out)
            return 0
        for stage in a.parts.split(","):
            if a.mode == "polpilot":
                cfgs = pol_configs(stage, scale=0.0)
                for c in cfgs:
                    c.chunk = 2
                meta = pol_execute(cfgs, a.out / "pilot", a.workers, a.seed + 7, f"pilot_{stage}", fit)
                est = sum(meta["summaries"][c.name]["cpu_seconds"] / 2 * c.analyses
                          for c in pol_configs(stage, 1.0))
                print(f"[pilot {stage}] estimated full cpu {est/3600:.2f} h", flush=True)
            else:
                cfgs = pol_configs(stage, a.scale)
                if a.only:
                    keep = set(a.only.split(","))
                    cfgs = [c for c in cfgs if c.name in keep]
                pol_execute(cfgs, a.out, a.workers, a.seed, stage + a.tag, fit)
        return 0
    if a.mode == "polfit":
        fits = {}
        for prior in ("aware", "flat", "marginal"):
            for delta in (0.03, 0.09):
                f = fit_arg(prior, delta)
                fits[f"{prior}_delta{delta:g}"] = f
                print(prior, delta, f["model"], "sse", round(f["sse"], 5), "rmse", round(f["rmse_te_q"], 4),
                      "unan", round(f["diag"]["te_frac_unanimous"], 3),
                      "unan_ok", round(f["diag"]["te_unanimous_correct"], 3),
                      "maj_ok", round(f["diag"]["te_majority_correct"], 3), flush=True)
        a.out.mkdir(parents=True, exist_ok=True)
        (a.out / "polfit.json").write_text(json.dumps(fits, indent=1))
        return 0
    if a.mode == "report":
        report(a.out)
        return 0
    self_check()
    for part in a.parts.split(","):
        if a.mode == "pilot":
            cfgs = configs_for(part, scale=0.0)  # 2 analyses per config
            for c in cfgs:
                c.chunk = 2
            meta = execute(cfgs, a.out / "pilot", a.workers, a.seed + 7, f"pilot_{part}")
            full = configs_for(part, 1.0)
            est = 0.0
            for c in full:
                s = meta["summaries"][c.name]
                est += s["cpu_seconds"] / s["n"] * c.analyses
            print(f"[pilot {part}] estimated full cpu {est/3600:.2f} h", flush=True)
        else:
            cfgs = configs_for(part, a.scale)
            if a.only:
                keep = set(a.only.split(","))
                cfgs = [c for c in cfgs if c.name in keep]
                if not cfgs:
                    continue
            execute(cfgs, a.out, a.workers, a.seed, part + a.tag)
    return 0


if __name__ == "__main__":
    sys.exit(main())
