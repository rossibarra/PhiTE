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


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("mode", choices=["pilot", "run", "report", "check"])
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
