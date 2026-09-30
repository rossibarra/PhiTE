#!/usr/bin/env python3
"""Calculate Phi-SFS for a focal A set and its matched SNP control sets.

Phi-SFS is the first Wasserstein distance between the projected, normalized,
unfolded site frequency spectrum of a focal set and that of one age-matched
SNP control set. README section 8 and PHI_SFS_WASSERSTEIN_CODING_PLAN.md carry
the full derivation.

Input assumptions, all of which are recorded in the output metadata:

* The VCF is biallelic. Multiallelic records are produced upstream only as
  separate biallelic records, so a comma in ALT is treated as an error rather
  than split here.
* The VCF FILTER column is ignored. The declared input is the already
  filtered preprocessing VCF, so every record at a requested coordinate is used.
* The VCF is **not** assumed to be polarized, and no REF or INFO annotation is
  consulted. Every site, TE and SNP alike, is polarized by the ARG: its weight
  is the posterior proportion of usable draws in which ALT is derived, and its
  contribution is the q-mixture of its two orientations. TE records are A/G
  biallelic sites in the ARG like any other, so the same ancestral table orients
  them, and a TE whose absence is derived is counted at the absence frequency.
  By default A, B0 and every null-left set are all posterior mixtures, so the
  observed and null comparisons are built identically.
  `--asymmetric-polarity-null` instead hard-orients A and every null-left set
  by one reproducible Bernoulli(q) draw per site, keeping B0 a mixture.
"""

from __future__ import annotations

import argparse
import csv
import datetime
import hashlib
import json
import math
import os
import shutil
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import NamedTuple, Sequence

import numpy as np

from .build_ancestral_states import SCHEMA_VERSION as ANCESTRAL_SCHEMA_VERSION
from .build_ancestral_states import verify_table_arrays
from .release_provenance import software_provenance
from .sample_age_matched_controls import _load_target, _sha256_arrays
from .vcf_io import (  # noqa: F401 -- re-exported for existing callers
    COMPRESSED_SUFFIXES,
    _HashingStream,
    _decode_genotype,
    _open_vcf,
    drain,
)


SCHEMA_VERSION = "phi-sfs-wasserstein-v3"
ASYMMETRIC_NULL_DESIGN = (
    "bernoulli-q-hard-vs-posterior-mixture"
)
SYMMETRIC_NULL_DESIGN = (
    "posterior-mixture-vs-posterior-mixture"
)

# The only supported matched-control schema, and the per-replicate identifier
# array it publishes. The bootstrap-target matcher produces replicates with no
# chain structure, so it identifies them by replicate alone; inventing chain
# and sample columns would imply a within-chain correlation that does not
# exist. "No chain structure" is the whole claim -- these replicates are not
# statistically independent, since they share the observed TE sample and the
# interval store. An earlier swap sampler identified replicates by chain and
# sample instead; `calculate()` now accepts only this schema, so that mapping
# is gone rather than carried as dead code.
SUPPORTED_MATCH_SCHEMA = "bootstrap-target-matches-v1"
MATCH_IDENTIFIER_NAMES = ("replicate_id",)

PROJECTION_SIZE = 20
RETAINED_BINS = np.arange(1, PROJECTION_SIZE, dtype=np.int64)
PROGRESS_RECORDS = 5_000_000

_UNSET = object()
_GENOTYPE_ERRORS = {
    "invalid": "invalid GT {gt!r} at {chrom}:{position}",
    "non-biallelic": "non-biallelic GT {gt!r} at {chrom}:{position}",
    "heterozygous": "heterozygous GT {gt!r} at inbred site {chrom}:{position}",
}


@dataclass(frozen=True)
class SiteCount:
    """One site's observed counts, with polarity carried rather than applied.

    `alt` and `callable` come from the genotypes and are polarity-independent.
    `p_alt_derived` is the probability that ALT is the derived allele: the ARG's
    posterior proportion at every site, TE or SNP, in A or in B. Orientation is
    applied when spectra are summed, not here, because that is the only step
    that depends on it.
    """

    alt: int
    callable: int
    p_alt_derived: float


class PhiResult(NamedTuple):
    """Wasserstein Phi-SFS and the arrays that determine its direction."""

    value: float
    cdf_a: np.ndarray
    cdf_b: np.ndarray
    cdf_residual: np.ndarray
    bin_residual: np.ndarray
    mean_daf_difference: float


class PhiCalibration(NamedTuple):
    """An observed Phi-SFS calibrated against neutral-null distances."""

    observed: float
    null: np.ndarray
    null_mean: float
    null_sd: float
    z_score: float
    p_value: float
    exceedances: int
    null_z_scores: np.ndarray


def hypergeometric_projection(k: int, n: int, m: int = PROJECTION_SIZE) -> np.ndarray:
    """Return the expected derived-count distribution after projection to m.

    A site observed with `k` derived alleles among `n` callable inbred
    individuals contributes probability mass to projected bin `j` equal to

        h_j(k, n) = C(k, j) * C(n - k, m - j) / C(n, m),   j = 0, ..., m,

    the hypergeometric probability of drawing `j` derived alleles in a sample
    of `m` drawn without replacement. This is the exact expectation over all
    subsamples, not a random downsampling draw, so the result is deterministic.

    The returned vector covers bins 0 through m inclusive and sums to one.
    Sites with `n < m` cannot be projected and are rejected here; callers drop
    them before reaching this function. Probabilities are evaluated in log
    space via `lgamma` so that large `n` stays finite: the identity holds to
    full double precision at n = 2e6.
    """
    if not isinstance(k, (int, np.integer)) or not isinstance(n, (int, np.integer)):
        raise TypeError("k and n must be integers")
    k, n = int(k), int(n)
    if n < m:
        raise ValueError(f"cannot project n={n} observations to m={m}")
    if k < 0 or k > n:
        raise ValueError(f"derived count k={k} must satisfy 0 <= k <= n={n}")
    result = np.zeros(m + 1, dtype=np.float64)
    lower = max(0, m - (n - k))
    upper = min(m, k)
    denominator = math.lgamma(n + 1) - math.lgamma(m + 1) - math.lgamma(n - m + 1)
    for j in range(lower, upper + 1):
        log_probability = (
            math.lgamma(k + 1) - math.lgamma(j + 1) - math.lgamma(k - j + 1)
            + math.lgamma(n - k + 1)
            - math.lgamma(m - j + 1)
            - math.lgamma(n - k - (m - j) + 1)
            - denominator
        )
        result[j] = math.exp(log_probability)
    total = float(result.sum())
    if not math.isfinite(total) or total <= 0:
        raise RuntimeError(f"invalid hypergeometric projection for k={k}, n={n}, m={m}")
    result /= total
    return result


def project_sites(
    counts: dict[tuple[str, int], SiteCount],
) -> tuple[dict[tuple[str, int], int], np.ndarray, np.ndarray]:
    """Project every eligible site, evaluating each distinct (k, n) pair once.

    Sites with fewer than `PROJECTION_SIZE` callable individuals are ineligible
    and are absent from the returned coordinate map; that is the `n >= 20`
    filter, and it is expected behaviour rather than an error.

    Each returned row holds bins 1 through 19 of one site's projection and is
    deliberately **not** renormalized after the endpoint bins are removed, so
    it sums to `1 - h_0 - h_20` rather than to one. A site whose derived count
    is very likely to project to 0 or 20 therefore contributes proportionally
    less polymorphic mass, which is the intended weighting. The companion
    endpoint array holds the excluded `h_0 + h_20` mass for diagnostics.

    Each projection is a polarity-weighted mixture of the site's two possible
    orientations (see below). Because reversing a projection swaps `h_0` with
    `h_20`, the retained mass `1 - h_0 - h_20` is invariant under polarity, so
    the weighting redistributes mass within a site and never reweights sites
    against each other.

    Returns the coordinate-to-row map, an (n_distinct, 19) projection matrix,
    and the aligned (n_distinct,) endpoint masses.
    """
    rows: dict[tuple[str, int], int] = {}
    distinct: dict[SiteCount, int] = {}
    retained: list[np.ndarray] = []
    endpoints: list[float] = []
    for coordinate, item in counts.items():
        if item.callable < PROJECTION_SIZE:
            continue
        row = distinct.get(item)
        if row is None:
            # Polarity is a probability, not a decision, so a site contributes
            # the posterior mean of its two orientations:
            #     p * h(k, n) + (1 - p) * h(n - k, n)
            # and h(n - k, n) is exactly h(k, n) reversed, so one projection
            # serves both. The mixture is linear in p and therefore unbiased at
            # any number of posterior draws, and a site the ARG cannot polarize
            # (p near 0.5) contributes a near-symmetric, effectively folded
            # shape rather than a confidently wrong one.
            vector = hypergeometric_projection(item.alt, item.callable)
            p = float(item.p_alt_derived)
            if not 0.0 <= p <= 1.0:
                raise ValueError(f"polarity weight out of range at {coordinate}: {p}")
            mixed = p * vector + (1.0 - p) * vector[::-1]
            row = len(retained)
            distinct[item] = row
            retained.append(mixed[1:PROJECTION_SIZE])
            endpoints.append(float(mixed[0] + mixed[PROJECTION_SIZE]))
        rows[coordinate] = row
    projections = (
        np.asarray(retained, dtype=np.float64) if retained
        else np.zeros((0, PROJECTION_SIZE - 1), dtype=np.float64)
    )
    return rows, projections, np.asarray(endpoints, dtype=np.float64)


def _site_uniform(seed: int, coordinate: tuple[str, int]) -> float:
    """Return a stable coordinate-keyed U[0,1) variate."""
    chrom, position = coordinate
    payload = f"phi-sfs-bernoulli-q-v1\0{seed}\0{chrom}\0{position}".encode()
    value = int.from_bytes(hashlib.sha256(payload).digest()[:8], "big")
    return value / 2**64


def project_sites_bernoulli_q(
    counts: dict[tuple[str, int], SiteCount],
    *,
    seed: int,
) -> tuple[dict[tuple[str, int], int], np.ndarray, np.ndarray]:
    """Project sites after one reproducible hard orientation drawn from q.

    For q = P(ALT derived | ARG), ALT is declared derived when a
    coordinate-keyed U[0,1) draw is below q; otherwise REF is declared derived.
    The expected hard projection is therefore the posterior-mixture projection.
    """
    rows: dict[tuple[str, int], int] = {}
    distinct: dict[tuple[int, int, bool], int] = {}
    retained: list[np.ndarray] = []
    endpoints: list[float] = []
    for coordinate, item in counts.items():
        if item.callable < PROJECTION_SIZE:
            continue
        q = float(item.p_alt_derived)
        if not 0.0 <= q <= 1.0:
            raise ValueError(f"polarity weight out of range at {coordinate}: {q}")
        alt_derived = _site_uniform(seed, coordinate) < q
        key = (item.alt, item.callable, alt_derived)
        row = distinct.get(key)
        if row is None:
            vector = hypergeometric_projection(item.alt, item.callable)
            if not alt_derived:
                vector = vector[::-1]
            row = len(retained)
            distinct[key] = row
            retained.append(vector[1:PROJECTION_SIZE])
            endpoints.append(float(vector[0] + vector[PROJECTION_SIZE]))
        rows[coordinate] = row
    projections = (
        np.asarray(retained, dtype=np.float64) if retained
        else np.zeros((0, PROJECTION_SIZE - 1), dtype=np.float64)
    )
    return rows, projections, np.asarray(endpoints, dtype=np.float64)


def accumulate_spectrum(
    coordinates: Sequence[tuple[str, int]],
    rows: dict[tuple[str, int], int],
    projections: np.ndarray,
    endpoints: np.ndarray,
) -> tuple[np.ndarray, float, int]:
    """Sum the projections of the eligible members of one site set.

    Sites are gathered by distinct (k, n) row and weighted by how often that
    row occurs, so a set of N sites costs one length-19 matrix product rather
    than N array additions. Returns the unnormalized bins 1 through 19, the
    excluded endpoint mass, and the number of eligible sites.
    """
    indices = np.fromiter(
        (rows[coordinate] for coordinate in coordinates if coordinate in rows),
        dtype=np.int64,
    )
    weights = np.bincount(indices, minlength=projections.shape[0]).astype(np.float64)
    return weights @ projections, float(weights @ endpoints), int(indices.size)


def normalized_spectrum(raw: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Normalize one accumulated spectrum over bins 1 through 19.

    Normalization happens once per set, after every site has contributed, so a
    site's weight in the final spectrum is proportional to its retained
    projection mass. Normalizing site by site instead would give every eligible
    site equal weight and would discard the endpoint-mass information.

    Normalization also discards the set's absolute scale, so differences in
    eligible-site count and total retained mass become invisible in the final
    spectrum. Callers must report the retained fractions separately.
    """
    values = np.asarray(raw, dtype=np.float64)
    if values.shape != (PROJECTION_SIZE - 1,):
        raise ValueError("spectrum must contain bins 1 through 19")
    total = float(values.sum())
    if not math.isfinite(total) or total <= 0:
        raise ValueError("spectrum has zero retained mass in bins 1 through 19")
    return values, values / total


def phi_sfs(
    a: np.ndarray,
    b: np.ndarray,
    *,
    daf: np.ndarray | None = None,
) -> PhiResult:
    """Return the one-dimensional Wasserstein distance between two SFS arrays.

    The inputs are probability masses on the same strictly increasing derived
    allele-frequency grid. On a discrete one-dimensional grid, Wasserstein-1
    is the area between their CDFs. The default grid is the retained projected
    counts 1 through 19 divided by the projection size of 20.

    The distance is unsigned. ``mean_daf_difference`` and the oriented
    residual arrays use ``a - b`` and retain information about direction.
    """
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    if a.ndim != 1 or b.ndim != 1 or a.shape != b.shape:
        raise ValueError("normalized spectra must be one-dimensional and equal length")
    if a.size < 2:
        raise ValueError("normalized spectra must contain at least two bins")
    if not np.all(np.isfinite(a)) or not np.all(np.isfinite(b)):
        raise ValueError("normalized spectra must be finite")
    if np.any(a < 0.0) or np.any(b < 0.0):
        raise ValueError("normalized spectra must be nonnegative")
    if not np.isclose(a.sum(), 1.0, rtol=0.0, atol=1e-8) or not np.isclose(
        b.sum(), 1.0, rtol=0.0, atol=1e-8
    ):
        raise ValueError("both spectra must be normalized")

    if daf is None:
        if a.size != RETAINED_BINS.size:
            raise ValueError(
                "daf is required unless spectra use the default bins 1 through 19"
            )
        grid = RETAINED_BINS.astype(np.float64) / PROJECTION_SIZE
    else:
        grid = np.asarray(daf, dtype=np.float64)
        if grid.ndim != 1 or grid.shape != a.shape:
            raise ValueError("daf must be one-dimensional and match the spectra")
        if not np.all(np.isfinite(grid)):
            raise ValueError("daf must be finite")
        if np.any(np.diff(grid) <= 0.0):
            raise ValueError("daf must be strictly increasing")

    cdf_a = np.cumsum(a)
    cdf_b = np.cumsum(b)
    cdf_residual = cdf_a - cdf_b
    bin_residual = a - b
    value = float(np.sum(np.abs(cdf_residual[:-1]) * np.diff(grid)))
    mean_daf_difference = float(np.dot(bin_residual, grid))
    return PhiResult(
        value,
        cdf_a,
        cdf_b,
        cdf_residual,
        bin_residual,
        mean_daf_difference,
    )


def calibrate_phi(observed: float, null: np.ndarray) -> PhiCalibration:
    """Standardize an observed distance and compute its Monte Carlo P-value.

    The Z-score uses the null sample standard deviation (``ddof=1``). The
    one-sided P-value counts null distances greater than or equal to the
    observed distance and applies the add-one correction to numerator and
    denominator.
    """
    observed_array = np.asarray(observed, dtype=np.float64)
    if observed_array.ndim != 0:
        raise ValueError("observed distance must be a scalar")
    observed_value = float(observed_array)
    if not math.isfinite(observed_value) or observed_value < 0.0:
        raise ValueError("observed distance must be finite and nonnegative")

    null_values = np.asarray(null, dtype=np.float64)
    if null_values.ndim != 1:
        raise ValueError("null distances must be one-dimensional")
    if null_values.size < 2:
        raise ValueError("at least two null distances are required")
    if not np.all(np.isfinite(null_values)) or np.any(null_values < 0.0):
        raise ValueError("null distances must be finite and nonnegative")

    null_mean = float(np.mean(null_values))
    null_sd = float(np.std(null_values, ddof=1))
    if null_sd == 0.0:
        raise ValueError("null distances have zero sample standard deviation")
    z_score = (observed_value - null_mean) / null_sd
    exceedances = int(np.count_nonzero(null_values >= observed_value))
    p_value = (1.0 + exceedances) / (null_values.size + 1.0)
    null_z_scores = (null_values - null_mean) / null_sd
    return PhiCalibration(
        observed_value,
        null_values.copy(),
        null_mean,
        null_sd,
        float(z_score),
        float(p_value),
        exceedances,
        null_z_scores,
    )


def _checked_table_array(path: Path, shape: tuple[int, ...]) -> np.ndarray:
    """Load an ancestral-table array, refusing a wrong shape or a signed dtype.

    Both arrays are counts indexed by store row, so a shape mismatch means the
    table does not describe this store and a signed or floating dtype means it
    was not written by the builder. Either would otherwise surface as silently
    misaligned polarity rather than an error.
    """
    array = np.load(path, mmap_mode="r", allow_pickle=False)
    if array.shape != shape:
        raise ValueError(f"{path.name} has shape {array.shape}, expected {shape}")
    if array.dtype.kind != "u":
        raise ValueError(f"{path.name} has dtype {array.dtype}, expected unsigned")
    return array


class PolarityResolver:
    """Resolves P(ALT is derived) per site from the ARG ancestral table.

    TE and SNP sites are resolved the same way. A TE insertion is biologically
    derived, but taking that as known while the SNP controls carry the ARG's
    uncertain polarity makes the observed and null comparisons differ in
    construction, and in simulation that asymmetry alone produced large excess
    rejection. So a TE is oriented by the ARG exactly as a SNP is; TE records
    are A/G sites (A = absence) and a draw calling G ancestral makes absence the
    derived allele. `te_coordinates` is used only to count the two kinds apart.

    The weight is the posterior proportion of ARG draws calling REF ancestral --
    equivalently, ALT derived -- among the draws that named one of the two
    observed alleles. Conditioning that way rather than on the raw present-draw
    count is what lets every requested site carry a weight without an
    intersection across draws or a fallback rule, and it discards draws naming a
    third base, which cannot orient the site.

    The proportion is used as reported. Against TE ground truth the ARG is only
    about 91% correct where all its draws agree, so this weight is somewhat
    overconfident and SNP spectra come out sharper than the ARG's measured
    accuracy warrants. That is a deliberate, recorded choice, not an oversight.

    Resolution happens during the VCF scan because the weight depends on which
    allele is ALT, which only the scan knows.
    """

    BASES = "ACGT"

    def __init__(
        self,
        te_coordinates: set[tuple[str, int]],
        *,
        store_positions: np.ndarray,
        chromosome_offsets: dict[str, int],
        ancestral_counts: np.ndarray,
        present_draw_count: np.ndarray,
    ) -> None:
        self._te = te_coordinates
        self._positions = np.asarray(store_positions)
        self._offsets = chromosome_offsets
        self._counts = ancestral_counts
        self._present = present_draw_count
        self.te_sites = 0
        self.te_usable_draws = 0
        self.te_unusable_draws = 0
        self.control_sites = 0
        self.control_usable_draws = 0
        self.control_unusable_draws = 0

    def __call__(self, chrom: str, position: int, ref: str, alt: str) -> float:
        offset = self._offsets.get(chrom)
        if offset is None:
            raise ValueError(f"no chromosome offset for {chrom!r}")
        target = float(offset + position)
        index = int(np.searchsorted(self._positions, target))
        if index >= self._positions.size or self._positions[index] != target:
            raise ValueError(
                f"site {chrom}:{position} is absent from the ancestral "
                "table; it cannot be polarized"
            )
        if ref not in self.BASES or alt not in self.BASES:
            raise ValueError(f"non-ACGT allele at {chrom}:{position}: {ref}/{alt}")
        # The table counts how many draws called each base ANCESTRAL, so ALT is
        # derived exactly when REF is the ancestral call. Condition on draws that
        # picked one of the two observed alleles: a draw naming a third base
        # cannot orient this site and is not evidence either way.
        row = self._counts[index]
        if row.sum() > self._present[index]:
            raise ValueError(
                f"ancestral table is inconsistent at {chrom}:{position}: "
                f"{int(row.sum())} ancestral calls across {int(self._present[index])} "
                "draws"
            )
        ref_calls = float(row[self.BASES.index(ref)])
        alt_calls = float(row[self.BASES.index(alt)])
        oriented = ref_calls + alt_calls
        if oriented <= 0:
            raise ValueError(
                f"no posterior draw at {chrom}:{position} calls either observed "
                f"allele ({ref}/{alt}) ancestral; it cannot be polarized"
            )
        unusable = int(self._present[index] - oriented)
        if (chrom, position) in self._te:
            self.te_sites += 1
            self.te_usable_draws += int(oriented)
            self.te_unusable_draws += unusable
        else:
            self.control_sites += 1
            self.control_usable_draws += int(oriented)
            self.control_unusable_draws += unusable
        return ref_calls / oriented


def read_site_counts(
    vcf: Path,
    coordinates: set[tuple[str, int]],
    *,
    polarity: PolarityResolver,
    heterozygous: str,
    progress: bool = True,
) -> tuple[dict[tuple[str, int], SiteCount], str]:
    """Read callable and polarized derived counts for the requested coordinates.

    Returns the per-coordinate counts and the SHA-256 of the VCF bytes, which
    is accumulated during this single pass rather than in a second full read.

    Only CHROM and POS are parsed for records that are not requested, because
    splitting every sample column of every record dominates the scan on a
    sample-rich VCF.
    """
    found: dict[tuple[str, int], SiteCount] = {}
    genotypes: dict[str, int | None | str] = {}
    handle, hashing, buffered = _open_vcf(vcf)
    next_report = PROGRESS_RECORDS
    records = 0
    try:
        for line_number, raw in enumerate(handle, 1):
            if raw.startswith("#"):
                continue
            records += 1
            if progress and records >= next_report:
                print(f"  scanned {records:,} VCF records", flush=True)
                next_report += PROGRESS_RECORDS
            try:
                chrom, position_text, rest = raw.split("\t", 2)
            except ValueError:
                raise ValueError(f"{vcf}:{line_number}: malformed VCF record") from None
            try:
                position = int(position_text)
            except ValueError as error:
                raise ValueError(f"{vcf}:{line_number}: invalid POS") from error
            coordinate = (chrom, position)
            if coordinate not in coordinates:
                continue
            if coordinate in found:
                raise ValueError(f"duplicate VCF record at {chrom}:{position}")
            fields = rest.rstrip("\n").split("\t")
            if len(fields) < 8:
                raise ValueError(f"{vcf}:{line_number}: expected VCF samples and FORMAT")
            ref, alt, formats = fields[1], fields[2], fields[6]
            if "," in alt:
                raise ValueError(f"multiallelic VCF record at {chrom}:{position}")
            weight = polarity(chrom, position, ref, alt)
            format_fields = formats.split(":")
            if "GT" not in format_fields:
                raise ValueError(f"VCF record lacks GT at {chrom}:{position}")
            gt_index = format_fields.index("GT")
            alt_count = 0
            callable_count = 0
            for sample in fields[7:]:
                if gt_index:
                    parts = sample.split(":")
                    gt = parts[gt_index] if gt_index < len(parts) else "."
                else:
                    end = sample.find(":")
                    gt = sample if end < 0 else sample[:end]
                allele = genotypes.get(gt, _UNSET)
                if allele is _UNSET:
                    allele = _decode_genotype(gt, heterozygous)
                    genotypes[gt] = allele
                if allele is None:
                    continue
                if allele.__class__ is str:
                    raise ValueError(_GENOTYPE_ERRORS[allele].format(
                        gt=gt, chrom=chrom, position=position,
                    ))
                alt_count += allele
                callable_count += 1
            found[coordinate] = SiteCount(
                alt=alt_count, callable=callable_count, p_alt_derived=weight,
            )
        drain(buffered)
    finally:
        handle.close()
    if progress:
        print(f"  scanned {records:,} VCF records", flush=True)
    return found, hashing.digest.hexdigest()


def _load_integers(path: Path, label: str) -> np.ndarray:
    """Load an integer array, refusing to coerce a non-integer dtype.

    Casting first would silently truncate: a position or row index stored as
    2.9 would become 2 and resolve to the wrong site rather than failing.
    """
    values = np.load(path, allow_pickle=False)
    if not np.issubdtype(values.dtype, np.integer):
        raise ValueError(f"{label} must be an integer array, not {values.dtype}")
    return values.astype(np.int64)


def _load_coordinates(target: Path, matches: Path, schema: str):
    """Load and cross-validate the target and matched-control site arrays.

    Returns the target's chromosomes, positions, and row indices; the matched
    sets' chromosomes, positions, and row indices; and the per-replicate
    identifier arrays for `schema` (a name-to-array mapping, carried through
    to the outputs unchanged). Row indices are returned here, rather than
    reloaded by the caller, because this function already loads and validates
    them.
    """
    if schema != SUPPORTED_MATCH_SCHEMA:
        raise ValueError(
            f"unsupported matched-control schema_version {schema!r}; "
            f"expected {SUPPORTED_MATCH_SCHEMA!r}"
        )
    te_chromosomes = np.load(target / "te_chromosomes.npy", allow_pickle=False).astype(str)
    labels = np.load(matches / "chromosome_labels.npy", allow_pickle=False).astype(str)
    te_positions = _load_integers(target / "te_positions.npy", "target positions")
    te_rows = _load_integers(target / "te_row_indices.npy", "target row indices")
    positions = _load_integers(matches / "positions.npy", "matched positions")
    codes = _load_integers(matches / "chromosome_codes.npy", "matched chromosome codes")
    rows = _load_integers(matches / "row_indices.npy", "matched row indices")
    identifiers = {
        name: _load_integers(matches / f"{name}.npy", f"{name} array")
        for name in MATCH_IDENTIFIER_NAMES
    }
    if te_chromosomes.shape != te_positions.shape or te_chromosomes.ndim != 1:
        raise ValueError("target chromosome and position arrays are not aligned 1-D arrays")
    if te_rows.shape != te_positions.shape:
        raise ValueError("target row indices do not align with target positions")
    if positions.shape != codes.shape or positions.ndim != 2:
        raise ValueError("matched chromosome codes and positions are not aligned 2-D arrays")
    if rows.shape != positions.shape:
        raise ValueError("matched row indices do not align with matched positions")
    if np.any(te_rows < 0) or np.any(rows < 0):
        raise ValueError("row indices must be non-negative")
    for name, values in identifiers.items():
        if values.shape != (positions.shape[0],):
            raise ValueError(f"{name} array does not align with matched sets")
    if np.any(codes < 0) or np.any(codes >= labels.size):
        raise ValueError("matched chromosome code is out of range")
    ordered = np.sort(rows, axis=1)
    if ordered.shape[1] > 1 and np.any(np.diff(ordered, axis=1) == 0):
        raise ValueError("a matched control set contains duplicate control rows")
    match_chromosomes = labels[codes]
    return (
        te_chromosomes, te_positions, te_rows,
        match_chromosomes, positions, rows,
        identifiers,
    )


def _json(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def _validate_provenance(target: Path, matches: Path) -> tuple[dict, dict, str]:
    """Require that the matched bundle was built from this exact target.

    The store hashes alone cannot establish this, because every target built
    from one SNP store shares them. The matcher records `target_digest`, a hash
    over the target's row indices, mean CDF, age grid, and acceptance
    threshold, so recomputing it from the target directory is what actually
    binds the two bundles together. The digest is computed with the matcher's
    own loader and hash helper so the two cannot drift apart.

    The store hashes must be *present* in both bundles and must agree, which
    rejects a hand-built or truncated bundle that carries no store identity at
    all. They are not required to be non-null: the dense store records neither
    digest, and every other step in this pipeline compares them only when a
    value exists, so demanding a non-null digest here would make this the one
    step that rejects a bundle the rest of the pipeline produced and accepted.
    """
    target_rows, target_cdf, age_bins, threshold, target_meta = _load_target(target)
    match_meta = _json(matches / "metadata.json")
    if match_meta.get("complete") is not True:
        raise ValueError(f"matched-control bundle is not marked complete: {matches}")
    for key in ("source_store_content_sha256", "source_catalog_sha256"):
        if key not in target_meta or key not in match_meta:
            raise ValueError(
                f"target and matched-control metadata must both record {key}"
            )
        if target_meta[key] != match_meta[key]:
            raise ValueError(f"target and matched-control {key} values differ")
    expected = match_meta.get("target_digest")
    if not expected:
        raise ValueError(f"matched-control metadata records no target_digest: {matches}")
    digest = _sha256_arrays(
        target_rows, target_cdf, age_bins, np.asarray([threshold], dtype=np.float64)
    )
    if digest != expected:
        raise ValueError(
            "matched-control bundle was built for a different target: "
            f"it records target_digest {expected}, but {target} hashes to {digest}"
        )
    return target_meta, match_meta, digest


def _validate_target_authority(target_meta: dict, match_meta: dict, a_type: str) -> None:
    """Require the target, not `-A`, to be the authority on A's variant type.

    Mirrors `bootstrap_target_matcher._validate_inputs` (wording only; that
    helper is private to the matcher and is not imported here). Without this,
    `--a-type` is trusted blindly and the output metadata mislabels A's variant
    type, which `phi_contrast` then groups on. These checks must pass before
    `a_polarity_rule` is written to the output metadata or the VCF is scanned.

    A target built with the TE polarity mask is refused for either type: its
    at-least-50%-derived filter and agreeing-draw ages condition A on the ARG's
    polarity in a way no SNP control set is, which is the asymmetry the
    posterior-polarity design removes.
    """
    target_a_type = target_meta.get("a_type")
    if target_a_type not in ("TE", "SNP"):
        raise ValueError(
            "target metadata does not declare a valid a_type; rebuild the "
            "target with normalize_tes.te_age_target --a-type TE or SNP"
        )
    if target_a_type != a_type:
        raise ValueError(
            f"--a-type {a_type} disagrees with target metadata a_type={target_a_type}"
        )
    match_a_type = match_meta.get("a_type")
    if match_a_type is not None and match_a_type != a_type:
        raise ValueError(
            f"--a-type {a_type} disagrees with matched-control metadata "
            f"a_type={match_a_type}"
        )
    if target_meta.get("te_polarity") is not None:
        raise ValueError(
            "target was built with the TE polarity mask; A is polarized by the "
            "ARG posterior like its SNP controls, so rebuild the target and "
            "matches without --te-polarity-mask or --max-flipped-fraction"
        )
    eligibility = target_meta.get("vcf_eligibility")
    if not isinstance(eligibility, dict) or not eligibility.get("mask"):
        raise ValueError(
            "age matching requires a final target built with --vcf-eligibility "
            "so A was filtered before its size and age CDF were fixed"
        )


def _validate_vcf_identity(target_meta: dict, match_meta: dict, *, heterozygous: str) -> dict:
    """Require the target and matched-control eligibility identities to agree.

    `target_meta["vcf_eligibility"]["identity"]` and
    `match_meta["vcf_eligibility_identity"]` are the same dict by contract --
    the matcher copies it from the candidate report the eligibility mask was
    built against -- so any difference means the two bundles rest on different
    VCF eligibility decisions and must not be combined. The caller still owes
    a check that the VCF actually scanned hashes to `identity["vcf_sha256"]`;
    that can only happen after the scan, so this returns the identity rather
    than checking the digest itself.
    """
    eligibility = target_meta.get("vcf_eligibility")
    identity = eligibility.get("identity") if isinstance(eligibility, dict) else None
    if not isinstance(identity, dict):
        raise ValueError(
            "target metadata vcf_eligibility record carries no identity; "
            "rebuild the target from an eligibility mask that records one"
        )
    match_identity = match_meta.get("vcf_eligibility_identity")
    if match_identity != identity:
        raise ValueError(
            "target and matched-control vcf_eligibility_identity values differ; "
            "the two bundles were built against different VCF eligibility masks"
        )
    if identity.get("heterozygous") != heterozygous:
        raise ValueError(
            f"--heterozygous {heterozygous} disagrees with the eligibility "
            f"mask's policy heterozygous={identity.get('heterozygous')!r}"
        )
    if identity.get("min_callable") != PROJECTION_SIZE:
        raise ValueError(
            f"eligibility mask min_callable={identity.get('min_callable')!r} "
            f"must equal the projection size {PROJECTION_SIZE}"
        )
    return identity


def calculate(args: argparse.Namespace) -> None:
    if args.min_null_replicates < 2:
        raise ValueError("--min-null-replicates must be at least 2 for a Z-score")
    if args.reference_sensitivity < 0:
        raise ValueError("--reference-sensitivity must be nonnegative")
    capped = args.max_null_replicates is not None
    if capped:
        if args.max_null_replicates < 2:
            raise ValueError("--max-null-replicates must be at least 2 for a Z-score")
        if args.min_null_replicates > args.max_null_replicates:
            raise ValueError("--min-null-replicates cannot exceed --max-null-replicates")
        # A capped run defines B0 as the first set of the seeded permutation and
        # the nulls as the next N; a named reference or alternative references
        # would each need their own rule, so neither is combined with a cap.
        if args.reference_replicate is not None:
            raise ValueError("--max-null-replicates cannot be combined with --reference-replicate")
        if args.reference_sensitivity:
            raise ValueError("--max-null-replicates cannot be combined with --reference-sensitivity")
    target_meta, match_meta, target_digest = _validate_provenance(args.target, args.matches)
    _validate_target_authority(target_meta, match_meta, args.a_type)
    vcf_eligibility_identity = _validate_vcf_identity(
        target_meta, match_meta, heterozygous=args.heterozygous
    )
    match_schema = match_meta.get("schema_version")
    if match_schema != "bootstrap-target-matches-v1":
        raise ValueError(
            "Wasserstein Phi-SFS requires bootstrap-target-matches-v1; "
            f"received {match_schema!r}"
        )
    if match_meta.get("phi_sfs_selection_blind") is not True:
        raise ValueError("matched controls must be selected without allele-frequency information")
    match_config = match_meta.get("config")
    if not isinstance(match_config, dict) or match_config.get("disjoint_replicates") is not True:
        raise ValueError("matched controls must be generated with --disjoint-replicates")
    if match_meta.get("maximum_control_reuse") != 1:
        raise ValueError("matched-control metadata must report maximum_control_reuse equal to 1")

    te_chrom, te_pos, target_rows, snp_chrom, snp_pos, match_rows, identifiers = _load_coordinates(
        args.target, args.matches, match_schema
    )
    identifier_names = list(identifiers)
    replicate_ids = identifiers["replicate_id"]
    if np.unique(replicate_ids).size != replicate_ids.size:
        raise ValueError("matched-control replicate_id values must be unique")
    qc_pass = np.load(args.matches / "qc_pass.npy", allow_pickle=False)
    if qc_pass.dtype.kind != "b" or qc_pass.shape != replicate_ids.shape:
        raise ValueError("qc_pass.npy must be a boolean array aligned with matched sets")
    replicate_order = np.argsort(replicate_ids, kind="stable")
    passing = replicate_order[qc_pass[replicate_order]]

    # B0 is drawn uniformly from the QC-passing sets with a seed derived from
    # --reference-seed and the target digest, unless --reference-replicate
    # names it. Replicate 0 is matched first, from the undepleted pool, so it
    # is not a typical set; a seeded draw keeps the choice prespecified,
    # reproducible and SFS-blind. That makes the choice of reference uniform
    # among accepted sets; it does not by itself make sequentially matched
    # sets exchangeable (CODE_REVIEW_ROUND12.md, finding 7).
    # The same permutation supplies the alternative references for
    # reference sensitivity, so they are prespecified in the same way.
    reference_seed = int.from_bytes(hashlib.sha256(
        f"phi-sfs-reference:{args.reference_seed}:{target_digest}".encode()
    ).digest()[:8], "little")
    order = passing[np.random.default_rng(reference_seed).permutation(passing.size)]
    if args.reference_replicate is None:
        if order.size == 0:
            raise ValueError("no matched set passes QC, so no reference can be drawn")
        reference_source_index = int(order[0])
        reference_rule = (
            f"uniform draw from QC-passing sets, seed from --reference-seed "
            f"{args.reference_seed} and target_digest"
        )
    else:
        reference_hits = np.flatnonzero(replicate_ids == args.reference_replicate)
        if reference_hits.size != 1:
            raise ValueError(
                f"reference replicate ID {args.reference_replicate} is absent from matches"
            )
        reference_source_index = int(reference_hits[0])
        if not bool(qc_pass[reference_source_index]):
            raise ValueError(
                f"reference replicate ID {args.reference_replicate} failed matching QC"
            )
        reference_rule = "prespecified --reference-replicate"
    alternatives = order[order != reference_source_index]

    # By default the nulls are every QC-passing set other than B0, so R is
    # whatever the matcher delivered, subject to a floor fixed before the run.
    # With --max-null-replicates N they are instead the next N sets of the same
    # seeded permutation, which fixes R so that (R + 1) * alpha can be an
    # integer and the nominal test has exact attainable size. QC is computed
    # from age matching alone, so neither choice ever looks at an SFS. The
    # add-one P-value (1 + exceedances) / (R + 1) is valid for any such R.
    if capped:
        if order.size < args.max_null_replicates + 1:
            raise ValueError(
                f"only {int(order.size)} matched sets pass QC; "
                f"--max-null-replicates {args.max_null_replicates} requires "
                f"{args.max_null_replicates + 1} (B0 plus the nulls)"
            )
        null_source_indices = order[1:args.max_null_replicates + 1]
        unused_source_indices = order[args.max_null_replicates + 1:]
    else:
        null_source_indices = passing[passing != reference_source_index]
        unused_source_indices = order[:0]
    null_count = int(null_source_indices.size)
    if null_count < args.min_null_replicates:
        raise ValueError(
            f"only {null_count} QC-passing null replicates remain after reserving "
            f"B0; --min-null-replicates requires {args.min_null_replicates}"
        )
    accepted_source_indices = np.sort(
        np.append(null_source_indices, reference_source_index)
    )
    reference_index = int(np.flatnonzero(
        accepted_source_indices == reference_source_index
    )[0])
    null_indices = np.delete(np.arange(accepted_source_indices.size), reference_index)

    # Reference sensitivity reruns the calibration with the next N sets of the
    # same seeded permutation as B0. Each alternative uses every other
    # QC-passing set as its nulls, so the primary B0 becomes a null there, and
    # all of them draw on the same scanned sets as the primary analysis.
    sensitivity_n = int(args.reference_sensitivity)
    if alternatives.size < sensitivity_n:
        raise ValueError(
            f"only {int(alternatives.size)} other QC-passing sets exist; "
            f"--reference-sensitivity {sensitivity_n} requires that many"
        )
    sensitivity_source_indices = alternatives[:sensitivity_n]
    sensitivity_accepted_sets: list[np.ndarray] = [
        np.sort(passing) for _ in sensitivity_source_indices.tolist()
    ]
    # Every set any analysis (primary or sensitivity) uses must be scanned, so
    # the union -- not just the primary accepted set -- determines `requested`
    # below. The primary computation still reads only its own accepted
    # sub-array afterwards, so requesting sensitivity cannot perturb it.
    used_source_indices = np.unique(
        np.concatenate([accepted_source_indices, *sensitivity_accepted_sets])
    )

    if np.unique(match_rows).size != match_rows.size:
        raise ValueError("disjoint matched bundle contains a control used more than once")
    if np.intersect1d(target_rows, match_rows).size:
        raise ValueError("matched B controls must exclude every row in focal set A")
    reuse_row_indices = _load_integers(args.matches / "reuse_row_indices.npy", "reuse row indices")
    reuse_counts = _load_integers(args.matches / "reuse_counts.npy", "reuse counts")
    if reuse_row_indices.ndim != 1 or reuse_counts.shape != reuse_row_indices.shape:
        raise ValueError("reuse_row_indices.npy and reuse_counts.npy must be aligned 1-D arrays")
    if np.any(np.diff(reuse_row_indices) <= 0):
        raise ValueError("reuse_row_indices.npy must be strictly increasing")
    if reuse_counts.size == 0 or int(reuse_counts.max()) != 1:
        raise ValueError("reuse_counts.npy must verify maximum control reuse equal to 1")

    def max_control_reuse(rows: np.ndarray) -> int:
        """Look up one matched set's rows and return their largest reuse count.

        `reuse_row_indices.npy`/`reuse_counts.npy` record, for every control
        row used in any published set, how many times it is used across the
        whole bundle. In a valid disjoint bundle that is 1 everywhere, so this
        is 1 for every set; the lookup is exact rather than assumed.
        """
        positions = np.searchsorted(reuse_row_indices, rows)
        if np.any(positions >= reuse_row_indices.size) or np.any(
            reuse_row_indices[positions] != rows
        ):
            raise ValueError("reuse_row_indices.npy does not cover every control row in use")
        return int(reuse_counts[positions].max())

    def aligned_match_array(name: str) -> np.ndarray:
        values = np.load(args.matches / f"{name}.npy", allow_pickle=False)
        if values.shape[:1] != replicate_ids.shape:
            raise ValueError(f"{name}.npy does not align with matched sets")
        return values

    match_to_bootstrap = aligned_match_array("match_to_bootstrap_w1").astype(np.float64)
    match_error_ratio = aligned_match_array("matching_error_ratio").astype(np.float64)
    bootstrap_to_observed = aligned_match_array("bootstrap_to_observed_w1").astype(np.float64)
    bootstrap_seeds = aligned_match_array("bootstrap_seeds")
    bootstrap_counts = aligned_match_array("bootstrap_counts")

    te_coordinates = list(zip(te_chrom.tolist(), te_pos.tolist()))
    all_snp_coordinates = [
        list(zip(snp_chrom[row].tolist(), snp_pos[row].tolist()))
        for row in range(snp_pos.shape[0])
    ]
    if not all_snp_coordinates:
        raise ValueError("matched-control bundle contains no sets")
    if not te_coordinates:
        raise ValueError("target contains no A sites")
    site_count = len(te_coordinates)
    if snp_pos.shape[1] != site_count:
        raise ValueError(
            f"matched sets contain {snp_pos.shape[1]} sites but A contains {site_count}; "
            "rebuild the target and controls with the shared eligibility mask"
        )
    if bootstrap_counts.shape != (replicate_ids.size, site_count):
        raise ValueError("bootstrap_counts.npy does not align with matched sets and M")
    if bootstrap_counts.dtype.kind not in "iu" or np.any(bootstrap_counts < 0):
        raise ValueError("bootstrap_counts.npy must contain nonnegative integers")
    if np.any(bootstrap_counts.sum(axis=1) != site_count):
        raise ValueError("every bootstrap count vector must sum to M")
    if bootstrap_seeds.dtype.kind not in "iu":
        raise ValueError("bootstrap_seeds.npy must contain integers")
    for name, values in (
        ("match_to_bootstrap_w1", match_to_bootstrap),
        ("matching_error_ratio", match_error_ratio),
        ("bootstrap_to_observed_w1", bootstrap_to_observed),
    ):
        selected = values[accepted_source_indices]
        if not np.all(np.isfinite(selected)) or np.any(selected < 0.0):
            raise ValueError(f"QC-passing {name}.npy values must be finite and nonnegative")
    # Scanned coordinates cover every set any analysis uses: the primary
    # accepted set and, when sensitivity is requested, every alternative
    # reference's own accepted set. Indexing by `used_source_indices` below
    # then recovers just the primary sub-array, so this union cannot change
    # the primary result.
    used_snp_coordinates = [all_snp_coordinates[index] for index in used_source_indices]
    selected_identifiers = {
        name: values[accepted_source_indices] for name, values in identifiers.items()
    }

    requested = set(te_coordinates)
    for row in used_snp_coordinates:
        requested.update(row)
    print(
        f"Scanning {args.vcf} for {len(requested):,} requested sites "
        f"across {len(used_snp_coordinates)} matched sets",
        flush=True,
    )
    table = Path(args.ancestral_table)
    store_meta = json.loads((table / "metadata.json").read_text(encoding="utf-8"))
    if store_meta.get("schema_version") != ANCESTRAL_SCHEMA_VERSION:
        raise ValueError(
            f"ancestral table at {table} has schema "
            f"{store_meta.get('schema_version')!r}, expected "
            f"{ANCESTRAL_SCHEMA_VERSION!r}; rebuild it with the current "
            "normalize_tes.build_ancestral_states"
        )
    if list(store_meta.get("bases", [])) != ["A", "C", "G", "T"]:
        raise ValueError(
            f"ancestral table declares bases {store_meta.get('bases')!r}; the "
            "count columns are read positionally and must be A, C, G, T"
        )
    # Bind the table to the same store the target and matches were built from.
    # Without this a table from another store with overlapping row coordinates
    # produces a complete, plausible result with the wrong control polarity --
    # a silent scientific error rather than a crash. The digest is the identity;
    # the recorded path supplies the row-coordinate catalog and is therefore
    # required at analysis time as well as recorded for provenance.
    table_digest = store_meta.get("store_content_sha256")
    expected_digest = (
        target_meta.get("source_store_content_sha256")
        or match_meta.get("source_store_content_sha256")
    )
    if expected_digest is None:
        raise ValueError(
            "target and matches record no store content digest, so the ancestral "
            "table cannot be bound to them; rebuild them from a store carrying "
            "content_sha256"
        )
    if table_digest != expected_digest:
        raise ValueError(
            f"ancestral table was built from a different interval store: it "
            f"records {table_digest!r}, the target and matches record "
            f"{expected_digest!r}"
        )
    source_store = Path(store_meta["store"])
    source_metadata = json.loads(
        (source_store / "metadata.json").read_text(encoding="utf-8")
    )
    if source_metadata.get("content_sha256") != table_digest:
        raise ValueError(
            "the ancestral table's recorded source-store path no longer carries "
            "the content digest used to build the table"
        )
    store_positions = np.load(
        source_store / "positions.npy", mmap_mode="r", allow_pickle=False)
    declared_rows = store_meta.get("store_rows")
    if declared_rows is not None and int(declared_rows) != store_positions.size:
        raise ValueError(
            f"ancestral table declares {declared_rows:,} store rows but its store "
            f"has {store_positions.size:,}"
        )
    chromosome_offsets = {
        entry["chrom"]: int(entry["offset"])
        for entry in source_metadata["chromosomes"]
    }
    # Validated for shape and dtype even though the resolver conditions on the
    # two observed alleles rather than this total; a malformed table should fail
    # here rather than in a later reader.
    if not store_meta.get("complete"):
        raise ValueError(
            f"ancestral table at {table} is not marked complete; it may be a "
            "partially written or interrupted build"
        )
    present_counts = _checked_table_array(
        table / "present_draw_count.npy", (store_positions.size,))
    ancestral_counts = _checked_table_array(
        table / "ancestral_counts.npy", (store_positions.size, 4))
    # Every site's polarity is read from these arrays, so bind their content to
    # the metadata before trusting either.
    verify_table_arrays(table, store_meta, ancestral_counts, present_counts)
    polarity = PolarityResolver(
        set(te_coordinates) if args.a_type == "TE" else set(),
        store_positions=store_positions,
        chromosome_offsets=chromosome_offsets,
        ancestral_counts=ancestral_counts,
        present_draw_count=present_counts,
    )
    counts, vcf_sha256 = read_site_counts(
        args.vcf,
        requested,
        polarity=polarity,
        heterozygous=args.heterozygous,
        progress=not args.quiet,
    )
    print(
        f"polarity: {polarity.te_sites:,} TE sites and "
        f"{polarity.control_sites:,} SNP sites from the ancestral table",
        flush=True,
    )
    if vcf_sha256 != vcf_eligibility_identity["vcf_sha256"]:
        raise ValueError(
            "the VCF differs from the one that defined eligibility: scanned "
            f"{args.vcf} hashes to {vcf_sha256}, but the eligibility identity "
            f"records {vcf_eligibility_identity['vcf_sha256']}"
        )
    missing = sorted(requested.difference(counts))
    if missing:
        preview = ", ".join(f"{chrom}:{pos}" for chrom, pos in missing[:10])
        raise ValueError(f"{len(missing)} requested sites are absent from the VCF: {preview}")

    site_rows, projections, endpoints = project_sites(counts)
    hard_site_rows = hard_projections = hard_endpoints = None
    if args.asymmetric_polarity_null:
        hard_site_rows, hard_projections, hard_endpoints = project_sites_bernoulli_q(
            counts, seed=args.polarity_imputation_seed,
        )
    print(
        f"Projected {projections.shape[0]:,} distinct (k, n) pairs "
        f"covering {len(site_rows):,} eligible sites",
        flush=True,
    )

    # A is polarized exactly as the null-left sets are, whatever its type, so
    # the observed and null distances differ only in which sites they hold.
    a_uses_hard_imputation = args.asymmetric_polarity_null
    a_counts, a_endpoint, a_eligible = accumulate_spectrum(
        te_coordinates,
        hard_site_rows if a_uses_hard_imputation else site_rows,
        hard_projections if a_uses_hard_imputation else projections,
        hard_endpoints if a_uses_hard_imputation else endpoints,
    )
    if a_eligible != site_count:
        raise ValueError(
            f"A retains {a_eligible} of {site_count} sites after callability filtering; "
            "rebuild the target and controls with the shared eligibility mask"
        )
    a_raw, a_normalized = normalized_spectrum(a_counts)

    b_raw_all = np.empty((len(used_snp_coordinates), PROJECTION_SIZE - 1), dtype=np.float64)
    b_normalized_all = np.empty_like(b_raw_all)
    b_endpoints_all = np.empty(len(used_snp_coordinates), dtype=np.float64)
    b_hard_raw_all = b_hard_normalized_all = b_hard_endpoints_all = None
    if args.asymmetric_polarity_null:
        b_hard_raw_all = np.empty_like(b_raw_all)
        b_hard_normalized_all = np.empty_like(b_raw_all)
        b_hard_endpoints_all = np.empty_like(b_endpoints_all)

    for position, coordinates in enumerate(used_snp_coordinates):
        counts_vector, endpoint, eligible = accumulate_spectrum(
            coordinates, site_rows, projections, endpoints
        )
        if eligible != site_count:
            replicate_id = int(replicate_ids[used_source_indices[position]])
            raise ValueError(
                f"B replicate {replicate_id} retains {eligible} of {site_count} sites "
                "after callability filtering; rebuild the target and controls with "
                "the shared eligibility mask"
            )
        raw, normalized = normalized_spectrum(counts_vector)
        b_raw_all[position] = raw
        b_normalized_all[position] = normalized
        b_endpoints_all[position] = endpoint
        if args.asymmetric_polarity_null:
            hard_counts, hard_endpoint, hard_eligible = accumulate_spectrum(
                coordinates, hard_site_rows, hard_projections, hard_endpoints
            )
            if hard_eligible != site_count:
                raise RuntimeError("Bernoulli-q projection changed SNP eligibility")
            hard_raw, hard_normalized = normalized_spectrum(hard_counts)
            b_hard_raw_all[position] = hard_raw
            b_hard_normalized_all[position] = hard_normalized
            b_hard_endpoints_all[position] = hard_endpoint

    daf = RETAINED_BINS.astype(np.float64) / PROJECTION_SIZE
    b_cdf_all = np.cumsum(b_normalized_all, axis=1)
    b_hard_cdf_all = (
        np.cumsum(b_hard_normalized_all, axis=1)
        if args.asymmetric_polarity_null else None
    )

    # The primary result is read out of the shared scan above by position, so
    # it is identical whether or not sensitivity pulled extra sets into that
    # scan: every array below depends only on this sub-selection.
    accepted_positions = np.searchsorted(used_source_indices, accepted_source_indices)
    b_raw = b_raw_all[accepted_positions]
    b_normalized = b_normalized_all[accepted_positions]
    b_endpoints = b_endpoints_all[accepted_positions]
    b_cdf = b_cdf_all[accepted_positions]
    b_hard_raw = b_hard_normalized = b_hard_endpoints = b_hard_cdf = None
    if args.asymmetric_polarity_null:
        b_hard_raw = b_hard_raw_all[accepted_positions]
        b_hard_normalized = b_hard_normalized_all[accepted_positions]
        b_hard_endpoints = b_hard_endpoints_all[accepted_positions]
        b_hard_cdf = b_hard_cdf_all[accepted_positions]
    reference_sfs = b_normalized[reference_index]
    reference_cdf = b_cdf[reference_index]
    observed_result = phi_sfs(a_normalized, reference_sfs, daf=daf)
    null_cdf = b_hard_cdf if args.asymmetric_polarity_null else b_cdf
    null_sfs = b_hard_normalized if args.asymmetric_polarity_null else b_normalized
    null_phi = (
        np.abs(null_cdf[null_indices, :-1] - reference_cdf[:-1])
        * np.diff(daf)
    ).sum(axis=1)
    calibration = calibrate_phi(observed_result.value, null_phi)
    null_mean_daf_difference = (
        (null_sfs[null_indices] - reference_sfs) @ daf
    )

    # Reference sensitivity: rerun the same phi_sfs/calibrate_phi pair once
    # per alternative B0, each using every other QC-passing set as its nulls
    # (see the selection above), reading out of the same shared scan.
    sensitivity_reference_ids: list[int] = []
    sensitivity_observed_phi_sfs: list[float] = []
    sensitivity_z_scores: list[float] = []
    sensitivity_p_values: list[float] = []
    for alt_source, alt_accepted in zip(
        sensitivity_source_indices.tolist(), sensitivity_accepted_sets
    ):
        alt_positions = np.searchsorted(used_source_indices, alt_accepted)
        alt_normalized = b_normalized_all[alt_positions]
        alt_cdf = b_cdf_all[alt_positions]
        alt_reference_local = int(np.flatnonzero(alt_accepted == alt_source)[0])
        alt_null_local = np.delete(np.arange(alt_accepted.size), alt_reference_local)
        alt_reference_sfs = alt_normalized[alt_reference_local]
        alt_reference_cdf = alt_cdf[alt_reference_local]
        alt_result = phi_sfs(a_normalized, alt_reference_sfs, daf=daf)
        alt_null_cdf = (
            b_hard_cdf_all[alt_positions]
            if args.asymmetric_polarity_null else alt_cdf
        )
        alt_null_phi = (
            np.abs(alt_null_cdf[alt_null_local, :-1] - alt_reference_cdf[:-1])
            * np.diff(daf)
        ).sum(axis=1)
        alt_calibration = calibrate_phi(alt_result.value, alt_null_phi)
        sensitivity_reference_ids.append(int(replicate_ids[alt_source]))
        sensitivity_observed_phi_sfs.append(alt_result.value)
        sensitivity_z_scores.append(alt_calibration.z_score)
        sensitivity_p_values.append(alt_calibration.p_value)

    selected_rows = match_rows[accepted_source_indices]
    reference_rows = selected_rows[reference_index]
    overlaps = np.asarray([
        np.intersect1d(selected_rows[index], reference_rows).size
        for index in null_indices
    ], dtype=np.int64)
    if np.any(overlaps != 0):
        raise ValueError("a null matched set overlaps B0 despite disjoint mode")

    selected_source = accepted_source_indices
    reference_source = int(selected_source[reference_index])
    reference_id = int(selected_identifiers["replicate_id"][reference_index])
    a_retained_mass = float(a_raw.sum())
    reference_retained_mass = float(b_raw[reference_index].sum())
    reference_max_control_reuse = max_control_reuse(reference_rows)
    mixture_polarity_rule = (
        "posterior q*h(k,n) + (1-q)*h(n-k,n) over usable ARG draws"
    )
    hard_polarity_rule = "coordinate-keyed hard Bernoulli(q) orientation"
    a_polarity_rule = (
        hard_polarity_rule if args.asymmetric_polarity_null else mixture_polarity_rule
    )
    null_left_polarity_rule = (
        hard_polarity_rule if args.asymmetric_polarity_null else mixture_polarity_rule
    )
    null_polarity_design = (
        ASYMMETRIC_NULL_DESIGN
        if args.asymmetric_polarity_null else SYMMETRIC_NULL_DESIGN
    )

    comparison_rows: list[dict[str, object]] = [{
        "role": "observed",
        "left_id": "A",
        "left_type": args.a_type,
        "right_id": reference_id,
        "right_type": args.b_type,
        "left_polarity_rule": a_polarity_rule,
        "right_polarity_rule": mixture_polarity_rule,
        "phi_sfs": observed_result.value,
        "null_z_score": "",
        "mean_daf_difference": observed_result.mean_daf_difference,
        "left_sites": site_count,
        "right_sites": site_count,
        "left_retained_mass": a_retained_mass,
        "right_retained_mass": reference_retained_mass,
        "left_endpoint_mass": a_endpoint,
        "right_endpoint_mass": b_endpoints[reference_index],
        "left_matching_qc_pass": "",
        "left_bootstrap_to_observed_w1": "",
        "left_match_to_bootstrap_w1": "",
        "left_matching_error_ratio": "",
        "right_matching_qc_pass": True,
        "right_bootstrap_to_observed_w1": bootstrap_to_observed[reference_source],
        "right_match_to_bootstrap_w1": match_to_bootstrap[reference_source],
        "right_matching_error_ratio": match_error_ratio[reference_source],
        "overlap_with_reference": "",
        "left_max_control_reuse": "",
        "right_max_control_reuse": reference_max_control_reuse,
    }]
    for null_offset, accepted_index in enumerate(null_indices):
        source_index = int(selected_source[accepted_index])
        comparison_rows.append({
            "role": "null",
            "left_id": int(selected_identifiers["replicate_id"][accepted_index]),
            "left_type": args.b_type,
            "right_id": reference_id,
            "right_type": args.b_type,
            "left_polarity_rule": null_left_polarity_rule,
            "right_polarity_rule": mixture_polarity_rule,
            "phi_sfs": calibration.null[null_offset],
            "null_z_score": calibration.null_z_scores[null_offset],
            "mean_daf_difference": null_mean_daf_difference[null_offset],
            "left_sites": site_count,
            "right_sites": site_count,
            "left_retained_mass": float(
                (b_hard_raw if args.asymmetric_polarity_null else b_raw)[accepted_index].sum()
            ),
            "right_retained_mass": reference_retained_mass,
            "left_endpoint_mass": (
                b_hard_endpoints if args.asymmetric_polarity_null else b_endpoints
            )[accepted_index],
            "right_endpoint_mass": b_endpoints[reference_index],
            "left_matching_qc_pass": True,
            "left_bootstrap_to_observed_w1": bootstrap_to_observed[source_index],
            "left_match_to_bootstrap_w1": match_to_bootstrap[source_index],
            "left_matching_error_ratio": match_error_ratio[source_index],
            "right_matching_qc_pass": True,
            "right_bootstrap_to_observed_w1": bootstrap_to_observed[reference_source],
            "right_match_to_bootstrap_w1": match_to_bootstrap[reference_source],
            "right_matching_error_ratio": match_error_ratio[reference_source],
            "overlap_with_reference": int(overlaps[null_offset]),
            "left_max_control_reuse": max_control_reuse(selected_rows[accepted_index]),
            "right_max_control_reuse": reference_max_control_reuse,
        })

    output = args.output
    if output.exists():
        raise FileExistsError(f"output already exists: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{output.name}.tmp.", dir=output.parent))
    try:
        arrays = {
            "bins.npy": RETAINED_BINS,
            "daf.npy": daf,
            "a_raw_sfs.npy": a_raw,
            "a_normalized_sfs.npy": a_normalized,
            "a_cdf.npy": observed_result.cdf_a,
            "b_raw_sfs.npy": b_raw,
            "b_normalized_sfs.npy": b_normalized,
            "b_cdf.npy": b_cdf,
            "reference_sfs.npy": reference_sfs,
            "reference_cdf.npy": reference_cdf,
            "observed_phi_sfs.npy": np.asarray(observed_result.value),
            "null_phi_sfs.npy": calibration.null,
            "null_z_scores.npy": calibration.null_z_scores,
            "observed_bin_residual.npy": observed_result.bin_residual,
            "observed_cdf_residual.npy": observed_result.cdf_residual,
            "reference_replicate_id.npy": np.asarray(reference_id, dtype=np.int64),
            "b_bootstrap_seeds.npy": bootstrap_seeds[selected_source],
            "b_bootstrap_counts.npy": bootstrap_counts[selected_source],
            **{
                f"b_{name}.npy": values
                for name, values in selected_identifiers.items()
            },
            **{
                f"null_{name}.npy": values[null_indices]
                for name, values in selected_identifiers.items()
            },
            "sensitivity_reference_ids.npy": np.asarray(
                sensitivity_reference_ids, dtype=np.int64
            ),
            "sensitivity_observed_phi_sfs.npy": np.asarray(
                sensitivity_observed_phi_sfs, dtype=np.float64
            ),
            "sensitivity_z_scores.npy": np.asarray(
                sensitivity_z_scores, dtype=np.float64
            ),
            "sensitivity_p_values.npy": np.asarray(
                sensitivity_p_values, dtype=np.float64
            ),
        }
        if args.asymmetric_polarity_null:
            arrays.update({
                "b_bernoulli_q_raw_sfs.npy": b_hard_raw,
                "b_bernoulli_q_normalized_sfs.npy": b_hard_normalized,
                "b_bernoulli_q_cdf.npy": b_hard_cdf,
            })
        for name, values in arrays.items():
            np.save(staging / name, values, allow_pickle=False)

        summary_row = {
            "focal_label": args.target.name,
            "target_digest": target_digest,
            "a_type": args.a_type,
            "b_type": args.b_type,
            "site_count_m": site_count,
            "null_replicates_r": null_count,
            "reference_replicate_id": reference_id,
            "observed_phi_sfs": calibration.observed,
            "null_mean": calibration.null_mean,
            "null_sample_sd": calibration.null_sd,
            "z_score": calibration.z_score,
            "exceedances": calibration.exceedances,
            "p_value": calibration.p_value,
            "minimum_attainable_p": 1.0 / (null_count + 1.0),
            "mean_daf_difference": observed_result.mean_daf_difference,
            "null_polarity_design": null_polarity_design,
            "polarity_imputation_seed": (
                args.polarity_imputation_seed if args.asymmetric_polarity_null else ""
            ),
            "reference_sensitivity_n": sensitivity_n,
            "reference_sensitivity_z_min": (
                min(sensitivity_z_scores) if sensitivity_z_scores else ""
            ),
            "reference_sensitivity_z_max": (
                max(sensitivity_z_scores) if sensitivity_z_scores else ""
            ),
            "reference_sensitivity_p_min": (
                min(sensitivity_p_values) if sensitivity_p_values else ""
            ),
            "reference_sensitivity_p_max": (
                max(sensitivity_p_values) if sensitivity_p_values else ""
            ),
        }
        with (staging / "summary.csv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(summary_row))
            writer.writeheader()
            writer.writerow(summary_row)
        with (staging / "comparisons.csv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(comparison_rows[0]))
            writer.writeheader()
            writer.writerows(comparison_rows)

        metadata = {
            "schema_version": SCHEMA_VERSION,
            "complete": True,
            "software": software_provenance(),
            "creation_command": " ".join(sys.argv),
            "creation_time_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
            "numpy_version": np.__version__,
            "projection_size": PROJECTION_SIZE,
            "retained_bins": [1, 19],
            "daf_grid": daf.tolist(),
            "site_projection_renormalized": False,
            "final_spectra_normalized": True,
            "phi_definition": "sum(abs(cdf_a[:-1] - cdf_b[:-1]) * diff(daf))",
            "phi_interpretation": (
                "one-dimensional Wasserstein distance between projected normalized SFS"
            ),
            "a_type": args.a_type,
            "b_type": args.b_type,
            "target": str(args.target.resolve()),
            "target_schema_version": target_meta.get("schema_version"),
            "matches": str(args.matches.resolve()),
            "matches_schema_version": match_schema,
            "replicate_identifiers": identifier_names,
            "target_digest": target_digest,
            "vcf": str(args.vcf.resolve()),
            "vcf_sha256": vcf_sha256,
            "vcf_eligibility_identity": vcf_eligibility_identity,
            "a_polarity_rule": a_polarity_rule,
            "b_polarity_rule": mixture_polarity_rule,
            "null_left_polarity_rule": null_left_polarity_rule,
            "null_polarity_design": null_polarity_design,
            "asymmetric_polarity_null": args.asymmetric_polarity_null,
            "polarity_imputation_seed": (
                args.polarity_imputation_seed if args.asymmetric_polarity_null else None
            ),
            "polarity_imputation_algorithm": (
                "sha256(phi-sfs-bernoulli-q-v1, seed, chromosome, position)"
                if args.asymmetric_polarity_null else None
            ),
            "ancestral_table": str(Path(args.ancestral_table).resolve()),
            "ancestral_table_schema_version": store_meta.get("schema_version"),
            "ancestral_table_array_sha256": store_meta.get("array_sha256"),
            "te_sites_polarized": polarity.te_sites,
            "te_usable_arg_draws": polarity.te_usable_draws,
            "te_unusable_arg_draws": polarity.te_unusable_draws,
            "te_polarity_source": "ARG posterior ancestral table, as for SNPs",
            "snp_sites_polarized": polarity.control_sites,
            "snp_usable_arg_draws": polarity.control_usable_draws,
            "snp_unusable_arg_draws": polarity.control_unusable_draws,
            "ancestral_case_policy": (
                "case-sensitive; a lowercase ancestral allele is rejected rather than folded"
            ),
            "heterozygous_policy": args.heterozygous,
            "biallelic_policy": (
                "records are assumed biallelic; a comma in ALT is rejected"
            ),
            "filter_policy": (
                "the VCF FILTER column is ignored; every record at a requested "
                "coordinate is used"
            ),
            "reference_selection_rule": reference_rule + ", chosen before the SFS scan",
            "reference_seed": args.reference_seed,
            "reference_replicate_id": reference_id,
            "reference_index_in_b_arrays": reference_index,
            "reference_bootstrap_seed": int(bootstrap_seeds[reference_source]),
            "reference_bootstrap_counts_array": "b_bootstrap_counts.npy",
            "minimum_null_replicates": args.min_null_replicates,
            "accepted_null_replicates": null_count,
            "null_selection_rule": (
                f"the {args.max_null_replicates} QC-passing sets following B0 in "
                "the seeded permutation"
                if capped else "every QC-passing non-reference set"
            ),
            "max_null_replicates": args.max_null_replicates,
            "qc_passing_sets": int(passing.size),
            "selection_permutation_algorithm": (
                "QC-passing sets ordered by replicate_id, permuted by "
                "numpy.random.default_rng(seed).permutation, seed = first 8 bytes "
                "(little-endian) of sha256('phi-sfs-reference:<reference_seed>:"
                "<target_digest>'); B0 is the first set"
            ),
            "selected_null_replicate_ids": sorted(
                int(value) for value in replicate_ids[null_source_indices]
            ),
            "unused_qc_passing_replicate_ids": sorted(
                int(value) for value in replicate_ids[unused_source_indices]
            ),
            "matched_sets_published": int(replicate_ids.size),
            "matched_sets_failing_qc": int(np.count_nonzero(~qc_pass)),
            "null_standard_deviation_ddof": 1,
            "p_value_tail_rule": "null distance >= observed distance",
            "p_value_formula": "(1 + exceedances) / (R + 1)",
            "observed_phi_sfs": calibration.observed,
            "null_mean": calibration.null_mean,
            "null_sample_sd": calibration.null_sd,
            "z_score": calibration.z_score,
            "p_value": calibration.p_value,
            "exceedances": calibration.exceedances,
            "matched_sets_passing_qc_including_reference": int(selected_source.size),
            "distinct_projections": int(projections.shape[0]),
            "equal_eligible_site_count": site_count,
            "a_input_sites": len(te_coordinates),
            "a_eligible_sites": a_eligible,
            "a_retained_mass": a_retained_mass,
            "a_endpoint_mass": a_endpoint,
            "a_retained_fraction": a_retained_mass / a_eligible,
            "a_endpoint_fraction": a_endpoint / a_eligible,
            "matching_qc_rule": "source qc_pass.npy; only passing sets analyzed",
            "matching_algorithm_same_for_reference_and_null": True,
            "matching_phi_sfs_selection_blind": True,
            "disjoint_replicates": True,
            "maximum_control_reuse": 1,
            "maximum_overlap_with_reference": int(overlaps.max(initial=0)),
            "bootstrap_to_observed_w1_range": [
                float(bootstrap_to_observed[selected_source].min()),
                float(bootstrap_to_observed[selected_source].max()),
            ],
            "matching_error_ratio_range": [
                float(match_error_ratio[selected_source].min()),
                float(match_error_ratio[selected_source].max()),
            ],
            "reference_sensitivity_run": sensitivity_n > 0,
            "reference_sensitivity_n": sensitivity_n,
            "reference_sensitivity_reference_ids": sensitivity_reference_ids,
            "reference_sensitivity_rule": (
                "alternative references are the next N sets of the seeded "
                "permutation of QC-passing sets that chose B0; each uses every "
                "other QC-passing set as its nulls, so the primary B0 becomes a null"
            ),
            "target_source_store_content_sha256": target_meta.get("source_store_content_sha256"),
            "matches_source_store_content_sha256": match_meta.get("source_store_content_sha256"),
        }
        with (staging / "metadata.json").open("w", encoding="utf-8") as handle:
            json.dump(metadata, handle, indent=2, sort_keys=True)
            handle.write("\n")
        os.replace(staging, output)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", type=Path, required=True,
                        help="TE target directory from normalize_tes.te_age_target")
    parser.add_argument("--matches", type=Path, required=True,
                        help="matched control sets from normalize_tes.bootstrap_target_matcher")
    parser.add_argument("--vcf", type=Path, required=True,
                        help="biallelic VCF of genotypes. Every requested site, TE "
                             "and control alike, must be present or the run stops")
    parser.add_argument("--output", type=Path, required=True,
                        help="destination directory for spectra and scores")
    parser.add_argument(
        "-A", "--a-type", choices=("TE", "SNP"), default="TE",
        help="variant type of focal set A (default: TE)",
    )
    parser.add_argument(
        "-B", "--b-type", choices=("SNP",), default="SNP",
        help="variant type of every matched control set B (default: SNP)",
    )
    parser.add_argument(
        "--reference-replicate", type=int, default=None,
        help="matched replicate ID to hold fixed as B0; by default B0 is drawn "
             "uniformly from the QC-passing sets (see --reference-seed)",
    )
    parser.add_argument(
        "--reference-seed", type=int, default=1002,
        help="seed, combined with the target digest, for drawing B0 and the "
             "reference-sensitivity alternatives (default: 1002)",
    )
    parser.add_argument(
        "--min-null-replicates", type=int, default=900,
        help="floor on R: every QC-passing non-reference set is a null, and "
             "the run fails if fewer than this many pass (default: 900)",
    )
    parser.add_argument(
        "--max-null-replicates", type=int, default=None,
        help="use exactly N nulls: the N QC-passing sets that follow B0 in the "
             "seeded permutation. Fails if fewer than N + 1 sets pass QC. "
             "Default: every QC-passing non-reference set",
    )
    parser.add_argument(
        "--reference-sensitivity", type=int, default=0,
        help="repeat calibration with N alternative references: the next N "
             "sets of the seeded permutation that chose B0 "
             "(default: 0, no sensitivity run)",
    )
    parser.add_argument(
        "--ancestral-table", type=Path, required=True,
        help="directory written by normalize_tes.build_ancestral_states, giving "
             "every site, TE and SNP, its posterior polarity",
    )
    parser.add_argument(
        "--asymmetric-polarity-null", action=argparse.BooleanOptionalAction,
        default=False,
        help="hard-polarize A and each null-left set with one Bernoulli(q) draw "
             "per site while keeping B0 the posterior q-mixture. The default, "
             "--no-asymmetric-polarity-null, keeps A, B0 and every null-left "
             "set as posterior q-mixtures",
    )
    parser.add_argument(
        "--polarity-imputation-seed", type=int, default=2001,
        help="seed for coordinate-keyed Bernoulli(q) hard orientations "
             "(default: 2001; used only with --asymmetric-polarity-null)",
    )
    parser.add_argument("--heterozygous", choices=("error", "missing"), default="error",
                        help="how to treat a heterozygous call in these inbred "
                             "lines: stop, or count the genotype as missing")
    parser.add_argument(
        "--quiet", action="store_true",
        help="suppress periodic VCF scan progress",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    calculate(args)
    print(f"Wrote {args.output}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
