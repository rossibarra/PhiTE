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
  consulted. TE sites are polarized by biology -- an insertion is the derived
  state -- and SNPs in either A or B by the ARG-derived table given to
  `--ancestral-table`, as a posterior-weighted mixture over the two observed
  alleles.
"""

from __future__ import annotations

import argparse
import csv
import datetime
import gzip
import hashlib
import io
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

from .release_provenance import software_provenance
from .sample_age_matched_controls import _load_target, _sha256_arrays


SCHEMA_VERSION = "phi-sfs-wasserstein-v1"

# Per-replicate identifier arrays published by each supported matched-control
# schema. The swap sampler saves ten correlated states from each of ten chains,
# so its replicates are identified by chain and position within that chain. The
# bootstrap-target matcher produces replicates with no chain structure, so it
# identifies them by replicate alone; inventing chain and sample columns would
# imply a within-chain correlation that does not exist. "No chain structure" is
# the whole claim -- these replicates are not statistically independent, since
# they share the observed TE sample and the interval store.
MATCH_IDENTIFIERS = {
    "swap-age-matched-controls-v1": ("chain_index", "sample_index"),
    "bootstrap-target-matches-v1": ("replicate_id",),
}

PROJECTION_SIZE = 20
RETAINED_BINS = np.arange(1, PROJECTION_SIZE, dtype=np.int64)
COMPRESSED_SUFFIXES = (".gz", ".bgz", ".bgzf")
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
    `p_alt_derived` is the probability that ALT is the derived allele: exactly 1
    at a TE site, where insertion is derived by biology, and the ARG's posterior
    proportion at a control SNP. Orientation is applied when spectra are summed,
    not here, because that is the only step that depends on it.
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


class _HashingStream(io.RawIOBase):
    """Raw byte stream that digests everything read through it."""

    def __init__(self, handle):
        self._handle = handle
        self.digest = hashlib.sha256()

    def readable(self) -> bool:
        return True

    def readinto(self, buffer) -> int:
        read = self._handle.readinto(buffer)
        if read:
            self.digest.update(memoryview(buffer)[:read])
        return read

    def close(self) -> None:
        try:
            self._handle.close()
        finally:
            super().close()


def _open_vcf(path: Path):
    """Open a VCF as text over a hashing stream, so one pass yields both.

    Returns the text handle, the hashing stream, and the buffered byte stream,
    so that the caller can drain any bytes the text layer did not consume
    before reading the digest.
    """
    hashing = _HashingStream(path.open("rb"))
    buffered = io.BufferedReader(hashing, buffer_size=1 << 20)
    compressed = path.suffix.lower() in COMPRESSED_SUFFIXES
    stream = gzip.GzipFile(fileobj=buffered) if compressed else buffered
    return io.TextIOWrapper(stream, encoding="utf-8"), hashing, buffered



def _decode_genotype(gt: str, heterozygous: str) -> int | None | str:
    """Return one individual's allele, None when not callable, or an error tag.

    Each inbred individual contributes a single observed allele, so haploid and
    homozygous diploid calls are accepted and any missing allele makes the
    whole individual uncallable. Results are cached by the caller because
    genotype strings are drawn from a very small alphabet.
    """
    alleles = gt.replace("|", "/").split("/")
    if not alleles or any(allele == "." for allele in alleles):
        return None
    values = []
    for allele in alleles:
        if not allele.isdigit():
            return "invalid"
        value = int(allele)
        if value not in (0, 1):
            return "non-biallelic"
        values.append(value)
    if len(set(values)) > 1:
        return None if heterozygous == "missing" else "heterozygous"
    return values[0]


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
    """Resolves P(ALT is derived) per site, from two sources.

    TE sites are polarized by biology: a TE insertion is the derived state, and
    the genotyping convention encodes presence as ALT, so the weight is exactly
    1. The rare exception -- a TE that fixed and was later removed by a deletion
    -- is not modelled; it accounts for about 3% of TE sites and identifying it
    would need an independent outgroup.

    Control SNPs cannot be polarized that way. Their weight is the posterior
    proportion of ARG draws calling REF ancestral -- equivalently, ALT derived --
    among the draws that named one of the two observed alleles. Conditioning that
    way rather than on the raw present-draw count is what lets every requested
    site carry a weight without an intersection across draws or a fallback rule,
    and it discards draws naming a third base, which cannot orient the site.

    The proportion is used as reported. Against TE ground truth the ARG is only
    about 91% correct where all its draws agree, so this weight is somewhat
    overconfident and control spectra come out sharper than the ARG's measured
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
        self.control_sites = 0
        self.control_usable_draws = 0
        self.control_unusable_draws = 0

    def __call__(self, chrom: str, position: int, ref: str, alt: str) -> float:
        if (chrom, position) in self._te:
            self.te_sites += 1
            return 1.0
        offset = self._offsets.get(chrom)
        if offset is None:
            raise ValueError(f"no chromosome offset for {chrom!r}")
        target = float(offset + position)
        index = int(np.searchsorted(self._positions, target))
        if index >= self._positions.size or self._positions[index] != target:
            raise ValueError(
                f"control site {chrom}:{position} is absent from the ancestral "
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
        self.control_sites += 1
        self.control_usable_draws += int(oriented)
        self.control_unusable_draws += int(self._present[index] - row.sum())
        self.control_unusable_draws += int(row.sum() - oriented)
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
        while buffered.read(1 << 20):
            pass
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

    The per-replicate identifier arrays depend on the matched-control schema;
    see MATCH_IDENTIFIERS. They are returned as a name-to-array mapping and
    carried through to the outputs unchanged.
    """
    names = MATCH_IDENTIFIERS.get(schema)
    if names is None:
        supported = ", ".join(sorted(MATCH_IDENTIFIERS))
        raise ValueError(
            f"unsupported matched-control schema_version {schema!r}; "
            f"expected one of: {supported}"
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
        for name in names
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
    return te_chromosomes, te_positions, match_chromosomes, positions, identifiers


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


def calculate(args: argparse.Namespace) -> None:
    if args.null_replicates < 2:
        raise ValueError("--null-replicates must be at least 2 for a Z-score")
    target_meta, match_meta, target_digest = _validate_provenance(args.target, args.matches)
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

    te_chrom, te_pos, snp_chrom, snp_pos, identifiers = _load_coordinates(
        args.target, args.matches, match_schema
    )
    identifier_names = list(identifiers)
    replicate_ids = identifiers["replicate_id"]
    if np.unique(replicate_ids).size != replicate_ids.size:
        raise ValueError("matched-control replicate_id values must be unique")
    reference_hits = np.flatnonzero(replicate_ids == args.reference_replicate)
    if reference_hits.size != 1:
        raise ValueError(
            f"reference replicate ID {args.reference_replicate} is absent from matches"
        )

    qc_pass = np.load(args.matches / "qc_pass.npy", allow_pickle=False)
    if qc_pass.dtype.kind != "b" or qc_pass.shape != replicate_ids.shape:
        raise ValueError("qc_pass.npy must be a boolean array aligned with matched sets")
    reference_source_index = int(reference_hits[0])
    if not bool(qc_pass[reference_source_index]):
        raise ValueError(f"reference replicate ID {args.reference_replicate} failed matching QC")
    # R is fixed so every category is calibrated with equal precision. The
    # matcher publishes spare sets; the nulls are the first R QC-passing
    # non-reference sets in replicate_id order. QC is age-only, so this choice
    # is made before and independently of any SFS. Unused spares are the
    # latest, most pool-depleted replicates.
    replicate_order = np.argsort(replicate_ids, kind="stable")
    passing_nulls = replicate_order[
        qc_pass[replicate_order] & (replicate_order != reference_source_index)
    ]
    available_null_count = int(passing_nulls.size)
    null_count = int(args.null_replicates)
    if available_null_count < null_count:
        raise ValueError(
            f"only {available_null_count} QC-passing null replicates remain after "
            f"reserving B0; --null-replicates requires {null_count}. Publish more "
            "spare matched sets rather than lowering R after the fact"
        )
    null_source_indices = passing_nulls[:null_count]
    unused_passing_ids = replicate_ids[passing_nulls[null_count:]]
    accepted_source_indices = np.sort(
        np.append(null_source_indices, reference_source_index)
    )
    reference_index = int(np.flatnonzero(
        accepted_source_indices == reference_source_index
    )[0])
    null_indices = np.delete(np.arange(accepted_source_indices.size), reference_index)

    match_rows = _load_integers(args.matches / "row_indices.npy", "matched row indices")
    if match_rows.shape != snp_pos.shape:
        raise ValueError("matched row indices do not align with matched positions")
    if np.unique(match_rows).size != match_rows.size:
        raise ValueError("disjoint matched bundle contains a control used more than once")
    target_rows = _load_integers(args.target / "te_row_indices.npy", "target row indices")
    if np.intersect1d(target_rows, match_rows).size:
        raise ValueError("matched B controls must exclude every row in focal set A")
    reuse_counts = _load_integers(args.matches / "reuse_counts.npy", "reuse counts")
    if reuse_counts.ndim != 1 or reuse_counts.size == 0 or int(reuse_counts.max()) != 1:
        raise ValueError("reuse_counts.npy must verify maximum control reuse equal to 1")

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
    snp_coordinates = [all_snp_coordinates[index] for index in accepted_source_indices]
    selected_identifiers = {
        name: values[accepted_source_indices] for name, values in identifiers.items()
    }

    requested = set(te_coordinates)
    for row in snp_coordinates:
        requested.update(row)
    print(
        f"Scanning {args.vcf} for {len(requested):,} requested sites "
        f"across {len(snp_coordinates)} matched sets",
        flush=True,
    )
    table = Path(args.ancestral_table)
    store_meta = json.loads((table / "metadata.json").read_text(encoding="utf-8"))
    if store_meta.get("schema_version") != "ancestral-state-counts-v1":
        raise ValueError(f"unexpected ancestral table schema in {table}")
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
    polarity = PolarityResolver(
        set(te_coordinates) if args.a_type == "TE" else set(),
        store_positions=store_positions,
        chromosome_offsets=chromosome_offsets,
        ancestral_counts=_checked_table_array(
            table / "ancestral_counts.npy", (store_positions.size, 4)),
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
        f"polarity: {polarity.te_sites:,} TE sites from biology, "
        f"{polarity.control_sites:,} SNP sites from the ancestral table",
        flush=True,
    )
    missing = sorted(requested.difference(counts))
    if missing:
        preview = ", ".join(f"{chrom}:{pos}" for chrom, pos in missing[:10])
        raise ValueError(f"{len(missing)} requested sites are absent from the VCF: {preview}")

    site_rows, projections, endpoints = project_sites(counts)
    print(
        f"Projected {projections.shape[0]:,} distinct (k, n) pairs "
        f"covering {len(site_rows):,} eligible sites",
        flush=True,
    )

    a_counts, a_endpoint, a_eligible = accumulate_spectrum(
        te_coordinates, site_rows, projections, endpoints
    )
    if a_eligible != site_count:
        raise ValueError(
            f"A retains {a_eligible} of {site_count} sites after callability filtering; "
            "rebuild the target and controls with the shared eligibility mask"
        )
    a_raw, a_normalized = normalized_spectrum(a_counts)

    b_raw = np.empty((len(snp_coordinates), PROJECTION_SIZE - 1), dtype=np.float64)
    b_normalized = np.empty_like(b_raw)
    b_endpoints = np.empty(len(snp_coordinates), dtype=np.float64)

    for replicate, coordinates in enumerate(snp_coordinates):
        counts_vector, endpoint, eligible = accumulate_spectrum(
            coordinates, site_rows, projections, endpoints
        )
        if eligible != site_count:
            replicate_id = int(selected_identifiers["replicate_id"][replicate])
            raise ValueError(
                f"B replicate {replicate_id} retains {eligible} of {site_count} sites "
                "after callability filtering; rebuild the target and controls with "
                "the shared eligibility mask"
            )
        raw, normalized = normalized_spectrum(counts_vector)
        b_raw[replicate] = raw
        b_normalized[replicate] = normalized
        b_endpoints[replicate] = endpoint

    daf = RETAINED_BINS.astype(np.float64) / PROJECTION_SIZE
    b_cdf = np.cumsum(b_normalized, axis=1)
    reference_sfs = b_normalized[reference_index]
    reference_cdf = b_cdf[reference_index]
    observed_result = phi_sfs(a_normalized, reference_sfs, daf=daf)
    null_phi = (
        np.abs(b_cdf[null_indices, :-1] - reference_cdf[:-1])
        * np.diff(daf)
    ).sum(axis=1)
    calibration = calibrate_phi(observed_result.value, null_phi)
    null_mean_daf_difference = (
        (b_normalized[null_indices] - reference_sfs) @ daf
    )

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

    comparison_rows: list[dict[str, object]] = [{
        "role": "observed",
        "left_id": "A",
        "left_type": args.a_type,
        "right_id": reference_id,
        "right_type": args.b_type,
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
    }]
    for null_offset, accepted_index in enumerate(null_indices):
        source_index = int(selected_source[accepted_index])
        comparison_rows.append({
            "role": "null",
            "left_id": int(selected_identifiers["replicate_id"][accepted_index]),
            "left_type": args.b_type,
            "right_id": reference_id,
            "right_type": args.b_type,
            "phi_sfs": calibration.null[null_offset],
            "null_z_score": calibration.null_z_scores[null_offset],
            "mean_daf_difference": null_mean_daf_difference[null_offset],
            "left_sites": site_count,
            "right_sites": site_count,
            "left_retained_mass": float(b_raw[accepted_index].sum()),
            "right_retained_mass": reference_retained_mass,
            "left_endpoint_mass": b_endpoints[accepted_index],
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
        }
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
            "a_polarity_rule": (
                "insertion presence is derived after upstream at-least-50%-derived "
                "retention" if args.a_type == "TE" else
                "posterior q*h(k,n) + (1-q)*h(n-k,n) over usable ARG draws"
            ),
            "b_polarity_rule": (
                "posterior q*h(k,n) + (1-q)*h(n-k,n) over usable ARG draws"
            ),
            "ancestral_table": str(Path(args.ancestral_table).resolve()),
            "ancestral_table_schema_version": store_meta.get("schema_version"),
            "te_sites_polarized": polarity.te_sites,
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
            "reference_selection_rule": "prespecified --reference-replicate before SFS scan",
            "reference_replicate_id": reference_id,
            "reference_index_in_b_arrays": reference_index,
            "reference_bootstrap_seed": int(bootstrap_seeds[reference_source]),
            "reference_bootstrap_counts_array": "b_bootstrap_counts.npy",
            "requested_null_replicates": null_count,
            "accepted_null_replicates": null_count,
            "null_selection_rule": (
                "first R QC-passing non-reference sets in replicate_id order"
            ),
            "matched_sets_published": int(replicate_ids.size),
            "matched_sets_failing_qc": int(np.count_nonzero(~qc_pass)),
            "qc_passing_null_sets_available": available_null_count,
            "unused_qc_passing_replicate_ids": unused_passing_ids.tolist(),
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
            "reference_sensitivity_run": False,
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
        "--reference-replicate", type=int, default=0,
        help="prespecified matched replicate ID to hold fixed as B0 (default: 0)",
    )
    parser.add_argument(
        "--null-replicates", type=int, default=1000,
        help="exact number R of null B_i sets: the first R QC-passing "
             "non-reference sets in replicate_id order; fails if fewer pass "
             "(default: 1000)",
    )
    parser.add_argument(
        "--ancestral-table", type=Path, required=True,
        help="directory written by normalize_tes.build_ancestral_states, giving SNP "
             "sites their posterior polarity; A uses it when --a-type SNP",
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
