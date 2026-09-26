#!/usr/bin/env python3
"""Build a reusable pre-matching VCF callability mask.

The age matcher must choose TE targets and SNP controls from sites that the
SFS can actually use.  This module scans the analysis VCF once, applies the
polarity-independent genotype rules, and publishes the interval-store rows
with at least ``m`` callable inbred individuals.  Polarity deliberately does
not appear here: TE majority support and the SNP posterior mixture are applied
by their respective upstream/downstream components.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .release_provenance import software_provenance
from .snp_age_store import open_snp_age_store, store_schema
from .vcf_io import (
    COMPRESSED_SUFFIXES,
    _HashingStream,
    _open_vcf,
    decode_inbred_genotype,
    drain,
)


SCHEMA_VERSION = "vcf-eligibility-v1"
DEFAULT_MIN_CALLABLE = 20


def _sha256_array(values: np.ndarray) -> str:
    array = np.ascontiguousarray(values)
    digest = hashlib.sha256()
    digest.update(str(array.dtype).encode("utf-8"))
    digest.update(str(array.shape).encode("utf-8"))
    digest.update(array.tobytes())
    return digest.hexdigest()


@dataclass(frozen=True)
class EligibilityResult:
    """Eligible store rows and their polarity-independent VCF counts."""

    rows: np.ndarray
    alt_counts: np.ndarray
    callable_counts: np.ndarray
    report: dict
    snp_rows: np.ndarray | None = None
    p_alt_derived: np.ndarray | None = None


@dataclass(frozen=True)
class LoadedEligibility:
    """Authenticated per-row inputs for one variant type and either A/B role."""

    rows: np.ndarray
    alt_counts: np.ndarray
    callable_counts: np.ndarray
    p_alt_derived: np.ndarray | None
    metadata: dict
    identity: dict


_MISSING = object()


def eligibility_identity(metadata: dict) -> dict:
    """Return the content-addressed identity of a loaded eligibility mask.

    Callers (the target and candidate-row reports) record this dict verbatim
    so that a downstream consumer can require an exact match against it,
    binding the artifact by what it contains rather than by where it was
    resolved from. Raises ``ValueError`` if any of the required fields is
    absent from ``metadata`` (or, for the array digests, from its nested
    ``array_sha256``).
    """
    array_hashes = metadata.get("array_sha256")
    if not isinstance(array_hashes, dict):
        array_hashes = {}
    fields = {
        "vcf_sha256": metadata.get("vcf_sha256", _MISSING),
        "heterozygous": metadata.get("heterozygous", _MISSING),
        "min_callable": metadata.get("min_callable", _MISSING),
        "store_content_sha256": metadata.get("store_content_sha256", _MISSING),
        "row_indices_sha256": array_hashes.get("row_indices", _MISSING),
        "snp_row_indices_sha256": array_hashes.get("snp_row_indices", _MISSING),
        "p_alt_derived_sha256": array_hashes.get("p_alt_derived", _MISSING),
    }
    missing = [key for key, value in fields.items() if value is _MISSING or value is None]
    if missing:
        raise ValueError(
            "eligibility metadata is missing required identity field(s): "
            + ", ".join(missing)
        )
    return fields


def _chromosome_map(store: object) -> dict[str, tuple[int, int]]:
    chromosomes = getattr(store, "chromosomes", None)
    if chromosomes is None:
        chromosomes = (getattr(store, "metadata", {}) or {}).get("chromosomes")
    if not chromosomes:
        raise ValueError("store has no chromosome coordinate metadata")
    return {
        str(entry["chrom"]): (int(entry["offset"]), int(entry["length"]))
        for entry in chromosomes
    }


def scan_vcf(
    vcf: Path,
    store: object,
    *,
    min_callable: int = DEFAULT_MIN_CALLABLE,
    heterozygous: str = "error",
    ancestral_counts: np.ndarray | None = None,
    present_draw_count: np.ndarray | None = None,
) -> EligibilityResult:
    """Return store rows passing the shared VCF genotype/callability rules.

    Records absent from the interval-store catalog are irrelevant to matching
    and are not genotype-decoded.  A malformed potential candidate is excluded
    with a reason count instead of aborting the complete scan; downstream code
    then never requests that record from the stricter SFS reader.
    """
    if min_callable < 1:
        raise ValueError("min_callable must be positive")
    if heterozygous not in ("error", "missing"):
        raise ValueError("heterozygous must be 'error' or 'missing'")

    catalog = np.asarray(store.positions, dtype=np.float64)
    if catalog.ndim != 1 or np.any(catalog[1:] <= catalog[:-1]):
        raise ValueError("store positions must be a strictly increasing 1-D array")
    store_eligible = np.asarray(store.eligible)
    if store_eligible.dtype != bool or store_eligible.shape != catalog.shape:
        raise ValueError("store eligible mask must be boolean and align with positions")
    chromosomes = _chromosome_map(store)
    if (ancestral_counts is None) != (present_draw_count is None):
        raise ValueError("ancestral_counts and present_draw_count must be supplied together")
    if ancestral_counts is not None:
        ancestral_counts = np.asarray(ancestral_counts)
        present_draw_count = np.asarray(present_draw_count)
        if ancestral_counts.shape != (catalog.size, 4):
            raise ValueError("ancestral_counts must have shape (store rows, 4)")
        if present_draw_count.shape != (catalog.size,):
            raise ValueError("present_draw_count must align with store rows")

    rows: list[int] = []
    alt_counts: list[int] = []
    callable_counts: list[int] = []
    snp_rows: list[int] = []
    q_values: list[float] = []
    seen: set[int] = set()
    reasons: dict[str, int] = {}
    snp_reasons: dict[str, int] = {}
    records = relevant = 0
    genotype_cache: dict[str, int | None | str] = {}
    handle, hashing, buffered = _open_vcf(Path(vcf))
    try:
        for line_number, raw in enumerate(handle, 1):
            if raw.startswith("#"):
                continue
            records += 1
            # Only CHROM and POS are parsed for a record that turns out not to
            # be in the catalog. Splitting every sample column of every record
            # dominates scan time on a sample-rich VCF; a catalog lookup on two
            # fields is cheap and rejects most records before that cost.
            try:
                chrom, position_text, rest = raw.split("\t", 2)
            except ValueError:
                raise ValueError(f"{vcf}:{line_number}: malformed VCF record") from None
            try:
                position = int(position_text)
            except ValueError as error:
                raise ValueError(f"{vcf}:{line_number}: invalid POS") from error
            chrom_info = chromosomes.get(chrom)
            if chrom_info is None:
                continue
            offset, length = chrom_info
            if position < 1 or position > length:
                continue
            global_position = float(offset + position)
            row = int(np.searchsorted(catalog, global_position))
            if row >= catalog.size or catalog[row] != global_position:
                continue
            relevant += 1
            if row in seen:
                raise ValueError(f"duplicate VCF record at {chrom}:{position}")
            seen.add(row)
            if not store_eligible[row]:
                reasons["store_ineligible"] = reasons.get("store_ineligible", 0) + 1
                continue
            # ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, then sample columns.
            rest_fields = rest.rstrip("\n").split("\t")
            if len(rest_fields) < 8:
                reasons["missing_samples_or_format"] = (
                    reasons.get("missing_samples_or_format", 0) + 1
                )
                continue
            ref, alt, formats = rest_fields[1], rest_fields[2], rest_fields[6]
            if "," in alt:
                reasons["multiallelic_record"] = (
                    reasons.get("multiallelic_record", 0) + 1
                )
                continue
            format_fields = formats.split(":")
            if "GT" not in format_fields:
                reasons["missing_gt"] = reasons.get("missing_gt", 0) + 1
                continue
            gt_index = format_fields.index("GT")
            alt_count = callable_count = 0
            failure: str | None = None
            for sample in rest_fields[7:]:
                parts = sample.split(":")
                gt = parts[gt_index] if gt_index < len(parts) else "."
                allele = genotype_cache.get(gt)
                if gt not in genotype_cache:
                    allele = decode_inbred_genotype(gt, heterozygous)
                    genotype_cache[gt] = allele
                if allele is None:
                    continue
                if isinstance(allele, str):
                    failure = allele
                    break
                alt_count += allele
                callable_count += 1
            if failure is not None:
                reasons[failure] = reasons.get(failure, 0) + 1
                continue
            if callable_count < min_callable:
                reasons["insufficient_callability"] = (
                    reasons.get("insufficient_callability", 0) + 1
                )
                continue
            rows.append(row)
            alt_counts.append(alt_count)
            callable_counts.append(callable_count)
            if ancestral_counts is not None:
                # These two reasons remove a row only from the SNP-orientable
                # subset; the row stays callable and remains in row_indices.npy,
                # so they are counted separately from `reasons` above.
                if ref not in "ACGT" or alt not in "ACGT":
                    snp_reasons["snp_non_acgt_alleles"] = (
                        snp_reasons.get("snp_non_acgt_alleles", 0) + 1
                    )
                    continue
                counts = ancestral_counts[row]
                present = int(present_draw_count[row])
                if int(counts.sum()) > present:
                    raise ValueError(
                        f"ancestral table is inconsistent at {chrom}:{position}: "
                        f"{int(counts.sum())} calls across {present} present draws"
                    )
                ref_calls = int(counts["ACGT".index(ref)])
                alt_calls = int(counts["ACGT".index(alt)])
                oriented = ref_calls + alt_calls
                if oriented == 0:
                    snp_reasons["snp_no_usable_orientation"] = (
                        snp_reasons.get("snp_no_usable_orientation", 0) + 1
                    )
                    continue
                snp_rows.append(row)
                q_values.append(ref_calls / oriented)
        drain(buffered)
    finally:
        handle.close()
    digest = hashing.digest.hexdigest()

    order = np.argsort(rows, kind="stable")
    row_array = np.asarray(rows, dtype=np.int64)[order]
    alt_array = np.asarray(alt_counts, dtype=np.uint32)[order]
    callable_array = np.asarray(callable_counts, dtype=np.uint32)[order]
    snp_order = np.argsort(snp_rows, kind="stable")
    snp_array = np.asarray(snp_rows, dtype=np.int64)[snp_order]
    q_array = np.asarray(q_values, dtype=np.float64)[snp_order]
    report = {
        "schema_version": SCHEMA_VERSION,
        "vcf": str(Path(vcf).resolve()),
        "vcf_sha256": digest,
        "min_callable": int(min_callable),
        "heterozygous": heterozygous,
        "polarity_filtering": "none",
        "snp_orientation": "full posterior q; never thresholded",
        "applicable_set_roles": ["A", "B"],
        "applicable_variant_types": ["TE", "SNP"],
        "vcf_records": records,
        "store_catalog_records": relevant,
        "eligible_rows": int(row_array.size),
        "snp_orientable_rows": int(snp_array.size),
        # `excluded_by_reason` accounts only for rows dropped from the
        # callable mask: store_catalog_records - sum(excluded_by_reason) ==
        # eligible_rows. `snp_excluded_by_reason` accounts for rows that stay
        # callable (and remain in row_indices.npy) but are dropped from the
        # SNP-orientable subset: eligible_rows - sum(snp_excluded_by_reason)
        # == snp_orientable_rows.
        "excluded_by_reason": reasons,
        "snp_excluded_by_reason": snp_reasons,
        "store_content_sha256": (getattr(store, "metadata", {}) or {}).get(
            "content_sha256"
        ),
        "store_catalog_sha256": (getattr(store, "metadata", {}) or {}).get(
            "catalog_sha256"
        ),
    }
    return EligibilityResult(
        row_array, alt_array, callable_array, report, snp_array, q_array
    )


def load_ancestral_table(
    table_dir: Path, store: object
) -> tuple[np.ndarray, np.ndarray, dict]:
    """Authenticate and load the ARG ancestral-state counts for SNP q."""
    metadata = json.loads((table_dir / "metadata.json").read_text(encoding="utf-8"))
    if metadata.get("schema_version") != "ancestral-state-counts-v1":
        raise SystemExit(f"{table_dir}: unsupported ancestral-table schema")
    if not metadata.get("complete"):
        raise SystemExit(f"{table_dir}: ancestral table is incomplete")
    if metadata.get("bases") != ["A", "C", "G", "T"]:
        raise SystemExit(f"{table_dir}: ancestral columns are not ordered A,C,G,T")
    expected = (getattr(store, "metadata", {}) or {}).get("content_sha256")
    if not expected or metadata.get("store_content_sha256") != expected:
        raise SystemExit(f"{table_dir}: ancestral table does not authenticate against this store")
    n_rows = int(np.asarray(store.positions).size)
    if int(metadata.get("store_rows", -1)) != n_rows:
        raise SystemExit(f"{table_dir}: ancestral-table row count disagrees with the store")
    counts = np.load(table_dir / "ancestral_counts.npy", mmap_mode="r", allow_pickle=False)
    present = np.load(
        table_dir / "present_draw_count.npy", mmap_mode="r", allow_pickle=False
    )
    if counts.shape != (n_rows, 4) or counts.dtype.kind != "u":
        raise SystemExit(f"{table_dir}: ancestral_counts.npy has invalid shape or dtype")
    if present.shape != (n_rows,) or present.dtype.kind != "u":
        raise SystemExit(f"{table_dir}: present_draw_count.npy has invalid shape or dtype")
    if np.any(counts.sum(axis=1, dtype=np.uint64) > present):
        raise SystemExit(f"{table_dir}: ancestral counts exceed present-draw counts")
    return counts, present, metadata


def publish(output: Path, result: EligibilityResult, metadata: dict) -> None:
    """Atomically publish an authenticated eligibility-mask directory."""
    arrays = (result.rows, result.alt_counts, result.callable_counts)
    if any(np.asarray(values).ndim != 1 for values in arrays):
        raise ValueError("eligibility arrays must be one-dimensional")
    if len({np.asarray(values).size for values in arrays}) != 1:
        raise ValueError("eligibility arrays must have equal length")
    if result.rows.size != int(result.report.get("eligible_rows", -1)):
        raise ValueError("eligible row count disagrees with report")
    snp_rows = (
        np.empty(0, dtype=np.int64)
        if result.snp_rows is None else np.asarray(result.snp_rows)
    )
    q_values = (
        np.empty(0, dtype=np.float64)
        if result.p_alt_derived is None else np.asarray(result.p_alt_derived)
    )
    if snp_rows.ndim != 1 or q_values.shape != snp_rows.shape:
        raise ValueError("SNP rows and q values must be aligned 1-D arrays")
    if snp_rows.size != int(result.report.get("snp_orientable_rows", snp_rows.size)):
        raise ValueError("SNP-orientable row count disagrees with report")
    if output.exists():
        raise SystemExit(f"output already exists: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    staging = output.with_name(f".{output.name}.staging.{os.getpid()}")
    if staging.exists():
        shutil.rmtree(staging)
    staging.mkdir()
    try:
        np.save(staging / "row_indices.npy", result.rows, allow_pickle=False)
        np.save(staging / "alt_counts.npy", result.alt_counts, allow_pickle=False)
        np.save(staging / "callable_counts.npy", result.callable_counts, allow_pickle=False)
        np.save(staging / "snp_row_indices.npy", snp_rows, allow_pickle=False)
        np.save(staging / "p_alt_derived.npy", q_values, allow_pickle=False)
        published_report = {
            **result.report,
            "eligible_rows": int(result.rows.size),
            "snp_orientable_rows": int(snp_rows.size),
            "array_sha256": {
                "row_indices": _sha256_array(result.rows),
                "alt_counts": _sha256_array(result.alt_counts),
                "callable_counts": _sha256_array(result.callable_counts),
                "snp_row_indices": _sha256_array(snp_rows),
                "p_alt_derived": _sha256_array(q_values),
            },
        }
        (staging / "metadata.json").write_text(
            json.dumps(
                {**metadata, **published_report, "complete": True},
                indent=2,
                sort_keys=True,
            ) + "\n",
            encoding="utf-8",
        )
        os.replace(staging, output)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise


def load_eligibility(
    mask_dir: Path,
    store: object,
    *,
    variant_type: str,
    expected_min_callable: int | None = None,
) -> LoadedEligibility:
    """Load authenticated TE-callable or SNP-callable-and-orientable rows."""
    if variant_type not in ("TE", "SNP"):
        raise ValueError("variant_type must be TE or SNP")
    metadata = json.loads((mask_dir / "metadata.json").read_text(encoding="utf-8"))
    if metadata.get("schema_version") != SCHEMA_VERSION or not metadata.get("complete"):
        raise SystemExit(f"{mask_dir}: incomplete or unsupported VCF eligibility mask")
    store_metadata = getattr(store, "metadata", {}) or {}
    for key in ("store_content_sha256", "store_catalog_sha256"):
        expected = store_metadata.get(key.removeprefix("store_"))
        observed = metadata.get(key)
        if not expected or not observed or expected != observed:
            raise SystemExit(f"{mask_dir}: {key} does not authenticate against this store")
    if expected_min_callable is not None and metadata.get("min_callable") != expected_min_callable:
        raise SystemExit(
            f"{mask_dir}: min_callable={metadata.get('min_callable')}, "
            f"expected {expected_min_callable}"
        )
    rows = np.load(mask_dir / "row_indices.npy", allow_pickle=False)
    alt_counts = np.load(mask_dir / "alt_counts.npy", allow_pickle=False)
    callable_counts = np.load(mask_dir / "callable_counts.npy", allow_pickle=False)
    expected_hashes = metadata.get("array_sha256")
    if not isinstance(expected_hashes, dict):
        raise SystemExit(f"{mask_dir}: eligibility artifact lacks array digests")
    for name, values in (
        ("row_indices", rows),
        ("alt_counts", alt_counts),
        ("callable_counts", callable_counts),
    ):
        if expected_hashes.get(name) != _sha256_array(values):
            raise SystemExit(f"{mask_dir}: {name}.npy fails its content digest")
    if rows.ndim != 1 or not np.issubdtype(rows.dtype, np.integer):
        raise SystemExit(f"{mask_dir}: row_indices.npy must be a 1-D integer array")
    rows = rows.astype(np.int64, copy=False)
    n_rows = int(np.asarray(store.positions).size)
    if rows.size and (rows[0] < 0 or rows[-1] >= n_rows or np.any(rows[1:] <= rows[:-1])):
        raise SystemExit(f"{mask_dir}: row indices must be sorted, unique, and within the store")
    if rows.size != int(metadata.get("eligible_rows", -1)):
        raise SystemExit(f"{mask_dir}: row count disagrees with metadata")
    for name, values in (("alt_counts", alt_counts), ("callable_counts", callable_counts)):
        if (
            values.ndim != 1
            or values.size != rows.size
            or not np.issubdtype(values.dtype, np.unsignedinteger)
        ):
            raise SystemExit(
                f"{mask_dir}: {name}.npy must be an aligned unsigned 1-D array"
            )
    minimum = int(metadata.get("min_callable", -1))
    if np.any(callable_counts < minimum) or np.any(alt_counts > callable_counts):
        raise SystemExit(f"{mask_dir}: stored allele/callability counts are inconsistent")
    identity = eligibility_identity(metadata)
    if variant_type == "TE":
        return LoadedEligibility(rows, alt_counts, callable_counts, None, metadata, identity)

    snp_rows = np.load(mask_dir / "snp_row_indices.npy", allow_pickle=False)
    q_values = np.load(mask_dir / "p_alt_derived.npy", allow_pickle=False)
    for name, values in (("snp_row_indices", snp_rows), ("p_alt_derived", q_values)):
        if expected_hashes.get(name) != _sha256_array(values):
            raise SystemExit(f"{mask_dir}: {name}.npy fails its content digest")
    if (
        snp_rows.ndim != 1
        or not np.issubdtype(snp_rows.dtype, np.integer)
        or q_values.shape != snp_rows.shape
        or q_values.dtype.kind != "f"
        or np.any(~np.isfinite(q_values))
        or np.any((q_values < 0.0) | (q_values > 1.0))
    ):
        raise SystemExit(f"{mask_dir}: SNP row/q arrays are invalid")
    snp_rows = snp_rows.astype(np.int64, copy=False)
    if snp_rows.size > 1 and np.any(snp_rows[1:] <= snp_rows[:-1]):
        raise SystemExit(f"{mask_dir}: SNP rows must be sorted and unique")
    indices = np.searchsorted(rows, snp_rows)
    if (
        snp_rows.size != int(metadata.get("snp_orientable_rows", -1))
        or np.any(indices >= rows.size)
        or np.any(rows[indices] != snp_rows)
    ):
        raise SystemExit(f"{mask_dir}: SNP rows are not a subset of callable rows")
    return LoadedEligibility(
        snp_rows,
        alt_counts[indices],
        callable_counts[indices],
        q_values,
        metadata,
        identity,
    )


def load_eligible_rows(
    mask_dir: Path,
    store: object,
    *,
    variant_type: str,
    expected_min_callable: int = DEFAULT_MIN_CALLABLE,
) -> np.ndarray:
    """Compatibility convenience returning only rows for an A/B variant type."""
    return load_eligibility(
        mask_dir,
        store,
        variant_type=variant_type,
        expected_min_callable=expected_min_callable,
    ).rows


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", type=Path, required=True)
    parser.add_argument("--store", type=Path, required=True)
    parser.add_argument("--ancestral-table", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--min-callable", type=int, default=DEFAULT_MIN_CALLABLE)
    parser.add_argument("--heterozygous", choices=("error", "missing"), default="error")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    store = open_snp_age_store(args.store)
    ancestral_counts, present_draw_count, ancestral_metadata = load_ancestral_table(
        args.ancestral_table, store
    )
    result = scan_vcf(
        args.vcf,
        store,
        min_callable=args.min_callable,
        heterozygous=args.heterozygous,
        ancestral_counts=ancestral_counts,
        present_draw_count=present_draw_count,
    )
    publish(args.output, result, {
        "store": str(args.store.resolve()),
        "store_schema": store_schema(store),
        "ancestral_table": str(args.ancestral_table.resolve()),
        "ancestral_table_store_content_sha256": ancestral_metadata.get(
            "store_content_sha256"
        ),
        "software": software_provenance(),
    })
    print(f"VCF records       {result.report['vcf_records']:,}")
    print(f"store records     {result.report['store_catalog_records']:,}")
    print(f"eligible (n >= {args.min_callable}) {result.rows.size:,}")
    print(f"SNP orientable     {result.report['snp_orientable_rows']:,}")
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
