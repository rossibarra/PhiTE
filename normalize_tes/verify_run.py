"""Verify that a published Phi-SFS run meets the production acceptance gates."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import tempfile
from pathlib import Path
from typing import Any

import numpy as np

from .bootstrap_target_matcher import _candidate_array_digest
from .phi_sfs import SCHEMA_VERSION as PHI_SCHEMA_VERSION
from .phi_sfs import SYMMETRIC_NULL_DESIGN
from .release_provenance import PROJECT_VERSION, software_provenance
from .sample_age_matched_controls import _load_target, _sha256_arrays


MATCH_SCHEMA_VERSION = "bootstrap-target-matches-v1"
BLOCK_SIZE = 100


def _json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError as error:
        raise ValueError(f"required file is missing: {path}") from error
    if not isinstance(value, dict):
        raise ValueError(f"expected a JSON object: {path}")
    return value


def _integer_array(path: Path, label: str) -> np.ndarray:
    values = np.load(path, allow_pickle=False)
    if not np.issubdtype(values.dtype, np.integer):
        raise ValueError(f"{label} must have an integer dtype, not {values.dtype}")
    return values.astype(np.int64, copy=False)


def _average_ranks(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(values.size, dtype=np.float64)
    ranks[order] = np.arange(1, values.size + 1, dtype=np.float64)
    _, inverse = np.unique(values, return_inverse=True)
    sums = np.bincount(inverse, weights=ranks)
    counts = np.bincount(inverse)
    return (sums / counts)[inverse]


def _spearman(x: np.ndarray, y: np.ndarray) -> float:
    ranked_x = _average_ranks(x)
    ranked_y = _average_ranks(y)
    if np.ptp(ranked_x) == 0 or np.ptp(ranked_y) == 0:
        return 0.0
    return float(np.corrcoef(ranked_x, ranked_y)[0, 1])


def _smd(early: np.ndarray, late: np.ndarray) -> float:
    pooled = math.sqrt((float(early.var(ddof=1)) + float(late.var(ddof=1))) / 2)
    difference = float(late.mean() - early.mean())
    if pooled == 0:
        return 0.0 if difference == 0 else math.copysign(math.inf, difference)
    return difference / pooled


def _normalize(rows: np.ndarray) -> np.ndarray:
    totals = rows.sum(axis=-1, keepdims=True)
    if np.any(~np.isfinite(totals)) or np.any(totals <= 0):
        raise ValueError("Phi-SFS control spectra must have positive finite mass")
    return rows / totals


def _phi_to_pooled(raw: np.ndarray) -> np.ndarray:
    normalized = _normalize(np.asarray(raw, dtype=np.float64))
    pooled = _normalize(np.asarray(raw, dtype=np.float64).sum(axis=0, keepdims=True))[0]
    pooled_cdf = np.cumsum(pooled)
    cdf = np.cumsum(normalized, axis=1)
    # Production retains projected counts 1..19, spaced 1/20 apart.
    return np.sum(np.abs(cdf[:, :-1] - pooled_cdf[:-1]) / 20.0, axis=1)


class Checks:
    def __init__(self) -> None:
        self.criteria: list[dict[str, Any]] = []

    def add(self, criterion: str, passed: bool, value: Any, expected: str) -> None:
        if isinstance(value, np.generic):
            value = value.item()
        self.criteria.append({
            "criterion": criterion,
            "value": value,
            "expected": expected,
            "pass": bool(passed),
        })

    @property
    def passed(self) -> bool:
        return all(row["pass"] for row in self.criteria)


def verify(
    *,
    candidate_rows: Path,
    target: Path,
    matches: Path,
    phi: Path,
    a_type: str,
    expected_version: str,
    expected_commit: str,
    expected_replicates: int = 500,
    min_qc_passes: int = 451,
    min_resolved_fraction: float = 0.70,
) -> dict[str, Any]:
    """Return a complete production-verification report without writing files."""
    checks = Checks()
    candidate_report_path = candidate_rows.with_suffix(candidate_rows.suffix + ".json")
    candidate_meta = _json(candidate_report_path)
    target_rows, target_cdf, age_bins, threshold, target_meta = _load_target(target)
    match_meta = _json(matches / "metadata.json")
    phi_meta = _json(phi / "metadata.json")

    artifacts = {
        "candidate rows": candidate_meta,
        "target": target_meta,
        "matches": match_meta,
        "Phi-SFS": phi_meta,
    }
    for label, metadata in artifacts.items():
        software = metadata.get("software")
        software = software if isinstance(software, dict) else {}
        checks.add(
            f"{label}: software release",
            software.get("version") == expected_version,
            software.get("version"),
            f"== {expected_version}",
        )
        checks.add(
            f"{label}: Git commit",
            software.get("git_commit") == expected_commit,
            software.get("git_commit"),
            f"== {expected_commit}",
        )
        checks.add(
            f"{label}: clean checkout",
            software.get("git_dirty") is False,
            software.get("git_dirty"),
            "is false",
        )

    candidates = np.load(candidate_rows, allow_pickle=False)
    checks.add(
        "candidate rows: one-dimensional integer array",
        candidates.ndim == 1 and np.issubdtype(candidates.dtype, np.integer),
        f"shape={candidates.shape}, dtype={candidates.dtype}",
        "1-D integer array",
    )
    candidate_ordered = (
        candidates.ndim == 1
        and candidates.size > 0
        and np.all(np.diff(candidates.astype(np.int64, copy=False)) > 0)
    )
    checks.add(
        "candidate rows: nonempty, sorted, and unique",
        candidate_ordered,
        int(candidates.size),
        "strictly increasing with at least one row",
    )
    checks.add(
        "candidate rows: recorded count",
        candidate_meta.get("candidate_rows") == int(candidates.size),
        candidate_meta.get("candidate_rows"),
        f"== {int(candidates.size)}",
    )
    candidate_content_digest = _candidate_array_digest(candidates)
    checks.add(
        "candidate rows: content digest",
        candidate_meta.get("candidate_rows_sha256") == candidate_content_digest,
        candidate_meta.get("candidate_rows_sha256"),
        f"== {candidate_content_digest}",
    )
    match_candidate_digest = _sha256_arrays(np.asarray(candidates))
    checks.add(
        "candidate rows: matcher digest",
        match_meta.get("candidate_rows_digest") == match_candidate_digest,
        match_meta.get("candidate_rows_digest"),
        f"== {match_candidate_digest}",
    )
    resolution_entries = [
        *candidate_meta.get("inclusion_lists", []),
        *candidate_meta.get("exclusion_lists", []),
    ]
    fractions = [entry.get("resolved_fraction") for entry in resolution_entries]
    fractions_valid = bool(fractions) and all(
        isinstance(value, (int, float)) and value >= min_resolved_fraction
        for value in fractions
    )
    checks.add(
        "candidate rows: position-list resolution",
        fractions_valid,
        (
            min(fractions)
            if fractions and all(isinstance(x, (int, float)) for x in fractions)
            else None
        ),
        f">= {min_resolved_fraction}",
    )
    checks.add(
        "candidate rows: restricted to VCF eligibility",
        candidate_meta.get("universe") == "restricted"
        and isinstance(candidate_meta.get("vcf_eligibility_identity"), dict),
        candidate_meta.get("universe"),
        "restricted with an eligibility identity",
    )

    store_pairs = (
        ("content", candidate_meta.get("store_content_sha256"),
         target_meta.get("source_store_content_sha256"),
         match_meta.get("source_store_content_sha256")),
        ("catalog", candidate_meta.get("store_catalog_sha256"),
         target_meta.get("source_catalog_sha256"),
         match_meta.get("source_catalog_sha256")),
    )
    for label, candidate_value, target_value, match_value in store_pairs:
        values = [candidate_value, target_value, match_value]
        checks.add(
            f"provenance: store {label} identity",
            candidate_value is not None and len(set(values)) == 1,
            values,
            "candidate, target, and matches agree",
        )

    target_eligibility = (target_meta.get("vcf_eligibility") or {}).get("identity")
    eligibility_values = [
        candidate_meta.get("vcf_eligibility_identity"),
        target_eligibility,
        match_meta.get("vcf_eligibility_identity"),
        phi_meta.get("vcf_eligibility_identity"),
    ]
    checks.add(
        "provenance: VCF eligibility identity",
        isinstance(target_eligibility, dict)
        and all(value == target_eligibility for value in eligibility_values),
        (
            "agree"
            if len({json.dumps(v, sort_keys=True) for v in eligibility_values}) == 1
            else "differ"
        ),
        "candidate, target, matches, and Phi-SFS agree",
    )

    target_digest = _sha256_arrays(
        target_rows, target_cdf, age_bins, np.asarray([threshold], dtype=np.float64)
    )
    checks.add(
        "provenance: target digest in matches",
        match_meta.get("target_digest") == target_digest,
        match_meta.get("target_digest"),
        f"== {target_digest}",
    )
    checks.add(
        "provenance: target digest in Phi-SFS",
        phi_meta.get("target_digest") == target_digest,
        phi_meta.get("target_digest"),
        f"== {target_digest}",
    )
    checks.add(
        "target: focal type",
        target_meta.get("a_type") == a_type and match_meta.get("a_type") == a_type
        and phi_meta.get("a_type") == a_type,
        [target_meta.get("a_type"), match_meta.get("a_type"), phi_meta.get("a_type")],
        f"all == {a_type}",
    )
    checks.add(
        "target: no TE polarity mask",
        target_meta.get("te_polarity") is None,
        target_meta.get("te_polarity"),
        "is null",
    )

    match_rows = _integer_array(matches / "row_indices.npy", "matched row indices")
    replicate_ids = _integer_array(matches / "replicate_id.npy", "replicate IDs")
    qc = np.load(matches / "qc_pass.npy", allow_pickle=False)
    if match_rows.ndim != 2:
        raise ValueError("matched row indices must be a two-dimensional array")
    if replicate_ids.shape != (match_rows.shape[0],):
        raise ValueError("replicate IDs do not align with matched sets")
    if qc.dtype.kind != "b" or qc.shape != replicate_ids.shape:
        raise ValueError("qc_pass.npy must be a boolean array aligned with matched sets")
    checks.add(
        "matches: complete supported bundle",
        match_meta.get("complete") is True
        and match_meta.get("schema_version") == MATCH_SCHEMA_VERSION,
        [match_meta.get("complete"), match_meta.get("schema_version")],
        f"complete {MATCH_SCHEMA_VERSION}",
    )
    checks.add(
        "matches: published set count",
        match_rows.shape[0] == expected_replicates
        and match_meta.get("replicates") == expected_replicates,
        [match_rows.shape[0], match_meta.get("replicates")],
        f"both == {expected_replicates}",
    )
    checks.add(
        "matches: unique ordered replicate IDs",
        np.array_equal(replicate_ids, np.arange(expected_replicates)),
        "0..N-1" if np.array_equal(replicate_ids, np.arange(replicate_ids.size)) else "invalid",
        f"0..{expected_replicates - 1}",
    )
    config = match_meta.get("config") if isinstance(match_meta.get("config"), dict) else {}
    checks.add(
        "matches: disjoint mode",
        config.get("disjoint_replicates") is True,
        config.get("disjoint_replicates"),
        "is true",
    )
    duplicate_sets = sum(np.unique(row).size != row.size for row in match_rows)
    flattened_rows = np.sort(match_rows, axis=None)
    globally_disjoint = flattened_rows.size > 0 and not np.any(
        np.diff(flattened_rows) == 0
    )
    reuse_row_indices = _integer_array(
        matches / "reuse_row_indices.npy", "reuse row indices"
    )
    reuse_counts = _integer_array(matches / "reuse_counts.npy", "reuse counts")
    maximum_reuse = int(reuse_counts.max(initial=0))
    checks.add(
        "matches: maximum control reuse",
        maximum_reuse == 1 and match_meta.get("maximum_control_reuse") == 1,
        [maximum_reuse, match_meta.get("maximum_control_reuse")],
        "both == 1",
    )
    checks.add(
        "matches: duplicate controls within sets",
        duplicate_sets == 0,
        int(duplicate_sets),
        "== 0",
    )
    target_overlap = int(np.intersect1d(target_rows, flattened_rows).size)
    checks.add(
        "matches: controls exclude the focal set",
        target_overlap == 0,
        target_overlap,
        "== 0 overlapping rows",
    )
    if candidate_ordered:
        candidate_locations = np.searchsorted(candidates, flattened_rows)
        candidates_cover_matches = bool(
            np.all(candidate_locations < candidates.size)
            and np.array_equal(candidates[candidate_locations], flattened_rows)
        )
    else:
        candidates_cover_matches = False
    checks.add(
        "matches: every control belongs to the candidate universe",
        candidates_cover_matches,
        int(flattened_rows.size),
        "all matched rows present",
    )
    checks.add(
        "matches: controls are globally disjoint",
        globally_disjoint
        and np.array_equal(reuse_row_indices, flattened_rows)
        and reuse_counts.shape == reuse_row_indices.shape
        and np.all(reuse_counts == 1),
        int(flattened_rows.size),
        "all rows unique and exactly represented by reuse arrays",
    )
    qc_passes = int(qc.sum())
    checks.add(
        "matches: QC-passing sets",
        qc_passes >= min_qc_passes and match_meta.get("qc_passes") == qc_passes,
        [qc_passes, match_meta.get("qc_passes")],
        f"actual >= {min_qc_passes} and metadata agrees",
    )
    blocks: list[dict[str, Any]] = []
    for block in range(replicate_ids.size // BLOCK_SIZE):
        start = block * BLOCK_SIZE
        stop = start + BLOCK_SIZE
        blocks.append({
            "block": block,
            "first_replicate": start,
            "last_replicate": stop - 1,
            "qc_passes": int(qc[start:stop].sum()),
        })
    fewest_block = min(block["qc_passes"] for block in blocks) if blocks else 0
    qc_fall = blocks[0]["qc_passes"] - blocks[-1]["qc_passes"] if blocks else math.inf
    checks.add(
        "matches: fewest QC passes per complete block of 100",
        fewest_block >= 80,
        fewest_block,
        ">= 80",
    )
    checks.add(
        "matches: QC pass-rate fall from first to last block",
        qc_fall <= 10,
        qc_fall,
        "<= 10 percentage points",
    )

    set_size = int(target_rows.size)
    checks.add(
        "site count: target and every matched set",
        match_rows.shape[1] == set_size and match_meta.get("set_size") == set_size,
        [set_size, match_rows.shape[1], match_meta.get("set_size")],
        "all equal",
    )
    checks.add(
        "Phi-SFS: complete supported result",
        phi_meta.get("complete") is True
        and phi_meta.get("schema_version") == PHI_SCHEMA_VERSION,
        [phi_meta.get("complete"), phi_meta.get("schema_version")],
        f"complete {PHI_SCHEMA_VERSION}",
    )
    checks.add(
        "Phi-SFS: posterior-mixture null design",
        phi_meta.get("null_polarity_design") == SYMMETRIC_NULL_DESIGN,
        phi_meta.get("null_polarity_design"),
        f"== {SYMMETRIC_NULL_DESIGN}",
    )
    checks.add(
        "Phi-SFS: equal eligible site count",
        phi_meta.get("equal_eligible_site_count") == set_size,
        phi_meta.get("equal_eligible_site_count"),
        f"== {set_size}",
    )
    checks.add(
        "Phi-SFS: matched set count and disjointness",
        phi_meta.get("matched_sets_published") == expected_replicates
        and phi_meta.get("disjoint_replicates") is True
        and phi_meta.get("maximum_control_reuse") == 1,
        [phi_meta.get("matched_sets_published"), phi_meta.get("disjoint_replicates"),
         phi_meta.get("maximum_control_reuse")],
        f"{expected_replicates}, true, 1",
    )
    phi_ids = _integer_array(phi / "b_replicate_id.npy", "Phi-SFS replicate IDs")
    if phi_ids.ndim != 1:
        raise ValueError("b_replicate_id.npy must be one-dimensional")
    accepted_ids = {
        int(phi_meta.get("reference_replicate_id")),
        *(int(value) for value in phi_meta.get("selected_null_replicate_ids", [])),
    }
    id_to_index = {int(value): index for index, value in enumerate(replicate_ids)}
    ids_known = all(value in id_to_index for value in accepted_ids)
    accepted_pass_qc = ids_known and all(qc[id_to_index[value]] for value in accepted_ids)
    checks.add(
        "Phi-SFS: every analyzed set passed matching QC",
        accepted_pass_qc
        and set(phi_ids.tolist()) == accepted_ids
        and phi_ids.size == len(accepted_ids)
        and np.unique(phi_ids).size == phi_ids.size,
        len(accepted_ids),
        "all IDs are QC-passing and match b_replicate_id.npy",
    )
    null_count = len(accepted_ids) - 1
    checks.add(
        "Phi-SFS: accepted null count",
        null_count >= min_qc_passes - 1
        and phi_meta.get("accepted_null_replicates") == null_count,
        [null_count, phi_meta.get("accepted_null_replicates")],
        f"actual >= {min_qc_passes - 1} and metadata agrees",
    )

    raw = np.load(phi / "b_raw_sfs.npy", allow_pickle=False)
    if (
        raw.ndim != 2
        or raw.shape != (phi_ids.size, 19)
        or np.any(~np.isfinite(raw))
        or np.any(raw < 0)
    ):
        raise ValueError(
            "b_raw_sfs.npy must be a finite nonnegative (analyzed sets x 19) "
            "array aligned with b_replicate_id.npy"
        )
    order = np.argsort(phi_ids, kind="stable")
    ordered_ids = phi_ids[order]
    distances = _phi_to_pooled(raw[order])
    quarter = max(2, ordered_ids.size // 4)
    rho = _spearman(ordered_ids.astype(np.float64), distances)
    effect = _smd(distances[:quarter], distances[-quarter:])
    checks.add(
        "depletion: absolute Spearman rho",
        math.isfinite(rho) and abs(rho) < 0.1,
        round(rho, 6),
        "< 0.1",
    )
    checks.add(
        "depletion: absolute early-vs-late SMD",
        math.isfinite(effect) and abs(effect) < 0.2,
        round(effect, 6),
        "< 0.2",
    )

    return {
        "schema_version": "phite-run-verification-v1",
        "pass": checks.passed,
        "software": software_provenance(),
        "inputs": {
            "candidate_rows": str(candidate_rows.resolve()),
            "target": str(target.resolve()),
            "matches": str(matches.resolve()),
            "phi": str(phi.resolve()),
        },
        "expected": {
            "a_type": a_type,
            "version": expected_version,
            "commit": expected_commit,
            "replicates": expected_replicates,
            "min_qc_passes": min_qc_passes,
            "min_resolved_fraction": min_resolved_fraction,
        },
        "criteria": checks.criteria,
        "blocks": blocks,
        "diagnostics": {
            "set_size": set_size,
            "qc_passes": qc_passes,
            "analyzed_sets": int(phi_ids.size),
            "depletion_quarter_size": quarter,
            "depletion_rho": rho,
            "depletion_smd": effect,
        },
    }


def _publish(output: Path, report: dict[str, Any]) -> None:
    if output.exists():
        raise FileExistsError(f"output already exists: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{output.name}.tmp.", dir=output.parent))
    try:
        (staging / "report.json").write_text(
            json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        with (staging / "criteria.csv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=["criterion", "value", "expected", "pass"]
            )
            writer.writeheader()
            writer.writerows(report["criteria"])
        with (staging / "blocks.csv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=["block", "first_replicate", "last_replicate", "qc_passes"],
            )
            writer.writeheader()
            writer.writerows(report["blocks"])
        os.replace(staging, output)
    except BaseException:
        for path in staging.iterdir():
            path.unlink()
        staging.rmdir()
        raise


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    current_commit = software_provenance().get("git_commit")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate-rows", type=Path, required=True)
    parser.add_argument("--target", type=Path, required=True)
    parser.add_argument("--matches", type=Path, required=True)
    parser.add_argument("--phi", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("-A", "--a-type", choices=("TE", "SNP"), default="TE")
    parser.add_argument("--expected-version", default=PROJECT_VERSION)
    parser.add_argument(
        "--expected-commit", default=current_commit,
        help="commit expected in every artifact (default: current checkout HEAD)",
    )
    parser.add_argument("--expected-replicates", type=int, default=500)
    parser.add_argument("--min-qc-passes", type=int, default=451)
    parser.add_argument("--min-resolved-fraction", type=float, default=0.70)
    args = parser.parse_args(argv)
    if not args.expected_commit:
        parser.error("--expected-commit is required outside a Git checkout")
    if args.expected_replicates < 2:
        parser.error("--expected-replicates must be at least 2")
    if not 1 <= args.min_qc_passes <= args.expected_replicates:
        parser.error("--min-qc-passes must be between 1 and --expected-replicates")
    if not 0 <= args.min_resolved_fraction <= 1:
        parser.error("--min-resolved-fraction must be between 0 and 1")
    return args


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    report = verify(
        candidate_rows=args.candidate_rows,
        target=args.target,
        matches=args.matches,
        phi=args.phi,
        a_type=args.a_type,
        expected_version=args.expected_version,
        expected_commit=args.expected_commit,
        expected_replicates=args.expected_replicates,
        min_qc_passes=args.min_qc_passes,
        min_resolved_fraction=args.min_resolved_fraction,
    )
    _publish(args.output, report)
    for row in report["criteria"]:
        print(
            f"{str(row['pass']):5}  {row['criterion']}: "
            f"{row['value']} ({row['expected']})"
        )
    print(f"wrote {args.output / 'criteria.csv'}")
    print(f"wrote {args.output / 'report.json'}")
    return 0 if report["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
