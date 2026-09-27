#!/usr/bin/env python3
"""Formal contrasts between Phi-SFS result categories.

`normalize_tes.phi_sfs` calibrates one focal category (A versus its matched
neutral SNP controls) and publishes a Z-score and a null distribution of
standardized null distances (`null_z_scores.npy`). This module compares that
Z-score across two or more already-published categories -- for example TE
subfamilies, or a TE category against a SNP-versus-SNP negative control -- so
that a plotted difference in Z between categories is a tested claim rather
than a description. PHI_SFS_WASSERSTEIN_CODING_PLAN.md section 8 gives the
statistical design this module implements.

Each category's `null_z_scores.npy` entries are already standardized (mean
zero, sample SD one) null distances for *that* category alone; they carry no
shared replicate identity across categories; two categories are matched
independently, with independent bootstrap draws and independent matched
control sets, so index i of one category's null vector and index i of
another's describe unrelated random draws. There is therefore no natural
pairing between them. This module fixes one arbitrary, prespecified,
deterministic pairing per ordered category pair and repeats it under
independently drawn pairings as a sensitivity check (`--pairing-repeats`).

Because the R null replicates within one category are exchangeable (nothing
distinguishes replicate i from replicate j once QC has passed), any bijection
between the two categories' R null indices is an equally valid finite-sample
null realization; the choice of pairing is a matter of reproducibility, not
of correctness, which is what the sensitivity repeats are for.
"""

from __future__ import annotations

import argparse
import csv
import datetime
import hashlib
import json
import math
import os
import re
import shutil
import sys
import tempfile
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import NamedTuple, Sequence

import numpy as np

from .release_provenance import software_provenance


SCHEMA_VERSION = "phi-contrast-v1"
REQUIRED_RESULT_SCHEMA = "phi-sfs-wasserstein-v2"
DEFAULT_SEED = 1002
DEFAULT_PAIRING_REPEATS = 100
_LABEL_SAFE = re.compile(r"[^A-Za-z0-9_.-]")


# --------------------------------------------------------------- pure math


class ContrastResult(NamedTuple):
    """One pairing's contrast between two categories' Z-scores."""

    delta_obs: float
    null_delta: np.ndarray
    pairing: np.ndarray
    p_value: float
    exceedances: int


def derive_seed(base_seed: int, label1: str, label2: str, repeat: int) -> int:
    """Return a deterministic seed for one pairing repeat of an unordered pair.

    The seed depends on the pair of labels only through their sorted order, so
    it is identical regardless of which label is passed first -- the pairing
    itself (and therefore the null distribution it produces) does not depend
    on argument order. `contrast()` is still order-sensitive in the sense that
    swapping which array is "1" and which is "2" changes which side gets
    permuted; `paired_contrast()` below is the order-independent wrapper that
    most callers want.
    """
    if not isinstance(base_seed, (int, np.integer)):
        raise TypeError("base_seed must be an integer")
    if not isinstance(repeat, (int, np.integer)) or repeat < 0:
        raise ValueError("repeat must be a nonnegative integer")
    lo, hi = sorted((str(label1), str(label2)))
    digest = hashlib.sha256(
        f"{int(base_seed)}|{lo}|{hi}|{int(repeat)}".encode("utf-8")
    ).digest()
    return int.from_bytes(digest[:8], "big")


def contrast(
    z1: float,
    z2: float,
    null1: np.ndarray,
    null2: np.ndarray,
    rng: np.random.Generator,
) -> ContrastResult:
    """Contrast two categories' observed Z-scores against a random null pairing.

    ``delta_obs = z1 - z2``. With ``R = min(R1, R2)``, the null pairing draws
    ``pairing = rng.permutation(R1)[:R]`` from the caller's generator and sets
    ``null_delta[i] = null1[pairing[i]] - null2[i]`` -- a prespecified
    deterministic random pairing, since the two categories' null replicates
    carry no shared identity. When the second category has more replicates
    than the first, a random subset of R of them is drawn from the same
    generator (``rng.permutation(R2)[:R]``) and paired in that order. With
    equal R this reduces exactly to a permutation of the first side. The two-sided add-one
    P-value counts null contrasts with ``abs(null_delta) >= abs(delta_obs)``,
    so an exact tie counts as an exceedance, and applies the add-one correction
    to numerator and denominator so the minimum attainable P-value is
    ``1 / (R + 1)``, never zero.

    This function is order-sensitive: calling it with the two categories
    swapped does not, in general, simply negate the result, because the
    pairing formula treats its first and second arguments asymmetrically (only
    the first side's index is permuted). Use `paired_contrast` for a contrast
    between two labelled categories that is guaranteed order-independent up to
    an overall sign.
    """
    z1 = float(z1)
    z2 = float(z2)
    null1 = np.asarray(null1, dtype=np.float64)
    null2 = np.asarray(null2, dtype=np.float64)
    if not math.isfinite(z1) or not math.isfinite(z2):
        raise ValueError("observed Z-scores must be finite")
    if null1.ndim != 1 or null2.ndim != 1:
        raise ValueError("null Z-score vectors must be one-dimensional")
    if null1.size < 1 or null2.size < 1:
        raise ValueError("at least one null replicate is required")
    if not np.all(np.isfinite(null1)) or not np.all(np.isfinite(null2)):
        raise ValueError("null Z-score vectors must be finite")

    r = min(null1.size, null2.size)
    delta_obs = z1 - z2
    pairing = rng.permutation(null1.size)[:r]
    second = rng.permutation(null2.size)[:r] if null2.size > r else np.arange(r)
    null_delta = null1[pairing] - null2[second]
    exceedances = int(np.count_nonzero(np.abs(null_delta) >= abs(delta_obs)))
    p_value = (1.0 + exceedances) / (r + 1.0)
    return ContrastResult(delta_obs, null_delta, pairing, p_value, exceedances)


def paired_contrast(
    label1: str,
    z1: float,
    null1: np.ndarray,
    label2: str,
    z2: float,
    null2: np.ndarray,
    *,
    seed: int,
    repeat: int,
) -> ContrastResult:
    """Order-independent-seed wrapper around `contrast` for a labelled pair.

    The pairing is always drawn in label-sorted order (using `derive_seed`,
    which ignores argument order) and then sign-flipped, if necessary, to
    match the caller's requested (label1, label2) order. Negating `delta_obs`
    and every entry of `null_delta` together leaves every ``abs(...)`` value,
    and therefore the P-value, exactly unchanged; only the reported sign of
    the contrast changes with reporting order. So: swapping (label1, z1, null1)
    with (label2, z2, null2) flips the sign of `delta_obs` and of every
    `null_delta` entry, and leaves `p_value` and `exceedances` exactly the
    same -- this is a deliberate design choice (see module docstring), not a
    general property of `contrast`.
    """
    seed_value = derive_seed(seed, label1, label2, repeat)
    rng = np.random.default_rng(seed_value)
    lo, hi = sorted((str(label1), str(label2)))
    if (str(label1), str(label2)) == (lo, hi):
        result = contrast(z1, z2, null1, null2, rng)
        return result
    result = contrast(z2, z1, null2, null1, rng)
    return ContrastResult(
        delta_obs=-result.delta_obs,
        null_delta=-result.null_delta,
        pairing=result.pairing,
        p_value=result.p_value,
        exceedances=result.exceedances,
    )


def holm(p_values: Sequence[float]) -> np.ndarray:
    """Holm step-down multiplicity-adjusted P-values, in the input's order.

    Ascending-sorted raw P-values are multiplied by ``n, n-1, ..., 1``, forced
    to be nondecreasing by a running maximum (the standard Holm enforcement so
    that a less significant test is never assigned a smaller adjusted
    P-value than a more significant one), and capped at 1.
    """
    p = np.asarray(p_values, dtype=np.float64)
    if p.ndim != 1:
        raise ValueError("p_values must be one-dimensional")
    n = p.size
    if n == 0:
        return np.zeros(0, dtype=np.float64)
    if not np.all(np.isfinite(p)) or np.any(p < 0.0) or np.any(p > 1.0):
        raise ValueError("p_values must be finite and lie in [0, 1]")
    order = np.argsort(p, kind="stable")
    sorted_p = p[order]
    multipliers = (n - np.arange(n)).astype(np.float64)
    raw = sorted_p * multipliers
    enforced = np.minimum(np.maximum.accumulate(raw), 1.0)
    adjusted = np.empty(n, dtype=np.float64)
    adjusted[order] = enforced
    return adjusted


def benjamini_hochberg(p_values: Sequence[float]) -> np.ndarray:
    """Benjamini-Hochberg step-up FDR-adjusted P-values, in the input's order.

    Ascending-sorted raw P-values are multiplied by ``n / rank``, forced to be
    nonincreasing by a running minimum taken from the largest rank down (the
    standard BH enforcement), and capped at 1.
    """
    p = np.asarray(p_values, dtype=np.float64)
    if p.ndim != 1:
        raise ValueError("p_values must be one-dimensional")
    n = p.size
    if n == 0:
        return np.zeros(0, dtype=np.float64)
    if not np.all(np.isfinite(p)) or np.any(p < 0.0) or np.any(p > 1.0):
        raise ValueError("p_values must be finite and lie in [0, 1]")
    order = np.argsort(p, kind="stable")
    sorted_p = p[order]
    ranks = np.arange(1, n + 1, dtype=np.float64)
    raw = sorted_p * n / ranks
    enforced = np.minimum(np.minimum.accumulate(raw[::-1])[::-1], 1.0)
    adjusted = np.empty(n, dtype=np.float64)
    adjusted[order] = enforced
    return adjusted


# ------------------------------------------------------------------- I/O


@dataclass(frozen=True)
class CategoryResult:
    """One loaded phi-sfs-wasserstein-v2 result, as needed for a contrast."""

    label: str
    directory: Path
    a_type: str
    b_type: str
    accepted_null_replicates: int
    z_score: float
    null_z_scores: np.ndarray
    target_digest: object
    schema_version: str
    null_polarity_design: str


def _json(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def _read_summary_z_score(directory: Path) -> float | None:
    summary_path = directory / "summary.csv"
    if not summary_path.exists():
        return None
    with summary_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 1 or "z_score" not in rows[0]:
        raise ValueError(f"{summary_path} must contain exactly one row with a z_score column")
    try:
        return float(rows[0]["z_score"])
    except (TypeError, ValueError) as error:
        raise ValueError(f"{summary_path} has a non-numeric z_score") from error


def load_category(label: str, directory: Path) -> CategoryResult:
    """Load and validate one phi-sfs-wasserstein-v2 result directory."""
    metadata_path = directory / "metadata.json"
    if not metadata_path.exists():
        raise ValueError(f"{label} ({directory}): no metadata.json found")
    metadata = _json(metadata_path)

    schema = metadata.get("schema_version")
    if schema != REQUIRED_RESULT_SCHEMA:
        raise ValueError(
            f"{label} ({directory}): unsupported schema_version {schema!r}; "
            f"expected {REQUIRED_RESULT_SCHEMA!r}"
        )
    if metadata.get("complete") is not True:
        raise ValueError(f"{label} ({directory}): result is not marked complete")

    a_type = metadata.get("a_type")
    b_type = metadata.get("b_type")
    if a_type not in ("TE", "SNP"):
        raise ValueError(f"{label} ({directory}): metadata a_type must be 'TE' or 'SNP', got {a_type!r}")
    if b_type not in ("SNP",):
        raise ValueError(f"{label} ({directory}): metadata b_type must be 'SNP', got {b_type!r}")

    null_polarity_design = metadata.get("null_polarity_design")
    if not isinstance(null_polarity_design, str) or not null_polarity_design:
        raise ValueError(
            f"{label} ({directory}): metadata null_polarity_design must be a "
            "non-empty string"
        )

    r = metadata.get("accepted_null_replicates")
    if not isinstance(r, (int, np.integer)) or r < 1:
        raise ValueError(
            f"{label} ({directory}): metadata accepted_null_replicates must be a "
            f"positive integer, got {r!r}"
        )
    r = int(r)

    z_score = metadata.get("z_score")
    if not isinstance(z_score, (int, float)) or not math.isfinite(float(z_score)):
        raise ValueError(f"{label} ({directory}): metadata z_score must be a finite number")
    z_score = float(z_score)

    summary_z = _read_summary_z_score(directory)
    if summary_z is not None and not math.isclose(summary_z, z_score, rel_tol=0.0, abs_tol=1e-9):
        raise ValueError(
            f"{label} ({directory}): metadata z_score {z_score!r} disagrees with "
            f"summary.csv z_score {summary_z!r}"
        )

    null_path = directory / "null_z_scores.npy"
    if not null_path.exists():
        raise ValueError(f"{label} ({directory}): no null_z_scores.npy found")
    null_z = np.load(null_path, allow_pickle=False).astype(np.float64)
    if null_z.ndim != 1 or null_z.size != r:
        raise ValueError(
            f"{label} ({directory}): null_z_scores.npy has {null_z.size} entries, "
            f"but metadata declares accepted_null_replicates={r}"
        )
    if not np.all(np.isfinite(null_z)):
        raise ValueError(f"{label} ({directory}): null_z_scores.npy must be finite")

    return CategoryResult(
        label=label,
        directory=directory,
        a_type=a_type,
        b_type=b_type,
        accepted_null_replicates=r,
        z_score=z_score,
        null_z_scores=null_z,
        target_digest=metadata.get("target_digest"),
        schema_version=schema,
        null_polarity_design=null_polarity_design,
    )


def _sanitize_label(label: str) -> str:
    sanitized = _LABEL_SAFE.sub("_", label)
    return sanitized or "_"


def calculate(args: argparse.Namespace) -> None:
    if len(args.result) < 2:
        raise ValueError("at least two --result categories are required for a contrast")
    labels_seen: dict[str, Path] = {}
    entries: list[tuple[str, Path]] = []
    for raw in args.result:
        if "=" not in raw:
            raise ValueError(f"--result must be LABEL=DIR, got {raw!r}")
        label, _, directory = raw.partition("=")
        label = label.strip()
        if not label:
            raise ValueError(f"--result has an empty label: {raw!r}")
        if label in labels_seen:
            raise ValueError(
                f"duplicate --result label {label!r}: {labels_seen[label]} and {directory}"
            )
        path = Path(directory)
        labels_seen[label] = path
        entries.append((label, path))

    if args.pairing_repeats < 1:
        raise ValueError("--pairing-repeats must be at least 1")

    categories = [load_category(label, directory) for label, directory in entries]

    b_types = {category.b_type for category in categories}
    if len(b_types) != 1:
        detail = ", ".join(f"{c.label}={c.b_type}" for c in categories)
        raise ValueError(f"all --result categories must share one b_type; got {detail}")

    null_polarity_designs = {
        category.null_polarity_design for category in categories
    }
    if len(null_polarity_designs) != 1:
        detail = ", ".join(
            f"{category.label}={category.null_polarity_design}"
            for category in categories
        )
        raise ValueError(
            "all --result categories must share one null_polarity_design; "
            f"got {detail}"
        )

    # Categories may carry different R (every QC-passing null is used); each
    # pair is contrasted over min(R1, R2) randomly paired replicates.

    a_types = {category.a_type for category in categories}
    a_type_consistent = len(a_types) == 1
    if not a_type_consistent:
        detail = ", ".join(f"{c.label}={c.a_type}" for c in categories)
        print(
            f"warning: --result categories do not share one a_type; comparing across "
            f"variant types: {detail}",
            flush=True,
        )

    # Labels are already unique (checked above), so a collision here can only
    # come from two distinct labels sanitizing to the same filename fragment.
    sanitized_names: dict[str, str] = {}
    for category in categories:
        sanitized = _sanitize_label(category.label)
        if sanitized in sanitized_names.values():
            raise ValueError(
                f"label {category.label!r} sanitizes to a filename that collides with "
                "another label; rename one of the --result labels"
            )
        sanitized_names[category.label] = sanitized

    pairs = list(combinations(range(len(categories)), 2))
    if not pairs:
        raise ValueError("no category pairs to contrast")

    rows: list[dict[str, object]] = []
    primary_null_arrays: dict[str, np.ndarray] = {}
    pair_seed_records: list[dict[str, object]] = []
    primary_p_values: list[float] = []

    for i, j in pairs:
        c1, c2 = categories[i], categories[j]
        p_values: list[float] = []
        primary_delta_obs = None
        primary_null_delta = None
        for repeat in range(args.pairing_repeats):
            result = paired_contrast(
                c1.label, c1.z_score, c1.null_z_scores,
                c2.label, c2.z_score, c2.null_z_scores,
                seed=args.seed, repeat=repeat,
            )
            p_values.append(result.p_value)
            if repeat == 0:
                primary_delta_obs = result.delta_obs
                primary_null_delta = result.null_delta
        p_array = np.asarray(p_values, dtype=np.float64)
        primary_p = float(p_array[0])
        primary_p_values.append(primary_p)

        filename = f"null_contrasts_{sanitized_names[c1.label]}__{sanitized_names[c2.label]}.npy"
        if filename in primary_null_arrays:
            raise ValueError(f"sanitized pair filename collision: {filename}")
        primary_null_arrays[filename] = primary_null_delta

        rows.append({
            "label1": c1.label,
            "label2": c2.label,
            "a_type1": c1.a_type,
            "a_type2": c2.a_type,
            "b_type": c1.b_type,
            "r1": c1.accepted_null_replicates,
            "r2": c2.accepted_null_replicates,
            "r": min(c1.accepted_null_replicates, c2.accepted_null_replicates),
            "z1": c1.z_score,
            "z2": c2.z_score,
            "delta_obs": float(primary_delta_obs),
            "p": primary_p,
            "p_min": float(p_array.min()),
            "p_median": float(np.median(p_array)),
            "p_max": float(p_array.max()),
            "null_contrasts_file": filename,
        })
        pair_seed_records.append({
            "label1": c1.label,
            "label2": c2.label,
            "seed_repeat0": derive_seed(args.seed, c1.label, c2.label, 0),
        })

    holm_p = holm(primary_p_values)
    bh_p = benjamini_hochberg(primary_p_values)
    for row, holm_value, bh_value in zip(rows, holm_p, bh_p):
        row["holm_p"] = float(holm_value)
        row["bh_p"] = float(bh_value)

    output = args.output
    if output.exists():
        raise FileExistsError(f"output already exists: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{output.name}.tmp.", dir=output.parent))
    try:
        for filename, values in primary_null_arrays.items():
            np.save(staging / filename, values, allow_pickle=False)

        fieldnames = list(rows[0])
        with (staging / "contrasts.csv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)

        metadata = {
            "schema_version": SCHEMA_VERSION,
            "complete": True,
            "software": software_provenance(),
            "creation_command": " ".join(sys.argv),
            "creation_time_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
            "numpy_version": np.__version__,
            "seed": int(args.seed),
            "pairing_repeats": int(args.pairing_repeats),
            "primary_repeat_index": 0,
            "b_type": next(iter(b_types)),
            "null_polarity_design": next(iter(null_polarity_designs)),
            "accepted_null_replicates": {
                category.label: category.accepted_null_replicates
                for category in categories
            },
            "a_type_consistent": a_type_consistent,
            "delta_obs_formula": "Z_{A1} - Z_{A2}",
            "null_delta_formula": (
                "null1[pairing] - null2[second] over R = min(R1, R2) pairs, where "
                "pairing = rng.permutation(R1)[:R], second = rng.permutation(R2)[:R] "
                "when R2 > R (else all of null2), and rng "
                "is seeded from derive_seed(seed, label1, label2, repeat); evaluated "
                "in label-sorted order and sign-flipped to match the reporting order "
                "of (label1, label2), so P is exactly invariant to which label is "
                "given first"
            ),
            "seed_derivation": (
                "int.from_bytes(sha256(f'{seed}|{sorted_label_lo}|{sorted_label_hi}"
                "|{repeat}'.encode()).digest()[:8], 'big')"
            ),
            "p_value_tail_rule": "abs(null_delta) >= abs(delta_obs), two-sided",
            "p_value_formula": "(1 + exceedances) / (R + 1), R = min(R1, R2)",
            "primary_p_is_repeat": 0,
            "sensitivity_p_summary": "min/median/max of P over --pairing-repeats independent pairings",
            "multiplicity_correction": {
                "holm": "step-down over all pairs' primary P-values",
                "benjamini_hochberg": "step-up FDR over all pairs' primary P-values",
            },
            "pair_seeds": pair_seed_records,
            "inputs": [
                {
                    "label": category.label,
                    "directory": str(category.directory.resolve()),
                    "schema_version": category.schema_version,
                    "a_type": category.a_type,
                    "b_type": category.b_type,
                    "accepted_null_replicates": category.accepted_null_replicates,
                    "z_score": category.z_score,
                    "target_digest": category.target_digest,
                    "null_polarity_design": category.null_polarity_design,
                }
                for category in categories
            ],
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
    parser.add_argument(
        "--result", action="append", required=True, metavar="LABEL=DIR",
        help="a phi-sfs-wasserstein-v2 result directory as LABEL=DIR; repeat for "
             "each category (at least two required)",
    )
    parser.add_argument("--output", type=Path, required=True,
                         help="destination directory for the contrast tables; must not exist")
    parser.add_argument(
        "--seed", type=int, default=DEFAULT_SEED,
        help=f"base seed for the deterministic null pairing (default: {DEFAULT_SEED})",
    )
    parser.add_argument(
        "--pairing-repeats", type=int, default=DEFAULT_PAIRING_REPEATS,
        help="independent pairing repeats per pair for a P sensitivity range; "
             f"repeat 0 is the reported primary result (default: {DEFAULT_PAIRING_REPEATS})",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    calculate(args)
    print(f"Wrote {args.output}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
