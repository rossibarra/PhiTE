#!/usr/bin/env python3
"""Sample one V6 negative-control focal SNP set from the candidate universe.

Validation item V6 (docs/REMAINING_VALIDATION_PROPOSAL.md) runs 300 tests, each
with a focal set of exactly M distinct SNPs drawn independently from the full
genome-wide candidate universe, without replacement within a set. Focal sets
are not made disjoint from one another; the overlap is reported afterwards.

Writes a two-column position file (chromosome, native 1-based VCF position) that
`normalize_tes.te_age_target --te-positions` and
`normalize_tes.build_candidate_rows --exclude-positions` read, plus a JSON
sidecar recording the seed and the sampled store rows.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from normalize_tes.release_provenance import software_provenance
from normalize_tes.snp_age_store import open_snp_age_store


def sample(candidate_rows: np.ndarray, m: int, seed: int) -> np.ndarray:
    if np.unique(candidate_rows).size != candidate_rows.size:
        raise ValueError("candidate rows contain duplicates")
    if m > candidate_rows.size:
        raise ValueError(f"cannot sample {m} of {candidate_rows.size} candidates")
    picked = np.random.default_rng(seed).choice(candidate_rows, size=m, replace=False)
    return np.sort(picked)


def native_positions(store: object, rows: np.ndarray) -> tuple[list[str], np.ndarray]:
    """Map store rows to (chromosome, native position) via the store's offsets."""
    positions = np.asarray(store.positions)[rows]
    chroms = sorted(store.metadata["chromosomes"], key=lambda c: int(c["offset"]))
    offsets = np.array([int(c["offset"]) for c in chroms], dtype=np.float64)
    which = np.searchsorted(offsets, positions, side="right") - 1
    if np.any(which < 0):
        raise ValueError("a sampled row lies before the first chromosome offset")
    native = positions - offsets[which]
    for index in np.unique(which):
        length = int(chroms[index]["length"])
        if np.any(native[which == index] > length) or np.any(native[which == index] < 1):
            raise ValueError(f"a sampled row falls outside chromosome {chroms[index]['chrom']}")
    if np.any(native != np.round(native)):
        raise ValueError("a sampled store position is not an integer coordinate")
    return [str(chroms[i]["chrom"]) for i in which], native.astype(np.int64)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--store", type=Path, required=True)
    p.add_argument("--candidate-rows", type=Path, required=True)
    p.add_argument("--m", type=int, default=4000)
    p.add_argument("--seed", type=int, required=True)
    p.add_argument("--output", type=Path, required=True, help="new position file")
    args = p.parse_args(argv)
    if args.output.exists():
        raise SystemExit(f"output already exists: {args.output}")
    candidates = np.load(args.candidate_rows, allow_pickle=False).astype(np.int64)
    rows = sample(candidates, args.m, args.seed)
    chroms, native = native_positions(open_snp_age_store(args.store), rows)
    staged = args.output.with_name(f".{args.output.name}.tmp")
    staged.write_text("".join(f"{c}\t{n}\n" for c, n in zip(chroms, native)))
    args.output.with_suffix(args.output.suffix + ".json").write_text(json.dumps({
        "m": args.m, "seed": args.seed,
        "candidate_rows": str(args.candidate_rows.resolve()),
        "candidate_rows_sha256": hashlib.sha256(candidates.tobytes()).hexdigest(),
        "sampled_rows": rows.tolist(),
        "software": software_provenance(),
    }) + "\n")
    staged.replace(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
