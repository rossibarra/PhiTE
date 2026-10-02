#!/usr/bin/env python3
"""Keep the first N sets of a disjoint bootstrap-target match bundle.

Matching is sequential: each replicate's seeds derive from the global seed, the
target digest, the replicate and the restart (bootstrap_target_matcher.derive_seed),
and set k is matched from the pool left by sets 0..k-1. So the first N sets of an
R-set run are the sets an N-set run would produce, and trimming stands in for
rerunning the matcher.

Every array whose first axis is the replicate axis is cut to its first N rows,
and so is every CSV row whose `replicate` is below N. The reuse table is
recomputed from the kept row indices. Metadata counts are recomputed, and
`trimmed_from` records the source bundle, its metadata digest, N and the reason.
Arrays without a replicate axis (age bins, target CDF, chromosome labels) are
copied unchanged.

The output directory is built in a staging directory beside it and renamed into
place, so it exists only when complete.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import tempfile
from pathlib import Path

import numpy as np


def trim(source: Path, output: Path, first: int, reason: str) -> None:
    meta_bytes = (source / "metadata.json").read_bytes()
    meta = json.loads(meta_bytes)
    if not meta.get("complete"):
        raise SystemExit(f"{source} is not a complete bundle")
    ids = np.load(source / "replicate_id.npy")
    total = ids.size
    if not np.array_equal(ids, np.arange(total)):
        raise SystemExit("replicate IDs are not 0..R-1 in order")
    if not 0 < first < total:
        raise SystemExit(f"--first must be between 1 and {total - 1}")
    if output.exists():
        raise SystemExit(f"output already exists: {output}")

    staging = Path(tempfile.mkdtemp(prefix=f".{output.name}.tmp.", dir=output.parent))
    try:
        for path in sorted(source.glob("*.npy")):
            if path.name in ("reuse_row_indices.npy", "reuse_counts.npy"):
                continue
            values = np.load(path, mmap_mode="r")
            if values.ndim and values.shape[0] == total:
                values = values[:first]
            np.save(staging / path.name, np.ascontiguousarray(values), allow_pickle=False)

        rows = np.load(staging / "row_indices.npy")
        kept, uses = np.unique(rows, return_counts=True)
        np.save(staging / "reuse_row_indices.npy", kept.astype(np.int64), allow_pickle=False)
        np.save(staging / "reuse_counts.npy", uses.astype(np.uint16), allow_pickle=False)

        for path in sorted(source.glob("*.csv")):
            with path.open(newline="") as src, (staging / path.name).open("w", newline="") as out:
                reader = csv.DictReader(src)
                writer = csv.DictWriter(out, fieldnames=reader.fieldnames)
                writer.writeheader()
                writer.writerows(r for r in reader if int(r["replicate"]) < first)

        qc = np.load(staging / "qc_pass.npy")
        meta.update({
            "replicates": first,
            "qc_passes": int(qc.sum()),
            "qc_failures": int((~qc).sum()),
            "maximum_control_reuse": int(uses.max()),
            "unique_controls_across_sets": int(kept.size),
            "trimmed_from": {
                "source": str(source.resolve()),
                "source_metadata_sha256": hashlib.sha256(meta_bytes).hexdigest(),
                "source_replicates": int(total),
                "first_replicates_kept": first,
                "reason": reason,
            },
        })
        meta["config"] = dict(meta["config"], replicates=first)
        with (staging / "metadata.json").open("w") as handle:
            json.dump(meta, handle, indent=2)
            handle.write("\n")
        os.replace(staging, output)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    print(f"kept replicates 0-{first - 1} of {total}: {int(qc.sum())} pass QC, "
          f"{kept.size} unique controls, maximum reuse {int(uses.max())} -> {output}")


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--matches", type=Path, required=True, help="source match bundle")
    p.add_argument("--first", type=int, required=True, help="number of leading sets to keep")
    p.add_argument("--reason", required=True, help="why the bundle is trimmed (recorded)")
    p.add_argument("--output", type=Path, required=True, help="new bundle directory")
    a = p.parse_args(argv)
    trim(a.matches, a.output, a.first, a.reason)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
