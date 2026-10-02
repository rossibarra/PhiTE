#!/usr/bin/env python3
"""V4 report: blocking criteria for full production matching under depletion.

Validation item V4 (docs/REMAINING_VALIDATION_PROPOSAL.md) is judged on the
matched-control bundle, plus, for the SFS-drift criterion, the Phi-SFS output
run on that bundle (V5), which holds one spectrum per QC-passing set.

Definitions the proposal leaves open are fixed here:
  blocks        replicates 0-99, 100-199, ...; only complete blocks of 100 count
  QC fall       pass rate of the first block minus that of the last complete block
  SFS distance  phi (as in tools.sim_polarity_arms.phi_many) between each
                published set's normalized spectrum and the pooled spectrum, the
                normalized sum of every published set's raw spectrum
  rho           Spearman correlation with replicate ID (average ranks for ties)
  early/late    first and last quarter of the published sets in replicate order,
                the convention tools.sim_polarity_arms uses for depletion order
  SMD           (late mean - early mean) / sqrt((early var + late var) / 2),
                sample variances; it passes when |SMD| < 0.2

Without --phi the drift criterion is reported as pending.

Outputs (new directory):
  criteria.csv     one row per blocking criterion: value, threshold, pass
  blocks.csv       QC passes and W1 / error-ratio medians per block of 100
  report.json      everything above plus the inputs used
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from tools.sim_polarity_arms import normalize, phi_many

BLOCK = 100


def average_ranks(x: np.ndarray) -> np.ndarray:
    order = np.argsort(x, kind="mergesort")
    ranks = np.empty(x.size, dtype=np.float64)
    ranks[order] = np.arange(1, x.size + 1)
    values, inverse = np.unique(x, return_inverse=True)
    sums = np.bincount(inverse, weights=ranks)
    counts = np.bincount(inverse)
    return (sums / counts)[inverse]


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    return float(np.corrcoef(average_ranks(x), average_ranks(y))[0, 1])


def smd(early: np.ndarray, late: np.ndarray) -> float:
    pooled = np.sqrt((early.var(ddof=1) + late.var(ddof=1)) / 2)
    return float((late.mean() - early.mean()) / pooled)


def matching_criteria(matches: Path, min_qc: int) -> tuple[list[dict], list[dict], dict]:
    meta = json.loads((matches / "metadata.json").read_text())
    if not meta.get("complete"):
        raise SystemExit(f"{matches} is not a complete bundle")
    rows = np.load(matches / "row_indices.npy")
    ids = np.load(matches / "replicate_id.npy")
    qc = np.load(matches / "qc_pass.npy").astype(bool)
    w1 = np.load(matches / "match_to_bootstrap_w1.npy")
    ratio = np.load(matches / "matching_error_ratio.npy")
    if not np.array_equal(ids, np.arange(ids.size)):
        raise SystemExit("replicate IDs are not 0..N-1 in order")

    duplicate_sets = int(sum(np.unique(r).size != r.size for r in rows))
    _, uses = np.unique(rows, return_counts=True)
    blocks = []
    for b in range(ids.size // BLOCK):
        s = slice(b * BLOCK, (b + 1) * BLOCK)
        blocks.append({
            "block": b, "first_replicate": b * BLOCK, "last_replicate": (b + 1) * BLOCK - 1,
            "qc_passes": int(qc[s].sum()),
            "median_match_to_bootstrap_w1": float(np.median(w1[s])),
            "median_matching_error_ratio": float(np.median(ratio[s])),
        })
    passes = [b["qc_passes"] for b in blocks]
    fall = (passes[0] - passes[-1]) / BLOCK
    criteria = [
        {"criterion": "maximum control reuse", "value": int(uses.max()),
         "threshold": "== 1", "pass": int(uses.max()) == 1},
        {"criterion": "sets containing a duplicate control", "value": duplicate_sets,
         "threshold": "== 0", "pass": duplicate_sets == 0},
        {"criterion": "QC-passing sets", "value": int(qc.sum()),
         "threshold": f">= {min_qc}", "pass": int(qc.sum()) >= min_qc},
        {"criterion": "fewest QC passes in a complete block of 100", "value": min(passes),
         "threshold": ">= 80", "pass": min(passes) >= 80},
        {"criterion": "QC pass-rate fall, first to last block (points)",
         "value": round(100 * fall, 2), "threshold": "<= 10", "pass": fall <= 0.10},
    ]
    diagnostics = {
        "rho_match_to_bootstrap_w1": spearman(ids.astype(float), w1),
        "rho_matching_error_ratio": spearman(ids.astype(float), ratio),
        "median_w1_shift_first_to_last_block":
            blocks[-1]["median_match_to_bootstrap_w1"] - blocks[0]["median_match_to_bootstrap_w1"],
        "median_error_ratio_shift_first_to_last_block":
            blocks[-1]["median_matching_error_ratio"] - blocks[0]["median_matching_error_ratio"],
        "metadata_qc_passes": meta["qc_passes"],
        "metadata_maximum_control_reuse": meta["maximum_control_reuse"],
        "replicates": int(ids.size), "set_size": int(rows.shape[1]),
    }
    return criteria, blocks, diagnostics


def drift_criteria(phi: Path) -> tuple[list[dict], dict]:
    meta = json.loads((phi / "metadata.json").read_text())
    if not meta.get("complete"):
        raise SystemExit(f"{phi} is not a complete Phi-SFS output")
    ids = np.load(phi / "b_replicate_id.npy")
    raw = np.load(phi / "b_raw_sfs.npy")
    order = np.argsort(ids)
    ids, raw = ids[order], raw[order]
    distance = phi_many(normalize(raw.sum(axis=0)), normalize(raw))
    quarter = max(1, ids.size // 4)
    rho = spearman(ids.astype(float), distance)
    effect = smd(distance[:quarter], distance[-quarter:])
    criteria = [
        {"criterion": "SFS distance to pooled spectrum: |rho| with replicate ID",
         "value": round(rho, 4), "threshold": "|rho| < 0.1", "pass": abs(rho) < 0.1},
        {"criterion": "SFS distance to pooled spectrum: early-vs-late |SMD|",
         "value": round(effect, 4), "threshold": "|SMD| < 0.2", "pass": abs(effect) < 0.2},
    ]
    return criteria, {"published_sets": int(ids.size), "quarter_size": quarter,
                      "early_mean_distance": float(distance[:quarter].mean()),
                      "late_mean_distance": float(distance[-quarter:].mean())}


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--matches", type=Path, required=True, help="V4 matched-control bundle")
    p.add_argument("--phi", type=Path, help="Phi-SFS output on that bundle (V5)")
    p.add_argument("--min-qc-passes", type=int, default=901,
                   help="QC-passing sets required (901 of 1,001 in the frozen V4)")
    p.add_argument("--output", type=Path, required=True, help="new output directory")
    return p.parse_args(argv)


def run(a: argparse.Namespace) -> None:
    if a.output.exists():
        raise SystemExit(f"output already exists: {a.output}")
    criteria, blocks, diagnostics = matching_criteria(a.matches, a.min_qc_passes)
    if a.phi is None:
        criteria.append({"criterion": "SFS drift (needs --phi)", "value": "",
                         "threshold": "", "pass": "pending"})
    else:
        drift, drift_info = drift_criteria(a.phi)
        criteria += drift
        diagnostics.update(drift_info)
    a.output.mkdir(parents=True)
    write_csv(a.output / "criteria.csv", criteria)
    write_csv(a.output / "blocks.csv", blocks)
    report = {"matches": str(a.matches), "phi": str(a.phi) if a.phi else None,
              "min_qc_passes": a.min_qc_passes,
              "criteria": criteria, "blocks": blocks, "diagnostics": diagnostics}
    (a.output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    for c in criteria:
        print(f"{str(c['pass']):8} {c['criterion']}: {c['value']} ({c['threshold']})")
    print(json.dumps(diagnostics, indent=2))


def main(argv=None) -> int:
    run(parse_args(argv))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
