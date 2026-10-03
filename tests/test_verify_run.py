import json
from pathlib import Path

import numpy as np

from normalize_tes.bootstrap_target_matcher import _candidate_array_digest
from normalize_tes.phi_sfs import SYMMETRIC_NULL_DESIGN
from normalize_tes.sample_age_matched_controls import _sha256_arrays
from normalize_tes.verify_run import _phi_to_pooled, main, verify
from tools.sim_polarity_arms import normalize, phi_many


COMMIT = "1" * 40
STORE_CONTENT = "2" * 64
STORE_CATALOG = "3" * 64
ELIGIBILITY = {
    "vcf_sha256": "4" * 64,
    "heterozygous": "error",
    "min_callable": 20,
    "store_content_sha256": STORE_CONTENT,
}


def _software():
    return {
        "name": "PhiTE",
        "version": "0.9.0",
        "git_commit": COMMIT,
        "git_dirty": False,
    }


def _write_fixture(root: Path) -> dict[str, Path]:
    candidate = root / "candidate_rows.npy"
    target = root / "target"
    matches = root / "matches"
    phi = root / "phi"
    target.mkdir()
    matches.mkdir()
    phi.mkdir()

    candidates = np.arange(10_000, 13_000, dtype=np.int64)
    np.save(candidate, candidates)
    candidate.with_suffix(".npy.json").write_text(json.dumps({
        "software": _software(),
        "store_content_sha256": STORE_CONTENT,
        "store_catalog_sha256": STORE_CATALOG,
        "candidate_rows": int(candidates.size),
        "candidate_rows_sha256": _candidate_array_digest(candidates),
        "universe": "restricted",
        "vcf_eligibility_identity": ELIGIBILITY,
        "inclusion_lists": [{"resolved_fraction": 0.99}],
        "exclusion_lists": [{"resolved_fraction": 0.98}],
        "min_resolved_fraction": 0.70,
    }))

    target_rows = np.arange(4, dtype=np.int64)
    target_cdf = np.asarray([0.2, 0.6, 1.0])
    age_bins = np.asarray([1.0, 2.0, 3.0])
    threshold = 5.0
    np.save(target / "te_row_indices.npy", target_rows)
    np.save(target / "target_cdf.npy", target_cdf)
    np.save(target / "age_bins.npy", age_bins)
    target_meta = {
        "software": _software(),
        "source_store_content_sha256": STORE_CONTENT,
        "source_catalog_sha256": STORE_CATALOG,
        "wasserstein_threshold_generations": threshold,
        "a_type": "TE",
        "te_polarity": None,
        "vcf_eligibility": {"identity": ELIGIBILITY},
    }
    (target / "metadata.json").write_text(json.dumps(target_meta))
    target_digest = _sha256_arrays(
        target_rows, target_cdf, age_bins, np.asarray([threshold], dtype=np.float64)
    )

    replicate_ids = np.arange(500, dtype=np.int64)
    match_rows = candidates[:2_000].reshape(500, 4)
    np.save(matches / "row_indices.npy", match_rows)
    np.save(matches / "replicate_id.npy", replicate_ids)
    np.save(matches / "qc_pass.npy", np.ones(500, dtype=bool))
    np.save(matches / "reuse_row_indices.npy", np.sort(match_rows, axis=None))
    np.save(matches / "reuse_counts.npy", np.ones(match_rows.size, dtype=np.int64))
    match_meta = {
        "software": _software(),
        "complete": True,
        "schema_version": "bootstrap-target-matches-v1",
        "source_store_content_sha256": STORE_CONTENT,
        "source_catalog_sha256": STORE_CATALOG,
        "target_digest": target_digest,
        "candidate_rows_digest": _sha256_arrays(candidates),
        "a_type": "TE",
        "vcf_eligibility_identity": ELIGIBILITY,
        "replicates": 500,
        "set_size": 4,
        "qc_passes": 500,
        "maximum_control_reuse": 1,
        "config": {"disjoint_replicates": True},
    }
    (matches / "metadata.json").write_text(json.dumps(match_meta))

    selected_nulls = list(range(1, 500))
    phi_meta = {
        "software": _software(),
        "complete": True,
        "schema_version": "phi-sfs-wasserstein-v3",
        "target_digest": target_digest,
        "a_type": "TE",
        "vcf_eligibility_identity": ELIGIBILITY,
        "equal_eligible_site_count": 4,
        "null_polarity_design": SYMMETRIC_NULL_DESIGN,
        "matched_sets_published": 500,
        "disjoint_replicates": True,
        "maximum_control_reuse": 1,
        "reference_replicate_id": 0,
        "selected_null_replicate_ids": selected_nulls,
        "accepted_null_replicates": 499,
    }
    (phi / "metadata.json").write_text(json.dumps(phi_meta))
    np.save(phi / "b_replicate_id.npy", replicate_ids)
    spectrum = np.arange(1, 20, dtype=np.float64)
    np.save(phi / "b_raw_sfs.npy", np.tile(spectrum, (500, 1)))
    return {"candidate": candidate, "target": target, "matches": matches, "phi": phi}


def _verify(paths):
    return verify(
        candidate_rows=paths["candidate"],
        target=paths["target"],
        matches=paths["matches"],
        phi=paths["phi"],
        a_type="TE",
        expected_version="0.9.0",
        expected_commit=COMMIT,
    )


def test_verify_accepts_a_complete_production_run(tmp_path):
    report = _verify(_write_fixture(tmp_path))
    assert report["pass"] is True
    assert all(row["pass"] for row in report["criteria"])
    assert report["diagnostics"]["analyzed_sets"] == 500


def test_depletion_distance_matches_the_prespecified_v4_calculation():
    raw = np.random.default_rng(7).uniform(0.1, 2.0, size=(12, 19))
    expected = phi_many(normalize(raw.sum(axis=0)), normalize(raw))
    assert np.allclose(_phi_to_pooled(raw), expected)


def test_verify_rejects_analyzed_set_that_failed_qc(tmp_path):
    paths = _write_fixture(tmp_path)
    qc = np.load(paths["matches"] / "qc_pass.npy")
    qc[0] = False
    np.save(paths["matches"] / "qc_pass.npy", qc)
    metadata_path = paths["matches"] / "metadata.json"
    metadata = json.loads(metadata_path.read_text())
    metadata["qc_passes"] = 499
    metadata_path.write_text(json.dumps(metadata))

    report = _verify(paths)
    criterion = next(
        row for row in report["criteria"]
        if row["criterion"] == "Phi-SFS: every analyzed set passed matching QC"
    )
    assert criterion["pass"] is False
    assert report["pass"] is False


def test_cli_publishes_report_and_returns_failure_status(tmp_path):
    paths = _write_fixture(tmp_path)
    metadata_path = paths["phi"] / "metadata.json"
    metadata = json.loads(metadata_path.read_text())
    metadata["null_polarity_design"] = "wrong"
    metadata_path.write_text(json.dumps(metadata))
    output = tmp_path / "verification"

    status = main([
        "--candidate-rows", str(paths["candidate"]),
        "--target", str(paths["target"]),
        "--matches", str(paths["matches"]),
        "--phi", str(paths["phi"]),
        "--output", str(output),
        "--expected-commit", COMMIT,
    ])

    assert status == 1
    assert (output / "criteria.csv").is_file()
    assert json.loads((output / "report.json").read_text())["pass"] is False
