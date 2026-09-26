import csv
import json
import sys
from pathlib import Path

import numpy as np
import pytest

from normalize_tes.bootstrap_target_matcher import (
    OptimizerConfig,
    _candidate_array_digest,
    bootstrap_cdf,
    bootstrap_counts,
    derive_seed,
    disjoint_stratum_capacity,
    log_search_grid,
    main,
    median_age_strata,
    optimize_restart,
    parse_args,
    validate_restart_result,
)
from normalize_tes.snp_age_store import open_snp_age_store
from normalize_tes.snp_interval_dataset import INTERVAL_SCHEMA_VERSION, pack_status
from normalize_tes.swap_control_sampler import analysis_points, eligible_candidates, search_grid
from test_swap_control_sampler import _interval_store, _target as _base_target


ELIGIBILITY_IDENTITY = {
    "vcf_sha256": "v" * 64,
    "heterozygous": "error",
    "min_callable": 20,
    "store_content_sha256": "a" * 64,
    "row_indices_sha256": "r" * 64,
    "snp_row_indices_sha256": "s" * 64,
    "p_alt_derived_sha256": "p" * 64,
}


def _target(path, store_path, **kwargs):
    """The shared fixture target, plus the eligibility identity it now records."""
    target = _base_target(path, store_path, **kwargs)
    metadata_path = target / "metadata.json"
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    metadata["vcf_eligibility"]["identity"] = dict(ELIGIBILITY_IDENTITY)
    metadata_path.write_text(json.dumps(metadata), encoding="utf-8")
    return target


def _strata_store(path, below):
    """An interval store with one 10-generation interval per SNP.

    Like `_interval_store`, but the lower endpoints are chosen by the caller
    so a test can place candidates in particular target age strata.
    """
    below = np.asarray(below, dtype=np.float64)
    count = below.size
    above = below + 10
    arrays = {
        "positions": np.arange(1, count + 1, dtype=np.float64),
        "offsets": np.arange(count + 1, dtype=np.uint64),
        "below": below,
        "above": above,
        "draw_id": np.zeros(count, dtype=np.uint8),
        "status": pack_status(np.full((1, count), 2, dtype=np.uint8)),
        "present_draw_count": np.ones(count, dtype=np.uint32),
        "missing_draw_count": np.zeros(count, dtype=np.uint32),
        "usable_draw_count": np.ones(count, dtype=np.uint32),
        "usable_interval_count": np.ones(count, dtype=np.uint32),
        "skipped_root_count": np.zeros(count, dtype=np.uint32),
    }
    path.mkdir()
    for name, values in arrays.items():
        np.save(path / f"{name}.npy", values)
    (path / "metadata.json").write_text(json.dumps({
        "schema_version": INTERVAL_SCHEMA_VERSION,
        "n_snps": count,
        "n_intervals": count,
        "n_posterior_draws": 1,
        "maximum_above": float(above.max()),
        "endpoint_dtype": "float64",
        "minimum_usable_draws": 1,
        "arrays": {
            name: {"dtype": value.dtype.name, "shape": list(value.shape)}
            for name, value in arrays.items()
        },
        "chromosomes": [{"chrom": "1", "offset": 0, "length": 100}],
        "catalog_sha256": "fixture-catalog",
        "content_sha256": "a" * 64,
    }), encoding="utf-8")
    return path


def _candidate_file(tmp_path, rows, *, identity=None):
    """A candidate array with the provenance report build_candidate_rows writes."""
    path = tmp_path / "candidates.npy"
    array = np.asarray(rows, dtype=np.int64)
    np.save(path, array)
    report = {
        "store_content_sha256": "a" * 64,
        "store_catalog_sha256": "fixture-catalog",
        "candidate_rows": int(array.size),
        "candidate_rows_sha256": _candidate_array_digest(array),
        "vcf_eligibility_identity": (
            dict(ELIGIBILITY_IDENTITY) if identity is None else identity
        ),
    }
    path.with_suffix(".npy.json").write_text(json.dumps(report), encoding="utf-8")
    return path


def test_default_replicate_count_provides_reference_and_null_sets():
    assert OptimizerConfig().replicates == 1201
    args = parse_args([
        "--store", "store", "--target", "target", "--all-eligible",
        "--output", "output", "--disjoint-replicates",
    ])
    assert args.replicates == 1201
    assert args.disjoint_replicates
    assert args.a_type == "TE"


def test_a_type_must_match_the_filtered_target(tmp_path):
    store = _interval_store(tmp_path / "store")
    target = _target(tmp_path / "target", store)
    with pytest.raises(ValueError, match="disagrees with target metadata"):
        _run_matcher(store, target, tmp_path / "output", "-A", "SNP")


def test_snp_a_uses_snp_filtered_target_without_te_polarity(tmp_path):
    store = _interval_store(tmp_path / "store")
    target = _target(tmp_path / "target", store)
    metadata_path = target / "metadata.json"
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    metadata["a_type"] = "SNP"
    metadata["te_polarity"] = None
    metadata_path.write_text(json.dumps(metadata), encoding="utf-8")
    (target / "te_keep_draws.npy").unlink()

    output = tmp_path / "output"
    assert _run_matcher(store, target, output, "-A", "SNP") == 0
    published = json.loads((output / "metadata.json").read_text(encoding="utf-8"))
    assert published["a_type"] == "SNP"


def test_te_a_requires_inclusive_half_derived_filter(tmp_path):
    store = _interval_store(tmp_path / "store")
    target = _target(tmp_path / "target", store)
    metadata_path = target / "metadata.json"
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    metadata["te_polarity"]["max_flipped_fraction"] = 0.49
    metadata_path.write_text(json.dumps(metadata), encoding="utf-8")
    with pytest.raises(ValueError, match="max_flipped_fraction=0.5"):
        _run_matcher(store, target, tmp_path / "output")


def test_disjoint_capacity_fails_before_creating_work_state(tmp_path):
    store = _interval_store(tmp_path / "store")
    target = _target(tmp_path / "target", store)
    output = tmp_path / "output"
    work = tmp_path / "work"
    with pytest.raises(ValueError, match=r"requires 6 unique candidate SNPs"):
        _run_matcher(
            store, target, output,
            "--replicates", "3", "--disjoint-replicates",
            "--work-dir", str(work),
        )
    assert not output.exists()
    assert not work.exists()


# Two target TEs with medians near 5 and 17 generations put the target's
# quota in strata 1 and 2 (one site each); candidates with lower endpoints 1
# and 13 fall in those same strata.
_YOUNG, _OLD = 1, 13


def test_median_age_strata_is_shared_and_chunk_invariant(tmp_path):
    store_path = _strata_store(tmp_path / "store", [0, 12] + [_YOUNG] * 6 + [_OLD] * 6)
    target = _target(tmp_path / "target", store_path)
    store = open_snp_age_store(store_path)
    boundaries = np.load(target / "interval_boundary_ages.npy")
    quotas = np.load(target / "interval_quotas.npy")
    candidates = np.arange(2, 14)
    strata = median_age_strata(store, candidates, boundaries, quotas.size)
    np.testing.assert_array_equal(
        strata, median_age_strata(store, candidates, boundaries, quotas.size,
                                  chunk_rows=1, block_rows=1),
    )
    np.testing.assert_array_equal(strata, [1] * 6 + [2] * 6)
    counts, capacity = disjoint_stratum_capacity(store, candidates, boundaries, quotas)
    np.testing.assert_array_equal(quotas, [0, 1, 1, 0])
    np.testing.assert_array_equal(counts, [0, 6, 6, 0])
    np.testing.assert_array_equal(capacity, [np.inf, 6.0, 6.0, np.inf])


def test_disjoint_stratum_capacity_passes_when_every_stratum_suffices(tmp_path):
    store = _strata_store(tmp_path / "store", [0, 12] + [_YOUNG] * 6 + [_OLD] * 6)
    target = _target(tmp_path / "target", store)
    candidates = _candidate_file(tmp_path, np.arange(2, 14))
    output = tmp_path / "output"
    assert _run_matcher(
        store, target, output, "--replicates", "3", "--disjoint-replicates",
        candidates=candidates,
    ) == 0
    metadata = json.loads((output / "metadata.json").read_text(encoding="utf-8"))
    assert metadata["stratum_quotas"] == [0, 1, 1, 0]
    assert metadata["disjoint_stratum_candidates"] == [0, 6, 6, 0]
    assert metadata["disjoint_stratum_capacity_sets"] == [None, 6.0, 6.0, None]
    assert metadata["maximum_control_reuse"] == 1
    assert metadata["vcf_eligibility_identity"] == ELIGIBILITY_IDENTITY


def test_disjoint_capacity_fails_when_one_stratum_is_thin(tmp_path):
    """Ten candidates hold 3 x 2 in total, but only two are young."""
    store = _strata_store(tmp_path / "store", [0, 12] + [_YOUNG] * 2 + [_OLD] * 8)
    target = _target(tmp_path / "target", store)
    candidates = _candidate_file(tmp_path, np.arange(2, 12))
    output, work = tmp_path / "output", tmp_path / "work"
    with pytest.raises(ValueError) as caught:
        _run_matcher(
            store, target, output, "--replicates", "3", "--disjoint-replicates",
            "--work-dir", str(work), candidates=candidates,
        )
    message = str(caught.value)
    assert "1 of 4 strata fall short" in message
    assert "stratum 1: 2 candidates / quota 1 = 2.0 sets" in message
    assert "stratum 2" not in message
    assert "will not be silently reused" in message
    assert not output.exists()
    assert not work.exists()


def test_candidate_and_target_eligibility_identities_must_agree(tmp_path):
    store = _strata_store(tmp_path / "store", [0, 12] + [_YOUNG] * 6 + [_OLD] * 6)
    target = _target(tmp_path / "target", store)
    other = {**ELIGIBILITY_IDENTITY, "vcf_sha256": "w" * 64}
    candidates = _candidate_file(tmp_path, np.arange(2, 14), identity=other)
    output = tmp_path / "output"
    with pytest.raises(ValueError, match=r"different VCF-eligibility artifacts "
                                         r"\(differing: vcf_sha256\)"):
        _run_matcher(store, target, output, candidates=candidates)
    assert not output.exists()


def test_candidate_report_without_eligibility_identity_is_rejected(tmp_path):
    store = _strata_store(tmp_path / "store", [0, 12] + [_YOUNG] * 6 + [_OLD] * 6)
    target = _target(tmp_path / "target", store)
    incomplete = {k: v for k, v in ELIGIBILITY_IDENTITY.items() if k != "heterozygous"}
    candidates = _candidate_file(tmp_path, np.arange(2, 14), identity=incomplete)
    with pytest.raises(ValueError, match="vcf_eligibility_identity lacks heterozygous"):
        _run_matcher(store, target, tmp_path / "output", candidates=candidates)


def test_target_without_eligibility_identity_is_rejected(tmp_path):
    store = _interval_store(tmp_path / "store")
    target = _base_target(tmp_path / "target", store)
    with pytest.raises(ValueError, match="vcf_eligibility.identity is missing"):
        _run_matcher(store, target, tmp_path / "output")


def test_bootstrap_counts_and_cdf_are_reproducible():
    rows = np.array([[0.0, 0.5, 1.0], [0.2, 0.8, 1.0]])
    first = bootstrap_counts(2, np.random.default_rng(7))
    second = bootstrap_counts(2, np.random.default_rng(7))
    np.testing.assert_array_equal(first, second)
    assert first.sum() == 2
    np.testing.assert_allclose(
        bootstrap_cdf(np.array([1, 1]), rows), rows.mean(axis=0)
    )
    with pytest.raises(ValueError, match="aligned"):
        bootstrap_cdf(np.array([1]), rows)


def test_seed_derivation_is_stable():
    assert derive_seed(3, "x", 2) == derive_seed(3, "x", 2)
    assert derive_seed(3, "x", 2, 0) != derive_seed(3, "x", 2, 1)
    # The bootstrap seed and every restart seed of the same replicate differ,
    # so restarts are independent stratified draws rather than repeats.
    assert derive_seed(3, "x", 2) != derive_seed(3, "x", 2, 0)
    assert derive_seed(3, "x", 2, 0) != derive_seed(3, "y", 2, 0)


def test_exact_optimizer_trace_is_monotone_and_certified(tmp_path):
    store_path = _interval_store(tmp_path / "store")
    target_path = _target(tmp_path / "target", store_path)
    store = open_snp_age_store(store_path)
    target_rows = np.load(target_path / "te_row_indices.npy")
    target = np.load(target_path / "target_cdf.npy")
    ages = np.load(target_path / "age_bins.npy")
    candidates = eligible_candidates(store, target_rows)
    config = OptimizerConfig(
        replicates=1, restarts=1,
        min_epochs=2, max_epochs=4, patience=2,
        material_improvement_ratio=0, cdf_block_rows=2,
        search_bin_width=int(ages[1] - ages[0]),
        qc_max_ratio=1, qc_max_absolute=100,
    )
    coarse_ages, coarse_points = search_grid(
        float(store.metadata["maximum_above"]), config.search_bin_width
    )
    coarse_target = store.aggregate_cdf_at(
        target_rows, coarse_points, side="left", weighting="interval"
    )
    result = optimize_restart(
        store, candidates, np.array([4, 5]), target, target, ages, 10,
        coarse_target=coarse_target,
        coarse_ages=coarse_ages,
        coarse_points=coarse_points,
        seed=19, config=config,
    )
    best = np.array([record["best_distance"] for record in result.trace])
    assert np.all(np.diff(best) <= 1e-12)
    certified = store.aggregate_cdf_at(
        result.rows, analysis_points(ages), side="left", weighting="interval"
    )
    np.testing.assert_allclose(result.cdf, certified)
    assert result.best_distance <= result.initial_distance
    assert np.unique(result.rows).size == result.rows.size
    validate_restart_result(
        result,
        store=store,
        candidates=candidates,
        bootstrap_target=target,
        observed_target=target,
        age_bins=ages,
        expected_seed=19,
    )


def _run_matcher(store, target, output, *extra, candidates=None):
    universe = (
        ["--all-eligible"] if candidates is None
        else ["--candidate-rows", str(candidates)]
    )
    return main([
        "--store", str(store), "--target", str(target),
        *universe, "--output", str(output),
        "--replicates", "2", "--restarts", "2",
        "--min-epochs", "2", "--max-epochs", "3", "--patience", "1",
        "--material-improvement-ratio", "0", "--cdf-block-rows", "2",
        "--search-bin-width", "10",
        "--qc-max-ratio", "10", "--qc-max-absolute", "1000", "--seed", "23",
        *extra,
    ])


def test_bootstrap_bundle_is_consumable_by_phi_sfs(tmp_path):
    """The whole point of the bundle is to be read by the Phi-SFS step.

    Rounds 7 found the two ways this silently failed: a target_digest computed
    over different arrays than phi_sfs recomputes, and missing per-replicate
    identifier arrays. Both were invisible to every other test.
    """
    from normalize_tes import phi_sfs
    from normalize_tes.snp_age_store import open_snp_age_store

    store = _interval_store(tmp_path / "store")
    target = _target(tmp_path / "target", store)
    output = tmp_path / "bootstrap_matches"
    assert _run_matcher(store, target, output) == 0

    # The digest must certify against the target directory ...
    target_meta, match_meta, digest = phi_sfs._validate_provenance(target, output)
    assert match_meta["target_digest"] == digest
    # ... and the identifiers must load under this bundle's own schema.
    assert match_meta["schema_version"] == "bootstrap-target-matches-v1"
    opened = open_snp_age_store(store)
    chromosomes, positions = opened.rows_to_native(
        np.load(target / "te_row_indices.npy")
    )
    np.save(target / "te_chromosomes.npy", chromosomes)
    np.save(target / "te_positions.npy", positions)
    *_, identifiers = phi_sfs._load_coordinates(
        target, output, match_meta["schema_version"]
    )
    assert list(identifiers) == ["replicate_id"]
    assert identifiers["replicate_id"].tolist() == [0, 1]
    # Chain/sample identifiers must NOT be invented for independent replicates.
    assert not (output / "chain_index.npy").exists()
    assert not (output / "sample_index.npy").exists()


def test_output_inside_work_dir_is_rejected(tmp_path):
    store = _interval_store(tmp_path / "store")
    target = _target(tmp_path / "target", store)
    work = tmp_path / "work"
    with pytest.raises(ValueError, match="must not be inside"):
        _run_matcher(store, target, work / "nested", "--work-dir", str(work))
    with pytest.raises(ValueError, match="different paths"):
        _run_matcher(store, target, work, "--work-dir", str(work))


def test_resume_rejects_a_changed_implementation(tmp_path):
    store = _interval_store(tmp_path / "store")
    target = _target(tmp_path / "target", store)
    work = tmp_path / "work"
    assert _run_matcher(
        store, target, tmp_path / "out",
        "--work-dir", str(work), "--keep-work",
    ) == 0
    identity = json.loads((work / "identity.json").read_text())
    assert identity["software"]["name"] == "PhiTE"
    assert identity["numpy_version"]
    identity["software"]["git_commit"] = "0" * 40
    (work / "identity.json").write_text(json.dumps(identity, indent=2, sort_keys=True))
    with pytest.raises(ValueError, match="parameters or provenance differ"):
        _run_matcher(
            store, target, tmp_path / "out2",
            "--work-dir", str(work), "--resume",
        )


def test_resume_reuses_completed_replicate_bundles(tmp_path):
    store = _interval_store(tmp_path / "store")
    target = _target(tmp_path / "target", store)
    work = tmp_path / "work"
    assert _run_matcher(
        store, target, tmp_path / "first",
        "--work-dir", str(work), "--keep-work",
    ) == 0
    first = np.load(tmp_path / "first" / "row_indices.npy")
    assert _run_matcher(
        store, target, tmp_path / "second",
        "--work-dir", str(work), "--resume",
    ) == 0
    np.testing.assert_array_equal(
        first, np.load(tmp_path / "second" / "row_indices.npy")
    )


def test_source_store_provenance_is_the_store_not_the_repository(tmp_path):
    store = _interval_store(tmp_path / "store")
    target = _target(tmp_path / "target", store)
    output = tmp_path / "out"
    assert _run_matcher(store, target, output) == 0
    metadata = json.loads((output / "metadata.json").read_text())
    assert metadata["source_store"] == str(store.resolve())
    # The seed-library path is gone, so the bundle must not claim one.
    assert metadata["initialisation"].startswith("stratified draw")
    assert "seed_sets" not in metadata and "seed_sets_digest" not in metadata


def test_bootstrap_cdf_accumulates_in_float64(tmp_path):
    """float32 accumulation over many TE rows displaces the bootstrap target."""
    rng = np.random.default_rng(0)
    n_sites, grid = 40_000, 64
    rows32 = rng.random((n_sites, grid)).astype(np.float32)
    rows32.sort(axis=1)
    counts = bootstrap_counts(n_sites, np.random.default_rng(5))
    reference = (
        counts.astype(np.float64) @ rows32.astype(np.float64) / counts.sum()
    )
    np.testing.assert_allclose(
        bootstrap_cdf(counts, rows32), reference, rtol=0, atol=1e-12
    )


def test_cli_writes_aligned_atomic_bundle(tmp_path):
    store = _interval_store(tmp_path / "store")
    target = _target(tmp_path / "target", store)
    output = tmp_path / "bootstrap_matches"
    assert main([
        "--store", str(store),
        "--target", str(target),
        "--all-eligible",
        "--output", str(output),
        "--replicates", "2",
        "--restarts", "2",
        "--min-epochs", "2",
        "--max-epochs", "3",
        "--patience", "1",
        "--material-improvement-ratio", "0",
        "--cdf-block-rows", "2",
        # The fixture store spans 27 generations, so the coarse log screen
        # needs a search width it can put at least two points inside.
        "--search-bin-width", "10",
        "--qc-max-ratio", "10",
        "--qc-max-absolute", "1000",
        "--seed", "23",
    ]) == 0
    rows = np.load(output / "row_indices.npy")
    assert rows.shape == (2, 2)
    assert np.load(output / "bootstrap_counts.npy").shape == (2, 2)
    assert np.load(output / "bootstrap_target_cdfs.npy").shape == (2, 4)
    assert np.load(output / "restart_best_rows.npy").shape == (2, 2, 2)
    assert np.all(np.load(output / "triangle_ok.npy"))
    with (output / "replicates.csv").open() as handle:
        assert len(list(csv.DictReader(handle))) == 2
    metadata = json.loads((output / "metadata.json").read_text())
    assert metadata["complete"] is True
    assert metadata["replicates"] == 2
    assert metadata["restarts_per_replicate"] == 2
    assert metadata["qc_interpretation"].startswith("optimizer convergence")
    assert not any(path.name.startswith(f".{output.name}.tmp") for path in tmp_path.iterdir())


def test_log_search_grid_resolves_the_young_end_without_costing_more_points():
    """The coarse screen has to see the region the log-age metric prices.

    A uniform 20,000-generation screen puts the whole lower quartile of a
    production TE age distribution inside its first cell, so under log-age
    weighting the optimizer would be blind to exactly the ages the metric is
    meant to price. The geometric screen keeps the young end at full exact
    resolution and is no more expensive than the uniform one.
    """
    ages = np.arange(36_746, dtype=np.float64) * 1_000.0
    points = ages + 500.0
    linear_points = int(36_745_000 // 20_000) + 1
    coarse_ages, coarse_points = log_search_grid(
        ages, points, linear_points, 1_000.0
    )
    assert coarse_ages.size <= linear_points
    assert np.all(np.diff(coarse_ages) > 0)
    assert np.all(np.isin(coarse_ages, ages))
    np.testing.assert_allclose(coarse_points, coarse_ages + 500.0)
    assert coarse_ages[0] == ages[0] and coarse_ages[-1] == ages[-1]
    # Full exact resolution across the lower three quartiles of a production
    # in-gene TE age distribution (q75 is about 80,000 generations).
    fine = np.flatnonzero(np.diff(coarse_ages) > 1_000)
    assert coarse_ages[fine[0]] > 100_000
    with pytest.raises(ValueError, match="at least two"):
        log_search_grid(ages, points, 1, 1_000.0)


def test_resume_rejects_a_different_dirty_source_state(tmp_path, monkeypatch):
    """Two dirty edits on one commit must not share a resume identity.

    `software_provenance` records the HEAD commit and a dirty flag, which are
    identical for any two sets of uncommitted edits. Without hashing the loaded
    modules a long job could resume across an implementation change and mix
    replicate bundles from two versions of the code.
    """
    from normalize_tes import release_provenance

    first = release_provenance.loaded_source_digest()
    assert first["sha256"] and first["modules"]

    # Simulate an edit to one loaded module by pointing the digest at a copy of
    # the repo whose matcher source differs by a single byte.
    shadow = tmp_path / "shadow"
    shadow.mkdir()
    root = Path(release_provenance.__file__).resolve().parent.parent
    for name in first["modules"]:
        destination = shadow / name
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes((root / name).read_bytes())
    target = shadow / "normalize_tes/bootstrap_target_matcher.py"
    target.write_bytes(target.read_bytes() + b"\n# one byte of drift\n")

    real_modules = {}
    for name in first["modules"]:
        stem = Path(name).stem
        module_name = "normalize_tes" if stem == "__init__" else f"normalize_tes.{stem}"
        if module_name in sys.modules:
            real_modules[name] = sys.modules[module_name]
    for name, module in real_modules.items():
        monkeypatch.setattr(module, "__file__", str(shadow / name), raising=False)
    second = release_provenance.loaded_source_digest(shadow)

    assert second["modules"] == first["modules"]
    assert second["sha256"] != first["sha256"], (
        "a one-byte change to a loaded module must change the resume identity"
    )
