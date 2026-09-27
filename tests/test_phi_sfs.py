import csv
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from normalize_tes import phi_sfs as phi_sfs_module
from normalize_tes.phi_sfs import (
    PROJECTION_SIZE,
    RETAINED_BINS,
    SiteCount,
    accumulate_spectrum,
    calibrate_phi,
    hypergeometric_projection,
    main,
    normalized_spectrum,
    phi_sfs,
    project_sites,
    project_sites_bernoulli_q,
)
from normalize_tes.sample_age_matched_controls import _sha256_arrays


# ---------------------------------------------------------------- projection


def test_projection_probability_and_point_mass():
    projected = hypergeometric_projection(7, 20)
    assert projected.sum() == pytest.approx(1)
    assert projected[7] == pytest.approx(1)
    assert np.count_nonzero(projected > 1e-12) == 1


def test_projection_known_case_and_validation():
    projected = hypergeometric_projection(1, 21)
    assert projected[0] == pytest.approx(1 / 21)
    assert projected[1] == pytest.approx(20 / 21)
    with pytest.raises(ValueError, match="cannot project"):
        hypergeometric_projection(3, 19)
    with pytest.raises(ValueError, match="0 <= k"):
        hypergeometric_projection(22, 21)


def test_projection_matches_random_subsampling():
    """Cross-check the closed form against the definition it stands for."""
    k, n, draws = 7, 30, 40_000
    rng = np.random.default_rng(0)
    alleles = np.zeros(n, dtype=np.int64)
    alleles[:k] = 1
    chosen = rng.random((draws, n)).argsort(axis=1)[:, :20]
    empirical = np.bincount(alleles[chosen].sum(axis=1), minlength=21) / draws
    assert np.allclose(empirical, hypergeometric_projection(k, n), atol=0.005)


def test_projection_stable_at_large_n():
    projected = hypergeometric_projection(1, 2_000_000)
    assert projected.sum() == pytest.approx(1)
    assert projected[0] == pytest.approx(1 - 1e-5, abs=1e-9)


# --------------------------------------------------------- site projection


def test_project_sites_filters_below_twenty_and_caches_pairs():
    counts = {
        ("chr1", 1): SiteCount(alt=1, callable=19, p_alt_derived=1.0),
        ("chr1", 2): SiteCount(alt=1, callable=20, p_alt_derived=1.0),
        ("chr1", 3): SiteCount(alt=1, callable=21, p_alt_derived=1.0),
        ("chr1", 4): SiteCount(alt=1, callable=21, p_alt_derived=1.0),
    }
    rows, projections, endpoints = project_sites(counts)
    assert ("chr1", 1) not in rows
    assert projections.shape == (2, 19)
    assert rows[("chr1", 3)] == rows[("chr1", 4)]
    assert projections[rows[("chr1", 2)]][0] == pytest.approx(1)
    assert endpoints[rows[("chr1", 2)]] == pytest.approx(0)
    assert projections[rows[("chr1", 3)]][0] == pytest.approx(20 / 21)
    assert endpoints[rows[("chr1", 3)]] == pytest.approx(1 / 21)


def test_sites_are_not_renormalized_after_endpoint_removal():
    counts = {
        ("chr1", 1): SiteCount(alt=1, callable=21, p_alt_derived=1.0),
        ("chr1", 2): SiteCount(alt=10, callable=20, p_alt_derived=1.0),
    }
    rows, projections, endpoints = project_sites(counts)
    coordinates = [("chr1", 1), ("chr1", 2)]
    raw_counts, endpoint, eligible = accumulate_spectrum(
        coordinates, rows, projections, endpoints
    )
    assert eligible == 2
    assert raw_counts.sum() == pytest.approx(20 / 21 + 1)
    assert endpoint == pytest.approx(1 / 21)
    assert raw_counts.sum() + endpoint == pytest.approx(eligible)
    raw, normalized = normalized_spectrum(raw_counts)
    assert normalized.sum() == pytest.approx(1)


def test_bernoulli_q_projection_is_hard_reproducible_and_has_expected_rate():
    counts = {
        ("chr1", i): SiteCount(alt=4, callable=20, p_alt_derived=0.3)
        for i in range(1, 10_001)
    }
    rows1, projected1, _ = project_sites_bernoulli_q(counts, seed=17)
    rows2, projected2, _ = project_sites_bernoulli_q(
        dict(reversed(counts.items())), seed=17
    )
    alt_derived = np.array([
        projected1[rows1[coordinate]][3] == 1 for coordinate in counts
    ])
    assert alt_derived.mean() == pytest.approx(0.3, abs=0.015)
    for coordinate in counts:
        np.testing.assert_array_equal(
            projected1[rows1[coordinate]], projected2[rows2[coordinate]]
        )


def test_accumulation_is_order_invariant_and_counts_repeats():
    counts = {
        ("chr1", 1): SiteCount(alt=3, callable=25, p_alt_derived=1.0),
        ("chr1", 2): SiteCount(alt=9, callable=25, p_alt_derived=1.0),
    }
    rows, projections, endpoints = project_sites(counts)
    forward = accumulate_spectrum([("chr1", 1), ("chr1", 2)], rows, projections, endpoints)
    reverse = accumulate_spectrum([("chr1", 2), ("chr1", 1)], rows, projections, endpoints)
    assert np.allclose(forward[0], reverse[0])
    repeated = accumulate_spectrum(
        [("chr1", 1), ("chr1", 1)], rows, projections, endpoints
    )
    assert repeated[2] == 2
    assert np.allclose(
        repeated[0],
        2 * accumulate_spectrum([("chr1", 1)], rows, projections, endpoints)[0],
    )


def test_zero_retained_mass_fails():
    counts = {("chr1", 1): SiteCount(alt=1, callable=10, p_alt_derived=1.0)}
    rows, projections, endpoints = project_sites(counts)
    raw_counts, _, eligible = accumulate_spectrum(
        [("chr1", 1)], rows, projections, endpoints
    )
    assert eligible == 0
    with pytest.raises(ValueError, match="zero retained mass"):
        normalized_spectrum(raw_counts)


# ------------------------------------------------------------------ the score


def test_phi_is_wasserstein_distance_and_symmetric():
    a = np.zeros(19)
    b = np.zeros(19)
    a[0], b[-1] = 1, 1
    result = phi_sfs(a, b)
    assert result.value == pytest.approx(18 / 20)
    assert result.mean_daf_difference == pytest.approx(-18 / 20)
    np.testing.assert_allclose(result.cdf_a, np.ones(19))
    np.testing.assert_allclose(result.cdf_b[:-1], np.zeros(18))
    np.testing.assert_allclose(result.cdf_residual, result.cdf_a - result.cdf_b)
    np.testing.assert_allclose(result.bin_residual, a - b)
    assert phi_sfs(a, a).value == pytest.approx(0)
    assert phi_sfs(b, a).value == pytest.approx(result.value)
    assert phi_sfs(b, a).mean_daf_difference == pytest.approx(
        -result.mean_daf_difference
    )


def test_phi_matches_cdf_area_on_an_irregular_grid():
    a = np.array([0.5, 0.0, 0.5])
    b = np.array([0.0, 1.0, 0.0])
    daf = np.array([0.1, 0.4, 0.9])
    # The CDF difference is +0.5 over [0.1, 0.4), then -0.5 over
    # [0.4, 0.9), so the total shaded area is 0.4.
    result = phi_sfs(a, b, daf=daf)
    assert result.value == pytest.approx(0.4)
    assert result.mean_daf_difference == pytest.approx(0.1)


def test_phi_random_spectra_are_symmetric_and_nonnegative():
    rng = np.random.default_rng(1)
    for _ in range(5):
        a = rng.random(19)
        a /= a.sum()
        b = rng.random(19)
        b /= b.sum()
        value = phi_sfs(a, b).value
        assert value == pytest.approx(phi_sfs(b, a).value)
        assert value >= 0


@pytest.mark.parametrize(
    ("a", "b", "daf", "message"),
    [
        (np.ones((1, 19)) / 19, np.ones(19) / 19, None, "one-dimensional"),
        (np.ones(18) / 18, np.ones(19) / 19, None, "equal length"),
        (np.r_[np.nan, np.ones(18)], np.ones(19) / 19, None, "finite"),
        (
            np.r_[-0.1, np.repeat(1.1 / 18, 18)],
            np.ones(19) / 19,
            None,
            "nonnegative",
        ),
        (np.full(19, 0.1), np.full(19, 1 / 19), None, "normalized"),
        (
            np.ones(3) / 3,
            np.ones(3) / 3,
            np.array([0.1, 0.1, 0.9]),
            "strictly increasing",
        ),
    ],
)
def test_phi_rejects_invalid_input(a, b, daf, message):
    with pytest.raises(ValueError, match=message):
        phi_sfs(a, b, daf=daf)


def test_phi_requires_explicit_grid_for_nondefault_spectra():
    with pytest.raises(ValueError, match="daf is required"):
        phi_sfs(np.ones(3) / 3, np.ones(3) / 3)


def test_phi_agrees_with_scipy_wasserstein_distance():
    """`phi_sfs` is exactly a discrete Wasserstein-1 distance, not a lookalike."""
    scipy_stats = pytest.importorskip("scipy.stats")
    daf = RETAINED_BINS.astype(np.float64) / PROJECTION_SIZE
    rng = np.random.default_rng(7)
    for _ in range(10):
        a = rng.random(daf.size)
        a /= a.sum()
        b = rng.random(daf.size)
        b /= b.sum()
        expected = scipy_stats.wasserstein_distance(daf, daf, a, b)
        assert phi_sfs(a, b, daf=daf).value == pytest.approx(expected)


# ----------------------------------------------------------- null calibration


def test_calibrate_phi_uses_sample_sd_and_add_one_p_value():
    result = calibrate_phi(0.4, np.array([0.1, 0.2, 0.3]))
    assert result.observed == pytest.approx(0.4)
    np.testing.assert_allclose(result.null, [0.1, 0.2, 0.3])
    assert result.null_mean == pytest.approx(0.2)
    assert result.null_sd == pytest.approx(0.1)
    assert result.z_score == pytest.approx(2.0)
    assert result.exceedances == 0
    assert result.p_value == pytest.approx(1 / 4)
    np.testing.assert_allclose(result.null_z_scores, [-1.0, 0.0, 1.0], atol=1e-15)
    assert result.null_z_scores.mean() == pytest.approx(0.0, abs=1e-15)
    assert result.null_z_scores.std(ddof=1) == pytest.approx(1.0)

    permuted = calibrate_phi(0.4, np.array([0.3, 0.1, 0.2]))
    assert permuted.null_mean == pytest.approx(result.null_mean)
    assert permuted.null_sd == pytest.approx(result.null_sd)
    assert permuted.z_score == pytest.approx(result.z_score)
    assert permuted.p_value == pytest.approx(result.p_value)


def test_calibrate_phi_counts_ties_as_exceedances():
    result = calibrate_phi(0.2, np.array([0.1, 0.2, 0.3]))
    assert result.exceedances == 2
    assert result.p_value == pytest.approx(3 / 4)


@pytest.mark.parametrize(
    ("observed", "null", "message"),
    [
        (-0.1, np.array([0.1, 0.2]), "observed.*nonnegative"),
        (np.nan, np.array([0.1, 0.2]), "observed.*finite"),
        (0.1, np.array([[0.1, 0.2]]), "one-dimensional"),
        (0.1, np.array([0.1]), "at least two"),
        (0.1, np.array([0.1, np.inf]), "null.*finite"),
        (0.1, np.array([0.1, -0.2]), "null.*nonnegative"),
        (0.1, np.array([0.2, 0.2]), "zero sample standard deviation"),
    ],
)
def test_calibrate_phi_rejects_invalid_input(observed, null, message):
    with pytest.raises(ValueError, match=message):
        calibrate_phi(observed, null)


# -------------------------------------------------------------- the fixtures


def _write_bundle(
    root: Path,
    *,
    positions=None,
    row_indices=None,
    target_digest=None,
    a_type="TE",
    include_te_polarity=None,
    max_flipped_fraction=0.5,
    heterozygous="error",
):
    """Write a target and matched-control bundle that pass provenance checks.

    Defaults declare a valid TE target: `a_type` "TE", a `te_polarity` record
    at the required `max_flipped_fraction=0.5`, and a `vcf_eligibility` record
    carrying an `identity` dict. The matched-control metadata mirrors `a_type`
    and carries the identical `vcf_eligibility_identity` dict, as the matcher
    is contracted to copy it from the candidate report. `identity["vcf_sha256"]`
    is a placeholder here because no VCF exists yet at bundle-construction
    time; `_run` (or a direct call to `_sync_vcf_identity`) fills in the real
    digest once the fixture's VCF is written.
    """
    target = root / "target"
    matches = root / "matches"
    target.mkdir()
    matches.mkdir()

    if include_te_polarity is None:
        include_te_polarity = a_type == "TE"

    te_rows = np.array([0, 1], dtype=np.int64)
    cdf = np.array([0.5, 1.0], dtype=np.float64)
    ages = np.array([100.0, 200.0], dtype=np.float64)
    threshold = 12.5
    np.save(target / "te_chromosomes.npy", np.array(["chr1", "chr1"]), allow_pickle=False)
    np.save(target / "te_positions.npy", np.array([10, 20]), allow_pickle=False)
    np.save(target / "te_row_indices.npy", te_rows, allow_pickle=False)
    np.save(target / "target_cdf.npy", cdf, allow_pickle=False)
    np.save(target / "age_bins.npy", ages, allow_pickle=False)

    if positions is None:
        positions = np.array([[30, 40], [50, 60], [70, 80]])
    if row_indices is None:
        row_indices = np.array([[2, 3], [4, 5], [6, 7]], dtype=np.int64)
    np.save(matches / "positions.npy", np.asarray(positions), allow_pickle=False)
    np.save(matches / "row_indices.npy", np.asarray(row_indices), allow_pickle=False)
    np.save(matches / "chromosome_codes.npy",
            np.zeros(np.shape(positions), dtype=np.uint16), allow_pickle=False)
    np.save(matches / "chromosome_labels.npy", np.array(["chr1"]), allow_pickle=False)
    replicate_count = np.shape(positions)[0]
    np.save(matches / "replicate_id.npy", np.arange(replicate_count), allow_pickle=False)
    np.save(matches / "qc_pass.npy", np.ones(replicate_count, dtype=bool), allow_pickle=False)
    np.save(matches / "match_to_bootstrap_w1.npy",
            np.linspace(0.01, 0.03, replicate_count), allow_pickle=False)
    np.save(matches / "matching_error_ratio.npy",
            np.linspace(0.1, 0.3, replicate_count), allow_pickle=False)
    np.save(matches / "bootstrap_to_observed_w1.npy",
            np.linspace(0.1, 0.2, replicate_count), allow_pickle=False)
    np.save(matches / "bootstrap_seeds.npy",
            np.arange(100, 100 + replicate_count, dtype=np.uint64), allow_pickle=False)
    np.save(matches / "bootstrap_counts.npy",
            np.ones((replicate_count, te_rows.size), dtype=np.uint32), allow_pickle=False)
    unique_rows = np.unique(np.asarray(row_indices, dtype=np.int64))
    np.save(matches / "reuse_row_indices.npy", unique_rows, allow_pickle=False)
    np.save(matches / "reuse_counts.npy",
            np.ones(unique_rows.shape, dtype=np.uint16), allow_pickle=False)

    digest = _sha256_arrays(
        te_rows, cdf, ages, np.asarray([threshold], dtype=np.float64)
    )
    identity = {
        "vcf_sha256": "0" * 64,   # placeholder; _sync_vcf_identity fills this in
        "heterozygous": heterozygous,
        "min_callable": 20,
        "store_content_sha256": "store",
        "row_indices_sha256": "rows",
        "snp_row_indices_sha256": "snp_rows",
        "p_alt_derived_sha256": "polarity",
    }
    vcf_eligibility = {"mask": str((root / "eligibility").resolve()), "identity": identity}
    target_metadata = {
        "source_store_content_sha256": "store",
        "source_catalog_sha256": "catalog",
        "wasserstein_threshold_generations": threshold,
        "a_type": a_type,
        "vcf_eligibility": vcf_eligibility,
    }
    if include_te_polarity:
        target_metadata["te_polarity"] = {"max_flipped_fraction": max_flipped_fraction}
    (target / "metadata.json").write_text(json.dumps(target_metadata))
    (matches / "metadata.json").write_text(json.dumps({
        "schema_version": "bootstrap-target-matches-v1",
        "source_store_content_sha256": "store",
        "source_catalog_sha256": "catalog",
        "complete": True,
        "target_digest": target_digest if target_digest is not None else digest,
        "phi_sfs_selection_blind": True,
        "maximum_control_reuse": 1,
        "config": {"disjoint_replicates": True},
        "a_type": a_type,
        "vcf_eligibility_identity": dict(identity),
    }))
    return target, matches


def _sync_vcf_identity(target: Path, matches: Path, vcf: Path) -> None:
    """Patch the identity's `vcf_sha256` to match this fixture's actual VCF.

    `_write_bundle` writes metadata before the VCF exists, so the digest it
    records is a placeholder. `_run` calls this automatically so the ordinary
    end-to-end tests get a consistent, passing identity without each test
    having to know about it. A test that wants to exercise the digest-mismatch
    check calls `_run(..., sync_identity=False)` after pinning the identity to
    a different VCF's digest itself.
    """
    digest = hashlib.sha256(Path(vcf).read_bytes()).hexdigest()
    target_meta_path = target / "metadata.json"
    target_meta = json.loads(target_meta_path.read_text())
    target_meta["vcf_eligibility"]["identity"]["vcf_sha256"] = digest
    target_meta_path.write_text(json.dumps(target_meta))

    matches_meta_path = matches / "metadata.json"
    matches_meta = json.loads(matches_meta_path.read_text())
    matches_meta["vcf_eligibility_identity"]["vcf_sha256"] = digest
    matches_meta_path.write_text(json.dumps(matches_meta))


def _record(position: int, derived: int, *, callable_count: int = 20, info: str = "."):
    """One biallelic haploid record with an exact derived count."""
    calls = (
        ["1"] * derived
        + ["0"] * (callable_count - derived)
        + ["."] * (20 - callable_count)
    )
    return f"chr1\t{position}\t.\tA\tG\t.\tPASS\t{info}\tGT\t" + "\t".join(calls)


def _vcf_text(info: str = "."):
    records = [
        _record(10, 4, info=info),
        _record(20, 8, info=info),
        _record(30, 4, info=info),
        _record(40, 12, info=info),
        _record(50, 8, info=info),
        _record(60, 12, info=info),
        _record(70, 1, info=info),
        _record(80, 19, info=info),
    ]
    header = (
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        + "\t".join(f"s{index}" for index in range(20))
        + "\n"
    )
    return header + "\n".join(records) + "\n"


def _ancestral_table(tmp_path, positions=(10, 20, 30, 40, 50, 60, 70, 80)):
    """A table asserting REF (`A`) is ancestral at every site, unanimously.

    Every record in `_vcf_text` is `A`/`G`, so an ancestral call of `A` makes the
    ALT allele derived with weight exactly 1. That reproduces the semantics these
    tests were written against, when polarity came from the REF column, so their
    hand-calculated expectations still hold. Tests that need a genuinely
    uncertain site build their own table.
    """
    store = tmp_path / "fake_store"
    store.mkdir(exist_ok=True)
    np.save(store / "positions.npy", np.asarray(positions, dtype=np.float64))
    (store / "metadata.json").write_text(json.dumps({
        "chromosomes": [{"chrom": "chr1", "length": 1000, "offset": 0}],
        "content_sha256": "store",
    }), encoding="utf-8")

    table = tmp_path / "ancestral"
    table.mkdir(exist_ok=True)
    counts = np.zeros((len(positions), 4), dtype=np.uint16)
    counts[:, 0] = 75                      # column 0 is "A"
    np.save(table / "ancestral_counts.npy", counts)
    np.save(table / "present_draw_count.npy",
            np.full(len(positions), 75, dtype=np.uint16))
    (table / "metadata.json").write_text(json.dumps({
        "schema_version": "ancestral-state-counts-v1",
        "bases": ["A", "C", "G", "T"],
        "store": str(store),
        "store_content_sha256": "store",   # matches the bundle fixture's digest
        "store_rows": len(positions),
        "complete": True,
    }), encoding="utf-8")
    return table


def _run(target, matches, vcf, output, *extra, sync_identity=True, reference=0):
    """Run the CLI with a floor of 2 nulls.

    Hand-calculated tests fix B0 at replicate 0 through `reference`; pass
    `reference=None` to exercise the default seeded draw of B0.
    """
    if sync_identity:
        _sync_vcf_identity(target, matches, Path(vcf))
    argv = [
        "--target", str(target), "--matches", str(matches),
        "--vcf", str(vcf), "--output", str(output),
        "--min-null-replicates", "2", *extra,
    ]
    if reference is not None and "--reference-replicate" not in extra:
        argv += ["--reference-replicate", str(reference)]
    if "--ancestral-table" not in argv:
        argv += ["--ancestral-table", str(_ancestral_table(Path(output).parent))]
    return main(argv)


# ------------------------------------------------------------- end to end


def test_end_to_end_matches_hand_calculation(tmp_path):
    """Every site has n = 20, so each projects to a point mass at its own k.

    A is k = 4 and 8. B0 is k = 4 and 12, so half the mass moves four
    DAF bins and observed Phi is 0.5 * (4 / 20) = 0.1. The two null sets
    give distances 0.1 and 0.25 from B0.
    """
    target, matches = _write_bundle(tmp_path)
    vcf = tmp_path / "sites.vcf"
    output = tmp_path / "phi"
    vcf.write_text(_vcf_text())
    assert _run(target, matches, vcf, output) == 0

    bins = np.load(output / "bins.npy")
    assert bins.tolist() == list(range(1, 20))

    a = np.load(output / "a_normalized_sfs.npy")
    assert a[3] == pytest.approx(0.5)   # bin 4
    assert a[7] == pytest.approx(0.5)   # bin 8
    assert a.sum() == pytest.approx(1)

    b = np.load(output / "b_normalized_sfs.npy")
    assert b[0][3] == pytest.approx(0.5)
    assert b[0][11] == pytest.approx(0.5)   # bin 12
    assert b[1][7] == pytest.approx(0.5)
    assert b[1][11] == pytest.approx(0.5)

    assert np.load(output / "observed_phi_sfs.npy").item() == pytest.approx(0.1)
    assert np.load(output / "null_phi_sfs.npy").tolist() == pytest.approx([0.1, 0.25])
    assert np.load(output / "reference_replicate_id.npy").item() == 0
    assert np.load(output / "null_replicate_id.npy").tolist() == [1, 2]

    residual = np.load(output / "observed_bin_residual.npy")
    assert residual[7] == pytest.approx(0.5)
    assert residual[11] == pytest.approx(-0.5)

    summary = (output / "summary.csv").read_text().splitlines()
    assert len(summary) == 2
    comparisons = (output / "comparisons.csv").read_text().splitlines()
    assert len(comparisons) == 4


def test_end_to_end_metadata_and_diagnostics(tmp_path):
    target, matches = _write_bundle(tmp_path)
    vcf = tmp_path / "sites.vcf"
    output = tmp_path / "phi"
    vcf.write_text(_vcf_text())
    assert _run(target, matches, vcf, output) == 0

    metadata = json.loads((output / "metadata.json").read_text())
    assert metadata["complete"] is True
    assert metadata["schema_version"] == "phi-sfs-wasserstein-v1"
    assert metadata["accepted_null_replicates"] == 2
    assert metadata["a_eligible_sites"] == 2
    assert metadata["equal_eligible_site_count"] == 2
    assert metadata["a_type"] == "TE"
    assert metadata["b_type"] == "SNP"
    assert metadata["reference_replicate_id"] == 0
    assert metadata["maximum_control_reuse"] == 1
    assert metadata["maximum_overlap_with_reference"] == 0
    # Every eligible site has n = 20, so no mass reaches bins 0 or 20.
    assert metadata["a_retained_fraction"] == pytest.approx(1)
    assert metadata["a_endpoint_fraction"] == pytest.approx(0)
    assert metadata["distinct_projections"] == 5
    assert metadata["software"]["name"] == "PhiTE"
    assert metadata["creation_command"]
    assert metadata["creation_time_utc"]
    assert metadata["numpy_version"]
    assert len(metadata["vcf_sha256"]) == 64
    assert metadata["target_digest"] == json.loads(
        (matches / "metadata.json").read_text()
    )["target_digest"]
    assert metadata["vcf_eligibility_identity"] == json.loads(
        (target / "metadata.json").read_text()
    )["vcf_eligibility"]["identity"]

    summary_header, summary_values = (
        line.split(",") for line in (output / "summary.csv").read_text().splitlines()
    )
    summary = dict(zip(summary_header, summary_values))
    assert float(summary["observed_phi_sfs"]) == pytest.approx(0.1)
    assert float(summary["p_value"]) == pytest.approx(1.0)


def test_comparisons_reuse_columns_are_present_and_equal_one(tmp_path):
    """A valid disjoint bundle has every control used exactly once, everywhere."""
    target, matches = _write_bundle(tmp_path)
    vcf = tmp_path / "sites.vcf"
    vcf.write_text(_vcf_text())
    output = tmp_path / "phi"
    assert _run(target, matches, vcf, output) == 0

    with (output / "comparisons.csv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert rows[0]["role"] == "observed"
    assert rows[0]["left_max_control_reuse"] == ""
    assert rows[0]["right_max_control_reuse"] == "1"
    for row in rows[1:]:
        assert row["role"] == "null"
        assert row["left_max_control_reuse"] == "1"
        assert row["right_max_control_reuse"] == "1"


def test_vcf_sha256_matches_a_direct_digest(tmp_path):
    target, matches = _write_bundle(tmp_path)
    vcf = tmp_path / "sites.vcf"
    vcf.write_text(_vcf_text())
    assert _run(target, matches, vcf, tmp_path / "phi") == 0
    metadata = json.loads((tmp_path / "phi" / "metadata.json").read_text())
    assert metadata["vcf_sha256"] == hashlib.sha256(vcf.read_bytes()).hexdigest()


def test_compressed_input_is_read_and_hashed(tmp_path):
    target, matches = _write_bundle(tmp_path)
    vcf = tmp_path / "sites.vcf.bgz"
    vcf.write_bytes(gzip.compress(_vcf_text().encode()))
    assert _run(target, matches, vcf, tmp_path / "phi") == 0
    metadata = json.loads((tmp_path / "phi" / "metadata.json").read_text())
    assert metadata["vcf_sha256"] == hashlib.sha256(vcf.read_bytes()).hexdigest()
    assert np.load(tmp_path / "phi" / "observed_phi_sfs.npy").item() == pytest.approx(0.1)



def test_heterozygous_calls_fail_by_default(tmp_path):
    target, matches = _write_bundle(tmp_path)
    vcf = tmp_path / "sites.vcf"
    vcf.write_text(_vcf_text().replace("\t1\t", "\t0/1\t", 1))
    with pytest.raises(ValueError, match="heterozygous"):
        _run(target, matches, vcf, tmp_path / "phi")


def test_heterozygous_missing_policy_drops_the_individual(tmp_path):
    # The eligibility mask's recorded policy must agree with --heterozygous.
    target, matches = _write_bundle(tmp_path, heterozygous="missing")
    vcf = tmp_path / "sites.vcf"
    # Site 10 loses one derived individual: k = 3 among n = 19, so it is dropped.
    vcf.write_text(_vcf_text().replace("\t1\t", "\t0/1\t", 1))
    with pytest.raises(ValueError, match="shared eligibility mask"):
        _run(target, matches, vcf, tmp_path / "phi",
             "--heterozygous", "missing")



def test_missing_site_is_reported(tmp_path):
    target, matches = _write_bundle(tmp_path)
    vcf = tmp_path / "sites.vcf"
    vcf.write_text(_vcf_text().replace(_record(30, 4) + "\n", ""))
    with pytest.raises(ValueError, match="absent from the VCF"):
        _run(target, matches, vcf, tmp_path / "phi")


def test_existing_output_is_never_overwritten(tmp_path):
    target, matches = _write_bundle(tmp_path)
    vcf = tmp_path / "sites.vcf"
    vcf.write_text(_vcf_text())
    output = tmp_path / "phi"
    output.mkdir()
    with pytest.raises(FileExistsError):
        _run(target, matches, vcf, output)


def test_calculate_is_byte_reproducible(tmp_path):
    target, matches = _write_bundle(tmp_path)
    vcf = tmp_path / "sites.vcf"
    vcf.write_text(_vcf_text())
    out1 = tmp_path / "phi1"
    out2 = tmp_path / "phi2"
    assert _run(target, matches, vcf, out1) == 0
    assert _run(target, matches, vcf, out2) == 0

    names1 = sorted(path.name for path in out1.glob("*.npy"))
    names2 = sorted(path.name for path in out2.glob("*.npy"))
    assert names1 == names2
    for name in names1:
        assert (out1 / name).read_bytes() == (out2 / name).read_bytes(), name


def test_injected_publish_failure_leaves_no_output_or_staging(tmp_path, monkeypatch):
    target, matches = _write_bundle(tmp_path)
    vcf = tmp_path / "sites.vcf"
    vcf.write_text(_vcf_text())
    output = tmp_path / "phi"

    def failing_dump(*args, **kwargs):
        raise RuntimeError("injected failure")

    # metadata.json is written last, after every .npy file and both CSVs, so
    # this exercises cleanup of a staging directory that already holds output.
    monkeypatch.setattr(phi_sfs_module.json, "dump", failing_dump)

    with pytest.raises(RuntimeError, match="injected failure"):
        _run(target, matches, vcf, output)

    assert not output.exists()
    assert list(output.parent.glob(f".{output.name}.tmp.*")) == []


def test_end_to_end_snp_a_uses_posterior_polarity(tmp_path):
    target, matches = _write_bundle(tmp_path, a_type="SNP")
    vcf = tmp_path / "sites.vcf"
    vcf.write_text(_vcf_text())
    table = _ancestral_table(tmp_path)
    counts = np.load(table / "ancestral_counts.npy")
    counts[:, 0] = 100
    counts[:2, 0] = 25
    counts[:2, 2] = 75
    np.save(table / "ancestral_counts.npy", counts)
    np.save(
        table / "present_draw_count.npy",
        np.full(counts.shape[0], 100, dtype=np.uint16),
    )

    output = tmp_path / "phi"
    assert _run(
        target, matches, vcf, output,
        "-A", "SNP", "-B", "SNP", "--ancestral-table", str(table),
    ) == 0
    a = np.load(output / "a_normalized_sfs.npy")
    assert a[3] == pytest.approx(0.125)    # q * site 10 at bin 4
    assert a[15] == pytest.approx(0.375)  # (1-q) * site 10 at bin 16
    assert a[7] == pytest.approx(0.125)   # q * site 20 at bin 8
    assert a[11] == pytest.approx(0.375)  # (1-q) * site 20 at bin 12
    metadata = json.loads((output / "metadata.json").read_text())
    assert metadata["a_type"] == "SNP"
    assert metadata["b_type"] == "SNP"
    assert metadata["te_sites_polarized"] == 0
    header, values = (
        line.split(",") for line in (output / "summary.csv").read_text().splitlines()
    )
    summary = dict(zip(header, values))
    assert summary["a_type"] == "SNP"
    assert summary["b_type"] == "SNP"


def test_asymmetric_null_uses_hard_snp_left_and_mixture_reference(tmp_path):
    target, matches = _write_bundle(tmp_path, a_type="SNP")
    vcf = tmp_path / "sites.vcf"
    vcf.write_text(_vcf_text())
    table = _ancestral_table(tmp_path)
    counts = np.load(table / "ancestral_counts.npy")
    counts[:, 0] = 25
    counts[:, 2] = 75
    np.save(table / "ancestral_counts.npy", counts)
    np.save(
        table / "present_draw_count.npy",
        np.full(counts.shape[0], 100, dtype=np.uint16),
    )

    output = tmp_path / "phi"
    assert _run(
        target, matches, vcf, output,
        "-A", "SNP", "-B", "SNP", "--ancestral-table", str(table),
        "--asymmetric-polarity-null", "--polarity-imputation-seed", "17",
    ) == 0
    a = np.load(output / "a_normalized_sfs.npy")
    b_mix = np.load(output / "b_normalized_sfs.npy")
    b_hard = np.load(output / "b_bernoulli_q_normalized_sfs.npy")
    assert np.count_nonzero(a) == 2
    assert not np.array_equal(b_mix, b_hard)
    assert np.load(output / "observed_phi_sfs.npy").item() == pytest.approx(
        phi_sfs(a, b_mix[0]).value
    )
    assert np.load(output / "null_phi_sfs.npy").tolist() == pytest.approx([
        phi_sfs(b_hard[1], b_mix[0]).value,
        phi_sfs(b_hard[2], b_mix[0]).value,
    ])
    metadata = json.loads((output / "metadata.json").read_text())
    assert metadata["asymmetric_polarity_null"] is True
    assert metadata["polarity_imputation_seed"] == 17
    assert metadata["a_polarity_rule"] == "coordinate-keyed hard Bernoulli(q) orientation"
    assert metadata["null_left_polarity_rule"] == (
        "coordinate-keyed hard Bernoulli(q) orientation"
    )


def test_explicit_reference_id_controls_b0_and_null_identities(tmp_path):
    target, matches = _write_bundle(tmp_path)
    vcf = tmp_path / "sites.vcf"
    vcf.write_text(_vcf_text())
    output = tmp_path / "phi"
    assert _run(
        target, matches, vcf, output, "--reference-replicate", "1"
    ) == 0
    assert np.load(output / "reference_replicate_id.npy").item() == 1
    assert np.load(output / "null_replicate_id.npy").tolist() == [0, 2]
    assert np.load(output / "observed_phi_sfs.npy").item() == pytest.approx(0.2)


# ------------------------------------------------------------- bundle checks


def test_matches_built_for_another_target_are_rejected(tmp_path):
    target, matches = _write_bundle(tmp_path, target_digest="0" * 64)
    vcf = tmp_path / "sites.vcf"
    vcf.write_text(_vcf_text())
    with pytest.raises(ValueError, match="built for a different target"):
        _run(target, matches, vcf, tmp_path / "phi")


def test_incomplete_matched_bundle_is_rejected(tmp_path):
    target, matches = _write_bundle(tmp_path)
    metadata = json.loads((matches / "metadata.json").read_text())
    metadata["complete"] = False
    (matches / "metadata.json").write_text(json.dumps(metadata))
    vcf = tmp_path / "sites.vcf"
    vcf.write_text(_vcf_text())
    with pytest.raises(ValueError, match="not marked complete"):
        _run(target, matches, vcf, tmp_path / "phi")


def test_nonbootstrap_match_schema_is_rejected(tmp_path):
    target, matches = _write_bundle(tmp_path)
    metadata = json.loads((matches / "metadata.json").read_text())
    metadata["schema_version"] = "swap-age-matched-controls-v1"
    (matches / "metadata.json").write_text(json.dumps(metadata))
    vcf = tmp_path / "sites.vcf"
    vcf.write_text(_vcf_text())
    with pytest.raises(ValueError, match="requires bootstrap-target-matches-v1"):
        _run(target, matches, vcf, tmp_path / "phi")


def test_nondisjoint_bundle_is_rejected(tmp_path):
    target, matches = _write_bundle(tmp_path)
    metadata = json.loads((matches / "metadata.json").read_text())
    metadata["config"]["disjoint_replicates"] = False
    (matches / "metadata.json").write_text(json.dumps(metadata))
    vcf = tmp_path / "sites.vcf"
    vcf.write_text(_vcf_text())
    with pytest.raises(ValueError, match="--disjoint-replicates"):
        _run(target, matches, vcf, tmp_path / "phi")


def test_global_control_reuse_is_rejected(tmp_path):
    target, matches = _write_bundle(tmp_path)
    rows = np.load(matches / "row_indices.npy")
    rows[1, 0] = rows[0, 0]
    np.save(matches / "row_indices.npy", rows)
    vcf = tmp_path / "sites.vcf"
    vcf.write_text(_vcf_text())
    with pytest.raises(ValueError, match="control used more than once"):
        _run(target, matches, vcf, tmp_path / "phi")


def test_focal_rows_are_excluded_from_controls(tmp_path):
    target, matches = _write_bundle(tmp_path)
    rows = np.load(matches / "row_indices.npy")
    rows[0, 0] = 0
    np.save(matches / "row_indices.npy", rows)
    vcf = tmp_path / "sites.vcf"
    vcf.write_text(_vcf_text())
    with pytest.raises(ValueError, match="exclude every row in focal set A"):
        _run(target, matches, vcf, tmp_path / "phi")


def test_reference_must_pass_qc(tmp_path):
    target, matches = _write_bundle(tmp_path)
    qc = np.load(matches / "qc_pass.npy")
    qc[0] = False
    np.save(matches / "qc_pass.npy", qc)
    vcf = tmp_path / "sites.vcf"
    vcf.write_text(_vcf_text())
    with pytest.raises(ValueError, match="reference replicate ID 0 failed"):
        _run(target, matches, vcf, tmp_path / "phi")


def test_minimum_null_count_is_enforced(tmp_path):
    target, matches = _write_bundle(tmp_path)
    vcf = tmp_path / "sites.vcf"
    vcf.write_text(_vcf_text())
    with pytest.raises(ValueError, match="only 2 QC-passing null replicates"):
        _run(
            target, matches, vcf, tmp_path / "phi",
            "--min-null-replicates", "3",
        )


def _spare_set_bundle(tmp_path):
    """Four matched sets: B0 plus up to three nulls."""
    target, matches = _write_bundle(
        tmp_path,
        positions=np.array([[30, 40], [50, 60], [70, 80], [90, 100]]),
        row_indices=np.array([[2, 3], [4, 5], [6, 7], [8, 9]], dtype=np.int64),
    )
    vcf = tmp_path / "sites.vcf"
    vcf.write_text(_vcf_text() + _record(90, 2) + "\n" + _record(100, 6) + "\n")
    table = _ancestral_table(
        tmp_path, positions=(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
    )
    return target, matches, vcf, table


def test_every_qc_passing_set_other_than_b0_is_a_null(tmp_path):
    target, matches, vcf, table = _spare_set_bundle(tmp_path)
    output = tmp_path / "phi"
    assert _run(target, matches, vcf, output, "--ancestral-table", str(table)) == 0
    assert np.load(output / "null_replicate_id.npy").tolist() == [1, 2, 3]
    assert np.load(output / "b_replicate_id.npy").tolist() == [0, 1, 2, 3]
    metadata = json.loads((output / "metadata.json").read_text())
    assert metadata["accepted_null_replicates"] == 3
    assert metadata["minimum_null_replicates"] == 2
    assert metadata["matched_sets_published"] == 4
    assert metadata["matched_sets_failing_qc"] == 0
    assert metadata["null_selection_rule"] == "every QC-passing non-reference set"
    summary = dict(zip(*(
        line.split(",") for line in (output / "summary.csv").read_text().splitlines()
    )))
    assert summary["null_replicates_r"] == "3"
    assert float(summary["minimum_attainable_p"]) == pytest.approx(1 / 4)


def test_a_failed_set_simply_reduces_r(tmp_path):
    target, matches, vcf, table = _spare_set_bundle(tmp_path)
    qc = np.load(matches / "qc_pass.npy")
    qc[1] = False
    np.save(matches / "qc_pass.npy", qc)
    output = tmp_path / "phi"
    assert _run(target, matches, vcf, output, "--ancestral-table", str(table)) == 0
    assert np.load(output / "null_replicate_id.npy").tolist() == [2, 3]
    assert np.load(output / "b_replicate_id.npy").tolist() == [0, 2, 3]
    metadata = json.loads((output / "metadata.json").read_text())
    assert metadata["matched_sets_failing_qc"] == 1
    assert metadata["accepted_null_replicates"] == 2


def test_default_b0_is_a_reproducible_seeded_draw(tmp_path):
    target, matches, vcf, table = _spare_set_bundle(tmp_path)
    first, second = tmp_path / "phi1", tmp_path / "phi2"
    for output in (first, second):
        assert _run(
            target, matches, vcf, output, "--ancestral-table", str(table),
            reference=None,
        ) == 0
    reference = np.load(first / "reference_replicate_id.npy").item()
    assert np.load(second / "reference_replicate_id.npy").item() == reference
    assert reference in (0, 1, 2, 3)
    nulls = np.load(first / "null_replicate_id.npy").tolist()
    assert sorted(nulls + [reference]) == [0, 1, 2, 3]
    metadata = json.loads((first / "metadata.json").read_text())
    assert metadata["reference_seed"] == 1002
    assert metadata["reference_selection_rule"].startswith("uniform draw from QC-passing sets")


def test_reference_seed_changes_the_drawn_b0(tmp_path):
    """Some seed among a handful must pick a different B0 from four sets."""
    target, matches, vcf, table = _spare_set_bundle(tmp_path)
    drawn = set()
    for seed in range(8):
        output = tmp_path / f"phi{seed}"
        assert _run(
            target, matches, vcf, output, "--ancestral-table", str(table),
            "--reference-seed", str(seed), reference=None,
        ) == 0
        drawn.add(np.load(output / "reference_replicate_id.npy").item())
    assert len(drawn) > 1


def test_a_drawn_b0_never_fails_qc(tmp_path):
    target, matches, vcf, table = _spare_set_bundle(tmp_path)
    qc = np.load(matches / "qc_pass.npy")
    qc[0] = False
    np.save(matches / "qc_pass.npy", qc)
    for seed in range(6):
        output = tmp_path / f"phi{seed}"
        assert _run(
            target, matches, vcf, output, "--ancestral-table", str(table),
            "--reference-seed", str(seed), reference=None,
        ) == 0
        assert np.load(output / "reference_replicate_id.npy").item() in (1, 2, 3)


def _sensitivity_bundle(tmp_path):
    """Five matched sets (ids 0-4).

    A sensitivity rerun's own null set can then reach a set the primary run
    never touches, exercising the requirement that every set any analysis
    uses -- not just the primary accepted set -- gets scanned.
    """
    target, matches = _write_bundle(
        tmp_path,
        positions=np.array(
            [[30, 40], [50, 60], [70, 80], [90, 100], [110, 120]]
        ),
        row_indices=np.array(
            [[2, 3], [4, 5], [6, 7], [8, 9], [10, 11]], dtype=np.int64
        ),
    )
    vcf = tmp_path / "sites.vcf"
    vcf.write_text(
        _vcf_text()
        + _record(90, 2) + "\n" + _record(100, 6) + "\n"
        + _record(110, 5) + "\n" + _record(120, 9) + "\n"
    )
    table = _ancestral_table(
        tmp_path,
        positions=(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120),
    )
    return target, matches, vcf, table


def test_reference_sensitivity_selection_and_primary_arrays_unchanged(tmp_path):
    """With reference id 2, sensitivity N=2 draws two other QC-passing sets.

    The alternatives come from the same seeded permutation that would draw B0,
    so they are prespecified, distinct, never the primary reference, and
    reproducible. Requesting sensitivity must not perturb any primary array.
    """
    target, matches, vcf, table = _sensitivity_bundle(tmp_path)

    out_plain = tmp_path / "phi_plain"
    assert _run(
        target, matches, vcf, out_plain,
        "--reference-replicate", "2", "--ancestral-table", str(table),
    ) == 0

    out_sensitivity = tmp_path / "phi_sensitivity"
    assert _run(
        target, matches, vcf, out_sensitivity,
        "--reference-replicate", "2", "--reference-sensitivity", "2",
        "--ancestral-table", str(table),
    ) == 0

    for name in (
        "observed_phi_sfs.npy", "null_phi_sfs.npy", "null_z_scores.npy",
        "a_normalized_sfs.npy", "b_normalized_sfs.npy", "b_raw_sfs.npy",
        "reference_replicate_id.npy", "null_replicate_id.npy",
        "b_replicate_id.npy",
    ):
        np.testing.assert_array_equal(
            np.load(out_plain / name), np.load(out_sensitivity / name)
        )

    def _summary(directory):
        header, values = (
            line.split(",") for line in (directory / "summary.csv").read_text().splitlines()
        )
        return dict(zip(header, values))

    plain_summary = _summary(out_plain)
    sensitivity_summary = _summary(out_sensitivity)
    for key in ("observed_phi_sfs", "null_mean", "null_sample_sd", "z_score", "p_value"):
        assert plain_summary[key] == sensitivity_summary[key]

    ids = np.load(out_sensitivity / "sensitivity_reference_ids.npy")
    assert len(set(ids.tolist())) == 2
    assert set(ids.tolist()) <= {0, 1, 3, 4}
    rerun = tmp_path / "phi_sensitivity_rerun"
    assert _run(
        target, matches, vcf, rerun,
        "--reference-replicate", "2", "--reference-sensitivity", "2",
        "--ancestral-table", str(table),
    ) == 0
    np.testing.assert_array_equal(np.load(rerun / "sensitivity_reference_ids.npy"), ids)
    z_scores = np.load(out_sensitivity / "sensitivity_z_scores.npy")
    p_values = np.load(out_sensitivity / "sensitivity_p_values.npy")
    observed = np.load(out_sensitivity / "sensitivity_observed_phi_sfs.npy")
    assert z_scores.shape == (2,)
    assert p_values.shape == (2,)
    assert observed.shape == (2,)

    metadata = json.loads((out_sensitivity / "metadata.json").read_text())
    assert metadata["reference_sensitivity_run"] is True
    assert metadata["reference_sensitivity_n"] == 2
    assert metadata["reference_sensitivity_reference_ids"] == ids.tolist()

    plain_metadata = json.loads((out_plain / "metadata.json").read_text())
    assert plain_metadata["reference_sensitivity_run"] is False
    assert plain_metadata["reference_sensitivity_n"] == 0

    header, values = (
        line.split(",")
        for line in (out_sensitivity / "summary.csv").read_text().splitlines()
    )
    summary = dict(zip(header, values))
    assert summary["reference_sensitivity_n"] == "2"
    assert float(summary["reference_sensitivity_z_min"]) == pytest.approx(
        float(min(z_scores))
    )
    assert float(summary["reference_sensitivity_z_max"]) == pytest.approx(
        float(max(z_scores))
    )
    assert float(summary["reference_sensitivity_p_min"]) == pytest.approx(
        float(min(p_values))
    )
    assert float(summary["reference_sensitivity_p_max"]) == pytest.approx(
        float(max(p_values))
    )


def test_reference_sensitivity_fails_if_too_few_alternatives(tmp_path):
    """Five sets leave four alternatives to reference id 3, not five."""
    target, matches, vcf, table = _sensitivity_bundle(tmp_path)
    with pytest.raises(ValueError, match="only 4 other QC-passing sets exist"):
        _run(
            target, matches, vcf, tmp_path / "phi",
            "--reference-replicate", "3", "--reference-sensitivity", "5",
            "--ancestral-table", str(table),
        )


def test_reference_sensitivity_defaults_to_zero_and_publishes_empty_arrays(tmp_path):
    target, matches = _write_bundle(tmp_path)
    vcf = tmp_path / "sites.vcf"
    vcf.write_text(_vcf_text())
    output = tmp_path / "phi"
    assert _run(target, matches, vcf, output) == 0

    header, values = (
        line.split(",") for line in (output / "summary.csv").read_text().splitlines()
    )
    summary = dict(zip(header, values))
    assert summary["reference_sensitivity_n"] == "0"
    assert summary["reference_sensitivity_z_min"] == ""
    assert summary["reference_sensitivity_z_max"] == ""
    assert summary["reference_sensitivity_p_min"] == ""
    assert summary["reference_sensitivity_p_max"] == ""

    for name in (
        "sensitivity_reference_ids.npy", "sensitivity_observed_phi_sfs.npy",
        "sensitivity_z_scores.npy", "sensitivity_p_values.npy",
    ):
        assert np.load(output / name).shape == (0,)

    metadata = json.loads((output / "metadata.json").read_text())
    assert metadata["reference_sensitivity_run"] is False
    assert metadata["reference_sensitivity_n"] == 0
    assert metadata["reference_sensitivity_reference_ids"] == []


def test_equal_eligible_site_count_is_enforced(tmp_path):
    target, matches = _write_bundle(tmp_path)
    vcf = tmp_path / "sites.vcf"
    vcf.write_text(_vcf_text().replace(_record(50, 8), _record(50, 5, callable_count=10)))
    with pytest.raises(ValueError, match="B replicate 1 retains 1 of 2"):
        _run(target, matches, vcf, tmp_path / "phi")


def test_absent_store_provenance_is_rejected(tmp_path):
    target, matches = _write_bundle(tmp_path)
    metadata = json.loads((matches / "metadata.json").read_text())
    del metadata["source_catalog_sha256"]
    (matches / "metadata.json").write_text(json.dumps(metadata))
    vcf = tmp_path / "sites.vcf"
    vcf.write_text(_vcf_text())
    with pytest.raises(ValueError, match="must both record"):
        _run(target, matches, vcf, tmp_path / "phi")


def test_null_store_provenance_blocks_ancestral_binding(tmp_path):
    """Null store digests are tolerated everywhere except table authentication.

    A dense store records no content digest, and `target_digest` alone binds the
    target to its matched sets, so null provenance was previously harmless. The
    ancestral table has no other identity to bind against: accepting null would
    accept a table from any store, which is the silent wrong-polarity failure the
    check exists to prevent. So it is fatal here and only here.
    """
    target, matches = _write_bundle(tmp_path)
    for path in (target / "metadata.json", matches / "metadata.json"):
        metadata = json.loads(path.read_text())
        metadata["source_store_content_sha256"] = None
        metadata["source_catalog_sha256"] = None
        path.write_text(json.dumps(metadata))
    vcf = tmp_path / "sites.vcf"
    vcf.write_text(_vcf_text())
    with pytest.raises(ValueError, match="no store content digest"):
        _run(target, matches, vcf, tmp_path / "phi")


def test_ancestral_table_from_another_store_is_rejected(tmp_path):
    """The digest is the identity, so a table from elsewhere must not be used.

    Row coordinates overlap between stores, so a foreign table yields a
    complete, plausible result with the wrong control polarity rather than an
    error. Only the digest catches it.
    """
    target, matches = _write_bundle(tmp_path)
    vcf = tmp_path / "sites.vcf"
    vcf.write_text(_vcf_text())
    other = tmp_path / "other"
    other.mkdir()
    table = _ancestral_table(other)
    metadata = json.loads((table / "metadata.json").read_text())
    metadata["store_content_sha256"] = "a-different-store"
    (table / "metadata.json").write_text(json.dumps(metadata))
    with pytest.raises(ValueError, match="different interval store"):
        _run(target, matches, vcf, tmp_path / "phi",
             "--ancestral-table", str(table))


def test_ancestral_table_source_path_replacement_is_rejected(tmp_path):
    """The coordinate catalog loaded by path must still be the authenticated one."""
    target, matches = _write_bundle(tmp_path)
    vcf = tmp_path / "sites.vcf"
    vcf.write_text(_vcf_text())
    table = _ancestral_table(tmp_path)
    metadata_path = Path(json.loads(
        (table / "metadata.json").read_text()
    )["store"]) / "metadata.json"
    metadata = json.loads(metadata_path.read_text())
    metadata["content_sha256"] = "replacement-store"
    metadata_path.write_text(json.dumps(metadata))
    with pytest.raises(ValueError, match="recorded source-store path"):
        _run(target, matches, vcf, tmp_path / "phi",
             "--ancestral-table", str(table))


def test_null_store_provenance_still_requires_a_matching_target(tmp_path):
    target, matches = _write_bundle(tmp_path, target_digest="0" * 64)
    for path in (target / "metadata.json", matches / "metadata.json"):
        metadata = json.loads(path.read_text())
        metadata["source_store_content_sha256"] = None
        metadata["source_catalog_sha256"] = None
        path.write_text(json.dumps(metadata))
    vcf = tmp_path / "sites.vcf"
    vcf.write_text(_vcf_text())
    with pytest.raises(ValueError, match="built for a different target"):
        _run(target, matches, vcf, tmp_path / "phi")


def test_mismatched_store_provenance_is_rejected(tmp_path):
    target, matches = _write_bundle(tmp_path)
    metadata = json.loads((matches / "metadata.json").read_text())
    metadata["source_store_content_sha256"] = "other"
    (matches / "metadata.json").write_text(json.dumps(metadata))
    vcf = tmp_path / "sites.vcf"
    vcf.write_text(_vcf_text())
    with pytest.raises(ValueError, match="values differ"):
        _run(target, matches, vcf, tmp_path / "phi")


def test_endpoint_fraction_is_reported_when_n_exceeds_twenty(tmp_path):
    """With 21 callable individuals, real mass falls in bins 0 and 20."""
    target, matches = _write_bundle(tmp_path)
    vcf = tmp_path / "sites.vcf"
    header = (
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        + "\t".join(f"s{index}" for index in range(21))
        + "\n"
    )
    records = []
    for position, derived in (
        (10, 1), (20, 8), (30, 4), (40, 12),
        (50, 8), (60, 12), (70, 1), (80, 20),
    ):
        calls = ["1"] * derived + ["0"] * (21 - derived)
        records.append(
            f"chr1\t{position}\t.\tA\tG\t.\tPASS\t.\tGT\t" + "\t".join(calls)
        )
    vcf.write_text(header + "\n".join(records) + "\n")
    assert _run(target, matches, vcf, tmp_path / "phi") == 0

    metadata = json.loads((tmp_path / "phi" / "metadata.json").read_text())
    assert metadata["a_eligible_sites"] == 2
    assert metadata["a_endpoint_fraction"] > 0
    # Site 10 is k = 1 of n = 21, so it loses exactly 1/21 to bin 0; site 20
    # loses nothing. The two fractions must together account for every site.
    assert metadata["a_endpoint_fraction"] == pytest.approx((1 / 21) / 2)
    assert (
        metadata["a_retained_fraction"] + metadata["a_endpoint_fraction"]
        == pytest.approx(1)
    )


def test_non_integer_row_indices_are_rejected(tmp_path):
    target, matches = _write_bundle(tmp_path)
    np.save(matches / "row_indices.npy", np.array([[2.9, 3.0], [3.0, 4.0]]),
            allow_pickle=False)
    vcf = tmp_path / "sites.vcf"
    vcf.write_text(_vcf_text())
    with pytest.raises(ValueError, match="must be an integer array"):
        _run(target, matches, vcf, tmp_path / "phi")


def test_negative_row_indices_are_rejected(tmp_path):
    target, matches = _write_bundle(
        tmp_path, row_indices=np.array([[-1, 3], [4, 5], [6, 7]])
    )
    vcf = tmp_path / "sites.vcf"
    vcf.write_text(_vcf_text())
    with pytest.raises(ValueError, match="non-negative"):
        _run(target, matches, vcf, tmp_path / "phi")


def test_non_integer_positions_are_rejected(tmp_path):
    target, matches = _write_bundle(tmp_path)
    np.save(
        matches / "positions.npy",
        np.array([[30.0, 40.0], [50.0, 60.0], [70.0, 80.0]]),
        allow_pickle=False,
    )
    vcf = tmp_path / "sites.vcf"
    vcf.write_text(_vcf_text())
    with pytest.raises(ValueError, match="must be an integer array"):
        _run(target, matches, vcf, tmp_path / "phi")


def test_misaligned_row_indices_are_rejected(tmp_path):
    target, matches = _write_bundle(tmp_path)
    np.save(matches / "row_indices.npy", np.array([[2, 3, 9], [3, 4, 9]]),
            allow_pickle=False)
    vcf = tmp_path / "sites.vcf"
    vcf.write_text(_vcf_text())
    with pytest.raises(ValueError, match="do not align"):
        _run(target, matches, vcf, tmp_path / "phi")


def test_duplicate_controls_within_a_set_are_rejected(tmp_path):
    target, matches = _write_bundle(
        tmp_path, row_indices=np.array([[2, 2], [4, 5], [6, 7]])
    )
    vcf = tmp_path / "sites.vcf"
    vcf.write_text(_vcf_text())
    with pytest.raises(ValueError, match="duplicate control rows"):
        _run(target, matches, vcf, tmp_path / "phi")


def test_table_calling_alt_ancestral_reverses_polarization(tmp_path):
    """Flipping the table mirrors the control spectra and leaves the TE one alone.

    Every record is `A`/`G`. For a control SNP, a table naming `A` ancestral
    makes the derived count its ALT count, and naming `G` ancestral makes it
    `n - alt`, so the two control spectra must be mirror images. TE sites are
    polarized by biology and never consult the table, so the TE spectrum must be
    identical across the two runs -- which is the asymmetry the two-arm design
    exists to produce.
    """
    target, matches = _write_bundle(tmp_path)
    vcf = tmp_path / "v.vcf"
    vcf.write_text(_vcf_text(), encoding="utf-8")

    fwd = tmp_path / "fwd"; fwd.mkdir()
    forward = _ancestral_table(fwd)
    out_f = tmp_path / "out_forward"
    assert _run(target, matches, vcf, out_f, "--ancestral-table", str(forward)) == 0

    rev = tmp_path / "rev"; rev.mkdir()
    reversed_table = _ancestral_table(rev)
    counts = np.load(reversed_table / "ancestral_counts.npy")
    counts[:, 0] = 0                       # not A
    counts[:, 2] = 75                      # G, the ALT allele, is ancestral
    np.save(reversed_table / "ancestral_counts.npy", counts)
    out_r = tmp_path / "out_reversed"
    assert _run(target, matches, vcf, out_r, "--ancestral-table",
                str(reversed_table)) == 0

    # The TE spectrum is polarized by biology and must be untouched by the table.
    np.testing.assert_allclose(
        np.load(out_f / "a_normalized_sfs.npy"),
        np.load(out_r / "a_normalized_sfs.npy"), atol=1e-12,
    )
    # The control spectra are polarized by the table, so they must mirror.
    a = np.load(out_f / "b_normalized_sfs.npy")
    b = np.load(out_r / "b_normalized_sfs.npy")
    np.testing.assert_allclose(a, b[:, ::-1], atol=1e-12)


def test_site_the_table_cannot_orient_is_rejected(tmp_path):
    """A draw naming a third base cannot orient the site, and none may remain.

    Conditioning on draws that named one of the two observed alleles is what
    makes the weight meaningful. If no draw named either, there is no weight to
    compute and guessing one would invent polarity the ARG never supplied.
    """
    target, matches = _write_bundle(tmp_path)
    vcf = tmp_path / "v.vcf"
    vcf.write_text(_vcf_text(), encoding="utf-8")
    bad = tmp_path / "bad"; bad.mkdir()
    table = _ancestral_table(bad)
    counts = np.load(table / "ancestral_counts.npy")
    counts[:] = 0
    counts[:, 1] = 75                      # every draw says C, neither A nor G
    np.save(table / "ancestral_counts.npy", counts)
    with pytest.raises(ValueError, match="cannot be polarized"):
        _run(target, matches, vcf, tmp_path / "out", "--ancestral-table", str(table))


# --------------------------------------- target authority and VCF identity


def test_snp_target_run_with_a_type_te_is_rejected(tmp_path):
    """A SNP target run with `-A TE` must fail, not silently mispolarize."""
    target, matches = _write_bundle(tmp_path, a_type="SNP")
    vcf = tmp_path / "sites.vcf"
    vcf.write_text(_vcf_text())
    with pytest.raises(ValueError, match="disagrees with target metadata a_type"):
        _run(target, matches, vcf, tmp_path / "phi", "-A", "TE")


def test_te_target_without_te_polarity_is_rejected(tmp_path):
    target, matches = _write_bundle(tmp_path, include_te_polarity=False)
    vcf = tmp_path / "sites.vcf"
    vcf.write_text(_vcf_text())
    with pytest.raises(ValueError, match="--te-polarity-mask"):
        _run(target, matches, vcf, tmp_path / "phi")


def test_te_target_with_wrong_max_flipped_fraction_is_rejected(tmp_path):
    target, matches = _write_bundle(tmp_path, max_flipped_fraction=0.6)
    vcf = tmp_path / "sites.vcf"
    vcf.write_text(_vcf_text())
    with pytest.raises(ValueError, match="max_flipped_fraction=0.5"):
        _run(target, matches, vcf, tmp_path / "phi")


def test_target_and_match_vcf_eligibility_identity_mismatch_is_rejected(tmp_path):
    target, matches = _write_bundle(tmp_path)
    metadata = json.loads((target / "metadata.json").read_text())
    metadata["vcf_eligibility"]["identity"]["row_indices_sha256"] = "different-rows"
    (target / "metadata.json").write_text(json.dumps(metadata))
    vcf = tmp_path / "sites.vcf"
    vcf.write_text(_vcf_text())
    with pytest.raises(ValueError, match="vcf_eligibility_identity values differ"):
        _run(target, matches, vcf, tmp_path / "phi")


def test_heterozygous_policy_mismatch_is_rejected(tmp_path):
    """The eligibility mask's policy must match `--heterozygous`, not just A/B."""
    target, matches = _write_bundle(tmp_path, heterozygous="missing")
    vcf = tmp_path / "sites.vcf"
    vcf.write_text(_vcf_text())
    with pytest.raises(
        ValueError, match="disagrees with the eligibility mask's policy"
    ):
        _run(target, matches, vcf, tmp_path / "phi")


def test_vcf_digest_mismatch_is_rejected(tmp_path):
    """A VCF with different bytes but identical sites must still be rejected.

    Every requested site is still callable and identical, so nothing else in
    `calculate()` would catch this; only the recorded `vcf_sha256` can.
    """
    target, matches = _write_bundle(tmp_path)
    vcf = tmp_path / "sites.vcf"
    vcf.write_text(_vcf_text())
    _sync_vcf_identity(target, matches, vcf)   # pin identity to this VCF's digest

    other_vcf = tmp_path / "other.vcf"
    other_vcf.write_text(
        _vcf_text().replace(
            "##fileformat=VCFv4.2\n", "##fileformat=VCFv4.2\n##note=repacked\n"
        )
    )
    with pytest.raises(
        ValueError, match="the VCF differs from the one that defined eligibility"
    ):
        _run(target, matches, other_vcf, tmp_path / "phi", sync_identity=False)
