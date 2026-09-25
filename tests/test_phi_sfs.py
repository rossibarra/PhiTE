import gzip
import json
from pathlib import Path

import numpy as np
import pytest

from normalize_tes.phi_sfs import (
    SiteCount,
    accumulate_spectrum,
    calibrate_phi,
    hypergeometric_projection,
    main,
    normalized_spectrum,
    phi_sfs,
    project_sites,
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


def _write_bundle(root: Path, *, positions=None, row_indices=None, target_digest=None):
    """Write a target and matched-control bundle that pass provenance checks."""
    target = root / "target"
    matches = root / "matches"
    target.mkdir()
    matches.mkdir()

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
    np.save(matches / "reuse_counts.npy",
            np.ones(np.size(row_indices), dtype=np.uint16), allow_pickle=False)

    digest = _sha256_arrays(
        te_rows, cdf, ages, np.asarray([threshold], dtype=np.float64)
    )
    (target / "metadata.json").write_text(json.dumps({
        "source_store_content_sha256": "store",
        "source_catalog_sha256": "catalog",
        "wasserstein_threshold_generations": threshold,
    }))
    (matches / "metadata.json").write_text(json.dumps({
        "schema_version": "bootstrap-target-matches-v1",
        "source_store_content_sha256": "store",
        "source_catalog_sha256": "catalog",
        "complete": True,
        "target_digest": target_digest if target_digest is not None else digest,
        "phi_sfs_selection_blind": True,
        "maximum_control_reuse": 1,
        "config": {"disjoint_replicates": True},
    }))
    return target, matches


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


def _run(target, matches, vcf, output, *extra):
    argv = [
        "--target", str(target), "--matches", str(matches),
        "--vcf", str(vcf), "--output", str(output),
        "--null-replicates", "2", *extra,
    ]
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

    summary_header, summary_values = (
        line.split(",") for line in (output / "summary.csv").read_text().splitlines()
    )
    summary = dict(zip(summary_header, summary_values))
    assert float(summary["observed_phi_sfs"]) == pytest.approx(0.1)
    assert float(summary["p_value"]) == pytest.approx(1.0)


def test_vcf_sha256_matches_a_direct_digest(tmp_path):
    import hashlib

    target, matches = _write_bundle(tmp_path)
    vcf = tmp_path / "sites.vcf"
    vcf.write_text(_vcf_text())
    assert _run(target, matches, vcf, tmp_path / "phi") == 0
    metadata = json.loads((tmp_path / "phi" / "metadata.json").read_text())
    assert metadata["vcf_sha256"] == hashlib.sha256(vcf.read_bytes()).hexdigest()


def test_compressed_input_is_read_and_hashed(tmp_path):
    import hashlib

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
    target, matches = _write_bundle(tmp_path)
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


def test_end_to_end_snp_a_uses_posterior_polarity(tmp_path):
    target, matches = _write_bundle(tmp_path)
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
            "--null-replicates", "3",
        )


def _spare_set_bundle(tmp_path):
    """Four matched sets, so R = 2 leaves one spare set."""
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


def test_nulls_are_the_first_r_passing_sets_in_replicate_order(tmp_path):
    target, matches, vcf, table = _spare_set_bundle(tmp_path)
    output = tmp_path / "phi"
    assert _run(target, matches, vcf, output, "--ancestral-table", str(table)) == 0
    assert np.load(output / "null_replicate_id.npy").tolist() == [1, 2]
    assert np.load(output / "b_replicate_id.npy").tolist() == [0, 1, 2]
    metadata = json.loads((output / "metadata.json").read_text())
    assert metadata["requested_null_replicates"] == 2
    assert metadata["matched_sets_published"] == 4
    assert metadata["matched_sets_failing_qc"] == 0
    assert metadata["qc_passing_null_sets_available"] == 3
    assert metadata["unused_qc_passing_replicate_ids"] == [3]


def test_a_failed_set_is_replaced_by_the_next_spare(tmp_path):
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
    assert metadata["unused_qc_passing_replicate_ids"] == []


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
