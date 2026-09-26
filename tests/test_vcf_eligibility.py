import gzip
import hashlib
import json
from types import SimpleNamespace

import numpy as np
import pytest

from normalize_tes import build_candidate_rows
from normalize_tes import te_age_target
from normalize_tes.vcf_eligibility import (
    EligibilityResult,
    eligibility_identity,
    load_ancestral_table,
    load_eligibility,
    load_eligible_rows,
    publish,
    scan_vcf,
)


def _store():
    return SimpleNamespace(
        positions=np.array([1.0, 2.0, 3.0, 4.0, 5.0]),
        eligible=np.array([True, True, True, True, False]),
        chromosomes=({"chrom": "chr1", "offset": 0, "length": 10},),
        metadata={"content_sha256": "content", "catalog_sha256": "catalog"},
    )


def _record(position, genotypes, *, alt="G", fmt="GT"):
    return "\t".join([
        "chr1", str(position), ".", "A", alt, ".", "PASS", ".", fmt,
        *genotypes,
    ])


def _write_vcf(path, records, *, compressed=False):
    text = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsamples\n"
    text += "\n".join(records) + "\n"
    if compressed:
        path.write_bytes(gzip.compress(text.encode()))
    else:
        path.write_text(text, encoding="utf-8")
    return path


def _ancestry(*, unoriented=()):
    counts = np.zeros((5, 4), dtype=np.uint16)
    counts[:, 0] = 3  # REF=A ancestral => ALT derived, q=1
    counts[list(unoriented)] = 0
    return counts, counts.sum(axis=1, dtype=np.uint16)


def _write_ancestral_table(path, *, digest="content"):
    path.mkdir()
    counts, present = _ancestry()
    np.save(path / "ancestral_counts.npy", counts)
    np.save(path / "present_draw_count.npy", present)
    (path / "metadata.json").write_text(json.dumps({
        "schema_version": "ancestral-state-counts-v1",
        "bases": ["A", "C", "G", "T"],
        "store_content_sha256": digest,
        "store_rows": 5,
        "complete": True,
    }), encoding="utf-8")
    return path


def test_scan_enforces_callability_and_nonpolarity_rules(tmp_path):
    called20 = ["0"] * 19 + ["1"]
    called19 = ["0"] * 19 + ["."]
    vcf = _write_vcf(tmp_path / "sites.vcf", [
        _record(1, called20),
        _record(2, called19),
        _record(3, called20, alt="G,T"),
        _record(4, ["0"] * 19 + ["0/1"]),
        _record(5, called20),  # store-ineligible independently of callability
        _record(9, called20),  # absent from the store catalog
    ])
    counts, present = _ancestry()
    result = scan_vcf(
        vcf, _store(), min_callable=20,
        ancestral_counts=counts, present_draw_count=present,
    )
    assert result.rows.tolist() == [0]
    assert result.alt_counts.tolist() == [1]
    assert result.callable_counts.tolist() == [20]
    assert result.report["excluded_by_reason"] == {
        "insufficient_callability": 1,
        "multiallelic_record": 1,
        "heterozygous_genotype": 1,
        "store_ineligible": 1,
    }
    assert result.snp_rows.tolist() == [0]
    assert result.p_alt_derived.tolist() == [1.0]
    assert result.report["polarity_filtering"] == "none"
    assert result.report["applicable_set_roles"] == ["A", "B"]
    assert result.report["applicable_variant_types"] == ["TE", "SNP"]


def test_ancestral_table_is_authenticated_against_store(tmp_path):
    table = _write_ancestral_table(tmp_path / "ancestral")
    counts, present, metadata = load_ancestral_table(table, _store())
    assert counts.shape == (5, 4)
    assert present.shape == (5,)
    assert metadata["store_content_sha256"] == "content"

    other = _store()
    other.metadata["content_sha256"] = "other"
    with pytest.raises(SystemExit, match="authenticate"):
        load_ancestral_table(table, other)


def test_heterozygous_missing_policy_matches_callability_semantics(tmp_path):
    vcf = _write_vcf(
        tmp_path / "sites.vcf.gz",
        [_record(1, ["0"] * 20 + ["0/1"])],
        compressed=True,
    )
    counts, present = _ancestry()
    result = scan_vcf(
        vcf, _store(), min_callable=20, heterozygous="missing",
        ancestral_counts=counts, present_draw_count=present,
    )
    assert result.rows.tolist() == [0]
    assert result.callable_counts.tolist() == [20]
    assert result.report["vcf_sha256"] == hashlib.sha256(vcf.read_bytes()).hexdigest()


def test_published_mask_authenticates_store_and_projection_size(tmp_path):
    store = _store()
    result = EligibilityResult(
        rows=np.array([0, 2], dtype=np.int64),
        alt_counts=np.array([1, 2], dtype=np.uint32),
        callable_counts=np.array([20, 21], dtype=np.uint32),
        report={
            "schema_version": "vcf-eligibility-v1",
            "eligible_rows": 2,
            "min_callable": 20,
            "vcf_sha256": "vcf-fake-hash",
            "heterozygous": "error",
            "store_content_sha256": "content",
            "store_catalog_sha256": "catalog",
        },
        snp_rows=np.array([0, 2], dtype=np.int64),
        p_alt_derived=np.array([0.25, 0.75]),
    )
    mask = tmp_path / "mask"
    publish(mask, result, {})
    np.testing.assert_array_equal(
        load_eligible_rows(
            mask, store, variant_type="TE", expected_min_callable=20
        ), [0, 2]
    )
    loaded_snp = load_eligibility(
        mask, store, variant_type="SNP", expected_min_callable=20
    )
    np.testing.assert_array_equal(loaded_snp.rows, [0, 2])
    np.testing.assert_allclose(loaded_snp.p_alt_derived, [0.25, 0.75])
    with pytest.raises(SystemExit, match="expected 21"):
        load_eligible_rows(mask, store, variant_type="TE", expected_min_callable=21)
    other = _store()
    other.metadata["content_sha256"] = "other"
    with pytest.raises(SystemExit, match="authenticate"):
        load_eligible_rows(mask, other, variant_type="TE")

    # Finding 4: the loaded object exposes the same content-addressed identity
    # a caller would get by hashing the published metadata directly, so two
    # consumers of the same mask are guaranteed to agree without comparing paths.
    published_metadata = json.loads((mask / "metadata.json").read_text())
    expected_identity = eligibility_identity(published_metadata)
    assert set(expected_identity) == {
        "vcf_sha256",
        "heterozygous",
        "min_callable",
        "store_content_sha256",
        "row_indices_sha256",
        "snp_row_indices_sha256",
        "p_alt_derived_sha256",
    }
    assert expected_identity["vcf_sha256"] == "vcf-fake-hash"
    assert expected_identity["heterozygous"] == "error"
    assert expected_identity["min_callable"] == 20
    assert expected_identity["store_content_sha256"] == "content"
    loaded_te = load_eligibility(mask, store, variant_type="TE", expected_min_callable=20)
    assert loaded_te.identity == expected_identity
    assert loaded_snp.identity == expected_identity


def test_eligibility_identity_requires_every_field():
    complete = {
        "vcf_sha256": "vcf-hash",
        "heterozygous": "missing",
        "min_callable": 20,
        "store_content_sha256": "store-hash",
        "array_sha256": {
            "row_indices": "row-hash",
            "alt_counts": "alt-hash",
            "callable_counts": "callable-hash",
            "snp_row_indices": "snp-row-hash",
            "p_alt_derived": "q-hash",
        },
    }
    assert eligibility_identity(complete) == {
        "vcf_sha256": "vcf-hash",
        "heterozygous": "missing",
        "min_callable": 20,
        "store_content_sha256": "store-hash",
        "row_indices_sha256": "row-hash",
        "snp_row_indices_sha256": "snp-row-hash",
        "p_alt_derived_sha256": "q-hash",
    }
    for key in ("vcf_sha256", "heterozygous", "min_callable", "store_content_sha256"):
        incomplete = {k: v for k, v in complete.items() if k != key}
        with pytest.raises(ValueError, match=key):
            eligibility_identity(incomplete)
    for key in ("row_indices", "snp_row_indices", "p_alt_derived"):
        incomplete_arrays = {**complete, "array_sha256": {
            k: v for k, v in complete["array_sha256"].items() if k != key
        }}
        with pytest.raises(ValueError, match=f"{key}_sha256"):
            eligibility_identity(incomplete_arrays)
    with pytest.raises(ValueError, match="row_indices_sha256"):
        eligibility_identity({k: v for k, v in complete.items() if k != "array_sha256"})


def test_candidate_builder_excludes_snp_a_target_from_vcf_eligible_b_pool(
    tmp_path, monkeypatch
):
    store = _store()
    mask = tmp_path / "mask"
    result = EligibilityResult(
        rows=np.array([0, 2, 3], dtype=np.int64),
        alt_counts=np.zeros(3, dtype=np.uint32),
        callable_counts=np.full(3, 20, dtype=np.uint32),
        report={
            "schema_version": "vcf-eligibility-v1",
            "eligible_rows": 3,
            "min_callable": 20,
            "vcf_sha256": "vcf-fake-hash",
            "heterozygous": "error",
            "store_content_sha256": "content",
            "store_catalog_sha256": "catalog",
        },
        snp_rows=np.array([0, 2, 3], dtype=np.int64),
        p_alt_derived=np.array([0.1, 0.5, 0.9]),
    )
    publish(mask, result, {})
    monkeypatch.setattr(build_candidate_rows, "open_snp_age_store", lambda _: store)
    monkeypatch.setattr(build_candidate_rows, "store_schema", lambda _: "test-store")
    monkeypatch.setattr(build_candidate_rows, "software_provenance", lambda: {})
    monkeypatch.setattr(
        build_candidate_rows,
        "_resolve_lists",
        # The exclusion list is the focal A set. Its variant type is
        # irrelevant: an A=SNP row must be removed from the B=SNP universe in
        # exactly the same way as an A=TE row.
        lambda _store, paths, minimum_fraction, kind: (
            np.array([0, 1, 2, 3], dtype=np.int64), []
        ) if kind == "inclusion" else (np.array([3], dtype=np.int64), []),
    )
    output = tmp_path / "candidates.npy"
    exclude = tmp_path / "exclude.txt"
    include = tmp_path / "include.txt"
    assert build_candidate_rows.main([
        "--store", str(tmp_path / "store"),
        "--exclude-positions", str(exclude),
        "--include-positions", str(include),
        "--vcf-eligibility", str(mask),
        "--output", str(output),
    ]) == 0
    np.testing.assert_array_equal(np.load(output), [0, 2])
    report = json.loads(output.with_suffix(".npy.json").read_text())
    assert report["universe_rows"] == 3
    assert report["vcf_eligibility"] == str(mask.resolve())
    # Finding 4: the report binds the mask by content, not only by path.
    published_metadata = json.loads((mask / "metadata.json").read_text())
    assert report["vcf_eligibility_identity"] == eligibility_identity(published_metadata)


def test_candidate_report_identity_is_null_without_vcf_eligibility(tmp_path, monkeypatch):
    store = _store()
    monkeypatch.setattr(build_candidate_rows, "open_snp_age_store", lambda _: store)
    monkeypatch.setattr(build_candidate_rows, "store_schema", lambda _: "test-store")
    monkeypatch.setattr(build_candidate_rows, "software_provenance", lambda: {})
    monkeypatch.setattr(
        build_candidate_rows,
        "_resolve_lists",
        lambda _store, paths, minimum_fraction, kind: (
            np.array([3], dtype=np.int64), []
        ),
    )
    output = tmp_path / "candidates.npy"
    exclude = tmp_path / "exclude.txt"
    assert build_candidate_rows.main([
        "--store", str(tmp_path / "store"),
        "--exclude-positions", str(exclude),
        "--output", str(output),
    ]) == 0
    report = json.loads(output.with_suffix(".npy.json").read_text())
    assert report["vcf_eligibility"] is None
    assert report["vcf_eligibility_identity"] is None


def test_loader_rejects_tampered_count_arrays(tmp_path):
    store = _store()
    result = EligibilityResult(
        rows=np.array([0], dtype=np.int64),
        alt_counts=np.array([1], dtype=np.uint32),
        callable_counts=np.array([20], dtype=np.uint32),
        report={
            "schema_version": "vcf-eligibility-v1",
            "eligible_rows": 1,
            "min_callable": 20,
            "store_content_sha256": "content",
            "store_catalog_sha256": "catalog",
        },
        snp_rows=np.array([0], dtype=np.int64),
        p_alt_derived=np.array([0.5]),
    )
    mask = tmp_path / "mask"
    publish(mask, result, {})
    np.save(mask / "callable_counts.npy", np.array([19], dtype=np.uint32))
    with pytest.raises(SystemExit, match="content digest"):
        load_eligible_rows(mask, store, variant_type="TE")


def test_target_eligibility_mask_preserves_a_site_order_and_alignment():
    target_rows = np.array([4, 1, 3, 2], dtype=np.int64)
    keep = te_age_target.target_eligibility_mask(
        target_rows, np.array([2, 3, 4], dtype=np.int64)
    )
    assert keep.tolist() == [True, False, True, True]
    # The same mask can be applied to target coordinates and an already
    # polarity-filtered draw matrix without changing their alignment.
    coordinates = np.array([40, 10, 30, 20])
    draw_mask = np.arange(12).reshape(4, 3)
    np.testing.assert_array_equal(coordinates[keep], [40, 30, 20])
    np.testing.assert_array_equal(draw_mask[keep], draw_mask[[0, 2, 3]])


def test_snp_a_requires_eligibility_artifact(monkeypatch):
    monkeypatch.setattr(te_age_target, "open_snp_age_store", lambda _: _store())
    monkeypatch.setattr(
        te_age_target,
        "load_native_position_list",
        lambda _: (np.array(["chr1"]), np.array([1])),
    )
    resolution = SimpleNamespace(
        included_global_positions=np.array([1.0]),
        included_chromosomes=np.array(["chr1"]),
        included_native_positions=np.array([1]),
        included_rows=np.array([0]),
    )
    monkeypatch.setattr(
        te_age_target, "resolve_native_position_requests", lambda *args, **kwargs: resolution
    )
    with pytest.raises(SystemExit, match="requires --vcf-eligibility"):
        te_age_target.main([
            "--store", "store", "--te-positions", "target", "--output", "out",
            "--a-type", "SNP",
        ])


def test_snp_subset_excludes_zero_orientation_but_keeps_full_q(tmp_path):
    called = ["0"] * 19 + ["1"]
    vcf = _write_vcf(tmp_path / "sites.vcf", [
        _record(1, called),
        _record(2, called),
        _record(3, called),
    ])
    counts = np.zeros((5, 4), dtype=np.uint16)
    counts[0, 0] = 1                         # q=1
    counts[1, 2] = 3                         # ALT=G ancestral, q=0
    counts[2, 0] = 1; counts[2, 2] = 3       # q=1/4
    present = counts.sum(axis=1, dtype=np.uint16)
    result = scan_vcf(
        vcf, _store(), ancestral_counts=counts, present_draw_count=present
    )
    assert result.rows.tolist() == [0, 1, 2]
    assert result.snp_rows.tolist() == [0, 1, 2]
    np.testing.assert_allclose(result.p_alt_derived, [1.0, 0.0, 0.25])

    counts[1] = 0
    present[1] = 0
    result = scan_vcf(
        vcf, _store(), ancestral_counts=counts, present_draw_count=present
    )
    assert result.rows.tolist() == [0, 1, 2]
    assert result.snp_rows.tolist() == [0, 2]
    assert result.report["excluded_by_reason"]["snp_no_usable_orientation"] == 1
