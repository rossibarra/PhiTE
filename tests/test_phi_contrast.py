import csv
import json
from pathlib import Path

import numpy as np
import pytest

from normalize_tes.phi_contrast import (
    REQUIRED_RESULT_SCHEMA,
    SCHEMA_VERSION,
    benjamini_hochberg,
    calculate,
    contrast,
    derive_seed,
    holm,
    load_category,
    main,
    paired_contrast,
    parse_args,
)


# --------------------------------------------------------------- contrast()


def test_contrast_hand_calculated_and_ties_count_as_exceedances():
    # null1 is constant, so null1[pairing] - null2 is independent of the
    # permutation actually drawn: it is always 2 - null2 = [3, 2, 1, 0].
    z1, z2 = 4.0, 2.0  # delta_obs = 2.0
    null1 = np.array([2.0, 2.0, 2.0, 2.0])
    null2 = np.array([-1.0, 0.0, 1.0, 2.0])
    result = contrast(z1, z2, null1, null2, np.random.default_rng(0))
    assert result.delta_obs == pytest.approx(2.0)
    assert np.allclose(np.sort(result.null_delta), [0.0, 1.0, 2.0, 3.0])
    # abs(null_delta) >= abs(2.0): values 3 and 2 (the exact tie) qualify.
    assert result.exceedances == 2
    assert result.p_value == pytest.approx(3 / 5)
    assert set(result.pairing.tolist()) == {0, 1, 2, 3}


def test_contrast_minimum_p_is_one_over_r_plus_one():
    null1 = np.array([0.0, 0.0, 0.0])
    null2 = np.array([0.0, 0.1, -0.1])
    result = contrast(100.0, 0.0, null1, null2, np.random.default_rng(1))
    assert result.exceedances == 0
    assert result.p_value == pytest.approx(1 / 4)


def test_contrast_with_unequal_r_uses_min_r_pairs():
    """Five and three null replicates give three pairs and P over 3 + 1."""
    null1 = np.array([0.0, 0.1, 0.2, 0.3, 0.4])
    null2 = np.array([0.0, 0.0, 0.0])
    result = contrast(100.0, 0.0, null1, null2, np.random.default_rng(4))
    assert result.null_delta.shape == (3,)
    assert len(set(result.pairing.tolist())) == 3
    assert result.p_value == pytest.approx(1 / 4)
    swapped = contrast(0.0, 100.0, null2, null1, np.random.default_rng(4))
    assert swapped.null_delta.shape == (3,)
    assert swapped.p_value == pytest.approx(1 / 4)


def test_contrast_with_equal_r_is_a_permutation_of_the_first_side():
    null1 = np.array([0.5, -0.5, 1.5, 2.0])
    null2 = np.array([0.1, 0.2, 0.3, 0.4])
    result = contrast(1.0, 0.0, null1, null2, np.random.default_rng(9))
    expected = np.random.default_rng(9).permutation(4)
    np.testing.assert_array_equal(result.pairing, expected)
    np.testing.assert_allclose(result.null_delta, null1[expected] - null2)


def test_contrast_rejects_empty_and_nonfinite():
    with pytest.raises(ValueError, match="at least one"):
        contrast(1.0, 0.0, np.array([0.0]), np.array([]), np.random.default_rng(0))
    with pytest.raises(ValueError, match="finite"):
        contrast(float("nan"), 0.0, np.array([0.0]), np.array([0.0]), np.random.default_rng(0))
    with pytest.raises(ValueError, match="finite"):
        contrast(1.0, 0.0, np.array([np.inf]), np.array([0.0]), np.random.default_rng(0))


# ------------------------------------------------------------- derive_seed


def test_derive_seed_is_order_independent_and_deterministic():
    a = derive_seed(1002, "alpha", "beta", 0)
    b = derive_seed(1002, "beta", "alpha", 0)
    assert a == b
    # Different repeat or different base seed must (almost certainly) differ.
    assert derive_seed(1002, "alpha", "beta", 1) != a
    assert derive_seed(7, "alpha", "beta", 0) != a
    # Deterministic: recomputing gives the same value.
    assert derive_seed(1002, "alpha", "beta", 0) == a


# ---------------------------------------------------------- paired_contrast


def test_paired_contrast_deterministic_given_seed():
    null1 = np.array([0.1, -0.2, 0.3, 0.0])
    null2 = np.array([0.2, 0.1, -0.1, 0.05])
    first = paired_contrast("A", 1.5, null1, "B", 0.5, null2, seed=1002, repeat=0)
    second = paired_contrast("A", 1.5, null1, "B", 0.5, null2, seed=1002, repeat=0)
    assert first.delta_obs == second.delta_obs
    assert np.array_equal(first.null_delta, second.null_delta)
    assert first.p_value == second.p_value


def test_paired_contrast_swap_flips_delta_sign_but_keeps_p():
    null1 = np.array([0.1, -0.2, 0.3, 0.0, 0.4])
    null2 = np.array([0.2, 0.1, -0.1, 0.05, -0.3])
    forward = paired_contrast("cat1", 2.0, null1, "cat2", 0.5, null2, seed=1002, repeat=0)
    backward = paired_contrast("cat2", 0.5, null2, "cat1", 2.0, null1, seed=1002, repeat=0)
    assert forward.delta_obs == pytest.approx(1.5)
    assert backward.delta_obs == pytest.approx(-1.5)
    assert np.allclose(backward.null_delta, -forward.null_delta)
    assert backward.p_value == pytest.approx(forward.p_value)
    assert backward.exceedances == forward.exceedances


# ---------------------------------------------------------------- holm/BH


def test_holm_matches_hand_worked_example():
    # Textbook example: sorted p already ascending, so output order == input order.
    p = [0.01, 0.02, 0.03, 0.5]
    adjusted = holm(p)
    assert adjusted == pytest.approx([0.04, 0.06, 0.06, 0.5])


def test_holm_respects_input_order_when_unsorted():
    sorted_p = [0.01, 0.02, 0.03, 0.5]
    expected_sorted_adjusted = [0.04, 0.06, 0.06, 0.5]
    permutation = [2, 0, 3, 1]  # p[permutation] reorders the sorted list
    shuffled = [sorted_p[i] for i in permutation]
    expected = [expected_sorted_adjusted[i] for i in permutation]
    assert holm(shuffled) == pytest.approx(expected)


def test_benjamini_hochberg_matches_hand_worked_example():
    p = [0.01, 0.02, 0.03, 0.5]
    adjusted = benjamini_hochberg(p)
    assert adjusted == pytest.approx([0.04, 0.04, 0.04, 0.5])


def test_benjamini_hochberg_respects_input_order_when_unsorted():
    sorted_p = [0.01, 0.02, 0.03, 0.5]
    expected_sorted_adjusted = [0.04, 0.04, 0.04, 0.5]
    permutation = [3, 1, 0, 2]
    shuffled = [sorted_p[i] for i in permutation]
    expected = [expected_sorted_adjusted[i] for i in permutation]
    assert benjamini_hochberg(shuffled) == pytest.approx(expected)


def test_holm_and_bh_reject_out_of_range_p():
    with pytest.raises(ValueError, match="\\[0, 1\\]"):
        holm([0.1, 1.5])
    with pytest.raises(ValueError, match="\\[0, 1\\]"):
        benjamini_hochberg([-0.1, 0.5])


# ------------------------------------------------------------------ helpers


def _write_result(
    directory: Path,
    *,
    a_type: str = "TE",
    b_type: str = "SNP",
    z_score: float,
    null_z_scores: np.ndarray,
    target_digest: str = "deadbeef",
    complete: bool = True,
    schema_version: str = REQUIRED_RESULT_SCHEMA,
    null_polarity_design: str = "bernoulli-q-hard-vs-posterior-mixture",
    write_summary: bool = True,
) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    r = int(null_z_scores.size)
    np.save(directory / "null_z_scores.npy", null_z_scores, allow_pickle=False)
    metadata = {
        "schema_version": schema_version,
        "complete": complete,
        "a_type": a_type,
        "b_type": b_type,
        "accepted_null_replicates": r,
        "z_score": z_score,
        "target_digest": target_digest,
        "null_polarity_design": null_polarity_design,
    }
    (directory / "metadata.json").write_text(json.dumps(metadata), encoding="utf-8")
    if write_summary:
        with (directory / "summary.csv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=["z_score"])
            writer.writeheader()
            writer.writerow({"z_score": z_score})
    return directory


def _rng_null(seed: int, r: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    values = rng.normal(size=r)
    return (values - values.mean()) / values.std(ddof=1)


# ---------------------------------------------------------------- load_category


def test_load_category_rejects_wrong_schema(tmp_path):
    directory = _write_result(
        tmp_path / "cat", z_score=1.0, null_z_scores=_rng_null(1, 5),
        schema_version="phi-sfs-v1",
    )
    with pytest.raises(ValueError, match="unsupported schema_version"):
        load_category("cat", directory)


def test_load_category_rejects_incomplete(tmp_path):
    directory = _write_result(
        tmp_path / "cat", z_score=1.0, null_z_scores=_rng_null(1, 5), complete=False,
    )
    with pytest.raises(ValueError, match="not marked complete"):
        load_category("cat", directory)


def test_load_category_requires_null_polarity_design(tmp_path):
    directory = _write_result(
        tmp_path / "cat", z_score=1.0, null_z_scores=_rng_null(1, 5),
    )
    metadata = json.loads((directory / "metadata.json").read_text())
    del metadata["null_polarity_design"]
    (directory / "metadata.json").write_text(json.dumps(metadata), encoding="utf-8")
    with pytest.raises(ValueError, match="null_polarity_design"):
        load_category("cat", directory)


def test_load_category_rejects_mismatched_summary_z_score(tmp_path):
    directory = _write_result(
        tmp_path / "cat", z_score=1.0, null_z_scores=_rng_null(1, 5),
    )
    with (directory / "summary.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["z_score"])
        writer.writeheader()
        writer.writerow({"z_score": 99.0})
    with pytest.raises(ValueError, match="disagrees with"):
        load_category("cat", directory)


def test_load_category_rejects_null_array_length_mismatch(tmp_path):
    directory = _write_result(
        tmp_path / "cat", z_score=1.0, null_z_scores=_rng_null(1, 5), write_summary=False,
    )
    metadata = json.loads((directory / "metadata.json").read_text())
    metadata["accepted_null_replicates"] = 7
    (directory / "metadata.json").write_text(json.dumps(metadata), encoding="utf-8")
    with pytest.raises(ValueError, match="accepted_null_replicates"):
        load_category("cat", directory)


# ---------------------------------------------------------------- CLI / end-to-end


def _build_two_categories(tmp_path, r=20, seed_a=1, seed_b=2, z_a=2.5, z_b=0.3):
    a_dir = _write_result(
        tmp_path / "phi_sfs_A", a_type="TE", z_score=z_a, null_z_scores=_rng_null(seed_a, r),
    )
    b_dir = _write_result(
        tmp_path / "phi_sfs_B", a_type="TE", z_score=z_b, null_z_scores=_rng_null(seed_b, r),
    )
    return a_dir, b_dir


def test_cli_end_to_end_writes_expected_outputs(tmp_path):
    a_dir, b_dir = _build_two_categories(tmp_path)
    output = tmp_path / "contrast_out"
    argv = [
        "--result", f"A={a_dir}",
        "--result", f"B={b_dir}",
        "--output", str(output),
        "--pairing-repeats", "10",
    ]
    assert main(argv) == 0

    metadata = json.loads((output / "metadata.json").read_text())
    assert metadata["schema_version"] == SCHEMA_VERSION
    assert metadata["complete"] is True
    assert metadata["accepted_null_replicates"] == {"A": 20, "B": 20}
    assert metadata["null_polarity_design"] == (
        "bernoulli-q-hard-vs-posterior-mixture"
    )
    assert metadata["pairing_repeats"] == 10
    assert len(metadata["inputs"]) == 2

    with (output / "contrasts.csv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 1
    row = rows[0]
    # --result order is preserved, so the pair is reported as (A, B).
    assert (row["label1"], row["label2"]) == ("A", "B")
    assert float(row["z1"]) == pytest.approx(2.5)
    assert float(row["z2"]) == pytest.approx(0.3)
    assert float(row["delta_obs"]) == pytest.approx(2.2)
    assert 0.0 <= float(row["p"]) <= 1.0
    assert float(row["p_min"]) <= float(row["p"]) <= float(row["p_max"])
    assert float(row["holm_p"]) >= float(row["p"])
    assert float(row["bh_p"]) >= float(row["p"])

    null_file = output / row["null_contrasts_file"]
    assert null_file.exists()
    null_values = np.load(null_file, allow_pickle=False)
    assert null_values.shape == (20,)


def test_cli_contrasts_categories_with_different_r(tmp_path):
    a_dir = _write_result(tmp_path / "A", z_score=1.0, null_z_scores=_rng_null(1, 10))
    b_dir = _write_result(tmp_path / "B", z_score=1.0, null_z_scores=_rng_null(2, 20))
    output = tmp_path / "out"
    argv = ["--result", f"A={a_dir}", "--result", f"B={b_dir}", "--output", str(output)]
    assert main(argv) == 0
    header, values = (
        line.split(",") for line in (output / "contrasts.csv").read_text().splitlines()
    )
    row = dict(zip(header, values))
    assert (row["r1"], row["r2"], row["r"]) == ("10", "20", "10")
    assert float(row["p"]) >= 1 / 11
    null_file = output / row["null_contrasts_file"]
    assert np.load(null_file).shape == (10,)


def test_cli_rejects_wrong_schema(tmp_path):
    a_dir = _write_result(tmp_path / "A", z_score=1.0, null_z_scores=_rng_null(1, 10))
    b_dir = _write_result(
        tmp_path / "B", z_score=1.0, null_z_scores=_rng_null(2, 10),
        schema_version="phi-sfs-v1",
    )
    output = tmp_path / "out"
    argv = ["--result", f"A={a_dir}", "--result", f"B={b_dir}", "--output", str(output)]
    with pytest.raises(ValueError, match="unsupported schema_version"):
        main(argv)
    assert not output.exists()


def test_cli_rejects_incomplete_result(tmp_path):
    a_dir = _write_result(tmp_path / "A", z_score=1.0, null_z_scores=_rng_null(1, 10))
    b_dir = _write_result(
        tmp_path / "B", z_score=1.0, null_z_scores=_rng_null(2, 10), complete=False,
    )
    output = tmp_path / "out"
    argv = ["--result", f"A={a_dir}", "--result", f"B={b_dir}", "--output", str(output)]
    with pytest.raises(ValueError, match="not marked complete"):
        main(argv)
    assert not output.exists()


def test_cli_rejects_mixed_null_polarity_designs(tmp_path):
    a_dir = _write_result(
        tmp_path / "A", z_score=1.0, null_z_scores=_rng_null(1, 10),
    )
    b_dir = _write_result(
        tmp_path / "B", z_score=1.0, null_z_scores=_rng_null(2, 10),
        null_polarity_design="posterior-mixture-vs-posterior-mixture",
    )
    output = tmp_path / "out"
    argv = ["--result", f"A={a_dir}", "--result", f"B={b_dir}", "--output", str(output)]
    with pytest.raises(ValueError, match="share one null_polarity_design"):
        main(argv)
    assert not output.exists()


def test_cli_refuses_to_overwrite(tmp_path):
    a_dir, b_dir = _build_two_categories(tmp_path)
    output = tmp_path / "out"
    argv = ["--result", f"A={a_dir}", "--result", f"B={b_dir}", "--output", str(output)]
    assert main(argv) == 0
    with pytest.raises(FileExistsError):
        main(argv)


def test_cli_atomic_cleanup_on_injected_failure(tmp_path, monkeypatch):
    a_dir, b_dir = _build_two_categories(tmp_path)
    output = tmp_path / "out"
    argv = ["--result", f"A={a_dir}", "--result", f"B={b_dir}", "--output", str(output)]

    import normalize_tes.phi_contrast as phi_contrast

    def _boom(*args, **kwargs):
        raise RuntimeError("injected failure")

    monkeypatch.setattr(phi_contrast.json, "dump", _boom)
    with pytest.raises(RuntimeError, match="injected failure"):
        main(argv)

    assert not output.exists()
    # No leftover staging directories beside the intended output.
    leftovers = [p for p in output.parent.iterdir() if p.name.startswith(f".{output.name}.tmp.")]
    assert leftovers == []


def test_cli_warns_on_a_type_mismatch(tmp_path, capsys):
    a_dir = _write_result(tmp_path / "A", a_type="TE", z_score=1.0, null_z_scores=_rng_null(1, 10))
    b_dir = _write_result(tmp_path / "B", a_type="SNP", z_score=1.0, null_z_scores=_rng_null(2, 10))
    output = tmp_path / "out"
    argv = ["--result", f"A={a_dir}", "--result", f"B={b_dir}", "--output", str(output)]
    assert main(argv) == 0
    captured = capsys.readouterr()
    assert "a_type" in captured.out
    metadata = json.loads((output / "metadata.json").read_text())
    assert metadata["a_type_consistent"] is False


def test_cli_rejects_fewer_than_two_results(tmp_path):
    a_dir, _ = _build_two_categories(tmp_path)
    output = tmp_path / "out"
    argv = ["--result", f"A={a_dir}", "--output", str(output)]
    with pytest.raises(ValueError, match="at least two"):
        main(argv)


def test_cli_rejects_duplicate_labels(tmp_path):
    a_dir, b_dir = _build_two_categories(tmp_path)
    output = tmp_path / "out"
    argv = ["--result", f"A={a_dir}", "--result", f"A={b_dir}", "--output", str(output)]
    with pytest.raises(ValueError, match="duplicate"):
        main(argv)


def test_parse_args_defaults():
    args = parse_args(["--result", "A=x", "--result", "B=y", "--output", "z"])
    assert args.seed == 1002
    assert args.pairing_repeats == 100
