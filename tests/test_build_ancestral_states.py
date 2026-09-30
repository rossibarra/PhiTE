import json
from types import SimpleNamespace

import numpy as np
import pytest

import normalize_tes.build_ancestral_states as builder
from ancestral_table_helpers import stamp_ancestral_table


def _store(digest="store-digest", inputs=None):
    return SimpleNamespace(
        positions=np.arange(3, dtype=np.float64),
        metadata={
            "content_sha256": digest,
            "chromosomes": [],
            "sequence_length": 0,
            "inputs": inputs if inputs is not None else [],
        },
    )


def _part(path, draws):
    path.mkdir()
    np.save(path / "ancestral_counts.npy", np.zeros((3, 4), dtype=np.uint16))
    np.save(path / "present_draw_count.npy", np.zeros(3, dtype=np.uint16))
    (path / "metadata.json").write_text(json.dumps({
        "complete": True,
        "store_content_sha256": "store-digest",
        "bases": ["A", "C", "G", "T"],
        "draws": draws,
    }))
    stamp_ancestral_table(path)


def test_merge_requires_expected_draw_count(tmp_path, monkeypatch):
    monkeypatch.setattr(builder, "open_snp_age_store", lambda _: _store())
    part = tmp_path / "part"
    _part(part, [{"path": str(tmp_path / "a.tsz"), "sites": 4}])
    with pytest.raises(SystemExit, match="requires --expect-draws"):
        builder.main([
            "--store", str(tmp_path / "store"), "--output", str(tmp_path / "out"),
            "--merge", str(part),
        ])


def test_merge_rejects_duplicate_draw_inside_one_part(tmp_path, monkeypatch):
    inputs = [{"draw_id": 0, "path": str(tmp_path / "a.tsz")}]
    monkeypatch.setattr(builder, "open_snp_age_store",
                        lambda _: _store(inputs=inputs))
    draw = {"path": str(tmp_path / "a.tsz"), "sites": 4}
    part = tmp_path / "part"
    _part(part, [draw, draw])
    with pytest.raises(SystemExit, match="same draw more than once"):
        builder.main([
            "--store", str(tmp_path / "store"), "--output", str(tmp_path / "out"),
            "--merge", str(part), "--expect-draws", "1",
        ])


def test_merge_detects_same_path_even_if_site_metadata_differs(tmp_path, monkeypatch):
    inputs = [{"draw_id": 0, "path": str(tmp_path / "a.tsz")}]
    monkeypatch.setattr(builder, "open_snp_age_store",
                        lambda _: _store(inputs=inputs))
    draw_path = str(tmp_path / "a.tsz")
    first = tmp_path / "first"
    second = tmp_path / "second"
    _part(first, [{"path": draw_path, "sites": 4}])
    _part(second, [{"path": draw_path, "sites": 999}])
    with pytest.raises(SystemExit, match="already counted"):
        builder.main([
            "--store", str(tmp_path / "store"), "--output", str(tmp_path / "out"),
            "--merge", str(first), str(second), "--expect-draws", "2",
        ])


def test_merge_rejects_tree_arguments(tmp_path, monkeypatch):
    monkeypatch.setattr(builder, "open_snp_age_store", lambda _: _store())
    part = tmp_path / "part"
    _part(part, [{"path": str(tmp_path / "a.tsz"), "sites": 4}])
    with pytest.raises(SystemExit, match="cannot be combined"):
        builder.main([
            "--store", str(tmp_path / "store"), "--output", str(tmp_path / "out"),
            str(tmp_path / "b.tsz"), "--merge", str(part), "--expect-draws", "1",
        ])


def test_builder_rejects_store_without_content_identity(tmp_path, monkeypatch):
    monkeypatch.setattr(builder, "open_snp_age_store", lambda _: _store(None))
    with pytest.raises(SystemExit, match="no content_sha256"):
        builder.main([
            "--store", str(tmp_path / "store"), "--output", str(tmp_path / "out"),
            str(tmp_path / "draw.tsz"),
        ])


def _counted_part(path, draw_path, value):
    path.mkdir()
    counts = np.zeros((3, 4), dtype=np.uint16)
    counts[:, 0] = value
    np.save(path / "ancestral_counts.npy", counts)
    np.save(path / "present_draw_count.npy", np.full(3, value, dtype=np.uint16))
    (path / "metadata.json").write_text(json.dumps({
        "complete": True,
        "store_content_sha256": "store-digest",
        "bases": ["A", "C", "G", "T"],
        "draws": [{"path": draw_path, "sites": 3}],
    }))
    return stamp_ancestral_table(path)


def _two_draw_setup(tmp_path, monkeypatch):
    inputs = [{"draw_id": 0, "path": str(tmp_path / "a.tsz")},
              {"draw_id": 1, "path": str(tmp_path / "b.tsz")}]
    monkeypatch.setattr(builder, "open_snp_age_store",
                        lambda _: _store(inputs=inputs))
    first = _counted_part(tmp_path / "first", str(tmp_path / "a.tsz"), 1)
    second = _counted_part(tmp_path / "second", str(tmp_path / "b.tsz"), 1)
    return first, second


def _merge(tmp_path, *parts, expect=2):
    return builder.main([
        "--store", str(tmp_path / "store"), "--output", str(tmp_path / "out"),
        "--merge", *map(str, parts), "--expect-draws", str(expect),
    ])


def test_merge_publishes_v2_table_with_array_digests(tmp_path, monkeypatch):
    first, second = _two_draw_setup(tmp_path, monkeypatch)
    assert _merge(tmp_path, first, second) == 0
    out = tmp_path / "out"
    metadata = json.loads((out / "metadata.json").read_text())
    assert metadata["schema_version"] == builder.SCHEMA_VERSION == "ancestral-state-counts-v2"
    counts = np.load(out / "ancestral_counts.npy")
    present = np.load(out / "present_draw_count.npy")
    assert counts[:, 0].tolist() == [2, 2, 2]
    builder.verify_table_arrays(out, metadata, counts, present)


def test_merge_rejects_part_whose_array_no_longer_matches_its_digest(tmp_path, monkeypatch):
    first, second = _two_draw_setup(tmp_path, monkeypatch)
    counts = np.load(second / "ancestral_counts.npy")
    counts[0, 0] += 1
    np.save(second / "ancestral_counts.npy", counts)
    with pytest.raises(SystemExit, match="does not match its recorded digest"):
        _merge(tmp_path, first, second)


def test_merge_rejects_v1_part(tmp_path, monkeypatch):
    first, second = _two_draw_setup(tmp_path, monkeypatch)
    metadata = json.loads((second / "metadata.json").read_text())
    metadata["schema_version"] = "ancestral-state-counts-v1"
    del metadata["array_sha256"]
    (second / "metadata.json").write_text(json.dumps(metadata))
    with pytest.raises(SystemExit, match="is not an ancestral-state-counts-v2 table"):
        _merge(tmp_path, first, second)


def test_merge_rejects_part_from_another_store(tmp_path, monkeypatch):
    first, second = _two_draw_setup(tmp_path, monkeypatch)
    metadata = json.loads((second / "metadata.json").read_text())
    metadata["store_content_sha256"] = "other-store"
    (second / "metadata.json").write_text(json.dumps(metadata))
    with pytest.raises(SystemExit, match="different interval store"):
        _merge(tmp_path, first, second)


def test_merge_rejects_incomplete_part(tmp_path, monkeypatch):
    first, second = _two_draw_setup(tmp_path, monkeypatch)
    metadata = json.loads((second / "metadata.json").read_text())
    metadata["complete"] = False
    (second / "metadata.json").write_text(json.dumps(metadata))
    with pytest.raises(SystemExit, match="incomplete"):
        _merge(tmp_path, first, second)


def test_merge_rejects_a_missing_draw(tmp_path, monkeypatch):
    first, _ = _two_draw_setup(tmp_path, monkeypatch)
    with pytest.raises(SystemExit, match="do not cover exactly"):
        _merge(tmp_path, first, expect=1)


def test_merge_rejects_swapped_arrays(tmp_path, monkeypatch):
    first, second = _two_draw_setup(tmp_path, monkeypatch)
    counts = (second / "ancestral_counts.npy").read_bytes()
    present = (second / "present_draw_count.npy").read_bytes()
    (second / "ancestral_counts.npy").write_bytes(present)
    (second / "present_draw_count.npy").write_bytes(counts)
    with pytest.raises(SystemExit, match="has shape"):
        _merge(tmp_path, first, second)
