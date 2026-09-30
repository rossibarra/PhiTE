import numpy as np
import pytest

from tools.sim_polarity_arms import (
    bootstrap_matched_sets,
    replicate_bootstrap_interval,
)


def test_bootstrap_matched_sets_reuse_preserves_size_and_excludes_focal():
    focal = np.array([0, 1, 2, 3])
    pool = np.arange(4, 44)
    bins = np.arange(44) % 2
    sets = bootstrap_matched_sets(
        focal, pool, bins, 10, np.random.default_rng(2), 2, disjoint=False,
    )
    assert sets.shape == (10, 4)
    assert not np.intersect1d(focal, sets).size
    assert all(np.unique(row).size == row.size for row in sets)


def test_bootstrap_matched_sets_depletion_is_globally_disjoint():
    focal = np.array([0, 1, 2, 3])
    pool = np.arange(4, 44)
    bins = np.arange(44) % 2
    sets = bootstrap_matched_sets(
        focal, pool, bins, 4, np.random.default_rng(3), 2, disjoint=True,
    )
    assert sets.shape == (4, 4)
    assert np.unique(sets).size == sets.size
    assert not np.intersect1d(focal, sets).size


def test_bootstrap_matched_sets_returns_completed_prefix_on_depletion():
    focal = np.array([0, 1, 2, 3])
    pool = np.arange(4, 12)
    bins = np.zeros(12, dtype=int)
    sets = bootstrap_matched_sets(
        focal, pool, bins, 5, np.random.default_rng(4), 1, disjoint=True,
    )
    assert sets.shape == (2, 4)


def test_bootstrap_matched_sets_accepts_prespecified_bootstrap_targets():
    focal = np.array([0, 1, 2, 3])
    pool = np.arange(4, 44)
    bins = np.arange(44) % 2
    targets = np.array([[1, 3], [3, 1], [2, 2]])
    sets = bootstrap_matched_sets(
        focal, pool, bins, 3, np.random.default_rng(5), 2,
        disjoint=False, target_counts=targets,
    )
    observed = np.stack([np.bincount(bins[row], minlength=2) for row in sets])
    np.testing.assert_array_equal(observed, targets)
    with pytest.raises(ValueError, match="sum to focal size"):
        bootstrap_matched_sets(
            focal, pool, bins, 1, np.random.default_rng(5), 2,
            disjoint=False, target_counts=np.array([[1, 1]]),
        )


def test_replicate_bootstrap_interval_resamples_replicates():
    lo, hi = replicate_bootstrap_interval(
        np.array([0.0, 0.5, 1.0]), seed=8, draws=10_000,
    )
    assert 0.0 <= lo < 0.5 < hi <= 1.0
    with pytest.raises(ValueError, match="nonempty finite"):
        replicate_bootstrap_interval(np.array([]), seed=1)
