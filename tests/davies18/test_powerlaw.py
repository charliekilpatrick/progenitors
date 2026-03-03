"""Unit tests for progenitors.davies18.powerlaw."""
import numpy as np
import pytest

from progenitors.davies18 import powerlaw


def test_randomp_shape():
    rng = np.random.default_rng(42)
    L = powerlaw.randomp(-1.5, 100, range_=(1e4, 1e5), rng=rng)
    assert L.shape == (100,)
    assert np.all(L >= 1e4)
    assert np.all(L <= 1e5)


def test_randomp_gamma_minus_one():
    rng = np.random.default_rng(0)
    L = powerlaw.randomp(-1.0, 500, range_=(10.0, 100.0), rng=rng)
    assert np.all(L >= 10.0)
    assert np.all(L <= 100.0)
    # Log-uniform in L
    logL = np.log10(L)
    assert np.min(logL) >= np.log10(10) and np.max(logL) <= np.log10(100)


def test_randomp_reproducible():
    rng = np.random.default_rng(123)
    L1 = powerlaw.randomp(-1.35, 50, range_=(1e4, 1e6), rng=rng)
    rng2 = np.random.default_rng(123)
    L2 = powerlaw.randomp(-1.35, 50, range_=(1e4, 1e6), rng=rng2)
    np.testing.assert_array_almost_equal(L1, L2)


def test_randomp_range_required():
    with pytest.raises(ValueError):
        powerlaw.randomp(-1.0, 10, range_=(None, 100.0))
    with pytest.raises(ValueError):
        powerlaw.randomp(-1.0, 10, range_=(10.0, None))
