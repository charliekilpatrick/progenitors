"""Unit tests for progenitors.davies18.distribution (mass-distribution modeling)."""
import numpy as np
import pytest

from progenitors.davies18 import distribution


def test_generate_distribution():
    masses, probs = distribution.generate_distribution(10.0, 15.0)
    assert len(masses) == int((15 - 10) / 0.1)
    assert len(probs) == len(masses)
    assert np.all(probs > 0)
    assert np.isclose(probs.sum(), 1.0)
    assert masses[0] == 10.0
    assert masses[-1] == 15.0


def test_generate_sample():
    masses = np.array([1.0, 2.0, 3.0])
    probs = np.array([0.2, 0.5, 0.3])
    rng = np.random.default_rng(42)
    out = distribution.generate_sample(masses, probs, 100, rng=rng)
    assert out.shape == (100,)
    assert np.all(np.isin(out, masses))


def test_calculate_chi2():
    inp = np.array([10.0, 12.0, 14.0])
    sim = np.array([10.1, 11.9, 14.2])
    lims = np.array([1.0, 1.0, 1.0])
    chi2 = distribution.calculate_chi2(inp, sim, lims)
    assert chi2 >= 0 and np.isfinite(chi2)
    assert distribution.calculate_chi2(inp, inp, lims) == 0.0


def test_run_mass_simulation_smoke():
    result = distribution.run_mass_simulation(
        input_file="input_masses.txt",
        n_trials=50,
        seed=0,
        verbose=False,
    )
    assert "vals" in result and "chi" in result
    assert result["vals"].shape[0] == 50
    assert result["chi"].shape[0] == 50
