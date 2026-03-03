"""Unit tests for progenitors.davies18.m_hi_lo_fit."""
import numpy as np
import pytest

from progenitors.davies18 import m_hi_lo_fit


def test_m_hi_lo_fit_bounds():
    # At x=0 should be near m_lo, at x=1 near m_hi
    m_lo, m_hi = 8.0, 20.0
    m0 = m_hi_lo_fit.m_hi_lo_fit(0.0, m_lo, m_hi)
    m1 = m_hi_lo_fit.m_hi_lo_fit(1.0, m_lo, m_hi)
    assert m0 >= m_lo and m0 <= m_hi
    assert m1 >= m_lo and m1 <= m_hi
    assert m0 < m1


def test_m_hi_lo_fit_array():
    x = np.linspace(0, 1, 11)
    m = m_hi_lo_fit.m_hi_lo_fit(x, 8.0, 20.0)
    assert m.shape == x.shape
    assert np.all(np.diff(m) >= 0)  # monotonically increasing


def test_fit_m_hi_lo_roundtrip():
    # Data from exact power-law with M_lo=8, M_hi=25
    nn = np.linspace(0, 1, 30)
    masses = m_hi_lo_fit.m_hi_lo_fit(nn, 8.0, 25.0)
    (m_lo, m_hi), popt, _ = m_hi_lo_fit.fit_m_hi_lo(nn, masses, p0=(7.0, 26.0))
    assert 6 <= m_lo <= 11
    assert 22 <= m_hi <= 28


def test_fit_m_hi_lo_returns_tuple():
    nn = np.array([0.0, 0.5, 1.0])
    masses = np.array([8.5, 12.0, 19.0])
    (m_lo, m_hi), popt, pcov = m_hi_lo_fit.fit_m_hi_lo(nn, masses)
    assert isinstance(m_lo, (int, float))
    assert isinstance(m_hi, (int, float))
    assert popt.shape == (2,)
    assert pcov.shape == (2, 2)
