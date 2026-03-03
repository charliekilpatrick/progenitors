"""Unit tests for progenitors.davies18.histutils."""
import numpy as np
import pytest

from progenitors.davies18 import histutils


def test_opthist_basic():
    x = np.array([1.0, 1.2, 1.5, 2.0, 2.5])
    h = histutils.opthist(x, 1.0, 3.0, 0.5)
    assert len(h) == 4
    assert h.sum() == len(x)
    assert np.all(h >= 0)


def test_opthist_matches_outxhist_length():
    binmin, binmax, binsize = 4.0, 5.8, 0.01
    centers = histutils.outxhist(binmin, binmax, binsize)
    data = np.random.uniform(binmin, binmax - 1e-6, 1000)
    h = histutils.opthist(data, binmin, binmax, binsize)
    assert len(h) == len(centers)


def test_outxhist():
    x = histutils.outxhist(4.0, 5.8, 0.01)
    assert len(x) == 180
    assert 4.0 <= x[0] < 4.1
    assert 5.7 <= x[-1] <= 5.8


def test_outxhist_empty():
    x = histutils.outxhist(1.0, 1.0, 0.1)
    assert len(x) == 0


def test_wmean():
    vals = np.array([1.0, 2.0, 3.0])
    w = np.array([1.0, 1.0, 1.0])
    assert np.isclose(histutils.wmean(vals, w), 2.0)
    w2 = np.array([0.0, 1.0, 0.0])
    assert np.isclose(histutils.wmean(vals, w2), 2.0)


def test_wmean_zero_sum_weights():
    vals = np.array([1.0, 2.0])
    w = np.array([0.0, 0.0])
    out = histutils.wmean(vals, w)
    assert np.isnan(out)
