"""Unit tests for progenitors.davies18.prog_mc."""
import numpy as np
import pytest

from progenitors.davies18 import prog_mc


def test_run_prog_mc_returns_dict():
    out = prog_mc.run_prog_mc(
        obsfile="IIPprog_obsdata_2019.csv",
        MLfile="M-L_STARS.txt",
        ntrials=50,
        seed=42,
    )
    assert "mlo_arr" in out
    assert "mhi_arr" in out
    assert "hist2d" in out
    assert "minax" in out
    assert "maxax" in out
    assert "lconf" in out
    assert "mmin_best" in out
    assert "mmax_best" in out
    assert "snname" in out


def test_run_prog_mc_shapes():
    out = prog_mc.run_prog_mc(
        obsfile="IIPprog_obsdata_2019.csv",
        MLfile="M-L_STARS.txt",
        ntrials=100,
        seed=0,
    )
    assert len(out["mlo_arr"]) == 100
    assert len(out["mhi_arr"]) == 100
    assert out["hist2d"].ndim == 2
    assert len(out["minax"]) == out["hist2d"].shape[0]
    assert len(out["maxax"]) == out["hist2d"].shape[1]
    assert out["masses_fit"].shape[0] == len(out["snname"])


def test_run_prog_mc_reproducible():
    out1 = prog_mc.run_prog_mc(obsfile="IIPprog_obsdata_2019.csv", MLfile="M-L_STARS.txt", ntrials=100, seed=99)
    out2 = prog_mc.run_prog_mc(obsfile="IIPprog_obsdata_2019.csv", MLfile="M-L_STARS.txt", ntrials=100, seed=99)
    np.testing.assert_array_almost_equal(out1["mlo_arr"], out2["mlo_arr"])
    np.testing.assert_array_almost_equal(out1["mhi_arr"], out2["mhi_arr"])
