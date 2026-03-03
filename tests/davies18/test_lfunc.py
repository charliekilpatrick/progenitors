"""Unit tests for progenitors.davies18.lfunc."""
import numpy as np
import pytest

from progenitors.davies18 import lfunc
from progenitors.davies18 import io_utils
from progenitors.davies18 import histutils


def test_create_Lobs_returns_shapes():
    parr, nsn, obs_str = lfunc.create_Lobs(
        obsfile="IIPprog_obsdata_2022.csv",
        NTRIALS=30,
        seed=42,
    )
    assert parr.ndim == 2
    assert parr.shape[1] == nsn
    assert len(obs_str) == nsn
    lhist = histutils.outxhist(4.0, 5.8, 0.01)
    assert parr.shape[0] == len(lhist)
    assert np.all(parr >= 0)
    # Detected SNe: normalized histograms sum to 1; upper-limit-only columns can sum to 0
    sums = parr.sum(axis=0)
    assert np.all(sums <= 1.0 + 1e-6)
    assert np.all(sums[sums > 0] >= 1.0 - 1e-6)


def test_create_modelgrid_small():
    df = io_utils.read_obs_csv("IIPprog_obsdata_2022.csv")
    obs_str = io_utils.obs_df_to_struct(df)
    LPMOD = lfunc.create_modelgrid(
        obs_str,
        NSN_OBS=min(10, len(obs_str)),
        NTRIALS=20,
        NLLO=2,
        NLHI=2,
        NGAM=2,
        seed=42,
    )
    assert "Llo" in LPMOD
    assert "Lhi" in LPMOD
    assert "GAMMA_L" in LPMOD
    assert "results" in LPMOD
    assert len(LPMOD["Llo"]) == 2
    assert len(LPMOD["Lhi"]) == 2
    assert len(LPMOD["GAMMA_L"]) == 2
    # results is 3D list [nllo][nlhi][ngam] of arrays (nl, nsn)
    r00 = LPMOD["results"][0][0][0]
    assert r00.ndim == 2
    assert np.all(r00 >= 0)


def test_compare_obs_mod():
    LMIN, LMAX = 4.0, 5.8
    parr, nsn, obs_str = lfunc.create_Lobs(
        obsfile="IIPprog_obsdata_2022.csv",
        NTRIALS=50,
        seed=42,
    )
    LPMOD = lfunc.create_modelgrid(
        obs_str,
        LMIN=LMIN,
        LMAX=LMAX,
        NSN_OBS=nsn,
        NTRIALS=100,
        NLLO=3,
        NLHI=3,
        NGAM=2,
        seed=42,
    )
    best = lfunc.compare_obs_mod(
        parr,
        LPMOD,
        LMIN=LMIN,
        LMAX=LMAX,
        NSN_OBS=nsn,
        out_tex=None,
        out_contour=None,
    )
    assert "llo" in best
    assert "lhi" in best
    assert "gam" in best
    assert len(best["llo"]) == 3
    assert len(best["lhi"]) == 3
    assert len(best["gam"]) == 3


def test_mainproc_small():
    best, LPOBS, LPMOD = lfunc.mainproc(
        obsfile="IIPprog_obsdata_2022.csv",
        NTRIALS_OBS=40,
        NTRIALS_MOD=50,
        NLLO=2,
        NLHI=2,
        NGAM=2,
        seed=42,
    )
    assert "llo" in best
    assert "lhi" in best
    assert LPOBS.ndim == 2
    assert "Llo" in LPMOD
