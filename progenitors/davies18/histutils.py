"""
Histogram utilities matching IDL histutils.pro / Lfunc_obs_v2.pro.

- opthist: histogram with fixed binmin, binmax, binsize
- outxhist: bin centers for that grid
- wmean: weighted mean
"""
import numpy as np


def opthist(var, binmin, binmax, binsize):
    """
    Histogram with fixed bin edges. Number of bins = (binmax - binmin) / binsize
    so that it can match outxhist(binmin, binmax, binsize) length when binmax is
    LMAX (then n = (LMAX-LMIN)/LSTEP). Use binmax = LMAX - LSTEP for IDL-compatible
    179 bins with LMAX=5.8, LSTEP=0.01.

    Parameters
    ----------
    var : array-like
        Data to histogram.
    binmin, binmax : float
        Range (max is exclusive for last bin).
    binsize : float
        Bin width.

    Returns
    -------
    ndarray
        Counts in each bin.
    """
    var = np.asarray(var, dtype=float)
    n = (binmax - binmin) / binsize
    nbins = int(round(n))
    if nbins <= 0:
        return np.array([], dtype=float)
    edges = binmin + np.arange(nbins + 1, dtype=float) * binsize
    outhist, _ = np.histogram(var, bins=edges)
    return outhist.astype(float)


def outxhist(binmin, binmax, binsize):
    """
    Bin center positions for the grid used by opthist.
    Matches IDL: findgen((binmax-binmin)/binsize)/((binmax-binmin)/binsize)*(binmax-binmin)+binmin+binsize/2

    Parameters
    ----------
    binmin, binmax : float
        Range.
    binsize : float
        Bin width.

    Returns
    -------
    ndarray
        Bin centers.
    """
    n = (binmax - binmin) / binsize
    n = int(round(n))
    if n <= 0:
        return np.array([])
    xarray = np.arange(n, dtype=float) / n * (binmax - binmin) + binmin + binsize / 2.0
    return xarray


def wmean(vals, weights):
    """
    Weighted mean: sum(vals * weights) / sum(weights).

    Parameters
    ----------
    vals : array-like
        Values.
    weights : array-like
        Weights (same length as vals). NaNs are ignored in the sum.

    Returns
    -------
    float
        Weighted mean, or np.nan if sum(weights) is zero.
    """
    vals = np.asarray(vals, dtype=float)
    weights = np.asarray(weights, dtype=float)
    dw = np.nansum(weights)
    if dw == 0:
        return np.nan
    return np.nansum(vals * weights) / dw
