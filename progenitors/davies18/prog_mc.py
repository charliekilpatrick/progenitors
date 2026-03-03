"""
Progenitor mass MC pipeline: L -> mass via M-L relation, fit M_lo/M_hi per trial.
Uses data in davies18/data/. Matches IDL prog_MC.pro logic.
"""
import os
import numpy as np
from . import io_utils
from . import m_hi_lo_fit

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(_THIS_DIR, "data")
MABS_SUN = 4.74


def _lum_to_mass(lum, mass_mod, lfin_mod):
    log_mass = np.interp(lum, lfin_mod, np.log10(mass_mod))
    mass = 10.0 ** log_mass
    out = (lum > np.max(lfin_mod)) | (lum < np.min(lfin_mod))
    if np.any(out):
        mass[out] = 10.0 ** (-1.2417 + 0.47868 * lum[out])
    return mass


def run_prog_mc(
    obsfile="IIPprog_obsdata_2019.csv",
    MLfile="M-L_STARS.txt",
    ntrials=10000,
    include_uppers=True,
    seed=42,
):
    rng = np.random.default_rng(seed)
    df = io_utils.read_obs_csv(obsfile)
    mass_mod, lfin_mod = io_utils.read_ml_relation(MLfile)
    notes = df["notes"].fillna("").astype(str)
    iupper = np.where(notes.str.strip().str.lower() == "upper lim")[0]
    idetect = np.setdiff1d(np.arange(len(df)), iupper)
    nsn = len(df)
    dmod = df["dmod"].values
    edmod = df["edmod"].values
    mag = df["mag"].values
    emag = df["emag"].values
    alam = df["alam"].values
    ealam = df["ealam"].values
    BClam = df["BClam"].values
    eBClam = df["eBClam"].values
    L_us = df["L_us"].values
    dL_us = df["dL_us"].values
    L_s15 = df["L_s15"].values if "L_s15" in df.columns else L_us
    dL_s15 = df["dL_s15"].values if "dL_s15" in df.columns else dL_us
    snname = df["SN name"].values if "SN name" in df.columns else df["snname"].values

    dmod_t = dmod[:, None] + rng.normal(0, 1, (nsn, ntrials)) * edmod[:, None]
    BC_t = BClam[:, None] + (rng.uniform(0, 1, (nsn, ntrials)) - 0.5) * 2.0 * eBClam[:, None]
    A_t = np.maximum(alam[:, None] + rng.normal(0, 1, (nsn, ntrials)) * ealam[:, None], 0.0)
    mag_t = mag[:, None] + rng.normal(0, 1, (nsn, ntrials)) * emag[:, None]
    lums_t = (mag_t - dmod_t - A_t + BC_t - MABS_SUN) / (-2.5)
    lums_t = np.where(
        np.isfinite(L_us[:, None]),
        L_us[:, None] + rng.normal(0, 1, (nsn, ntrials)) * dL_us[:, None],
        lums_t,
    )
    masses_t = np.zeros((nsn, ntrials))
    for t in range(ntrials):
        masses_t[:, t] = _lum_to_mass(lums_t[:, t], mass_mod, lfin_mod)

    mlo_arr = np.zeros(ntrials)
    mhi_arr = np.zeros(ntrials)
    nn0 = np.arange(nsn) / max(nsn - 1.0, 1)
    for t in range(ntrials):
        masses = masses_t[:, t].copy()
        order = np.argsort(masses)
        masses = masses[order]
        if include_uppers:
            masses[iupper] = np.nan
        i_to_fit = np.where(np.isfinite(masses))[0]
        nfin = len(i_to_fit)
        if nfin < 2:
            mlo_arr[t], mhi_arr[t] = 8.0, 20.0
            continue
        masses_tofit = masses[i_to_fit]
        nn_tofit = nn0[i_to_fit] / np.nanmax(nn0[i_to_fit])
        (mlo, mhi), _, _ = m_hi_lo_fit.fit_m_hi_lo(nn_tofit, masses_tofit, p0=(8.0, 20.0))
        mlo_arr[t] = mlo
        mhi_arr[t] = mhi

    min1, max1 = 5.5, 10.5
    min2, max2 = 14.0, 40.0
    bin1, bin2 = 0.25, 0.75
    hist2d, xe, ye = np.histogram2d(
        mlo_arr, mhi_arr,
        bins=[np.arange(min1, max1 + bin1, bin1), np.arange(min2, max2 + bin2, bin2)],
    )
    hist2d = hist2d / np.sum(hist2d)
    minax = (xe[:-1] + xe[1:]) / 2
    maxax = (ye[:-1] + ye[1:]) / 2
    iorder = np.argsort(hist2d.ravel())[::-1]
    cumul = np.cumsum(hist2d.ravel()[iorder])
    pconf = np.interp([0.68, 0.95, 0.997], cumul, np.arange(len(cumul))).astype(int)
    lconf = hist2d.ravel()[iorder[pconf]]
    mmax_prof = np.nanmax(hist2d, axis=1)
    mmin_prof = np.nanmax(hist2d, axis=0)
    mmin_best = float(minax[np.argmax(mmax_prof)])
    mmax_best = float(maxax[np.argmax(mmin_prof)])
    masses_fit = np.column_stack([
        np.median(masses_t, axis=1),
        np.percentile(masses_t, 84, axis=1) - np.median(masses_t, axis=1),
        np.median(masses_t, axis=1) - np.percentile(masses_t, 16, axis=1),
    ])
    return {
        "mlo_arr": mlo_arr,
        "mhi_arr": mhi_arr,
        "hist2d": hist2d,
        "minax": minax,
        "maxax": maxax,
        "lconf": lconf,
        "mmin_best": mmin_best,
        "mmax_best": mmax_best,
        "masses_t": masses_t,
        "masses_fit": masses_fit,
        "snname": snname,
        "idetect": idetect,
        "iupper": iupper,
    }
