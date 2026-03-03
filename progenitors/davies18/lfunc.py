"""
Luminosity function pipeline: create_Lobs, create_modelgrid, compare_obs_mod.
Matches IDL Lfunc_obs_v2.pro logic for algorithm parity.
"""
import os
import numpy as np
from . import histutils
from . import powerlaw
from . import io_utils

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(_THIS_DIR, "data")
LEGACY_DIR = os.path.join(_THIS_DIR, "legacy")


def create_Lobs(
    obsfile="IIPprog_obsdata_2022.csv",
    NTRIALS=10000,
    LSTEP=0.01,
    LMIN=4.0,
    LMAX=5.8,
    seed=42,
    outpath=None,
):
    """
    Create observed P-distribution of L (log L) from CSV.
    Returns (LPOBS, NSN_OBS, obs_str), and optionally saves LPOBS to outpath.
    """
    rng = np.random.default_rng(seed)
    df = io_utils.read_obs_csv(obsfile)
    notes = df["notes"].fillna("").astype(str)
    L_us = df["L_us"].values
    dL_us = df["dL_us"].values
    # IDL: iupper = where(notes eq 'upper lim'); iuse = where(L_us le max(L_us[idetect]))
    iupper = np.where(notes.str.strip().str.lower() == "upper lim")[0]
    idetect = np.setdiff1d(np.arange(len(df)), iupper)
    if len(idetect) == 0:
        max_L = np.nanmax(L_us)
    else:
        max_L = np.nanmax(L_us[idetect])
    iuse = np.where(L_us <= max_L)[0]
    nsn = len(iuse)
    L_us = L_us[iuse]
    dL_us = dL_us[iuse]
    df_use = df.iloc[iuse].reset_index(drop=True)
    obs_str = io_utils.obs_df_to_struct(df_use)
    iupper_in_use = np.array([i for i in range(nsn) if iuse[i] in iupper], dtype=int)

    lhist = histutils.outxhist(LMIN, LMAX, LSTEP)
    nl = len(lhist)
    parr = np.zeros((nl, nsn))
    tarr = np.zeros((NTRIALS, nsn))

    for t in range(NTRIALS):
        Li = L_us + rng.normal(0, 1, nsn) * dL_us
        Li_nan = Li.copy()
        Li_nan[iupper_in_use] = np.nan
        order = np.argsort(np.where(np.isnan(Li_nan), np.inf, Li_nan))
        tarr[t, :] = Li_nan[order]

    for i in range(nsn):
        cur = tarr[:, i]
        cur = cur[np.isfinite(cur)]
        # Use LMAX (not LMAX-LSTEP) so bin count matches outxhist(LMIN,LMAX,LSTEP)
        nhist = histutils.opthist(cur, LMIN, LMAX, LSTEP)
        nhist = nhist.astype(float) / max(1.0, np.sum(nhist))
        parr[:, i] = nhist

    if outpath is None:
        outpath = os.path.join(LEGACY_DIR, "LPOBS.npz")
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    np.savez(outpath, LPOBS=parr, lhist=lhist, nsn=nsn)
    return parr, nsn, obs_str


def create_modelgrid(
    obs_str,
    LSTEP=0.01,
    LMIN=4.45,
    LMAX=5.28,
    NSN_OBS=24,
    NTRIALS=10000,
    NLLO=20,
    NLHI=19,
    NGAM=10,
    LLO_RANGE=(4.1, 4.5),
    LHI_RANGE=(5.1, 5.7),
    GAM_RANGE=(0.0, 2.0),
    TAPER_LLO=0.0,
    TAPER_LHI=0.0,
    seed=42,
    outname="LPMOD",
):
    """
    Build grid of model probability arrays for (L_lo, L_hi, gamma).
    Returns dict with Llo, Lhi, GAMMA_L, results (list of parr arrays).
    """
    rng = np.random.default_rng(seed)
    nsn_real = len(obs_str)
    nsn = NSN_OBS
    llo_min, llo_max = LLO_RANGE
    lhi_min, lhi_max = LHI_RANGE
    gam_min, gam_max = GAM_RANGE

    llo = np.linspace(llo_min, llo_max, NLLO)
    lhi = np.linspace(lhi_min, lhi_max, NLHI)
    pl = np.linspace(gam_min, gam_max, NGAM) if NGAM > 1 else np.array([-1.7])
    lhmin, lhmax, lhstep = LMIN, LMAX, LSTEP
    lxhist = histutils.outxhist(lhmin, lhmax, lhstep)
    nl = len(lxhist)
    dL_interp = np.interp(
        np.linspace(0, 1, nsn),
        np.linspace(0, 1, nsn_real),
        np.array([o["dL_us"] for o in obs_str]),
    )

    results = []
    for xx in range(len(llo)):
        for yy in range(len(lhi)):
            for pp in range(len(pl)):
                LRAN_MIN = 10 ** llo[xx]
                LRAN_MAX = 10 ** lhi[yy]
                PLCUR = pl[pp]
                if TAPER_LLO != 0:
                    LRAN_MIN = 10 ** (llo[xx] - TAPER_LLO)
                if TAPER_LHI != 0:
                    LRAN_MAX = 10 ** (lhi[yy] + TAPER_LHI)
                lums_mast = powerlaw.randomp(PLCUR, int(1e6), range_=(LRAN_MIN, LRAN_MAX), rng=rng)
                lums_mast = np.log10(lums_mast)
                larray = np.zeros((nsn, NTRIALS))
                for t in range(NTRIALS):
                    it = rng.integers(0, len(lums_mast), size=nsn)
                    lums_t = np.sort(lums_mast[it])
                    lums_t = lums_t + rng.normal(0, 1, nsn) * dL_interp
                    lums_t = np.sort(lums_t)
                    larray[:, t] = lums_t
                parr = np.zeros((nl, nsn))
                for i in range(nsn):
                    lyhist = histutils.opthist(larray[i, :], lhmin, lhmax, lhstep)
                    parr[:, i] = lyhist / float(NTRIALS)
                results.append(parr)
    nllo, nlhi, ngam = len(llo), len(lhi), len(pl)
    idx = 0
    results_grid = []
    for _ in range(nllo):
        row = []
        for _ in range(nlhi):
            col = []
            for _ in range(ngam):
                col.append(results[idx])
                idx += 1
            row.append(col)
        results_grid.append(row)
    return {
        "Llo": llo,
        "Lhi": lhi,
        "GAMMA_L": pl,
        "results": results_grid,
    }


def compare_obs_mod(
    LPOBS,
    LPMOD,
    LSTEP=0.01,
    LMIN=4.0,
    LMAX=5.8,
    NSN_OBS=None,
    silent=False,
    out_tex=None,
    out_contour=None,
):
    """
    Compare observed and model distributions; find best-fit (L_lo, L_hi, gamma).
    Returns BESTFITPARS dict and writes LFUNC_results.tex and LUM_CONTOUR.sav (or .npz).
    """
    Llo = LPMOD["Llo"]
    Lhi = LPMOD["Lhi"]
    GAML = LPMOD["GAMMA_L"]
    mod_probs = LPMOD["results"]
    lhist = histutils.outxhist(LMIN, LMAX, LSTEP)
    nl = len(lhist)
    nsn = NSN_OBS if NSN_OBS is not None else LPOBS.shape[1]
    nlo, nhi, ngam = len(Llo), len(Lhi), len(GAML)

    loglikarr = np.full((nlo, nhi, ngam), np.nan)
    probarr = np.full((nlo, nhi, ngam), np.nan)
    for xx in range(nlo):
        for yy in range(nhi):
            for gg in range(ngam):
                parr_mod = mod_probs[xx][yy][gg]
                with np.errstate(divide="ignore", invalid="ignore"):
                    loglik_i = np.nansum(np.log(parr_mod + 1e-300) + np.log(LPOBS + 1e-300))
                loglikarr[xx, yy, gg] = loglik_i
                prob_prof = np.nansum(LPOBS * parr_mod, axis=0)
                prob_prof = prob_prof[np.isfinite(prob_prof) & (prob_prof > 0)]
                prob_i = np.sum(np.log(prob_prof)) if len(prob_prof) else -np.inf
                probarr[xx, yy, gg] = prob_i

    chisqarr = -2.0 * probarr
    probarr_n = np.exp(probarr - np.nanmax(probarr))
    probarr_n = probarr_n / (np.nansum(probarr_n) * 1.2)
    imax = np.nanargmin(chisqarr)
    xyg = np.unravel_index(imax, chisqarr.shape)
    llo_best_i, lhi_best_i, gam_best_i = xyg
    minchisq = np.nanmin(chisqarr)
    lsig = minchisq + np.array([3.53, 7.81, 14.16])
    prob2d = np.nanmin(chisqarr, axis=2)
    chisq2d = np.nanmin(chisqarr, axis=2)
    i95 = np.where(chisqarr < lsig[0])
    xyg_95 = np.array(np.unravel_index(i95[0], chisqarr.shape))
    llo_lo = np.min(Llo[xyg_95[0]])
    llo_hi = np.max(Llo[xyg_95[0]])
    lhi_lo = np.min(Lhi[xyg_95[1]])
    lhi_hi = np.max(Lhi[xyg_95[1]])
    gam_lo = np.min(GAML[xyg_95[2]])
    gam_hi = np.max(GAML[xyg_95[2]])
    llo_m = np.outer(Llo, np.ones(nhi))
    lhi_m = np.outer(np.ones(nlo), Lhi)
    weights = np.exp(-0.5 * chisq2d)
    i1sig = np.where(chisq2d < lsig[0])
    llo_best_w = histutils.wmean(llo_m[i1sig], weights[i1sig])
    lhi_best_w = histutils.wmean(lhi_m[i1sig], weights[i1sig])
    llo_bf = np.array([llo_best_w, llo_hi - llo_best_w, llo_best_w - llo_lo])
    lhi_bf = np.array([lhi_best_w, lhi_hi - lhi_best_w, lhi_best_w - lhi_lo])
    gam_bf = np.array([GAML[gam_best_i], gam_hi - GAML[gam_best_i], GAML[gam_best_i] - gam_lo])
    BESTFITPARS = {"llo": llo_bf, "lhi": lhi_bf, "gam": gam_bf}

    if out_tex is None:
        out_tex = os.path.join(LEGACY_DIR, "LFUNC_results.tex")
    with open(out_tex, "w") as f:
        vars_ = [
            r"$\log(L_{\rm lo}/L_\odot)$",
            r"$\log(L_{\rm hi}/L_\odot)$",
            r"$\Gamma_{\rm L}$",
        ]
        for i, (var, v) in enumerate(zip(vars_, [llo_bf, lhi_bf, gam_bf])):
            f.write(f"{var} & = & $ {v[0]:.2f}^{{+{v[1]:.2f}}}_{{-{v[2]:.2f}}}$ \\\\\n")

    if out_contour is None:
        out_contour = os.path.join(LEGACY_DIR, "LUM_CONTOUR.npz")
    im_c = -chisq2d + np.nanmax(chisq2d)
    np.savez(out_contour, im_c=im_c, llo=Llo, lhi=Lhi, lev_c=lsig, chisq2d=chisq2d)
    return BESTFITPARS


def mainproc(
    obsfile="IIPprog_obsdata_2022.csv",
    NTRIALS_OBS=100000,
    NTRIALS_MOD=100000,
    LSTEP=0.01,
    LMIN=4.0,
    LMAX=5.8,
    NSN_OBS=25,
    NLLO=30,
    NLHI=29,
    NGAM=23,
    LLO_RANGE=(4.0, 4.6),
    LHI_RANGE=(5.0, 5.7),
    GAM_RANGE=(0.5, -2.2),
    seed=42,
):
    """Run full luminosity-function pipeline (create_Lobs -> create_modelgrid -> compare_obs_mod)."""
    if os.path.isfile(obsfile):
        obs_path = obsfile
    else:
        obs_path = io_utils.get_data_path(os.path.basename(obsfile))
    if not os.path.isfile(obs_path):
        obs_path = os.path.join(DATA_DIR, obsfile)
    LPOBS, nsn, obs_str = create_Lobs(
        obsfile=obs_path,
        NTRIALS=NTRIALS_OBS,
        LSTEP=LSTEP,
        LMIN=LMIN,
        LMAX=LMAX,
        seed=seed,
    )
    LPMOD = create_modelgrid(
        obs_str,
        LSTEP=LSTEP,
        LMIN=LMIN,
        LMAX=LMAX,
        NSN_OBS=nsn,
        NTRIALS=NTRIALS_MOD,
        NLLO=NLLO,
        NLHI=NLHI,
        NGAM=NGAM,
        LLO_RANGE=LLO_RANGE,
        LHI_RANGE=LHI_RANGE,
        GAM_RANGE=GAM_RANGE,
        seed=seed,
    )
    BESTFITPARS = compare_obs_mod(
        LPOBS,
        LPMOD,
        LSTEP=LSTEP,
        LMIN=LMIN,
        LMAX=LMAX,
        NSN_OBS=nsn,
    )
    return BESTFITPARS, LPOBS, LPMOD
