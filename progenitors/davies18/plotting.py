"""
Matplotlib plotting routines that produce EPS figures similar to IDL output.
Figures: luminosity contour (L_lo vs L_hi), Mlo-Mhi contour, lum comparison, mass spectrum.
"""
import os
import numpy as np


def _plt():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
FIGS_DIR = os.path.join(_THIS_DIR, "figs")
LEGACY_DIR = os.path.join(_THIS_DIR, "legacy")


def _ensure_figs_dir():
    os.makedirs(FIGS_DIR, exist_ok=True)


def plot_lum_contour(im_c, llo, lhi, lev_c=None, outpath=None):
    """
    Plot luminosity-function contour (L_lo vs L_hi) with chi^2 levels.
    im_c: 2D array (nllo, nlhi) to display; lev_c: contour levels (e.g. 1, 2, 3 sigma).
    """
    plt = _plt()
    _ensure_figs_dir()
    if outpath is None:
        outpath = os.path.join(FIGS_DIR, "LUM_CONTOUR.eps")
    fig, ax = plt.subplots(1, 1, figsize=(6, 5))
    Llo, Lhi = np.meshgrid(llo, lhi, indexing="ij")
    if lev_c is not None:
        ax.contour(Lhi, Llo, im_c, levels=np.sort(lev_c)[::-1], colors=["yellow", "orange", "red"], linewidths=1.5)
    ax.contourf(Lhi, Llo, im_c, levels=20, cmap="viridis_r", alpha=0.7)
    ax.set_xlabel(r"$\log(L_{\rm hi}/L_\odot)$")
    ax.set_ylabel(r"$\log(L_{\rm lo}/L_\odot)$")
    ax.set_title("Luminosity function fit")
    fig.tight_layout()
    fig.savefig(outpath, format="eps", bbox_inches="tight")
    plt.close(fig)


def plot_mlo_mhi_contour(
    hist2d,
    minax,
    maxax,
    lconf,
    mmin_best,
    mmax_best,
    title="with upper limits",
    outpath=None,
):
    """M_min vs M_max 2D histogram with confidence contours and best-fit (IDL Mlo-Mhi style)."""
    plt = _plt()
    _ensure_figs_dir()
    if outpath is None:
        outpath = os.path.join(FIGS_DIR, "Mlo-Mhi.eps")
    fig, ax = plt.subplots(1, 1, figsize=(6, 5))
    Mlo, Mhi = np.meshgrid(minax, maxax, indexing="ij")
    ax.contourf(Mhi, Mlo, hist2d, levels=20, cmap="viridis_r", alpha=0.8)
    ax.contour(Mhi, Mlo, hist2d, levels=np.sort(lconf)[::-1], colors=["yellow", "orange", "red"], linewidths=1)
    ax.plot(mmax_best, mmin_best, "o", color="white", markersize=10, markeredgecolor="black")
    ax.set_xlabel(r"$M_{\rm max}/M_\odot$")
    ax.set_ylabel(r"$M_{\rm min}/M_\odot$")
    ax.set_title(title)
    ax.text(0.95, 0.95, r"$M_{\rm min}=%.1f$, $M_{\rm max}=%.1f$" % (mmin_best, mmax_best),
            transform=ax.transAxes, ha="right", va="top", fontsize=10)
    fig.tight_layout()
    fig.savefig(outpath, format="eps", bbox_inches="tight")
    plt.close(fig)


def plot_lum_comp(L_us, L_s15, dL_us, dL_s15, outpath=None):
    """Luminosity comparison: L_us vs L_s15 with error bars (1:1 line)."""
    plt = _plt()
    _ensure_figs_dir()
    if outpath is None:
        outpath = os.path.join(FIGS_DIR, "lum_comp.eps")
    fig, ax = plt.subplots(1, 1, figsize=(6, 5))
    ax.errorbar(L_us, L_s15, xerr=dL_us, yerr=dL_s15, fmt="o", capsize=2)
    lims = [min(np.nanmin(L_us), np.nanmin(L_s15)), max(np.nanmax(L_us), np.nanmax(L_s15))]
    ax.plot(lims, lims, "k--", alpha=0.7)
    ax.set_xlabel(r"$L/L_\odot$ (this work)")
    ax.set_ylabel(r"$L_{\rm S15}/L_\odot$")
    ax.set_title("Luminosity comparison")
    ax.set_aspect("equal")
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    fig.tight_layout()
    fig.savefig(outpath, format="eps", bbox_inches="tight")
    plt.close(fig)


def plot_mass_spec(masses_fit, snname, mso, mmin_best, mmax_best, outpath=None):
    """Mass spectrum: progenitor mass vs SN name with best-fit IMF curve."""
    plt = _plt()
    from . import m_hi_lo_fit
    _ensure_figs_dir()
    if outpath is None:
        outpath = os.path.join(FIGS_DIR, "mass_spec.eps")
    nsn = len(snname)
    yax = np.arange(nsn)
    yax0 = np.arange(nsn) / max(nsn - 1.0, 1)
    xax = m_hi_lo_fit.m_hi_lo_fit(yax0, mmin_best, mmax_best)
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    ax.errorbar(
        masses_fit[mso, 0],
        yax,
        xerr=[masses_fit[mso, 2], masses_fit[mso, 1]],
        fmt="o",
        capsize=2,
    )
    ax.plot(xax, yax, "b-", lw=2, label=r"$M_{\rm min}=%.1f$, $M_{\rm max}=%.1f$" % (mmin_best, mmax_best))
    ax.set_yticks(yax)
    ax.set_yticklabels([str(snname[i]) for i in mso], fontsize=8)
    ax.set_xlabel(r"$M_{\rm init}/M_\odot$")
    ax.set_ylabel("SN progenitor")
    ax.set_ylim(-0.5, nsn - 0.5)
    ax.legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(outpath, format="eps", bbox_inches="tight")
    plt.close(fig)


def plot_all_from_compare(contour_npz_path=None):
    """After compare_obs_mod: plot luminosity contour from saved LUM_CONTOUR.npz."""
    path = contour_npz_path or os.path.join(LEGACY_DIR, "LUM_CONTOUR.npz")
    if os.path.isfile(path):
        d = np.load(path)
        plot_lum_contour(
            d["im_c"], d["llo"], d["lhi"], lev_c=d["lev_c"],
            outpath=os.path.join(FIGS_DIR, "LUM_CONTOUR.eps"),
        )
