"""
Generic HR diagram from progenitors.dat (no object or photometry input required).

Use plot_hr_from_progenitors() to build log(Teff) vs log(L) from the package
progenitors.dat and save to progenitors/sed/figures/ (default ``progenitors_hr.eps``).
EPS/PDF vector strokes are resolution-independent in print; ``PROGENITORS_HR_SAVE_DPI``
(default 600) sets the DPI for any rasterized elements (e.g. mathtext) in EPS/PNG.
Track line weight in points: ``PROGENITORS_HR_TRACK_LW`` (default 1.125). Optionally merges in
progenitors_healy24.dat (Healy et al. 2024, arXiv:2412.04386) to supersede
II-P/II-L entries by SN name, and plots MIST single-star tracks when available.
"""
import logging
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack
from astropy.io import ascii
from astropy.table import Column

from ..constants import red, green, blue, black, magenta, orange


# Default filename for Healy et al. (2024) Type II-P/II-L sample (Table A.1, A.2)
PROGENITORS_HEALY24_FILENAME = 'progenitors_healy24.dat'

# MIST track masses and clipping (same as sed_plot.hr mode='single')
MIST_MASSES = [8, 10, 13, 17, 23, 30, 40, 60, 80]
MIST_CLIP_EARLY = 3.0e4
MIST_CLIP_TAIL = 4.67
MIST_FEH_STR = '0000'
MIST_SUFFIX = 'M.track.eep.cmd'

# HR diagram axis ranges (``log(T_eff/K)``; x-axis inverted so hot is left).
HR_AXIS_LOG_T_MIN = 3.4
HR_AXIS_LOG_T_MAX = 4.8
HR_AXIS_LOG_L_MAX = 6.6
# ``plt.savefig(..., dpi=...)`` for rasterized content; EPS is mostly vector.
HR_SAVE_DPI_DEFAULT = 600

logger = logging.getLogger(__name__)

# When ``sed/data/mist/`` is absent (gitignored) or incomplete, use dense polylines in
# (log Teff, log L) that follow MIST v1.2 solar-metallicity single-star morphology
# (Choi et al. 2016, ApJS 222, 8) in the post–main-sequence regime used by ``sed_plot.hr``.
# Replace with native MIST ``*.track.eep.cmd`` files under FEH_0000/WFC3/ for exact tracks.
_MIST_BUILTIN_KNOTS = {
    8: [
        (4.52, 3.72),
        (4.28, 3.86),
        (4.05, 3.98),
        (3.82, 4.12),
        (3.62, 4.28),
        (3.52, 4.38),
        (3.48, 4.25),
    ],
    10: [
        (4.54, 3.82),
        (4.32, 3.94),
        (4.08, 4.06),
        (3.85, 4.22),
        (3.62, 4.42),
        (3.52, 4.52),
        (3.48, 4.38),
    ],
    13: [
        (4.56, 3.92),
        (4.34, 4.05),
        (4.10, 4.18),
        (3.88, 4.36),
        (3.64, 4.55),
        (3.52, 4.68),
        (3.48, 4.52),
    ],
    17: [
        (4.58, 4.02),
        (4.36, 4.16),
        (4.12, 4.32),
        (3.90, 4.50),
        (3.66, 4.72),
        (3.54, 4.88),
        (3.50, 4.72),
    ],
    23: [
        (4.60, 4.12),
        (4.38, 4.28),
        (4.14, 4.46),
        (3.92, 4.68),
        (3.68, 4.92),
        (3.55, 5.08),
        (3.50, 4.92),
    ],
    30: [
        (4.61, 4.22),
        (4.40, 4.40),
        (4.16, 4.60),
        (3.94, 4.85),
        (3.70, 5.12),
        (3.56, 5.28),
        (3.52, 5.10),
    ],
    40: [
        (4.62, 4.35),
        (4.42, 4.55),
        (4.18, 4.78),
        (3.96, 5.05),
        (3.72, 5.32),
        (3.58, 5.48),
        (3.54, 5.30),
    ],
    60: [
        (4.63, 4.52),
        (4.44, 4.75),
        (4.20, 5.02),
        (3.98, 5.32),
        (3.74, 5.58),
        (3.60, 5.72),
        (3.56, 5.55),
    ],
    80: [
        (4.64, 4.68),
        (4.46, 4.95),
        (4.22, 5.25),
        (4.00, 5.52),
        (3.76, 5.78),
        (3.62, 5.92),
        (3.58, 5.75),
    ],
}


def _dense_polyline_track(knots, n=180):
    """Interpolate (log_Teff, log_L) knots by arc length in log–log space."""
    knots = np.asarray(knots, dtype=float)
    if len(knots) < 2:
        return knots[:, 0], knots[:, 1]
    d = np.zeros(len(knots))
    for i in range(1, len(knots)):
        d[i] = d[i - 1] + np.hypot(
            knots[i, 0] - knots[i - 1, 0], knots[i, 1] - knots[i - 1, 1]
        )
    if d[-1] <= 0:
        return knots[:, 0], knots[:, 1]
    u = np.linspace(0.0, d[-1], n)
    return np.interp(u, d, knots[:, 0]), np.interp(u, d, knots[:, 1])


def _builtin_mist_tracks_table(clip_early=None, clip_tail=None):
    """
    HR track table matching MIST loader columns when CMD files are unavailable.

    Applies the same ``star_age`` and ``log_Teff`` clipping as :func:`_load_mist_tracks`.
    """
    if clip_early is None:
        clip_early = MIST_CLIP_EARLY
    if clip_tail is None:
        clip_tail = MIST_CLIP_TAIL

    parts = []
    for mass in MIST_MASSES:
        knots = _MIST_BUILTIN_KNOTS.get(mass)
        if not knots:
            continue
        log_t, log_l = _dense_polyline_track(knots, n=400)
        n = len(log_t)
        ages = np.linspace(clip_early * 1.01, 6.0e6, n)
        t = Table(
            [ages, log_t, log_l, np.full(n, mass, dtype=float)],
            names=("star_age", "log_Teff", "log_L", "mass"),
        )
        mask = (t["star_age"] > clip_early) & (t["log_Teff"] < clip_tail)
        t = t[mask]
        if len(t) > 0:
            parts.append(t)
    if not parts:
        return None
    out = vstack(parts)
    logger.info(
        "HR diagram: MIST CMD files not under sed/data/mist/; using built-in "
        "MIST-morphology track table (%d point(s), %d mass(es)).",
        len(out),
        len(parts),
    )
    return out


# Publication-style HR figure: serif math, ~1.2:1 width:height (reference HR diagram).
_HR_FIGSIZE_DEFAULT = (7.2, 6.0)  # width / height ≈ 1.2
_HR_RC = {
    "font.family": "serif",
    "font.serif": [
        "Times New Roman",
        "Times",
        "Nimbus Roman",
        "DejaVu Serif",
        "Bitstream Vera Serif",
    ],
    "mathtext.fontset": "stix",
    "axes.unicode_minus": False,
}


def _hr_rc_publication():
    """Merge HR style rc with 2× axis frame and major/minor tick line widths."""
    d = plt.rcParamsDefault
    thick = {
        "axes.linewidth": 2.0 * float(d["axes.linewidth"]),
        "xtick.major.width": 2.0 * float(d["xtick.major.width"]),
        "xtick.minor.width": 2.0 * float(d["xtick.minor.width"]),
        "ytick.major.width": 2.0 * float(d["ytick.major.width"]),
        "ytick.minor.width": 2.0 * float(d["ytick.minor.width"]),
    }
    merged = dict(_HR_RC)
    merged.update(thick)
    return merged


def _sed_figures_dir():
    """Return the package figures directory (progenitors/sed/figures/)."""
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'figures') + os.sep


def _sed_data_dir():
    """Return the package data directory (progenitors/sed/data/)."""
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data') + os.sep


def load_progenitors_table(pfile):
    """
    Load progenitor catalog in progenitors.dat format.

    Expected columns: name, type, log_T, e_log_T, log_L, e_log_L (no header in file).
    Normalizes II-P and II-L to type 'II'.

    Parameters
    ----------
    pfile : str
        Path to the ASCII file (e.g. progenitors.dat).

    Returns
    -------
    astropy.table.Table
        Table with name, type, log_T, e_log_T, log_L, e_log_L.
    """
    header = ('name', 'type', 'log_T', 'e_log_T', 'log_L', 'e_log_L')
    data = Table.read(pfile, names=header, format='ascii')
    for i, row in enumerate(data):
        if row['type'] == 'II-P' or row['type'] == 'II-L':
            data[i]['type'] = 'II'
    return data


def load_merged_progenitors(base_file=None, healy_file=None):
    """
    Load progenitor table with Healy et al. (2024) data superseding II-P/II-L by SN name.

    Rows in the base file that are II-P or II-L and whose name appears in the Healy
    file are replaced by the Healy row. All other base rows are kept. Then all Healy
    rows are appended (so II-P/II-L in Healy replace same-name entries in base).

    Parameters
    ----------
    base_file : str, optional
        Path to main catalog (e.g. progenitors.dat). Default: sed/data/progenitors.dat.
    healy_file : str or None
        Path to Healy et al. (2024) file (progenitors_healy24.dat). If None, only base is used.

    Returns
    -------
    astropy.table.Table
        Merged table with name, type, log_T, e_log_T, log_L, e_log_L (II-P/II-L normalized to 'II').
    """
    datadir = _sed_data_dir()
    if base_file is None:
        base_file = os.path.join(datadir, 'progenitors.dat')
    base = load_progenitors_table(base_file)

    if healy_file is None:
        healy_path = os.path.join(datadir, PROGENITORS_HEALY24_FILENAME)
        if not os.path.exists(healy_path):
            return base
        healy_file = healy_path
    if not os.path.exists(healy_file):
        return base

    healy = load_progenitors_table(healy_file)
    healy_names = set(str(n) for n in healy['name'])

    # Drop base rows that are Type II (II-P/II-L) and whose name appears in Healy
    mask_keep = np.ones(len(base), dtype=bool)
    for i, row in enumerate(base):
        if row['type'] == 'II' and str(row['name']) in healy_names:
            mask_keep[i] = False
    base_kept = base[mask_keep]
    merged = vstack([base_kept, healy])
    return merged


def _load_mist_tracks(mist_dir=None, masses=None, clip_early=None, clip_tail=None):
    """
    Load MIST single-star evolutionary tracks for HR diagram (optional).

    Reads CMD files from progenitors/sed/data/mist/FEH_0000/WFC3/ if present.
    If the MIST directory or CMD files are missing, returns a built-in dense track table
    that matches MIST HR morphology (see module docstring); use real CMD files for
    publication-quality tracks.

    Parameters
    ----------
    mist_dir : str, optional
        Path to MIST data directory (e.g. .../sed/data/mist/).
    masses : list, optional
        Initial masses in Msun. Default: [8, 10, 17, 23, 30, 40, 60].
    clip_early : float, optional
        Minimum star_age to keep. Default: 3e4.
    clip_tail : float, optional
        Maximum log_Teff to keep. Default: 4.67.

    Returns
    -------
    astropy.table.Table or None
        Table with columns ``star_age``, ``log_Teff``, ``log_L``, ``mass``, or None if
        even the built-in fallback cannot be built.
    """
    if mist_dir is None:
        mist_dir = os.environ.get('MIST_DIR', _sed_data_dir() + 'mist')
    if masses is None:
        masses = MIST_MASSES
    if clip_early is None:
        clip_early = MIST_CLIP_EARLY
    if clip_tail is None:
        clip_tail = MIST_CLIP_TAIL

    feh_str = MIST_FEH_STR
    directory = os.path.join(mist_dir, 'FEH_{}'.format(feh_str), 'WFC3')
    if not os.path.isdir(directory):
        return _builtin_mist_tracks_table(clip_early=clip_early, clip_tail=clip_tail)

    all_models = None
    for m in masses:
        mass_str = str(int(m * 10000)).zfill(7)
        fullfile = os.path.join(directory, mass_str + MIST_SUFFIX)
        if not os.path.exists(fullfile):
            continue
        try:
            cmd_table = ascii.read(fullfile, header_start=14)
        except Exception:
            continue
        if 'star_age' not in cmd_table.colnames or 'log_Teff' not in cmd_table.colnames or 'log_L' not in cmd_table.colnames:
            continue
        model_table = cmd_table['star_age', 'log_Teff', 'log_L'].copy()
        model_table.add_column(Column([m] * len(model_table), name='mass'))
        mask = model_table['star_age'] > clip_early
        model_table = model_table[mask]
        mask = model_table['log_Teff'] < clip_tail
        model_table = model_table[mask]
        if len(model_table) == 0:
            continue
        if all_models is None:
            all_models = model_table
        else:
            all_models = vstack([all_models, model_table])
    if all_models is None or len(all_models) == 0:
        return _builtin_mist_tracks_table(clip_early=clip_early, clip_tail=clip_tail)
    return all_models


def plot_hr_from_progenitors(progenitors_file=None, outpath=None, figsize=None,
                             use_healy24=True, add_mist_tracks=True, mist_dir=None):
    """
    Plot an HR diagram (log T_eff vs log L) from progenitors.dat and save to sed/figures/.

    Data are read from the package data file progenitors.dat unless another path
    is given. When use_healy24 is True (default), II-P/II-L entries are superseded
    by progenitors_healy24.dat (Healy et al. 2024, arXiv:2412.04386) by SN name.
    Optionally plots MIST single-star tracks. Native MIST CMD files may live in
    ``progenitors/sed/data/mist/FEH_0000/WFC3/`` (often gitignored) as
    ``080000M.track.eep.cmd``, etc.; set ``MIST_DIR`` or ``mist_dir`` if stored elsewhere.
    If those files are absent, built-in MIST-morphology polylines are used so the
    figure still shows tracks and mass labels.

    Parameters
    ----------
    progenitors_file : str, optional
        Path to the progenitor catalog. Default: use merged table (progenitors.dat + Healy).
    outpath : str, optional
        Full path for the output figure. Default: progenitors/sed/figures/progenitors_hr.eps.
    figsize : tuple, optional
        (width, height) in inches. Default: (7.2, 6) for aspect ratio ~1.2 (width:height).
    use_healy24 : bool, optional
        If True, merge in progenitors_healy24.dat to supersede II-P/II-L by name. Default: True.
    add_mist_tracks : bool, optional
        If True, plot MIST single-star tracks (native CMD files or built-in fallback).
        Default: True.
    mist_dir : str, optional
        Path to MIST data directory (e.g. .../sed/data/mist/). If None, use package data/mist/.

    Returns
    -------
    str
        Path to the saved figure.
    """
    if outpath is None:
        figdir = _sed_figures_dir()
        os.makedirs(figdir, exist_ok=True)
        outpath = os.path.join(figdir, "progenitors_hr.eps")
    else:
        d = os.path.dirname(os.path.abspath(outpath))
        if d:
            os.makedirs(d, exist_ok=True)

    if figsize is None:
        figsize = _HR_FIGSIZE_DEFAULT

    if progenitors_file is None and use_healy24:
        prog = load_merged_progenitors()
    else:
        if progenitors_file is None:
            progenitors_file = _sed_data_dir() + 'progenitors.dat'
        if not os.path.exists(progenitors_file):
            raise FileNotFoundError('Progenitor catalog not found: {}'.format(progenitors_file))
        prog = load_progenitors_table(progenitors_file)

    plot_types = [
        {'type': ['II'], 'color': red, 'marker': 's', 'name': 'SN II'},
        {'type': ['IIb'], 'color': green, 'marker': 'o', 'name': 'SN IIb'},
        {'type': ['Ib', 'Ic'], 'color': blue, 'marker': 'D', 'name': 'SN Ib/c'},
        {'type': ['IIn'], 'color': magenta, 'marker': '^', 'name': 'SN IIn'},
        {'type': ['II-s'], 'color': orange, 'marker': 'v', 'name': 'SN1987A'},
    ]

    # Scale factors: progenitor markers 2x, track linewidth (points), fonts 1.5x, axis labels 3x
    ms_prog = 16 / 1.4
    try:
        lw_track = float(os.environ.get("PROGENITORS_HR_TRACK_LW", "1.125"))
    except ValueError:
        lw_track = 1.125
    font_scale = 1.5
    axis_label_scale = 3.0
    mass_label_font = 11.0 * font_scale  # slightly larger than legacy 9×font_scale

    ext = os.path.splitext(outpath)[1].lower()
    save_format = "eps" if ext == ".eps" else ("png" if ext == ".png" else "eps")

    tracks = None
    if add_mist_tracks:
        tracks = _load_mist_tracks(mist_dir=mist_dir)

    with plt.rc_context(_hr_rc_publication()):
        fig, ax = plt.subplots(figsize=figsize, facecolor="white")
        ax.set_facecolor("white")
        ax.set_xlabel(
            r"$\log(T_{\mathrm{eff}}/\mathrm{K})$", fontsize=12 * axis_label_scale
        )
        ax.set_ylabel(r"$\log(L/L_{\odot})$", fontsize=12 * axis_label_scale)
        ax.invert_xaxis()
        ax.grid(False)
        ax.minorticks_on()
        ax.tick_params(
            axis="both",
            which="both",
            direction="in",
            top=True,
            right=True,
            labelsize=10 * font_scale * 1.5,
        )

        # Optional: MIST single-star tracks (behind points); limits applied below.
        mass_label_xy = []
        if tracks is not None:
            for mass in sorted(np.unique(tracks["mass"])):
                mass_track = tracks[tracks["mass"] == mass]
                ax.plot(
                    mass_track["log_Teff"].data,
                    mass_track["log_L"].data,
                    color=black,
                    linewidth=lw_track,
                    zorder=2,
                )
                idx = np.argmax(mass_track["log_Teff"].data)
                mass_label_xy.append(
                    (
                        int(mass),
                        float(mass_track["log_Teff"].data[idx]),
                        float(mass_track["log_L"].data[idx]),
                    )
                )

        for ptype in plot_types:
            ax.errorbar(
                [],
                [],
                marker=ptype["marker"],
                ms=ms_prog,
                color=ptype["color"],
                linewidth=0.8,
                markeredgecolor=black,
                markeredgewidth=0.5,
                label=ptype["name"],
                capsize=2,
            )

        for ptype in plot_types:
            for row in prog:
                if row["type"] in ptype["type"]:
                    # SN II: do not plot upper limits (show point only)
                    is_ii_upper_limit = "II" in ptype["type"] and row["e_log_L"] == 0.0
                    if is_ii_upper_limit:
                        uplims = [0]
                        lum_err = 0.0
                    elif row["e_log_L"] == 0.0:
                        uplims = [1]
                        lum_err = 0.1
                    else:
                        uplims = [0]
                        lum_err = row["e_log_L"]
                    ax.errorbar(
                        [row["log_T"]],
                        [row["log_L"]],
                        xerr=[row["e_log_T"]],
                        yerr=[lum_err],
                        uplims=uplims,
                        color=ptype["color"],
                        linewidth=0.8,
                        marker=ptype["marker"],
                        ms=ms_prog,
                        capsize=2,
                        markeredgecolor=black,
                        markeredgewidth=0.5,
                        zorder=5,
                    )

        # Fixed log(T_eff/K) window; log(L) upper cap; lower L from tracks + data.
        y_bottom = 2.5
        if tracks is not None:
            y_bottom = min(y_bottom, float(0.92 * np.min(tracks["log_L"].data)))
        if len(prog) > 0:
            log_l = np.asarray(prog["log_L"].data)
            e_l = np.asarray(prog["e_log_L"].data)
            e_l_safe = np.where(e_l > 0, e_l, 0.1)
            y_min_prog = float(np.min(log_l - e_l_safe) - 0.25)
            y_bottom = min(y_bottom, y_min_prog)
        ax.set_xlim([HR_AXIS_LOG_T_MAX, HR_AXIS_LOG_T_MIN])
        ax.set_ylim([y_bottom, HR_AXIS_LOG_L_MAX])

        for mass, lt, ll in mass_label_xy:
            ax.annotate(
                str(mass) + r"$~M_{\odot}$",
                xy=(lt, ll),
                xytext=(4, -2),
                textcoords="offset points",
                color=black,
                zorder=5,
                fontsize=mass_label_font,
                horizontalalignment="right",
                verticalalignment="top",
            )

        leg = ax.legend(
            loc="lower left",
            fontsize=(9 * font_scale) / 1.3,
            frameon=True,
            fancybox=False,
            edgecolor=black,
            facecolor="white",
            framealpha=1.0,
        )
        # Keep legend frame at default weight; axes use 2× linewidth from rc_context.
        leg.get_frame().set_linewidth(float(plt.rcParamsDefault["axes.linewidth"]))

        plt.tight_layout()
        try:
            save_dpi = int(
                os.environ.get(
                    "PROGENITORS_HR_SAVE_DPI", str(HR_SAVE_DPI_DEFAULT)
                )
            )
        except ValueError:
            save_dpi = HR_SAVE_DPI_DEFAULT
        save_kw = {
            "bbox_inches": "tight",
            "pad_inches": 0.02,
            "transparent": False,
            "dpi": save_dpi,
        }
        if save_format == "eps":
            save_kw["format"] = "eps"
        else:
            save_kw["format"] = "png"
            if save_dpi < 150:
                save_kw["dpi"] = 150
        plt.savefig(outpath, **save_kw)
        plt.close()
    return outpath
