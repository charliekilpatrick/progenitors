"""
Generic HR diagram from progenitors.dat (no object or photometry input required).

Use plot_hr_from_progenitors() to build log(Teff) vs log(L) from the package
progenitors.dat and save to progenitors/sed/figures/. Optionally merges in
progenitors_healy24.dat (Healy et al. 2024, arXiv:2412.04386) to supersede
II-P/II-L entries by SN name, and plots MIST single-star tracks when available.
"""
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
    Returns None if the directory or files are missing (no dependency on full sed_fitter).

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
        Table with columns log_Teff, log_L, mass (and star_age), or None.
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
        return None

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
    return all_models


def plot_hr_from_progenitors(progenitors_file=None, outpath=None, figsize=None,
                             use_healy24=True, add_mist_tracks=True, mist_dir=None):
    """
    Plot an HR diagram (log T_eff vs log L) from progenitors.dat and save to sed/figures/.

    Data are read from the package data file progenitors.dat unless another path
    is given. When use_healy24 is True (default), II-P/II-L entries are superseded
    by progenitors_healy24.dat (Healy et al. 2024, arXiv:2412.04386) by SN name.
    Optionally plots MIST single-star tracks when data/mist/ is available. MIST data
    must be in ``progenitors/sed/data/mist/FEH_0000/WFC3/`` with files like
    ``080000M.track.eep.cmd`` (8 Msun), or set the ``MIST_DIR`` env var / pass
    ``mist_dir`` to point to your MIST data directory.

    Parameters
    ----------
    progenitors_file : str, optional
        Path to the progenitor catalog. Default: use merged table (progenitors.dat + Healy).
    outpath : str, optional
        Full path for the output figure. Default: progenitors/sed/figures/progenitors_hr.png.
    figsize : tuple, optional
        (width, height) in inches. Default: (8, 6).
    use_healy24 : bool, optional
        If True, merge in progenitors_healy24.dat to supersede II-P/II-L by name. Default: True.
    add_mist_tracks : bool, optional
        If True, plot MIST single-star tracks when data/mist/ exists. Default: True.
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
        outpath = os.path.join(figdir, 'progenitors_hr.png')
    else:
        d = os.path.dirname(os.path.abspath(outpath))
        if d:
            os.makedirs(d, exist_ok=True)

    if figsize is None:
        figsize = (8 * 1.3, 6)  # 30% wider than (8, 6)

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

    # Scale factors: progenitor markers 2x, MIST linewidth 0.5x, fonts 1.5x, axis labels 3x
    ms_prog = 16 / 1.4
    lw_mist = 0.2
    font_scale = 1.5
    axis_label_scale = 3.0

    fig, ax = plt.subplots(figsize=figsize)
    ax.set_xlabel(r'$\log(T_{\mathrm{eff}}/\mathrm{K})$', fontsize=12 * axis_label_scale)
    ax.set_ylabel(r'$\log(L/L_{\odot})$', fontsize=12 * axis_label_scale)
    ax.invert_xaxis()

    # Optional: MIST single-star tracks (behind points), styled as in original plotting.py hr()
    tracks = None
    if add_mist_tracks:
        tracks = _load_mist_tracks(mist_dir=mist_dir)
    figscale = figsize[0] if isinstance(figsize, (list, tuple)) else figsize
    if tracks is not None:
        ylim_t = [0.92 * np.min(tracks['log_L'].data), 1.02 * np.max(tracks['log_L'].data)]
        xlim_t = [1.045 * np.max(tracks['log_Teff'].data), 0.975 * np.min(tracks['log_Teff'].data)]
        ax.set_ylim(ylim_t)
        ax.set_xlim(xlim_t)
        for mass in sorted(np.unique(tracks['mass'])):
            mass_track = tracks[tracks['mass'] == mass]
            ax.plot(mass_track['log_Teff'].data, mass_track['log_L'].data,
                    color=black, linewidth=lw_mist * figscale, zorder=2)
            # Mass label at hottest point (as in original sed_plot.hr)
            idx = np.argmax(mass_track['log_Teff'].data)
            label_teff = mass_track['log_Teff'].data[idx]
            label_logL = mass_track['log_L'].data[idx]
            textpos = [label_teff + 0.0065 * figscale * (xlim_t[0] - xlim_t[1]),
                       label_logL - 0.004 * figscale * (ylim_t[1] - ylim_t[0])]
            ax.text(textpos[0], textpos[1], str(int(mass)) + r'$~M_{\odot}$',
                    color=black, zorder=5, fontsize=9 * font_scale,
                    horizontalalignment='center', verticalalignment='center',
                    bbox=dict(facecolor='white', edgecolor=black))

    for ptype in plot_types:
        ax.errorbar([], [], marker=ptype['marker'], ms=ms_prog, color=ptype['color'],
                    linewidth=0.8, markeredgecolor=black, markeredgewidth=0.5,
                    label=ptype['name'], capsize=2)

    for ptype in plot_types:
        for row in prog:
            if row['type'] in ptype['type']:
                # SN II: do not plot upper limits (show point only)
                is_ii_upper_limit = ('II' in ptype['type'] and row['e_log_L'] == 0.0)
                if is_ii_upper_limit:
                    uplims = [0]
                    lum_err = 0.0
                elif row['e_log_L'] == 0.0:
                    uplims = [1]
                    lum_err = 0.1
                else:
                    uplims = [0]
                    lum_err = row['e_log_L']
                ax.errorbar([row['log_T']], [row['log_L']],
                            xerr=[row['e_log_T']], yerr=[lum_err],
                            uplims=uplims, color=ptype['color'],
                            linewidth=0.8, marker=ptype['marker'], ms=ms_prog,
                            capsize=2, markeredgecolor=black, markeredgewidth=0.5, zorder=5)

    # Set axis limits so all data points (and labels) are within the frame
    if len(prog) > 0:
        log_t = np.asarray(prog['log_T'].data)
        log_l = np.asarray(prog['log_L'].data)
        e_t = np.asarray(prog['e_log_T'].data)
        e_l = np.asarray(prog['e_log_L'].data)
        e_l_safe = np.where(e_l > 0, e_l, 0.1)
        x_min_prog = np.min(log_t - e_t) - 0.08
        x_max_prog = np.max(log_t + e_t) + 0.08
        y_min_prog = np.min(log_l - e_l_safe) - 0.25
        y_max_prog = np.max(log_l + e_l_safe) + 0.3  # extra for SN1987A label
        if tracks is not None:
            # x is inverted: lims are [high, low]; expand to include prog
            ax.set_xlim([max(xlim_t[0], x_max_prog), min(xlim_t[1], x_min_prog)])
            ax.set_ylim([min(ylim_t[0], y_min_prog), max(ylim_t[1], y_max_prog)])
        else:
            ax.set_xlim([x_max_prog, x_min_prog])  # x inverted
            ax.set_ylim([y_min_prog, y_max_prog])

    ax.legend(loc='lower left', fontsize=(9 * font_scale) / 1.3)
    tick_label_scale = 1.5
    ax.tick_params(axis='both', labelsize=10 * font_scale * tick_label_scale)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(outpath, format='png', dpi=150)
    plt.close()
    return outpath
