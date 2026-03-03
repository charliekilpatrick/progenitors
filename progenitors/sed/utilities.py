"""
SED utilities: dust maps, coordinate parsing, number formatting.
"""
from astropy.coordinates import SkyCoord
from astropy.table import Column
from astropy import units as u
import numpy as np
import dustmaps.sfd
import os


def is_number(val):
    """
    Check if a value can be converted to a float.

    Parameters
    ----------
    val : any
        Value to test.

    Returns
    -------
    bool
        True if float(val) succeeds, False otherwise.
    """
    try:
        float(val)
        return True
    except (ValueError, TypeError):
        return False


def import_dustmap():
    """
    Ensure SFD dust map data is available; fetch if missing.

    Checks for SFD_dust_4096_ngp.fits and SFD_dust_4096_sgp.fits in the
    dustmaps data directory and calls dustmaps.sfd.fetch() if needed.
    """
    for pole in ['ngp', 'sgp']:
        datadir = dustmaps.sfd.data_dir()
        file = os.path.join(datadir, 'sfd',
            'SFD_dust_4096_{}.fits'.format(pole))
        if not os.path.exists(file):
            dustmaps.sfd.fetch()
            break


def get_sfd():
    """
    Return an SFD dust map query object.

    Returns
    -------
    dustmaps.sfd.SFDQuery
        Query object for E(B-V) at (l, b) or coordinates.
    """
    return dustmaps.sfd.SFDQuery()


def round_to_n(f, n):
    """
    Round a float to n significant digits.

    Parameters
    ----------
    f : float
        Number to round.
    n : int
        Number of significant figures.

    Returns
    -------
    float
        Rounded value.
    """
    code = '%.{0}g'.format(int(n))
    val = float(code % f)
    return val


def parse_coord(ra, dec):
    """
    Parse RA/Dec into SkyCoord (degrees or sexagesimal).

    Parameters
    ----------
    ra : str, float, array-like, or Column
        Right ascension (e.g. '12:30:00' or 187.5).
    dec : str, float, array-like, or Column
        Declination (e.g. '+30:00:00' or 30.0).

    Returns
    -------
    SkyCoord or ndarray of SkyCoord
        Single coordinate or array of coordinates.
    """
    def _check(r, d):
        if ':' in str(r) and ':' in str(d):
            return SkyCoord(r, d, unit=(u.hour, u.deg))
        try:
            return SkyCoord(float(r), float(d), unit=(u.deg, u.deg))
        except (ValueError, TypeError):
            return None

    if (isinstance(ra, (list, np.ndarray, Column)) and
            isinstance(dec, (list, np.ndarray, Column))):
        return np.array([_check(r, d) for r, d in zip(ra, dec)])
    return _check(ra, dec)
