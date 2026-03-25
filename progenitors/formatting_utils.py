"""
Shared numeric and coordinate helpers used across the package.

Kept free of heavy optional dependencies (no astroquery, no pipeline I/O).
"""
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Column


def is_number(val):
    """
    Return True if *val* can be interpreted as a finite float.

    None and NaN are not considered numbers.
    """
    if val is None:
        return False
    try:
        x = float(val)
    except (ValueError, TypeError):
        return False
    if np.isnan(x):
        return False
    return True


def parse_coord(ra, dec):
    """
    Parse RA/Dec into ``SkyCoord`` (ICRS), degrees or sexagesimal strings.

    If both inputs are array-like (including ``astropy.table.Column``), returns
    a 1-D ``numpy.ndarray`` of ``SkyCoord``. On failure for a single pair,
    returns ``None``.
    """
    def _single(r, d):
        sr, sd = str(r), str(d)
        if ':' in sr and ':' in sd:
            try:
                return SkyCoord(r, d, unit=(u.hourangle, u.deg), frame='icrs')
            except (ValueError, TypeError):
                return None
        try:
            return SkyCoord(float(r), float(d), unit=(u.deg, u.deg), frame='icrs')
        except (ValueError, TypeError):
            return None

    if (isinstance(ra, (list, np.ndarray, Column)) and
            isinstance(dec, (list, np.ndarray, Column))):
        return np.array([_single(r, d) for r, d in zip(ra, dec)])
    return _single(ra, dec)


def round_to_decimal_places(value, n_places):
    """
    Format *value* as a string with *n_places* digits after the decimal point.

    Matches the legacy ``progenitors.util.round_to_n`` behavior (fixed decimals,
    not significant figures).
    """
    code = '%7.{0}f'.format(int(n_places))
    return (code % float(value)).strip()


def round_to_significant_figures(value, n_sigfig):
    """
    Round *value* to *n_sigfig* significant figures; return a float.

    Matches the legacy ``progenitors.sed.utilities.round_to_n`` behavior.
    """
    code = '%.{0}g'.format(int(n_sigfig))
    return float(code % float(value))
