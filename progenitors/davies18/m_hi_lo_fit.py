"""
IMF-style mass fit (IDL m_hi_lo_fit.pro).

Cumulative mass at fraction x from power-law CDF between m_lo and m_hi
with slope GAMMA = -1.35. Used for progenitor mass distribution fitting.
"""
import numpy as np
from scipy.optimize import curve_fit

GAMMA = -1.35


def m_hi_lo_fit(x, m_lo, m_hi):
    """
    Mass at cumulative fraction x for power-law IMF between m_lo and m_hi.

    Parameters
    ----------
    x : array-like
        Cumulative fraction(s) in [0, 1].
    m_lo, m_hi : float
        Lower and upper mass bounds (solar masses).

    Returns
    -------
    ndarray or float
        Mass(es) at the given cumulative fraction(s).
    """
    x = np.asarray(x, dtype=float)
    a0, a1 = float(m_lo), float(m_hi)
    a0g = a0 ** GAMMA
    a1g = a1 ** GAMMA
    norm = 1.0 / (a1g - a0g)
    inner = np.clip(x / norm + a0g, 1e-300, np.inf)
    expo = (1.0 / GAMMA) * np.log10(inner)
    return 10.0 ** expo


def fit_m_hi_lo(nn_tofit, masses_tofit, p0=(8.0, 20.0)):
    """
    Fit m_lo and m_hi to mass vs cumulative fraction with scipy.optimize.curve_fit.

    Parameters
    ----------
    nn_tofit : array-like
        Cumulative fraction values (e.g. normalized rank).
    masses_tofit : array-like
        Sorted masses at those fractions.
    p0 : tuple of float, optional
        Initial guess (m_lo, m_hi). Default (8.0, 20.0).

    Returns
    -------
    tuple
        ((m_lo, m_hi), popt, pcov) from curve_fit.
    """
    nn_tofit = np.asarray(nn_tofit, dtype=float)
    masses_tofit = np.asarray(masses_tofit, dtype=float)
    popt, pcov = curve_fit(
        lambda x, m_lo, m_hi: m_hi_lo_fit(x, m_lo, m_hi),
        nn_tofit,
        masses_tofit,
        p0=p0,
        bounds=([0.1, 1.0], [200, 200]),
    )
    return (float(popt[0]), float(popt[1])), popt, pcov
