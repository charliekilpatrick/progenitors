"""
Power-law sampling matching IDL randomp.

IDL: randomp, lums, gamma, n, range=[Lmin, Lmax]
Samples L (linear scale) from pdf ∝ L^gamma between Lmin and Lmax.
Inverse CDF: for gamma != -1, CDF(L) = (L^(gamma+1) - Lmin^(gamma+1)) / (Lmax^(gamma+1) - Lmin^(gamma+1)).
So L = (u * (Lmax^(g+1) - Lmin^(g+1)) + Lmin^(g+1))^(1/(g+1)).
For gamma == -1: CDF = log(L/Lmin)/log(Lmax/Lmin), so L = Lmin * (Lmax/Lmin)^u.
"""
import numpy as np


def randomp(gamma, n, range_=(None, None), rng=None):
    """
    Sample n values from power-law pdf ∝ L^gamma on [Lmin, Lmax].

    Parameters
    ----------
    gamma : float
        Power-law index (pdf ∝ L^gamma).
    n : int
        Number of samples.
    range_ : tuple of (float, float)
        (Lmin, Lmax) in linear luminosity units.
    rng : np.random.Generator, optional
        Random generator for reproducibility.

    Returns
    -------
    ndarray
        Shape (n,) of luminosities in linear space.
    """
    if rng is None:
        rng = np.random.default_rng()
    Lmin, Lmax = range_
    if Lmin is None or Lmax is None:
        raise ValueError("range_ must be (Lmin, Lmax)")
    Lmin, Lmax = float(Lmin), float(Lmax)
    u = rng.uniform(0, 1, size=int(n))
    u = np.clip(u, 1e-300, 1 - 1e-300)
    g = float(gamma)
    if np.isclose(g, -1.0):
        L = Lmin * (Lmax / Lmin) ** u
    else:
        a = g + 1.0
        Lmin_a = Lmin ** a
        Lmax_a = Lmax ** a
        L = (u * (Lmax_a - Lmin_a) + Lmin_a) ** (1.0 / a)
    return L
