"""
Mass-distribution and χ² utilities for progenitor analysis.

Uses the same IMF slope (γ = -1.35) as m_hi_lo_fit. Provides discrete
IMF grid sampling, weighted sampling, and χ² comparison for mass
distributions. For the full luminosity-function and progenitor mass MC
pipelines, use lfunc and prog_mc.
"""
import os
import numpy as np

from . import m_hi_lo_fit

GAMMA_IMF = m_hi_lo_fit.GAMMA  # -1.35


def generate_distribution(m_min, m_max, step=0.1):
    """
    Discrete IMF mass grid and normalized weights (pdf ∝ m^γ, γ = -1.35).

    Parameters
    ----------
    m_min, m_max : float
        Mass range (solar masses).
    step : float
        Grid step.

    Returns
    -------
    masses : ndarray
        1D mass grid.
    probs : ndarray
        Normalized probabilities (sum 1).
    """
    n = max(1, int(round((m_max - m_min) / step)))
    masses = np.linspace(m_min, m_max, n)
    probs = np.power(masses, GAMMA_IMF)
    probs = np.maximum(probs, 1e-300)
    probs = probs / np.sum(probs)
    return masses, probs


def generate_sample(masses, probs, size, rng=None):
    """
    Draw `size` masses from the discrete distribution (masses, probs).

    Parameters
    ----------
    masses, probs : array-like
        Same length; probs should sum to 1.
    size : int
        Number of samples.
    rng : np.random.Generator, optional
        Random generator.

    Returns
    -------
    ndarray
        Shape (size,) of sampled masses.
    """
    if rng is None:
        rng = np.random.default_rng()
    masses = np.asarray(masses)
    probs = np.asarray(probs, dtype=float)
    probs = probs / np.sum(probs)
    return rng.choice(masses, size=int(size), p=probs)


def calculate_chi2(input_masses, sim_masses, lims):
    """
    χ² between observed and simulated masses with per-point weights.

    chi2 = sum(lims * (input_masses - sim_masses)^2) / len(input_masses)^2

    Parameters
    ----------
    input_masses, sim_masses : array-like
        Same length.
    lims : array-like
        Weights (e.g. 0/1 for upper limits).

    Returns
    -------
    float
    """
    input_masses = np.asarray(input_masses, dtype=float)
    sim_masses = np.asarray(sim_masses, dtype=float)
    lims = np.asarray(lims, dtype=float)
    n = len(input_masses)
    return float(np.sum(lims * (input_masses - sim_masses) ** 2) / (n * n))


def run_mass_simulation(input_file="input_masses.txt", n_trials=100000, seed=None, verbose=True):
    """
    Run mass-distribution MC: load masses/lims from file, sample (m_min, m_max),
    generate IMF samples, compute χ². Requires davies18/data/input_masses.txt
    (two columns: mass, lim).

    Parameters
    ----------
    input_file : str
        Filename in davies18/data/ (or path).
    n_trials : int
        Number of MC trials.
    seed : int, optional
        Random seed.
    verbose : bool
        Print progress and sigma intervals.

    Returns
    -------
    dict
        vals: (n_trials, 2) (m_min, m_max); chi: (n_trials,) normalized χ².
    """
    from . import io_utils
    path = input_file if os.path.isfile(input_file) else io_utils.get_data_path(input_file)
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Need {path} to run mass simulation")
    rng = np.random.default_rng(seed)
    masses, lims = np.loadtxt(path, dtype=float, unpack=True)
    idx = np.argsort(masses)
    masses = masses[idx]
    lims = lims[idx]
    vals = []
    chi = []
    for i in range(n_trials):
        if verbose and i > 0 and i % 10000 == 0:
            print(i)
        m_min = float(rng.uniform(6.0, 10.0))
        m_max = float(rng.uniform(15.0, 35.0))
        mask = masses < m_max
        inp_m = masses[mask]
        inp_lims = lims[mask]
        vals.append((m_min, m_max))
        sim_m, prob = generate_distribution(m_min, m_max)
        samp = generate_sample(sim_m, prob, len(inp_m), rng=rng)
        samp = np.sort(samp)
        chi.append(calculate_chi2(inp_m, samp, inp_lims))
    chi = np.array(chi)
    vals = np.array(vals)
    chi = chi / np.min(chi)
    if verbose:
        for s in [1.0, 2.3, 3.5]:
            mask = chi < 1.0 + s
            m_min = vals[mask, 0]
            m_max = vals[mask, 1]
            print(s, np.min(m_max), np.max(m_max), np.median(m_max))
    return {"vals": vals, "chi": chi}
