"""
Synphot compatibility layer for progenitors.sed.

Provides a uniform API on top of synphot for spectra, bandpasses, and
observations. Public API: setref, get_spectral_data_root, ArraySpectrum,
BlackBody, Observation, ObsBandpass, ArrayBandpass, FileSpectrum.
Uses synphot; HST/JWST filters require stsynphot.
"""
import os
import numpy as np
from astropy import units as u

from synphot import SourceSpectrum, SpectralElement, Observation as SynphotObservation
from synphot.models import Empirical1D, BlackBodyNorm1D
from synphot import units as synphot_units

# Optional: HST/JWST bandpasses
try:
    import stsynphot as stsyn
    _HAS_STSYN = True
except ImportError:
    stsyn = None
    _HAS_STSYN = False


def setref(area=None):
    """
    No-op for compatibility with legacy code.

    synphot passes area to countrate() when needed; not used for magnitude
    calculations.

    Parameters
    ----------
    area : float, optional
        Telescope area (ignored).
    """
    pass


def get_spectral_data_root():
    """
    Root directory for spectral data (e.g. Pickles, CDBS grid).

    Prefers SYNPATH (synphot convention); falls back to PYSYN_CDBS
    (stsynphot/legacy). Returns empty string if neither is set.

    Returns
    -------
    str
    """
    return os.environ.get("SYNPATH", os.environ.get("PYSYN_CDBS", ""))


class _SpectrumWrapper:
    """Wrap synphot SourceSpectrum with .wave, .flux, .name."""
    def __init__(self, sp, name=''):
        self._sp = sp
        self._name = name
        self._wave = None
        self._flux = None

    @property
    def waveset(self):
        return self._sp.waveset

    @property
    def wave(self):
        if self._wave is None:
            ws = self._sp.waveset
            self._wave = np.asarray(ws.to(u.AA).value if hasattr(ws, 'to') else ws)
        return self._wave

    @property
    def flux(self):
        if self._flux is None:
            f = self._sp(self._sp.waveset)
            self._flux = np.asarray(f.value if hasattr(f, 'value') else f)
        return self._flux

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, value):
        self._name = value

    def __mul__(self, other):
        if np.isscalar(other):
            return _SpectrumWrapper(self._sp * other, self._name)
        return NotImplemented

    def __rmul__(self, other):
        if np.isscalar(other):
            return _SpectrumWrapper(other * self._sp, self._name)
        return NotImplemented

    def convert(self, unit):
        """Convert spectrum to given flux unit (e.g. 'flam', 'fnu')."""
        if unit == 'flam':
            flux_unit = u.erg / u.s / u.cm**2 / u.AA
        elif unit == 'fnu':
            flux_unit = u.erg / u.s / u.cm**2 / u.Hz
        elif unit == 'count' or unit == 'counts':
            flux_unit = synphot_units.PHOTLAM
        else:
            flux_unit = u.Unit(unit) if isinstance(unit, str) else unit
        wave = self._sp.waveset
        f = self._sp(wave)
        flux_val = np.asarray(f.value if hasattr(f, 'value') else f)
        if unit == 'fnu' and hasattr(synphot_units, 'convert_flux'):
            flux_val = synphot_units.convert_flux(
                wave, flux_val, flux_unit
            ).value
        elif hasattr(f, 'to'):
            try:
                flux_val = f.to(flux_unit).value
            except Exception:
                pass
        self._sp = SourceSpectrum(
            Empirical1D, points=wave, lookup_table=flux_val * flux_unit
        )
        self._wave = None
        self._flux = None


def ArraySpectrum(wave, flux, waveunits='angstrom', fluxunits='flam', name=''):
    """
    Build a spectrum from wavelength and flux arrays (synphot backend).

    Returns a wrapper with .wave, .flux, .name.

    Parameters
    ----------
    wave : array-like
        Wavelength values.
    flux : array-like
        Flux values (same length as wave).
    waveunits : str, optional
        Wavelength unit (e.g. 'angstrom', 'nm'). Default is 'angstrom'.
    fluxunits : str, optional
        Flux unit: 'flam', 'fnu', 'count'/'counts', or astropy unit string.
        Default is 'flam'.
    name : str, optional
        Label for the spectrum.

    Returns
    -------
    _SpectrumWrapper
        Object with .wave, .flux, .name and .convert(unit) method.
    """
    wave = np.asarray(wave)
    flux = np.asarray(flux)
    if waveunits == 'angstrom' or waveunits == 'Angstrom':
        wave_u = wave * u.AA
    else:
        wave_u = wave * u.Unit(waveunits)
    if fluxunits == 'flam':
        flux_u = flux * (u.erg / u.s / u.cm**2 / u.AA)
    elif fluxunits == 'fnu':
        flux_u = flux * (u.erg / u.s / u.cm**2 / u.Hz)
    elif fluxunits == 'count' or fluxunits == 'counts':
        # Dimensionless (e.g. extinction curve) or PHOTLAM
        flux_u = flux * u.dimensionless_unscaled
    else:
        flux_u = flux * u.Unit(fluxunits)
    sp = SourceSpectrum(Empirical1D, points=wave_u, lookup_table=flux_u)
    return _SpectrumWrapper(sp, name=name)


def BlackBody(temp):
    """
    Blackbody spectrum at given temperature (synphot backend).

    Returns wrapper; use .convert('fnu') then scale.

    Parameters
    ----------
    temp : float
        Temperature in Kelvin.

    Returns
    -------
    _SpectrumWrapper
        Normalized blackbody with .wave, .flux.
    """
    sp = SourceSpectrum(BlackBodyNorm1D, temperature=temp)
    return _SpectrumWrapper(sp, name='')


def _unwrap_spectrum(sp):
    return sp._sp if isinstance(sp, _SpectrumWrapper) else sp


def _taper_bp(bp):
    if hasattr(bp, 'taper'):
        return bp.taper()
    return bp


class _ObservationWrapper:
    """Wrap synphot Observation to provide .effstim('abmag') returning float."""
    def __init__(self, obs):
        self._obs = obs

    def effstim(self, flux_unit):
        val = self._obs.effstim(flux_unit=flux_unit)
        return float(val.value) if hasattr(val, 'value') else float(val)


def Observation(spectrum, bandpass, binset=None, force='taper', **kwargs):
    """
    Observe a spectrum through a bandpass (synphot backend).

    Use .effstim('abmag') or .effstim('vegamag').

    Parameters
    ----------
    spectrum : _SpectrumWrapper or SourceSpectrum
        Source spectrum.
    bandpass : SpectralElement or array bandpass
        Filter bandpass.
    binset : array-like, optional
        Wavelength set for integration.
    force : str, optional
        If 'taper', taper bandpass edges. Default is 'taper'.

    Returns
    -------
    _ObservationWrapper
        Object with .effstim(flux_unit) returning a float (e.g. magnitude).
    """
    sp = _unwrap_spectrum(spectrum)
    bp = _taper_bp(bandpass) if force == 'taper' else bandpass
    obs = SynphotObservation(sp, bp, binset=binset)
    return _ObservationWrapper(obs)


def ObsBandpass(name):
    """
    Load a named filter bandpass.

    Uses synphot for standard names; stsynphot for HST/JWST (e.g. wfc3,uvis1,f438w).

    Parameters
    ----------
    name : str
        Filter name (e.g. 'johnson,V', 'wfc3,uvis1,f438w').

    Returns
    -------
    SpectralElement
        Bandpass transmission curve.

    Raises
    ------
    ImportError
        If HST/JWST filter requested and stsynphot is not installed.
    """
    name_lower = name.lower()
    if any(x in name_lower for x in ('wfc3', 'acs', 'wfpc2', 'jwst', 'nircam', 'miri')):
        if not _HAS_STSYN:
            raise ImportError(
                "HST/JWST filter '%s' requires stsynphot. "
                "Install with: pip install stsynphot" % name
            )
        return stsyn.band(name)
    try:
        return SpectralElement.from_filter(name)
    except Exception:
        if not _HAS_STSYN:
            raise ImportError(
                "Filter '%s' may require stsynphot. "
                "Install with: pip install stsynphot" % name
            )
        return stsyn.band(name)


def ArrayBandpass(wave, trans, name='', waveunits='Angstrom'):
    """
    Build a bandpass from wavelength and transmission arrays.

    Parameters
    ----------
    wave : array-like
        Wavelength values.
    trans : array-like
        Transmission (0-1) at each wavelength.
    name : str, optional
        Unused; for API compatibility.
    waveunits : str, optional
        Wavelength unit. Default is 'Angstrom'.

    Returns
    -------
    SpectralElement
        Bandpass for use with Observation.
    """
    wave = np.asarray(wave)
    trans = np.asarray(trans)
    if waveunits == 'Angstrom' or waveunits == 'angstrom':
        wave_u = wave * u.AA
    else:
        wave_u = wave * u.Unit(waveunits)
    return SpectralElement(Empirical1D, points=wave_u, lookup_table=trans)


def FileSpectrum(filename):
    """
    Load a spectrum from a file (synphot backend; e.g. FITS).

    Parameters
    ----------
    filename : str
        Path to spectrum file.

    Returns
    -------
    _SpectrumWrapper
        Object with .wave, .flux, .name.

    Raises
    ------
    OSError, FileNotFoundError
        If file cannot be opened.
    """
    sp = SourceSpectrum.from_file(filename)
    return _SpectrumWrapper(sp, name='')
