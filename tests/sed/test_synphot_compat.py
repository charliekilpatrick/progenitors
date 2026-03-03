"""Unit tests for progenitors.sed.synphot_compat (synphot backend)."""
import os
import numpy as np
import pytest

pytest.importorskip("synphot", reason="synphot required")
pytest.importorskip("astropy", reason="astropy required")


def test_setref():
    """setref is a no-op; should not raise."""
    from progenitors.sed import synphot_compat as S
    S.setref()
    S.setref(area=25.0)


def test_array_spectrum_basic():
    """ArraySpectrum from wave/flux arrays returns wrapper with .wave, .flux."""
    from progenitors.sed.synphot_compat import ArraySpectrum
    wave = np.linspace(3000, 8000, 100)
    flux = np.ones(100) * 1e-15
    sp = ArraySpectrum(wave, flux, waveunits='angstrom', fluxunits='flam', name='test')
    assert hasattr(sp, 'wave')
    assert hasattr(sp, 'flux')
    assert hasattr(sp, 'name')
    assert sp.wave is not None and sp.flux is not None
    assert len(sp.wave) == len(sp.flux)
    assert np.all(np.isfinite(sp.flux)) and np.all(sp.flux >= 0)
    assert sp.name == 'test'


def test_array_spectrum_fnu():
    """ArraySpectrum with fluxunits='fnu'."""
    from progenitors.sed.synphot_compat import ArraySpectrum
    wave = np.linspace(4000, 7000, 50)
    flux = np.ones(50) * 1e-20
    sp = ArraySpectrum(wave, flux, fluxunits='fnu')
    assert sp.wave is not None
    assert sp.flux is not None
    assert len(sp.flux) == len(sp.wave)


def test_blackbody():
    """BlackBody returns wrapper with .wave, .flux."""
    from progenitors.sed.synphot_compat import BlackBody
    bb = BlackBody(5000)
    assert hasattr(bb, 'wave')
    assert hasattr(bb, 'flux')
    assert bb.flux is not None
    assert bb.wave is not None
    assert np.all(bb.flux >= 0)


def test_spectrum_wrapper_mult_rmul():
    """Wrapper supports scalar multiply (left and right)."""
    from progenitors.sed.synphot_compat import ArraySpectrum
    wave = np.linspace(4000, 6000, 20)
    flux = np.ones(20)
    sp = ArraySpectrum(wave, flux, fluxunits='flam')
    sp2 = sp * 2.0
    sp3 = 3.0 * sp
    np.testing.assert_allclose(sp2.flux, 2.0 * sp.flux)
    np.testing.assert_allclose(sp3.flux, 3.0 * sp.flux)


def test_spectrum_wrapper_convert_flam():
    """convert('flam') keeps spectrum in flam units."""
    from progenitors.sed.synphot_compat import ArraySpectrum
    wave = np.linspace(4000, 6000, 30)
    flux = np.ones(30) * 1e-15
    sp = ArraySpectrum(wave, flux, fluxunits='flam')
    sp.convert('flam')
    assert sp.flux is not None
    assert len(sp.flux) == len(sp.wave)


def test_spectrum_wrapper_convert_fnu():
    """convert('fnu') changes to fnu."""
    from progenitors.sed.synphot_compat import ArraySpectrum
    wave = np.linspace(4000, 6000, 30)
    flux = np.ones(30) * 1e-15
    sp = ArraySpectrum(wave, flux, fluxunits='flam')
    sp.convert('fnu')
    assert sp.flux is not None
    assert len(sp.flux) == len(sp.wave)


def test_observation_abmag():
    """Observation with spectrum and bandpass gives effstim('abmag') float."""
    from progenitors.sed.synphot_compat import (
        ArraySpectrum, ArrayBandpass, Observation
    )
    # Use ArrayBandpass so we don't depend on johnson,V / stsynphot
    wave = np.linspace(3000, 10000, 200)
    trans = np.ones(200)
    trans[:20] = 0
    trans[180:] = 0
    bp = ArrayBandpass(wave, trans)
    flux = np.ones(200) * 1e-15
    sp = ArraySpectrum(wave, flux, fluxunits='flam')
    obs = Observation(sp, bp)
    mag = obs.effstim('abmag')
    assert isinstance(mag, (float, np.floating))
    assert np.isfinite(mag)


def test_array_bandpass():
    """ArrayBandpass from wave/trans arrays."""
    from progenitors.sed.synphot_compat import ArrayBandpass
    wave = np.linspace(3000, 8000, 100)
    trans = np.ones(100)
    trans[:20] = 0
    trans[80:] = 0
    bp = ArrayBandpass(wave, trans, waveunits='Angstrom')
    assert bp is not None
    # Can use in Observation
    from progenitors.sed.synphot_compat import ArraySpectrum, Observation
    sp = ArraySpectrum(wave, np.ones(100) * 1e-15, fluxunits='flam')
    obs = Observation(sp, bp)
    val = obs.effstim('abmag')
    assert np.isfinite(val)


def test_unwrap_spectrum():
    """_unwrap_spectrum returns internal synphot SourceSpectrum for wrapper."""
    from synphot import SourceSpectrum
    from progenitors.sed.synphot_compat import ArraySpectrum, _unwrap_spectrum
    sp = ArraySpectrum([4000, 5000], [1e-15, 1e-15], fluxunits='flam')
    inner = _unwrap_spectrum(sp)
    assert inner is not sp
    assert not isinstance(inner, type(sp))
    # Backend must be synphot SourceSpectrum (not pysynphot)
    assert isinstance(inner, SourceSpectrum)


def test_file_spectrum_requires_file():
    """FileSpectrum raises if file does not exist."""
    from progenitors.sed.synphot_compat import FileSpectrum
    # synphot may raise FileNotFoundError or (when closing) UnboundLocalError
    with pytest.raises((OSError, FileNotFoundError, UnboundLocalError)):
        FileSpectrum('/nonexistent/path/spectrum.fits')


def test_no_pysynphot_import():
    """synphot_compat must not import pysynphot (refactor uses synphot only)."""
    import sys
    mods_before = set(sys.modules.keys())
    from progenitors.sed import synphot_compat  # noqa: F401
    new_mods = set(sys.modules.keys()) - mods_before
    pysyn_mods = [m for m in new_mods if 'pysynphot' in m]
    assert not pysyn_mods, "synphot_compat should not load pysynphot; loaded: %s" % pysyn_mods


def test_get_spectral_data_root():
    """get_spectral_data_root prefers SYNPATH, falls back to PYSYN_CDBS."""
    from progenitors.sed.synphot_compat import get_spectral_data_root
    orig_syn = os.environ.pop('SYNPATH', None)
    orig_pysyn = os.environ.pop('PYSYN_CDBS', None)
    try:
        assert get_spectral_data_root() == ''
        os.environ['PYSYN_CDBS'] = '/legacy/cdbs'
        assert get_spectral_data_root() == '/legacy/cdbs'
        os.environ['SYNPATH'] = '/synphot/data'
        assert get_spectral_data_root() == '/synphot/data'
    finally:
        if orig_syn is not None:
            os.environ['SYNPATH'] = orig_syn
        elif 'SYNPATH' in os.environ:
            del os.environ['SYNPATH']
        if orig_pysyn is not None:
            os.environ['PYSYN_CDBS'] = orig_pysyn
        elif 'PYSYN_CDBS' in os.environ:
            del os.environ['PYSYN_CDBS']


def test_observation_uses_synphot():
    """Observation result is from synphot (effstim returns float magnitude)."""
    from synphot import Observation as SynphotObservation
    from progenitors.sed.synphot_compat import (
        ArraySpectrum, ArrayBandpass, Observation, _unwrap_spectrum,
    )
    wave = np.linspace(3500, 9000, 150)
    trans = np.ones(150)
    trans[:15] = 0
    trans[135:] = 0
    bp = ArrayBandpass(wave, trans)
    sp = ArraySpectrum(wave, np.ones(150) * 1e-15, fluxunits='flam')
    obs = Observation(sp, bp)
    mag = obs.effstim('abmag')
    assert np.isfinite(mag)
    # Internal observation is synphot Observation
    assert hasattr(obs, '_obs')
    assert isinstance(obs._obs, SynphotObservation)
