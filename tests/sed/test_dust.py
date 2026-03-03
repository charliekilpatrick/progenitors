"""Unit tests for progenitors.sed.dust (rebin, get_dust, get_rsg, etc.)."""
import numpy as np
import pytest

# Dust module loads data at import; skip if data or synphot missing
pytest.importorskip("synphot", reason="synphot required for dust")
try:
    from progenitors.sed import dust
    _DUST_AVAILABLE = True
except Exception:
    _DUST_AVAILABLE = False


@pytest.mark.skipif(not _DUST_AVAILABLE, reason="dust data or deps not available")
class TestRebin:
    """Tests for rebin."""

    def test_rebin_shape(self):
        a = np.arange(100, dtype=float)
        out = dust.rebin(a, 10)
        assert out.shape == (10,)
        np.testing.assert_allclose(out[0], np.mean(a[:10]))
        np.testing.assert_allclose(out[-1], np.mean(a[90:]))

    def test_rebin_conserves_mean(self):
        a = np.ones(100)
        out = dust.rebin(a, 20)
        assert np.allclose(out, 1.0)


@pytest.mark.skipif(not _DUST_AVAILABLE, reason="dust data or deps not available")
class TestGetDust:
    """Tests for get_avg2, get_gvg2, get_avg10, get_avs2, get_avs10, get_dust."""

    def test_get_avg2_scalar_p(self):
        x = np.array([5000.0])  # wavelength in Angstrom
        p = 1.0  # single param for get_avg2
        a = dust.get_avg2(x, p)
        assert np.isscalar(a) or (isinstance(a, np.ndarray) and a.size >= 1)
        assert a >= 0 or np.all(np.isfinite(a))

    def test_get_gvg2(self):
        x = np.array([5000.0])
        p = np.array([1.0])  # get_gvg2 expects p[0]
        g = dust.get_gvg2(x, p)
        assert np.isscalar(g) or (isinstance(g, np.ndarray) and g.size >= 1)

    def test_get_dust_models(self):
        x = np.linspace(4000, 8000, 50)
        for model in ('g2', 's2', 'g10', 's10'):
            p = 0.5 if model == 'g2' else np.array([0.5])
            d = dust.get_dust(x, p, model=model)
            assert d.shape == x.shape
            assert np.all(d >= 0)

    def test_get_dust_negative_clipped(self):
        x = np.array([3000.0, 5000.0, 7000.0])
        d = dust.get_dust(x, 0.1, model='g2')
        assert np.all(d >= 0)


@pytest.mark.skipif(not _DUST_AVAILABLE, reason="dust data or deps not available")
class TestGetRsg:
    """Tests for get_rsg (uses module-level wavelength and rsg_10)."""

    def test_get_rsg_model_10(self):
        scale = 1.0  # L_sun
        temp = 3500.0
        flux = dust.get_rsg(scale, temp, model='10')
        assert flux.shape == dust.wavelength.shape
        assert np.all(flux >= 0)
        assert np.all(np.isfinite(flux))

    def test_get_rsg_temperature_in_range(self):
        # temp must be in dust.temp range for interpolation
        scale = 1.0
        flux = dust.get_rsg(scale, 3500.0, model='10')
        assert np.sum(flux) > 0


@pytest.mark.skipif(not _DUST_AVAILABLE, reason="dust data or deps not available")
class TestConstants:
    """Test module-level constants."""

    def test_flux_scale_positive(self):
        assert dust.FLUX_SCALE > 0

    def test_temp_array(self):
        assert hasattr(dust, 'temp')
        assert len(dust.temp) > 0
        assert np.all(dust.temp > 0)

    def test_wavelength_rebinned(self):
        assert dust.wavelength.shape == (7748,)
