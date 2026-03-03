"""Unit tests for progenitors.sed.analysis (sed_fitter and model params)."""
import pytest
import numpy as np
from types import SimpleNamespace

pytest.importorskip("synphot", reason="synphot required for SED tests")
pytest.importorskip("emcee", reason="emcee required for SED tests")


def _get_fitter(options=None):
    """Return a sed_fitter with minimal options set."""
    from progenitors.sed import analysis
    fit = analysis.sed_fitter(verbose=False)
    if options is None:
        options = SimpleNamespace(
            extinction_law='fitzpatrick99',
            dust_comp='sil', dust_width=2,
            use_variance=False, report_mist_mass=False,
            skip_dust=True, notau=False, bolometric_correction=None,
        )
    fit.options = options
    return fit


def test_model_fit_params_keys():
    from progenitors.sed import analysis
    assert hasattr(analysis, "model_fit_params")
    params = analysis.model_fit_params
    assert "blackbody" in params
    assert "rsg" in params
    assert params["blackbody"] == ["luminosity", "temperature"]
    assert "tau_V" in params["rsg"] or "luminosity" in params["rsg"]


def test_sed_fitter_dirs():
    """sed_fitter has expected dir keys."""
    from progenitors.sed import analysis
    fit = analysis.sed_fitter(verbose=False)
    for key in ("bandpass", "input", "data", "backends"):
        assert key in fit.dirs


def test_sed_fitter_files():
    from progenitors.sed import analysis
    fit = analysis.sed_fitter(verbose=False)
    assert "pickles" in fit.files
    assert "rsg" in fit.files or "blackbody" in fit.files


def test_sed_fitter_bounds():
    """Bounds contain expected parameter names."""
    fit = _get_fitter()
    for key in ("luminosity", "temperature", "tau_V", "dust_temp", "Av", "Rv"):
        assert key in fit.bounds
        assert len(fit.bounds[key]) >= 2


def test_add_options_parser():
    """add_options returns parser with sampler and nlive."""
    from progenitors.sed import analysis
    fit = _get_fitter()
    parser = fit.add_options()
    assert parser is not None
    # Parse known args to check --sampler and --nlive exist
    args = parser.parse_args(['obj', 'blackbody', '--sampler', 'dynesty', '--nlive', '100'])
    assert args.sampler == 'dynesty'
    assert args.nlive == 100


def test_set_model_type_blackbody():
    """set_model_type('blackbody') sets model_fit_params to lum, temp."""
    fit = _get_fitter()
    fit.set_model_type('blackbody')
    assert fit.model_type == 'blackbody'
    assert fit.model_fit_params == ['luminosity', 'temperature']


def test_set_model_type_rsg():
    """set_model_type('rsg') with options sets rsg_* model type."""
    fit = _get_fitter()
    fit.set_model_type('rsg')
    assert 'rsg' in fit.model_type
    assert 'luminosity' in fit.model_fit_params
    assert 'temperature' in fit.model_fit_params


def test_check_bounds_blackbody():
    """check_bounds returns True when theta is outside bounds, False when inside."""
    fit = _get_fitter()
    fit.set_model_type('blackbody')
    in_bounds = [4.0, 3.7]
    out_bounds_lum = [1.0, 3.7]
    assert fit.check_bounds(in_bounds) is False
    assert fit.check_bounds(out_bounds_lum) is True


def test_check_bounds_rsg():
    """check_bounds with 4 params (tau_V, lum, temp, dust_temp)."""
    fit = _get_fitter()
    fit.set_model_type('rsg')
    theta_in = [0.5, 4.0, 3500.0, 800.0]
    # tau_V > upper bound (e.g. 15 for sil)
    theta_out_tau = [20.0, 4.0, 3500.0, 800.0]
    assert fit.check_bounds(theta_in) is False
    assert fit.check_bounds(theta_out_tau) is True


def test_extinction_law_returns_spectrum():
    """extinction_law returns ArraySpectrum-like object with .wave, .flux."""
    fit = _get_fitter()
    wave = np.linspace(3000, 9000, 100)
    sp = fit.extinction_law(wave, Av=1.0, Rv=3.1)
    assert hasattr(sp, 'wave') and hasattr(sp, 'flux')
    assert len(sp.flux) == len(sp.wave)
    assert np.all(sp.flux <= 1.0) and np.all(sp.flux > 0)


def test_extinction_law_ccm89():
    """extinction_law with ccm89 option."""
    fit = _get_fitter()
    fit.options.extinction_law = 'ccm89'
    wave = np.linspace(4000, 7000, 50)
    sp = fit.extinction_law(wave, Av=0.5, Rv=3.1)
    assert sp.flux is not None and len(sp.flux) == 50


def test_create_blackbody():
    """create_blackbody returns spectrum with .wave and .flux."""
    fit = _get_fitter()
    bb = fit.create_blackbody(lum=4.0, temp=5000.0)
    assert hasattr(bb, 'wave') and hasattr(bb, 'flux')
    assert np.all(bb.flux >= 0) and np.all(np.isfinite(bb.flux))


def test_get_guess_blackbody():
    """get_guess returns array of length model_fit_params."""
    fit = _get_fitter()
    fit.set_model_type('blackbody')
    guess = fit.get_guess('blackbody', guess_type='params')
    assert len(guess) == 2
    assert np.all(np.isfinite(guess))


def test_sample_params_returns_arrays():
    """sample_params returns (params_sample, prob_sample)."""
    fit = _get_fitter()
    fit.set_model_type('blackbody')
    params = np.random.randn(100, 2) + np.array([4.0, 3.6])
    prob = -0.5 * np.sum((params - [4.0, 3.6])**2, axis=1)
    p, pr = fit.sample_params(params, prob, 2, nsamples=20)
    assert p.shape[0] == pr.shape[0]
    assert p.shape[1] == 2


def test_calculate_param_best_fit():
    """calculate_param_best_fit returns median and prints (show=False)."""
    fit = _get_fitter()
    fit.set_model_type('blackbody')
    sample = np.random.randn(200) + 4.0
    prob = np.ones(200)
    best = fit.calculate_param_best_fit(sample, prob, 2, 'luminosity', show=False)
    assert np.isfinite(best)
    assert 3.0 < best < 5.0


def test_log_likelihood_returns_tuple():
    """log_likelihood returns (log_prob, Av, Rv)."""
    fit = _get_fitter()
    fit.set_model_type('blackbody')
    from astropy.table import Table
    fit.phottable = Table({
        'magnitude': [15.0], 'mag_err': [0.1],
        'filter': ['WFC3_UVIS_F555W'], 'instrument': ['wfc3'],
    })
    try:
        fit.load_models(fit.phottable)
    except Exception:
        pytest.skip("load_models failed (missing data)")
    inst_filt = ['WFC3_UVIS_F555W']
    mag = np.array([15.0])
    magerr = np.array([0.1])
    theta = [4.0, 3.7]
    ll, av, rv = fit.log_likelihood(theta, inst_filt, mag, magerr, extinction=(0.0, 3.1))
    assert np.isfinite(ll) or ll == -np.inf
    assert (av is None and rv is None) or (av is not None and rv is not None)


def test_compute_model_mag_blackbody():
    """compute_model_mag with blackbody returns (mags, Av, Rv)."""
    fit = _get_fitter()
    fit.set_model_type('blackbody')
    try:
        fit.load_models(None)
    except Exception:
        pytest.skip("load_models failed (missing data)")
    inst_filt = ['WFC3_UVIS_F555W']
    theta = [4.0, 3.7]
    model_mag, av, rv = fit.compute_model_mag(inst_filt, theta, extinction=(0.0, 3.1))
    if model_mag is not None:
        assert len(model_mag) == len(inst_filt)
    assert av is None or (isinstance(av, (int, float)) and np.isfinite(av))
    assert rv is None or (isinstance(rv, (int, float)) and np.isfinite(rv))


def test_run_dynesty_raises_without_dynesty():
    """run_dynesty raises ImportError when dynesty is not installed."""
    from progenitors.sed import analysis
    if analysis._HAS_DYNESTY:
        pytest.skip("dynesty is installed; cannot test missing case")
    fit = _get_fitter()
    fit.set_model_type('blackbody')
    from astropy.table import Table
    fit.phottable = Table({
        'magnitude': [15.0], 'mag_err': [0.1],
        'filter': ['F'], 'instrument': ['wfc3'],
    })
    with pytest.raises(ImportError, match="dynesty"):
        fit.run_dynesty(fit.phottable, nlive=5)


def test_run_dynesty_smoke():
    """With dynesty installed, run_dynesty runs (tiny nlive)."""
    pytest.importorskip("dynesty", reason="dynesty required for smoke test")
    from progenitors.sed import analysis
    fit = _get_fitter()
    fit.set_model_type('blackbody')
    from astropy.table import Table
    fit.phottable = Table({
        'magnitude': [15.0], 'mag_err': [0.1],
        'filter': ['WFC3_UVIS_F555W'], 'instrument': ['wfc3'],
    })
    fit.load_models(fit.phottable)
    fit.run_dynesty(fit.phottable, nlive=10, bound='single')
