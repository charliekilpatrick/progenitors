"""Unit tests for progenitors.sed.plotting (truncate, sed_plot, axis_titles)."""
import pytest

pytest.importorskip("matplotlib", reason="plotting requires matplotlib")


def test_truncate():
    """truncate rounds to given decimals."""
    from progenitors.sed.plotting import truncate
    assert truncate(1.2345, 0) == 1.0
    assert truncate(1.2345, 1) == 1.2
    assert truncate(1.2345, 2) == 1.23
    assert truncate(1.999, 0) == 1.0


def test_sed_plot_init():
    """sed_plot initializes with fitter and axis_titles."""
    from progenitors.sed.plotting import sed_plot
    p = sed_plot()
    assert hasattr(p, 'fit')
    assert hasattr(p, 'axis_titles')
    assert hasattr(p, 'plot_types')
    assert 'wavelength' in p.axis_titles
    assert 'flux' in p.axis_titles
    assert p.plot_types.get('sed') is not None
    assert p.plot_types.get('hr') is not None


def test_sed_plot_axis_titles_keys():
    """Axis titles include expected keys for SED/HR/CMD plots."""
    from progenitors.sed.plotting import sed_plot
    p = sed_plot()
    for key in ('wavelength', 'flux', 'log_T', 'log_L', 'A_V', 'Teff'):
        assert key in p.axis_titles
