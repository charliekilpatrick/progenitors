"""Unit tests for progenitors.davies18.plotting (optional matplotlib)."""
import os
import numpy as np
import pytest

matplotlib = pytest.importorskip("matplotlib", reason="matplotlib optional for plotting")

from progenitors.davies18 import plotting


def test_plot_lum_contour_smoke():
    llo = np.linspace(4.1, 4.5, 5)
    lhi = np.linspace(5.1, 5.6, 6)
    im_c = np.random.rand(5, 6).cumsum(axis=0).cumsum(axis=1)
    lev_c = np.array([1.0, 2.0, 3.0])
    with open(os.devnull, "w") as devnull:
        pass
    out = os.path.join(os.path.dirname(__file__), "..", "tmp_LUM_test.eps")
    try:
        plotting.plot_lum_contour(im_c, llo, lhi, lev_c=lev_c, outpath=out)
        assert os.path.isfile(out)
    finally:
        if os.path.isfile(out):
            os.remove(out)


def test_plot_all_from_compare_no_file():
    # Should not raise if path does not exist (just returns)
    plotting.plot_all_from_compare(contour_npz_path="/nonexistent/LUM_CONTOUR.npz")


def test_plot_all_from_compare_with_file():
    # If legacy LUM_CONTOUR.npz exists, should produce fig
    from progenitors.davies18 import io_utils
    leg = os.path.join(os.path.dirname(io_utils.__file__), "legacy", "LUM_CONTOUR.npz")
    if not os.path.isfile(leg):
        pytest.skip("LUM_CONTOUR.npz not found; run lfunc.mainproc first")
    plotting.plot_all_from_compare(contour_npz_path=leg)
    figs_dir = os.path.join(os.path.dirname(io_utils.__file__), "figs")
    assert os.path.isdir(figs_dir)
    eps = os.path.join(figs_dir, "LUM_CONTOUR.eps")
    assert os.path.isfile(eps)
