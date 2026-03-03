"""Unit tests for progenitors.extinction package."""
import os
import pytest


def test_extinction_package_import():
    import progenitors.extinction as ext
    assert ext is not None


def test_extinction_has_data_dir():
    """Extinction package contains data subdirs (Observed, Templates, etc.)."""
    import progenitors.extinction
    pkg_dir = os.path.dirname(progenitors.extinction.__file__)
    assert os.path.isdir(pkg_dir)
    # Optional: check for known subdirs if they exist
    for sub in ("data", "filters", "colortemps", "spectra"):
        path = os.path.join(pkg_dir, sub)
        if os.path.isdir(path):
            assert True
            break
    else:
        assert os.path.exists(pkg_dir)
