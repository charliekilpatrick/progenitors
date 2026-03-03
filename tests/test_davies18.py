"""Smoke tests for progenitors.davies18 package. See tests/davies18/ for full unit tests."""
import os
import pytest


def test_davies18_package_import():
    import progenitors.davies18 as d18
    assert d18 is not None
    assert hasattr(d18, "lfunc")
    assert hasattr(d18, "prog_mc")
    assert hasattr(d18, "io_utils")


def test_davies18_has_data_dir():
    """Davies18 package has data/ and legacy/ layout."""
    import progenitors.davies18
    pkg_dir = os.path.dirname(progenitors.davies18.__file__)
    assert os.path.isdir(os.path.join(pkg_dir, "data"))
    assert os.path.isdir(os.path.join(pkg_dir, "legacy"))
