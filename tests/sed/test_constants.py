"""Unit tests for progenitors.sed.constants."""
import pytest

pytest.importorskip("astropy")


def test_color():
    from progenitors.sed.constants import color
    c = color(255, 0, 0)
    assert len(c) == 4
    assert c[0] == 1.0
    assert c[1] == 0.0
    assert c[2] == 0.0


def test_named_colors():
    from progenitors.sed.constants import black, red, blue, green
    assert len(black) == 4
    assert red[0] == 1.0
    assert blue[2] == 1.0


def test_palette():
    from progenitors.sed.constants import palette
    assert len(palette) >= 1


def test_physical_constants():
    from progenitors.sed.constants import SOLAR_LUM_CGS, BB_SCALE, PICKLES_CONST
    assert SOLAR_LUM_CGS > 0
    assert BB_SCALE > 0
    assert PICKLES_CONST > 0


def test_alternate_names():
    from progenitors.sed.constants import alternate_names
    assert "magnitude" in alternate_names
    assert "mag" in alternate_names["magnitude"]
