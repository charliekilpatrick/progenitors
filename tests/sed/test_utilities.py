"""Unit tests for progenitors.sed.utilities."""
import numpy as np
import pytest
from astropy.table import Column


def test_is_number():
    from progenitors.sed.utilities import is_number
    assert is_number(1.0) is True
    assert is_number("2.5") is True
    assert is_number("x") is False


def test_round_to_n():
    from progenitors.sed.utilities import round_to_n
    # sed.utilities round_to_n returns float with %.Ng formatting
    val = round_to_n(123.456, 4)
    assert isinstance(val, (int, float))
    assert val > 0


def test_parse_coord_deg():
    from progenitors.sed.utilities import parse_coord
    c = parse_coord(10.0, 20.0)
    assert c is not None
    assert abs(c.ra.deg - 10.0) < 0.01


def test_parse_coord_sexagesimal():
    from progenitors.sed.utilities import parse_coord
    c = parse_coord("00:00:00", "+00:00:00")
    assert c is not None


def test_parse_coord_array():
    from progenitors.sed.utilities import parse_coord
    ra = [10.0, 20.0]
    dec = [5.0, 10.0]
    coords = parse_coord(ra, dec)
    assert len(coords) == 2
