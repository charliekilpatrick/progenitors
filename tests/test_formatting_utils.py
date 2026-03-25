"""Tests for shared formatting_utils (no full util import required)."""
import numpy as np
import pytest
from astropy.table import Column

from progenitors.formatting_utils import (
    is_number,
    parse_coord,
    round_to_decimal_places,
    round_to_significant_figures,
)


def test_is_number_unified():
    assert is_number(1.0) is True
    assert is_number("1.5") is True
    assert is_number(0) is True
    assert is_number("abc") is False
    assert is_number(None) is False
    assert is_number(np.nan) is False


def test_round_to_decimal_places():
    assert round_to_decimal_places(123.456, 2) == "123.46"
    assert round_to_decimal_places(0.00123, 2).replace(" ", "") == "0.00"


def test_round_to_significant_figures():
    val = round_to_significant_figures(123.456, 4)
    assert isinstance(val, (int, float))
    assert val > 0


def test_parse_coord_deg():
    c = parse_coord(10.5, 20.0)
    assert c is not None
    assert abs(c.ra.deg - 10.5) < 0.01


def test_parse_coord_array():
    ra = [10.0, 20.0]
    dec = [5.0, 10.0]
    coords = parse_coord(ra, dec)
    assert len(coords) == 2


def test_parse_coord_column():
    ra = Column([10.0, 20.0])
    dec = Column([5.0, 10.0])
    coords = parse_coord(ra, dec)
    assert len(coords) == 2
