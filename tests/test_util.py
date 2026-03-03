"""Unit tests for util module (root)."""
import os
import numpy as np
import pytest
from astropy.table import Table

# Set minimal env so util can be imported (it reads PROGENITORS_SHEET, etc.)
for key in ("PROGENITORS_SHEET", "PROGENITORS_SHEET_TEST", "TNS_API_KEY", "TNS_BOT_NAME",
            "TNS_BOT_ID", "YSE_USER", "YSE_PASSWORD", "ADS_AUTHCODE", "GMAIL_LOGIN",
            "GMAIL_PASSWORD", "MY_EMAIL"):
    os.environ.setdefault(key, "test-dummy")


def test_is_number():
    from progenitors.util import is_number
    assert is_number(1.0) is True
    assert is_number("1.5") is True
    assert is_number(0) is True
    assert is_number("abc") is False
    assert is_number(None) is False
    assert is_number(np.nan) is False


def test_check_dict():
    from progenitors.util import check_dict
    d = {"a": {"b": {"c": 1}}}
    assert check_dict(d, ["a", "b", "c"]) == 1
    assert check_dict(d, ["a", "b"]) == {"c": 1}
    assert check_dict(d, ["a", "x"]) is None
    assert check_dict(d, ["z"]) is None
    assert check_dict(None, ["a"]) is None


def test_round_to_n():
    from progenitors.util import round_to_n
    assert round_to_n(123.456, 2) == "123.46"
    assert round_to_n(0.00123, 2).replace(" ", "") == "0.00"


def test_parse_coord_deg():
    from progenitors.util import parse_coord
    coord = parse_coord(10.5, 20.0)
    assert coord is not None
    assert abs(coord.ra.deg - 10.5) < 0.01
    assert abs(coord.dec.deg - 20.0) < 0.01


def test_parse_coord_sexagesimal():
    from progenitors.util import parse_coord
    coord = parse_coord("00:00:00", "+00:00:00")
    assert coord is not None
    assert abs(coord.ra.deg) < 0.01
    assert abs(coord.dec.deg) < 0.01


def test_parse_coord_array():
    from progenitors.util import parse_coord
    ra = [10.0, 20.0]
    dec = [5.0, 10.0]
    coords = parse_coord(ra, dec)
    assert coords is not None
    assert len(coords) == 2


def test_format_date():
    from progenitors.util import format_date
    out = format_date("2020-01-15 12:30:00")
    assert "2020" in out
    assert "01" in out
    assert "15" in out


def test_get_tns_header():
    from progenitors.util import get_tns_header
    # Just check it returns a dict with User-Agent
    h = get_tns_header()
    assert isinstance(h, dict)
    assert "User-Agent" in h
    assert "tns_marker" in h["User-Agent"]


def test_get_coord_from_row():
    from progenitors.util import get_coord
    row = Table({"RA": ["01:36:48.161"], "Dec": ["+15:45:31.00"]})
    coord = get_coord(row[0])
    assert coord is not None


def test_get_coord_dict_row():
    from progenitors.util import get_coord
    class DictRow:
        def keys(self):
            return ["RA", "Dec"]
        def __getitem__(self, k):
            return "01:36:48.161" if k == "RA" else "+15:45:31.00"
    coord = get_coord(DictRow())
    assert coord is not None


def test_get_classification_mask_minimal_table():
    from progenitors.util import get_classification_mask
    t = Table({
        "Classification": ["SN II"],
        "Name": ["SN2020test"],
        "Discovery Date": ["2020-01-01"],
        "Redshift": [""],
    })
    t.meta = {}
    mask = get_classification_mask(t)
    assert mask is not None
    assert len(mask) == 1
    assert mask.dtype == bool
