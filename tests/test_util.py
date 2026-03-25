"""Unit tests for util module (root)."""
import copy
import io
import os
import zipfile

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


def _sample_tns_public_objects_zip_bytes():
    csv = (
        "2025-03-24 00:00:00\n"
        "objid,name_prefix,name,ra,declination\n"
        "1,SN,SN2025a,10.0,20.0\n"
    )
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as zf:
        zf.writestr("tns_public_objects.csv", csv.encode("utf-8"))
    return buf.getvalue()


def test_parse_tns_public_objects_zip_bytes():
    from progenitors.util import parse_tns_public_objects_zip_bytes

    t = parse_tns_public_objects_zip_bytes(_sample_tns_public_objects_zip_bytes())
    assert len(t) == 1
    assert "name" in t.colnames
    assert t["name"][0] == "SN2025a"
    assert t["objid"][0] == 1


def test_parse_tns_public_objects_zip_bytes_no_csv_member():
    from progenitors.util import parse_tns_public_objects_zip_bytes

    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as zf:
        zf.writestr("readme.txt", b"x")
    with pytest.raises(ValueError, match="No .csv"):
        parse_tns_public_objects_zip_bytes(buf.getvalue())


def test_fetch_tns_public_objects_zip_bytes_requires_api_key():
    import progenitors.util as u
    from progenitors.util import fetch_tns_public_objects_zip_bytes

    saved = copy.deepcopy(u.params["tns"])
    try:
        u.params["tns"] = {"api_key": "", "bot_id": "1", "bot_name": "bot"}
        with pytest.raises(ValueError, match="TNS_API_KEY"):
            fetch_tns_public_objects_zip_bytes()
    finally:
        u.params["tns"] = saved


def test_fetch_tns_public_objects_zip_bytes_requires_bot_marker_fields():
    import progenitors.util as u
    from progenitors.util import fetch_tns_public_objects_zip_bytes

    saved = copy.deepcopy(u.params["tns"])
    try:
        u.params["tns"] = {"api_key": "k", "bot_id": "", "bot_name": "bot"}
        with pytest.raises(ValueError, match="TNS_BOT_ID"):
            fetch_tns_public_objects_zip_bytes()
    finally:
        u.params["tns"] = saved


def test_normalize_spectral_classification():
    from progenitors.util import normalize_spectral_classification

    assert normalize_spectral_classification("") == ""
    assert normalize_spectral_classification("Ia") == "SN Ia"
    assert normalize_spectral_classification("II P") == "SN II-P"


def test_apply_tns_public_catalog_to_sndata_name_only():
    from astropy.table import Table

    from progenitors.util import apply_tns_public_catalog_to_sndata

    tns = Table(
        rows=[
            ("SN 2024AA", "Ia", ""),
            ("AT 2024BB", "II P", "AT2024BB, SN 2024BB"),
        ],
        names=("name", "type", "internal_names"),
    )
    sheet = Table(
        {
            "Name": ["2024AA", "2024BB", "nomatch"],
            "RA": ["01:00:00", "12:00:00.0", "00:00:01"],
            "Dec": ["+10:00:00", "+00:00:00", "+89:00:00"],
            "Classification": ["Candidate", "Candidate", "Candidate"],
        }
    )
    sheet.meta = {}
    sndata = {"Type Ia": sheet}
    apply_tns_public_catalog_to_sndata(sndata, tns)

    assert sheet["Classification"][0] == "SN Ia"
    assert sheet["Classification"][1] == "SN II-P"
    assert sheet["Classification"][2] == "Candidate"
    assert sheet.meta["2024AA"]["tns"]["object_type"]["name"] == "Ia"


def test_prepare_tns_public_crossmatch_missing_column():
    from astropy.table import Table

    from progenitors.util import prepare_tns_public_crossmatch

    bad = Table({"name": ["a"]})
    with pytest.raises(ValueError, match="missing"):
        prepare_tns_public_crossmatch(bad)


@pytest.mark.remote
def test_fetch_tns_public_objects_table_live():
    """Download real TNS archive when credentials are set (not dummy test values)."""
    key = (os.environ.get("TNS_API_KEY") or "").strip()
    bot = (os.environ.get("TNS_BOT_NAME") or "").strip()
    bid = (os.environ.get("TNS_BOT_ID") or "").strip()
    if not key or not bot or not bid:
        pytest.skip("TNS_API_KEY, TNS_BOT_NAME, TNS_BOT_ID not all set")
    if key == "test-dummy" or bot == "test-dummy":
        pytest.skip("TNS env vars are placeholder test-dummy values")

    from progenitors.util import fetch_tns_public_objects_table

    t = fetch_tns_public_objects_table(timeout=600)
    cols = {c.lower() for c in t.colnames}
    assert "name" in cols
    assert "objid" in cols
    assert len(t) > 1000
