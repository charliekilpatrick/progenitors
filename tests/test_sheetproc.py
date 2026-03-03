"""Unit tests for sheetproc module."""
import os
import pytest

# Allow tests to run without real credentials
if "PROGENITORS_SHEET" not in os.environ:
    os.environ.setdefault("PROGENITORS_SHEET", "test-dummy")
if "YSE_USER" not in os.environ:
    os.environ.setdefault("YSE_USER", "test")
if "YSE_PASSWORD" not in os.environ:
    os.environ.setdefault("YSE_PASSWORD", "test")


def test_transpose():
    from progenitors.sheets.sheetproc import transpose
    data = [[1, 2, 3], [4, 5, 6]]
    out = transpose(data)
    assert out == [[1, 4], [2, 5], [3, 6]]


def test_params_cols():
    from progenitors.sheets.sheetproc import params
    assert "cols" in params
    assert "Name" in params["cols"]
    assert "SHEET" in params


def test_convert_table_to_lists():
    from astropy.table import Table
    from progenitors.sheets.sheetproc import convert_table_to_lists, params
    t = Table({"Name": ["SN1", "SN2"], "RA": ["01:00:00", "02:00:00"], "Dec": ["+00:00", "+01:00"]})
    out = convert_table_to_lists(t)
    assert out[0] == list(params["cols"])  # header row
    assert len(out) == 3  # header + 2 rows
    assert "SN1" in str(out[1])


def test_load_save_metadata_roundtrip(tmp_path):
    from astropy.table import Table
    from progenitors.sheets.sheetproc import load_metadata_from_file, save_metadata_to_file
    all_data = {"Type Ia": Table({"Name": ["SN1"], "RA": ["00:00"], "Dec": ["+00"]})}
    all_data["Type Ia"].meta = {"SN1": {"test": 1}}
    save_metadata_to_file(all_data, str(tmp_path))
    loaded = {"Type Ia": Table({"Name": ["SN1"], "RA": ["00:00"], "Dec": ["+00"]})}
    loaded = load_metadata_from_file(loaded, str(tmp_path))
    assert loaded["Type Ia"].meta.get("SN1", {}).get("test") == 1
