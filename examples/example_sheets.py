"""Example: test progenitors.sheets (no live Google Sheets). Run from repo root."""
from astropy.table import Table
from progenitors.sheets.sheetproc import (
    transpose,
    params,
    convert_table_to_lists,
    load_metadata_from_file,
    save_metadata_to_file,
)
import tempfile
import os

# transpose
data = [[1, 2, 3], [4, 5, 6]]
assert transpose(data) == [[1, 4], [2, 5], [3, 6]]
# params
assert "cols" in params and "Name" in params["cols"]
# convert_table_to_lists
t = Table({"Name": ["SN1", "SN2"], "RA": ["01:00:00", "02:00:00"], "Dec": ["+00:00", "+01:00"]})
out = convert_table_to_lists(t)
assert out[0] == list(params["cols"]) and len(out) == 3
# load/save metadata roundtrip
all_data = {"Type Ia": Table({"Name": ["SN1"], "RA": ["00:00"], "Dec": ["+00"]})}
all_data["Type Ia"].meta = {"SN1": {"test": 1}}
with tempfile.TemporaryDirectory() as tmp:
    save_metadata_to_file(all_data, tmp)
    loaded = {"Type Ia": Table({"Name": ["SN1"], "RA": ["00:00"], "Dec": ["+00"]})}
    loaded = load_metadata_from_file(loaded, tmp)
    assert loaded["Type Ia"].meta.get("SN1", {}).get("test") == 1
print("progenitors.sheets OK")
