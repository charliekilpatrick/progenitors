"""Example: test progenitors.util (no network). Run from repo root: python examples/example_util.py"""
from progenitors.util import (
    is_number,
    check_dict,
    round_to_n,
    parse_coord,
    format_date,
    get_tns_header,
)

# is_number
assert is_number(1.0) and is_number("1.5") and not is_number("abc")
# check_dict
d = {"a": {"b": {"c": 1}}}
assert check_dict(d, ["a", "b", "c"]) == 1 and check_dict(d, ["a", "x"]) is None
# round_to_n
assert round_to_n(123.456, 2) == "123.46"
# parse_coord (deg and sexagesimal)
c = parse_coord(10.5, 20.0)
assert c is not None and abs(c.ra.deg - 10.5) < 0.01
c2 = parse_coord("00:00:00", "+00:00:00")
assert c2 is not None
# format_date
assert "2020" in format_date("2020-01-15 12:30:00")
# get_tns_header (no request sent)
h = get_tns_header()
assert isinstance(h, dict) and "User-Agent" in h
print("progenitors.util OK")
