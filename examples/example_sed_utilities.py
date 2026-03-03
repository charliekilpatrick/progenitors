"""Example: test progenitors.sed.utilities. Run from repo root: python examples/example_sed_utilities.py"""
from progenitors.sed.utilities import is_number, round_to_n, parse_coord

# is_number
assert is_number(1.0) and is_number("2.5") and not is_number("x")
# round_to_n
val = round_to_n(123.456, 4)
assert isinstance(val, (int, float)) and val > 0
# parse_coord (deg and sexagesimal)
c = parse_coord(10.0, 20.0)
assert c is not None and abs(c.ra.deg - 10.0) < 0.01
c2 = parse_coord("00:00:00", "+00:00:00")
assert c2 is not None
# array
ra, dec = [10.0, 20.0], [5.0, 10.0]
coords = parse_coord(ra, dec)
assert len(coords) == 2
print("progenitors.sed.utilities OK")
