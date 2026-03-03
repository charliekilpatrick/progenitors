"""Example: test progenitors.sed.constants. Run from repo root: python examples/example_sed_constants.py"""
from progenitors.sed.constants import (
    color,
    black,
    red,
    blue,
    palette,
    SOLAR_LUM_CGS,
    BB_SCALE,
    PICKLES_CONST,
    alternate_names,
)

# color
c = color(255, 0, 0)
assert len(c) == 4 and c[0] == 1.0 and c[1] == 0.0
# named colors
assert len(black) == 4 and red[0] == 1.0 and blue[2] == 1.0
# palette
assert len(palette) >= 1
# physical constants
assert SOLAR_LUM_CGS > 0 and BB_SCALE > 0 and PICKLES_CONST > 0
# alternate_names
assert "magnitude" in alternate_names and "mag" in alternate_names["magnitude"]
print("progenitors.sed.constants OK")
