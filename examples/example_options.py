"""Example: test progenitors.options. Run from repo root: python examples/example_options.py"""
import sys
from progenitors.options import parse_coord, parse_two_floats, parse_arguments, message

# parse_coord
c = parse_coord("10.5", "20.0")
assert c is not None and abs(c.ra.deg - 10.5) < 0.01
# parse_two_floats
a, b = parse_two_floats("1.0 2.0")
assert a == 1.0 and b == 2.0
# parse_arguments (use minimal argv)
sys.argv = ["progenitors"]
args = parse_arguments()
assert args.redo is False and args.alert is False
# message prints to stdout
message("example_options check")
print("progenitors.options OK")
