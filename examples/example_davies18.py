"""Example: test progenitors.davies18 package. Run from repo root: python examples/example_davies18.py"""
import os
import progenitors.davies18 as d18

assert d18 is not None
pkg_dir = os.path.dirname(d18.__file__)
assert os.path.isdir(pkg_dir)
files = os.listdir(pkg_dir)
csv = [f for f in files if f.endswith(".csv")]
pro = [f for f in files if f.endswith(".pro")]
assert len(csv) + len(pro) >= 1, "davies18 should contain .csv or .pro files"
print("progenitors.davies18 OK")
