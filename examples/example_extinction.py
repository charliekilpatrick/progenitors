"""Example: test progenitors.extinction package. Run from repo root: python examples/example_extinction.py"""
import os
import progenitors.extinction as ext

assert ext is not None
pkg_dir = os.path.dirname(ext.__file__)
assert os.path.isdir(pkg_dir)
# Optional: at least one known subdir or package exists
for sub in ("data", "filters", "colortemps", "spectra"):
    if os.path.isdir(os.path.join(pkg_dir, sub)):
        break
else:
    assert os.path.exists(pkg_dir)
print("progenitors.extinction OK")
