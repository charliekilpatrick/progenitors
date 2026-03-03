"""Example: test progenitors package and entry points. Run from repo root: python examples/example_package.py"""
from progenitors import params, parse_arguments, message
from progenitors.__main__ import main
from progenitors.pipeline import main as pipeline_main, do_trim_metadata

assert "SHEET" in params or "metadata" in params
assert "target" in params
assert callable(parse_arguments) and callable(message)
assert callable(main) and callable(pipeline_main) and callable(do_trim_metadata)
print("progenitors package OK")
