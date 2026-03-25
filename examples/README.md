# Examples

Scripts to exercise main modules. Run from the **repository root** with the package installed (`pip install -e .` or `pip install -e ".[dev]"`). No live credentials or network unless noted. See the main [README.md](../README.md) for installation and package layout.

```bash
python examples/example_extinction.py
# or run all:
for f in examples/example_*.py; do python "$f" || exit 1; done
```

| Script | Module | Description |
|--------|--------|-------------|
| `example_sheets.py` | `progenitors.sheets` | `transpose`, `params`, load/save metadata |
| `example_extinction.py` | `progenitors.extinction` | Package import |
| `example_davies18.py` | `progenitors.davies18` | Package import and data paths |
| `example_simulation.py` | `progenitors.davies18.distribution` | `generate_distribution`, `generate_sample`, `calculate_chi2` |
