# Examples

Scripts to exercise main modules. Run from **repository root** with the package installed (`pip install -e .`). No live credentials or network unless noted.

```bash
python examples/example_util.py
# or run all:
for f in examples/example_*.py; do python "$f" || exit 1; done
```

| Script | Module | Description |
|--------|--------|-------------|
| `example_util.py` | `progenitors.util` | `is_number`, `parse_coord`, `round_to_n`, etc. |
| `example_options.py` | `progenitors.options` | `parse_arguments`, `message` |
| `example_env_config.py` | `progenitors.env_config` | `get_env`, `validate_env` |
| `example_sheets.py` | `progenitors.sheets` | `transpose`, `params`, load/save metadata |
| `example_sed_constants.py` | `progenitors.sed.constants` | Colors, physical constants |
| `example_sed_utilities.py` | `progenitors.sed.utilities` | `is_number`, `round_to_n`, `parse_coord` |
| `example_extinction.py` | `progenitors.extinction` | Package import |
| `example_davies18.py` | `progenitors.davies18` | Package import and data paths |
| `example_simulation.py` | `progenitors.davies18.distribution` | `generate_distribution`, `generate_sample`, `calculate_chi2` |
| `example_package.py` | `progenitors` | `params`, `parse_arguments`, pipeline main |
