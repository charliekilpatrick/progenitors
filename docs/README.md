# Building the documentation

From the **repository root** (see main [README.md](../README.md) for installation):

```bash
pip install -e ".[dev]"   # includes sphinx, sphinx-rtd-theme
sphinx-build -b html docs docs/_build
# open docs/_build/index.html
```

To include API docs for `progenitors.sed.analysis` and `progenitors.sed.dust`, install optional deps so those modules import: `pip install -e ".[sed_fit]"` (emcee, dynesty). Dust module also needs scipy and data files in `progenitors/sed/data/dust/`.

Or from inside `docs/`:

```bash
make html
# open _build/html/index.html
```

To use the Read the Docs theme, set in `docs/conf.py`:

```python
html_theme = "sphinx_rtd_theme"
```

and install `sphinx-rtd-theme`.
