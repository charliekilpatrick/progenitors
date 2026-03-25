# Progenitors

[![CI](https://github.com/charliekilpatrick/progenitors/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/charliekilpatrick/progenitors/actions/workflows/ci.yml) [![Documentation](https://github.com/charliekilpatrick/progenitors/actions/workflows/docs.yml/badge.svg?branch=main)](https://charliekilpatrick.github.io/progenitors/)

Methods for downloading, organizing, and analyzing data for progenitor systems of explosive transients (e.g., supernovae). Includes a Google Sheets–driven pipeline, SED fitting from photometry, and the Davies et al. (2018) progenitor analysis (luminosity function and mass MC) in Python.

**Version:** `progenitors.__version__` (single source: `progenitors._version`).  
**Python:** 3.12+

---

## In this repository


| Document                                                                                             | Description                                                                                            |
| ---------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------ |
| **[examples/README.md](examples/README.md)**                                                         | Example scripts for each major module; run from repo root.                                             |
| **[docs/README.md](docs/README.md)**                                                                 | How to build Sphinx documentation locally.                                                             |
| **[progenitors/davies18/README.md](progenitors/davies18/README.md)**                                 | Davies et al. (2018) citation, data layout, and Python usage (lfunc, prog_mc, distribution, plotting). |
| **[progenitors/davies18/README_PYTHON_REFACTOR.md](progenitors/davies18/README_PYTHON_REFACTOR.md)** | Refactor notes and validation for the IDL → Python port.                                               |
| **[CONTRIBUTING.md](CONTRIBUTING.md)**                                                                 | Submitting issues, contact information, and how to contribute.                                          |


---

## Installation

From the **repository root**:

```bash
pip install -e .
pip install -e ".[dev]"   # pytest, sphinx, sphinx-rtd-theme
```

**Optional extras:** `[sed]`, `[sed_hst]` (stsynphot), `[sed_fit]` (emcee, dynesty).

**Conda:** `conda env create -f environment.yml && conda activate progenitors && pip install -e .`

**Verify:** `progenitors --help` and `pytest tests/ -v`

---

## Configuration

Credentials and API keys use **environment variables** only. See `progenitors.env_config` and `progenitors.settings.env_config`.


| Purpose      | Variables                                                         |
| ------------ | ----------------------------------------------------------------- |
| Google Sheet | `PROGENITORS_SHEET`                                               |
| OAuth        | `credentials.json` at repo root or `PROGENITORS_CREDENTIALS_PATH` |
| YSE-PZ       | `YSE_USER`, `YSE_PASSWORD`                                        |
| TNS          | `TNS_API_KEY`, `TNS_BOT_NAME`, `TNS_BOT_ID`                       |
| ADS          | `ADS_AUTHCODE`                                                    |
| Email alerts | `GMAIL_LOGIN`, `GMAIL_PASSWORD`, `MY_EMAIL`                       |
| Pipeline log level | `PROGENITORS_LOG_LEVEL` (`DEBUG`, `INFO`, `WARNING`, `ERROR`; default `INFO`) or CLI `--log-level`. Set `NO_COLOR=1` to disable ANSI colors. |
| Metadata compression | `PROGENITORS_METADATA_ZSTD_LEVEL` — zstd level for `.pkl.zst` cache (default `1` = fast; use `3`–`6` for smaller files). |
| Metadata fsync | `PROGENITORS_METADATA_FSYNC=0` — skip `fsync` after metadata writes (faster on cloud folders; slightly less crash-safe). |
| Metadata cache directory | `PROGENITORS_METADATA_DIR` — override cache path (use a local SSD path for faster load/save than cloud-synced folders). |
| Corrupt main cache | `PROGENITORS_METADATA_QUARANTINE_CORRUPT=0` — do not rename a bad `*.pkl.zst` to `*.corrupt` after loading from backup (default is to quarantine). |
| Google Sheets writes | `PROGENITORS_SHEETS_WRITE_INTERVAL` — minimum seconds between write calls (default `1.05`, ~60/min quota). Retries with backoff on HTTP 429/503. |


Example: `cp progenitors/scripts/set_env_example.sh my_env.sh`, edit, then `source my_env.sh`.

Pipeline logs go to **stderr** in a compact form (`HH:MM:SS LEVEL module:message`). On a TTY, only the level word (**INFO**, **WARNING**, **ERROR**, **CRITICAL**) is colored (unless `NO_COLOR` is set). Per-object metadata progress is **DEBUG** only; **INFO** summarizes each metadata pass and sheet step.

---

## Package layout


| Path                                   | Description                                                                                                                                                                                                                                                                                                                                  |
| -------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **progenitors/**                       | Main package.                                                                                                                                                                                                                                                                                                                                |
| **pipeline**, **options**              | CLI `progenitors`: sync Sheets, fetch metadata (TNS, YSE, NED, ADS, OSC, MAST), optional alerts.                                                                                                                                                                                                                                             |
| **sheets**                             | Google Sheets I/O, column definitions.                                                                                                                                                                                                                                                                                                       |
| **util**, **env_config**, **settings** | Helpers, env validation, pipeline/sheet/SED settings.                                                                                                                                                                                                                                                                                        |
| **sed**                                | SED fitting (analysis, dust, constants, synphot_compat).                                                                                                                                                                                                                                                                                     |
| **davies18**                           | Davies et al. (2018): lfunc (luminosity function), prog_mc (mass MC), **distribution** (mass-distribution and χ²: `generate_distribution`, `generate_sample`, `calculate_chi2`, `run_mass_simulation`), plotting. Data in `davies18/data/`, IDL in `davies18/legacy/`. See [progenitors/davies18/README.md](progenitors/davies18/README.md). |
| **extinction**                         | Extinction data and placeholders.                                                                                                                                                                                                                                                                                                            |
| **scripts**                            | Env examples, SED helpers.                                                                                                                                                                                                                                                                                                                   |


---

## Usage

**Google Sheets pipeline — entry point:** After install, the console script runs `progenitors.__main__:main`, which parses CLI flags in `progenitors.options` and calls `progenitors.pipeline.main` (download sheet → merge YSE → fetch metadata → optional reclassification → upload). Equivalent invocations:

```bash
python -m progenitors
progenitors
```

**Pipeline:** Set env vars and `credentials.json`, then:

```bash
progenitors
progenitors --redo --metadata ned,distance,tns,yse,ads,hst
progenitors --update-classification --alert
```

Refresh TNS object/classification metadata, re-sort sheets by classification, trim stale metadata keys, and email on new candidates (similar to legacy `progenitors.py --update-tns-class --update-classification --trim-metadata --alert`):

```bash
progenitors --update-tns-class --update-classification --trim-metadata --alert
```

(`--update-tns-class` ensures `tns` is included in the metadata types for this run; the default metadata set already includes `tns`, so this flag matters most when you pass a custom `--metadata` list that omitted it.)

**Local metadata cache:** zstd-compressed pickle (`.pkl.zst`) under [`progenitors/metadata/`](progenitors/metadata/README.md) by default, or `PROGENITORS_METADATA_DIR` if set; legacy `.pkl` is still read until a save migrates it. Gitignored.

**Python (sheets and metadata):**

```python
from progenitors import params, parse_arguments
from progenitors.sheets import download_progenitor_data, load_metadata_from_file
from progenitors.env_config import validate_env
validate_env(["sheets", "yse"])
sndata = download_progenitor_data(params["SHEET"])
sndata = load_metadata_from_file(sndata, params["metadata"])
```

**Mass-distribution and Davies et al. (2018):**

```python
from progenitors.davies18 import distribution, lfunc, prog_mc

masses, probs = distribution.generate_distribution(10.0, 15.0)
samples = distribution.generate_sample(masses, probs, 100)
chi2 = distribution.calculate_chi2(obs_masses, sim_masses, lims)

BESTFITPARS, LPOBS, LPMOD = lfunc.mainproc(NTRIALS_OBS=1000, NTRIALS_MOD=1000, seed=42)
result = prog_mc.run_prog_mc(obsfile="IIPprog_obsdata_2019.csv", MLfile="M-L_STARS.txt", ntrials=1000, seed=42)
```

Details, citation, and data layout: [progenitors/davies18/README.md](progenitors/davies18/README.md).

---

## Examples and tests

**Examples:** From repo root, run scripts under `examples/` (package must be installed). Full list and descriptions: [examples/README.md](examples/README.md).

```bash
python examples/example_extinction.py
# or run all:
for f in examples/example_*.py; do python "$f" || exit 1; done
```

**Tests:** `pytest`, or e.g. `pytest tests/davies18/ -v`, `pytest --cov=progenitors --cov-report=term-missing`. Markers: `slow`, `remote` (deselect with `-m "not slow and not remote"` for fast local runs).

---

## Documentation

Build Sphinx docs from repo root: `pip install -e ".[dev]"` then `sphinx-build -b html docs docs/_build`. Full instructions and theme options: [docs/README.md](docs/README.md).

---

## Issues, contributing, and contact

- **Submitting issues:** Use [GitHub Issues](https://github.com/YOUR_ORG/progenitors/issues) for bugs, feature requests, or documentation improvements. Include steps to reproduce (for bugs), your environment (OS, Python version), and avoid posting credentials. See [CONTRIBUTING.md](CONTRIBUTING.md) for detailed guidelines.
- **Contributing:** Fork the repo, create a branch, install with `pip install -e ".[dev]"`, run tests (`pytest -m "not slow and not remote" -v`), and open a pull request. Keep changes focused and update docs as needed. Full instructions: [CONTRIBUTING.md](CONTRIBUTING.md).
- **Contact:** **Charlie Kilpatrick** — [ckilpatrick@northwestern.edu](mailto:ckilpatrick@northwestern.edu). For sensitive or security-related matters, contact by email rather than opening a public issue.

Replace `YOUR_ORG/progenitors` in the Issues link with your repository URL if different.