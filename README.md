# Progenitors

Methods for downloading, organizing, and analyzing data for progenitor systems of explosive transients (e.g., supernovae). Includes a Google Sheets–driven pipeline, SED fitting from photometry, and the Davies et al. (2018) progenitor analysis (luminosity function and mass MC) in Python.

**Version:** `progenitors.__version__` (single source: `progenitors._version`).  
**Python:** 3.12+

---

## Installation

From repository root:

```bash
pip install -e .
pip install -e ".[dev]"   # pytest, sphinx
```

Optional: `[sed]`, `[sed_hst]` (stsynphot), `[sed_fit]` (emcee, dynesty).

Conda: `conda env create -f environment.yml && conda activate progenitors && pip install -e .`

**Check:** `progenitors --help` and `pytest tests/ -v`

---

## Configuration

Credentials and API keys come from **environment variables** only. See `progenitors.env_config` and `progenitors.settings.env_config`.

| Purpose      | Variables |
|-------------|-----------|
| Google Sheet| `PROGENITORS_SHEET` |
| OAuth       | `credentials.json` at repo root or `PROGENITORS_CREDENTIALS_PATH` |
| YSE-PZ      | `YSE_USER`, `YSE_PASSWORD` |
| TNS         | `TNS_API_KEY`, `TNS_BOT_NAME`, `TNS_BOT_ID` |
| ADS         | `ADS_AUTHCODE` |
| Email alerts| `GMAIL_LOGIN`, `GMAIL_PASSWORD`, `MY_EMAIL` |

Example: `cp progenitors/scripts/set_env_example.sh my_env.sh` then edit and `source my_env.sh`.

---

## Package layout

| Path | Description |
|------|-------------|
| **progenitors/** | Main package. |
| **pipeline**, **options** | CLI `progenitors`: sync Sheets, fetch metadata (TNS, YSE, NED, ADS, OSC, MAST), optional alerts. |
| **sheets** | Google Sheets I/O, column definitions. |
| **util**, **env_config**, **settings** | Helpers, env validation, pipeline/sheet/SED settings. |
| **sed** | SED fitting (analysis, dust, constants, synphot_compat). |
| **davies18** | Davies et al. (2018): lfunc (luminosity function), prog_mc (mass MC), **distribution** (mass-distribution and χ²: `generate_distribution`, `generate_sample`, `calculate_chi2`, `run_mass_simulation`), plotting. Data in `davies18/data/`, IDL in `davies18/legacy/`. |
| **extinction** | Extinction data and placeholders. |
| **scripts** | Env examples, SED helpers. |

---

## Usage

**Pipeline:** Set env vars and `credentials.json`, then:

```bash
progenitors
progenitors --redo --metadata ned,distance,tns,yse,ads,hst
progenitors --update-classification --alert
```

**Python:**

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

See `progenitors/davies18/README.md` for citation and data layout.

---

## Examples and tests

**Examples** (from repo root): `python examples/example_util.py` etc. List in `examples/README.md`.

**Tests:** `pytest`, `pytest tests/davies18/ -v`, `pytest --cov=progenitors --cov-report=term-missing`.

---

## Documentation

`pip install -e ".[dev]"` then `sphinx-build -b html docs docs/_build`. See `docs/README.md`.

---

## Contact

Charlie Kilpatrick — ckilpatrick@northwestern.edu
