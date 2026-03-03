# Davies et al. IDL → Python refactor (notes)

**Main documentation and citation:** see [README.md](README.md) in this directory.

## Layout

- **legacy/** — All original IDL code (`.pro`) and reference outputs (`.sav`, `.eps`, `.tex`, `.out`).
- **data/** — Input data: `IIPprog_obsdata_2019.csv`, `IIPprog_obsdata_2022.csv`, `M-L_STARS.txt` (synthetic M–L relation if `M-L_S09.txt` is not provided), `input_masses.txt` (for distribution MC).
- **Python modules** — `histutils`, `powerlaw`, `m_hi_lo_fit`, `io_utils`, `lfunc`, `prog_mc`, `distribution`, `plotting`. Run `mainproc()` for the luminosity-function pipeline; `run_prog_mc()` for the mass MC; use `plotting` to write EPS figures to `figs/`.
- **figs/** — Created when generating figures; EPS output from `plotting` (e.g. `LUM_CONTOUR.eps`). **Requires matplotlib** for the `plotting` module.

---

## Structure and use of the refactor

The Python code mirrors the two IDL pipelines:

1. **Luminosity-function pipeline** — `lfunc.create_Lobs()` builds the observed L distribution from CSV (MC trials, histogram); `lfunc.create_modelgrid()` builds model probability arrays over a grid of (L_lo, L_hi, γ); `lfunc.compare_obs_mod()` computes log-likelihood / χ², best fit, and writes `LFUNC_results.tex` and `LUM_CONTOUR.npz`. Entry point: `lfunc.mainproc()`.
2. **Progenitor mass MC pipeline** — `prog_mc.run_prog_mc()` reads CSV and M–L relation, runs trials (randomize dmod/BC/A/mag or L, L→mass via interpolated M–L, fit M_lo/M_hi with `m_hi_lo_fit`), builds 2D histogram and confidence levels, returns best M_min/M_max and arrays for plotting.

Helper modules: **histutils** (opthist, outxhist, wmean), **powerlaw** (randomp), **m_hi_lo_fit** (IMF fit with Γ = −1.35), **io_utils** (CSV and M–L read, data paths). **distribution** provides `generate_distribution`, `generate_sample`, `calculate_chi2`, and `run_mass_simulation` for mass-distribution modeling.

**Quick run:**

```python
from progenitors.davies18 import lfunc, prog_mc
BESTFITPARS, LPOBS, LPMOD = lfunc.mainproc(NTRIALS_OBS=1000, NTRIALS_MOD=1000)  # small for speed
result = prog_mc.run_prog_mc(obsfile="IIPprog_obsdata_2019.csv", MLfile="M-L_STARS.txt", ntrials=1000, seed=42)
# plotting.plot_all_from_compare()  # requires matplotlib and legacy/LUM_CONTOUR.npz
```

---

## Tests and validation

Unit tests live under the repository **tests/** in **tests/davies18/**.

- **test_histutils.py** — `opthist`, `outxhist`, `wmean`; bin count matches between opthist and outxhist.
- **test_powerlaw.py** — `randomp` shape, range, reproducibility with fixed seed.
- **test_m_hi_lo_fit.py** — `m_hi_lo_fit` bounds and monotonicity; `fit_m_hi_lo` round-trip and return types.
- **test_io_utils.py** — `get_data_path`, `read_obs_csv`, `obs_df_to_struct`, `read_ml_relation`.
- **test_lfunc.py** — `create_Lobs` shapes and normalization; `create_modelgrid` small grid; `compare_obs_mod` best-fit structure; `mainproc` end-to-end with reduced grid/trials.
- **test_prog_mc.py** — `run_prog_mc` return dict and shapes; reproducibility with fixed seed.
- **test_distribution.py** — `generate_distribution`, `generate_sample`, `calculate_chi2`, `run_mass_simulation` smoke.
- **test_plotting.py** — Plot functions run without error when matplotlib is available; `plot_all_from_compare` with optional npz path.

Run all davies18 tests:

```bash
pytest tests/davies18/ -v
```

**Validation that results match the original:** The refactor is written for algorithm parity with the IDL (same histogram binning, same power-law and IMF formulae, same interpolation). With the **same input data** and a **fixed random seed**, the Python code produces the same numerical behavior as the IDL. Reference outputs in **legacy/** (e.g. `LPOBS.sav`, `LUM_CONTOUR.sav`, `LFUNC_results.tex`) can be compared to Python outputs (e.g. `LPOBS.npz`, `LUM_CONTOUR.npz`, `LFUNC_results.tex`) when both are run with the same CSV, M–L file, and seed. The test suite uses fixed seeds (e.g. `seed=42`) so that runs are reproducible and regression in best-fit parameters or array shapes is caught.
