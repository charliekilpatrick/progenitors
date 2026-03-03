"""
I/O for davies18: progenitor CSV and M–L relation files.

Data files live in the package data/ directory. Paths resolve via
get_data_path(); CSV columns are normalized to match IDL readcol names.
"""
import os
import pandas as pd
import numpy as np

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(_THIS_DIR, "data")
LEGACY_DIR = os.path.join(_THIS_DIR, "legacy")


def get_data_path(filename):
    """
    Resolve path to a file in data/ or legacy/.

    Parameters
    ----------
    filename : str
        Base name (e.g. 'IIPprog_obsdata_2022.csv').

    Returns
    -------
    str
        Full path if file exists in data/ or legacy/; otherwise
        data_dir/filename.
    """
    for d in (DATA_DIR, LEGACY_DIR, _THIS_DIR):
        p = os.path.join(d, filename)
        if os.path.isfile(p):
            return p
    return os.path.join(DATA_DIR, filename)


def read_obs_csv(filename="IIPprog_obsdata_2022.csv"):
    """
    Read progenitor observation CSV into a DataFrame with normalized column names.

    Uses get_data_path if filename is not an existing path. Column names are
    mapped to IDL-style names (snname, dmod, L_us, dL_us, etc.).

    Parameters
    ----------
    filename : str, optional
        CSV filename or path. Default 'IIPprog_obsdata_2022.csv'.

    Returns
    -------
    pandas.DataFrame
        Numeric columns coerced; snname, band, notes left as object.
    """
    path = filename if os.path.isfile(filename) else get_data_path(filename)
    df = pd.read_csv(path)
    # Normalize column names (first row may have spaces)
    df.columns = [c.strip() for c in df.columns]
    # Map to expected names (IDL form)
    col_map = {
        "SN name": "snname",
        "D(Mpc)": "dist",
        "eD": "edist",
        "DM": "dmod",
        "eDM": "edmod",
        "mag": "mag",
        "emag": "emag",
        "band": "band",
        "E(B-v)": "ebmv",
        "dE(B-V)": "eebmv",
        "Av": "av",
        "eAI": "eai",
        "notes": "notes",
        "A_lam": "alam",
        "aA": "ealam",
        "BC_lam": "BClam",
        "dBC": "eBClam",
        "L_s09": "L_s09",
        "dL_s09": "dL_s09",
        "Ls15": "L_s15",
        "dLs15": "dL_s15",
        "L_s15": "L_s15",
        "dL_s15": "dL_s15",
        "L_us": "L_us",
        "dL_us": "dL_us",
    }
    for old, new in col_map.items():
        if old in df.columns and new not in df.columns:
            df = df.rename(columns={old: new})
    # Coerce numeric; leave snname, band, notes as object
    for c in df.columns:
        if c in ("snname", "band", "notes"):
            continue
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def obs_df_to_struct(df, iuse=None):
    """
    Convert DataFrame to list of dicts for create_Lobs / create_modelgrid.

    Parameters
    ----------
    df : pandas.DataFrame
        Output of read_obs_csv (or similar columns).
    iuse : array-like of int, optional
        Row indices to keep. If None, all rows are used.

    Returns
    -------
    list of dict
        One dict per row with keys snname, dmod, edmod, mag, emag, alam,
        ealam, BClam, eBClam, L_us, dL_us, notes.
    """
    if iuse is not None:
        df = df.iloc[iuse].reset_index(drop=True)
    rows = []
    for i in range(len(df)):
        r = df.iloc[i]
        rows.append({
            "snname": r.get("snname", ""),
            "dmod": float(r.get("dmod", np.nan)),
            "edmod": float(r.get("edmod", np.nan)),
            "mag": float(r.get("mag", np.nan)),
            "emag": float(r.get("emag", np.nan)),
            "alam": float(r.get("alam", np.nan)),
            "ealam": float(r.get("ealam", np.nan)),
            "BClam": float(r.get("BClam", np.nan)),
            "eBClam": float(r.get("eBClam", np.nan)),
            "L_us": float(r.get("L_us", np.nan)),
            "dL_us": float(r.get("dL_us", np.nan)),
            "notes": str(r.get("notes", "")),
        })
    return rows


def read_ml_relation(filename="M-L_STARS.txt"):
    """
    Read M–L relation file (whitespace-separated: mass, log L).

    Parameters
    ----------
    filename : str, optional
        Filename in data/ or path. Default 'M-L_STARS.txt'.

    Returns
    -------
    tuple of ndarray
        (mass_mod, lfin_mod) as 1D float arrays.

    Raises
    ------
    FileNotFoundError
        If file is not found.
    """
    path = get_data_path(filename)
    if not os.path.isfile(path):
        raise FileNotFoundError(f"M–L file not found: {path}")
    df = pd.read_csv(path, sep=r"\s+", header=None, names=["mass_mod", "lfin_mod"], comment="#")
    return df["mass_mod"].values.astype(float), df["lfin_mod"].values.astype(float)
