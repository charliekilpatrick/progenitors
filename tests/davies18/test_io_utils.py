"""Unit tests for progenitors.davies18.io_utils."""
import os
import tempfile
import pytest

import pandas as pd
import numpy as np

from progenitors.davies18 import io_utils


def test_get_data_path_resolves_csv():
    p = io_utils.get_data_path("IIPprog_obsdata_2022.csv")
    assert "davies18" in p
    assert p.endswith(".csv") or "data" in p


def test_read_obs_csv_returns_dataframe():
    df = io_utils.read_obs_csv("IIPprog_obsdata_2022.csv")
    assert isinstance(df, pd.DataFrame)
    assert len(df) >= 1
    assert "L_us" in df.columns or "L_s09" in df.columns or "snname" in df.columns or "SN name" in df.columns


def test_obs_df_to_struct():
    df = io_utils.read_obs_csv("IIPprog_obsdata_2022.csv")
    rows = io_utils.obs_df_to_struct(df)
    assert len(rows) == len(df)
    for r in rows:
        assert "snname" in r
        assert "dmod" in r
        assert "L_us" in r


def test_obs_df_to_struct_iuse():
    df = io_utils.read_obs_csv("IIPprog_obsdata_2022.csv")
    iuse = [0, 1, 2]
    rows = io_utils.obs_df_to_struct(df, iuse=iuse)
    assert len(rows) == 3


def test_read_ml_relation():
    mass, logl = io_utils.read_ml_relation("M-L_STARS.txt")
    assert len(mass) == len(logl)
    assert len(mass) >= 1
    assert np.all(np.isfinite(mass))
    assert np.all(np.isfinite(logl))


def test_read_obs_csv_absolute_path():
    df = io_utils.read_obs_csv("IIPprog_obsdata_2022.csv")
    path = io_utils.get_data_path("IIPprog_obsdata_2022.csv")
    if os.path.isfile(path):
        df2 = io_utils.read_obs_csv(path)
        assert len(df2) == len(df)
