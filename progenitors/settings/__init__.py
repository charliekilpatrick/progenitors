"""
Pipeline and module settings.

Hard-coded configuration values are kept here so they are easy to find and change.
Each submodule documents which package modules use it.

| Settings module   | Used by |
|-------------------|---------|
| pipeline          | progenitors.pipeline, progenitors.__main__, progenitors.options |
| env_config        | progenitors.env_config |
| sheets            | progenitors.sheets.sheetproc |
| util              | progenitors.util |
| sed_analysis      | progenitors.sed.analysis (SED fit parameter names) |
"""
# sed_analysis.MODEL_FIT_PARAMS is used by progenitors.sed.analysis; import via from progenitors.settings import sed_analysis
from .pipeline import (
    DEFAULT_METADATA,
    DEFAULT_FEATURES,
    TRIM_METADATA_KEEP_KEYS,
    DEFAULT_YSE_SQL_QUERY,
)
from .env_config import ENV_BY_FEATURE, ENV_DESCRIPTIONS, ENV_EXAMPLE_SCRIPT
from .sheets import SHEET_COLUMNS, SCOPES
from .util import DM_HIERARCHY, EMAIL_FROM_ADDR, EMAIL_SMTP_SERVER, EMAIL_SMTP_PORT, EMAIL_SUBJECT, EMAIL_MSG_TEMPLATE
from . import sed_analysis

__all__ = [
    "DEFAULT_METADATA",
    "DEFAULT_FEATURES",
    "TRIM_METADATA_KEEP_KEYS",
    "DEFAULT_YSE_SQL_QUERY",
    "ENV_BY_FEATURE",
    "ENV_DESCRIPTIONS",
    "ENV_EXAMPLE_SCRIPT",
    "SHEET_COLUMNS",
    "SCOPES",
    "DM_HIERARCHY",
    "EMAIL_FROM_ADDR",
    "EMAIL_SMTP_SERVER",
    "EMAIL_SMTP_PORT",
    "EMAIL_SUBJECT",
    "EMAIL_MSG_TEMPLATE",
]
