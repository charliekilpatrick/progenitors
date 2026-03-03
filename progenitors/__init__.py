"""
Progenitors: data and analysis for progenitor systems of explosive transients.

Provides a Google Sheets–driven pipeline for progenitor targets, SED fitting,
and the Davies et al. (2018) progenitor analysis (luminosity function and
mass MC) in Python. Credentials and API keys are read from environment
variables; see progenitors.env_config.

Attributes
----------
__version__ : str
    Package version (from _version).
"""
from ._version import __version__
from .util import params
from .options import parse_arguments, message

__all__ = ["__version__", "params", "parse_arguments", "message"]
