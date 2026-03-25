"""
Environment variable configuration and validation.

All authentication (API keys, tokens, passwords) is read from environment
variables. Call validate_env() at pipeline start to fail fast with clear
messages if required variables are missing.

Settings: env var names and descriptions live in progenitors/settings/env_config.py
"""
import logging
import os
import sys

from .settings.env_config import (
    ENV_BY_FEATURE,
    ENV_DESCRIPTIONS,
    ENV_EXAMPLE_SCRIPT,
)

logger = logging.getLogger(__name__)


def get_env(key, default=""):
    """Return environment variable value or default. Use for building params."""
    return os.environ.get(key, default)


def validate_env(features, credentials_path=None):
    """
    Check that all environment variables required for the given features are set.

    Args:
        features: Iterable of feature names (e.g. 'sheets', 'yse', 'tns', 'ads', 'email').
        credentials_path: For feature 'sheets', path to credentials.json. If provided
            and file does not exist, validation fails.

    Exits with code 1 and prints missing variables if any are missing or empty.
    """
    required = []
    for f in features:
        required.extend(ENV_BY_FEATURE.get(f, []))
    required = list(dict.fromkeys(required))  # preserve order, no duplicates

    missing = [k for k in required if not (get_env(k) or "").strip()]
    if missing:
        lines = [
            "Missing required environment variables:",
            "",
        ]
        for k in missing:
            desc = ENV_DESCRIPTIONS.get(k, k)
            lines.append(f"  {k}")
            lines.append(f"    -> {desc}")
        lines.append("")
        lines.append("Set them in your shell or use the example script:")
        lines.append(f"  {ENV_EXAMPLE_SCRIPT}")
        lines.append("")
        logger.error("%s", "\n".join(lines))
        sys.exit(1)

    if "sheets" in features and credentials_path and not os.path.isfile(credentials_path):
        logger.error(
            "Google Sheets credentials file not found: %s. "
            "Download from Google Cloud Console and place at repo root, or set "
            "PROGENITORS_CREDENTIALS_PATH.",
            credentials_path,
        )
        sys.exit(1)
