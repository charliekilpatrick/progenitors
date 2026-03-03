"""
Environment variable names and descriptions for validation.

Used by: progenitors.env_config
Location: progenitors/settings/env_config.py
"""
# Feature name -> list of required environment variable names
ENV_BY_FEATURE = {
    "sheets": [
        "PROGENITORS_SHEET",
    ],
    "yse": [
        "YSE_USER",
        "YSE_PASSWORD",
    ],
    "tns": [
        "TNS_API_KEY",
        "TNS_BOT_NAME",
        "TNS_BOT_ID",
    ],
    "ads": [
        "ADS_AUTHCODE",
    ],
    "email": [
        "GMAIL_LOGIN",
        "GMAIL_PASSWORD",
        "MY_EMAIL",
    ],
}

# Human-readable descriptions for missing-variable error messages
ENV_DESCRIPTIONS = {
    "PROGENITORS_SHEET": "Google Sheet ID for progenitor data",
    "PROGENITORS_SHEET_TEST": "Optional test Google Sheet ID",
    "YSE_USER": "YSE-PZ portal username",
    "YSE_PASSWORD": "YSE-PZ portal password",
    "TNS_API_KEY": "TNS API key",
    "TNS_BOT_NAME": "TNS bot name",
    "TNS_BOT_ID": "TNS bot ID",
    "ADS_AUTHCODE": "NASA ADS API token",
    "GMAIL_LOGIN": "Gmail address for sending alerts",
    "GMAIL_PASSWORD": "Gmail app password or token",
    "MY_EMAIL": "Recipient email for alerts",
}

# Path shown in error message when env vars are missing
ENV_EXAMPLE_SCRIPT = "progenitors/scripts/set_env_example.sh"
