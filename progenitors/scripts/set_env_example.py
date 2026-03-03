#!/usr/bin/env python3
"""
Example script that populates environment variables with example/placeholder
values so the progenitors pipeline can be run (e.g. for testing env validation).

All real authentication must come from environment variables. This script sets
example values so that:
  - The pipeline's startup validation sees every variable as "set".
  - You can replace the placeholders with real values in your shell or a .env.

Usage:
  # Set example vars in current shell (then run pipeline)
  eval $(python -m progenitors.scripts.set_env_example)

  # Or export to a file and source it
  python -m progenitors.scripts.set_env_example > my_env.sh
  source my_env.sh
  progenitors
"""
import os

# Variables required by feature (see progenitors.env_config.ENV_BY_FEATURE)
EXAMPLE_VARS = {
    "PROGENITORS_SHEET": "your-google-sheet-id",
    "PROGENITORS_SHEET_TEST": "your-test-sheet-id",
    "YSE_USER": "your-yse-username",
    "YSE_PASSWORD": "your-yse-password",
    "TNS_API_KEY": "your-tns-api-key",
    "TNS_BOT_NAME": "your-bot-name",
    "TNS_BOT_ID": "your-bot-id",
    "ADS_AUTHCODE": "your-ads-api-token",
    "GMAIL_LOGIN": "your@gmail.com",
    "GMAIL_PASSWORD": "your-app-password",
    "MY_EMAIL": "recipient@example.com",
}


def main():
    import sys
    for key, value in EXAMPLE_VARS.items():
        safe = value.replace("'", "'\"'\"'")
        print(f"export {key}='{safe}'")


if __name__ == "__main__":
    main()
