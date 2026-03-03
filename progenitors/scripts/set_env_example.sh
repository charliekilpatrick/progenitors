#!/usr/bin/env bash
# Example: set environment variables for the progenitors pipeline.
# All authentication (API keys, tokens, passwords) is read from the environment.
# The pipeline validates required variables at startup and exits with a clear
# message if any are missing.
#
# Usage (from repo root, with package installed):
#   source progenitors/scripts/set_env_example.sh   # then run: progenitors
# Or copy this file, fill in real values, and source your copy:
#   cp progenitors/scripts/set_env_example.sh my_env.sh
#   # edit my_env.sh with your credentials
#   source my_env.sh && progenitors

# ----- Google Sheets (required for pipeline) -----
export PROGENITORS_SHEET="your-google-sheet-id"
# Optional: test sheet
# export PROGENITORS_SHEET_TEST="your-test-sheet-id"
# Optional: path to credentials.json (default: repo_root/credentials.json)
# export PROGENITORS_CREDENTIALS_PATH="/path/to/credentials.json"

# ----- YSE-PZ (required for pipeline) -----
export YSE_USER="your-yse-username"
export YSE_PASSWORD="your-yse-password"

# ----- TNS (required if --metadata includes tns) -----
export TNS_API_KEY="your-tns-api-key"
export TNS_BOT_NAME="your-bot-name"
export TNS_BOT_ID="your-bot-id"

# ----- NASA ADS (required if --metadata includes ads) -----
export ADS_AUTHCODE="your-ads-api-token"

# ----- Email alerts (required if using --alert) -----
export GMAIL_LOGIN="your@gmail.com"
export GMAIL_PASSWORD="your-app-password"
export MY_EMAIL="recipient@example.com"

# After sourcing, run the pipeline from repo root:
#   progenitors
#   progenitors --metadata osc,jwst,ned,distance,tns,yse,ads,hst
#   progenitors --alert
