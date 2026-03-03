"""
Google Sheets and table column configuration.

Used by: progenitors.sheets.sheetproc
Location: progenitors/settings/sheets.py
"""
# Column names for the progenitor Google Sheet (order matters for upload/download)
SHEET_COLUMNS = [
    "Name",
    "YSEPZ",
    "TNS",
    "RA",
    "Dec",
    "Classification",
    "Host",
    "NED",
    "Discovery Date",
    "HST (pre-explosion)",
    "HST (post-explosion)",
    "JWST (pre-explosion)",
    "JWST (post-explosion)",
    "Spectrum",
    "Distance",
    "Distance Method",
    "Ref. (Distance)",
    "Redshift",
    "Ref. (Discovery)",
    "Ref. (Classification)",
    "Post-Explosion",
]

# OAuth scope for Google Sheets API
SCOPES = ["https://www.googleapis.com/auth/spreadsheets"]
