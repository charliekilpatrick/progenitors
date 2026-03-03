"""
Pipeline default settings.

Used by: progenitors.pipeline, progenitors.__main__, progenitors.options
Location: progenitors/settings/pipeline.py
"""
# Default metadata types to run when --metadata is not given
DEFAULT_METADATA = (
    "osc",
    "jwst",
    "ned",
    "distance",
    "tns",
    "yse",
    "ads",
    "hst",
)

# Features required at startup when no --metadata override (Sheets + YSE always)
DEFAULT_FEATURES = ["sheets", "yse"]

# Metadata keys to keep when trimming (do_trim_metadata); others are removed if not in table
TRIM_METADATA_KEEP_KEYS = ("all_sndata", "curr_table", "mask", "jwst", "hst")

# Default YSE SQL query id when not overridden by --yse-sql-query
DEFAULT_YSE_SQL_QUERY = "160"
