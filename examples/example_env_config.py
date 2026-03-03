"""Example: test progenitors.env_config (no validation of real credentials). Run from repo root."""
from progenitors.env_config import get_env, ENV_BY_FEATURE, validate_env

# get_env
assert get_env("NONEXISTENT_XYZ") == ""
assert get_env("NONEXISTENT_XYZ", "default") == "default"
# ENV_BY_FEATURE
assert "sheets" in ENV_BY_FEATURE and "PROGENITORS_SHEET" in ENV_BY_FEATURE["sheets"]
assert "yse" in ENV_BY_FEATURE and "tns" in ENV_BY_FEATURE
# validate_env with no features required does not exit
validate_env([])
print("progenitors.env_config OK")
