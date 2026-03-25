"""Unit tests for progenitors.env_config."""
import logging
import os
import pytest


def test_get_env():
    from progenitors.env_config import get_env
    assert get_env("NONEXISTENT_VAR_XYZ") == ""
    assert get_env("NONEXISTENT_VAR_XYZ", "default") == "default"
    os.environ["_TEST_PROGENITORS_VAR"] = "foo"
    try:
        assert get_env("_TEST_PROGENITORS_VAR") == "foo"
    finally:
        del os.environ["_TEST_PROGENITORS_VAR"]


def test_alert_reclassification_sheet_keys():
    from progenitors.settings.pipeline import ALERT_RECLASSIFICATION_SHEET_KEYS

    assert "Type Ia" in ALERT_RECLASSIFICATION_SHEET_KEYS
    assert "Type II-P/II-L" in ALERT_RECLASSIFICATION_SHEET_KEYS
    assert len(ALERT_RECLASSIFICATION_SHEET_KEYS) == 5


def test_env_by_feature():
    from progenitors.env_config import ENV_BY_FEATURE
    assert "sheets" in ENV_BY_FEATURE
    assert "PROGENITORS_SHEET" in ENV_BY_FEATURE["sheets"]
    assert "yse" in ENV_BY_FEATURE
    assert "tns" in ENV_BY_FEATURE
    assert "ads" in ENV_BY_FEATURE
    assert "email" in ENV_BY_FEATURE


def test_validate_env_missing_exits(caplog, monkeypatch):
    from progenitors.env_config import validate_env

    monkeypatch.delenv("PROGENITORS_SHEET", raising=False)
    with caplog.at_level(logging.ERROR, logger="progenitors.env_config"):
        with pytest.raises(SystemExit):
            validate_env(["sheets"])
    blob = caplog.text
    assert "PROGENITORS_SHEET" in blob
    assert "Missing" in blob or "required" in blob.lower()


def test_validate_env_sheets_ok_when_set(tmp_path, monkeypatch):
    from progenitors.env_config import validate_env
    cred = tmp_path / "credentials.json"
    cred.write_text("{}")
    monkeypatch.setenv("PROGENITORS_SHEET", "test-sheet-id")
    # Should not exit when sheet id is set and credentials path exists
    validate_env(["sheets"], credentials_path=str(cred))


def test_validate_env_empty_features():
    from progenitors.env_config import validate_env
    # No features -> nothing required -> should not exit
    validate_env([])
    validate_env(["nonexistent_feature"])
