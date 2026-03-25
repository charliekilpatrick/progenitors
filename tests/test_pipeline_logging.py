"""Tests for pipeline logging setup."""
import logging
import re

import pytest


@pytest.fixture(autouse=True)
def _reset_progenitors_loggers():
    """Avoid idempotent handler skip affecting repeated tests."""
    log = logging.getLogger("progenitors")
    log.handlers.clear()
    log.setLevel(logging.WARNING)
    yield
    log.handlers.clear()


def test_setup_pipeline_logging_timestamp_and_plain_without_tty(monkeypatch, capsys):
    monkeypatch.setattr("sys.stderr.isatty", lambda: False)
    from progenitors.pipeline_logging import setup_pipeline_logging

    setup_pipeline_logging(level="INFO")
    logging.getLogger("progenitors.testlog").info("hello")
    err = capsys.readouterr().err
    assert "hello" in err
    assert "INFO" in err
    assert re.search(r"\d{2}:\d{2}:\d{2} INFO testlog:hello", err)


def test_setup_pipeline_logging_respects_no_color(monkeypatch, capsys):
    monkeypatch.setenv("NO_COLOR", "1")
    monkeypatch.setattr("sys.stderr.isatty", lambda: True)
    from progenitors.pipeline_logging import setup_pipeline_logging

    setup_pipeline_logging(level="INFO")
    logging.getLogger("progenitors.testlog2").warning("warn")
    err = capsys.readouterr().err
    assert "warn" in err
    assert "\033[" not in err
