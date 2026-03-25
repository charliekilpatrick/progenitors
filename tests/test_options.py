"""Unit tests for options module."""
import pytest


def test_is_number():
    from progenitors.options import is_number
    assert is_number(1.0) is True
    assert is_number("2.5") is True
    assert is_number("x") is False
    assert is_number(None) is False


def test_parse_coord_deg():
    from progenitors.options import parse_coord
    c = parse_coord("10.5", "20.0")
    assert c is not None
    assert abs(c.ra.deg - 10.5) < 0.01


def test_parse_coord_sexagesimal():
    from progenitors.options import parse_coord
    c = parse_coord("00:00:00", "00:00:00")
    assert c is not None


def test_parse_two_floats():
    from progenitors.options import parse_two_floats
    a, b = parse_two_floats("1.0 2.0")
    assert a == 1.0
    assert b == 2.0


def test_parse_two_floats_error():
    import argparse
    from progenitors.options import parse_two_floats
    with pytest.raises((argparse.ArgumentError, ValueError, TypeError)):
        parse_two_floats("1.0")


def test_parse_arguments_defaults(monkeypatch):
    import sys
    monkeypatch.setattr(sys, "argv", ["progenitors"])
    from progenitors.options import parse_arguments
    args = parse_arguments()
    assert args.redo is False
    assert args.alert is False
    assert args.update_classification is False
    assert args.update_tns_class is False
    assert args.log_level is None
    assert args.profile is False


def test_parse_arguments_profile(monkeypatch):
    import sys
    monkeypatch.setattr(sys, "argv", ["progenitors", "--profile"])
    from progenitors.options import parse_arguments
    args = parse_arguments()
    assert args.profile is True


def test_message(capsys):
    from progenitors.pipeline_logging import setup_pipeline_logging
    from progenitors.options import message

    setup_pipeline_logging("INFO")
    message("hello")
    _, err = capsys.readouterr()
    assert "hello" in err
