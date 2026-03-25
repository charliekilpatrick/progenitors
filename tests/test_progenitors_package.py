"""Tests for progenitors package public API and pipeline entry point."""
import pytest


def test_package_import():
    import progenitors
    assert hasattr(progenitors, "params")
    assert hasattr(progenitors, "parse_arguments")
    assert hasattr(progenitors, "message")
    assert hasattr(progenitors, "__version__")


def test_version_format():
    """Version string is non-empty (importlib, git describe, or fallback)."""
    from progenitors import __version__
    assert isinstance(__version__, str)
    assert len(__version__) >= 1


def test_params_from_init():
    from progenitors import params
    assert "SHEET" in params or "metadata" in params
    assert "target" in params


def test_pipeline_import():
    from progenitors.pipeline import main, do_trim_metadata
    assert callable(main)
    assert callable(do_trim_metadata)


def test_main_entry_point():
    from progenitors.__main__ import main
    assert callable(main)


