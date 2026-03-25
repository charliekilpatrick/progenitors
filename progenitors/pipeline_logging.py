"""
Central logging setup for the CLI pipeline and related modules.

Uses the stdlib ``logging`` package (do not name any module ``logging.py``).
Default: compact stderr lines; on a TTY only the level word (INFO / WARNING /
ERROR / CRITICAL) is colored. Set ``NO_COLOR`` to disable ANSI. Set
``PROGENITORS_LOG_LEVEL`` or pass ``--log-level`` on the ``progenitors`` CLI.
"""
from __future__ import annotations

import logging
import os
import sys
from typing import Optional

# Compact wall time only (saves columns vs full ISO date).
_DATE_FMT = "%H:%M:%S"

# Color only these level names (DEBUG left plain).
_COLORED_LEVELS = frozenset(("INFO", "WARNING", "ERROR", "CRITICAL"))
_LEVEL_COLORS = {
    logging.INFO: "\033[32m",
    logging.WARNING: "\033[33m",
    logging.ERROR: "\033[31m",
    logging.CRITICAL: "\033[31m",
}
_RESET = "\033[0m"


class _CompactColorFormatter(logging.Formatter):
    """Minimal spacing; ANSI color on the level token only."""

    def __init__(self, use_color: bool) -> None:
        super().__init__(datefmt=_DATE_FMT)
        self._use_color = use_color

    def format(self, record: logging.LogRecord) -> str:
        record.message = record.getMessage()
        t = self.formatTime(record, self.datefmt)
        name = record.name
        if name.startswith("progenitors."):
            name = name[len("progenitors.") :]
        elif name == "progenitors":
            name = "pkg"

        lvl = record.levelname
        if self._use_color and lvl in _COLORED_LEVELS:
            color = _LEVEL_COLORS.get(record.levelno, "")
            if color:
                lvl_s = f"{color}{lvl}{_RESET}"
            else:
                lvl_s = lvl
        else:
            lvl_s = lvl

        # One space between fields; no space before message after colon.
        line = f"{t} {lvl_s} {name}:{record.message}"
        if record.exc_info:
            if not record.exc_text:
                record.exc_text = self.formatException(record.exc_info)
        if record.exc_text:
            if line[-1:] != "\n":
                line += "\n"
            line += record.exc_text
        if record.stack_info:
            if line[-1:] != "\n":
                line += "\n"
            line += self.formatStack(record.stack_info)
        return line


def setup_pipeline_logging(level: Optional[str] = None) -> None:
    """
    Configure the ``progenitors`` logger tree (idempotent).

    Parameters
    ----------
    level : str, optional
        DEBUG, INFO, WARNING, or ERROR. If None, uses env ``PROGENITORS_LOG_LEVEL``
        or INFO.
    """
    name = (level or os.environ.get("PROGENITORS_LOG_LEVEL") or "INFO").upper()
    numeric = getattr(logging, name, logging.INFO)

    root = logging.getLogger("progenitors")
    if not root.handlers:
        h = logging.StreamHandler(sys.stderr)
        use_color = (
            hasattr(sys.stderr, "isatty")
            and sys.stderr.isatty()
            and not os.environ.get("NO_COLOR", "").strip()
            and os.environ.get("TERM", "") != "dumb"
        )
        h.setFormatter(_CompactColorFormatter(use_color=use_color))
        root.addHandler(h)
    root.setLevel(numeric)
    for h in root.handlers:
        h.setLevel(numeric)

    for noisy in ("urllib3", "googleapiclient", "google_auth_httplib2"):
        logging.getLogger(noisy).setLevel(
            logging.DEBUG if numeric <= logging.DEBUG else logging.WARNING
        )
