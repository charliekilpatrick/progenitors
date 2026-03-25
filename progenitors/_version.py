"""Package version: resolved at import time (installed dist, git, or fallback)."""

from __future__ import annotations

import os
import re
import subprocess


def _normalize_git_describe_to_pep440(raw: str) -> str:
    """
    Convert ``git describe`` output to a PEP 440 string setuptools/pip accept.

    Examples
    --------
    ``v1.1.0-42-g001adbd-dirty`` -> ``1.1.0.dev42+g001adbd.dirty``
    ``v1.1.0`` -> ``1.1.0``
    ``abc1234`` (no tags) -> ``0.1.0+abc1234``
    """
    s = raw.strip()
    if not s:
        return "0.1.0+unknown"
    s = re.sub(r"^[vV]", "", s)

    m = re.match(
        r"^(\d+\.\d+(?:\.\d+)?)-(\d+)-g([0-9a-f]+)(-dirty)?$",
        s,
    )
    if m:
        base, n, sha, dirty = m.group(1), m.group(2), m.group(3), m.group(4)
        loc = "g" + sha + (".dirty" if dirty else "")
        return f"{base}.dev{n}+{loc}"

    m2 = re.match(r"^(\d+\.\d+(?:\.\d+)?)-dirty$", s)
    if m2:
        return f"{m2.group(1)}+dirty"

    if re.fullmatch(r"\d+\.\d+(?:\.\d+)?", s):
        return s

    if re.fullmatch(r"[0-9a-f]{7,40}", s):
        return f"0.1.0+g{s}"

    if re.fullmatch(r"[0-9a-f]{7,40}-dirty", s):
        sha = s[:-6]
        return f"0.1.0+g{sha}.dirty"

    safe = re.sub(r"[^a-zA-Z0-9.]+", "_", s)[:48]
    return f"0.1.0+local.{safe}"


def _compute_version() -> str:
    try:
        from importlib.metadata import PackageNotFoundError, version

        try:
            return version("progenitors")
        except PackageNotFoundError:
            pass
    except ImportError:
        pass

    pkg_dir = os.path.dirname(os.path.abspath(__file__))
    repo_root = os.path.dirname(pkg_dir)
    try:
        proc = subprocess.run(
            ["git", "describe", "--tags", "--dirty", "--always"],
            cwd=repo_root,
            capture_output=True,
            text=True,
            timeout=8,
            check=False,
        )
        if proc.returncode == 0 and proc.stdout.strip():
            return _normalize_git_describe_to_pep440(proc.stdout.strip())
    except (OSError, subprocess.SubprocessError):
        pass

    return "0.1.0+unknown"


__version__ = _compute_version()
