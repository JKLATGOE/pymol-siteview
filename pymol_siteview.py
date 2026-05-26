#!/usr/bin/env python
"""Compatibility wrapper for running from a source checkout with PyMOL."""

from __future__ import annotations

import sys
from pathlib import Path


try:
    from pymol_siteview.cli import main
except ModuleNotFoundError:
    src_dir = Path(__file__).resolve().parent / "src"
    if src_dir.exists():
        sys.path.insert(0, str(src_dir))
    from pymol_siteview.cli import main


if __name__ == "__main__":
    raise SystemExit(main())
