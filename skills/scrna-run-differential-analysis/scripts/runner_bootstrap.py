#!/usr/bin/env python3
"""Locate the shared runtime from source checkout or an installed skill bundle."""

import os
import sys
from pathlib import Path


def run(skill, script_file):
    script = Path(script_file).resolve()
    candidates = []
    configured = os.environ.get("SCRNA_SKILLS_HOME")
    if configured:
        candidates.append(Path(configured) / "toolkit" / "python")
    candidates.extend(parent / "toolkit" / "python" for parent in script.parents)
    candidates.append(script.parent)
    for candidate in candidates:
        if (candidate / "scrna_runtime.py").is_file():
            sys.path.insert(0, str(candidate))
            from scrna_runtime import main
            raise SystemExit(main(skill))
    raise SystemExit("scrna_runtime.py not found; set SCRNA_SKILLS_HOME or install the bundled runtime")
