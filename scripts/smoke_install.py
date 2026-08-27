#!/usr/bin/env python3
"""Smoke-test isolated copies of all seven released skills."""

import argparse
import py_compile
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from install_skills import SKILLS


def check(command):
    result = subprocess.run(command, text=True, capture_output=True)
    if result.returncode:
        raise SystemExit("FAILED: " + " ".join(map(str, command)) + "\n" + result.stdout + result.stderr)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quick-validate", type=Path)
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    with tempfile.TemporaryDirectory(prefix="scrna-skills-") as temp:
        root = Path(temp)
        for name in SKILLS:
            installed = root / name
            shutil.copytree(repo / "skills" / name, installed)
            if args.quick_validate:
                check([sys.executable, str(args.quick_validate), str(installed)])
            for path in installed.rglob("*.py"):
                py_compile.compile(str(path), doraise=True)
            check([sys.executable, str(installed / "scripts" / "run.py"), "--help"])
            check([sys.executable, str(installed / "scripts" / "check_dependencies.py"), "--help"])
            print(f"PASS {name}")


if __name__ == "__main__":
    main()
