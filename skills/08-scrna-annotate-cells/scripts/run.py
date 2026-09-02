#!/usr/bin/env python3
from pathlib import Path
import sys
for parent in Path(__file__).resolve().parents:
    candidate = parent / "toolkit" / "python"
    if (candidate / "runner_bootstrap.py").is_file():
        sys.path.insert(0, str(candidate)); break
else:
    sys.path.insert(0, str(Path(__file__).resolve().parent))
from runner_bootstrap import run
run("08-scrna-annotate-cells", __file__)
