#!/usr/bin/env python3
"""Exercise tmux supervisor planning, execution, status, and duplicate refusal."""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import tempfile
import time
import uuid
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--runner",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "toolkit/python/run_in_tmux.py",
    )
    args = parser.parse_args()
    if shutil.which("tmux") is None:
        raise SystemExit("tmux is required for this release test")
    session = f"scrna-test-{uuid.uuid4().hex[:12]}"
    with tempfile.TemporaryDirectory(prefix="scrna-tmux-test-") as temporary:
        root = Path(temporary)
        log = root / "tmux.log"
        status = root / "tmux_status.json"
        base = [
            sys.executable,
            str(args.runner.resolve()),
            "--session",
            session,
            "--cwd",
            str(root),
            "--log",
            str(log),
            "--status",
            str(status),
            "--",
            sys.executable,
            "-c",
            "import time; print('tmux-test-ok', flush=True); time.sleep(2)",
        ]
        try:
            plan = subprocess.run(base, text=True, capture_output=True, check=True)
            if json.loads(plan.stdout).get("execute") is not False:
                raise RuntimeError("supervisor dry run did not report execute=false")
            subprocess.run(base[: base.index("--")] + ["--execute"] + base[base.index("--") :], check=True)
            duplicate = subprocess.run(
                base[: base.index("--")] + ["--execute"] + base[base.index("--") :],
                text=True,
                capture_output=True,
            )
            if duplicate.returncode == 0 or "duplicate active tmux session" not in duplicate.stderr:
                raise RuntimeError(f"duplicate session was not refused: {duplicate}")
            deadline = time.monotonic() + 10
            state = None
            while time.monotonic() < deadline:
                if status.is_file():
                    state = json.loads(status.read_text(encoding="utf-8"))
                    if state.get("status") in {"completed", "failed"}:
                        break
                time.sleep(0.1)
            if state is None or state.get("status") != "completed" or state.get("exit_code") != 0:
                raise RuntimeError(f"unexpected terminal status: {state}")
            if "tmux-test-ok" not in log.read_text(encoding="utf-8"):
                raise RuntimeError("command output missing from tmux log")
        finally:
            subprocess.run(["tmux", "kill-session", "-t", session], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    print("PASS tmux supervisor")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
