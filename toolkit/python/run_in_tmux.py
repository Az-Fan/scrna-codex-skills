#!/usr/bin/env python3
"""Launch a confirmed skill command in a detached, auditable tmux session."""
from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path


SESSION_RE = re.compile(r"^[A-Za-z0-9_.-]{1,64}$")


def timestamp() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def atomic_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--session", required=True, help="Unique tmux session name")
    parser.add_argument("--cwd", type=Path, required=True, help="Working directory for the command")
    parser.add_argument("--log", type=Path, required=True, help="Supervisor log (prefer tmux.log)")
    parser.add_argument("--status", type=Path, required=True, help="Machine-readable supervisor status JSON")
    parser.add_argument("--execute", action="store_true", help="Create the detached tmux session")
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("command", nargs=argparse.REMAINDER, help="Command following --")
    args = parser.parse_args()
    if args.command and args.command[0] == "--":
        args.command = args.command[1:]
    if not args.command:
        parser.error("a command is required after --")
    if not SESSION_RE.fullmatch(args.session):
        parser.error("--session must match [A-Za-z0-9_.-]{1,64}")
    args.cwd = args.cwd.expanduser().resolve()
    args.log = args.log.expanduser().resolve()
    args.status = args.status.expanduser().resolve()
    return args


def worker(args: argparse.Namespace) -> int:
    started = timestamp()
    state = {
        "schema_version": 1,
        "session": args.session,
        "status": "running",
        "started_at": started,
        "finished_at": None,
        "exit_code": None,
        "cwd": str(args.cwd),
        "command": args.command,
        "log": str(args.log),
    }
    atomic_json(args.status, state)
    args.log.parent.mkdir(parents=True, exist_ok=True)
    with args.log.open("a", encoding="utf-8") as handle:
        handle.write(f"[{started}] START session={args.session} command={shlex.join(args.command)}\n")
        handle.flush()
        try:
            result = subprocess.run(
                args.command,
                cwd=args.cwd,
                stdin=subprocess.DEVNULL,
                stdout=handle,
                stderr=subprocess.STDOUT,
                text=True,
            )
            exit_code = result.returncode
        except Exception as exc:  # Preserve a status record even when exec fails.
            handle.write(f"launcher error: {type(exc).__name__}: {exc}\n")
            exit_code = 125
        finished = timestamp()
        handle.write(f"[{finished}] END session={args.session} exit_code={exit_code}\n")
    state.update(
        status="completed" if exit_code == 0 else "failed",
        finished_at=finished,
        exit_code=exit_code,
    )
    atomic_json(args.status, state)
    return exit_code


def main() -> int:
    args = parse_args()
    if args.worker:
        return worker(args)
    if not args.cwd.is_dir():
        raise SystemExit(f"Working directory does not exist: {args.cwd}")
    plan = {
        "session": args.session,
        "cwd": str(args.cwd),
        "command": args.command,
        "log": str(args.log),
        "status_file": str(args.status),
        "execute": args.execute,
    }
    if not args.execute:
        print(json.dumps(plan, indent=2))
        return 0
    tmux = shutil.which("tmux")
    if tmux is None:
        raise SystemExit("tmux is not available on PATH")
    existing = subprocess.run(
        [tmux, "has-session", "-t", args.session],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    if existing.returncode == 0:
        raise SystemExit(f"Refusing duplicate active tmux session: {args.session}")
    args.log.parent.mkdir(parents=True, exist_ok=True)
    args.status.parent.mkdir(parents=True, exist_ok=True)
    worker_argv = [
        sys.executable,
        str(Path(__file__).resolve()),
        "--worker",
        "--session",
        args.session,
        "--cwd",
        str(args.cwd),
        "--log",
        str(args.log),
        "--status",
        str(args.status),
        "--",
        *args.command,
    ]
    launched = subprocess.run(
        [tmux, "new-session", "-d", "-s", args.session, shlex.join(worker_argv)],
        text=True,
        capture_output=True,
    )
    if launched.returncode:
        detail = launched.stderr.strip() or launched.stdout.strip()
        raise SystemExit(f"tmux launch failed: {detail}")
    plan["execute"] = True
    plan["status"] = "launched"
    print(json.dumps(plan, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
