#!/usr/bin/env python3
"""Build and install released skills from canonical toolkit sources."""
import argparse
import json
import shutil
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
MANIFEST = REPO / "release" / "runtime-manifest.json"

def load_manifest():
    return json.loads(MANIFEST.read_text(encoding="utf-8"))["skills"]

SKILLS = tuple(load_manifest())

def build_skill(name, destination):
    spec = load_manifest()[name]
    source = REPO / "skills" / name
    shutil.copytree(source, destination)
    scripts = destination / "scripts"
    scripts.mkdir(exist_ok=True)
    for filename in spec.get("python", []):
        shutil.copy2(REPO / "toolkit" / "python" / filename, scripts / filename)
    for filename in spec.get("r", []):
        shutil.copy2(REPO / "toolkit" / "R" / filename, scripts / filename)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target",type=Path,default=Path.home()/".codex"/"skills")
    parser.add_argument("--force",action="store_true")
    args=parser.parse_args(); target=args.target.expanduser().resolve(); target.mkdir(parents=True,exist_ok=True)
    for name in SKILLS:
        destination=target/name
        if destination.exists():
            if not args.force: raise SystemExit(f"Already exists: {destination}; rerun with --force to replace")
            shutil.rmtree(destination)
        build_skill(name,destination); print(f"installed {name} -> {destination}")
if __name__=="__main__": main()
