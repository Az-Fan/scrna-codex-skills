#!/usr/bin/env python3
"""Build and install released skills from canonical toolkit sources."""
import argparse
import json
import shutil
import tempfile
import datetime
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
MANIFEST = REPO / "release" / "runtime-manifest.json"

def load_manifest():
    return json.loads(MANIFEST.read_text(encoding="utf-8"))["skills"]

SKILLS = tuple(load_manifest())
LEGACY_SKILLS = {"03.1-scrna-apply-qc-filter": "04-scrna-apply-qc-filter"}

def build_skill(name, destination):
    spec = load_manifest()[name]
    source = REPO / "skills" / name
    shutil.copytree(source, destination, ignore=shutil.ignore_patterns("__pycache__", "*.pyc", "*.pyo"))
    scripts = destination / "scripts"
    scripts.mkdir(exist_ok=True)
    for filename in spec.get("python", []):
        shutil.copy2(REPO / "toolkit" / "python" / filename, scripts / filename)
    for filename in spec.get("r", []):
        shutil.copy2(REPO / "toolkit" / "R" / filename, scripts / filename)
    references = destination / "references"
    references.mkdir(exist_ok=True)
    for filename in spec.get("references", []):
        shutil.copy2(REPO / "toolkit" / "references" / filename, references / filename)

def install_skills(target, force=False):
    target = Path(target).expanduser().resolve()
    if target == REPO or target == REPO / "skills" or REPO / "skills" in target.parents:
        raise SystemExit("Refusing to install generated runtimes into canonical skills source")
    target.mkdir(parents=True, exist_ok=True)
    existing = [name for name in SKILLS if (target / name).exists()]
    if existing and not force:
        raise SystemExit(f"Already exists: {target / existing[0]}; rerun with --force to replace")
    legacy = [name for name in LEGACY_SKILLS if (target / name).exists()]
    backup = None
    moved = []
    published = []
    # 先构建全部 skill，再迁出旧目录；失败时恢复原安装。
    with tempfile.TemporaryDirectory(prefix=".scrna-install-", dir=target.parent) as temporary:
        staged = Path(temporary)
        for name in SKILLS:
            build_skill(name, staged / name)
        if existing or legacy:
            backup_root = target.parent / ("." + target.name + "-backups")
            backup_root.mkdir(exist_ok=True)
            backup = Path(tempfile.mkdtemp(prefix=datetime.datetime.now().strftime("%Y%m%dT%H%M%S-") , dir=backup_root))
        try:
            for name in existing + legacy:
                (target / name).rename(backup / name)
                moved.append(name)
            for name in SKILLS:
                (staged / name).rename(target / name)
                published.append(name)
        except BaseException:
            for name in published:
                (target / name).rename(staged / name)
            for name in reversed(moved):
                (backup / name).rename(target / name)
            raise
    for name in SKILLS:
        print(f"installed {name} -> {target / name}")
    if backup:
        print(f"Previous installation preserved outside discovery: {backup}")
    return backup


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target",type=Path,default=Path.home()/".codex"/"skills")
    parser.add_argument("--force",action="store_true")
    args=parser.parse_args()
    install_skills(args.target, args.force)
if __name__=="__main__": main()
