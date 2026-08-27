#!/usr/bin/env python3
"""Install the nine released skills into a Codex-compatible skill directory."""
import argparse, shutil
from pathlib import Path

SKILLS=("scrna-calculate-qc-metrics","scrna-review-qc","scrna-standardize-input","scrna-analyze-subset","scrna-annotate-cells","scrna-find-cluster-markers","scrna-benchmark-integration","scrna-score-programs","scrna-run-differential-analysis")

def main():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument("--target",type=Path,default=Path.home()/".codex"/"skills"); p.add_argument("--force",action="store_true"); a=p.parse_args(); repo=Path(__file__).resolve().parents[1]; target=a.target.expanduser().resolve(); target.mkdir(parents=True,exist_ok=True)
    for name in SKILLS:
        source,destination=repo/"skills"/name,target/name
        if destination.exists():
            if not a.force: raise SystemExit(f"Already exists: {destination}; rerun with --force to replace")
            shutil.rmtree(destination)
        shutil.copytree(source,destination); print(f"installed {name} -> {destination}")
if __name__=="__main__": main()
