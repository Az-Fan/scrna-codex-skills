#!/usr/bin/env python3
"""Build self-contained .skill archives without committing generated runtimes."""
import argparse, shutil, tempfile
from pathlib import Path
from install_skills import SKILLS, build_skill

def main():
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument("--output",type=Path,default=Path("dist")); args=parser.parse_args(); output=args.output.resolve(); output.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="scrna-package-") as temp:
        for name in SKILLS:
            built=Path(temp)/name; build_skill(name,built)
            archive=shutil.make_archive(str(output/name),"zip",root_dir=built.parent,base_dir=name)
            final=output/(name+".skill"); Path(archive).replace(final); print(final)
if __name__=="__main__": main()
