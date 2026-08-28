#!/usr/bin/env python3
"""Build and smoke-test isolated copies of all released skills."""
import argparse, py_compile, subprocess, sys, tempfile
from pathlib import Path
from install_skills import SKILLS, build_skill

def check(command):
    result=subprocess.run(command,text=True,capture_output=True)
    if result.returncode: raise SystemExit("FAILED: "+" ".join(map(str,command))+"\n"+result.stdout+result.stderr)

def main():
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument("--quick-validate",type=Path); args=parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="scrna-skills-") as temp:
        root=Path(temp)
        for name in SKILLS:
            installed=root/name; build_skill(name,installed)
            if args.quick_validate: check([sys.executable,str(args.quick_validate),str(installed)])
            for path in installed.rglob("*.py"): py_compile.compile(str(path),doraise=True)
            check([sys.executable,str(installed/"scripts"/"run.py"),"--help"])
            check([sys.executable,str(installed/"scripts"/"check_dependencies.py"),"--help"])
            tmux_runner=installed/"scripts"/"run_in_tmux.py"
            if tmux_runner.is_file(): check([sys.executable,str(tmux_runner),"--help"])
            print(f"PASS {name}")
if __name__=="__main__": main()
