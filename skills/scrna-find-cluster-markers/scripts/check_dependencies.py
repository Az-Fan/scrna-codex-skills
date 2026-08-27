#!/usr/bin/env python3
"""Check one released skill against its server-side pixi environment."""

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path


PROFILES = {
    "scrna-standardize-input": ("01-scrna-qc", ["Seurat", "SeuratObject", "Matrix", "jsonlite"], ["qs"]),
    "scrna-analyze-subset": ("02-annotation", ["Seurat", "SeuratObject", "Matrix", "jsonlite"], ["qs"]),
    "scrna-annotate-cells": ("02-annotation", ["Seurat", "SeuratObject", "Matrix", "jsonlite"], ["qs"]),
    "scrna-find-cluster-markers": ("02-annotation", ["Seurat", "SeuratObject", "Matrix", "jsonlite"], ["qs"]),
    "scrna-benchmark-integration": ("03-integration", ["Seurat", "SeuratObject", "Matrix", "jsonlite"], ["qs", "harmony"]),
    "scrna-run-differential-analysis": ("06-deg-analysis", ["Seurat", "Matrix", "jsonlite", "DESeq2"], ["qs", "clusterProfiler", "msigdbr"]),
    "scrna-score-programs": ("05-pathway_program", ["Seurat", "Matrix", "jsonlite"], ["qs", "VISION", "AUCell", "progeny"]),
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skill", choices=sorted(PROFILES))
    parser.add_argument("--pixi-root", default="~/projects/scrna_envs")
    parser.add_argument("--pixi-executable")
    args = parser.parse_args()
    skill = args.skill or Path(__file__).resolve().parents[1].name
    if skill not in PROFILES:
        raise SystemExit(f"Unknown released skill: {skill}")
    env_dir, required, optional = PROFILES[skill]
    packages = required + optional
    pixi = args.pixi_executable or shutil.which("pixi") or str(Path.home() / ".pixi/bin/pixi")
    manifest = Path(args.pixi_root).expanduser().resolve() / env_dir / "pixi.toml"
    report = {"skill": skill, "pixi_manifest": str(manifest), "pixi": pixi,
              "required_packages": {}, "optional_packages": {}, "compatible": False}
    if not Path(pixi).is_file() or not manifest.is_file():
        report["error"] = "pixi executable or manifest is missing"
        print(json.dumps(report, indent=2))
        return 2
    expression = (
        "p<-c(" + ",".join(json.dumps(x) for x in packages) + ");"
        "v<-vapply(p,function(x)if(requireNamespace(x,quietly=TRUE))"
        "as.character(packageVersion(x))else NA_character_,character(1));"
        "cat(jsonlite::toJSON(as.list(v),auto_unbox=TRUE,na='null'))"
    )
    environment_rscript = manifest.parent / ".pixi" / "envs" / "default" / "bin" / "Rscript"
    if environment_rscript.is_file():
        command = [str(environment_rscript), "-e", expression]
        report["runtime"] = str(environment_rscript)
    else:
        command = [pixi, "run", "--locked", "--manifest-path", str(manifest), "--", "Rscript", "-e", expression]
        report["runtime"] = "pixi run --locked"
    completed = subprocess.run(command, text=True, capture_output=True)
    if completed.returncode:
        report["error"] = completed.stderr.strip() or "pixi dependency probe failed"
    else:
        versions = json.loads(completed.stdout)
        report["required_packages"] = {name: versions[name] for name in required}
        report["optional_packages"] = {name: versions[name] for name in optional}
        missing = [name for name in required if versions[name] is None]
        report["compatible"] = not missing
        if missing:
            report["error"] = "Missing required R packages: " + ", ".join(missing)
        unavailable_optional = [name for name in optional if versions[name] is None]
        if unavailable_optional:
            report["note"] = "Optional branches unavailable: " + ", ".join(unavailable_optional)
    print(json.dumps(report, indent=2))
    return 0 if report["compatible"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
