#!/usr/bin/env python3
"""Check one released skill against its server-side pixi environment."""
import argparse, json, shutil, subprocess
from pathlib import Path

from scrna_runtime import resolved_rscript

PROFILES = {
    "02-scrna-calculate-qc-metrics": ("01-scrna-qc", ["Seurat", "SeuratObject", "Matrix", "jsonlite"], ["qs", "RANN", "S4Vectors"]),
    "03-scrna-review-qc": ("01-scrna-qc", ["Seurat", "SeuratObject", "Matrix", "jsonlite"], ["qs"]),
    "01-scrna-standardize-input": ("01-scrna-qc", ["Seurat", "SeuratObject", "Matrix", "jsonlite"], ["qs"]),
    "08-scrna-analyze-subset": ("02-annotation", ["Seurat", "SeuratObject", "Matrix", "jsonlite"], ["qs"]),
    "07-scrna-annotate-cells": ("02-annotation", ["Seurat", "SeuratObject", "Matrix", "jsonlite"], ["qs"]),
    "06-scrna-find-cluster-markers": ("02-annotation", ["Seurat", "SeuratObject", "Matrix", "jsonlite"], ["qs"]),
    "04-scrna-preprocess-and-cluster": ("03-integration", ["Seurat", "SeuratObject", "Matrix", "jsonlite"], ["qs", "harmony"]),
    "05-scrna-benchmark-integration": ("03-integration", ["Seurat", "SeuratObject", "Matrix", "jsonlite"], ["qs", "harmony"]),
    "09-scrna-score-programs": ("05-pathway_program", ["Seurat", "Matrix", "jsonlite"], ["qs", "VISION", "AUCell", "progeny"]),
    "10-scrna-run-differential-analysis": ("06-deg-analysis", ["Seurat", "Matrix", "jsonlite", "DESeq2"], ["qs", "clusterProfiler", "msigdbr"]),
}

def main():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument("--skill",choices=sorted(PROFILES)); p.add_argument("--pixi-root",default="~/projects/scrna_envs"); p.add_argument("--pixi-executable"); a=p.parse_args()
    skill=a.skill or Path(__file__).resolve().parents[1].name
    if skill not in PROFILES: raise SystemExit(f"Unknown released skill: {skill}")
    env_dir,required,optional=PROFILES[skill]; packages=required+optional
    pixi=a.pixi_executable or shutil.which("pixi") or str(Path.home()/".pixi/bin/pixi")
    manifest=Path(a.pixi_root).expanduser().resolve()/env_dir/"pixi.toml"
    report={"skill":skill,"pixi_manifest":str(manifest),"pixi":pixi,"required_packages":{},"optional_packages":{},"compatible":False}
    if not manifest.is_file(): report["error"]="pixi manifest is missing"; print(json.dumps(report,indent=2)); return 2
    expr="p<-c("+",".join(json.dumps(x) for x in packages)+");v<-vapply(p,function(x)if(requireNamespace(x,quietly=TRUE))as.character(packageVersion(x))else NA_character_,character(1));cat(jsonlite::toJSON(as.list(v),auto_unbox=TRUE,na='null'))"
    env_r=resolved_rscript(skill, {"runtime": {"pixi_root": str(Path(a.pixi_root).expanduser().resolve())}})
    if env_r is None:
        report["error"]="registered pixi Rscript is missing; system R fallback is disabled"
        print(json.dumps(report,indent=2)); return 2
    command=[str(env_r),"-e",expr]
    report["runtime"]=str(env_r)
    done=subprocess.run(command,text=True,capture_output=True)
    if done.returncode: report["error"]=done.stderr.strip() or "dependency probe failed"
    else:
        versions=json.loads(done.stdout); report["required_packages"]={x:versions[x] for x in required}; report["optional_packages"]={x:versions[x] for x in optional}
        missing=[x for x in required if versions[x] is None]; report["compatible"]=not missing
        if missing: report["error"]="Missing required R packages: "+", ".join(missing)
        unavailable=[x for x in optional if versions[x] is None]
        if unavailable: report["note"]="Optional branches unavailable: "+", ".join(unavailable)
    print(json.dumps(report,indent=2)); return 0 if report["compatible"] else 2

if __name__=="__main__": raise SystemExit(main())
