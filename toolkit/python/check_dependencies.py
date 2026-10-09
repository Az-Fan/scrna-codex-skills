#!/usr/bin/env python3
"""Check registered runtimes and the branches requested in a skill config."""
import argparse
import json
import subprocess
from pathlib import Path
from scrna_runtime import nested_get, environment_project, environment_record, integration_python_prefix, sccoda_python, cnmf_python

CORE = ["Seurat", "SeuratObject", "Matrix", "jsonlite"]
PROFILES = {
    "16-scrna-discover-programs": (CORE, ["qs"]),
    "01-scrna-standardize-input": (CORE, ["qs", "hdf5r"]),
    "02-scrna-calculate-qc-metrics": (CORE + ["ggplot2"], ["qs", "RANN", "S4Vectors", "celda", "SingleCellExperiment"]),
    "03-scrna-review-qc": (CORE + ["ggplot2"], ["qs"]),
    "04-scrna-apply-qc-filter": (CORE + ["digest"], ["qs"]),
    "05-scrna-benchmark-integration": (CORE + ["qs"], ["harmony"]),
    "06-scrna-preprocess-and-cluster": (CORE + ["ggplot2", "patchwork"], ["qs", "harmony", "clustree"]),
    "07-scrna-find-cluster-markers": (CORE + ["ggplot2", "patchwork"], ["qs"]),
    "08-scrna-annotate-cells": (CORE + ["qs", "ggplot2", "patchwork", "digest"], []),
    "09-scrna-export-subset": (CORE + ["qs"], []),
    "10-scrna-score-programs": (CORE, ["qs", "ggplot2", "VISION", "AUCell", "UCell", "progeny", "scMetabolism", "msigdbr"]),
    "11-scrna-run-differential-analysis": (CORE + ["DESeq2", "ggplot2"], ["qs", "clusterProfiler", "msigdbr", "apeglm"]),
    "12-scrna-run-pathway-enrichment": (["jsonlite", "clusterProfiler", "msigdbr", "ggplot2"], ["Seurat", "Matrix", "qs", "DESeq2"]),
    "13-scrna-test-cell-abundance": (CORE + ["ggplot2"], ["qs", "speckle", "miloR", "sccomp", "DCATS", "cmdstanr"]),
    "14-scrna-visualize-cell-composition": (["jsonlite", "ggplot2", "patchwork", "scales"], ["Seurat", "qs"]),
    "15-scrna-visualize-gene": (CORE + ["ggplot2", "patchwork"], ["qs", "ggrepel"]),
}


def package_requirements(skill, config):
    required, optional = (list(x) for x in PROFILES[skill])
    if skill == "16-scrna-discover-programs" and nested_get(config, "input.type") == "matrix":
        return [], []
    paths = [nested_get(config, key) for key in ("input.object", "input.path", "output.object_name")]
    if any(str(x or "").lower().endswith(".qs") for x in paths) or nested_get(config, "output.object_format") == "qs":
        required.append("qs")
    if skill in {"06-scrna-preprocess-and-cluster", "10-scrna-score-programs"} and (nested_get(config, "output.object_format") or "qs") == "qs":
        required.append("qs")
    if skill == "01-scrna-standardize-input" and (str(paths[1] or "").lower().endswith(".h5") or nested_get(config, "input.format") in {"h5", "10x_h5"}):
        required.append("hdf5r")
    if skill == "14-scrna-visualize-cell-composition" and nested_get(config, "input.object"):
        required.append("Seurat")
    if skill == "06-scrna-preprocess-and-cluster" and any(nested_get(x, "harmony.enabled") for x in config.get("scenarios", [])):
        required.append("harmony")
    if skill == "05-scrna-benchmark-integration" and any(str(x.get("name", "")).lower() == "harmony" for x in nested_get(config, "benchmark.methods") or []):
        required.append("harmony")
    if skill == "10-scrna-score-programs":
        packages = {"vision": "VISION", "aucell": "AUCell", "ucell": "UCell", "progeny": "progeny", "addmodulescore": "Seurat"}
        for task in config.get("tasks", []):
            package = packages.get(str(task.get("method", "")).lower())
            if package:
                required.append(package)
            source = nested_get(task, "gene_sets.source")
            if source in {"msigdb", "scmetabolism"}:
                required.append("msigdbr" if source == "msigdb" else "scMetabolism")
        if nested_get(config, "visualization.enabled") is not False:
            required.append("ggplot2")
    if skill == "13-scrna-test-cell-abundance":
        method_packages = {"propeller": "speckle", "milo": "miloR", "sccomp": "sccomp", "dcats": "DCATS"}
        methods = nested_get(config, "analysis.methods") or list(method_packages) + ["sccoda"]
        required.extend(method_packages[x] for x in methods if x in method_packages)
        if "sccomp" in methods:
            required.append("cmdstanr")
    return sorted(set(required)), sorted(set(optional) - set(required))


def probe_python(prefix, modules):
    expression = """import importlib.util, importlib.metadata, json
modules = MODULES
distributions = {"scib_metrics": "scib-metrics", "scvi": "scvi-tools", "sklearn": "scikit-learn"}
versions = {}
for name in modules:
    try:
        versions[name] = importlib.metadata.version(distributions.get(name, name))
    except importlib.metadata.PackageNotFoundError:
        versions[name] = None
print(json.dumps({"modules": {name: importlib.util.find_spec(name) is not None for name in modules}, "versions": versions}))
""".replace("MODULES", repr(modules))
    try:
        done = subprocess.run(list(prefix) + ["-c", expression], text=True, capture_output=True)
    except OSError as exc:
        return {}, str(exc), {}
    if done.returncode:
        return {}, done.stderr.strip() or "Python dependency probe failed", {}
    try:
        result = json.loads(done.stdout)
        present, versions = result["modules"], result["versions"]
    except ValueError:
        return {}, "Python dependency probe did not return JSON", {}
    missing = [name for name, available in present.items() if not available]
    error = ("Missing Python modules: " + ", ".join(missing)) if missing else None
    if versions.get("scib_metrics") == "0.5.7" and int((versions.get("pandas") or "0").split(".")[0]) >= 3:
        error = "scib-metrics 0.5.7 result tables require pandas < 3"
    return present, error, versions


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skill", choices=sorted(PROFILES))
    parser.add_argument("--config", type=Path, help="check the selected input/output formats and analysis branches")
    parser.add_argument("--pixi-root")
    parser.add_argument("--pixi-executable")
    args = parser.parse_args()
    skill = args.skill or Path(__file__).resolve().parents[1].name
    if skill not in PROFILES:
        raise SystemExit(f"Unknown released skill: {skill}")
    config = json.loads(args.config.read_text()) if args.config else {}
    if args.pixi_root:
        config.setdefault("runtime", {})["pixi_root"] = args.pixi_root
    project = environment_project(skill, config)
    required, optional = package_requirements(skill, config)
    report = {"skill": skill, "pixi_manifest": str(project / "pixi.toml"), "required_packages": {}, "optional_packages": {}, "compatible": False, "errors": []}
    report["environment"] = environment_record(skill, config)
    rscript = project / ".pixi/envs" / (nested_get(config, "pixi.environment") or "default") / "bin/Rscript"
    if skill == "16-scrna-discover-programs" and nested_get(config, "input.type") == "matrix":
        if not (project / "pixi.toml").is_file():
            report["errors"].append("Registered pixi manifest is missing")
    elif not (project / "pixi.toml").is_file() or not rscript.is_file():
        report["errors"].append("Registered pixi manifest or Rscript is missing; system R fallback is disabled")
    else:
        report["runtime"] = str(rscript)
        packages = required + optional
        expr = "p<-c(" + ",".join(json.dumps(x) for x in packages) + ");v<-vapply(p,function(x)if(requireNamespace(x,quietly=TRUE))as.character(packageVersion(x))else NA_character_,character(1));cat(jsonlite::toJSON(as.list(v),auto_unbox=TRUE,na='null'))"
        done = subprocess.run([str(rscript), "-e", expr], text=True, capture_output=True)
        if done.returncode:
            report["errors"].append(done.stderr.strip() or "R dependency probe failed")
        else:
            versions = json.loads(done.stdout)
            report["required_packages"] = {x: versions[x] for x in required}
            report["optional_packages"] = {x: versions[x] for x in optional}
            missing = [x for x in required if versions[x] is None]
            if missing:
                report["errors"].append("Missing required R packages: " + ", ".join(missing))
    if skill == "05-scrna-benchmark-integration":
        methods = [str(x.get("name", "")).lower() for x in nested_get(config, "benchmark.methods") or []]
        metrics = config.get("metrics", {})
        metric_requested = bool(metrics.get("batch_removal") or metrics.get("biological_conservation"))
        if not args.config or config.get("plots") or metric_requested or set(methods) & {"scvi", "scanvi", "bbknn"}:
            modules = ["numpy", "pandas", "scipy", "anndata", "scanpy"]
            if metric_requested or not args.config:
                modules.append("scib_metrics")
            if set(methods) & {"scvi", "scanvi"}:
                modules.append("scvi")
            if "bbknn" in methods:
                modules.append("bbknn")
            prefix = integration_python_prefix(config)
            report["python_argv_prefix"] = prefix
            present, error, versions = probe_python(prefix, modules)
            report["python_modules"] = present
            report["python_versions"] = versions
            if error:
                report["errors"].append(error)
    if skill == "13-scrna-test-cell-abundance":
        methods = nested_get(config, "analysis.methods") or ["sccomp", "sccoda"]
        if "sccomp" in methods:
            cmdstan = project / ".pixi/envs/default/bin/cmdstan"
            report["cmdstan_path"] = str(cmdstan)
            if not cmdstan.is_dir():
                report["errors"].append("Registered CmdStan installation is missing")
        if "sccoda" in methods:
            python = sccoda_python(config)
            report["sccoda_runtime"] = python
            present, error, versions = probe_python([python], ["pertpy", "anndata", "numpyro", "jax"])
            report["python_modules"] = present
            report["python_versions"] = versions
            if error:
                report["errors"].append(error)
    if skill == "16-scrna-discover-programs":
        python = cnmf_python(config)
        present, error, versions = probe_python([python], ["cnmf", "numpy", "pandas", "scipy", "anndata", "matplotlib", "sklearn"])
        report.update(cnmf_runtime=python, python_modules=present, python_versions=versions)
        if error:
            report["errors"].append(error)
    report["compatible"] = not report["errors"]
    if report["errors"]:
        report["error"] = "; ".join(report["errors"])
    unavailable = [x for x, version in report["optional_packages"].items() if version is None]
    if unavailable:
        report["note"] = "Optional branches unavailable: " + ", ".join(unavailable)
    if not args.config:
        report["scope"] = "baseline; use --config to check the selected formats and methods"
    print(json.dumps(report, indent=2))
    return 0 if report["compatible"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
