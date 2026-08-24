#!/usr/bin/env python3
"""Dependency-free runtime contract for the reusable scRNA-seq skills."""

import argparse
import datetime as dt
import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path


SPECS = {
    "standardize-scrna-input": {
        "required": ["project.id", "input.path", "input.format", "metadata.sample", "output_dir"],
        "artifacts": ["standardized_object", "samples_table", "cell_metadata", "field_mapping", "provenance"],
    },
    "review-scrna-qc": {
        "required": ["project.id", "input.object", "metadata.sample", "output_dir"],
        "artifacts": ["qc_summary", "qc_plots", "decision_table", "run_manifest"],
    },
    "annotate-scrna-cells": {
        "required": ["project.id", "input.object", "metadata.sample", "output_dir"],
        "artifacts": ["cluster_markers", "annotation_review", "cluster_plot", "run_manifest"],
    },
    "benchmark-scrna-integration": {
        "required": ["project.id", "input.object", "metadata.sample", "metadata.batch", "output_dir"],
        "artifacts": ["benchmark_table", "benchmark_object", "diagnostic_plots", "run_manifest"],
    },
    "analyze-scrna-subset": {
        "required": ["project.id", "input.object", "metadata.sample", "metadata.cell_type", "subset.include", "output_dir"],
        "artifacts": ["subset_counts", "subset_metadata", "subset_summary", "provenance"],
    },
    "run-scrna-differential-analysis": {
        "required": ["project.id", "input.object", "metadata.sample", "metadata.condition", "comparison.numerator", "comparison.denominator", "output_dir"],
        "artifacts": ["design_audit", "sample_cell_counts", "run_manifest"],
    },
}

DRIVERS = {
    "standardize-scrna-input": "standardize_input.R",
    "review-scrna-qc": "qc_report.R",
    "annotate-scrna-cells": "annotate_cells.R",
    "benchmark-scrna-integration": "integration_benchmark.R",
    "analyze-scrna-subset": "analyze_subset.R",
    "run-scrna-differential-analysis": "differential_analysis.R",
}


def nested_get(data, dotted):
    value = data
    for key in dotted.split("."):
        if not isinstance(value, dict) or key not in value:
            return None
        value = value[key]
    return value


def is_blank(value):
    return value is None or value == "" or value == []


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def default_argv(skill, config_path):
    rscript = shutil.which("Rscript")
    here = Path(__file__).resolve()
    drivers = [here.parents[1] / "R" / DRIVERS[skill], here.parent / DRIVERS[skill]]
    driver = next((path for path in drivers if path.is_file()), None)
    if rscript and driver is not None:
        return [rscript, str(driver), str(config_path.resolve())]
    return None


def validate(skill, config, config_path):
    errors, warnings = [], []
    spec = SPECS[skill]
    for field in spec["required"]:
        if is_blank(nested_get(config, field)):
            errors.append(f"missing required field: {field}")
    source = nested_get(config, "input.object") or nested_get(config, "input.path")
    if source:
        path = Path(os.path.expandvars(os.path.expanduser(str(source))))
        if not path.exists():
            errors.append(f"input does not exist in this execution context: {path}")
    numerator = nested_get(config, "comparison.numerator")
    denominator = nested_get(config, "comparison.denominator")
    if numerator is not None and numerator == denominator:
        errors.append("comparison numerator and denominator must differ")
    executor = config.get("executor", {})
    if executor and not isinstance(executor.get("argv", []), list):
        errors.append("executor.argv must be a JSON array, never a shell command string")
    if not executor and default_argv(skill, config_path) is None:
        warnings.append("default R executor is unavailable; dry-run remains available")
    return errors, warnings


def make_manifest(skill, config, config_path, errors, warnings):
    source = nested_get(config, "input.object") or nested_get(config, "input.path")
    source_path = Path(os.path.expandvars(os.path.expanduser(str(source)))) if source else None
    input_record = {"path": str(source_path) if source_path else None}
    if source_path and source_path.is_file():
        input_record.update({"bytes": source_path.stat().st_size, "sha256": sha256(source_path)})
    return {
        "schema_version": 1,
        "skill": skill,
        "created_at": dt.datetime.now(dt.timezone.utc).isoformat(),
        "config_path": str(config_path.resolve()),
        "project_id": nested_get(config, "project.id"),
        "input": input_record,
        "output_dir": nested_get(config, "output_dir"),
        "expected_artifacts": SPECS[skill]["artifacts"],
        "errors": errors,
        "warnings": warnings,
        "status": "blocked" if errors else "ready",
    }


def main(skill):
    if skill not in SPECS:
        raise SystemExit(f"unknown skill: {skill}")
    parser = argparse.ArgumentParser(description=f"Validate, plan, or execute {skill}")
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--execute", action="store_true")
    args = parser.parse_args()
    try:
        config = json.loads(args.config.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise SystemExit(f"invalid config: {exc}")
    errors, warnings = validate(skill, config, args.config)
    manifest = make_manifest(skill, config, args.config, errors, warnings)
    manifest_path = args.manifest or args.config.with_name(f"{skill}.manifest.json")
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    for message in warnings:
        print(f"WARNING: {message}", file=sys.stderr)
    if errors:
        for message in errors:
            print(f"ERROR: {message}", file=sys.stderr)
        return 2
    print(f"READY: {skill}; manifest={manifest_path}")
    if not args.execute:
        return 0
    argv = config.get("executor", {}).get("argv") or default_argv(skill, args.config)
    if not argv:
        print("ERROR: --execute requires executor.argv", file=sys.stderr)
        return 2
    executable = shutil.which(str(argv[0]))
    if executable is None:
        print(f"ERROR: executable not found: {argv[0]}", file=sys.stderr)
        return 2
    completed = subprocess.run([executable] + [str(x) for x in argv[1:]], check=False)
    return completed.returncode
