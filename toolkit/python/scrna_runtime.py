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
    "01-scrna-standardize-input": {
        "required": ["project.id", "input.path", "input.format", "metadata.sample", "output_dir"],
        "artifacts": ["standardized_object", "samples_table", "cell_metadata", "field_mapping", "provenance"],
    },
    "07-scrna-annotate-cells": {
        "required": ["project.id", "input.object", "metadata.sample", "output_dir"],
        "artifacts": ["cluster_markers", "annotation_review", "cluster_plot", "run_manifest"],
    },
    "06-scrna-find-cluster-markers": {
        "required": ["project.id", "input.object", "output_dir"],
        "artifacts": ["cluster_markers", "top_cluster_markers", "cluster_marker_summary", "marker_dotplot", "run_manifest"],
    },
    "05-scrna-benchmark-integration": {
        "required": ["project.id", "input.object", "metadata.sample", "metadata.batch_variables", "benchmark.methods", "metrics", "plots", "output_dir"],
        "artifacts": ["method_runs", "metric_results", "method_summary", "design_confounding", "selected_plots", "benchmark_object", "recommendation", "run_manifest"],
    },
    "08-scrna-analyze-subset": {
        "required": ["project.id", "input.object", "metadata.sample", "metadata.cell_type", "subset.include", "output_dir"],
        "artifacts": ["subset_counts", "subset_metadata", "subset_summary", "provenance"],
    },
    "10-scrna-run-differential-analysis": {
        "required": ["project.id", "output_dir"],
        "artifacts": ["design_audit", "task_status", "complete_results", "differential_plots", "optional_enrichment", "run_manifest"],
    },
    "04-scrna-preprocess-and-cluster": {
        "required": ["project.id", "input.object", "output_dir"],
        "artifacts": ["preprocessed_clustered_object", "scenario_summary", "cell_assignments", "cluster_sizes", "resolution_stability", "workflow_state", "run_manifest"],
    },
    "09-scrna-score-programs": {
        "required": ["project.id", "input.object", "tasks", "output_dir"],
        "artifacts": ["scored_object", "score_matrices", "signature_coverage", "assay_feature_mapping", "score_summaries", "diagnostic_plots", "run_manifest"],
    },
}

DRIVERS = {
    "01-scrna-standardize-input": "standardize_input.R",
    "07-scrna-annotate-cells": "annotate_cells.R",
    "06-scrna-find-cluster-markers": "find_cluster_markers.R",
    "05-scrna-benchmark-integration": "integration_benchmark.R",
    "08-scrna-analyze-subset": "analyze_subset.R",
    "10-scrna-run-differential-analysis": "differential_analysis.R",
    "09-scrna-score-programs": "score_programs.R",
    "04-scrna-preprocess-and-cluster": "preprocess_cluster.R",
}

ENV_PROFILES = {
    "02-scrna-calculate-qc-metrics": "01-scrna-qc",
    "03-scrna-review-qc": "01-scrna-qc",
    "01-scrna-standardize-input": "01-scrna-qc",
    "07-scrna-annotate-cells": "02-annotation",
    "06-scrna-find-cluster-markers": "02-annotation",
    "05-scrna-benchmark-integration": "03-integration",
    "08-scrna-analyze-subset": "02-annotation",
    "10-scrna-run-differential-analysis": "06-deg-analysis",
    "04-scrna-preprocess-and-cluster": "03-integration",
    "09-scrna-score-programs": "05-pathway_program",
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


def resolved_rscript(skill, config):
    pixi_root = nested_get(config, "runtime.pixi_root") or os.environ.get("SCRNA_PIXI_ROOT")
    root = Path(os.path.expandvars(os.path.expanduser(str(pixi_root or "~/projects/scrna_envs"))))
    profile = ENV_PROFILES.get(skill)
    if not profile:
        return None
    candidate = root / profile / ".pixi/envs/default/bin/Rscript"
    return candidate.resolve() if candidate.is_file() else None


def default_argv(skill, config_path, config):
    rscript = resolved_rscript(skill, config)
    here = Path(__file__).resolve()
    drivers = [here.parents[1] / "R" / DRIVERS[skill], here.parent / DRIVERS[skill]]
    driver = next((path for path in drivers if path.is_file()), None)
    if rscript and driver is not None:
        return [str(rscript), str(driver), str(config_path.resolve())]
    return None


def expected_artifacts(skill, config):
    if skill != "04-scrna-preprocess-and-cluster":
        return SPECS[skill]["artifacts"]
    action = nested_get(config, "workflow.action") or "run"
    if action == "finalize_resolution":
        return [
            "preprocessed_clustered_object", "scenario_summary", "cell_assignments",
            "cluster_sizes", "sample_cluster_counts", "umap_diagnostics",
            "workflow_state", "session_info", "run_log", "run_manifest_finalize",
        ]
    scenarios = config.get("scenarios") or []
    awaiting_review = any(
        nested_get(scenario, "clustering.mode") == "scan"
        and (nested_get(scenario, "clustering.selection") or "review") == "review"
        for scenario in scenarios if isinstance(scenario, dict)
    )
    artifacts = [
        "preprocessed_clustered_object", "scenario_summary",
        "scenario_cluster_similarity", "elbow", "workflow_state", "session_info", "run_log", "run_manifest_preprocess",
    ]
    if any(nested_get(scenario, "clustering.mode") == "scan" for scenario in scenarios if isinstance(scenario, dict)):
        artifacts.extend(["resolution_stability", "resolution_umap_grid", "optional_clustree"])
    if not awaiting_review:
        artifacts.extend(["cell_assignments", "cluster_sizes", "sample_cluster_counts", "umap_diagnostics"])
    return artifacts


def validate(skill, config, config_path):
    errors, warnings = [], []
    spec = SPECS[skill]
    for field in spec["required"]:
        if is_blank(nested_get(config, field)):
            errors.append(f"missing required field: {field}")
    stage = nested_get(config, "analysis.stage") or "differential"
    source = nested_get(config, "input.object") or nested_get(config, "input.differential_table") or nested_get(config, "enrichment.input_results") or nested_get(config, "input.path")
    if skill == "10-scrna-run-differential-analysis":
        if stage == "enrichment_only":
            tables = nested_get(config, "input.differential_tables")
            if not source and not tables:
                errors.append("enrichment_only requires input.differential_table or input.differential_tables")
            if tables is not None and (not isinstance(tables, list) or not tables):
                errors.append("input.differential_tables must be a non-empty array")
            for index, item in enumerate(tables or []):
                if not isinstance(item, dict) or is_blank(item.get("path")):
                    errors.append(f"differential table {index + 1} requires path")
                elif not Path(os.path.expandvars(os.path.expanduser(str(item["path"])))).exists():
                    errors.append(f"differential table does not exist in this execution context: {item['path']}")
        else:
            for field in ("input.object", "metadata.sample", "metadata.condition"):
                if is_blank(nested_get(config, field)):
                    errors.append(f"missing required field: {field}")
    if skill == "04-scrna-preprocess-and-cluster":
        action = nested_get(config, "workflow.action") or "run"
        qc_status = str(nested_get(config, "input.qc_status") or "").lower()
        if qc_status not in {"filtered", "unfiltered"}:
            errors.append("input.qc_status must be filtered or unfiltered")
        elif qc_status == "unfiltered":
            if nested_get(config, "input.allow_unfiltered") is not True:
                errors.append("unfiltered input requires input.allow_unfiltered=true")
            else:
                warnings.append("unfiltered input explicitly authorized; outputs are exploratory_unfiltered")
        if action == "finalize_resolution":
            for field in ("finalize.scenario", "finalize.resolution"):
                if is_blank(nested_get(config, field)):
                    errors.append(f"missing required field: {field}")
            source_value = nested_get(config, "input.object")
            output_value = nested_get(config, "output_dir")
            if source_value and output_value and nested_get(config, "finalize.allow_separate_output") is not True:
                source_parent = Path(os.path.expandvars(os.path.expanduser(str(source_value)))).resolve().parent
                output_path = Path(os.path.expandvars(os.path.expanduser(str(output_value)))).resolve()
                if source_parent != output_path:
                    errors.append("finalize output_dir must equal the scan object directory unless finalize.allow_separate_output=true")
        else:
            for field in ("input.assay", "metadata.sample", "scenarios"):
                if is_blank(nested_get(config, field)):
                    errors.append(f"missing required field: {field}")
            scenarios = config.get("scenarios")
            if not isinstance(scenarios, list) or not scenarios:
                errors.append("scenarios must be a non-empty array")
            else:
                for index, scenario in enumerate(scenarios):
                    if not isinstance(scenario, dict) or is_blank(scenario.get("name")):
                        errors.append(f"scenario {index + 1} requires a name")
    if source:
        path = Path(os.path.expandvars(os.path.expanduser(str(source))))
        if not path.exists():
            errors.append(f"input does not exist in this execution context: {path}")
    if skill == "10-scrna-run-differential-analysis" and stage != "enrichment_only":
        comparisons = config.get("comparisons")
        if comparisons is None:
            comparison = config.get("comparison")
            comparisons = [comparison] if isinstance(comparison, dict) else []
        if not isinstance(comparisons, list) or not comparisons:
            errors.append("provide a non-empty comparisons array or legacy comparison object")
        else:
            for index, comparison in enumerate(comparisons):
                if not isinstance(comparison, dict) or is_blank(comparison.get("numerator")) or is_blank(comparison.get("denominator")):
                    errors.append(f"comparison {index + 1} requires numerator and denominator")
                elif comparison["numerator"] == comparison["denominator"]:
                    errors.append(f"comparison {index + 1} numerator and denominator must differ")
    if skill == "09-scrna-score-programs":
        tasks = config.get("tasks")
        if not isinstance(tasks, list) or not tasks:
            errors.append("tasks must be a non-empty JSON array")
        else:
            names = []
            supported = {"vision", "aucell", "ucell", "addmodulescore", "progeny"}
            for index, task in enumerate(tasks):
                if not isinstance(task, dict) or is_blank(task.get("name")) or is_blank(task.get("method")):
                    errors.append(f"task {index + 1} requires name and method")
                    continue
                names.append(str(task["name"]))
                if str(task["method"]).lower() not in supported:
                    errors.append(f"task {index + 1} has unsupported method: {task['method']}")
                if str(task["method"]).lower() != "progeny" and not isinstance(task.get("gene_sets"), dict):
                    errors.append(f"task {index + 1} requires a gene_sets object")
            if len(names) != len(set(names)):
                errors.append("task names must be unique")
    if skill == "05-scrna-benchmark-integration":
        supported_methods = {"none", "harmony", "rpca", "scvi", "scanvi", "bbknn", "precomputed"}
        supported_batch_metrics = {"ilisi", "batch_asw", "pcr_comparison", "graph_connectivity", "kbet"}
        supported_bio_metrics = {"clisi", "label_asw", "isolated_labels", "nmi", "ari"}
        supported_plots = {
            "score_barplot", "score_heatmap", "ranking_plot", "metric_tradeoff",
            "umap_by_batch", "umap_by_sample", "umap_by_condition", "umap_by_label",
            "marker_dotplot", "program_retention",
        }
        batch_variables = nested_get(config, "metadata.batch_variables")
        if not isinstance(batch_variables, list) or not batch_variables or any(is_blank(x) for x in batch_variables):
            errors.append("metadata.batch_variables must be a non-empty array of column names")
        elif len(batch_variables) != len(set(batch_variables)):
            errors.append("metadata.batch_variables must not contain duplicates")
        methods = nested_get(config, "benchmark.methods")
        if not isinstance(methods, list) or not methods:
            errors.append("benchmark.methods must be a non-empty array")
        else:
            scenario_names = []
            for index, method in enumerate(methods):
                if not isinstance(method, dict) or is_blank(method.get("name")):
                    errors.append(f"benchmark method {index + 1} requires name")
                    continue
                name = str(method["name"]).lower()
                scenario_names.append(name)
                if name not in supported_methods:
                    errors.append(f"benchmark method {index + 1} has unsupported name: {name}")
                grid = method.get("parameter_grid", {})
                if not isinstance(grid, dict):
                    errors.append(f"benchmark method {index + 1} parameter_grid must be an object")
                elif any(not isinstance(values, list) or not values for values in grid.values()):
                    errors.append(f"benchmark method {index + 1} parameter_grid values must be non-empty arrays")
                if name == "precomputed" and is_blank(method.get("reduction")):
                    errors.append(f"benchmark method {index + 1} precomputed method requires reduction")
                numeric_positive = {
                    "harmony": {"theta"}, "rpca": {"k_anchor", "k_weight", "nfeatures"},
                    "scvi": {"n_latent", "n_layers", "max_epochs"},
                    "scanvi": {"n_latent", "n_layers", "max_epochs", "scanvi_max_epochs"},
                    "bbknn": {"neighbors_within_batch", "n_pcs", "n_top_genes"},
                }.get(name, set())
                for parameter in numeric_positive & set(grid):
                    if any(not isinstance(value, (int, float)) or value <= 0 for value in grid[parameter]):
                        errors.append(f"benchmark method {index + 1} parameter {parameter} must contain positive numbers")
            if "none" not in scenario_names:
                warnings.append("uncorrected method 'none' will be injected as the required baseline")
        metric_cfg = config.get("metrics", {})
        if not isinstance(metric_cfg, dict):
            errors.append("metrics must be an object")
        else:
            batch_metrics = metric_cfg.get("batch_removal", [])
            bio_metrics = metric_cfg.get("biological_conservation", [])
            if not isinstance(batch_metrics, list) or not isinstance(bio_metrics, list):
                errors.append("metric groups must be arrays")
            else:
                unknown = (set(batch_metrics) - supported_batch_metrics) | (set(bio_metrics) - supported_bio_metrics)
                if unknown:
                    errors.append("unsupported benchmark metrics: " + ", ".join(sorted(unknown)))
                labels = nested_get(config, "metadata.biological_labels") or []
                if bio_metrics and (not isinstance(labels, list) or not labels):
                    errors.append("biological-conservation metrics require metadata.biological_labels")
        plots = config.get("plots")
        if not isinstance(plots, list):
            errors.append("plots must be an array")
        else:
            unknown_plots = set(plots) - supported_plots
            if unknown_plots:
                errors.append("unsupported benchmark plots: " + ", ".join(sorted(unknown_plots)))
            if set(plots) & {"marker_dotplot", "program_retention"} and not config.get("gene_programs"):
                errors.append("marker_dotplot and program_retention require non-empty gene_programs")
        scoring = config.get("scoring", {})
        if scoring.get("enabled"):
            batch_weight = scoring.get("batch_weight")
            biology_weight = scoring.get("biology_weight")
            if not isinstance(batch_weight, (int, float)) or not isinstance(biology_weight, (int, float)):
                errors.append("enabled scoring requires numeric batch_weight and biology_weight")
            elif batch_weight < 0 or biology_weight < 0 or abs(batch_weight + biology_weight - 1) > 1e-9:
                errors.append("scoring weights must be non-negative and sum to 1")
        prefix = nested_get(config, "benchmark.python_argv_prefix")
        if prefix is not None and (not isinstance(prefix, list) or not prefix or any(is_blank(x) for x in prefix)):
            errors.append("benchmark.python_argv_prefix must be a non-empty argv array")
    executor = config.get("executor", {})
    if executor and not isinstance(executor.get("argv", []), list):
        errors.append("executor.argv must be a JSON array, never a shell command string")
    if not executor and default_argv(skill, config_path, config) is None:
        warnings.append("registered pixi R executor is unavailable; dry-run remains available and system R will not be used")
    return errors, warnings


def make_manifest(skill, config, config_path, errors, warnings):
    source = nested_get(config, "input.object") or nested_get(config, "input.differential_table") or nested_get(config, "enrichment.input_results") or nested_get(config, "input.path")
    source_path = Path(os.path.expandvars(os.path.expanduser(str(source)))) if source else None
    input_record = {"path": str(source_path) if source_path else None}
    if source_path and source_path.is_file():
        input_record.update({"bytes": source_path.stat().st_size, "sha256": sha256(source_path)})
    argv = config.get("executor", {}).get("argv") or default_argv(skill, config_path, config)
    return {
        "schema_version": 1,
        "skill": skill,
        "created_at": dt.datetime.now(dt.timezone.utc).isoformat(),
        "config_path": str(config_path.resolve()),
        "config_sha256": sha256(config_path),
        "project_id": nested_get(config, "project.id"),
        "input": input_record,
        "output_dir": nested_get(config, "output_dir"),
        "expected_artifacts": expected_artifacts(skill, config),
        "resolved_argv": [str(value) for value in argv] if argv else None,
        "resolved_rscript": str(resolved_rscript(skill, config)) if not config.get("executor") and resolved_rscript(skill, config) else None,
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
    action_suffix = ""
    if skill == "04-scrna-preprocess-and-cluster":
        action = nested_get(config, "workflow.action") or "run"
        action_suffix = ".finalize" if action == "finalize_resolution" else ".scan"
    manifest_path = args.manifest or args.config.with_name(f"{skill}{action_suffix}.manifest.json")
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
    argv = config.get("executor", {}).get("argv") or default_argv(skill, args.config, config)
    if not argv:
        print("ERROR: --execute requires executor.argv", file=sys.stderr)
        return 2
    executable = shutil.which(str(argv[0]))
    if executable is None:
        print(f"ERROR: executable not found: {argv[0]}", file=sys.stderr)
        return 2
    output_dir = Path(os.path.expandvars(os.path.expanduser(str(config["output_dir"]))))
    output_dir.mkdir(parents=True, exist_ok=True)
    log_path = output_dir / "run.log"
    started = dt.datetime.now(dt.timezone.utc).isoformat()
    command = [executable] + [str(x) for x in argv[1:]]
    with log_path.open("a", encoding="utf-8") as log:
        log.write(f"[{started}] START {' '.join(command)}\n")
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1)
        assert process.stdout is not None
        for line in process.stdout:
            sys.stdout.write(line)
            sys.stdout.flush()
            log.write(line)
            log.flush()
        returncode = process.wait()
        finished = dt.datetime.now(dt.timezone.utc).isoformat()
        log.write(f"[{finished}] EXIT {returncode}\n")
    return returncode
