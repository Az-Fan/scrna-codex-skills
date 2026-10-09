#!/usr/bin/env python3
"""Dependency-free runtime contract for the reusable scRNA-seq skills."""

import importlib.util
import argparse
import datetime as dt
import hashlib
import json
import os
import shutil
import subprocess
import sys
import uuid
import re
import itertools
from decimal import Decimal
from pathlib import Path


_delivery_spec = importlib.util.spec_from_file_location("scrna_result_delivery", Path(__file__).with_name("result_delivery.py"))
_delivery = importlib.util.module_from_spec(_delivery_spec)
_delivery_spec.loader.exec_module(_delivery)


SPECS = {
    "17-scrna-infer-grn": {
        "required": ["project.id", "metadata.sample", "output_dir", "inference.mode"],
        "artifacts": ["input_audit", "cell_membership", "inference_unit_audit", "feature_status", "resource_coverage", "adjacencies", "motif_enrichment", "regulons", "single_cell_activity", "sample_summary", "task_status", "workflow_state", "handoff", "run_manifest"],
    },
    "16-scrna-discover-programs": {
        "required": ["project.id", "metadata.sample", "output_dir", "cnmf.components"],
        "artifacts": ["input_audit", "k_selection_stats", "k_selection_plot", "reviewed_consensus_usage", "spectra", "top_genes", "sample_summaries", "program_review", "workflow_state", "run_manifest"],
    },
    "01-scrna-standardize-input": {
        "required": ["project.id", "input.path", "input.format", "metadata.sample", "output_dir"],
        "artifacts": ["standardized_object", "samples_table", "cell_metadata", "field_mapping", "provenance"],
    },
    "04-scrna-apply-qc-filter": {
        "required": ["project.id", "input.object", "input.decision_table", "metadata.sample", "decision.include_all_true", "decision.expected_retained_cells", "approval.status", "output_dir"],
        "artifacts": ["filtered_object", "cell_filter_decisions", "sample_retention", "optional_condition_retention", "decision_reason_counts", "approval_record", "session_info", "run_manifest"],
    },
    "08-scrna-annotate-cells": {
        "required": ["project.id", "input.object", "metadata.sample", "workflow.action", "output_dir"],
        "artifacts": ["cluster_markers_or_existing_marker_record", "annotation_review", "annotation_umaps", "optional_annotated_object", "annotation_tables", "run_manifest"],
    },
    "07-scrna-find-cluster-markers": {
        "required": ["project.id", "input.object", "output_dir"],
        "artifacts": ["cluster_markers", "top_cluster_markers", "cluster_marker_summary", "marker_dotplot", "run_manifest"],
    },
    "05-scrna-benchmark-integration": {
        "required": ["project.id", "input.object", "metadata.sample", "metadata.batch_variables", "benchmark.methods", "metrics", "output_dir"],
        "artifacts": ["method_runs", "metric_results", "method_summary", "design_confounding", "selected_plots", "benchmark_object", "recommendation", "recommendation_status", "run_manifest"],
    },
    "09-scrna-export-subset": {
        "required": ["project.id", "input.object", "metadata.sample", "metadata.cell_type", "subset.include", "output_dir"],
        "artifacts": ["subset_counts", "subset_metadata", "subset_summary", "provenance"],
    },
    "11-scrna-run-differential-analysis": {
        "required": ["project.id", "output_dir"],
        "artifacts": ["design_audit", "task_status", "complete_results", "differential_plots", "optional_enrichment", "run_manifest"],
    },
    "12-scrna-run-pathway-enrichment": {
        "required": ["project.id", "output_dir"],
        "artifacts": ["task_status", "complete_enrichment_results", "compact_enrichment_plots", "identifier_mapping", "run_manifest"],
    },
    "13-scrna-test-cell-abundance": {
        "required": ["project.id", "metadata.sample", "metadata.condition", "metadata.cell_type", "analysis.methods", "analysis.denominator.mode", "comparisons", "output_dir"],
        "artifacts": ["sample_cell_counts", "sample_cell_proportions", "design_audit", "task_status", "complete_method_results", "compact_diagnostic_plots", "method_concordance", "run_manifest"],
    },
    "14-scrna-visualize-cell-composition": {
        "required": ["project.id", "output_dir", "metadata.sample", "composition.variables", "composition.denominator.description"],
        "artifacts": ["composition_counts", "composition_proportions", "sample_coverage", "composition_audit", "plot_status", "figures", "run_manifest"],
    },
    "15-scrna-visualize-gene": {
        "required": ["project.id", "input.object", "genes", "metadata.sample", "expression.assay", "output_dir"],
        "artifacts": ["gene_status", "cell_expression_summary", "sample_expression", "target_gene_summary", "plot_status", "figures", "run_manifest"],
    },
    "06-scrna-preprocess-and-cluster": {
        "required": ["project.id", "input.object", "output_dir"],
        "artifacts": ["preprocessed_clustered_object", "scenario_summary", "cell_assignments", "cluster_sizes", "resolution_stability", "workflow_state", "run_manifest"],
    },
    "10-scrna-score-programs": {
        "required": ["project.id", "input.object", "tasks", "output_dir"],
        "artifacts": ["scored_object", "score_matrices", "signature_coverage", "assay_feature_mapping", "score_summaries", "diagnostic_plots", "run_manifest"],
    },
}

DRIVERS = {
    "17-scrna-infer-grn": "grn_pipeline.py",
    "16-scrna-discover-programs": "cnmf_discovery.py",
    "01-scrna-standardize-input": "standardize_input.R",
    "04-scrna-apply-qc-filter": "apply_qc_filter.R",
    "08-scrna-annotate-cells": "annotate_cells.R",
    "07-scrna-find-cluster-markers": "find_cluster_markers.R",
    "05-scrna-benchmark-integration": "integration_benchmark.R",
    "09-scrna-export-subset": "analyze_subset.R",
    "11-scrna-run-differential-analysis": "differential_analysis.R",
    "12-scrna-run-pathway-enrichment": "differential_analysis.R",
    "13-scrna-test-cell-abundance": "cell_abundance.R",
    "14-scrna-visualize-cell-composition": "visualize_cell_composition.R",
    "15-scrna-visualize-gene": "visualize_gene.R",
    "10-scrna-score-programs": "score_programs.R",
    "06-scrna-preprocess-and-cluster": "preprocess_cluster.R",
}

ENV_PROFILES = {
    "17-scrna-infer-grn": "04-grn",
    "16-scrna-discover-programs": "05-pathway_program",
    "02-scrna-calculate-qc-metrics": "01-scrna-qc",
    "03-scrna-review-qc": "01-scrna-qc",
    "01-scrna-standardize-input": "01-scrna-qc",
    "04-scrna-apply-qc-filter": "01-scrna-qc",
    "08-scrna-annotate-cells": "02-annotation",
    "07-scrna-find-cluster-markers": "02-annotation",
    "05-scrna-benchmark-integration": "03-integration",
    "09-scrna-export-subset": "02-annotation",
    "11-scrna-run-differential-analysis": "06-deg-analysis",
    "12-scrna-run-pathway-enrichment": "06-deg-analysis",
    "13-scrna-test-cell-abundance": "07-cell-abundance",
    "14-scrna-visualize-cell-composition": "02-annotation",
    "15-scrna-visualize-gene": "02-annotation",
    "06-scrna-preprocess-and-cluster": "03-integration",
    "10-scrna-score-programs": "05-pathway_program",
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


def environment_root(config):
    pixi_root = nested_get(config, "runtime.pixi_root") or os.environ.get("SCRNA_PIXI_ROOT")
    return Path(os.path.expandvars(os.path.expanduser(str(pixi_root or "~/projects/scrna_envs")))).resolve()


def environment_project(skill, config):
    configured = nested_get(config, "pixi.project") if skill in {"02-scrna-calculate-qc-metrics", "03-scrna-review-qc"} else None
    if configured:
        project = Path(os.path.expandvars(os.path.expanduser(str(configured)))).resolve()
        return project.parent if project.name == "pixi.toml" else project
    return environment_root(config) / ENV_PROFILES[skill]


def environment_record(skill, config):
    project = environment_project(skill, config)
    files = {}
    for name in ("pixi.toml", "pixi.lock", "supplemental-r.json"):
        path = project / name
        files[name] = {"path": str(path), "sha256": sha256(path) if path.is_file() else None}
    record = {"project": str(project), "files": files}
    if skill == "05-scrna-benchmark-integration":
        record["python_argv_prefix"] = integration_python_prefix(config)
    if skill == "13-scrna-test-cell-abundance":
        record["sccoda_python"] = sccoda_python(config)
    if skill == "16-scrna-discover-programs":
        record["cnmf_python"] = cnmf_python(config)
    if skill == "17-scrna-infer-grn":
        record["pyscenic_python"] = pyscenic_python(config)
    return record


def integration_python_prefix(config):
    prefix = list(nested_get(config, "benchmark.python_argv_prefix") or [str(environment_root(config) / "03-integration/.pixi/envs/scvi/bin/python")])
    prefix[0] = os.path.expandvars(os.path.expanduser(str(prefix[0])))
    return prefix


def sccoda_python(config):
    configured = nested_get(config, "runtime.sccoda_python") or str(environment_root(config) / "07-cell-abundance/.pixi/envs/sccoda/bin/python")
    return os.path.expandvars(os.path.expanduser(str(configured)))


def cnmf_python(config):
    configured = nested_get(config, "runtime.cnmf_python") or str(environment_root(config) / "05-pathway_program/.pixi/envs/cnmf/bin/python")
    return os.path.expandvars(os.path.expanduser(str(configured)))


def resolved_rscript(skill, config):
    profile = ENV_PROFILES.get(skill)
    if not profile:
        return None
    candidate = environment_project(skill, config) / ".pixi/envs/default/bin/Rscript"
    return candidate.resolve() if candidate.is_file() else None


def pyscenic_python(config):
    configured = nested_get(config, "runtime.pyscenic_python") or str(environment_root(config) / "04-grn/.pixi/envs/pyscenic/bin/python")
    return os.path.expandvars(os.path.expanduser(str(configured)))


def default_argv(skill, config_path, config):
    if skill == "17-scrna-infer-grn":
        driver = Path(__file__).with_name("grn_pipeline.py")
        python = pyscenic_python(config)
        return [python, str(driver), str(config_path.resolve())] if driver.is_file() and Path(python).is_file() else None
    if skill == "16-scrna-discover-programs":
        here = Path(__file__).resolve()
        driver = here.with_name("cnmf_discovery.py")
        python = cnmf_python(config)
        if driver.is_file() and Path(python).is_file():
            return [python, str(driver), str(config_path.resolve())]
        return None
    rscript = resolved_rscript(skill, config)
    here = Path(__file__).resolve()
    drivers = [here.parents[1] / "R" / DRIVERS[skill], here.parent / DRIVERS[skill]]
    driver = next((path for path in drivers if path.is_file()), None)
    if rscript and driver is not None:
        return [str(rscript), str(driver), str(config_path.resolve())]
    return None


def expected_artifacts(skill, config):
    if skill == "17-scrna-infer-grn" and (nested_get(config, "workflow.action") or "prepare") == "prepare":
        return ["input_audit", "cell_membership", "inference_unit_audit", "feature_status", "resource_coverage", "task_status", "workflow_state", "run_manifest"]
    if skill == "16-scrna-discover-programs":
        artifacts = ["input_audit", "feature_status", "cell_metadata", "k_selection_stats", "k_selection_plot", "workflow_state", "run_manifest"]
        if (nested_get(config, "workflow.action") or "discover") != "discover":
            artifacts.extend(["consensus_usage_raw_and_fractional", "gene_spectra", "top_genes", "program_gmt", "sample_summaries", "program_review", "task_status"])
        if nested_get(config, "interpretation.gmt"):
            artifacts.extend(["optional_program_enrichment", "gene_set_coverage", "enrichment_universe"])
        return artifacts
    if skill != "06-scrna-preprocess-and-cluster":
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
    if skill == "16-scrna-discover-programs":
        from cnmf_contract import validate_config
        errors.extend(validate_config(config, nested_get))
    if skill == "17-scrna-infer-grn":
        from grn_contract import validate_config
        errors.extend(validate_config(config, nested_get))
    stage = nested_get(config, "analysis.stage") or "differential"
    source = nested_get(config, "input.object") or nested_get(config, "input.counts_table") or nested_get(config, "input.differential_table") or nested_get(config, "enrichment.input_results") or nested_get(config, "input.path")
    if skill in {"11-scrna-run-differential-analysis", "12-scrna-run-pathway-enrichment"}:
        if skill == "12-scrna-run-pathway-enrichment" and stage != "enrichment_only":
            errors.append("12-scrna-run-pathway-enrichment requires analysis.stage=enrichment_only")
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
    if skill == "04-scrna-apply-qc-filter":
        if str(nested_get(config, "approval.status") or "").lower() != "approved":
            errors.append("approval.status must be exactly 'approved'")
        include = nested_get(config, "decision.include_all_true")
        if not isinstance(include, list) or not include or any(is_blank(x) for x in include):
            errors.append("decision.include_all_true must be a non-empty array")
        exclude = nested_get(config, "decision.exclude_any_true")
        if exclude is not None and (not isinstance(exclude, list) or any(is_blank(x) for x in exclude)):
            errors.append("decision.exclude_any_true must be an array")
        expected = nested_get(config, "decision.expected_retained_cells")
        if not isinstance(expected, int) or isinstance(expected, bool) or expected < 1:
            errors.append("decision.expected_retained_cells must be a positive integer")
        decision_table = nested_get(config, "input.decision_table")
        if decision_table and not Path(os.path.expandvars(os.path.expanduser(str(decision_table)))).exists():
            errors.append(f"decision table does not exist in this execution context: {decision_table}")
        object_name = str(nested_get(config, "output.object_name") or "filtered_object.rds")
        if not object_name.lower().endswith((".qs", ".rds")):
            errors.append("output.object_name must end in .qs or .rds")
        if (nested_get(config, "output.object_format") or "auto") not in {"auto", "qs", "rds"}:
            errors.append("output.object_format must be auto, qs or rds")
    if skill == "08-scrna-annotate-cells":
        action = str(nested_get(config, "workflow.action") or "")
        if action not in {"prepare_review", "apply_confirmed"}:
            errors.append("workflow.action must be prepare_review or apply_confirmed")
        if action == "prepare_review" and is_blank(nested_get(config, "metadata.cluster")) and nested_get(config, "clustering.compute_if_missing") is not True:
            errors.append("prepare_review requires metadata.cluster or clustering.compute_if_missing=true")
        if action == "apply_confirmed":
            for field in ("input.decisions", "metadata.reduction", "annotation.broad_column", "annotation.fine_column"):
                if is_blank(nested_get(config, field)):
                    errors.append(f"apply_confirmed requires {field}")
            decisions = nested_get(config, "input.decisions")
            if decisions and not Path(os.path.expandvars(os.path.expanduser(str(decisions)))).exists():
                errors.append(f"annotation decisions do not exist in this execution context: {decisions}")
            for field in ("input.review_record", "approval.review_run_id", "approval.decision_sha256", "approval.review_record_sha256", "metadata.cluster"):
                if is_blank(nested_get(config, field)):
                    errors.append(f"apply_confirmed requires {field}")
            if nested_get(config, "approval.status") != "approved":
                errors.append("annotation approval.status must be exactly 'approved'")
            if nested_get(config, "approval.source") != "human":
                errors.append("annotation approval.source must be exactly 'human'")
            approved_hash = nested_get(config, "approval.decision_sha256")
            if approved_hash is not None and (not isinstance(approved_hash, str) or not re.fullmatch(r"[0-9a-fA-F]{64}", approved_hash)):
                errors.append("annotation approval.decision_sha256 must be a SHA256 hex digest")
            review_hash = nested_get(config, "approval.review_record_sha256")
            if not isinstance(review_hash, str) or not re.fullmatch(r"[0-9a-fA-F]{64}", review_hash):
                errors.append("annotation approval.review_record_sha256 must be a SHA256 hex digest")
            review_value = nested_get(config, "input.review_record")
            if review_value:
                review_path = Path(os.path.expandvars(os.path.expanduser(str(review_value))))
                try:
                    review = json.loads(review_path.read_text(encoding="utf-8"))
                    if isinstance(review_hash, str) and sha256(review_path) != review_hash.lower():
                        errors.append("annotation review record SHA256 does not match the approved review version")
                    if not isinstance(review, dict) or review.get("schema_version") != 1 or review.get("kind") != "annotation_review":
                        errors.append("annotation review record has an unsupported schema or kind")
                    else:
                        if not review.get("review_run_id") or review.get("review_run_id") != nested_get(config, "approval.review_run_id"):
                            errors.append("annotation approval.review_run_id does not match the review record")
                        if review.get("cluster_column") != nested_get(config, "metadata.cluster"):
                            errors.append("annotation metadata.cluster does not match the reviewed cluster column")
                        bound_object = review.get("apply_object", {})
                        object_value = nested_get(config, "input.object")
                        object_path = Path(os.path.expandvars(os.path.expanduser(str(object_value)))) if object_value else None
                        if not isinstance(bound_object, dict) or not re.fullmatch(r"[0-9a-f]{64}", str(bound_object.get("sha256", ""))):
                            errors.append("annotation review record requires the reviewed apply_object SHA256")
                        elif object_path and object_path.is_file() and sha256(object_path) != bound_object["sha256"]:
                            errors.append("annotation input.object SHA256 does not match the reviewed apply object")
                        decision_path = Path(os.path.expandvars(os.path.expanduser(str(decisions)))) if decisions else None
                        if decision_path and decision_path.is_file():
                            if isinstance(approved_hash, str) and sha256(decision_path) != approved_hash.lower():
                                errors.append("annotation decisions SHA256 does not match the approved final decision version")
                            import csv
                            with decision_path.open(encoding="utf-8-sig", newline="") as stream:
                                decision_rows = list(csv.DictReader(stream, delimiter="\t"))
                            actual_clusters = [row.get("cluster") for row in decision_rows]
                            reviewed_clusters = review.get("cluster_ids")
                            if not isinstance(reviewed_clusters, list) or not reviewed_clusters or any(not isinstance(value, str) for value in reviewed_clusters):
                                errors.append("annotation review record requires string cluster_ids")
                            elif len(actual_clusters) != len(set(actual_clusters)) or set(actual_clusters) != set(reviewed_clusters):
                                errors.append("annotation decisions must cover exactly the reviewed clusters")
                except (OSError, ValueError, UnicodeError) as exc:
                    errors.append(f"annotation review binding cannot be read: {exc}")
            warnings.append("Annotation approval fields prevent content drift; they do not authenticate human authorization. Apply requires the user's explicit approval of this final decision version.")
    if skill == "06-scrna-preprocess-and-cluster":
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
    if skill == "11-scrna-run-differential-analysis" and stage != "enrichment_only":
        transform = nested_get(config, "analysis.pca_transform") or "vst"
        if transform not in ("vst", "log2_normalized"):
            errors.append("analysis.pca_transform must be vst or log2_normalized")
        analysis = config.get("analysis", {})
        if not isinstance(analysis, dict):
            errors.append("analysis must be an object")
            analysis = {}
        for field, minimum in (("min_samples_per_group", 2), ("min_cells_per_sample_population", 1),
                               ("min_total_count", 0), ("min_count_per_sample", 0), ("min_samples_expressed", 1)):
            if field in analysis:
                value = analysis[field]
                if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
                    errors.append(f"analysis.{field} must be an integer of at least {minimum}")
        for field in ("padj_threshold", "lfc_threshold"):
            if field in analysis:
                value = analysis[field]
                valid = (not isinstance(value, bool) and isinstance(value, (int, float))
                         and value == value and value not in (float("inf"), -float("inf")))
                if not valid or (not 0 < value < 1 if field == "padj_threshold" else value < 0):
                    errors.append(f"analysis.{field} must be a finite scalar " +
                                  ("between 0 and 1" if field == "padj_threshold" else "of at least 0"))
        if (nested_get(config, "analysis.method") or "pseudobulk_deseq2") == "pseudobulk_deseq2":
            if "design" in analysis:
                # Parse the complete supported grammar; operators inside quoted column
                # names are data, while ^, calls, dot expansion and interactions fail.
                design = analysis["design"]
                token = r"(?:`(?:[^`\\]|\\.)+`|(?:[A-Za-z]|\.(?![0-9]))[A-Za-z0-9._]*|1)"
                if not isinstance(design, str) or not re.fullmatch(r"\s*~\s*" + token + r"(?:\s*\+\s*" + token + r")*\s*", design):
                    errors.append("analysis.design supports only an intercept and explicit additive column main effects")
                else:
                    columns = []
                    for match in re.finditer(token, design[design.index("~") + 1:]):
                        column = match.group()
                        if column.startswith("`"):
                            column = re.sub(r"\\(.)", r"\1", column[1:-1])
                        elif column == "1":
                            continue
                        columns.append(column)
                    condition = nested_get(config, "metadata.condition")
                    declared = [condition] + (nested_get(config, "metadata.covariates") or [])
                    if condition not in columns or any(column not in declared or column == "." for column in columns):
                        errors.append("analysis.design requires metadata.condition as an independent main effect and only declared covariates")
            kind = nested_get(config, "analysis.counts_source.kind")
            if kind not in ("raw_umi", "raw_read"):
                errors.append("formal pseudobulk requires analysis.counts_source.kind = raw_umi or raw_read")
            assay = str(nested_get(config, "analysis.assay") or "RNA")
            if re.search(r"(^|[._-])(integrated|sct|harmony|corrected|scaled|normalized)([._-]|$)", assay, re.IGNORECASE):
                errors.append("formal pseudobulk requires an uncorrected raw-count assay")
            reference = nested_get(config, "analysis.counts_source.reference_object")
            if reference and not Path(os.path.expandvars(os.path.expanduser(str(reference)))).is_file():
                errors.append(f"raw-count reference does not exist in this execution context: {reference}")
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
    if skill in {"11-scrna-run-differential-analysis", "12-scrna-run-pathway-enrichment"}:
        universe_mode = nested_get(config, "enrichment.universe_mode") or "multiple_testing_eligible"
        if universe_mode not in ("multiple_testing_eligible", "tested"):
            errors.append("enrichment.universe_mode must be multiple_testing_eligible or tested")
    if skill == "13-scrna-test-cell-abundance":
        has_object = bool(nested_get(config, "input.object"))
        has_counts = bool(nested_get(config, "input.counts_table"))
        if has_object == has_counts:
            errors.append("provide exactly one of input.object or input.counts_table")
        methods = nested_get(config, "analysis.methods")
        supported = {"propeller", "sccomp", "sccoda", "milo", "dcats"}
        if not isinstance(methods, list) or not methods:
            errors.append("analysis.methods must be a non-empty array")
        else:
            unknown = {str(method).lower() for method in methods} - supported
            if unknown:
                errors.append("unsupported cell-abundance methods: " + ", ".join(sorted(unknown)))
            if len(methods) != len({str(method).lower() for method in methods}):
                errors.append("analysis.methods must not contain duplicates")
        comparisons = config.get("comparisons")
        if not isinstance(comparisons, list) or not comparisons:
            errors.append("comparisons must be a non-empty array")
        else:
            comparison_ids = []
            for index, comparison in enumerate(comparisons):
                if not isinstance(comparison, dict) or is_blank(comparison.get("id")) or is_blank(comparison.get("numerator")) or is_blank(comparison.get("denominator")):
                    errors.append(f"comparison {index + 1} requires id, numerator, and denominator")
                elif comparison["numerator"] == comparison["denominator"]:
                    errors.append(f"comparison {index + 1} numerator and denominator must differ")
                if isinstance(comparison, dict) and not is_blank(comparison.get("id")):
                    comparison_ids.append(str(comparison["id"]))
            if len(comparison_ids) != len(set(comparison_ids)):
                errors.append("comparison ids must be unique")
        denominator_mode = nested_get(config, "analysis.denominator.mode")
        if denominator_mode not in {"all_input_cells", "selected_cell_types"}:
            errors.append("analysis.denominator.mode must be all_input_cells or selected_cell_types")
        if denominator_mode == "selected_cell_types":
            include = nested_get(config, "analysis.denominator.include")
            if not isinstance(include, list) or not include or any(is_blank(value) for value in include):
                errors.append("selected_cell_types denominator requires a non-empty analysis.denominator.include array")
        if is_blank(nested_get(config, "analysis.denominator.description")):
            errors.append("analysis.denominator.description is required so relative abundance has an explicit interpretation")
        fdr = nested_get(config, "analysis.fdr")
        if fdr is not None and (not isinstance(fdr, (int, float)) or isinstance(fdr, bool) or not 0 < fdr < 1):
            errors.append("analysis.fdr must be a number between 0 and 1")
        min_samples = nested_get(config, "analysis.min_samples_per_group")
        if min_samples is not None and (not isinstance(min_samples, int) or isinstance(min_samples, bool) or min_samples < 2):
            errors.append("analysis.min_samples_per_group must be an integer of at least 2")
        min_cells = nested_get(config, "analysis.min_cells_per_sample")
        if min_cells is not None and (not isinstance(min_cells, int) or isinstance(min_cells, bool) or min_cells < 1):
            errors.append("analysis.min_cells_per_sample must be a positive integer")
        if (nested_get(config, "analysis.min_cells_policy") or "audit_only") not in {"audit_only", "exclude_samples", "stop"}:
            errors.append("analysis.min_cells_policy must be audit_only, exclude_samples or stop")
        covariate_types = nested_get(config, "metadata.covariate_types") or {}
        covariates = nested_get(config, "metadata.covariates") or []
        if not isinstance(covariate_types, dict) or any(key not in covariates or kind not in {"continuous", "categorical"} for key, kind in covariate_types.items()):
            errors.append("metadata.covariate_types must map configured covariates to continuous or categorical")
        method_names = {str(method).lower() for method in methods} if isinstance(methods, list) else set()
        if "propeller" in method_names:
            transform = nested_get(config, "method_options.propeller.transform")
            if transform is not None and transform not in {"logit", "asin"}:
                errors.append("method_options.propeller.transform must be logit or asin")
        if "sccoda" in method_names:
            references = nested_get(config, "method_options.sccoda.reference_cell_types")
            if not isinstance(references, list) or not references or any(is_blank(value) for value in references):
                errors.append("sccoda requires a non-empty method_options.sccoda.reference_cell_types array")
        if isinstance(methods, list) and "milo" in {str(method).lower() for method in methods}:
            if not nested_get(config, "input.object"):
                errors.append("milo requires input.object; an aggregated counts table is insufficient")
            if is_blank(nested_get(config, "method_options.milo.reduction")):
                errors.append("milo requires method_options.milo.reduction")
            for field in ("k", "d"):
                value = nested_get(config, f"method_options.milo.{field}")
                if value is not None and (not isinstance(value, int) or isinstance(value, bool) or value < 1):
                    errors.append(f"method_options.milo.{field} must be a positive integer")
            prop = nested_get(config, "method_options.milo.prop")
            if prop is not None and (not isinstance(prop, (int, float)) or isinstance(prop, bool) or not 0 < prop <= 1):
                errors.append("method_options.milo.prop must be in (0, 1]")
        similarity = nested_get(config, "method_options.dcats.similarity_matrix")
        if "dcats" in method_names and similarity:
            similarity_path = Path(os.path.expandvars(os.path.expanduser(str(similarity))))
            if not similarity_path.is_file():
                errors.append(f"DCATS similarity matrix does not exist: {similarity_path}")
    if skill == "14-scrna-visualize-cell-composition":
        has_object = bool(nested_get(config, "input.object"))
        has_counts = bool(nested_get(config, "input.counts_table"))
        if has_object == has_counts:
            errors.append("provide exactly one of input.object or input.counts_table")
        variables = nested_get(config, "composition.variables")
        if not isinstance(variables, list) or not variables:
            errors.append("composition.variables must be a non-empty array")
        else:
            seen = set()
            for index, variable in enumerate(variables):
                if not isinstance(variable, dict) or is_blank(variable.get("column")):
                    errors.append(f"composition variable {index + 1} requires column")
                elif variable["column"] in seen:
                    errors.append("composition.variables must not contain duplicate columns")
                else:
                    seen.add(variable["column"])
                kind = str(variable.get("kind", "metadata")).lower() if isinstance(variable, dict) else ""
                if kind not in {"cluster", "annotation", "state", "metadata"}:
                    errors.append(f"composition variable {index + 1} kind must be cluster, annotation, state, or metadata")
        mode = nested_get(config, "composition.denominator.mode") or "all_input_cells"
        if mode not in {"all_input_cells", "selected_parent", "selected_cell_types"}:
            errors.append("composition.denominator.mode must be all_input_cells, selected_parent, or selected_cell_types")
        if mode == "selected_parent" and is_blank(nested_get(config, "composition.parent_column")):
            errors.append("selected_parent requires composition.parent_column")
        if mode in {"selected_parent", "selected_cell_types"} and is_blank(nested_get(config, "composition.denominator.include")):
            errors.append("selected denominator requires composition.denominator.include")
        modes = config.get("modes", ["composition_summary"])
        if not isinstance(modes, list) or not modes or any(str(x) not in {"composition_summary", "batch_diagnostic", "hierarchical_composition"} for x in modes):
            errors.append("modes must contain composition_summary, batch_diagnostic, or hierarchical_composition")
        figure_format = nested_get(config, "plots.figure_format") or "png"
        if figure_format not in {"png", "pdf", "both"}:
            errors.append("plots.figure_format must be png, pdf, or both")
    if skill == "15-scrna-visualize-gene":
        selection = config.get("differential_selection", {})
        if not isinstance(selection, dict) or any(value is not None and (isinstance(value, (list, dict)) or is_blank(value)) for value in selection.values()):
            errors.append("differential_selection must map column names to single non-empty values")
        genes = config.get("genes")
        if not isinstance(genes, list) or not genes:
            errors.append("genes must be a non-empty array")
        else:
            symbols = []
            for index, gene in enumerate(genes):
                if isinstance(gene, str):
                    symbol = gene
                elif isinstance(gene, dict):
                    symbol = gene.get("symbol")
                else:
                    symbol = None
                if is_blank(symbol):
                    errors.append(f"gene {index + 1} requires symbol")
                else:
                    symbols.append(str(symbol))
            if len(symbols) != len(set(symbols)):
                errors.append("genes must not contain duplicate symbols")
        figure_format = nested_get(config, "plots.figure_format") or "png"
        if figure_format not in {"png", "pdf", "both"}:
            errors.append("plots.figure_format must be png, pdf, or both")
    if skill == "10-scrna-score-programs":
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
                if name == "none" and method.get("id", "none") != "none":
                    errors.append("The uncorrected baseline scenario ID must be none")
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
            if not errors:
                ids = []
                for method in methods:
                    name = str(method["name"]).lower()
                    base = method.get("id") or (name + "__reduction_" + str(method.get("reduction", "")) if name == "precomputed" else name)
                    grid = method.get("parameter_grid", {})
                    for values in itertools.product(*grid.values()):
                        def text(value):
                            if isinstance(value, bool): return str(value).lower()
                            if isinstance(value, (int, float)): return format(Decimal(format(value, ".15g")), "f")
                            return str(value)
                        suffix = "__".join(str(key) + "_" + text(value) for key, value in zip(grid, values))
                        ids.append(re.sub(r"[^A-Za-z0-9_.-]+", "_", str(base) + ("__" + suffix if suffix else "")).strip("_"))
                if "none" not in scenario_names: ids.append("none")
                if any(not value for value in ids) or len(ids) != len(set(ids)):
                    errors.append("Integration scenario IDs must be non-empty and unique after sanitization")
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
        warnings.append("registered pixi executor is unavailable; dry-run remains available and system interpreters will not be used")
    return errors, warnings


def make_manifest(skill, config, config_path, errors, warnings):
    source = nested_get(config, "input.object") or nested_get(config, "input.counts_table") or nested_get(config, "input.differential_table") or nested_get(config, "enrichment.input_results") or nested_get(config, "input.path")
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
        "environment": environment_record(skill, config),
        "errors": errors,
        "warnings": warnings,
        "status": "blocked" if errors else "ready",
    }


def run_manifest_name(skill, config):
    if skill == "06-scrna-preprocess-and-cluster":
        return "run_manifest_finalize.json" if nested_get(config, "workflow.action") == "finalize_resolution" else "run_manifest_preprocess.json"
    return "run_manifest.json"


def archive_previous_output(skill, config, config_path, output_dir, run_id):
    """Keep an earlier analysis intact outside the current result directory."""
    if skill not in {"11-scrna-run-differential-analysis", "12-scrna-run-pathway-enrichment", "13-scrna-test-cell-abundance"}:
        return None
    if not output_dir.exists() or not any(output_dir.iterdir()):
        return None
    def strings(value):
        if isinstance(value, dict):
            for child in value.values(): yield from strings(child)
        elif isinstance(value, list):
            for child in value: yield from strings(child)
        elif isinstance(value, str): yield value
    protected = [str(config_path)] + list(strings(config.get("input", {})))
    if nested_get(config, "enrichment.input_results"):
        protected.append(nested_get(config, "enrichment.input_results"))
    resolved_output = output_dir.resolve()
    for value in protected:
        candidate = Path(os.path.expandvars(os.path.expanduser(value))).resolve()
        if candidate.exists() and (candidate == resolved_output or resolved_output in candidate.parents):
            raise ValueError("Input/config is inside an existing output directory; use a separate output_dir to preserve it: " + str(candidate))
    archive_root = output_dir.parent / ("." + output_dir.name + "-previous-runs")
    archive_root.mkdir(exist_ok=True)
    archived = archive_root / run_id
    output_dir.rename(archived)
    archived_provenance = archived / "_provenance"
    archived_provenance.mkdir(exist_ok=True)
    (archived_provenance / "archive_record.json").write_text(json.dumps({
        "original_output_dir": str(resolved_output), "archived_output_dir": str(archived.resolve()),
        "archived_at": dt.datetime.now(dt.timezone.utc).isoformat(), "superseding_run_id": run_id,
        "historical_path_resolution": "Replace the original_output_dir prefix with archived_output_dir; prior manifests are preserved unchanged."
    }, indent=2) + "\n", encoding="utf-8")
    return str(archived)


def write_execution_manifest(path, record):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(record, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    temporary.replace(path)


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
    if skill == "06-scrna-preprocess-and-cluster":
        action = nested_get(config, "workflow.action") or "run"
        action_suffix = ".finalize" if action == "finalize_resolution" else ".scan"
    manifest_path = args.manifest or args.config.parent / "_provenance" / f"{skill}{action_suffix}.manifest.json"
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
    run_id = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%S") + "-" + uuid.uuid4().hex[:12]
    try:
        previous_output = archive_previous_output(skill, config, args.config, output_dir, run_id)
    except (OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    output_dir.mkdir(parents=True, exist_ok=True)
    if previous_output and output_dir.resolve() in manifest_path.resolve().parents:
        manifest_path.parent.mkdir(parents=True, exist_ok=True)
        manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    delivery_before = _delivery.snapshot(output_dir)
    technical_dir = output_dir / "_provenance"
    technical_dir.mkdir(parents=True, exist_ok=True)
    log_path = technical_dir / "run.log"
    started = dt.datetime.now(dt.timezone.utc).isoformat()
    command = [executable] + [str(x) for x in argv[1:]]
    execution_path = technical_dir / run_manifest_name(skill, config)
    execution = {"schema_version": 2, "skill": skill, "run_id": run_id, "project_id": nested_get(config, "project.id"),
        "started_at": started, "status": "running", "exit_status": None, "output_dir": str(output_dir),
        "config": {"path": str(args.config.resolve()), "sha256": sha256(args.config)}, "artifacts": [],
        "previous_output": previous_output, "executor_argv": command, "environment": manifest["environment"]}
    write_execution_manifest(execution_path, execution)
    returncode = 1
    process = None
    try:
        with log_path.open("a", encoding="utf-8") as log:
            log.write(f"[{started}] START {' '.join(command)}\n")
            child_env = os.environ.copy()
            child_env["SCRNA_ACTIVE_SKILL"] = skill
            child_env["SCRNA_RUN_ID"] = run_id
            child_env["SCRNA_PIXI_ROOT"] = str(environment_root(config))
            if not config.get("executor") and manifest["resolved_rscript"]:
                child_env["PATH"] = str(Path(manifest["resolved_rscript"]).parent) + os.pathsep + child_env.get("PATH", "")
            process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1, env=child_env)
            assert process.stdout is not None
            for line in process.stdout:
                sys.stdout.write(line)
                sys.stdout.flush()
                log.write(line)
                log.flush()
            returncode = process.wait()
            process.stdout.close()
    except KeyboardInterrupt:
        returncode = 130
        if process and process.poll() is None:
            process.terminate()
            try: process.wait(timeout=10)
            except subprocess.TimeoutExpired: process.kill(); process.wait()
    except OSError as exc:
        print(f"ERROR: executor failed: {exc}", file=sys.stderr)
    finally:
        finished = dt.datetime.now(dt.timezone.utc).isoformat()
        with log_path.open("a", encoding="utf-8") as log:
            log.write(f"[{finished}] EXIT {returncode}\n")
        try:
            record = json.loads(execution_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            record = execution
        record.update(run_id=run_id, status="completed" if returncode == 0 else "failed", exit_status=returncode,
            executor_started_at=started, executor_finished_at=finished, executor_argv=command, previous_output=previous_output,
            environment=manifest["environment"])
        returncode = _delivery.supervise(output_dir, delivery_before, skill, config, returncode, record, run_id)
        write_execution_manifest(execution_path, record)
    return returncode
