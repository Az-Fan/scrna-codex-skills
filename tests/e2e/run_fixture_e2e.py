#!/usr/bin/env python3
"""Run all released skills against the deterministic tiny Seurat fixture.

This is an execution test, not a biological validation.  It exercises the
installed, self-contained skill bundles in the project's existing pixi
environments and verifies their minimum output contracts.
"""
from __future__ import annotations

import argparse
import copy
import csv
import gzip
import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path


SKILL_ENVS = {
    "01-scrna-standardize-input": "01-scrna-qc",
    "02-scrna-calculate-qc-metrics": "01-scrna-qc",
    "03-scrna-review-qc": "01-scrna-qc",
    "04-scrna-apply-qc-filter": "01-scrna-qc",
    "05-scrna-benchmark-integration": "03-integration",
    "06-scrna-preprocess-and-cluster": "03-integration",
    "07-scrna-find-cluster-markers": "02-annotation",
    "08-scrna-annotate-cells": "02-annotation",
    "09-scrna-export-subset": "02-annotation",
    "10-scrna-score-programs": "05-pathway_program",
    "11-scrna-run-differential-analysis": "06-deg-analysis",
    "12-scrna-run-pathway-enrichment": "06-deg-analysis",
    "13-scrna-test-cell-abundance": "07-cell-abundance",
}

EXPECTED = {
    "01-scrna-standardize-input": ["cell_metadata.tsv", "samples.tsv", "field_mapping.json", "run_manifest.json"],
    "02-scrna-calculate-qc-metrics": ["qc_metrics_object.rds", "metadata.tsv.gz", "metric_status.tsv", "run_manifest.json"],
    "03-scrna-review-qc": ["qc_summary_by_sample.tsv", "threshold_review.tsv", "qc_atlas.pdf", "run_manifest.json"],
    "04-scrna-apply-qc-filter": ["filtered_object.rds", "cell_filter_decisions.tsv.gz", "filter_summary_by_sample.tsv", "filter_summary_by_condition.tsv", "filter_decision_counts.tsv", "approved_filter_record.tsv", "session_info.txt", "run_manifest.json"],
    "05-scrna-benchmark-integration": ["method_runs_r.tsv", "integration_benchmark_object.qs", "run_manifest.json"],
    "06-scrna-preprocess-and-cluster": ["preprocessed_clustered_object.qs", "cell_assignments.tsv", "standard_cluster_sizes.tsv", "standard_sample_cluster_counts.tsv", "standard_umap_diagnostics.pdf", "scenario_summary.tsv", "workflow_state.json", "session_info.txt", "run.log", "run_manifest_preprocess.json", "run_manifest_finalize.json"],
    "07-scrna-find-cluster-markers": ["cluster_markers.tsv", "top_cluster_markers.tsv", "cluster_marker_summary.tsv", "top_marker_dotplot.pdf", "run_manifest.json"],
    "08-scrna-annotate-cells": ["clustered_object.qs", "cluster_markers.tsv", "annotation_review.tsv", "cluster_umap.pdf", "cluster_sample_umap.pdf", "canonical_marker_dotplot.pdf", "annotated_object.qs", "cell_annotations.tsv", "annotation_summary.tsv", "annotated_umap.pdf", "cluster_sample_condition_umap.pdf", "session_info.txt", "run_manifest.json"],
    "09-scrna-export-subset": ["subset_object.qs", "subset_counts.mtx", "subset_metadata.tsv", "subset_summary.tsv", "features.tsv", "barcodes.tsv", "run_manifest.json"],
    "10-scrna-score-programs": ["signature_coverage.tsv", "assay_feature_mapping.tsv", "score_summary.tsv", "task_manifest.json", "run_manifest.json"],
    "11-scrna-run-differential-analysis": ["design_audit.tsv", "task_status.tsv", "all_comparisons.tsv", "sessionInfo.txt", "run_manifest.json"],
    "12-scrna-run-pathway-enrichment": ["task_status.tsv", "sessionInfo.txt", "run_manifest.json"],
    "13-scrna-test-cell-abundance": ["sample_cell_counts.tsv", "sample_cell_proportions.tsv", "design_audit.tsv", "cell_type_eligibility.tsv", "task_status.tsv", "all_method_results.tsv", "method_concordance.tsv", "sample_composition.pdf", "cell_type_proportions_by_condition.pdf", "sample_proportion_heatmap.pdf", "sessionInfo.txt", "run_manifest.json"],
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def call(argv, *, env=None, cwd=None):
    print("+", " ".join(map(str, argv)), flush=True)
    subprocess.run([str(x) for x in argv], check=True, env=env, cwd=cwd)


def write_json(path: Path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--env-root", type=Path, default=Path("/home/faz_laptop/projects/scrna_envs"))
    parser.add_argument("--pixi", type=Path, default=Path.home() / ".pixi/bin/pixi")
    parser.add_argument("--start-at", choices=list(SKILL_ENVS), help="resume at this skill without deleting earlier outputs")
    args = parser.parse_args()
    repo = args.repo.resolve()
    output_root = repo / "test-output/e2e"
    installed = repo / "test-output/installed"
    configs = output_root / "configs"
    fixture = repo / "tests/fixtures/tiny_scrna.rds"
    qc_fixture = repo / "tests/fixtures/tiny_scrna_multilayer.rds"
    if output_root.exists() and args.start_at is None:
        shutil.rmtree(output_root)
    output_root.mkdir(parents=True, exist_ok=True)

    annotation_r = args.env_root / "02-annotation/.pixi/envs/default/bin/Rscript"
    call([annotation_r, repo / "tests/fixtures/create_fixture.R", fixture, qc_fixture], cwd=repo)
    fixture_hash = sha256(fixture)
    qc_fixture_hash = sha256(qc_fixture)
    call([sys.executable, repo / "scripts/install_skills.py", "--target", installed, "--force"], cwd=repo)

    decision_table = configs / "approved_qc_decisions.tsv"
    decision_table.parent.mkdir(parents=True, exist_ok=True)
    with decision_table.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["cell_id", "in_scope", "passes_core_qc", "exclude_severe", "reason"])
        for index in range(1, 81):
            writer.writerow([f"cell_{index:03d}", "TRUE", "TRUE", "TRUE" if index == 80 else "FALSE", "approved_fixture_exclusion" if index == 80 else "retained"])

    qc_project = args.env_root / "01-scrna-qc"
    base = {"project": {"id": "tiny_fixture"}, "input": {"object": str(fixture)}}
    definitions = {
        "01-scrna-standardize-input": {**base, "input": {"path": str(fixture), "format": "auto"}, "metadata": {"sample": "sample_label"}},
        "02-scrna-calculate-qc-metrics": {
            "pixi": {"executable": str(args.pixi), "project": str(qc_project), "environment": "default"},
            "input": {"type": "seurat", "object": str(qc_fixture), "starsolo_dir": None, "gtf_file": None},
            "metadata": {"sample": "sample_label", "batch": "batch_id"}, "samples": [],
            "ambient_rna": {"method": "skip", "cluster_column": None},
            "parallel": {"workers": 2},
            "parameters": {"species": "mouse", "assay": "RNA", "min_genes_hq": 3, "mito_pattern": "^mt-", "doublet_artificial_fraction": 0.2, "seed": 1},
        },
        "03-scrna-review-qc": {
            "pixi": {"executable": str(args.pixi), "project": str(qc_project), "environment": "default"},
            "input": {"object": str(output_root / "02-scrna-calculate-qc-metrics/qc_metrics_object.rds")},
            "metadata": {"sample": "sample_label", "condition": "condition", "batch": "batch_id", "cluster": "seurat_clusters", "annotation": "cell_type"},
            "output": {"detail_level": "compact"},
            "thresholds": {},
        },
        "04-scrna-apply-qc-filter": {
            **base,
            "input": {"object": str(fixture), "decision_table": str(decision_table), "assay": "RNA"},
            "metadata": {"sample": "sample_label", "condition": "condition"},
            "decision": {"cell_id_column": "cell_id", "include_all_true": ["in_scope", "passes_core_qc"], "exclude_any_true": ["exclude_severe"], "reason_column": "reason", "expected_retained_cells": 79, "carry_columns": ["passes_core_qc", "exclude_severe"]},
            "approval": {"status": "approved", "approved_at": "2026-09-02", "note": "deterministic E2E fixture"},
            "output": {"object_name": "filtered_object.rds"},
        },
        "05-scrna-benchmark-integration": {**base, "input": {"object": str(fixture), "assay": "RNA"}, "metadata": {"sample": "sample_label", "batch_variables": ["batch_id"], "biological_labels": ["cell_type"], "condition": "condition"}, "benchmark": {"methods": [{"name": "none"}, {"name": "harmony", "parameter_grid": {"theta": [2]}}], "dims": [1, 2, 3, 4, 5], "neighbors": {"k": 10}, "seed": 1, "python_argv_prefix": [str(args.env_root / "03-integration/.pixi/envs/scvi/bin/python3.12")]}, "metrics": {"batch_removal": ["batch_asw"], "biological_conservation": ["label_asw", "nmi", "ari"]}, "plots": ["score_heatmap", "metric_tradeoff", "umap_by_batch", "umap_by_label"], "gene_programs": {}, "scoring": {"enabled": False, "batch_weight": 0.3, "biology_weight": 0.7}},
        "06-scrna-preprocess-and-cluster": {
            **base,
            "workflow": {"mode": "guided", "action": "run"},
            "input": {"object": str(fixture), "assay": "RNA", "qc_status": "filtered"},
            "metadata": {"sample": "sample_label", "condition": "condition", "batch": "batch_id"},
            "preprocessing": {"normalization_method": "LogNormalize", "scale_factor": 10000, "variable_feature_method": "vst", "n_variable_features": 30, "pca_npcs": 10, "dims": [1, 2, 3, 4, 5], "seed": 1},
            "cell_cycle": {"score": False, "species": "mouse"},
            "scenarios": [{"name": "standard", "regress_variables": [], "harmony": {"enabled": False}, "neighbors": {"k_param": 7}, "umap": {"n_neighbors": 9, "min_dist": 0.15, "metric": "euclidean", "method": "uwot"}, "clustering": {"mode": "scan", "resolutions": [0.2, 0.4], "selection": "review", "stability": {"n_subsamples": 2, "subsample_frac": 0.8, "min_clusters": 2}}}],
            "plots": {"group_by": ["sample_label", "condition", "batch_id"], "pt_size": 0.2, "preview_png": False},
            "output": {"object_format": "qs", "drop_scale_data": True},
        },
        "07-scrna-find-cluster-markers": {**base, "metadata": {"cluster": "seurat_clusters"}, "analysis": {"assay": "RNA", "test_use": "wilcox", "only_pos": True, "logfc_threshold": 0.1, "min_pct": 0.1, "min_diff_pct": -0.1, "return_thresh": 1.0, "max_cells_per_ident": None, "random_seed": 1, "join_layers": True, "normalize_if_missing": True}, "reporting": {"top_n": 10, "dotplot_top_n": 3}},
        "08-scrna-annotate-cells": {**base, "workflow": {"action": "prepare_review"}, "input": {"object": str(fixture), "markers": str(output_root / "07-scrna-find-cluster-markers/cluster_markers.tsv")}, "metadata": {"sample": "sample_label", "condition": "condition", "cluster": "seurat_clusters", "reduction": "umap"}, "clustering": {"compute_if_missing": False}, "markers": {"assay": "RNA", "canonical": {"endothelial": ["Kdr", "Pecam1", "Cdh5"], "fibroblast": ["Col1a1", "Col3a1", "Dcn"]}}},
        "09-scrna-export-subset": {**base, "metadata": {"sample": "sample_label", "cell_type": "cell_type"}, "subset": {"include": ["Endothelial"]}},
        "10-scrna-score-programs": {**base, "input": {"object": str(fixture), "assay": "RNA", "layer": "counts"}, "species": "mouse", "tasks": [{"name": "vascular_program", "method": "addmodulescore", "gene_sets": {"source": "inline", "sets": {"vascular": ["Kdr", "Pecam1", "Cdh5"]}}, "coverage": {"min_genes": 3, "min_fraction": 1.0, "on_insufficient": "error"}, "parameters": {"normalize_if_missing": True, "nbin": 4, "ctrl": 2}}], "summarize_by": ["sample_label", "condition", "cell_type"], "random_seed": 1, "cores": 1, "cache": {"enabled": False}, "output": {"object_format": "rds"}},
        "11-scrna-run-differential-analysis": {**base, "metadata": {"sample": "sample_label", "condition": "condition", "covariates": []}, "population": {"mode": "all", "include": [], "exclude": []}, "comparisons": [{"id": "stz_vs_control", "numerator": "stz", "denominator": "control"}], "analysis": {"stage": "differential", "method": "pseudobulk_deseq2", "assay": "RNA", "design": "~ condition", "min_cells_per_sample_population": 10, "min_samples_per_group": 2, "min_total_count": 1, "min_count_per_sample": 1, "min_samples_expressed": 2, "padj_threshold": 0.1, "lfc_threshold": 0.1, "lfc_shrink": False}, "plots": {"top_genes": 10}, "enrichment": {"enabled": False, "species": "mouse", "gene_id_type": "SYMBOL"}},
        "12-scrna-run-pathway-enrichment": {"project": {"id": "tiny_fixture_enrichment"}, "random_seed": 1, "input": {"differential_table": str(output_root / "11-scrna-run-differential-analysis/all_comparisons.tsv")}, "analysis": {"stage": "enrichment_only", "padj_threshold": 1.0, "lfc_threshold": 0.0}, "enrichment": {"enabled": True, "species": "mouse", "gene_id_type": "SYMBOL", "databases": ["GO_BP"], "min_input_genes": 1, "min_gene_set_size": 1, "max_gene_set_size": 500, "plot_top_terms": 5, "plot_label_width": 30, "plot_terms_per_page": 8}},
        "13-scrna-test-cell-abundance": {
            "project": {"id": "tiny_fixture_abundance"},
            "input": {"counts_table": str(repo / "tests/fixtures/cell_abundance_counts.tsv"), "count_column": "n_cells"},
            "metadata": {"sample": "sample_label", "condition": "condition", "cell_type": "cell_type", "covariates": []},
            "comparisons": [{"id": "case_vs_control", "numerator": "case", "denominator": "control"}],
            "analysis": {"methods": ["propeller", "dcats"], "denominator": {"mode": "all_input_cells", "description": "All annotated fixture populations"}, "min_samples_per_group": 3, "min_cells_per_sample": 20, "fdr": 0.1, "random_seed": 13},
            "method_options": {"propeller": {"transform": "logit", "robust": True, "trend": False}, "dcats": {"similarity_matrix": None, "reference_cell_types": []}},
            "runtime": {"pixi_root": str(args.env_root)},
        },
    }

    report = {"fixture_sha256_before": fixture_hash, "skills": {}}
    started = args.start_at is None
    for skill, config in definitions.items():
        if skill == args.start_at:
            started = True
        if not started:
            continue
        out = output_root / skill
        config["output_dir"] = str(out)
        config_path = configs / f"{skill}.json"
        manifest_path = configs / f"{skill}.manifest.json"
        write_json(config_path, config)
        launcher = installed / skill / "scripts/run.py"
        if skill in {"02-scrna-calculate-qc-metrics", "03-scrna-review-qc"}:
            command = [sys.executable, launcher, "--config", config_path]
            call(command)
            call(command + ["--execute"])
        else:
            env_dir = args.env_root / SKILL_ENVS[skill] / ".pixi/envs/default/bin"
            env = os.environ.copy()
            env["PATH"] = str(env_dir) + os.pathsep + env.get("PATH", "")
            command = [sys.executable, launcher, "--config", config_path, "--manifest", manifest_path]
            if skill == "06-scrna-preprocess-and-cluster":
                blocked_config = copy.deepcopy(config)
                blocked_config["input"]["qc_status"] = "unfiltered"
                blocked_config["input"].pop("allow_unfiltered", None)
                blocked_path = configs / f"{skill}.unfiltered-blocked.json"
                blocked_manifest = configs / f"{skill}.unfiltered-blocked.manifest.json"
                write_json(blocked_path, blocked_config)
                blocked = subprocess.run(
                    [sys.executable, launcher, "--config", blocked_path, "--manifest", blocked_manifest],
                    text=True, capture_output=True, env=env,
                )
                if blocked.returncode != 2 or "input.allow_unfiltered=true" not in blocked.stderr:
                    raise RuntimeError(f"{skill}: unfiltered input was not blocked correctly: rc={blocked.returncode}; stderr={blocked.stderr}")
            call(command, env=env)
            call(command + ["--execute"], env=env)
        if skill == "08-scrna-annotate-cells":
            decisions_path = configs / "confirmed_annotation_decisions.tsv"
            with decisions_path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(["cluster", "annotation_broad", "annotation_fine", "annotation_state", "decision", "confidence", "evidence", "conflicts", "sample_bias", "qc_flag"])
                writer.writerow(["0", "Endothelial", "Endothelial", "baseline", "confirmed", "high", "Kdr;Pecam1;Cdh5", "none", "none", "none"])
                writer.writerow(["1", "Stromal", "Fibroblast", "baseline", "confirmed", "high", "Col1a1;Col3a1;Dcn", "none", "none", "none"])
            apply_config = {
                "project": {"id": "tiny_fixture"},
                "workflow": {"action": "apply_confirmed"},
                "input": {"object": str(fixture), "decisions": str(decisions_path)},
                "metadata": {"sample": "sample_label", "condition": "condition", "cluster": "seurat_clusters", "reduction": "umap"},
                "annotation": {"broad_column": "annotation_broad", "fine_column": "annotation_fine", "state_column": "annotation_state", "decision_column": "decision", "confirmed_values": ["confirmed"]},
                "output_dir": str(out),
            }
            apply_path = configs / f"{skill}.apply.json"
            apply_manifest = configs / f"{skill}.apply.manifest.json"
            write_json(apply_path, apply_config)
            apply_command = [sys.executable, launcher, "--config", apply_path, "--manifest", apply_manifest]
            call(apply_command, env=env)
            call(apply_command + ["--execute"], env=env)
        if skill == "06-scrna-preprocess-and-cluster":
            scan_state = json.loads((out / "workflow_state.json").read_text(encoding="utf-8"))
            if scan_state.get("status") != "awaiting_resolution_confirmation":
                raise RuntimeError(f"{skill}: guided scan did not pause for confirmation: {scan_state}")
            if (out / "cell_assignments.tsv").exists():
                raise RuntimeError(f"{skill}: review scan wrote a meaningless final cell_assignments.tsv")
            scan_required = ["preprocessed_clustered_object.qs", "standard_resolution_stability.tsv", "standard_umap_clusters_by_resolution.png", "run_manifest_preprocess.json"]
            scan_missing = [name for name in scan_required if not (out / name).is_file() or not (out / name).stat().st_size]
            if scan_missing:
                raise RuntimeError(f"{skill}: scan artifacts missing/empty: {scan_missing}")
            finalize = {
                "project": {"id": "tiny_fixture"},
                "workflow": {"mode": "guided", "action": "finalize_resolution"},
                "input": {"object": str(out / "preprocessed_clustered_object.qs"), "assay": "RNA", "qc_status": "filtered"},
                "metadata": {"sample": "sample_label", "condition": "condition", "batch": "batch_id"},
                "finalize": {"scenario": "standard", "resolution": 0.4},
                "plots": {"group_by": ["sample_label", "condition", "batch_id"], "pt_size": 0.2, "preview_png": False},
                "output": {"object_format": "qs", "drop_scale_data": True},
                "output_dir": str(out),
            }
            finalize_config = configs / f"{skill}.finalize.json"
            finalize_manifest = configs / f"{skill}.finalize.manifest.json"
            write_json(finalize_config, finalize)
            finalize_command = [sys.executable, launcher, "--config", finalize_config, "--manifest", finalize_manifest]
            call(finalize_command, env=env)
            call(finalize_command + ["--execute"], env=env)
        missing = [name for name in EXPECTED[skill] if not (out / name).is_file() or (out / name).stat().st_size == 0]
        if skill == "01-scrna-standardize-input" and not any(path.is_file() and path.stat().st_size for path in out.glob("standardized_object.*")):
            missing.append("standardized_object.(qs|rds)")
        if missing:
            raise RuntimeError(f"{skill}: missing/empty required artifacts: {missing}")
        run_manifest_name = "run_manifest_finalize.json" if skill == "06-scrna-preprocess-and-cluster" else "run_manifest.json"
        manifest = json.loads((out / run_manifest_name).read_text(encoding="utf-8"))
        if manifest.get("skill") != skill:
            raise RuntimeError(f"{skill}: run manifest skill mismatch: {manifest.get('skill')!r}")
        if skill == "02-scrna-calculate-qc-metrics":
            with gzip.open(out / "metadata.tsv.gz", "rt", encoding="utf-8", newline="") as handle:
                metadata_rows = list(csv.DictReader(handle, delimiter="\t"))
            if len(metadata_rows) != 80:
                raise RuntimeError(f"{skill}: multi-layer input produced {len(metadata_rows)} cells instead of 80")
            with (out / "metric_status.tsv").open(encoding="utf-8", newline="") as handle:
                ambient_rows = [row for row in csv.DictReader(handle, delimiter="\t") if row["metric"] == "ambient_frac_decontx"]
            if len(ambient_rows) != 4 or any(row["status"] != "skipped" or row["reason"] != "disabled_by_config" for row in ambient_rows):
                raise RuntimeError(f"{skill}: configured DecontX skip was not recorded correctly: {ambient_rows}")
            if manifest.get("parallel_workers") != 2 or manifest.get("ambient_rna_method") != "skip":
                raise RuntimeError(f"{skill}: execution controls missing from run manifest")
        if skill == "03-scrna-review-qc":
            actual_files = {path.name for path in out.iterdir() if path.is_file()}
            expected_files = set(EXPECTED[skill])
            if actual_files != expected_files:
                raise RuntimeError(f"{skill}: compact output mismatch: expected {sorted(expected_files)}, got {sorted(actual_files)}")
            if manifest.get("output_detail_level") != "compact":
                raise RuntimeError(f"{skill}: compact detail level missing from run manifest")
        if skill == "06-scrna-preprocess-and-cluster":
            with (out / "scenario_summary.tsv").open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            if len(rows) != 1:
                raise RuntimeError(f"{skill}: expected one scenario summary row, got {len(rows)}")
            row = rows[0]
            expected_parameters = {"k_param": "7", "umap_n_neighbors": "9", "umap_min_dist": "0.15", "umap_metric": "euclidean", "umap_method": "uwot", "confirmed_resolution": "0.4", "status": "complete"}
            mismatched = {key: (row.get(key), value) for key, value in expected_parameters.items() if row.get(key) != value}
            if mismatched:
                raise RuntimeError(f"{skill}: resolved/finalized parameters mismatch: {mismatched}")
            if manifest.get("input", {}).get("sha256") is None or manifest.get("config", {}).get("sha256") is None:
                raise RuntimeError(f"{skill}: execution manifest lacks input/config SHA-256")
        if skill == "11-scrna-run-differential-analysis":
            with (out / "task_status.tsv").open(encoding="utf-8", newline="") as handle:
                statuses = [row["status"] for row in csv.DictReader(handle, delimiter="\t")]
            if not statuses or any(status != "completed" for status in statuses):
                raise RuntimeError(f"{skill}: not every E2E task completed: {statuses}")
        if skill == "13-scrna-test-cell-abundance":
            with (out / "task_status.tsv").open(encoding="utf-8", newline="") as handle:
                abundance_status = list(csv.DictReader(handle, delimiter="\t"))
            if {row["method"] for row in abundance_status} != {"propeller", "dcats"} or any(row["status"] != "completed" for row in abundance_status):
                raise RuntimeError(f"{skill}: abundance method execution mismatch: {abundance_status}")
        report["skills"][skill] = {"status": "passed", "artifacts": EXPECTED[skill]}

    report["fixture_sha256_after"] = sha256(fixture)
    if report["fixture_sha256_after"] != fixture_hash:
        raise RuntimeError("source fixture was modified by an executor")
    report["qc_fixture_sha256_before"] = qc_fixture_hash
    report["qc_fixture_sha256_after"] = sha256(qc_fixture)
    if report["qc_fixture_sha256_after"] != qc_fixture_hash:
        raise RuntimeError("multi-layer QC fixture was modified by an executor")
    report["status"] = "passed"
    write_json(output_root / "e2e-report.json", report)
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
