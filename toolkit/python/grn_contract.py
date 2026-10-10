"""Configuration checks for the core pySCENIC workflow."""
from figure_output import figure_format
import math
import os
from pathlib import Path

SKILL = "17-scrna-infer-grn"


def path_value(value):
    return Path(os.path.expandvars(os.path.expanduser(str(value)))).resolve()


def validate_config(config, get):
    errors = []
    try:
        figure_format(config)
    except ValueError as exc:
        errors.append(str(exc))
    action = get(config, "workflow.action") or "prepare"
    if action not in {"prepare", "infer", "run"}:
        errors.append("workflow.action must be prepare, infer or run")
    kind = get(config, "input.type") or "seurat"
    if kind not in {"seurat", "matrix"}:
        errors.append("input.type must be seurat or matrix")
    fields = ["input.object"] if kind == "seurat" else ["input.counts", "input.features", "input.barcodes", "input.metadata"]
    fields += ["resources.tf_list", "resources.motif_annotations"]
    for key in fields:
        value = get(config, key)
        if not isinstance(value, str) or not value.strip() or not path_value(value).is_file():
            errors.append(f"{key} must name an existing file")
    dbs = get(config, "resources.ranking_databases")
    if not isinstance(dbs, list) or not dbs or any(not isinstance(p, str) or not path_value(p).is_file() for p in dbs):
        errors.append("resources.ranking_databases must list existing cisTarget Feather files")
    elif len(set(map(str, map(path_value, dbs)))) != len(dbs):
        errors.append("resources.ranking_databases must not contain duplicate paths")
    output_value = get(config, "output_dir")
    if not isinstance(output_value, str) or not output_value.strip():
        errors.append("output_dir must name a separate output directory")
    else:
        output = path_value(output_value)
        sources = [get(config, key) for key in fields] + (dbs if isinstance(dbs, list) else [])
        if any(isinstance(source, str) and output in path_value(source).parents for source in sources):
            errors.append("Inputs/resources must be outside output_dir to preserve their source files")
    if get(config, "input.counts_source") != "raw_umi":
        errors.append("input.counts_source must be raw_umi")
    for role in ["species", "genome", "gene_identifier"]:
        value = get(config, "input." + role)
        if not isinstance(value, str) or not value.strip() or value != get(config, "resources." + role):
            errors.append(f"input.{role} and resources.{role} must explicitly match")
    assay = get(config, "input.assay")
    if assay is not None and (not isinstance(assay, str) or not assay.strip() or assay.lower() in {"sct", "integrated"}):
        errors.append("input.assay must name an uncorrected RNA assay")
    if (get(config, "input.orientation") or "genes_by_cells") not in {"genes_by_cells", "cells_by_genes"}:
        errors.append("input.orientation must be genes_by_cells or cells_by_genes")
    for key in ["sample", "cell_type", "condition", "batch"]:
        value = get(config, "metadata." + key)
        if (key == "sample" or value is not None) and (not isinstance(value, str) or not value.strip()):
            errors.append(f"metadata.{key} must name a non-empty metadata column")
    mode = get(config, "inference.mode")
    if mode not in {"single_cell", "metacell"}:
        errors.append("inference.mode must explicitly be single_cell or metacell")
    if mode == "metacell":
        column = get(config, "metacell.column")
        if not isinstance(column, str) or not column.strip():
            errors.append("metacell.column must name reviewed cell memberships; use skill 06 to generate clusters if needed")
    if action == "infer" and not isinstance(get(config, "workflow.review_reason"), str):
        errors.append("infer requires workflow.review_reason documenting input/metacell review")
    elif action == "infer" and not get(config, "workflow.review_reason").strip():
        errors.append("workflow.review_reason must be non-empty")
    if (get(config, "grn.method") or "grnboost2") != "grnboost2":
        errors.append("This skill implements GRNBoost2; grn.method must be grnboost2")
    for key, default, minimum in [("grn.workers", 1, 1), ("grn.seed", 777, 0), ("ctx.min_genes", 20, 2),
                                  ("ctx.rank_threshold", 1500, 1), ("regulons.min_genes", 20, 2),
                                  ("aucell.batch_size", 500, 1), ("input.feature_column", 1, 1),
                                  ("metacell.min_cells", 1, 1)]:
        value = get(config, key)
        if value is None:
            value = default
        if type(value) is not int or value < minimum:
            errors.append(f"{key} must be an integer >= {minimum}")
    for key, default, upper in [("ctx.auc_threshold", .05, 1), ("aucell.rank_fraction", .05, 1), ("ctx.nes_threshold", 3., None)]:
        value = get(config, key)
        value = default if value is None else value
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0 or (upper and not 0 < value <= upper):
            errors.append(f"{key} must be finite and non-negative" + (f" in (0,{upper}]" if upper else ""))
    override = get(config, "runtime.pyscenic_python")
    if override is not None and (not isinstance(override, str) or not override.strip()):
        errors.append("runtime.pyscenic_python must name an interpreter")
    return errors
