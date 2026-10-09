"""Dependency-free configuration contract for cNMF discovery."""
import math
import os
import re
from pathlib import Path

SKILL = "16-scrna-discover-programs"


def validate_config(config, get):
    errors = []
    action = get(config, "workflow.action") or "discover"
    if action not in {"discover", "consensus", "run"}:
        errors.append("workflow.action must be discover, consensus or run")
    kind = get(config, "input.type") or "seurat"
    if kind not in {"seurat", "matrix"}:
        errors.append("input.type must be seurat or matrix")
    fields = ["input.object"] if kind == "seurat" else ["input.counts", "input.features", "input.barcodes", "input.metadata"]
    for field in fields:
        value = get(config, field)
        if not isinstance(value, str) or not value.strip():
            errors.append(f"{field} is required")
        elif not Path(os.path.expandvars(value)).expanduser().is_file():
            errors.append(f"{field} does not exist: {value}")
    if get(config, "input.counts_source") != "raw_umi":
        errors.append("input.counts_source must be raw_umi; corrected or normalized values are unsupported")
    if kind == "seurat" and (get(config, "input.assay") or "RNA").lower() in {"integrated", "sct"}:
        errors.append("cNMF requires an uncorrected RNA raw-count assay")
    if (get(config, "input.orientation") or "genes_by_cells") not in {"genes_by_cells", "cells_by_genes"}:
        errors.append("input.orientation must be genes_by_cells or cells_by_genes")
    if kind == "matrix" and get(config, "metadata.reduction"):
        errors.append("matrix input has no embedding; omit metadata.reduction")
    name = get(config, "cnmf.name") or "programs"
    if not isinstance(name, str) or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", name):
        errors.append("cnmf.name must be a simple alphanumeric run name")
    ranks = get(config, "cnmf.components")
    valid_ranks = isinstance(ranks, list) and bool(ranks) and all(type(x) is int and x >= 2 for x in ranks)
    if not valid_ranks or len(ranks) != len(set(ranks)):
        errors.append("cnmf.components must contain unique integers >= 2")
    if action == "discover" and valid_ranks and len(ranks) < 2:
        errors.append("discover requires at least two candidate ranks")
    for field, default, minimum in [("cnmf.n_iter", 100, 2), ("cnmf.num_highvar_genes", 2000, 2),
                                    ("cnmf.workers", 1, 1), ("cnmf.seed", 14, 0),
                                    ("cnmf.max_nmf_iter", 1000, 1), ("reporting.top_n", 50, 1),
                                    ("input.feature_column", 1, 1), ("reporting.programs_per_page", 6, 1)]:
        value = get(config, field)
        value = default if value is None else value
        if type(value) is not int or value < minimum:
            errors.append(f"{field} must be an integer >= {minimum}")
    for field, default, upper in [("cnmf.density_threshold", 0.02, 2), ("cnmf.local_neighborhood_size", 0.3, 1)]:
        value = get(config, field)
        value = default if value is None else value
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or not 0 < value <= upper:
            errors.append(f"{field} must be finite and in (0, {upper}]")
    ks = get(config, "cnmf.consensus_k")
    if action in {"consensus", "run"}:
        if not isinstance(ks, list) or not ks or any(type(x) is not int for x in ks):
            errors.append("cnmf.consensus_k must explicitly list the ranks to export")
        elif len(set(ks)) != len(ks) or (valid_ranks and not set(ks).issubset(ranks)):
            errors.append("consensus_k must contain unique ranks present in components")
    reason = get(config, "workflow.selection_reason")
    if action == "consensus" and (not isinstance(reason, str) or not reason.strip()):
        errors.append("consensus requires workflow.selection_reason documenting the rank review")
    n_iter = get(config, "cnmf.n_iter")
    n_iter = 100 if n_iter is None else n_iter
    neighborhood = get(config, "cnmf.local_neighborhood_size")
    neighborhood = .3 if neighborhood is None else neighborhood
    if type(n_iter) is int and isinstance(neighborhood, (int, float)) and math.isfinite(neighborhood) and n_iter * neighborhood < 1:
        errors.append("n_iter * local_neighborhood_size must be >= 1 for density estimation")
    prefix = get(config, "runtime.cnmf_python")
    if prefix is not None and (not isinstance(prefix, str) or not prefix.strip()):
        errors.append("runtime.cnmf_python must be an interpreter path")
    for field in ("metadata.sample", "metadata.condition", "metadata.cell_type", "metadata.batch", "metadata.reduction"):
        value = get(config, field)
        if value is not None and (not isinstance(value, str) or not value.strip()):
            errors.append(f"{field} must be a non-empty column/reduction name when provided")
    gmt = get(config, "interpretation.gmt")
    if gmt is not None:
        if not isinstance(gmt, str) or not gmt.strip():
            errors.append("interpretation.gmt must be a non-empty file path when provided")
        elif not Path(os.path.expandvars(gmt)).expanduser().is_file():
            errors.append(f"interpretation.gmt does not exist: {gmt}")
    return errors
