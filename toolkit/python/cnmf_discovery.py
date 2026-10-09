#!/usr/bin/env python3
"""Run cNMF discovery and reviewed consensus on raw single-cell counts."""
import argparse
import importlib.metadata
import itertools
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

from scrna_runtime import nested_get as get, sha256, resolved_rscript
from cnmf_contract import SKILL, validate_config


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def path_value(value):
    return Path(os.path.expandvars(os.path.expanduser(str(value)))).resolve()


def input_records(config):
    fields = ["object"] if get(config, "input.type") != "matrix" else ["counts", "features", "barcodes", "metadata"]
    return {key: {"path": str(path_value(get(config, "input." + key))),
                  "sha256": sha256(path_value(get(config, "input." + key)))} for key in fields}


def preparation_parameters(config):
    return {key: get(config, "cnmf." + key) if get(config, "cnmf." + key) is not None else default
            for key, default in {"name": "programs", "components": None, "n_iter": 100,
                                 "num_highvar_genes": 2000, "seed": 14, "max_nmf_iter": 1000}.items()}


def input_contract(config):
    return {"input": config.get("input", {}), "metadata": config.get("metadata", {})}


def read_identifiers(path, column=1):
    frame = pd.read_csv(path, sep="\t", header=None, dtype=str, keep_default_na=False)
    if column > frame.shape[1]:
        raise ValueError(f"Identifier column {column} absent in {path}")
    ids = pd.Index(frame.iloc[:, column - 1])
    if ids.has_duplicates or any(not value.strip() for value in ids):
        raise ValueError(f"Empty or duplicate identifiers in {path}; resolve gene identifiers explicitly")
    return ids


def load_input(config, destination):
    kind = get(config, "input.type") or "seurat"
    if kind == "seurat":
        rscript = resolved_rscript(SKILL, config)
        if rscript is None:
            raise ValueError("Registered pathway-program Rscript is unavailable")
        here = Path(__file__).resolve()
        exporter = next((p for p in [here.parents[1] / "R/export_cnmf_input.R", here.parent / "export_cnmf_input.R"] if p.is_file()), None)
        if exporter is None:
            raise ValueError("Bundled export_cnmf_input.R is missing")
        subprocess.run([str(rscript), str(exporter), str(CONFIG_PATH), str(destination)], check=True)
        files = {x: destination / y for x, y in {"counts": "counts.mtx", "features": "features.tsv", "barcodes": "barcodes.tsv", "metadata": "metadata.tsv"}.items()}
        orientation, feature_column = "genes_by_cells", 1
    else:
        files = {x: path_value(get(config, "input." + x)) for x in ["counts", "features", "barcodes", "metadata"]}
        orientation = get(config, "input.orientation") or "genes_by_cells"
        feature_column = get(config, "input.feature_column") or 1
    genes = read_identifiers(files["features"], feature_column)
    cells = read_identifiers(files["barcodes"])
    counts = sparse.csr_matrix(mmread(files["counts"]), dtype=float)
    if orientation == "genes_by_cells":
        counts = counts.T.tocsr()
    counts.sum_duplicates()
    if counts.shape != (len(cells), len(genes)):
        raise ValueError("Matrix orientation/dimensions do not match features and barcodes")
    if not np.isfinite(counts.data).all() or np.any(counts.data < 0) or not np.allclose(counts.data, np.rint(counts.data), rtol=0, atol=1e-8):
        raise ValueError("cNMF input must contain finite, non-negative integer raw UMI counts")
    if np.any(np.asarray(counts.sum(axis=1)).ravel() <= 0):
        raise ValueError("Zero-count cells are present; review QC before program discovery")
    metadata = pd.read_csv(files["metadata"], sep="\t", index_col=0, dtype=str, keep_default_na=False)
    metadata.index = metadata.index.astype(str)
    if metadata.index.has_duplicates or set(metadata.index) != set(cells):
        raise ValueError("Metadata cell IDs must match barcodes exactly without duplicates")
    metadata = metadata.loc[cells]
    metadata.index.name = "cell_id"
    for role in ["sample", "condition", "cell_type", "batch"]:
        column = get(config, "metadata." + role)
        if column and (column not in metadata or metadata[column].str.strip().eq("").any()):
            raise ValueError(f"Missing/empty metadata values: {role} ({column})")
    sample = get(config, "metadata.sample")
    for role in ["condition", "batch"]:
        column = get(config, "metadata." + role)
        if column and metadata.groupby(sample)[column].nunique().gt(1).any():
            raise ValueError(f"Each sample must have one {role} value")
    expressed = np.asarray(counts.sum(axis=0)).ravel() > 0
    pd.DataFrame({"gene": genes, "status": np.where(expressed, "retained", "zero_counts")}).to_csv(destination / "feature_status.tsv", sep="\t", index=False)
    # 显式记录零表达基因；不静默删除细胞或合并重复基因。
    counts, genes = counts[:, expressed], genes[expressed]
    maximum = max(get(config, "cnmf.components"))
    hvg = get(config, "cnmf.num_highvar_genes") or 2000
    if hvg > len(genes):
        raise ValueError(f"num_highvar_genes={hvg} exceeds {len(genes)} expressed genes; choose it explicitly")
    if maximum >= min(len(cells), hvg):
        raise ValueError("Every candidate k must be smaller than both cell count and HVG count")
    genes.name = "gene"
    adata = anndata.AnnData(counts, obs=metadata, var=pd.DataFrame(index=genes))
    adata.write_h5ad(destination / "counts.h5ad", compression="gzip")
    metadata.index.name = "cell_id"
    metadata.to_csv(destination / "cell_metadata.tsv", sep="\t")
    embedding = destination / "embedding.tsv"
    if embedding.is_file():
        emb = pd.read_csv(embedding, sep="\t", index_col=0, dtype={"cell_id": str}, keep_default_na=False)
        if set(emb.index) != set(cells) or emb.index.has_duplicates or not np.isfinite(emb.to_numpy()).all():
            raise ValueError("Embedding must contain finite coordinates for every cell")
    return adata


def summarize_usage(usage, metadata, config):
    columns = list(dict.fromkeys(get(config, "metadata." + role) for role in ["sample", "condition", "cell_type", "batch"] if get(config, "metadata." + role)))
    if set(usage.index) != set(metadata.index) or usage.index.has_duplicates:
        raise ValueError("Consensus usage cell IDs do not match discovery input")
    if set(columns) & set(usage.columns):
        raise ValueError("Metadata column names collide with exported program IDs")
    joined = metadata.loc[usage.index, columns].join(usage)
    grouped = joined.groupby(columns, observed=True, sort=True)
    summary = grouped[list(usage.columns)].mean().reset_index()
    summary["n_cells"] = grouped.size().to_numpy()
    return summary


def plot_summary(summary, programs, config, destination):
    sample = get(config, "metadata.sample")
    population = get(config, "metadata.cell_type")
    labels = summary[sample].astype(str)
    if population:
        labels = labels + " | " + summary[population].astype(str)
    per_page = get(config, "reporting.programs_per_page") or 6
    for offset in range(0, len(programs), per_page):
        selected = programs[offset:offset + per_page]
        fig, ax = plt.subplots(figsize=(max(6, len(selected) * 1.1), max(3, len(labels) * .24)))
        plotted = ax.imshow(summary[selected].to_numpy(), aspect="auto", vmin=0, vmax=1, cmap="viridis")
        ax.set_xticks(range(len(selected)), selected, rotation=45, ha="right")
        ax.set_yticks(range(len(labels)), labels)
        ax.set_title("Mean fractional usage per sample and population")
        fig.colorbar(plotted, ax=ax, label="Fractional usage")
        fig.tight_layout()
        fig.savefig(destination / f"usage_sample_heatmap_page_{offset // per_page + 1}.png", dpi=300)
        plt.close(fig)


def plot_embedding(usage, embedding, config, destination):
    embedding = embedding.loc[usage.index]
    per_page = get(config, "reporting.programs_per_page") or 6
    for offset in range(0, usage.shape[1], per_page):
        selected = usage.columns[offset:offset + per_page]
        fig, axes = plt.subplots((len(selected) + 1) // 2, 2, figsize=(9, 3.7 * ((len(selected) + 1) // 2)), squeeze=False)
        for ax, program in zip(axes.flat, selected):
            points = ax.scatter(embedding.iloc[:, 0], embedding.iloc[:, 1], c=usage[program], s=2, cmap="viridis", vmin=0, vmax=1, rasterized=True)
            ax.set_title(program)
            ax.set_xlabel(embedding.columns[0]); ax.set_ylabel(embedding.columns[1])
            fig.colorbar(points, ax=ax, label="Fractional usage")
        for ax in list(axes.flat)[len(selected):]:
            ax.set_visible(False)
        fig.tight_layout()
        fig.savefig(destination / f"usage_embedding_page_{offset // per_page + 1}.png", dpi=300)
        plt.close(fig)


def enrich_programs(scores, top_n, gmt, destination):
    from scipy.stats import hypergeom
    universe = set(scores.index)
    sets = {}
    with path_value(gmt).open() as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3 or parts[0] in sets:
                raise ValueError("GMT rows need a unique set name, description, and genes")
            sets[parts[0]] = set(parts[2:])
    if not sets:
        raise ValueError("GMT is empty")
    rows, coverage = [], []
    for name, genes in sets.items():
        matched = genes & universe
        coverage.append({"gene_set": name, "input_genes": len(genes), "matched_genes": len(matched), "universe_size": len(universe), "status": "eligible" if matched else "no_overlap"})
    for program in scores:
        selected = set(scores[program][scores[program] > 0].nlargest(top_n).index)
        for name, genes in sets.items():
            matched = genes & universe
            if not matched:
                continue
            overlap = selected & matched
            p = hypergeom.sf(len(overlap) - 1, len(universe), len(matched), len(selected))
            rows.append({"program": program, "gene_set": name, "n_selected": len(selected), "set_size": len(matched), "universe_size": len(universe), "overlap": len(overlap), "overlap_genes": ";".join(sorted(overlap)), "p_value": p})
    table = pd.DataFrame(rows, columns=["program", "gene_set", "n_selected", "set_size", "universe_size", "overlap", "overlap_genes", "p_value"])
    if len(table):
        # 每个 k 的全部 program × gene set 检验共同校正，不先按显著性筛选。
        order = np.argsort(table.p_value.to_numpy())
        adjusted = np.minimum.accumulate((table.p_value.to_numpy()[order] * len(table) / np.arange(1, len(table) + 1))[::-1])[::-1]
        table["p_adjust"] = 1.0
        table.loc[order, "p_adjust"] = np.minimum(adjusted, 1)
    else:
        table["p_adjust"] = pd.Series(dtype=float)
    table.to_csv(destination / "program_enrichment.tsv", sep="\t", index=False)
    pd.DataFrame(coverage).to_csv(destination / "gene_set_coverage.tsv", sep="\t", index=False)
    (destination / "enrichment_universe.txt").write_text("\n".join(sorted(universe)) + "\n")


def export_consensus(model, k, config, metadata, output, embedding):
    threshold = get(config, "cnmf.density_threshold") or .02
    destination = output / f"k_{k}"
    if destination.exists():
        raise ValueError(f"Consensus output already exists: {destination}; use a new output run to preserve it")
    destination.mkdir()
    model.consensus(k=k, density_threshold=threshold,
                    local_neighborhood_size=get(config, "cnmf.local_neighborhood_size") or .3,
                    show_clustering=True, close_clustergram_fig=True)
    top_n = get(config, "reporting.top_n") or 50
    # Native NPZ preserves identifiers such as 001 and NA; load_results reads
    # text with pandas type/NA inference and can alter otherwise valid cell IDs.
    suffix = (k, str(threshold).replace(".", "_"))
    raw_usage = load_df_from_npz(model.paths["consensus_usages"] % suffix)
    scores = load_df_from_npz(model.paths["gene_spectra_score"] % suffix).T
    tpm = load_df_from_npz(model.paths["gene_spectra_tpm"] % suffix).T
    if raw_usage.shape[1] != k or list(scores.columns) != list(raw_usage.columns) or list(tpm.columns) != list(raw_usage.columns):
        raise ValueError("Unexpected cNMF result orientation or mismatched program IDs")
    programs = [f"cnmf{k}_{x}" for x in raw_usage.columns]
    pd.DataFrame({"upstream_program": raw_usage.columns, "program": programs}).to_csv(destination / "program_id_mapping.tsv", sep="\t", index=False)
    raw_usage.columns = scores.columns = tpm.columns = programs
    raw_usage.index = raw_usage.index.astype(str)
    scores.index = scores.index.astype(str)
    tpm.index = tpm.index.astype(str)
    if scores.index.has_duplicates or tpm.index.has_duplicates or not scores.index.equals(tpm.index):
        raise ValueError("Consensus gene spectra IDs disagree or contain duplicates")
    if not np.isfinite(raw_usage.to_numpy()).all() or (raw_usage.to_numpy() < 0).any() or raw_usage.sum(axis=1).le(0).any():
        raise ValueError("Consensus contains non-finite, negative, or zero-total usage")
    if not np.isfinite(scores.to_numpy()).all() or not np.isfinite(tpm.to_numpy()).all() or (tpm.to_numpy() < 0).any():
        raise ValueError("Consensus spectra contain non-finite values or negative TPM weights")
    usage = raw_usage.div(raw_usage.sum(axis=1), axis=0)
    for name, frame in [("usage_raw", raw_usage), ("usage", usage), ("spectra_scores", scores), ("spectra_tpm", tpm)]:
        frame.to_csv(destination / (name + ".tsv"), sep="\t", index_label="cell_id" if "usage" in name else "gene")
    ranked = []
    for program in programs:
        for rank, (gene, score) in enumerate(scores[program].nlargest(top_n).items(), 1):
            ranked.append({"program": program, "rank": rank, "gene": gene, "score": score, "tpm_weight": tpm.loc[gene, program]})
    pd.DataFrame(ranked).to_csv(destination / "top_genes.tsv", sep="\t", index=False)
    with (destination / "programs.gmt").open("w") as handle:
        for program in programs:
            genes = scores[program][scores[program] > 0].nlargest(top_n).index.astype(str)
            handle.write("\t".join([program, "positive_top_score_genes", *genes]) + "\n")
    summary = summarize_usage(usage, metadata, config)
    summary.to_csv(destination / "usage_summary_by_sample.tsv", sep="\t", index=False)
    pd.DataFrame({"program": programs, "candidate_label": "", "evidence": "", "sample_consistency": "", "qc_confounding": "", "decision": "pending"}).to_csv(destination / "program_review.tsv", sep="\t", index=False)
    plot_summary(summary, programs, config, destination)
    if embedding is not None:
        plot_embedding(usage, embedding, config, destination)
    if get(config, "interpretation.gmt"):
        enrich_programs(scores, top_n, get(config, "interpretation.gmt"), destination)
    for figure in Path(model.output_dir, model.name).glob("*clustering*k_%d.*" % k):
        if figure.suffix == ".png":
            shutil.copy2(figure, destination / figure.name)
    return usage


def consensus_registry(output, task_rows):
    """Retain all rank outcomes and compare every completed rank across calls."""
    status_path = output / "task_status.tsv"
    previous = pd.read_csv(status_path, sep="\t", keep_default_na=False).to_dict("records") if status_path.is_file() else []
    rows = previous + task_rows
    if len({row["k"] for row in rows}) != len(rows):
        raise ValueError("Duplicate consensus ranks in task registry")
    pd.DataFrame(rows).to_csv(status_path, sep="\t", index=False)
    usages = {}
    for row in rows:
        if row["status"] == "completed":
            k = int(row["k"])
            usages[k] = pd.read_csv(output / f"k_{k}/usage.tsv", sep="\t", index_col=0,
                                    dtype={"cell_id": str}, keep_default_na=False)
    return rows, usages


def execute(config):
    output = path_value(config["output_dir"])
    technical = output / "_provenance"
    technical.mkdir(parents=True, exist_ok=True)
    prepared = technical / "cnmf_input"
    native = technical / "cnmf"
    binding_path = technical / "discovery_record.json"
    action = get(config, "workflow.action") or "discover"
    params = preparation_parameters(config)
    versions = {x: importlib.metadata.version(x) for x in ["cnmf", "numpy", "pandas", "scipy", "anndata", "scikit-learn"]}
    write_json(technical / "python_session.json", {"executable": sys.executable, "versions": versions})
    sources = input_records(config)
    if get(config, "input.type") == "matrix" and get(config, "metadata.reduction"):
        raise ValueError("Matrix input has no embedding; omit metadata.reduction")
    if action == "consensus":
        binding = json.loads(binding_path.read_text())
        if binding.get("status") != "discovered" or binding["input"] != sources or binding["parameters"] != params or binding["input_contract"] != input_contract(config):
            raise ValueError("Consensus input/parameters differ from the completed discovery; use its original configuration")
        if binding["versions"] != versions:
            raise ValueError("Scientific package versions changed since discovery")
        for name, digest in binding["reusable_files"].items():
            if not (technical / name).is_file() or sha256(technical / name) != digest:
                raise ValueError(f"Discovery artifact changed or is missing: {name}")
        adata = anndata.read_h5ad(prepared / "counts.h5ad")
        model = cNMF(output_dir=str(native), name=params["name"])
    else:
        if prepared.exists() or native.exists() or binding_path.exists():
            raise ValueError("Discovery artifacts already exist; choose a new output_dir or use consensus on the completed run")
        prepared.mkdir(parents=True)
        adata = load_input(config, prepared)
        shutil.copy2(prepared / "feature_status.tsv", output / "feature_status.tsv")
        shutil.copy2(prepared / "cell_metadata.tsv", output / "cell_metadata.tsv")
        write_json(output / "input_audit.json", {"cells": adata.n_obs, "expressed_genes": adata.n_vars, "counts_source": "raw_umi", "matrix_orientation": "cells_by_genes", "input": sources, "parameters": params})
        model = cNMF(output_dir=str(native), name=params["name"])
        model.prepare(counts_fn=str(prepared / "counts.h5ad"), components=params["components"], n_iter=params["n_iter"], seed=params["seed"], num_highvar_genes=params["num_highvar_genes"], max_NMF_iter=params["max_nmf_iter"])
        workers = get(config, "cnmf.workers") or 1
        if workers == 1:
            model.factorize(worker_i=0, total_workers=1)
        else:
            model.factorize_multi_process(total_workers=workers)
        model.combine()
        model.k_selection_plot(close_fig=True)
        stats = load_df_from_npz(model.paths["k_selection_stats"])
        if len(stats) != len(params["components"]) or set(stats["k"]) != set(params["components"]) or not np.isfinite(stats[["silhouette", "prediction_error"]].to_numpy()).all():
            raise ValueError("Incomplete/non-finite k-selection diagnostics")
        stats["k"] = stats["k"].astype(int)
        stats["density_filter_applied"] = False
        stats.to_csv(output / "k_selection_stats.tsv", sep="\t", index=False)
        shutil.copy2(model.paths["k_selection_plot"], output / "k_selection.png")
        reusable = {str(p.relative_to(technical)): sha256(p) for folder in [prepared, native] for p in folder.rglob("*") if p.is_file()}
        write_json(binding_path, {"schema_version": 1, "status": "discovered", "input": sources, "parameters": params, "input_contract": input_contract(config), "versions": versions, "reusable_files": reusable})
    ks = [] if action == "discover" else get(config, "cnmf.consensus_k")
    if ks:
        existing = [k for k in ks if (output / f"k_{k}").exists()]
        if existing:
            raise ValueError(f"Consensus ranks already exported: {existing}; preserve them and choose new ranks")
    embedding = pd.read_csv(prepared / "embedding.tsv", sep="\t", index_col=0, dtype={"cell_id": str}, keep_default_na=False) if (prepared / "embedding.tsv").is_file() else None
    task_rows, usages = [], {}
    for k in ks:
        try:
            usages[k] = export_consensus(model, k, config, adata.obs, output, embedding)
            task_rows.append({"k": k, "status": "completed", "reason": ""})
        except Exception as exc:
            task_rows.append({"k": k, "status": "failed", "reason": str(exc)})
        task_rows[-1].update(density_threshold=get(config, "cnmf.density_threshold") or .02,
                             local_neighborhood_size=get(config, "cnmf.local_neighborhood_size") or .3,
                             selection_reason=get(config, "workflow.selection_reason") or "")
    if ks:
        task_rows, usages = consensus_registry(output, task_rows)
    cross = []
    for first, second in itertools.combinations(usages, 2):
        correlations = usages[first].join(usages[second]).corr(method="spearman")
        for a in usages[first]:
            for b in usages[second]:
                value = correlations.loc[a, b]
                cross.append({"k_a": first, "program_a": a, "k_b": second, "program_b": b, "spearman": value})
    if cross:
        pd.DataFrame(cross).to_csv(output / "cross_k_usage_correlations.tsv", sep="\t", index=False)
    failed = any(x["status"] == "failed" for x in task_rows)
    state = {"status": "failed" if failed else "awaiting_k_confirmation" if action == "discover" else "complete", "components": params["components"], "exported_k": list(usages), "selection_reason": get(config, "workflow.selection_reason"), "next_action": "Review k_selection.png and set action=consensus with consensus_k" if action == "discover" else "Review top genes, sample consistency and program_review.tsv"}
    write_json(technical / "workflow_state.json", state)
    if get(config, "interpretation.gmt"):
        write_json(technical / "gene_set_record.json", {"path": str(path_value(get(config, "interpretation.gmt"))), "sha256": sha256(path_value(get(config, "interpretation.gmt")))})
    manifest_path = technical / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.is_file() else {"skill": SKILL}
    manifest.update(inputs=sources, input=sources.get("object"), parameters=params, python_versions=versions,
                    workflow_action=action, exported_k=list(usages), selection_reason=get(config, "workflow.selection_reason"))
    manifest["artifacts"] = [{"path": str(p), "bytes": p.stat().st_size, "sha256": sha256(p)}
                             for p in output.rglob("*") if p.is_file() and technical not in p.parents and p.name != "RESULTS.md"]
    write_json(manifest_path, manifest)
    return 1 if failed else 0


def main():
    global CONFIG_PATH, np, pd, sparse, mmread, anndata, plt, cNMF, load_df_from_npz
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", type=Path)
    args = parser.parse_args()
    CONFIG_PATH = args.config.resolve()
    config = json.loads(CONFIG_PATH.read_text())
    errors = validate_config(config, get)
    if errors:
        raise ValueError("; ".join(errors))
    # 启动科学库之前设置缓存和线程，避免多个 worker 各自开启大量 BLAS 线程。
    technical = path_value(config["output_dir"]) / "_provenance"
    technical.mkdir(parents=True, exist_ok=True)
    os.environ["MPLCONFIGDIR"] = str(technical / "matplotlib")
    for variable in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"]:
        os.environ[variable] = "1"
    import numpy as np
    import pandas as pd
    from scipy import sparse
    from scipy.io import mmread
    import anndata
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from cnmf import cNMF
    from cnmf.cnmf import load_df_from_npz
    return execute(config)


if __name__ == "__main__":
    raise SystemExit(main())
