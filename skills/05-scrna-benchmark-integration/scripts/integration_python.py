#!/usr/bin/env python3
"""Run Python integration adapters and scIB benchmarking for the skill exchange bundle."""
import argparse
import itertools
import json
import math
import re
import sys
import traceback
from pathlib import Path

import numpy as np
import pandas as pd


def safe(value):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_")


def expand_methods(config):
    methods = config["benchmark"]["methods"]
    if not any(str(x["name"]).lower() == "none" for x in methods):
        methods = [{"name": "none"}] + methods
    scenarios = []
    for spec in methods:
        name = str(spec["name"]).lower()
        grid = spec.get("parameter_grid", {})
        keys = list(grid)
        rows = itertools.product(*(grid[k] for k in keys)) if keys else [()]
        for values in rows:
            params = dict(zip(keys, values))
            suffix = "__".join("{}_{}".format(k, params[k]) for k in keys)
            scenarios.append({"name": name, "id": safe(name + ("__" + suffix if suffix else "")), "params": params, "spec": spec})
    return scenarios


def read_exchange(exchange):
    import anndata as ad
    from scipy.io import mmread

    meta = pd.read_csv(exchange / "metadata.tsv", sep="\t", index_col="cell_id")
    cells = pd.read_csv(exchange / "barcodes.tsv", header=None)[0].astype(str).tolist()
    genes = pd.read_csv(exchange / "features.tsv", header=None)[0].astype(str).tolist()
    matrix = mmread(str(exchange / "counts.mtx")).tocsr()
    if matrix.shape != (len(cells), len(genes)):
        raise ValueError("counts.mtx shape does not match barcodes/features")
    if set(cells) != set(meta.index.astype(str)):
        raise ValueError("metadata and matrix cell IDs differ")
    meta.index = meta.index.astype(str)
    meta = meta.loc[cells].copy()
    adata = ad.AnnData(X=matrix, obs=meta, var=pd.DataFrame(index=pd.Index(genes, name="gene")))
    adata.var_names_make_unique()
    adata.layers["counts"] = adata.X.copy()
    return adata


def load_r_embeddings(adata, exchange):
    manifest_path = exchange / "embedding_manifest.tsv"
    if not manifest_path.exists() or manifest_path.stat().st_size == 0:
        return []
    manifest = pd.read_csv(manifest_path, sep="\t")
    rows = []
    for row in manifest.to_dict("records"):
        frame = pd.read_csv(row["path"], sep="\t", index_col="cell_id")
        missing = adata.obs_names.difference(frame.index.astype(str))
        if len(missing):
            raise ValueError("{} embedding misses {} cells".format(row["scenario"], len(missing)))
        key = "X__" + row["scenario"]
        adata.obsm[key] = frame.loc[adata.obs_names].to_numpy(dtype=np.float32)
        rows.append({"scenario": row["scenario"], "method": row["method"], "parameters": row.get("parameters", "{}"), "representation": key, "representation_type": row.get("representation_type", "embedding"), "status": "completed", "notes": "R adapter"})
    return rows


def merge_r_status(output, completed_embeddings):
    status_path = output / "method_runs_r.tsv"
    if not status_path.exists():
        return completed_embeddings
    completed = {row["scenario"]: row for row in completed_embeddings}
    merged = []
    for row in pd.read_csv(status_path, sep="\t", keep_default_na=False).to_dict("records"):
        if row["status"] == "completed" and row["scenario"] in completed:
            merged.append(completed[row["scenario"]])
        else:
            merged.append({"scenario": row["scenario"], "method": row["method"], "parameters": row.get("parameters", "{}"), "representation": row.get("representation", ""), "representation_type": "embedding", "status": row["status"], "notes": row.get("notes", "")})
    return merged


def train_python_methods(adata, scenarios, batch_key, label_key, seed):
    import scanpy as sc

    rows = []
    for scenario in [x for x in scenarios if x["name"] in {"scvi", "scanvi", "bbknn"}]:
        rep = "X__" + scenario["id"]
        status, notes, rep_type = "completed", "", "embedding"
        try:
            if scenario["name"] in {"scvi", "scanvi"}:
                import scvi
                scvi.settings.seed = seed
                work = adata.copy()
                scvi.model.SCVI.setup_anndata(work, layer="counts", batch_key=batch_key)
                model = scvi.model.SCVI(work, n_latent=int(scenario["params"].get("n_latent", 30)),
                    n_layers=int(scenario["params"].get("n_layers", 2)), gene_likelihood=scenario["params"].get("gene_likelihood", "nb"))
                model.train(max_epochs=int(scenario["params"].get("max_epochs", 200)), accelerator="cpu")
                if scenario["name"] == "scanvi":
                    if not label_key:
                        raise ValueError("scANVI requires metadata.biological_labels")
                    unlabeled = str(scenario["spec"].get("unlabeled_category", "Unknown"))
                    work.obs[label_key] = work.obs[label_key].astype(str).fillna(unlabeled).astype("category")
                    if unlabeled not in work.obs[label_key].cat.categories:
                        work.obs[label_key] = work.obs[label_key].cat.add_categories([unlabeled])
                    model = scvi.model.SCANVI.from_scvi_model(model, labels_key=label_key, unlabeled_category=unlabeled)
                    model.train(max_epochs=int(scenario["params"].get("scanvi_max_epochs", 100)), accelerator="cpu")
                adata.obsm[rep] = model.get_latent_representation().astype(np.float32)
            else:
                import bbknn
                work = adata.copy()
                sc.pp.normalize_total(work, target_sum=1e4); sc.pp.log1p(work)
                sc.pp.highly_variable_genes(work, n_top_genes=int(scenario["params"].get("n_top_genes", 2000)), subset=True)
                sc.pp.scale(work, max_value=10); sc.tl.pca(work, n_comps=int(scenario["params"].get("n_pcs", 30)), random_state=seed)
                bbknn.bbknn(work, batch_key=batch_key, neighbors_within_batch=int(scenario["params"].get("neighbors_within_batch", 3)), n_pcs=int(scenario["params"].get("n_pcs", 30)))
                sc.tl.umap(work, random_state=seed)
                adata.obsm[rep] = work.obsm["X_umap"].astype(np.float32)
                rep_type = "graph_umap"
        except Exception as exc:
            status = "missing_dependency" if isinstance(exc, (ImportError, ModuleNotFoundError)) else "failed"
            notes = "{}: {}".format(type(exc).__name__, exc)
        rows.append({"scenario": scenario["id"], "method": scenario["name"], "parameters": json.dumps(scenario["params"], sort_keys=True), "representation": rep, "representation_type": rep_type, "status": status, "notes": notes})
    return rows


def confounding_table(obs, batch_keys, condition_key):
    rows = []
    for batch in batch_keys:
        tab = pd.crosstab(obs[batch], obs[condition_key]) if condition_key else pd.DataFrame()
        if tab.empty:
            v, perfect = np.nan, False
        else:
            from scipy.stats import chi2_contingency
            chi2 = chi2_contingency(tab, correction=False)[0]
            n = tab.to_numpy().sum(); denom = max(1, min(tab.shape) - 1)
            v = math.sqrt(chi2 / (n * denom)) if n else np.nan
            perfect = bool(v > 0.999999)
        rows.append({"batch_variable": batch, "condition": condition_key or "", "cramers_v": v, "perfect_confounding": perfect, "n_batch_levels": obs[batch].nunique(), "n_condition_levels": obs[condition_key].nunique() if condition_key else np.nan})
    return pd.DataFrame(rows)


METRIC_MATCH = {
    "ilisi": ("ilisi",), "batch_asw": ("bras",), "pcr_comparison": ("pcr",),
    "graph_connectivity": ("graph", "connect"), "kbet": ("kbet",), "clisi": ("clisi",),
    "label_asw": ("silhouette", "label"), "isolated_labels": ("isolated",), "nmi": ("nmi",), "ari": ("ari",),
}
METRIC_COLUMNS = ["scenario", "batch_variable", "biological_label", "metric", "metric_group", "value", "status", "notes"]


def normalize_metric(value):
    return re.sub(r"[^a-z0-9]+", " ", str(value).lower()).strip()


def benchmark_metrics(adata, runs, config, output):
    requested_batch = config["metrics"].get("batch_removal", [])
    requested_bio = config["metrics"].get("biological_conservation", [])
    if not requested_batch and not requested_bio:
        return pd.DataFrame(columns=METRIC_COLUMNS)
    try:
        from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection
    except Exception as exc:
        rows = []
        for run in runs:
            if run["status"] == "completed":
                for batch in config["metadata"]["batch_variables"]:
                    for metric in requested_batch + requested_bio:
                        rows.append({"scenario": run["scenario"], "batch_variable": batch, "biological_label": "", "metric": metric, "metric_group": "batch_removal" if metric in requested_batch else "biological_conservation", "value": np.nan, "status": "missing_dependency", "notes": str(exc)})
        return pd.DataFrame(rows, columns=METRIC_COLUMNS)
    embedding_runs = [x for x in runs if x["status"] == "completed" and x.get("representation_type") == "embedding"]
    audit_rows = []
    for run in runs:
        if run in embedding_runs:
            continue
        status = "skipped_incompatible_representation" if run["status"] == "completed" else run["status"]
        note = "graph/UMAP representation is not scored as a latent embedding" if run["status"] == "completed" else run.get("notes", "method did not complete")
        for batch in config["metadata"]["batch_variables"]:
            for label in (config["metadata"].get("biological_labels") or [""]):
                for metric in requested_batch + requested_bio:
                    audit_rows.append({"scenario": run["scenario"], "batch_variable": batch, "biological_label": label, "metric": metric, "metric_group": "batch_removal" if metric in requested_batch else "biological_conservation", "value": np.nan, "status": status, "notes": note})
    keys = [x["representation"] for x in embedding_runs]
    lookup = {x["representation"]: x["scenario"] for x in embedding_runs}
    rows = list(audit_rows)
    for batch in config["metadata"]["batch_variables"]:
        for label in config["metadata"].get("biological_labels", []):
            try:
                bio_options = BioConservation(
                    isolated_labels="isolated_labels" in requested_bio,
                    nmi_ari_cluster_labels_leiden=False,
                    nmi_ari_cluster_labels_kmeans=bool({"nmi", "ari"} & set(requested_bio)),
                    silhouette_label="label_asw" in requested_bio,
                    clisi_knn="clisi" in requested_bio,
                )
                batch_options = BatchCorrection(
                    bras="batch_asw" in requested_batch,
                    ilisi_knn="ilisi" in requested_batch,
                    kbet_per_label="kbet" in requested_batch,
                    graph_connectivity="graph_connectivity" in requested_batch,
                    pcr_comparison="pcr_comparison" in requested_batch,
                )
                bench = Benchmarker(adata, batch_key=batch, label_key=label, embedding_obsm_keys=keys,
                    bio_conservation_metrics=bio_options, batch_correction_metrics=batch_options, n_jobs=-1)
                bench.benchmark()
                result = bench.get_results(min_max_scale=True)
                if "Embedding" in result.columns:
                    result = result.set_index("Embedding")
                if any("batch correction" in normalize_metric(x) for x in result.index):
                    result = result.T
                result = result.loc[result.index.intersection(keys)]
                for key, values in result.iterrows():
                    normalized = {column: normalize_metric(column) for column in result.columns}
                    for metric in requested_batch + requested_bio:
                        tokens = METRIC_MATCH[metric]
                        matches = [col for col, text in normalized.items() if all(token in text for token in tokens)]
                        if matches:
                            value = pd.to_numeric(pd.Series([values[matches[0]]]), errors="coerce").iloc[0]
                            status, note = "completed", "source_metric=" + str(matches[0])
                        else:
                            value, status, note = np.nan, "skipped_incompatible_representation", "metric absent from scIB result"
                        rows.append({"scenario": lookup[key], "batch_variable": batch, "biological_label": label, "metric": metric, "metric_group": "batch_removal" if metric in requested_batch else "biological_conservation", "value": value, "status": status, "notes": note})
            except Exception as exc:
                for key in keys:
                    for metric in requested_batch + requested_bio:
                        rows.append({"scenario": lookup[key], "batch_variable": batch, "biological_label": label, "metric": metric, "metric_group": "batch_removal" if metric in requested_batch else "biological_conservation", "value": np.nan, "status": "failed", "notes": "{}: {}".format(type(exc).__name__, exc)})
    return pd.DataFrame(rows, columns=METRIC_COLUMNS)


def summarize(metrics, config):
    complete = metrics[metrics["status"] == "completed"].copy()
    if complete.empty:
        return pd.DataFrame(), pd.DataFrame()
    summary = complete.groupby(["scenario", "metric_group"], observed=True)["value"].mean().unstack()
    summary = summary.rename(columns={"batch_removal": "batch_score", "biological_conservation": "biology_score"}).reset_index()
    for column in ("batch_score", "biology_score"):
        if column not in summary: summary[column] = np.nan
    scoring = config.get("scoring", {})
    summary["weighted_score"] = np.nan
    if scoring.get("enabled"):
        both = summary[["batch_score", "biology_score"]].notna().all(axis=1)
        summary.loc[both, "weighted_score"] = scoring.get("batch_weight", .3) * summary.loc[both, "batch_score"] + scoring.get("biology_weight", .7) * summary.loc[both, "biology_score"]
    summary["pareto_efficient"] = False
    valid = summary.dropna(subset=["batch_score", "biology_score"])
    for idx, row in valid.iterrows():
        dominated = ((valid["batch_score"] >= row["batch_score"]) & (valid["biology_score"] >= row["biology_score"]) & ((valid["batch_score"] > row["batch_score"]) | (valid["biology_score"] > row["biology_score"]))).any()
        summary.loc[idx, "pareto_efficient"] = not dominated
    ranking = summary.sort_values(["pareto_efficient", "weighted_score", "biology_score", "batch_score"], ascending=False, na_position="last")
    return summary, ranking


def recommendation_status(metrics, ranking, confounding, config):
    """Return a machine-readable decision state; never infer a winner from partial evidence."""
    requested = list(config["metrics"].get("batch_removal", [])) + list(
        config["metrics"].get("biological_conservation", [])
    )
    reasons = []
    if requested and metrics[metrics["status"] == "completed"].empty:
        reasons.append("no_requested_metric_completed")
    if ranking.empty:
        reasons.append("no_comparable_ranking")
    if not confounding.empty and confounding["perfect_confounding"].fillna(False).any():
        reasons.append("batch_condition_perfectly_confounded")
    state = "resolved" if not reasons else "unresolved"
    recommended = None
    if state == "resolved" and not ranking.empty:
        candidates = ranking[ranking["pareto_efficient"]]
        if len(candidates) == 1:
            recommended = str(candidates.iloc[0]["scenario"])
        else:
            state = "review_required"
            reasons.append("multiple_pareto_efficient_scenarios")
    return {
        "status": state,
        "recommended_scenario": recommended,
        "reasons": reasons,
        "requested_metric_count": len(requested),
        "completed_metric_rows": int((metrics["status"] == "completed").sum()) if not metrics.empty else 0,
    }


def make_plots(adata, runs, metrics, summary, config, output, seed):
    import matplotlib.pyplot as plt
    import scanpy as sc
    plots = set(config.get("plots", [])); figure_dir = output / "figures"; figure_dir.mkdir(exist_ok=True)
    colors = {"umap_by_batch": config["metadata"]["batch_variables"], "umap_by_sample": [config["metadata"]["sample"]], "umap_by_condition": [config["metadata"].get("condition")], "umap_by_label": config["metadata"].get("biological_labels", [])}
    for run in [x for x in runs if x["status"] == "completed"]:
        rep = run["representation"]; umap_key = "X_umap__" + run["scenario"]
        if run.get("representation_type") == "graph_umap":
            adata.obsm[umap_key] = adata.obsm[rep]
        else:
            sc.pp.neighbors(adata, use_rep=rep, key_added="neighbors__" + run["scenario"])
            sc.tl.umap(adata, neighbors_key="neighbors__" + run["scenario"], random_state=seed)
            adata.obsm[umap_key] = adata.obsm["X_umap"].copy()
        adata.obsm["X_plot"] = adata.obsm[umap_key]
        for plot, fields in colors.items():
            if plot not in plots: continue
            fields = [x for x in fields if x and x in adata.obs]
            if fields:
                fig = sc.pl.embedding(adata, basis="plot", color=fields, show=False, return_fig=True, title=[run["scenario"] + " | " + x for x in fields])
                fig.savefig(figure_dir / (plot + "__" + run["scenario"] + ".png"), dpi=200, bbox_inches="tight"); plt.close(fig)
    if not summary.empty and "metric_tradeoff" in plots:
        fig, ax = plt.subplots(figsize=(7, 6)); ax.scatter(summary["batch_score"], summary["biology_score"])
        for _, row in summary.iterrows(): ax.annotate(row["scenario"], (row["batch_score"], row["biology_score"]), fontsize=7)
        ax.set(xlabel="Batch removal", ylabel="Biological conservation"); fig.tight_layout(); fig.savefig(figure_dir / "metric_tradeoff.png", dpi=250); plt.close(fig)
    if not metrics.empty and "score_heatmap" in plots:
        table = metrics[metrics.status == "completed"].pivot_table(index="scenario", columns="metric", values="value", aggfunc="mean")
        if not table.empty:
            import seaborn as sns
            fig, ax = plt.subplots(figsize=(max(7, .8 * table.shape[1]), max(4, .4 * table.shape[0]))); sns.heatmap(table, annot=True, fmt=".2f", cmap="viridis", ax=ax); fig.tight_layout(); fig.savefig(figure_dir / "score_heatmap.png", dpi=250); plt.close(fig)
    if not summary.empty and "score_barplot" in plots:
        long = summary.melt(id_vars="scenario", value_vars=["batch_score", "biology_score"], var_name="score_type", value_name="score")
        import seaborn as sns
        fig, ax = plt.subplots(figsize=(max(8, .8 * len(summary)), 5)); sns.barplot(data=long, x="scenario", y="score", hue="score_type", ax=ax); ax.tick_params(axis="x", rotation=45); fig.tight_layout(); fig.savefig(figure_dir / "score_barplot.png", dpi=250); plt.close(fig)
    if not summary.empty and "ranking_plot" in plots:
        ranked = summary.sort_values("weighted_score" if summary.weighted_score.notna().any() else "biology_score")
        value = "weighted_score" if ranked.weighted_score.notna().any() else "biology_score"
        fig, ax = plt.subplots(figsize=(8, max(4, .35 * len(ranked)))); ax.barh(ranked.scenario, ranked[value]); ax.set_xlabel(value); fig.tight_layout(); fig.savefig(figure_dir / "ranking_plot.png", dpi=250); plt.close(fig)
    programs = config.get("gene_programs", {})
    if programs and ({"marker_dotplot", "program_retention"} & plots):
        expr = adata.copy(); expr.X = expr.layers["counts"].copy(); sc.pp.normalize_total(expr, target_sum=1e4); sc.pp.log1p(expr)
        genes = sorted({gene for values in programs.values() for gene in values if gene in expr.var_names})
        retention_rows = []
        for run in [x for x in runs if x["status"] == "completed"]:
            neighbor_key = "neighbors__" + run["scenario"]
            if run.get("representation_type") == "graph_umap":
                continue
            if neighbor_key not in adata.uns:
                sc.pp.neighbors(adata, use_rep=run["representation"], key_added=neighbor_key)
            cluster_key = "clusters__" + run["scenario"]
            sc.tl.leiden(adata, neighbors_key=neighbor_key, key_added=cluster_key, resolution=1.0)
            expr.obs[cluster_key] = adata.obs[cluster_key].astype(str).values
            if genes and "marker_dotplot" in plots:
                dot = sc.pl.dotplot(expr, genes, groupby=cluster_key, show=False, return_fig=True)
                dot.savefig(figure_dir / ("marker_dotplot__" + run["scenario"] + ".png"))
            for program, members in programs.items():
                found = [x for x in members if x in expr.var_names]
                if not found: continue
                score_key = "program__" + safe(program)
                sc.tl.score_genes(expr, found, score_name=score_key)
                grouped = expr.obs.groupby(cluster_key, observed=True)[score_key].mean()
                retention_rows.append({"scenario": run["scenario"], "program": program, "max_cluster_score": grouped.max(), "n_genes_found": len(found), "n_genes_requested": len(members)})
        if retention_rows:
            retention = pd.DataFrame(retention_rows)
            baseline = retention[retention.scenario == "none"].set_index("program")["max_cluster_score"]
            retention["retention_vs_none"] = retention.apply(lambda row: row.max_cluster_score / baseline.get(row.program, np.nan) if baseline.get(row.program, 0) != 0 else np.nan, axis=1)
            retention.to_csv(output / "program_retention.tsv", sep="\t", index=False)
            if "program_retention" in plots:
                import seaborn as sns
                table = retention.pivot(index="scenario", columns="program", values="retention_vs_none")
                fig, ax = plt.subplots(figsize=(max(7, .8 * table.shape[1]), max(4, .4 * table.shape[0]))); sns.heatmap(table, annot=True, fmt=".2f", center=1, cmap="vlag", ax=ax); fig.tight_layout(); fig.savefig(figure_dir / "program_retention.png", dpi=250); plt.close(fig)


def main():
    parser = argparse.ArgumentParser(); parser.add_argument("--config", required=True); parser.add_argument("--exchange", required=True)
    args = parser.parse_args(); config_path = Path(args.config); exchange = Path(args.exchange)
    config = json.loads(config_path.read_text(encoding="utf-8")); output = Path(config["output_dir"]).resolve(); output.mkdir(parents=True, exist_ok=True)
    seed = int(config.get("benchmark", {}).get("seed", 1234)); np.random.seed(seed)
    adata = read_exchange(exchange); scenarios = expand_methods(config)
    runs = merge_r_status(output, load_r_embeddings(adata, exchange))
    primary_batch = config["metadata"]["batch_variables"][0]
    labels = config["metadata"].get("biological_labels", []); primary_label = labels[0] if labels else None
    runs.extend(train_python_methods(adata, scenarios, primary_batch, primary_label, seed))
    pd.DataFrame(runs).to_csv(output / "method_runs.tsv", sep="\t", index=False)
    confounding = confounding_table(adata.obs, config["metadata"]["batch_variables"], config["metadata"].get("condition"))
    confounding.to_csv(output / "design_confounding.tsv", sep="\t", index=False)
    metrics = benchmark_metrics(adata, runs, config, output)
    metrics.to_csv(output / "metric_results_long.tsv", sep="\t", index=False)
    metrics[metrics.status != "completed"].to_csv(output / "skipped_metrics.tsv", sep="\t", index=False)
    summary, ranking = summarize(metrics, config); summary.to_csv(output / "method_summary.tsv", sep="\t", index=False); ranking.to_csv(output / "method_ranking.tsv", sep="\t", index=False)
    decision = recommendation_status(metrics, ranking, confounding, config)
    (output / "recommendation_status.json").write_text(json.dumps(decision, indent=2) + "\n", encoding="utf-8")
    make_plots(adata, runs, metrics, summary, config, output, seed)
    adata.write_h5ad(output / "benchmark_embeddings.h5ad", compression="gzip")
    lines = ["# Integration benchmark recommendation", "", "Do not select a method from UMAP alone.", ""]
    lines += ["Decision status: **{}**".format(decision["status"]), ""]
    if decision["reasons"]:
        lines += ["Reasons: " + ", ".join(decision["reasons"]), ""]
    if decision["status"] == "resolved" and not ranking.empty:
        try:
            ranking_text = ranking.to_markdown(index=False)
        except ImportError:
            ranking_text = "```text\n" + ranking.to_csv(sep="\t", index=False) + "```"
        lines += ["Pareto-efficient scenarios: " + ", ".join(ranking.loc[ranking.pareto_efficient, "scenario"].astype(str)), "", ranking_text]
    else: lines.append("No integration method is recommended. Resolve the recorded dependency, metric, or design limitation before selecting a correction scenario.")
    (output / "recommendation.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    try: main()
    except Exception:
        traceback.print_exc(); sys.exit(1)
