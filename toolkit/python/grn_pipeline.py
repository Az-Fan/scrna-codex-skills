#!/usr/bin/env python3
"""Prepare, infer and score GRN regulons, without downstream interpretation."""
import argparse
import importlib.metadata
import json
import os
from pathlib import Path
import subprocess
import sys

from figure_output import figure_format
from grn_contract import SKILL, path_value, validate_config
from scrna_runtime import nested_get as get, sha256, resolved_rscript


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def value(config, key, default):
    result = get(config, key)
    return default if result is None else result


def binding_config(config):
    bound = {key: config.get(key, {}) for key in ["input", "metadata", "resources", "inference", "metacell", "grn", "ctx", "regulons", "aucell"]}
    bound["figure_format"] = figure_format(config)
    return bound


def file_records(config):
    fields = ["object"] if value(config, "input.type", "seurat") == "seurat" else ["counts", "features", "barcodes", "metadata"]
    paths = {"input." + key: path_value(get(config, "input." + key)) for key in fields}
    paths.update({key: path_value(get(config, key)) for key in ["resources.tf_list", "resources.motif_annotations"]})
    paths.update({f"ranking_database_{i}": path_value(p) for i, p in enumerate(get(config, "resources.ranking_databases"))})
    return {key: {"path": str(path), "sha256": sha256(path)} for key, path in paths.items()}


def run_r(config, stage, technical):
    rscript = resolved_rscript(SKILL, config)
    here = Path(__file__).resolve()
    script = next((p for p in [here.parents[1] / "R/grn_input_and_activity.R", here.with_name("grn_input_and_activity.R")] if p.is_file()), None)
    if not rscript or not script:
        raise ValueError("Registered GRN Rscript or bundled R helper is missing")
    normalized = json.loads(json.dumps(config))
    normalized["output_dir"] = str(path_value(config["output_dir"]))
    for key in ["object", "counts", "features", "barcodes", "metadata"]:
        if get(config, "input." + key):
            normalized["input"][key] = str(path_value(get(config, "input." + key)))
    config_path = technical / "resolved_config.json"
    write_json(config_path, normalized)
    subprocess.run([str(rscript), str(script), str(config_path), stage], check=True)


def make_loom(prepared):
    from scipy.io import mmread
    import loompy
    import numpy as np
    counts = mmread(prepared / "inference_counts.mtx").tocsc()
    genes = (prepared / "inference_genes.tsv").read_text().splitlines()
    ids = (prepared / "inference_ids.tsv").read_text().splitlines()
    if counts.shape != (len(genes), len(ids)):
        raise ValueError("Prepared inference matrix dimensions disagree")
    path = prepared / "inference.loom"
    with loompy.new(str(path)) as ds:
        for offset in range(0, len(ids), 500):
            ds.add_columns(counts[:, offset:offset + 500].toarray(),
                           {"CellID": np.array(ids[offset:offset + 500])},
                           row_attrs={"Gene": np.array(genes)})
    return path, genes, ids


def audit_resources(config, genes, output):
    import pandas as pd
    from ctxcore.rnkdb import FeatherRankingDatabase
    declared = [line.strip() for line in path_value(get(config, "resources.tf_list")).read_text().splitlines() if line.strip()]
    if len(set(declared)) != len(declared):
        raise ValueError("TF list contains duplicate identifiers")
    matched = sorted(set(declared) & set(genes))
    pd.DataFrame({"tf": declared, "in_inference": [x in genes for x in declared]}).to_csv(output / "tf_coverage.tsv", sep="\t", index=False)
    if len(matched) < 2:
        raise ValueError("Fewer than two variable TFs match the inference matrix; check scope and identifiers")
    records = []
    for source in get(config, "resources.ranking_databases"):
        db = FeatherRankingDatabase(str(path_value(source)), name=path_value(source).stem)
        available = set(db.genes)
        overlap = set(genes) & available
        records.append({"database": str(path_value(source)), "database_genes": len(available),
                        "inference_genes": len(genes), "matched_genes": len(overlap), "coverage_fraction": len(overlap) / len(genes)})
        if not overlap:
            raise ValueError("cisTarget database has no overlap with inference genes")
        rank = value(config, "ctx.rank_threshold", 1500)
        if rank >= len(available):
            raise ValueError("ctx.rank_threshold must be smaller than each database's gene count")
    pd.DataFrame(records).to_csv(output / "database_coverage.tsv", sep="\t", index=False)
    return matched


def export_regulons(config, output):
    import pandas as pd
    from pyscenic.utils import load_motifs
    from pyscenic.transform import df2regulons
    motifs = load_motifs(str(output / "motif_enrichment.tsv"), sep="\t")
    sizes = motifs[("Enrichment", "TargetGenes")].map(len)
    audit = motifs.index.to_frame(index=False)
    audit["n_targets"] = sizes.to_numpy()
    audit["status"] = ["retained" if n else "empty_leading_edge" for n in sizes]
    audit.to_csv(output / "motif_row_status.tsv", sep="\t", index=False)
    # ctx 可输出没有 leading-edge 靶基因的行；记录排除理由后再聚合，
    # 保留完整原生 motif 表，不放宽富集阈值。
    regulons = df2regulons(motifs.loc[sizes > 0]) if sizes.gt(0).any() else []
    minimum = value(config, "regulons.min_genes", 20)
    universe = set((output / "_provenance/grn_input/genes.tsv").read_text().splitlines())
    summaries, targets = [], []
    names = set()
    with (output / "regulons.gmt").open("w") as gmt:
        for regulon in regulons:
            if regulon.name in names:
                raise ValueError("pySCENIC returned duplicated regulon names")
            names.add(regulon.name)
            # 保留 upstream 名称和上下文，不把所有网络强行改名为激活态。
            matched = [(gene, weight) for gene, weight in regulon.gene2weight.items() if gene in universe]
            retained = len(matched) >= minimum
            context = ";".join(sorted(map(str, regulon.context)))
            summaries.append({"regulon": regulon.name, "tf": regulon.transcription_factor,
                              "n_targets": len(regulon.genes), "matched_targets": len(matched),
                              "status": "retained" if retained else "below_min_genes", "context": context})
            if retained:
                gmt.write("\t".join([regulon.name, regulon.transcription_factor, *[gene for gene, _ in matched]]) + "\n")
                targets.extend({"regulon": regulon.name, "tf": regulon.transcription_factor, "target": gene, "weight": weight} for gene, weight in matched)
    pd.DataFrame(summaries, columns=["regulon", "tf", "n_targets", "matched_targets", "status", "context"]).to_csv(output / "regulon_summary.tsv", sep="\t", index=False)
    pd.DataFrame(targets, columns=["regulon", "tf", "target", "weight"]).to_csv(output / "regulon_targets.tsv", sep="\t", index=False)
    if not targets:
        raise ValueError("No regulons pass the declared size/coverage threshold; inspect motif_enrichment.tsv and regulon_summary.tsv")
    return sum(row["status"] == "retained" for row in summaries)


def execute(config):
    import pandas as pd
    output = path_value(config["output_dir"])
    technical = output / "_provenance"
    technical.mkdir(parents=True, exist_ok=True)
    prepared = technical / "grn_input"
    record_path = technical / "grn_prepare_record.json"
    action = value(config, "workflow.action", "prepare")
    versions = {p: importlib.metadata.version(p) for p in ["pyscenic", "arboreto", "ctxcore", "numpy", "pandas", "loompy", "scipy"]}
    sources = file_records(config)
    tasks = []
    current = "input_preparation"
    try:
        if action != "prepare" and any((output / name).exists() for name in ["adjacencies.tsv", "motif_enrichment.tsv", "regulon_activity.tsv"]):
            raise ValueError("GRN inference output already exists; preserve it and choose a fresh run")
        if action == "infer":
            record = json.loads(record_path.read_text())
            if record["inputs"] != sources or record["config"] != binding_config(config) or record["versions"] != versions:
                raise ValueError("GRN inputs/resources/parameters/versions differ from the prepared run")
            for relative, digest in record["prepared_files"].items():
                if not (prepared / relative).is_file() or sha256(prepared / relative) != digest:
                    raise ValueError("Prepared GRN artifact changed: " + relative)
        else:
            if prepared.exists() or record_path.exists():
                raise ValueError("GRN prepared input already exists; use infer or a new output_dir")
            run_r(config, "prepare", technical)
            _, genes, ids = make_loom(prepared)
            tfs = audit_resources(config, genes, output)
            write_json(output / "input_audit.json", {"inference_mode": get(config, "inference.mode"), "inference_units": len(ids),
                       "inference_genes": len(genes), "matched_tfs": len(tfs), "counts_source": "raw_umi", "inputs": sources})
            write_json(record_path, {"inputs": sources, "config": binding_config(config), "versions": versions,
                                     "prepared_files": {p.relative_to(prepared).as_posix(): sha256(p) for p in prepared.rglob("*") if p.is_file()}})
        tasks.append({"stage": current, "status": "completed", "reason": ""})
        if action != "prepare":
            workers = str(value(config, "grn.workers", 1))
            expression = str(prepared / "inference.loom")
            current = "grnboost2"
            command = [sys.executable, "-m", "pyscenic.cli.arboreto_with_multiprocessing", expression,
                       str(path_value(get(config, "resources.tf_list"))), "--method", "grnboost2", "--num_workers", workers,
                       "--seed", str(value(config, "grn.seed", 777)), "--output", str(output / "adjacencies.tsv")]
            write_json(technical / "grn_command.json", command)
            subprocess.run(command, check=True)
            adjacency = pd.read_csv(output / "adjacencies.tsv", sep="\t")
            import numpy as np
            if adjacency.empty or not {"TF", "target", "importance"}.issubset(adjacency) or not np.isfinite(adjacency.importance).all() or (adjacency.importance < 0).any():
                raise ValueError("GRNBoost2 produced empty/invalid adjacencies")
            tasks.append({"stage": current, "status": "completed", "reason": ""})
            current = "motif_pruning"
            command = [sys.executable, "-m", "pyscenic.cli.pyscenic", "ctx", str(output / "adjacencies.tsv"),
                       *[str(path_value(p)) for p in get(config, "resources.ranking_databases")],
                       "--annotations_fname", str(path_value(get(config, "resources.motif_annotations"))),
                       "--expression_mtx_fname", expression, "--output", str(output / "motif_enrichment.tsv"),
                       "--num_workers", workers, "--mode", "custom_multiprocessing"]
            for key, default in [("min_genes", 20), ("rank_threshold", 1500), ("auc_threshold", .05), ("nes_threshold", 3.)]:
                command += ["--" + key, str(value(config, "ctx." + key, default))]
            write_json(technical / "ctx_command.json", command)
            subprocess.run(command, check=True)
            n_regulons = export_regulons(config, output)
            tasks.append({"stage": current, "status": "completed", "reason": ""})
            current = "single_cell_aucell"
            run_r(config, "score", technical)
            tasks.append({"stage": current, "status": "completed", "reason": ""})
            write_json(output / "handoff.json", {"regulons": n_regulons, "activity": "regulon_activity.tsv", "targets": "regulon_targets.tsv",
                       "membership": "cell_membership.tsv", "metadata": "cell_metadata.tsv", "downstream_analysis": "not_run"})
        state = {"status": "awaiting_input_confirmation" if action == "prepare" else "complete",
                 "action": action, "review_reason": get(config, "workflow.review_reason"),
                 "next_action": "Review input, inference units and resource coverage; set action=infer" if action == "prepare" else "GRN core complete; select downstream analysis separately"}
    except Exception as exc:
        tasks.append({"stage": current, "status": "failed", "reason": str(exc)})
        state = {"status": "failed", "action": action, "failed_stage": current, "reason": str(exc)}
        raise
    finally:
        if tasks:
            pd.DataFrame(tasks).to_csv(output / "task_status.tsv", sep="\t", index=False)
            write_json(technical / "workflow_state.json", state)
            manifest_path = technical / "run_manifest.json"
            manifest = json.loads(manifest_path.read_text()) if manifest_path.is_file() else {"skill": SKILL}
            manifest.update(inputs=sources, parameters=binding_config(config), python_versions=versions, workflow_action=action,
                            artifacts=[{"path": str(p), "sha256": sha256(p)} for p in output.rglob("*") if p.is_file() and technical not in p.parents and p.name != "RESULTS.md"])
            write_json(manifest_path, manifest)
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", type=Path)
    args = parser.parse_args()
    config = json.loads(args.config.read_text())
    errors = validate_config(config, get)
    if errors:
        raise ValueError("; ".join(errors))
    technical = path_value(config["output_dir"]) / "_provenance"
    technical.mkdir(parents=True, exist_ok=True)
    for name in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"]:
        os.environ[name] = "1"
    os.environ["NUMBA_CACHE_DIR"] = str(technical / "numba_cache")
    os.environ["MPLCONFIGDIR"] = str(technical / "matplotlib")
    return execute(config)


if __name__ == "__main__":
    raise SystemExit(main())
