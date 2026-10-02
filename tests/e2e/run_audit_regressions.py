#!/usr/bin/env python3
"""Reproduce the 2026-10-02 data-contract failures in an isolated workspace."""
import argparse
import csv
import hashlib
import json
import subprocess
import sys
from pathlib import Path


def read_tsv(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--env-root", type=Path, default=Path.home() / "projects/scrna_envs")
    args = parser.parse_args()
    repo = args.repo.resolve()
    root = args.output_root.resolve()
    root.mkdir(parents=True, exist_ok=False)
    installed = root / "installed"
    subprocess.run([sys.executable, str(repo / "scripts/install_skills.py"), "--target", str(installed)], check=True, stdout=subprocess.DEVNULL)
    rscript = args.env_root / "02-annotation/.pixi/envs/default/bin/Rscript"
    fixture = root / "fixture.rds"
    multi = root / "multilayer.rds"
    subprocess.run([str(rscript), str(repo / "tests/fixtures/create_fixture.R"), str(fixture), str(multi)], check=True, stdout=subprocess.DEVNULL)
    normalized = root / "normalized.rds"
    conflicting = root / "conflicting.rds"
    prepare = root / "prepare.R"
    prepare.write_text('obj <- readRDS(' + json.dumps(str(fixture)) + '); obj <- Seurat::NormalizeData(obj, verbose=FALSE); saveRDS(obj,' + json.dumps(str(normalized)) + '); obj$condition[1] <- "stz"; saveRDS(obj,' + json.dumps(str(conflicting)) + ')\n')
    subprocess.run([str(rscript), str(prepare)], check=True, stdout=subprocess.DEVNULL)
    hashes = {path: hashlib.sha256(path.read_bytes()).hexdigest() for path in [fixture, multi, normalized, conflicting]}

    def run(name, config, label, error=None):
        config = dict(config, output_dir=str(root / label))
        path = root / (label + ".json")
        path.write_text(json.dumps(config, indent=2))
        done = subprocess.run([sys.executable, str(installed / name / "scripts/run.py"), "--config", str(path), "--execute"], text=True, capture_output=True)
        (root / (label + ".log")).write_text(done.stdout + done.stderr)
        if error:
            assert done.returncode and error in done.stdout + done.stderr, (label, done.stdout, done.stderr)
        else:
            assert done.returncode == 0, (label, done.stdout, done.stderr)
        print("PASS", label, flush=True)
        return root / label

    base = {"project": {"id": "audit_regression"}, "runtime": {"pixi_root": str(args.env_root)}}
    metadata = {"sample": "sample_label", "condition": "condition", "batch": "batch_id"}
    standard = dict(base, input={"path": str(fixture), "format": "auto"}, metadata=metadata, output={"object_format": "rds"})
    out = run("01-scrna-standardize-input", standard, "standardized")
    assert len(read_tsv(out / "samples.tsv")) == 4
    subprocess.run([sys.executable, str(installed / "01-scrna-standardize-input/scripts/validate_project.py"), str(out / "samples.tsv")], check=True)
    run("01-scrna-standardize-input", dict(standard, input={"path": str(conflicting), "format": "auto"}), "conflicting_samples", "exactly one condition and batch")
    assert not (root / "conflicting_samples/samples.tsv").exists()

    counts = root / "hierarchy.tsv"
    with counts.open("w") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "condition", "parent", "child", "n_cells"])
        for sample, condition in [("s1", "control"), ("s2", "case")]:
            for parent, child, n in [("vascular", "EC", 80), ("vascular", "Pericyte", 20), ("stromal", "Fibro", 25), ("stromal", "Immune", 75)]:
                writer.writerow([sample, condition, parent, child, n])
    hierarchy = dict(base, input={"counts_table": str(counts), "count_column": "n_cells"}, metadata={"sample": "sample", "condition": "condition"}, composition={"variables": [{"column": "child", "kind": "annotation"}], "parent_column": "parent", "denominator": {"mode": "selected_parent", "include": ["vascular", "stromal"], "description": "Within each parent"}}, modes=["hierarchical_composition"], plots={"umap": False})
    out = run("14-scrna-visualize-cell-composition", hierarchy, "hierarchy")
    assert len(list(out.glob("composition_overview_child_parent_*.png"))) == 2
    assert "Removed" not in (root / "hierarchy.log").read_text()
    assert all(int(row["denominator_cells"]) == 100 for row in read_tsv(out / "sample_coverage.tsv"))
    rows = read_tsv(out / "composition_proportions.tsv")
    for sample in ["s1", "s2"]:
        for parent in ["vascular", "stromal"]:
            assert abs(sum(float(row["proportion"]) for row in rows if row["sample"] == sample and row["parent"] == parent) - 1) < 1e-10

    gene = dict(base, input={"object": str(normalized)}, genes=[{"symbol": "Kdr"}], metadata=dict(metadata, population="cell_type"), expression={"assay": "RNA"}, plots={key: False for key in ["featureplot", "dotplot", "violin", "sample_expression", "sample_heatmap", "de_effect", "volcano_highlight"]})
    de = root / "combined_de.tsv"
    records = [["Kdr", 1, .01, "EC", "case_vs_control"], ["Kdr", -2, .02, "Fibroblast", "case_vs_control"]]
    for reverse in [False, True]:
        with de.open("w") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["gene", "log2FoldChange", "padj", "population", "comparison_id"])
            writer.writerows(records[::-1] if reverse else records)
        config = dict(gene, input=dict(gene["input"], differential_table=str(de)))
        run("15-scrna-visualize-gene", config, f"ambiguous_de_{reverse}", "ambiguous column")
        assert not (root / f"ambiguous_de_{reverse}/target_gene_summary.tsv").exists()
        out = run("15-scrna-visualize-gene", dict(config, differential_selection={"population": "EC", "comparison_id": "case_vs_control"}), f"selected_de_{reverse}")
        row = read_tsv(out / "target_gene_summary.tsv")[0]
        assert float(row["log2_fold_change"]) == 1 and row["population"] == "EC"

    program = dict(base, input={"object": str(multi), "assay": "RNA", "layer": "counts"}, species="mouse",
                   tasks=[{"name": "programs", "method": "addmodulescore", "gene_sets": {"source": "inline", "sets": {"vascular_program": ["Kdr", "Pecam1", "Cdh5"], "stromal_program": ["Col1a1", "Col3a1", "Dcn"]}}, "coverage": {"min_genes": 3, "min_fraction": 1.0, "on_insufficient": "error"}, "parameters": {"nbin": 4, "ctrl": 2}}],
                   summarize_by=["sample_label", "condition", "cell_type"], random_seed=1,
                   cache={"enabled": True}, output={"object_format": "rds"},
                   visualization={"enabled": True, "group_heatmap": {"x": "condition", "facet": "cell_type", "focus_levels": ["Endothelial"]}, "umap": {"enabled": True, "features_per_page": 1}})
    out = run("10-scrna-score-programs", program, "programs")
    task = json.loads((out / "_provenance/task_manifest.json").read_text())["programs"]
    assert task["n_cells"] == 80 and task["cache_hit"] is False
    assert len(list((out / "figures").glob("programs_umap_activity_page*.png"))) == 2
    assert (out / "figures/programs_group_heatmap_focused.png").is_file()
    run("10-scrna-score-programs", program, "programs")
    cached = json.loads((out / "_provenance/task_manifest.json").read_text())["programs"]
    assert cached["cache_hit"] is True and cached["cache_key"] == task["cache_key"]
    run("10-scrna-score-programs", dict(program, random_seed=2), "programs")
    changed = json.loads((out / "_provenance/task_manifest.json").read_text())["programs"]
    assert changed["cache_hit"] is False and changed["cache_key"] != task["cache_key"]
    for path, fingerprint in hashes.items():
        assert hashlib.sha256(path.read_bytes()).hexdigest() == fingerprint
    (root / "regression-report.json").write_text(json.dumps({"status": "passed", "source_fixtures_unchanged": True, "checks": ["sample_mapping", "parent_denominators", "DE_scope_and_order", "program_visualizations_and_cache"]}, indent=2))
    print("PASS all audit regressions")


if __name__ == "__main__":
    main()
