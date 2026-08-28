# Input, output, and usage

## Input

Required:

- A Seurat `.qs` or `.rds` object with non-negative integer raw counts in the chosen assay.
- `input.qc_status = filtered | unfiltered`. Unfiltered input requires `input.allow_unfiltered = true` and is labelled `exploratory_unfiltered`.
- `input.assay` (default `RNA`).
- `metadata.sample` (sample column that exists in the object).
- `scenarios`: one or more explicitly selected scenarios; exactly one is the normal case.
- `output_dir` (never overwritten by upstream steps).

Optional per scenario: `cell_cycle.score`, `regress_variables`, `harmony.enabled` + `harmony.group_by` + `harmony.theta`, `clustering.mode = fixed | scan`, `clustering.stability`, `umap.*`, `neighbors.*`.

## Output

Every run:

- `preprocessed_clustered_object.qs` (counts, data, reductions, graphs, cluster columns; `scale.data` dropped)
- `scenario_summary.tsv`
- `scenario_cluster_similarity.tsv`
- `<scenario>_elbow.png`
- `run_manifest_preprocess.json`, `workflow_state.json`, `session_info.txt`, and append-only `run.log`; finalization adds `run_manifest_finalize.json` without overwriting the preprocessing manifest

Fixed, automatically authorized, or finalized clustering adds:

- `cell_assignments.tsv`
- `<scenario>_cluster_sizes.tsv`
- `<scenario>_sample_cluster_counts.tsv`
- `<scenario>_umap_diagnostics.pdf`; an additional PNG is written only with `plots.preview_png = true`

Resolution scan adds:

- `<scenario>_umap_clusters_by_resolution.png`
- `<scenario>_clustree_resolution.png`
- `<scenario>_resolution_stability.tsv` and `.png`

## Usage

1. Copy `references/config.example.json` (batch, fixed resolution) or `references/config.guided.example.json` (guided, resolution scan with `selection: review`).
2. Dependency check: `python3 scripts/check_dependencies.py --pixi-root <registered-root>`.
3. Dry-run: `python3 scripts/run.py --config <config>.json` and inspect the exact resolved argv and action-aware artifact contract.
4. Execute: `python3 scripts/run.py --config <config>.json --execute`.
5. For a scan with `selection: review`, open the resolution UMAP/clustree/stability outputs, confirm a resolution, then run `references/config.finalize.example.json` in the same output directory. Finalize reuses the already-computed reductions and graphs, validates a temporary replacement object by reading it back, and atomically replaces the scan object.

The skill never auto-selects a resolution. The stability recommendation is recorded separately from the confirmed resolution.
