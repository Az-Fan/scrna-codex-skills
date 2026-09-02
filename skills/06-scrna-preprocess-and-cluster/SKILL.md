---
name: 06-scrna-preprocess-and-cluster
description: Preprocess and cluster an already filtered Seurat scRNA-seq object using one user-selected workflow or an explicitly requested comparison, with configurable normalization, cell-cycle scoring, covariate regression, fixed-resolution clustering, and optional Harmony. Use after QC filtering and before marker detection or annotation, or when comparing selected preprocessing choices without automatically choosing one.
---

# Preprocess and Cluster scRNA-seq

Use the existing project environment. Never install, create, or update an environment.

## Required conversation

Obtain the Seurat RDS/QS object, assay, sample field, QC status, and explicit output directory. Require `input.qc_status: filtered | unfiltered`; unfiltered input also requires `input.allow_unfiltered: true` and remains explicitly exploratory.

## Mandatory workflow selection gate

Before inspecting parameters deeply or creating any config, use the available interactive user-choice tool. In Wisp, call `ask_user` with exactly these five clickable options:

1. Standard workflow without regression or Harmony.
2. Regress cell-cycle scores.
3. Run Harmony using an explicit batch field.
4. Regress selected covariates and then run Harmony.
5. Compare two or more explicitly selected workflows.

Use the question: `Which preprocessing and clustering workflow should this run use?`

Do not set a timeout or automatic answer. Do not choose on the user's behalf. Do not create a config until the tool returns the user's selection. If no interactive choice tool is available, ask the same five-choice question in text and stop for the answer.

After selection, collect only the required follow-up:

- Option 1: no workflow-specific follow-up.
- Option 2: confirm whether existing cell-cycle score columns should be used or scores should be calculated, including species when calculation is needed.
- Option 3: ask for the batch field and optional Harmony theta.
- Option 4: ask for regression variables, batch field, and optional Harmony theta.
- Option 5: ask which two to four workflows from options 1-4 to compare, then collect the parameters required by those workflows.

Run only the selected workflow by default. Never turn one selected workflow into a multi-scenario comparison. Never infer a batch field from sample names or regression variables from available metadata.

Then ask for the execution mode:

- **Guided**: confirm regression, PCA dimensions, Harmony, resolution scan, and UMAP choices in sequence. Run the scan with `clustering.selection: review`, show the generated resolution UMAP grid, clustree, and stability recommendation, and stop. After the user confirms a resolution, run `workflow.action: finalize_resolution` against the scan object.
- **Batch**: collect all parameters once and run through with fixed resolution or `clustering.selection: recommended` when the user explicitly authorizes automatic use of the stability recommendation.

Do not describe the stability recommendation as a biologically optimal resolution. Keep `recommended_resolution` and `confirmed_resolution` distinct.

In guided mode, collect parameters in dependency order: normalization/HVG; cell-cycle scoring and regression; PCA count and downstream dimensions; optional Harmony fields and theta; neighbor `k.param`; resolution strategy; then UMAP `n.neighbors`, `min.dist`, metric, method, and seed. Explain that all resolution panels share one UMAP and that UMAP parameters do not determine graph clusters.

Before execution, inspect metadata, assays, raw counts, existing reductions, scenario names, and condition-batch confounding. Warn when regression may remove the biology under study or Harmony's batch field is confounded with condition.

## Execution

1. Create a project-local config from [references/config.example.json](references/config.example.json) for batch mode or [references/config.guided.example.json](references/config.guided.example.json) for guided scanning.
2. Run `python3 scripts/check_dependencies.py --pixi-root <registered-root>`, then dry-run `python3 scripts/run.py --config <config>` and verify that the manifest resolves the same exact pixi `Rscript`.
3. Execute only after confirmation with `python3 scripts/run.py --config <config> --execute`.
4. In guided mode, inspect `workflow_state.json`. When it reports `awaiting_resolution_confirmation`, present the resolution UMAP grid, clustree availability, stability table, and recommended value, then wait.
5. After explicit confirmation, create a config from [references/config.finalize.example.json](references/config.finalize.example.json) and finalize the selected candidate column in the same output directory without recomputing preprocessing. The validated replacement object atomically replaces the scan object.
6. When comparison mode was explicitly selected, compare scenarios using diagnostics and composition tables; do not select a winner automatically.

If execution may exceed 10 minutes, has uncertain duration, or could outlive the remote session, read [references/long-running-execution.md](references/long-running-execution.md) and launch the confirmed `--execute` command with `scripts/run_in_tmux.py`. Keep dependency checks, dry runs, and guided-mode confirmation pauses in the foreground. Use a new session name for finalization after a guided scan.

Read [references/preprocessing-contract.md](references/preprocessing-contract.md) when defining scenarios or interpreting outputs.
Read [references/input-output.md](references/input-output.md) when validating inputs or handing off artifacts.

## Guarantees

- Run exactly the selected scenario in normal mode.
- Run multiple scenarios only after the user explicitly selects comparison mode and names the workflows.
- Treat regression and Harmony as explicit choices, never implicit defaults.
- Preserve raw counts, cell identities, and the input object.
- Use separate reductions, graphs, UMAPs, and cluster columns for every scenario.
- Write only one confirmed final cluster column per scenario, even when multiple candidates were scanned.
- Permit a resolution scan, but never silently promote its recommendation in guided mode.
- Do not compute markers, annotate cells, delete clusters, or filter cells.
- Do not treat regression as a substitute for QC filtering or UMAP mixing as proof of successful integration.
