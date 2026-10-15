---
name: 18-scrna-prioritize-perturbation
description: Prioritize annotated scRNA-seq cell types by transcriptional separability between two conditions using the native pertpy implementation of Augur, balanced cell subsampling and random-forest cross-validation. Use to rank populations responding to disease, treatment or stimulation; not for relative cell abundance, causal perturbation strength, donor-level significance or pathway-specific SCPA analysis.
---

# Prioritize cell-type responses with Augur

Use the existing registered environment; do not install packages or alter environments during analysis. This skill uses `pertpy.Augur`, a Python implementation linked by the original Augur project, not the original R package. Never silently substitute a custom classifier score for Augur.

## Workflow

1. Confirm the annotated input, biological sample or donor identities, two conditions, cell-type labels, batch field if available, assay and population scope. `metadata.sample` must identify biological samples, not cells. Repeated donors should also have `metadata.donor`. Do not silently subset a population or infer technical batches from sample names.
2. Read [references/method-and-contract.md](references/method-and-contract.md). Adapt [references/config.example.json](references/config.example.json). Accept Seurat RDS/QS or a declared Matrix Market bundle of raw RNA UMI counts. Never use integrated or SCT values as raw counts. Select raw counts (default native-compatible input) or explicit in-memory `log1p_cp10k` normalization; record the choice and inspect depth/QC differences.
3. Run `scripts/check_dependencies.py --config <config>` and dry-run `scripts/run.py --config <config>`. Resolve Python from `07-cell-abundance/sccoda` and Seurat export R from the same registered profile, or use an explicit valid `runtime.augur_python`. Missing interpreters must fail rather than use system fallback.
4. Review cell-type eligibility, independent replicate coverage, sample imbalance and condition/batch confounding. Require at least two biological replicates per condition within each tested cell type. This minimum is an eligibility check, not proof of adequate power. Confounded whole-input designs stop; confounded or insufficient populations remain explicit untestable rows.
5. Execute the reviewed configuration with `--execute`. Native Augur performs variance-based feature selection, balanced subsampling per condition and random-forest cell-level CV. Record every fold, subsample, classifier feature importance and failed or untestable population. Do not silently lower cell-count thresholds or change feature selection after failure.
6. Deliver `RESULTS.md`, the full AUC ranking, subsample variation, sample/QC audits and `interpretation.json`. Inspect all plot pages. Use [references/long-running-execution.md](references/long-running-execution.md) and `scripts/run_in_tmux.py` for long runs.

## Interpretation

AUC measures how well expression distinguishes the two conditions within a cell type. Binary balanced classification has a chance baseline of 0.5; higher scores prioritize more separable populations under the recorded design. It is not an absolute biological perturbation magnitude, activation direction, damage score or causal effect. AUC differences between separate projects, preprocessing choices or contrasts are not automatically comparable.

Native Augur uses pooled **cell-level CV**, not donor-held-out evaluation. Cells from the same sample can appear in training and test folds. Sample/replicate checks and batch audits do not eliminate this limitation. Report rankings as exploratory; never call folds/subsamples biological replicates or turn their variation into donor confidence intervals or P values. Batch audit does not correct batches. Feature importance is not a DEG test or mechanistic evidence. Differential prioritization, permutations, multiclass/regression labels, donor-blocked validation, pathway SCPA and condition-level formal inference are outside this first version.

Use skill 13 for abundance changes, 11 for replicated expression inference, and 10/12 for predefined program scoring or enrichment. Augur answers a different question: which annotated populations are most transcriptionally distinguishable across the selected conditions.

Figure output: one format per run, PDF by default. Read [references/figure-output.md](references/figure-output.md); never export duplicate PDF and PNG copies. Technical input caches, versions, implementation-source checksum and execution records belong under `_provenance/`; keep complete scientific tables accessible.
