# Method and input/output contract

## Method

Native `pertpy.tools.Augur` performs per-cell-type variance feature selection and repeated fixed-size subsampling from each condition, followed by stratified cell CV using a random forest. Defaults are 50 subsamples, 20 cells per condition per subsample, 3 folds, feature fraction 0.5, variance quantile 0.5, LOESS span 0.75, 100 trees, seed 14 and one worker. Defaults are choices to review, not universal recommendations. A successful ranking uses the native mean ROC AUC and verifies agreement with exported fold scores.

Each population is passed separately with a fixed internal label, then mapped back using exact original cell-type IDs. This avoids reserved-key and prefix collisions (e.g. `EC` versus `EC_2`) in native result structures. No native cell-wise AUC annotations are interpreted as individual-cell responses. The installed implementation source checksum is recorded in addition to package versions.

Default input is raw RNA counts. Optional `preprocessing.normalization=log1p_cp10k` normalizes within each cell in memory; it is a wrapper preprocessing choice, not a native Augur default. Constant/all-zero feature failures and LOESS errors are reported rather than triggering an automatic alternative algorithm. Scope and preprocessing can change ranks. Inspect sample-level UMI depth and detected-gene differences before interpreting raw-count classification.

## Inputs and design

`input.type=seurat` requires `input.object` and an uncorrected `input.assay` (default RNA). Matrix input requires `counts`, `features`, `barcodes` and `metadata`; orientation defaults to `genes_by_cells`, with explicit `cells_by_genes` supported. Feature and barcode files each contain exactly one nonempty unique ID per line; metadata is a TSV with cell IDs in its first column. Seurat v5 raw layers are joined only in memory by the bundled exporter. Keep input files and the run config outside the fresh output directory so primary exports cannot overwrite their sources. Counts must be finite nonnegative integers with positive totals per cell; source objects and cell IDs are preserved.

Required metadata roles are `sample`, `condition` and `cell_type`. Optional `batch` and `donor` must be declared explicitly. Each sample has one condition and one value of any declared batch/donor field. The two contrast labels must be present and distinct. Other conditions are excluded only by the declared contrast, with an explicit cell-inclusion table. Zero-expression genes in this contrast are recorded, not silently conflated with missing genes.

Each tested type needs at least `augur.subsample_size` cells and `augur.min_replicates_per_condition` independent replicate IDs per condition (default/minimum two). If donor metadata exists, replicate counting uses donors; otherwise it uses samples. All absent/insufficient populations remain in eligibility and task tables. Eligibility does not make cell CV sample-independent. Declared batch/condition design rank is checked globally and within populations. No batch correction, regression, sample balancing or donor-level hypothesis testing is performed. Unequal sample contributions remain visible in `sample_cell_type_qc.tsv`.

## Outputs

- `auc_ranking.tsv`: full completed-population ordering, AUC and explicit cell-CV/exploratory interpretation.
- `fold_auc.tsv`, `subsample_auc.tsv`: all native fold scores and means per subsample. Their variability is not a donor confidence interval.
- `feature_importance.tsv`: complete classifier importance records; no causal or DEG significance interpretation.
- `cell_type_eligibility.tsv`, `task_status.tsv`: every population within the contrast, with reasons for untestable or failed tasks.
- `sample_design.tsv`, `model_design.tsv` (a batch-confounding rank-check matrix only, not fitted classifier covariate adjustment), `sample_cell_type_qc.tsv`, `cell_inclusion.tsv`, `feature_status.tsv`: design, coverage, depth, scope and gene audits.
- `auc_ranking.pdf`, `subsample_auc.pdf`: ranking against the 0.5 baseline and subsampling variation, with bounded pagination when needed. PNG is an alternative format, never an additional copy.
- `interpretation.json` and `interpretation_audit.tsv`: limitations and current status; the audit makes cell-CV limitations visible in `RESULTS.md`. `_provenance/run_manifest.json` binds source inputs/config, options, package versions and implementation-source SHA256. `_provenance/workflow_state.json` distinguishes completed, untestable and failed populations. All failed tasks or no completed population produce a nonzero exit; scientific outputs are preserved.

## Primary sources

- [Original Augur repository and Python implementation link](https://github.com/neurorestore/Augur)
- [Augur original paper](https://www.nature.com/articles/s41587-020-0605-1)
- [Augur protocol](https://www.nature.com/articles/s41596-021-00561-x)
- [pertpy Augur API](https://pertpy.readthedocs.io/en/stable/api/tools/pertpy.tools.Augur.html)

SCPA tests pathway expression distributions across conditions; it does not replace whole-transcriptome cell-type prioritization. It is not implemented by this skill.
