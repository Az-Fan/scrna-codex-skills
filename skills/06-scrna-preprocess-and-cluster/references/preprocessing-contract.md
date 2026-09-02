# Preprocessing contract

## Scenario boundary

Each scenario is an explicitly selected preprocessing choice. Normal execution contains one scenario; comparison mode contains two to four user-selected scenarios. Never add regression or Harmony implicitly.

## Input state

Declare `input.qc_status` as `filtered` or `unfiltered`. An unfiltered run requires `input.allow_unfiltered: true`; its object metadata, workflow state, and manifest remain labelled `exploratory_unfiltered`. Raw counts must exist and be non-negative integers.

## Resolved parameters

Pass `neighbors.k_param` to `FindNeighbors`. Pass `umap.n_neighbors`, `umap.min_dist`, `umap.metric`, and `umap.method` to `RunUMAP`. Record these resolved values in `scenario_summary.tsv`. UMAP settings affect visualization, not graph-cluster membership.

## Resolution selection

For scans, report stability, cluster count, minimum cluster size, maximum cluster fraction, counts of clusters below 50 and 100 cells, and ARI to the previous resolution. `recommended_resolution` is a computational stability candidate, not a biological optimum. Guided mode requires a distinct human-confirmed resolution.

## Finalization

Finalize in the scan output directory. Reuse the stored graph and candidate cluster column, generate the complete final artifact set, write a temporary object, verify that it can be read with unchanged dimensions and cell identities, and atomically replace the scan object. A separate finalize directory is blocked unless `finalize.allow_separate_output: true` is explicitly set.

## Environment

Dependency inspection and execution must resolve the same exact `Rscript` inside the selected existing pixi environment. Do not fall back to a system `Rscript` or create or modify an environment.
