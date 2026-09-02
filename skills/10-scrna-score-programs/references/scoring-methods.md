# Scoring methods and resources

## Method selection

| Method | Prefer for | Default layer | Important behavior |
|---|---|---|---|
| VISION | Hallmark, metabolism, and moderate collections of curated signatures | `counts` | Reproduces the sc06 convention: scale each cell to median UMI, then calculate signature scores. |
| AUCell | Rank-based activity robust to library-size differences | `counts` | Runtime and memory increase with cells and genes; tune `auc_max_rank` only deliberately. |
| UCell | Rank-based scoring for many cells with stable per-cell interpretation | `counts` | Requires UCell; scores are independent of dataset composition. |
| AddModuleScore | Small, familiar signatures where Seurat-compatible control-bin scoring is desired | `data` | Scores depend on the dataset's expression bins and control-gene sampling. Preserve seed, `nbin`, and `ctrl`. |
| PROGENy | Footprint-based signaling-pathway activity | `data` | Uses the PROGENy model rather than user-provided pathway members. Preserve organism, `top`, and package version. |

When `data` is requested but absent, normalize the selected assay in memory by default and record the resulting session; set `normalize_if_missing` to `false` to require a pre-existing normalized layer. This never changes the input file.

Do not compare numeric values from different methods as if they shared a scale. If the user requests a method comparison, compare within-signature ranks, correlations, group contrasts, and stability rather than raw magnitudes.

Run AUCell's current ranking and AUC calls with one core because its legacy `nCores` path can require an uninstalled `doMC` backend and is deprecated upstream. Treat the skill-level `cores` setting as applicable to methods that safely honor it.

The PROGENy package interface densifies expression internally; the executor does so explicitly for compatibility and defaults to `scale=false`, matching sc06. Confirm adequate RAM before applying PROGENy to a very large full object.

## Gene-set sources

- `inline`: accept a named JSON object whose values are gene-symbol arrays.
- `gmt`: read a local GMT and record its normalized path and MD5.
- `msigdb`: obtain the declared collection through `msigdbr`, cache it as GMT, and record collection and package version. For mouse Hallmark use `MH`; for human use `H`.
- `scmetabolism`: for human data, read the KEGG or REACTOME GMT installed with `scMetabolism`. For mouse, require a reviewed species-matched GMT instead of silently scoring human symbols.

Prefer cached, declared GMT files for strict long-term reproduction. Never silently substitute a collection, release, organism, or identifier namespace.

## Coverage policy

For every signature, record input genes, matched genes, fraction matched, status, and missing genes. `min_genes` protects against scores based on too few genes; `min_fraction` marks suspicious partial matches. By default, skip only `insufficient_genes` signatures and retain `low_fraction` signatures with a warning status.

Investigate systematic low coverage as a likely species, symbol-case, Ensembl/symbol, or gene-set-version mismatch. Do not automatically title-case symbols or perform ortholog conversion because those operations can create false matches; prepare a reviewed GMT when conversion is required.

## cNMF handoff

cNMF discovers programs rather than scoring declared sets. Export non-negative gene-level counts with cell and feature identifiers, test several ranks, retain stability/error diagnostics, and preserve spectra and usage matrices. Interpret programs using top genes, enrichment, sample consistency, and cellular distribution. Keep cNMF outputs separate from curated-score assays.
