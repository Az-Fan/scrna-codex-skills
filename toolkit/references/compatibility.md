# Environment compatibility

The repository uses existing pixi projects and never creates or repairs environments automatically.

| Skill | Registered pixi project |
|---|---|
| `01-scrna-standardize-input` | `01-scrna-qc` |
| `02-scrna-calculate-qc-metrics` | `01-scrna-qc` |
| `03-scrna-review-qc` | `01-scrna-qc` |
| `04-scrna-apply-qc-filter` | `01-scrna-qc` |
| `05-scrna-benchmark-integration` | `03-integration` |
| `06-scrna-preprocess-and-cluster` | `03-integration` |
| `07-scrna-find-cluster-markers` | `02-annotation` |
| `08-scrna-annotate-cells` | `02-annotation` |
| `09-scrna-export-subset` | `02-annotation` |
| `10-scrna-score-programs` | `05-pathway_program` |
| `11-scrna-run-differential-analysis` | `06-deg-analysis` |
| `12-scrna-run-pathway-enrichment` | `06-deg-analysis` |
| `13-scrna-test-cell-abundance` | `07-cell-abundance` (`default` for R methods; `sccoda` for pertpy/scCODA) |
| `14-scrna-visualize-cell-composition` | `02-annotation` |
| `15-scrna-visualize-gene` | `02-annotation` |

Use the registered probe or dependency checker to resolve the exact interpreter. A missing optional package is reported according to the skill contract; dependencies are never installed implicitly.

Run `scripts/check_dependencies.py --config /absolute/path/config.json` to check selected formats and methods. QS inputs or explicitly requested QS outputs require `qs`; skills 01 and 04 choose RDS when their automatic output format has no `qs` available. The integration check probes the configured Python interpreter and requires `scib_metrics` when metrics are requested. Without a config, checks cover the baseline and report optional branches separately.

The dependency report includes Python distribution versions and rejects the observed `scib-metrics==0.5.7` / `pandas>=3` result-table incompatibility. That scIB version needs `pandas<3`. Configure Conda-to-PyPI mapping for mirrored channels before resolving mixed environments, and reinstall previously overwritten packages when switching their provider. Run an actual benchmark after dependency repair; module presence alone does not prove metrics complete. scIB preparation uses the supplied uncorrected embedding and `benchmark.n_jobs` (default 1).
