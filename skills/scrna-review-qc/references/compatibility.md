# Compatibility

The client computer needs Codex and/or WispScience only. Analysis dependencies are resolved on the registered server where pixi runs.

Tested server baseline: Linux x86_64, Python 3.8.10, R 4.4.3, Seurat 5.4.0, and the project pixi locks. No GPU is required for the nine released executors.

| Skill | Server pixi environment |
|---|---|
| `scrna-calculate-qc-metrics` | `01-scrna-qc` |
| `scrna-review-qc` | `01-scrna-qc` |
| `scrna-standardize-input` | `01-scrna-qc` |
| `scrna-analyze-subset` | `02-annotation` |
| `scrna-annotate-cells` | `02-annotation` |
| `scrna-find-cluster-markers` | `02-annotation` |
| `scrna-benchmark-integration` | `03-integration` |
| `scrna-score-programs` | `05-pathway_program` |
| `scrna-run-differential-analysis` | `06-deg-analysis` |

Run `python scripts/check_dependencies.py --pixi-root ~/projects/scrna_envs` on the server before execution. The checker is read-only. Missing optional format or method packages disable only the corresponding branch.
