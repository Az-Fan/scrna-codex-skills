# Compatibility

The client computer needs Codex and/or WispScience only. Analysis dependencies are resolved on the registered server where pixi runs.

Tested server baseline: Linux x86_64, Python 3.8.10, R 4.4.3, Seurat 5.4.0, and the committed project pixi locks. No GPU is required for the seven core executors.

| Skill | Server pixi environment |
|---|---|
| `scrna-standardize-input` | `01-scrna-qc` |
| `scrna-analyze-subset` | `02-annotation` |
| `scrna-annotate-cells` | `02-annotation` |
| `scrna-find-cluster-markers` | `02-annotation` |
| `scrna-benchmark-integration` | `03-integration` |
| `scrna-score-programs` | `05-pathway_program` |
| `scrna-run-differential-analysis` | `06-deg-analysis` |

Run `python scripts/check_dependencies.py --pixi-root ~/projects/scrna_envs` on the server before execution. It is read-only and never creates environments or installs packages. Optional methods such as scVI, Harmony, VISION, AUCell, UCell, PROGENy, and enrichment backends remain conditional on the selected configuration and pixi feature.
