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

Use the registered probe or dependency checker to resolve the exact interpreter. A missing optional package is reported according to the skill contract; dependencies are never installed implicitly.
