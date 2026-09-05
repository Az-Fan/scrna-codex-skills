# Marker interpretation

## Review checklist

- Confirm that each cluster contains enough cells and occurs across the samples expected for the biological design.
- Prefer coherent marker programs over a single famous gene. Check both positive evidence and exclusion markers.
- Inspect whether top genes are dominated by mitochondrial, ribosomal, heat-shock, immediate-early, cell-cycle, hemoglobin, immunoglobulin, or ambient-RNA signals.
- Review neighboring clusters together. A marker can separate two clusters without being specific to a biological cell type.
- Investigate clusters with no passing markers, very small effect sizes, or nearly identical marker programs before annotation.

## Statistical scope

`FindAllMarkers` compares cells in each identity against cells in the remaining identities. Its adjusted P values do not account for biological replication, and the one-versus-rest contrast changes with cluster composition. Use these results to discover and rank annotation evidence. Do not report them as replicated condition-level differential-expression findings.

## Output contract

- `cluster_markers.tsv`: every marker returned by Seurat plus a within-cluster rank.
- `top_cluster_markers.tsv`: the first `reporting.top_n` markers per cluster using adjusted P value, effect size, and gene as deterministic tie breakers.
- `cluster_marker_summary.tsv`: cell and passing-marker counts for every input cluster, including clusters with zero markers.
- `top_marker_dotplot.png`: expression and detection overview for unique top genes.
- `_provenance/run_manifest.json`: configuration provenance and artifact paths.
