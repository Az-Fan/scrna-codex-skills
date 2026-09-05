# Fixed paper style v1

Implemented in skills 06, 07 and 08. Other families in figure-catalog.md remain inventoried, pending rollout; do not claim the entire suite has been restyled.

- Default output is white-background PNG at 300 dpi. `output.figure_format` accepts `png`, `pdf`, or `both`. `plots.preview_png` is obsolete for these figures; PNG is now the default primary output.
- Use the bundled `figure_style.R` helpers, not a newly selected template on each run. The theme uses sans-serif text, 11-point base typography, quiet axis lines, fixed sizing rules, and explicit legends.
- UMAP: at most four grouping panels per page, two columns. Keep embeddings unchanged. Group order follows the configured fields. Resolution scans keep every resolution and the same underlying embedding; colours do not imply biological equivalence of cluster IDs across resolutions.
- Dotplot: at most 30 genes by 25 clusters per page. Compute the complete DotPlot data before splitting; fixed expression colour limits and 0–100% size limits apply across every page. Retain every requested available gene and cluster.
- Optional FeaturePlots in 08: `figure_style.feature_genes` supplies the ordered gene list; `figure_style.feature_assay` selects the normalized assay. Default assay is the object's active assay. A normalized data layer is required; plotting does not silently normalize. Six genes per page, three columns; grey-to-red expression with a per-gene scale stated in the subtitle. This is descriptive expression, not a differential-expression test.
- A single page retains the original filename stem with `.png`; multiple pages use `_page01`, `_page02`, etc. Optional PDF uses the identical page stems. Do not leave old-format files in new output directories. Existing historical outputs are not automatically removed.
- Categorical colours are deterministic per label, independent of row order or subset composition. For a project-approved palette use named `figure_style.colors.<metadata_field>` maps and reuse the same config across stages. Arbitrary labels can yield similar hues; inspect high-cardinality legends and supply explicit colours when needed. Do not infer condition semantics or rename categories to obtain a colour.
- `_provenance/figure_colors.tsv` records resolved mappings; `_provenance/figure_status.tsv` records generated files and unavailable requested genes. These are technical records, not primary deliverables.

Example configuration addition (replace the example field names and genes with confirmed project values):

```json
{
  "output": {"figure_format": "png"},
  "figure_style": {
    "colors": {"condition": {"control": "#4477AA", "pah": "#D55E00"}},
    "feature_assay": "RNA",
    "feature_genes": ["PECAM1", "VWF"]
  }
}
```

The addition must be merged with existing output settings, not replace object-format or other fields. Do not add unrequested genes or run new statistical analyses. When an empty result or missing dependency prevents a plot, retain the scientific status and report the reason.
