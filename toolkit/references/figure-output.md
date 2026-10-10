# Single-format scientific figures

Every scRNA skill uses one scientific figure format per run. Default: `output.figure_format="pdf"`. Multi-page QC, abundance and cNMF families use PDF containers; other paginated families retain their page filenames with PDF extensions. Choose `"png"` explicitly only when standalone raster files are needed. `"both"`, arrays of formats and conflicting settings are rejected. A figure-free stage remains figure-free.

The legacy fields `plots.figure_format`, `enrichment.plot_format` and `visualization.figure_format` select that same run-wide format. Keep them consistent with `output.figure_format`, or use only the canonical field. For QC skill 03, `output.preview_png=true` is a compatibility alias selecting PNG instead of PDF; it never adds PNG copies beside PDF. Full diagnostics add audit tables, and export plots only in the selected format.

Preserve every eligible figure family, page, tested result and provenance record. Internal temporary raster layers embedded inside PDF and technical caches are not duplicate scientific exports. Never silently delete existing historical figures when revising a run; use a fresh output directory to avoid mixing old and new formats. Result handoff links only the current run's outputs.
