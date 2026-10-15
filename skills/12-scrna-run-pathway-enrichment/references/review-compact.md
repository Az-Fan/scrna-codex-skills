# Configuration

Required: `input_table`, `output_dir` (fresh directory), `population`, `comparison_id`.
Optional: `title`, `positive_label`, `negative_label`, `style` (`dot`, `lollipop`, `both`; default dot), `fdr_threshold` (0.05), `top_per_direction` (6), `matrix_per_database` (4), `max_jaccard` (0.6), `terms_per_page` (12), `figure_format` (`pdf` default or `png`). No statistical tests are run. PDF pages are consolidated into one book per figure family; PNG alternative retains numbered page images.

Input is skill 12 `enrichment_all_comparisons.tsv`: database, method (ORA/GSEA), direction (Up/Down for ORA), ID, Description, p.adjust, NES, Count, RichFactor, geneID, core_enrichment, population, comparison_id. RichFactor may be derived from Count and the numerator of BgRatio when absent. Positive/negative group labels describe GSEA enrichment; ORA uses literal Up/Down input genes.

The matrix selects representatives round-robin from significant GSEA-positive, GSEA-negative, ORA-Up and ORA-Down lists, sorted within each by FDR then stable pathway ID. All three method cells for each chosen pathway are fetched from complete results, including nonsignificant ones. Candidate and suppression tables record selection. This is a display budget, not an exhaustive biological summary.

All modes select one population/comparison. Supply task_id and/or source_table (matching the result column input_table) when combined tables contain multiple source tasks. Unresolved task ambiguity fails. Source IDs are read as strings, preserving leading zeros; scientific columns are validated numerically. Source tables and config must lie outside a fresh output_dir. Tables are never rewritten.
