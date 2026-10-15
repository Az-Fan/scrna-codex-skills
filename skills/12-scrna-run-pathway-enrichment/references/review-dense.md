# Dense review

Use a fresh output_dir and existing complete skill 12 and skill 11 tables. Config fields: enrichment_table, differential_table, mapping_table (skill 12 SYMBOL/ENTREZID mapping), output_dir, population, comparison_id, title, databases, top_per_direction, fdr_threshold, max_jaccard, max_genes. See config.review.example.json.

Select GSEA terms by FDR within database and NES direction; suppress redundant leading edges using Jaccard, retaining decisions. This mode does not include ORA and does not merge terms across databases. Information table shows database, pathway, |NES| bar, signed NES, FDR, leading-edge count/tested set size, and four leading-edge examples sorted by gene log2FC in the NES direction. Gene examples need not individually pass DE thresholds.

Resolve leading-edge Entrez IDs via the exact enrichment mapping and differential table; retain original_gene display symbols where supplied. If a mapped ID has multiple input genes, retain greatest absolute log2FC deterministically for display. Report unmapped IDs rather than fabricating genes.

Membership matrix reserves up to two available mapped core genes per pathway, then fills its column budget with shared core genes, breaking ties by absolute log2FC. Fail if the budget cannot retain the reserved genes. Colour shows gene log2FC only for leading-edge members; grey means absent from that leading edge, not zero expression or nonsignificance. Colour clipping is stated. This is descriptive membership, not causal evidence.

Keep selected pathways, complete core lists, chosen genes, membership cells and selection decisions as TSV. Generate exactly two PDFs and a provenance manifest. Inspect both plots before sharing; do not automatically promote a preview to default.

All modes select one population/comparison. Supply task_id and/or source_table (matching the result column input_table) when combined tables contain multiple source tasks. Unresolved task ambiguity fails. Source IDs are read as strings, preserving leading zeros; scientific columns are validated numerically. Source tables and config must lie outside a fresh output_dir. Tables are never rewritten.

This book mode supports PDF only; an explicit PNG request fails instead of silently exporting a second format. Compact dot/matrix mode supports either PDF or PNG.
