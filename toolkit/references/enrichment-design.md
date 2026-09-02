# Enrichment design

Run enrichment as an optional but comprehensive second stage. Confirm `enrichment.species` and `enrichment.gene_id_type` before execution.

## Default suite

When enrichment is enabled and `enrichment.databases` is omitted, run all of:

- GO Biological Process, Molecular Function, and Cellular Component.
- KEGG canonical pathways.
- Reactome pathways.
- MSigDB Hallmark gene sets.

Use the installed human or mouse organism annotation package for GO and identifier mapping. Use `msigdbr` gene sets for KEGG, Reactome, and Hallmark so the workflow is reproducible without live pathway API calls. Record package and database versions in `sessionInfo.txt`.

## ORA and GSEA

Run ORA separately for significant Up and Down genes in every requested database. Use genes with `tested=true` as the universe. Mark a direction `skipped_too_few_genes` when fewer than `enrichment.min_input_genes` mapped genes remain. Set pathway P- and Q-value cutoffs to 1 during computation so full tables are retained, then interpret adjusted P values from those tables.

Run GSEA for every database over all finite tested genes ranked by the DESeq2 Wald statistic when present, otherwise by signed log2 fold change. Use `eps=0`, retain full results before reporting thresholds, and apply configurable minimum and maximum gene-set sizes. Collapse duplicate Entrez mappings by retaining the statistic with greatest absolute magnitude.

Always retain `gene_id_mapping.tsv` and `enrichment_status.tsv`. A database failure, empty result, or missing optional resource must not discard DE results or successful enrichment from other databases. Review mapping yield before interpretation.

Write full database-specific tables plus a combined table. Keep figures readable by producing a fixed-size overview with at most three terms per database/method/direction and separate bounded-height detail PDFs for each database and method. Use `enrichment.plot_top_terms` for detail panels, `enrichment.plot_label_width` (default 55) to wrap display labels, and `enrichment.plot_terms_per_page` (default 20) to paginate detail and GSEA PDFs without altering full tables. Produce database-specific GSEA NES plots rather than one unbounded figure. Treat pathway adjusted P values as multiple-tested results. Use ORA input direction or GSEA NES for direction, and avoid treating overlapping pathway names as independent evidence.
