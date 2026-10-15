# GO semantic modules

Inputs: complete enrichment_table, mapping_table, differential_table; one population, comparison_id, database and matching GO ontology. Config also supplies title, output_dir (fresh), semantic_measure="Wang", similarity_threshold, display_neighbors (3), rows_per_page (<=30), random_seed. See config.review.example.json. R dependencies: GOSemSim, GO.db, the explicitly selected orgdb_package, igraph, ggplot2, jsonlite; ggrepel optional; gs merges variable-height PDF pages. Reuse the existing analysis environment.

Calculate pairwise Wang similarity over all unique significant GO terms, not previously selected Top terms. Fail on unresolved similarity rather than silently dropping IDs. Group GSEA/ORA and NES/input directions independently. Threshold the semantic graph, cluster with weighted Louvain resolution 1 and explicit seed reset per group. Canonical module IDs sort by term count, best term FDR and GO ID. Singletons remain distinct modules and count toward all totals.

The optional similarity_cache is an RDS list containing similarity, ontology, measure, orgdb_package, orgdb_version, go_db_version, gosemsim_version. Verify coverage and annotation/package versions before reuse. Record cache provenance. Never substitute hit-set or leading-edge Jaccard while calling the result GO semantic similarity.

Choose the similarity threshold after inspecting a small scan of module counts, singleton counts and largest module sizes. Require explicit similarity_threshold, random_seed and orgdb_package; verify database matches ontology (GO_BP/BP, GO_MF/MF or GO_CC/CC). Use the species annotation package matching the original analysis. No threshold scan is automatically performed; never claim one was run. Record any explicit scan without treating the chosen threshold as optimal or the modules as stable. Do not compare module IDs across groups as if they were matched themes.

Module representative = semantic medoid within that group, ties broken by FDR and GO ID. Representative labels are original GO descriptions, not curated theme annotations or cell identity claims. Median NES and best term FDR are descriptive; modules have no newly tested P value. Shared contribution genes and their term frequencies aid interpretation but do not determine these clusters.

The network retains every significant node. Its full threshold graph determines clusters and layout; only the union of strongest display_neighbors edges per node is rendered to reduce clutter. This union can give individual nodes more than display_neighbors visible connections. Export all edges with displayed flags, coordinates and module colour maps. Label multi-term modules with ID, size and representative description; singleton nodes remain visible and explicitly counted. Graph distance is not an effect-size metric.

Deliver module_overview.pdf, semantic_networks.pdf and module_term_atlas.pdf, with group_counts.tsv, module_summary.tsv and term_module_membership.tsv. The term atlas contains all significant rows, not Top N. Preserve source statistics, validate complete per-group coverage and module counts, and visually inspect at least one GSEA and one ORA network. Keep this as an optional preview unless the user adopts it.

All modes select one population/comparison. Supply task_id and/or source_table (matching the result column input_table) when combined tables contain multiple source tasks. Unresolved task ambiguity fails. Source IDs are read as strings, preserving leading zeros; scientific columns are validated numerically. Source tables and config must lie outside a fresh output_dir. Tables are never rewritten.

This book mode supports PDF only; an explicit PNG request fails instead of silently exporting a second format. Compact dot/matrix mode supports either PDF or PNG.

Method references: [GOSemSim Wang semantics](https://yulab-smu.top/biomedical-knowledge-mining-book/go-semantic-similarity.html), [igraph weighted Louvain](https://r.igraph.org/reference/cluster_louvain.html).
