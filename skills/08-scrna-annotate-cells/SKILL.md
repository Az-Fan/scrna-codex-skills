---
name: 08-scrna-annotate-cells
description: Prepare and apply cluster-level scRNA-seq annotations using existing cluster markers or an explicit fallback marker calculation, canonical-marker and UMAP review, and a required human-confirmed decision table. Use after clustering and marker discovery to create broad, fine, and optional state labels without conflating cluster IDs with biological annotations.
---

# Annotate scRNA Cells

Use two explicit actions. Never place candidate labels directly into the object.

## Prepare review

Set `workflow.action=prepare_review` and provide a confirmed cluster column and UMAP reduction. Prefer `input.markers` from `07-scrna-find-cluster-markers`; when it is absent, the executor preserves the former behavior and can calculate markers. Set `clustering.compute_if_missing=true` only to retain the legacy cluster-and-prepare mode.

The action writes the complete marker table, an annotation-review table, the clustered object, cluster/sample UMAPs, and an optional canonical-marker dot plot. Review sample-restricted clusters, conflicting lineage markers, QC states, doublets, and identity-versus-state distinctions.

## Apply confirmed decisions

Set `workflow.action=apply_confirmed` and provide the reviewed TSV using [references/config.apply.example.json](references/config.apply.example.json). Every object cluster must occur exactly once and every row must have an allowed confirmed decision. The action writes broad, fine, and optional state labels to a derivative object and always produces the complete cell-level annotation table, cluster summary, annotated UMAP, cluster/sample/condition audit UMAP, session information, and manifest.

## Guardrails

- Reuse the marker output from skill 07 instead of recomputing it when available.
- Preserve the legacy marker and clustering calculation only as explicit fallback modes.
- Keep cluster IDs, broad identity, fine identity, activation state, QC status, and exclusion decisions separate.
- Do not hard-code tissue-specific markers into execution logic.
- Do not delete cells during annotation.
- Do not accept partial, duplicated, provisional, or unconfirmed decision tables.

Read [references/annotation-review.md](references/annotation-review.md) before preparing or applying decisions. Dry-run `scripts/run.py --config <config>` first, then execute only the reviewed action with `--execute`.

If execution may exceed 10 minutes, read [references/long-running-execution.md](references/long-running-execution.md) and use `scripts/run_in_tmux.py` for the confirmed execution.
