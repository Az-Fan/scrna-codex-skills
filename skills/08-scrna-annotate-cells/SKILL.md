---
name: 08-scrna-annotate-cells
description: Prepare and apply cluster-level scRNA-seq annotations using existing cluster markers or an explicit fallback marker calculation, canonical-marker and UMAP review, and a required human-confirmed decision table. Use after clustering and marker discovery to create broad, fine, and optional state labels without conflating cluster IDs with biological annotations.
---

# Annotate scRNA Cells

Use two explicit actions. Never place candidate labels directly into the object.

## Prepare review

Set `workflow.action=prepare_review` and provide a confirmed cluster column and UMAP reduction. Prefer `input.markers` from `07-scrna-find-cluster-markers`; when it is absent, the executor preserves the former behavior and can calculate markers. Set `clustering.compute_if_missing=true` only to retain the legacy cluster-and-prepare mode.

The action writes `_provenance/annotation_review_record.json` with a unique review ID, source-object hash, prepared-object hash, template/marker hashes and cluster evidence. Apply must read the generated `clustered_object.qs`, whose bytes are bound in this record; the original source object may differ after preparation. The action writes the complete marker table, an annotation-review table, the clustered object, cluster/sample UMAPs, and an optional canonical-marker dot plot. Review sample-restricted clusters, conflicting lineage markers, QC states, doublets, and identity-versus-state distinctions.

## Apply confirmed decisions

Also bind the reviewed record bytes using mandatory `approval.review_record_sha256`; recompute and obtain renewed approval when the record changes.

Obtain the real user's explicit approval of the specific final decision TSV and reviewed run before applying labels. Present the completed table and its SHA256 for that approval; the agent must never mark its own proposal approved or substitute an earlier approval of another version. A pre-existing user approval is usable only when it covers this exact final decision version and review run.

Set `workflow.action=apply_confirmed`, `approval.status=approved`, `approval.source=human`, `approval.review_run_id` from the prepare record and `approval.decision_sha256` for the final approved TSV. Provide `input.review_record`, the bound prepared `input.object`, and the reviewed TSV using [references/config.apply.example.json](references/config.apply.example.json). Every object cluster must occur exactly once and every row must have an allowed confirmed decision. Broad and fine labels must be present and nonblank after trimming whitespace. The action writes broad, fine, and optional state labels to a derivative object and always produces the complete cell-level annotation table, cluster summary, annotated UMAP, cluster/sample/condition audit UMAP, session information, and manifest.

The template may be edited and its label columns may be renamed as configured. Approve the hash of the final edited TSV, not the preparation template hash. Python validation and direct R execution both reject changed decisions, another review ID, a different cluster column, or an object whose bytes differ from the reviewed prepared object. `_provenance/annotation_approval_record.json` preserves the binding and supplied approval reason/provenance. These hashes detect content drift; the `source=human` declaration does not prove that a human authorized the operation.

Older unbound review directories require a fresh `prepare_review`; no external or legacy bypass is supported.

## Guardrails

- Reuse the marker output from skill 07 instead of recomputing it when available.
- Preserve the legacy marker and clustering calculation only as explicit fallback modes.
- Keep cluster IDs, broad identity, fine identity, activation state, QC status, and exclusion decisions separate.
- Do not hard-code tissue-specific markers into execution logic.
- Do not delete cells during annotation.
- Do not accept partial, duplicated, provisional, or unconfirmed decision tables.

Read [references/annotation-review.md](references/annotation-review.md) before preparing or applying decisions. Dry-run `scripts/run.py --config <config>` first, then execute only the reviewed action with `--execute`.

If execution may exceed 10 minutes, read [references/long-running-execution.md](references/long-running-execution.md) and use `scripts/run_in_tmux.py` for the confirmed execution.

## Result organization

Keep primary figures, complete scientific tables, and review decisions directly accessible. Store execution manifests, session information, logs, and workflow state under `_provenance/`; do not list them as primary results. Read [references/output-layout.md](references/output-layout.md) when configuring outputs, locating legacy records, or adding custom plots and diagnostics.

## Fixed figures

Read [references/figure-style.md](references/figure-style.md) before plotting or configuring figure outputs. Use the bundled paper_v1 templates with PDF by default and deterministic pagination/colours. Preserve the existing figure families listed in [references/figure-catalog.md](references/figure-catalog.md); do not choose new chart styles on each run.

Result handoff: start with `<output_dir>/RESULTS.md` from the execute runner. It links scientific deliverables and reports current run status and limitations; follow [output-layout.md](references/output-layout.md) for retained results and technical records.

Figure output: select one format per run, default PDF; follow [figure-output.md](references/figure-output.md). Never automatically export both PDF and PNG.
