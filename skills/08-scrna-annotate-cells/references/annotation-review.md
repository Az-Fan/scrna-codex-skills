# Annotation review table

The preparation table uses `cluster`, `candidate_broad`, `candidate_fine`, `candidate_state`, `evidence`, `conflicts`, `sample_bias`, `qc_flag`, `confidence`, `decision`, and `notes`.

Keep positive evidence and conflicting/exclusion evidence explicit. Record sample restriction, possible doublets, low-quality states, and condition dominance without automatically deleting a cluster or converting an activation state into an anatomical identity.

Before `apply_confirmed`, fill the configured broad and fine label columns, optional state column, and decision column. Every observed cluster must occur exactly once; extra, missing, duplicated, provisional, or blank rows are blocking. Applying labels writes a derivative object and never overwrites cluster IDs.

## Version-bound approval

Bind the review record itself with `approval.review_record_sha256`, calculated from its final reviewed bytes. A record changed after approval requires renewed approval even when its run ID remains the same. This field is mandatory alongside the final decision TSV hash.

Preparation writes `_provenance/annotation_review_record.json` (schema version 1). It binds a unique `review_run_id` to the actual `clustered_object.qs` for apply, its SHA256, the cluster column/IDs/counts and assignment digest, the original source object, and the initially emitted review template and marker-table hashes. The source object may be normalized or clustered during preparation, so apply uses the prepared object, not the original file.

Finish editing the decision TSV before asking the user to approve its exact final contents and SHA256. It is valid to edit the original template; its preparation hash is historical context, not the hash of the approved final decisions. Record `approval.status=approved`, `approval.source=human`, `approval.review_run_id`, `approval.decision_sha256`, and the actual approval reason. `input.review_record` supplies the prepare record. Do not invent or self-authorize a human approval. The agent must retain the user's approval context; hashes and config declarations cannot establish its authenticity.

The apply executor records the final decision hash, review-record hash/ID, apply-object hash, source, reason and optional provenance in `_provenance/annotation_approval_record.json`. Cluster labels must still pass all semantic decision checks. Changes to the table, review run, cluster column or prepared object require the matching review and a new user approval of the final decision version. Old reviews without this record must be prepared again. Test fixtures may emulate the production approval fields only with explicit `approval.provenance=synthetic_automated_test_not_human_approval`; that test declaration is not genuine human authorization and is never a production bypass.
