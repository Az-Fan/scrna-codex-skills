# Annotation review table

The preparation table uses `cluster`, `candidate_broad`, `candidate_fine`, `candidate_state`, `evidence`, `conflicts`, `sample_bias`, `qc_flag`, `confidence`, `decision`, and `notes`.

Keep positive evidence and conflicting/exclusion evidence explicit. Record sample restriction, possible doublets, low-quality states, and condition dominance without automatically deleting a cluster or converting an activation state into an anatomical identity.

Before `apply_confirmed`, fill the configured broad and fine label columns, optional state column, and decision column. Every observed cluster must occur exactly once; extra, missing, duplicated, provisional, or blank rows are blocking. Applying labels writes a derivative object and never overwrites cluster IDs.
