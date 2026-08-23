---
name: infer-scrna-regulatory-networks
description: Build and interpret scRNA-seq gene regulatory networks and regulon activity using metacells, pySCENIC-compatible workflows, regulon scoring, RSS, differential regulon analysis, and target-to-TF screening. Use when investigating transcription-factor programs within cell types or conditions and when migrating the sc06 GRN workflow to another species or project.
---

# Infer scRNA Regulatory Networks

Treat regulons as model-derived hypotheses requiring database provenance and orthogonal evidence.

## Workflow

1. Define species, gene identifiers, target populations, sample structure, and regulatory question.
2. Verify the TF list, motif rankings, motif annotations, and genome/database compatibility.
3. Filter genes transparently and create metacells only when their grouping does not erase the contrast of interest.
4. Run network inference and motif pruning with recorded seeds, versions, and database paths.
5. Convert regulons to stable tabular/GMT representations and score activity.
6. Evaluate RSS, cell-type specificity, condition association, variance components, and target overlap.
7. Screen TFs for specified targets or gene lists and link candidates to pathways and DEG evidence.
8. Report database coverage, dropped genes, and robustness limitations.

## Guardrails

- Never reuse mouse motif resources for human data or vice versa.
- Do not interpret regulon activity as direct proof of TF activation.
- Avoid condition inference without biological-sample replication.
- Keep unsupervised discovery separate from hypothesis-driven ranking.

Read [references/grn-provenance.md](references/grn-provenance.md) before running inference.

