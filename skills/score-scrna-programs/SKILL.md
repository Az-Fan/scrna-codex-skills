---
name: score-scrna-programs
description: Quantify and interpret pathway activity, metabolism, curated gene sets, and de novo gene programs in scRNA-seq data using approaches such as PROGENy, Hallmark scoring, scMetabolism, module scores, and cNMF. Use when comparing cellular programs across cell types, samples, or conditions or selecting and interpreting cNMF programs.
---

# Score scRNA Programs

Match the method to the question and separate cell-level visualization from sample-level inference.

## Workflow

1. Define whether the target is a curated pathway, metabolism, a known signature, or de novo programs.
2. Confirm species, identifier mapping, assay/layer, population, and comparison design.
3. Use PROGENy or curated signatures for hypothesis-led activity; use cNMF for recurrent de novo programs.
4. For cNMF, export non-negative counts, test several ranks, assess stability, and retain rank-selection evidence.
5. Attach scores with namespaced fields and preserve the gene-set version and method parameters.
6. Summarize scores per sample and population before testing condition effects.
7. Interpret programs using top genes, enrichment, cell-type distribution, sample consistency, and condition association.

## Guardrails

- Do not call a program condition-specific from pooled cell-level significance alone.
- Do not choose cNMF rank only because one heatmap looks clean.
- Do not compare scores from different methods as though their numeric scales were identical.
- Preserve full loading and usage matrices.

Read [references/program-methods.md](references/program-methods.md) to choose a method.

