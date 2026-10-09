#!/usr/bin/env python3
"""Validate cNMF raw-count identity and rank-review contracts."""
import copy
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

os.environ.setdefault("MPLCONFIGDIR", "/tmp/scrna-cnmf-test-matplotlib")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "toolkit/python"))
from scrna_runtime import validate, default_argv
from cnmf_contract import SKILL, validate_config
import cnmf_discovery as discovery
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread, mmwrite
import anndata

discovery.np, discovery.pd, discovery.sparse = np, pd, sparse
discovery.mmread, discovery.anndata = mmread, anndata


class CNMFContracts(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.matrix = np.array([[1, 2, 3, 4], [4, 1, 2, 3], [2, 5, 1, 6], [0, 0, 0, 0], [1, 4, 3, 2]], dtype=float)
        mmwrite(self.root / "counts.mtx", self.matrix)
        (self.root / "features.tsv").write_text("g1\ng2\ng3\nzero\ng4\n")
        (self.root / "barcodes.tsv").write_text("c1\nc2\nc3\nc4\n")
        (self.root / "metadata.tsv").write_text("cell_id\tsample\tcondition\nc3\ts2\tcase\nc1\ts1\tcontrol\nc4\ts2\tcase\nc2\ts1\tcontrol\n")
        self.config = {"project": {"id": "fixture"}, "input": {"type": "matrix", "counts_source": "raw_umi", **{x: str(self.root / (x + (".mtx" if x == "counts" else ".tsv"))) for x in ["counts", "features", "barcodes", "metadata"]}}, "metadata": {"sample": "sample", "condition": "condition"}, "cnmf": {"components": [2, 3], "num_highvar_genes": 4, "n_iter": 8}, "output_dir": str(self.root / "output")}
        self.destination = self.root / "prepared"
        self.destination.mkdir()

    def test_shuffled_metadata_and_zero_features(self):
        adata = discovery.load_input(self.config, self.destination)
        self.assertEqual(list(adata.obs_names), ["c1", "c2", "c3", "c4"])
        self.assertEqual(list(adata.var_names), ["g1", "g2", "g3", "g4"])
        np.testing.assert_array_equal(adata.X.toarray(), self.matrix[[0, 1, 2, 4]].T)
        self.assertEqual(adata.obs.loc["c1", "condition"], "control")
        audit = pd.read_csv(self.destination / "feature_status.tsv", sep="\t")
        self.assertEqual(audit.loc[audit.gene == "zero", "status"].iloc[0], "zero_counts")
        self.assertEqual(anndata.read_h5ad(self.destination / "counts.h5ad").shape, (4, 4))

    def test_explicit_transposed_orientation(self):
        mmwrite(self.root / "counts.mtx", self.matrix.T)
        self.config["input"]["orientation"] = "cells_by_genes"
        adata = discovery.load_input(self.config, self.destination)
        np.testing.assert_array_equal(adata.X.toarray(), self.matrix[[0, 1, 2, 4]].T)

    def test_bad_counts_and_zero_cells(self):
        for value in [-1, .5, float("nan")]:
            broken = self.matrix.copy(); broken[0, 0] = value
            mmwrite(self.root / "counts.mtx", broken)
            with self.assertRaisesRegex(ValueError, "integer raw UMI"):
                discovery.load_input(self.config, self.destination)
        broken = self.matrix.copy(); broken[:, 0] = 0
        mmwrite(self.root / "counts.mtx", broken)
        with self.assertRaisesRegex(ValueError, "Zero-count cells"):
            discovery.load_input(self.config, self.destination)

    def test_duplicate_genes_and_missing_cells(self):
        (self.root / "features.tsv").write_text("g1\ng1\ng3\nzero\ng4\n")
        with self.assertRaisesRegex(ValueError, "duplicate identifiers"):
            discovery.load_input(self.config, self.destination)
        (self.root / "features.tsv").write_text("g1\ng2\ng3\nzero\ng4\n")
        p = self.root / "metadata.tsv"; p.write_text(p.read_text().replace("c3", "unexpected"))
        with self.assertRaisesRegex(ValueError, "must match barcodes exactly"):
            discovery.load_input(self.config, self.destination)

    def test_review_and_interpreter_override(self):
        cfg = copy.deepcopy(self.config); cfg["workflow"] = {"action": "consensus"}
        errors, _ = validate(SKILL, cfg, self.root / "config.json")
        self.assertTrue(any("consensus_k" in x for x in errors))
        self.assertTrue(any("selection_reason" in x for x in errors))
        cfg["workflow"]["selection_reason"] = "reviewed ranks"
        cfg["cnmf"]["consensus_k"] = [2]
        cfg["runtime"] = {"cnmf_python": "/missing/explicit/python"}
        self.assertIsNone(default_argv(SKILL, self.root / "config.json", cfg))
        cfg["runtime"]["cnmf_python"] = sys.executable
        self.assertEqual(default_argv(SKILL, self.root / "config.json", cfg)[0], sys.executable)
        self.assertEqual(validate(SKILL, cfg, self.root / "config.json")[0], [])

    def test_sample_summary_uses_cell_ids(self):
        usage = pd.DataFrame({"p1": [.1, .9, .2, .8], "p2": [.9, .1, .8, .2]}, index=["c4", "c2", "c3", "c1"])
        metadata = pd.read_csv(self.root / "metadata.tsv", sep="\t", index_col=0)
        summary = discovery.summarize_usage(usage, metadata, self.config)
        row = summary.loc[summary["sample"] == "s1"].iloc[0]
        self.assertAlmostEqual(row.p1, .85)
        self.assertEqual(row.n_cells, 2)

    def test_standalone_requires_sample_and_handles_invalid_assay(self):
        cfg = copy.deepcopy(self.config)
        del cfg["metadata"]["sample"]
        cfg["input"]["assay"] = 4
        errors = validate_config(cfg, discovery.get)
        self.assertIn("metadata.sample is required", errors)
        self.assertIn("input.assay must be a non-empty assay name", errors)

    def test_native_consensus_preserves_identifiers_and_values(self):
        from cnmf import cNMF
        from cnmf.cnmf import load_df_from_npz, save_df_to_npz
        model = cNMF(output_dir=str(self.root / "native"), name="fixture")
        cells, genes = ["001", "01", "NA", "nan"], ["001", "01", "NA"]
        raw = pd.DataFrame([[1., 3.], [2., 1.], [3., 2.], [4., 1.]], index=cells, columns=[1, 2])
        scores = pd.DataFrame([[2., -1., 1.], [-1., 3., 2.]], index=[1, 2], columns=genes)
        tpm = scores.abs() * 10
        for key, frame in [("consensus_usages", raw), ("gene_spectra_score", scores), ("gene_spectra_tpm", tpm)]:
            save_df_to_npz(frame, model.paths[key] % (2, "2"))
        cfg = copy.deepcopy(self.config)
        cfg["cnmf"]["density_threshold"] = 2
        metadata = pd.DataFrame({"sample": ["s1", "s1", "s2", "s2"], "condition": ["control", "control", "case", "case"]}, index=cells)
        with patch.object(discovery, "load_df_from_npz", load_df_from_npz, create=True), patch.object(model, "consensus"), patch.object(discovery, "plot_summary"):
            usage = discovery.export_consensus(model, 2, cfg, metadata, self.root, None)
        self.assertEqual(list(usage.index), cells)
        np.testing.assert_allclose(usage.to_numpy(), raw.div(raw.sum(axis=1), axis=0))
        exported = pd.read_csv(self.root / "k_2/spectra_scores.tsv", sep="\t", index_col=0, dtype={"gene": str}, keep_default_na=False)
        self.assertEqual(list(exported.index), genes)
        np.testing.assert_array_equal(exported.to_numpy(), scores.T.to_numpy())

    def test_consensus_registry_accumulates_successes_and_failures(self):
        output = self.root / "output"
        output.mkdir()
        for k in [2, 3]:
            destination = output / f"k_{k}"
            destination.mkdir()
            pd.DataFrame({f"cnmf{k}_1": [.2, .8]}, index=["001", "NA"]).to_csv(destination / "usage.tsv", sep="\t", index_label="cell_id")
        rows, usages = discovery.consensus_registry(output, [{"k": 2, "status": "completed", "reason": ""}, {"k": 4, "status": "failed", "reason": "density filtering"}])
        rows, usages = discovery.consensus_registry(output, [{"k": 3, "status": "completed", "reason": ""}])
        self.assertEqual(set(usages), {2, 3})
        self.assertEqual(list(usages[2].index), ["001", "NA"])
        self.assertEqual([row for row in rows if row["status"] == "failed"][0]["reason"], "density filtering")
        self.assertEqual(len(pd.read_csv(output / "task_status.tsv", sep="\t")), 3)

    def test_ora_keeps_complete_tests_and_audits_unmatched_sets(self):
        gmt = self.root / "sets.gmt"; gmt.write_text("A\tdescription\tg1\nghost\tdescription\tunknown\nB\tdescription\tg4\n")
        scores = pd.DataFrame({"p1": [5, -1, -2, -3], "p2": [-1, -2, -3, 7], "p3": [-1, -2, -3, -4]}, index=["g1", "g2", "g3", "g4"])
        discovery.enrich_programs(scores, 3, gmt, self.destination)
        table = pd.read_csv(self.destination / "program_enrichment.tsv", sep="\t")
        self.assertEqual(len(table), 6)
        self.assertEqual(set(table.n_selected), {0, 1})
        self.assertEqual(set(table.universe_size), {4})
        self.assertEqual(len(table.loc[table.overlap == 0]), 4)
        self.assertTrue(table.loc[table.program == "p3", "p_value"].eq(1).all())
        self.assertTrue((table.p_adjust >= table.p_value).all())
        coverage = pd.read_csv(self.destination / "gene_set_coverage.tsv", sep="\t")
        self.assertEqual(coverage.loc[coverage.gene_set == "ghost", "status"].iloc[0], "no_overlap")


if __name__ == "__main__":
    unittest.main()
