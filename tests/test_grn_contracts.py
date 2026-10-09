"""GRN resource, stage, membership and preservation contracts."""
import copy
import csv
import gzip
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "toolkit/python"))
from grn_contract import SKILL, validate_config
from scrna_runtime import nested_get as get, default_argv, expected_artifacts, resolved_rscript
from grn_pipeline import binding_config, file_records, execute, export_regulons


class GRNContracts(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.config = json.loads((ROOT / "skills/17-scrna-infer-grn/references/config.example.json").read_text())
        for field in ["object"]:
            p = self.root / field; p.write_text("fixture")
            self.config["input"][field] = str(p)
        for field in ["tf_list", "motif_annotations"]:
            p = self.root / field; p.write_text("fixture")
            self.config["resources"][field] = str(p)
        p = self.root / "ranking.feather"; p.write_text("fixture")
        self.config["resources"]["ranking_databases"] = [str(p)]
        self.config["output_dir"] = str(self.root / "output")

    def test_valid_prepare_and_reviewed_infer(self):
        self.assertEqual(validate_config(self.config, get), [])
        self.config["workflow"]["action"] = "infer"
        self.assertTrue(any("review_reason" in x for x in validate_config(self.config, get)))
        self.config["workflow"]["review_reason"] = "reviewed memberships"
        self.assertEqual(validate_config(self.config, get), [])

    def test_resource_mismatch_and_duplicate_databases(self):
        for role in ["species", "genome", "gene_identifier"]:
            cfg = copy.deepcopy(self.config)
            cfg["resources"][role] = "other"
            self.assertTrue(any(role in x for x in validate_config(cfg, get)))
        self.config["resources"]["ranking_databases"] *= 2
        self.assertTrue(any("duplicate paths" in x for x in validate_config(self.config, get)))

    def test_output_cannot_contain_source_files(self):
        self.config["output_dir"] = str(self.root)
        self.assertTrue(any("outside output_dir" in x for x in validate_config(self.config, get)))

    def test_memberships_and_raw_counts_must_be_declared(self):
        cfg = copy.deepcopy(self.config)
        del cfg["metacell"]["column"]
        cfg["input"]["counts_source"] = "normalized"
        errors = validate_config(cfg, get)
        self.assertTrue(any("metacell.column" in x for x in errors))
        self.assertTrue(any("raw_umi" in x for x in errors))
        cfg["inference"]["mode"] = "single_cell"
        self.assertFalse(any("metacell.column" in x for x in validate_config(cfg, get)))

    def test_binding_ignores_stage_but_binds_aggregation_and_scoring(self):
        infer = copy.deepcopy(self.config)
        infer["workflow"] = {"action": "infer", "review_reason": "reviewed"}
        self.assertEqual(binding_config(self.config), binding_config(infer))
        infer["aucell"]["batch_size"] += 1
        self.assertNotEqual(binding_config(self.config), binding_config(infer))
        before = file_records(self.config)
        Path(self.config["resources"]["tf_list"]).write_text("changed")
        self.assertNotEqual(before, file_records(self.config))

    def test_missing_explicit_runtime_has_no_fallback(self):
        self.config["runtime"]["pyscenic_python"] = "/missing/python"
        self.assertIsNone(default_argv(SKILL, self.root / "config.json", self.config))
        self.config["runtime"]["pyscenic_python"] = sys.executable
        self.assertEqual(default_argv(SKILL, self.root / "config.json", self.config)[0], sys.executable)
        self.assertNotIn("single_cell_activity", expected_artifacts(SKILL, self.config))

    def test_run_preserves_existing_inference_output(self):
        output = Path(self.config["output_dir"])
        output.mkdir()
        adjacency = output / "adjacencies.tsv"
        adjacency.write_text("previous network\n")
        self.config["workflow"]["action"] = "run"
        with patch("grn_pipeline.importlib.metadata.version", return_value="fixture"):
            with self.assertRaisesRegex(ValueError, "already exists"):
                execute(self.config)
        self.assertEqual(adjacency.read_text(), "previous network\n")
        self.assertFalse((output / "_provenance/grn_input").exists())

    def test_infer_rejects_changed_prepared_artifact(self):
        technical = Path(self.config["output_dir"]) / "_provenance"
        prepared = technical / "grn_input"
        prepared.mkdir(parents=True)
        self.config["workflow"] = {"action": "infer", "review_reason": "reviewed"}
        record = {"inputs": file_records(self.config), "config": binding_config(self.config),
                  "versions": {p: "fixture" for p in ["pyscenic", "arboreto", "ctxcore", "numpy", "pandas", "loompy", "scipy"]},
                  "prepared_files": {"inference.loom": "original_sha"}}
        (technical / "grn_prepare_record.json").write_text(json.dumps(record))
        (prepared / "inference.loom").write_text("tampered")
        with patch("grn_pipeline.importlib.metadata.version", return_value="fixture"):
            with self.assertRaisesRegex(ValueError, "artifact changed"):
                execute(self.config)
        self.assertFalse((technical.parent / "adjacencies.tsv").exists())

    def test_empty_motif_leading_edges_are_audited_before_aggregation(self):
        import pandas as pd
        output = Path(self.config["output_dir"])
        prepared = output / "_provenance/grn_input"
        prepared.mkdir(parents=True)
        (prepared / "genes.tsv").write_text("g1\ng2\n")
        motifs = pd.DataFrame({("Enrichment", "TargetGenes"): [[], [("g1", 1.), ("g2", 2.)]]},
                              index=pd.MultiIndex.from_tuples([("TF", "empty"), ("TF", "valid")], names=["TF", "MotifID"]))
        regulon = SimpleNamespace(name="TF(+)", transcription_factor="TF", genes=["g1", "g2"], gene2weight={"g1": 1., "g2": 2.}, context={"activating"})
        def transform(frame):
            self.assertEqual(list(frame.index), [("TF", "valid")])
            return [regulon]
        with patch.dict(sys.modules, {"pyscenic.utils": SimpleNamespace(load_motifs=lambda *args, **kwargs: motifs),
                                     "pyscenic.transform": SimpleNamespace(df2regulons=transform)}):
            self.assertEqual(export_regulons({"regulons": {"min_genes": 2}}, output), 1)
        audit = pd.read_csv(output / "motif_row_status.tsv", sep="\t")
        self.assertEqual(list(audit.status), ["empty_leading_edge", "retained"])
        self.assertEqual(len(pd.read_csv(output / "regulon_targets.tsv", sep="\t")), 2)

    def test_matrix_metacells_preserve_ids_and_sum_within_samples(self):
        import numpy as np
        from scipy.io import mmwrite, mmread
        from scipy import sparse
        cfg = copy.deepcopy(self.config)
        cfg["input"] = {"type": "matrix", "counts_source": "raw_umi", "species": "mouse", "genome": "mm10", "gene_identifier": "MGI_symbol", "orientation": "cells_by_genes"}
        cfg["metadata"] = {"sample": "sample", "cell_type": "cell_type", "condition": "condition"}
        cfg["metacell"]["column"] = "members"
        matrix = np.arange(1, 49).reshape(8, 6)
        mmwrite(self.root / "counts.mtx", sparse.csr_matrix(matrix))
        with gzip.open(self.root / "counts.mtx.gz", "wb") as handle:
            handle.write((self.root / "counts.mtx").read_bytes())
        cells = ["001", "01", "NA", "nan", "c5", "c6", "c7", "c8"]
        (self.root / "features.tsv").write_text("\n".join(f"g{i}" for i in range(6)) + "\n")
        (self.root / "barcodes.tsv").write_text("\n".join(cells) + "\n")
        with (self.root / "metadata.tsv").open("w") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["cell_id", "sample", "cell_type", "condition", "members"])
            for i in reversed(range(8)):
                writer.writerow([cells[i], "s1" if i < 4 else "s2", "EC", "control" if i < 4 else "case", "a" if i % 4 < 2 else "b"])
        for key in ["features", "barcodes", "metadata"]:
            cfg["input"][key] = str(self.root / (key + ".tsv"))
        cfg["input"]["counts"] = str(self.root / "counts.mtx.gz")
        rscript = resolved_rscript(SKILL, cfg)
        if not rscript:
            self.skipTest("GRN R runtime unavailable")
        config_path = self.root / "matrix.json"
        config_path.write_text(json.dumps(cfg))
        done = subprocess.run([str(rscript), str(ROOT / "toolkit/R/grn_input_and_activity.R"), str(config_path), "prepare"], capture_output=True, text=True)
        self.assertEqual(done.returncode, 0, done.stderr)
        out = Path(cfg["output_dir"])
        with (out / "cell_membership.tsv").open() as handle:
            memberships = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual([row["cell_id"] for row in memberships], cells)
        with (out / "inference_unit_audit.tsv").open() as handle:
            units = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(units), 4)
        self.assertEqual({row["n_cells"] for row in units}, {"2"})
        self.assertEqual({row["sample_purity"] for row in units}, {"1"})
        actual = mmread(out / "_provenance/grn_input/inference_counts.mtx").toarray()
        expected = np.column_stack([matrix[i:i + 2].sum(axis=0) for i in range(0, 8, 2)])
        np.testing.assert_array_equal(actual, expected)
        # 错误计数应阻止聚合，不静默修正。
        cfg["output_dir"] = str(self.root / "negative")
        matrix[0, 0] = -1
        mmwrite(self.root / "negative.mtx", sparse.csr_matrix(matrix))
        cfg["input"]["counts"] = str(self.root / "negative.mtx")
        config_path.write_text(json.dumps(cfg))
        failed = subprocess.run([str(rscript), str(ROOT / "toolkit/R/grn_input_and_activity.R"), str(config_path), "prepare"], capture_output=True, text=True)
        self.assertNotEqual(failed.returncode, 0)
        self.assertIn("non-negative integer", failed.stderr)

    def test_seurat_factor_purity_and_dense_counts(self):
        import numpy as np
        from scipy.io import mmread
        cfg = copy.deepcopy(self.config)
        cfg["metadata"] = {"sample": "sample", "cell_type": "cell_type", "condition": "condition", "batch": "batch"}
        cfg["metacell"] = {"column": "members", "min_cells": 1}
        rscript = resolved_rscript(SKILL, cfg)
        if not rscript:
            self.skipTest("GRN R runtime unavailable")
        creator = self.root / "create.R"
        creator.write_text('''args <- commandArgs(TRUE)
library(Seurat)
x <- matrix(1:48, 6, 8, dimnames=list(paste0("g", 1:6), paste0("c", 1:8)))
obj <- CreateSeuratObject(x)
obj$sample <- factor(rep(c("s1", "s2"), each=4), levels=c("unused", "s1", "s2"))
obj$cell_type <- factor(rep("EC", 8), levels=c("unused", "EC"))
obj$condition <- factor(rep(c("control", "case"), each=4))
obj$batch <- factor(rep("batch1", 8), levels=c("unused", "batch1"))
obj$members <- factor(rep(c("a", "a", "b", "b"), 2))
if (args[[2]] == "dense") LayerData(obj, assay="RNA", layer="counts") <- x
stopifnot(validObject(obj))
saveRDS(obj, args[[1]])
''')
        for representation in ["sparse", "dense"]:
            with self.subTest(representation=representation):
                obj = self.root / (representation + ".rds")
                created = subprocess.run([str(rscript), str(creator), str(obj), representation], capture_output=True, text=True)
                self.assertEqual(created.returncode, 0, created.stderr)
                cfg["input"]["object"] = str(obj)
                cfg["output_dir"] = str(self.root / representation)
                config_path = self.root / (representation + ".json")
                config_path.write_text(json.dumps(cfg))
                done = subprocess.run([str(rscript), str(ROOT / "toolkit/R/grn_input_and_activity.R"), str(config_path), "prepare"], capture_output=True, text=True)
                self.assertEqual(done.returncode, 0, done.stderr)
                out = Path(cfg["output_dir"])
                with (out / "inference_unit_audit.tsv").open() as handle:
                    units = list(csv.DictReader(handle, delimiter="\t"))
                self.assertEqual(len(units), 4)
                for column in ["sample", "cell_type", "condition", "batch"]:
                    self.assertEqual({row[column + "_purity"] for row in units}, {"1"})
                matrix = np.arange(1, 49).reshape(6, 8, order="F")
                expected = np.column_stack([matrix[:, i:i + 2].sum(axis=1) for i in range(0, 8, 2)])
                np.testing.assert_array_equal(mmread(out / "_provenance/grn_input/inference_counts.mtx").toarray(), expected)


if __name__ == "__main__":
    unittest.main()
