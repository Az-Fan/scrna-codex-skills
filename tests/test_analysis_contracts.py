#!/usr/bin/env python3
import importlib.util
import json
import tempfile
import unittest
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
INTEGRATION = ROOT / "skills/05-scrna-benchmark-integration/scripts/integration_python.py"
SPEC = importlib.util.spec_from_file_location("integration_python", INTEGRATION)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)
RUNTIME_SPEC = importlib.util.spec_from_file_location("scrna_runtime", ROOT / "toolkit/python/scrna_runtime.py")
RUNTIME = importlib.util.module_from_spec(RUNTIME_SPEC)
RUNTIME_SPEC.loader.exec_module(RUNTIME)
SCCODA_SPEC = importlib.util.spec_from_file_location("sccoda_adapter", ROOT / "toolkit/python/cell_abundance_sccoda.py")
SCCODA = importlib.util.module_from_spec(SCCODA_SPEC)
SCCODA_SPEC.loader.exec_module(SCCODA)


class RecommendationContractTests(unittest.TestCase):
    def setUp(self):
        self.config = {
            "metrics": {
                "batch_removal": ["ilisi"],
                "biological_conservation": ["label_asw"],
            }
        }
        self.ranking = pd.DataFrame(
            [{"scenario": "none", "pareto_efficient": True, "biology_score": 0.8, "batch_score": 0.2}]
        )

    def test_missing_requested_metrics_is_unresolved(self):
        metrics = pd.DataFrame(
            [{"status": "missing_dependency", "metric": "ilisi", "value": np.nan}]
        )
        confounding = pd.DataFrame([{"perfect_confounding": False}])
        decision = MODULE.recommendation_status(metrics, self.ranking, confounding, self.config)
        self.assertEqual(decision["status"], "unresolved")
        self.assertIsNone(decision["recommended_scenario"])
        self.assertIn("no_requested_metric_completed", decision["reasons"])

    def test_perfect_confounding_blocks_recommendation(self):
        metrics = pd.DataFrame([{"status": "completed", "metric": "ilisi", "value": 0.5}])
        confounding = pd.DataFrame([{"perfect_confounding": True}])
        decision = MODULE.recommendation_status(metrics, self.ranking, confounding, self.config)
        self.assertEqual(decision["status"], "unresolved")
        self.assertIn("batch_condition_perfectly_confounded", decision["reasons"])

    def test_single_pareto_scenario_can_resolve(self):
        metrics = pd.DataFrame([{"scenario": "none", "batch_variable": "batch", "biological_label": "cell_type", "status": "completed", "metric": metric, "value": 0.5} for metric in ["ilisi", "label_asw"]])
        confounding = pd.DataFrame([{"perfect_confounding": False}])
        decision = MODULE.recommendation_status(metrics, self.ranking, confounding, self.config)
        self.assertEqual(decision["status"], "resolved")
        self.assertEqual(decision["recommended_scenario"], "none")

    def test_malformed_completed_rows_do_not_resolve_a_ranking(self):
        metrics = pd.DataFrame([{"status":"completed","metric":"ilisi","value":.5}])
        decision = MODULE.recommendation_status(metrics,self.ranking,pd.DataFrame(),self.config)
        self.assertEqual(decision["status"],"unresolved")

    def test_partial_high_scores_do_not_create_a_winner(self):
        config = {"metrics": {"batch_removal": ["ilisi"], "biological_conservation": ["label_asw", "nmi"]}}
        rows = []
        for scenario, value in [("partial", .99), ("none", .5)]:
            for metric in ["ilisi", "label_asw", "nmi"]:
                failed = scenario == "partial" and metric == "nmi"
                rows.append(dict(scenario=scenario, batch_variable="batch", biological_label="cell_type", metric=metric,
                    metric_group="batch_removal" if metric == "ilisi" else "biological_conservation",
                    value=np.nan if failed else value, status="failed" if failed else "completed"))
        metrics = pd.DataFrame(rows)
        _, ranking = MODULE.summarize(metrics, config)
        decision = MODULE.recommendation_status(metrics, ranking, pd.DataFrame(), config)
        self.assertEqual(decision["status"], "unresolved")
        self.assertIsNone(decision["recommended_scenario"])
        self.assertFalse(ranking.set_index("scenario").loc["partial", "pareto_efficient"])
        weighted_config = dict(config,scoring={"enabled":True,"batch_weight":.3,"biology_weight":.7})
        weighted,_ = MODULE.summarize(metrics,weighted_config)
        self.assertTrue(pd.isna(weighted.set_index("scenario").loc["partial","weighted_score"]))

    def test_an_entire_missing_label_metric_grid_blocks_recommendation(self):
        config = dict(self.config, metadata={"batch_variables": ["batch"], "biological_labels": ["broad", "fine"]})
        metrics = pd.DataFrame([dict(scenario="none", batch_variable="batch", biological_label="broad", metric=metric,
            status="completed", value=.5) for metric in ["ilisi", "label_asw"]])
        decision = MODULE.recommendation_status(metrics, self.ranking, pd.DataFrame(), config)
        self.assertEqual(decision["status"], "unresolved")
        self.assertIn("incomplete_requested_metric_evidence", decision["reasons"])

    def test_precomputed_reductions_have_distinct_ids_and_duplicate_grids_are_rejected(self):
        config = {"benchmark": {"methods": [{"name": "precomputed", "reduction": "pca"}, {"name": "precomputed", "reduction": "harmony"}]}}
        scenarios = MODULE.expand_methods(config)
        self.assertEqual(len({x["id"] for x in scenarios}), 3)
        with self.assertRaises(ValueError):
            MODULE.expand_methods({"benchmark": {"methods": [{"name": "harmony", "parameter_grid": {"theta": [2, 2.0]}}]}})

    def test_validation_allows_an_explicit_empty_plot_selection(self):
        with tempfile.TemporaryDirectory() as temporary:
            obj = Path(temporary)/"object.rds"; obj.write_bytes(b"fixture")
            config = {"project":{"id":"test"},"input":{"object":str(obj)},"metadata":{"sample":"sample","batch_variables":["batch"]},
                "benchmark":{"methods":[{"name":"none"}]},"metrics":{"batch_removal":[],"biological_conservation":[]},"plots":[],"output_dir":str(Path(temporary)/"out")}
            errors,_ = RUNTIME.validate("05-scrna-benchmark-integration",config,Path(temporary)/"config.json")
            self.assertEqual(errors,[])


class ExecutorAuditTests(unittest.TestCase):
    def test_failure_updates_manifest_and_archives_old_results(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            obj = root / "object.rds"; obj.write_bytes(b"fixture")
            output = root / "output"; output.mkdir()
            (output / "all_comparisons.tsv").write_text("old result\n")
            config = {"project": {"id": "test"}, "input": {"object": str(obj)},
                "metadata": {"sample": "sample", "condition": "condition"},
                "comparison": {"numerator": "case", "denominator": "control"}, "output_dir": str(output),
                "analysis": {"counts_source": {"kind": "raw_umi"}},
                "executor": {"argv": [sys.executable, "-c", "import json,sys; from pathlib import Path; p=Path(sys.argv[1]); p.write_text(json.dumps({'exit_status':0,'skill':'11-scrna-run-differential-analysis'})); sys.exit(3)", str(output / "_provenance/run_manifest.json")]}}
            config_path = root / "config.json"; config_path.write_text(json.dumps(config))
            done = subprocess.run([sys.executable, str(ROOT / "skills/11-scrna-run-differential-analysis/scripts/run.py"), "--config", str(config_path), "--execute"], capture_output=True, text=True)
            self.assertEqual(done.returncode, 3, done.stdout + done.stderr)
            manifest = json.loads((output / "_provenance/run_manifest.json").read_text())
            self.assertEqual(manifest["exit_status"], 3)
            self.assertEqual(manifest["status"], "failed")
            self.assertFalse((output / "all_comparisons.tsv").exists())
            self.assertEqual((Path(manifest["previous_output"]) / "all_comparisons.tsv").read_text(), "old result\n")

    def test_archiving_does_not_move_an_input_in_the_output_directory(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "output"; output.mkdir()
            obj = output / "object.rds"; obj.write_bytes(b"source data")
            with self.assertRaises(ValueError):
                RUNTIME.archive_previous_output("11-scrna-run-differential-analysis", {"input": {"object": str(obj)}}, Path(temporary) / "config.json", output, "run")
            self.assertEqual(obj.read_bytes(), b"source data")


class CovariateHandoffTests(unittest.TestCase):
    def test_numeric_category_and_continuous_covariate_keep_their_r_types(self):
        metadata = pd.DataFrame({"age":[41,42,43,44],"donor":[101,102,101,102]})
        observed = SCCODA.apply_covariate_types(metadata,["age","donor"],{"age":"continuous","donor":"categorical"})
        self.assertTrue(pd.api.types.is_numeric_dtype(observed.age))
        self.assertIsInstance(observed.donor.dtype,pd.CategoricalDtype)
        self.assertTrue(pd.api.types.is_numeric_dtype(metadata.donor))

    def test_sample_and_declared_category_strings_survive_tsv_roundtrip(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary)/"metadata.tsv"
            path.write_text("sample\tdonor\tage\n001\t01\t41\n1\t1\t42\n")
            observed = SCCODA.read_sample_table(path,{"donor":str})
            self.assertEqual(observed.index.tolist(),["001","1"])
            self.assertEqual(observed.donor.tolist(),["01","1"])


class V3WorkflowContractTests(unittest.TestCase):
    def test_qc_filter_requires_explicit_approval(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            obj = root / "object.rds"; obj.write_bytes(b"fixture")
            decisions = root / "decisions.tsv"; decisions.write_text("cell_id\tkeep\ncell_1\tTRUE\n", encoding="utf-8")
            config_path = root / "config.json"; config_path.write_text("{}", encoding="utf-8")
            config = {
                "project": {"id": "test"}, "input": {"object": str(obj), "decision_table": str(decisions)},
                "metadata": {"sample": "sample"},
                "decision": {"include_all_true": ["keep"], "expected_retained_cells": 1},
                "approval": {"status": "pending"}, "output_dir": str(root / "out"),
            }
            errors, _ = RUNTIME.validate("04-scrna-apply-qc-filter", config, config_path)
            self.assertIn("approval.status must be exactly 'approved'", errors)

    def test_annotation_apply_requires_decision_contract(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            obj = root / "object.rds"; obj.write_bytes(b"fixture")
            config_path = root / "config.json"; config_path.write_text("{}", encoding="utf-8")
            config = {
                "project": {"id": "test"}, "workflow": {"action": "apply_confirmed"},
                "input": {"object": str(obj)}, "metadata": {"sample": "sample", "cluster": "cluster"},
                "output_dir": str(root / "out"),
            }
            errors, _ = RUNTIME.validate("08-scrna-annotate-cells", config, config_path)
            self.assertTrue(any("apply_confirmed requires input.decisions" in error for error in errors))

    def test_enrichment_entry_rejects_differential_stage(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            table = root / "de.tsv"; table.write_text("gene\tlog2FoldChange\nA\t1\n", encoding="utf-8")
            config_path = root / "config.json"; config_path.write_text("{}", encoding="utf-8")
            config = {
                "project": {"id": "test"}, "input": {"differential_table": str(table)},
                "analysis": {"stage": "differential"}, "output_dir": str(root / "out"),
            }
            errors, _ = RUNTIME.validate("12-scrna-run-pathway-enrichment", config, config_path)
            self.assertIn("12-scrna-run-pathway-enrichment requires analysis.stage=enrichment_only", errors)

    def test_formal_de_requires_raw_count_declaration(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            obj = root / "object.rds"; obj.write_bytes(b"fixture")
            cfg_path = root / "config.json"; cfg_path.write_text("{}")
            config = {"project": {"id": "test"}, "input": {"object": str(obj)},
                      "metadata": {"sample": "sample", "condition": "condition"},
                      "output_dir": str(root / "out"), "analysis": {},
                      "comparisons": [{"numerator": "B", "denominator": "A"}]}
            skill = "11-scrna-run-differential-analysis"
            errors, _ = RUNTIME.validate(skill, config, cfg_path)
            self.assertTrue(any("counts_source.kind" in error for error in errors))
            config["analysis"]["counts_source"] = {"kind": "raw_umi"}
            errors, _ = RUNTIME.validate(skill, config, cfg_path)
            self.assertFalse(errors)
            for assay in ["integrated", "SCT", "harmony", "RNA.corrected"]:
                config["analysis"]["assay"] = assay
                errors, _ = RUNTIME.validate(skill, config, cfg_path)
                self.assertTrue(any("uncorrected" in error for error in errors))
            config["analysis"] = {"method": "seurat_wilcox"}
            errors, _ = RUNTIME.validate(skill, config, cfg_path)
            self.assertFalse(errors)

    def test_de_rejects_invalid_alpha_and_background(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            obj = root / "object.rds"; obj.write_bytes(b"fixture")
            cfg_path = root / "config.json"; cfg_path.write_text("{}")
            config = {"project": {"id": "test"}, "input": {"object": str(obj)},
                      "metadata": {"sample": "sample", "condition": "condition"},
                      "output_dir": str(root / "out"), "analysis": {"counts_source": {"kind": "raw_read"}, "padj_threshold": 1},
                      "comparisons": [{"numerator": "B", "denominator": "A"}], "enrichment": {"universe_mode": "all_expressed"}}
            errors, _ = RUNTIME.validate("11-scrna-run-differential-analysis", config, cfg_path)
            self.assertTrue(any("padj_threshold" in error for error in errors))
            self.assertTrue(any("universe_mode" in error for error in errors))

    def test_cell_abundance_requires_explicit_denominator(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            table = root / "counts.tsv"
            table.write_text("sample\tcondition\tcell_type\tn_cells\ns1\tA\tEC\t10\n", encoding="utf-8")
            config_path = root / "config.json"; config_path.write_text("{}", encoding="utf-8")
            config = {
                "project": {"id": "test"}, "input": {"counts_table": str(table)},
                "metadata": {"sample": "sample", "condition": "condition", "cell_type": "cell_type"},
                "comparisons": [{"id": "B_vs_A", "numerator": "B", "denominator": "A"}],
                "analysis": {"methods": ["propeller"], "denominator": {"mode": "all_input_cells"}},
                "output_dir": str(root / "out"),
            }
            errors, _ = RUNTIME.validate("13-scrna-test-cell-abundance", config, config_path)
            self.assertTrue(any("analysis.denominator.description" in error for error in errors))

    def test_milo_rejects_aggregated_counts_input(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            table = root / "counts.tsv"
            table.write_text("sample\tcondition\tcell_type\tn_cells\ns1\tA\tEC\t10\n", encoding="utf-8")
            config_path = root / "config.json"; config_path.write_text("{}", encoding="utf-8")
            config = {
                "project": {"id": "test"}, "input": {"counts_table": str(table)},
                "metadata": {"sample": "sample", "condition": "condition", "cell_type": "cell_type"},
                "comparisons": [{"id": "B_vs_A", "numerator": "B", "denominator": "A"}],
                "analysis": {"methods": ["milo"], "denominator": {"mode": "all_input_cells", "description": "all cells"}},
                "method_options": {"milo": {"reduction": "pca"}},
                "output_dir": str(root / "out"),
            }
            errors, _ = RUNTIME.validate("13-scrna-test-cell-abundance", config, config_path)
            self.assertTrue(any("milo requires input.object" in error.lower() for error in errors))


if __name__ == "__main__":
    unittest.main()
