#!/usr/bin/env python3
import importlib.util
import json
import tempfile
import unittest
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
        metrics = pd.DataFrame([{"status": "completed", "metric": "ilisi", "value": 0.5}])
        confounding = pd.DataFrame([{"perfect_confounding": False}])
        decision = MODULE.recommendation_status(metrics, self.ranking, confounding, self.config)
        self.assertEqual(decision["status"], "resolved")
        self.assertEqual(decision["recommended_scenario"], "none")


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


if __name__ == "__main__":
    unittest.main()
