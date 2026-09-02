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


if __name__ == "__main__":
    unittest.main()
