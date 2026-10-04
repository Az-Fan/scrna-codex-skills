import copy
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("runtime", ROOT / "toolkit/python/scrna_runtime.py")
RUNTIME = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(RUNTIME)


class AnnotationApprovalTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.obj = self.root / "object.qs"
        self.obj.write_text("prepared object")
        self.table = self.root / "decisions.tsv"
        self.table.write_text("cluster\tbroad\tfine\tdecision\n0\tEndothelial\tEndothelial\tconfirmed\n")
        self.record = self.root / "review.json"
        self.record.write_text(json.dumps({"schema_version": 1, "kind": "annotation_review",
            "review_run_id": "review-1", "cluster_column": "cluster", "cluster_ids": ["0"],
            "apply_object": {"sha256": RUNTIME.sha256(self.obj)}}))
        self.config = {"project": {"id": "approval-fixture"}, "workflow": {"action": "apply_confirmed"},
            "input": {"object": str(self.obj), "decisions": str(self.table), "review_record": str(self.record)},
            "metadata": {"sample": "sample", "cluster": "cluster", "reduction": "umap"},
            "annotation": {"broad_column": "broad", "fine_column": "fine"},
            "approval": {"status": "approved", "source": "human", "review_run_id": "review-1",
                         "decision_sha256": RUNTIME.sha256(self.table), "review_record_sha256": RUNTIME.sha256(self.record),
                         "provenance": "synthetic_automated_test_not_human_approval"},
            "output_dir": str(self.root / "results")}

    def errors(self, config=None):
        return RUNTIME.validate("08-scrna-annotate-cells", config or self.config, self.root / "config.json")[0]

    def test_bound_declaration_passes_and_missing_approval_fails(self):
        self.assertEqual(self.errors(), [])
        config = copy.deepcopy(self.config)
        del config["approval"]
        self.assertTrue(self.errors(config))
        config = copy.deepcopy(self.config)
        config["approval"]["source"] = "agent"
        self.assertTrue(any("source" in error for error in self.errors(config)))

    def test_decision_and_input_object_changes_invalidate_approval(self):
        self.table.write_text(self.table.read_text().replace("Endothelial", "Fibroblast"))
        self.assertTrue(any("decisions SHA256" in error for error in self.errors()))
        self.config["approval"]["decision_sha256"] = RUNTIME.sha256(self.table)
        self.assertEqual(self.errors(), [])
        self.obj.write_text("different object")
        self.assertTrue(any("input.object SHA256" in error for error in self.errors()))

    def test_review_content_and_run_identity_are_bound(self):
        config = copy.deepcopy(self.config)
        config["approval"]["review_run_id"] = "different-review"
        self.assertTrue(any("review_run_id" in error for error in self.errors(config)))
        content = json.loads(self.record.read_text())
        content["note"] = "edited after approval"
        self.record.write_text(json.dumps(content))
        self.assertTrue(any("review record SHA256" in error for error in self.errors()))

    def test_cluster_coverage_and_malformed_record_are_rejected(self):
        self.table.write_text(self.table.read_text() + "1\tTcell\tTcell\tconfirmed\n")
        self.config["approval"]["decision_sha256"] = RUNTIME.sha256(self.table)
        self.assertTrue(any("reviewed clusters" in error for error in self.errors()))
        self.record.write_text("[]")
        self.config["approval"]["review_record_sha256"] = RUNTIME.sha256(self.record)
        self.assertTrue(any("schema" in error for error in self.errors()))


if __name__ == "__main__":
    unittest.main()
