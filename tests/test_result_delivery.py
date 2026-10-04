#!/usr/bin/env python3
"""Verify result claims, stage inheritance, and portable navigation."""
import importlib.util
import json
import re
import shutil
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from urllib.parse import unquote

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("delivery", ROOT / "toolkit/python/result_delivery.py")
DELIVERY = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(DELIVERY)


def put(root, name, text="scientific result"):
    path = root / name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


def check_links(test, root):
    for index in root.rglob("RESULTS.md"):
        for target in re.findall(r"\]\(([^)]+)\)", index.read_text()):
            resolved = (index.parent / unquote(target)).resolve()
            test.assertTrue(resolved.exists(), (index, target))
            test.assertTrue(resolved == root.resolve() or root.resolve() in resolved.parents)


class DeliveryTests(unittest.TestCase):
    def test_current_run_excludes_stale_artifacts_and_status(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            put(root, "old.pdf")
            put(root, "task_status.tsv", "task_id\tstatus\nold\tfailed\n")
            before = DELIVERY.snapshot(root)
            put(root, "current.tsv")
            record = DELIVERY.finish(root, before, "07-markers", {}, 0)
            self.assertEqual(record["status"], "completed")
            self.assertEqual(record["current_files"], ["current.tsv"])
            index = (root / "RESULTS.md").read_text()
            self.assertNotIn("[old.pdf]", index)
            self.assertNotIn("status=failed", index)
            self.assertIn("2 unchanged files", index)
            check_links(self, root)

    def test_partial_comparison_and_relocated_links(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "results"
            put(root, "task_status.tsv", "task_id\tstatus\tlow_replication_warning\nok\tcompleted\tTRUE\nbad\tfailed\tFALSE\n")
            put(root, "comparisons/ok/all_genes.tsv")
            put(root, "comparisons/ok/figure space.png")
            put(root, "comparisons/bad/ERROR.txt", "invalid design")
            record = DELIVERY.finish(root, {}, "11-de", {"project": {"id": "test study"}}, 0)
            self.assertEqual(record["status"], "partial")
            index = (root / "RESULTS.md").read_text()
            self.assertIn("low_replication_warning", index)
            self.assertIn("comparisons/ok/RESULTS.md", index)
            self.assertNotIn("comparisons/ok/all_genes.tsv", index)
            self.assertIn("Status: **failed**", (root / "comparisons/bad/RESULTS.md").read_text())
            self.assertFalse((root / "comparisons/ok/_provenance").exists())
            moved = Path(directory) / "copied results"
            shutil.copytree(root, moved)
            check_links(self, moved)

    def test_skips_empty_and_missing_genes_are_visible_without_false_failure(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            put(root, "plot_status.tsv", "family\tstatus\treason\nDE\tskipped\tnot configured\n")
            put(root, "gene_status.tsv", "symbol\tfound_in_object\nabsent\tFALSE\n")
            put(root, "_provenance/metric_status.tsv", "metric\tstatus\treason\nY\tskipped\tmissing annotation\n")
            put(root, "task_status.tsv", "task_id\tstatus\nempty\tempty\n")
            record = DELIVERY.finish(root, {}, "15-genes", {}, 0)
            self.assertEqual(record["status"], "completed")
            index = (root / "RESULTS.md").read_text()
            for word in ("not configured", "absent", "missing annotation", "no eligible results"):
                self.assertIn(word, index)
            check_links(self, root)

    def test_awaiting_review_and_confirmed_stage_inheritance(self):
        for skill, prepare, apply in (("06-cluster", "run", "finalize_resolution"), ("08-annotation", "prepare_review", "apply_confirmed")):
            with self.subTest(skill=skill), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                put(root, "review.png")
                put(root, "threshold_review.tsv")
                if skill.startswith("06"):
                    put(root, "_provenance/workflow_state.json", json.dumps({"status": "awaiting_resolution_confirmation"}))
                config = {"workflow": {"action": prepare}}
                self.assertEqual(DELIVERY.finish(root, {}, skill, config, 0)["status"], "awaiting confirmation")
                # An unrelated old run file is not part of the registered stage.
                put(root, "unrelated.tsv")
                before = DELIVERY.snapshot(root)
                put(root, "final_object.rds")
                if skill.startswith("06"):
                    put(root, "_provenance/workflow_state.json", json.dumps({"status": "complete"}))
                record = DELIVERY.finish(root, before, skill, {"workflow": {"action": apply}}, 0)
                self.assertEqual(record["status"], "completed")
                self.assertEqual(record["retained_stage_files"], ["review.png", "threshold_review.tsv"])
                self.assertNotIn("unrelated.tsv", record["retained_stage_files"])
                self.assertIn("retained from the preceding review stage", (root / "RESULTS.md").read_text())
                check_links(self, root)

    def test_failed_executor_never_claims_success(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            record = DELIVERY.finish(root, {}, "07-markers", {}, 7)
            self.assertEqual(record["status"], "failed")
            self.assertIn("No scientific deliverables", (root / "RESULTS.md").read_text())
            check_links(self, root)

    def test_legacy_inheritance_whitelist_and_unverified_label_persist(self):
        for skill, action, approved_name in (("08-annotation", "apply_confirmed", "annotation_review.tsv"),
                                              ("06-cluster", "finalize_resolution", "standard_resolution_stability.png")):
            with self.subTest(skill=skill), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                put(root, approved_name)
                put(root, "old_umap.pdf")
                put(root, "old_marker.tsv")
                put(root, "old_annotation.tsv")
                before = DELIVERY.snapshot(root)
                put(root, "final_object.qs")
                config = {"workflow": {"action": action}}
                record = DELIVERY.finish(root, before, skill, config, 0)
                self.assertEqual(record["retained_stage_files"], [approved_name])
                self.assertEqual(record["legacy_unverified_files"], [approved_name])
                index = (root / "RESULTS.md").read_text()
                self.assertIn("stage provenance unverified", index)
                self.assertNotIn("retained from the preceding review stage", index)
                self.assertNotIn("[old_umap.pdf]", index)
                before = DELIVERY.snapshot(root)
                DELIVERY.finish(root, before, skill, config, 0)
                self.assertIn("stage provenance unverified", (root / "RESULTS.md").read_text())
                check_links(self, root)

    def test_corrupt_or_other_skill_registry_does_not_claim_arbitrary_files(self):
        for content in ("[]", "not json", json.dumps({"skill": "other", "current_files": ["old.pdf"]}),
                        json.dumps({"skill": "08-annotation", "current_files": "old.pdf"}),
                        json.dumps({"skill": "08-annotation", "current_files": ["../old.pdf"]})):
            with self.subTest(content=content), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                put(root, "old.pdf")
                put(root, "_provenance/result_delivery.json", content)
                before = DELIVERY.snapshot(root)
                record = DELIVERY.finish(root, before, "08-annotation", {"workflow": {"action": "apply_confirmed"}}, 0)
                self.assertEqual(record["retained_stage_files"], [])

    def test_unresolved_recommendation_and_compressed_decisions(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            put(root, "recommendation_status.json", json.dumps({"status": "review_required", "reasons": ["incomplete evidence"]}))
            self.assertEqual(DELIVERY.finish(root, {}, "05-integration", {}, 0)["status"], "partial")
            self.assertIn("incomplete evidence", (root / "RESULTS.md").read_text())
            self.assertEqual(DELIVERY.category("cell_filter_decisions.tsv.gz"), "Review and decisions")
            self.assertEqual(DELIVERY.category("metadata.tsv.gz"), "Complete tables and summaries")
            self.assertEqual(DELIVERY.category("counts.mtx.gz"), "Downstream objects and matrices")

    def test_delivery_failure_keeps_executor_outcome_and_fails_wrapper(self):
        with tempfile.TemporaryDirectory() as directory:
            record = {"exit_status": 0, "status": "completed"}
            with patch.object(DELIVERY, "finish", side_effect=OSError("disk full")):
                code = DELIVERY.supervise(Path(directory), {}, "07-markers", {}, 0, record)
            self.assertEqual(code, 1)
            self.assertEqual(record["executor_exit_status"], 0)
            self.assertEqual(record["delivery_status"], "failed")
            self.assertIn("disk full", record["delivery_error"])

    def test_mutable_entrypoint_does_not_claim_old_child_indexes(self):
        with tempfile.TemporaryDirectory() as directory:
            record = {}
            root = Path(directory)
            put(root, "comparisons/old/RESULTS.md", "old index")
            before = DELIVERY.snapshot(root)
            put(root, "results.tsv")
            self.assertEqual(DELIVERY.supervise(root, before, "07-markers", {}, 0, record), 0)
            self.assertEqual(record["delivery_status"], "completed")
            artifact = record["artifacts"][0]
            self.assertEqual(len(record["artifacts"]), 1)
            self.assertIsNone(artifact["sha256"])
            self.assertIsNone(artifact["bytes"])
            self.assertEqual(Path(artifact["path"]).name, "RESULTS.md")
            self.assertTrue(artifact["mutable"])


if __name__ == "__main__":
    unittest.main()
