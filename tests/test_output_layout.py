#!/usr/bin/env python3
"""Exercise technical output placement without running a biological analysis."""
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("runtime", ROOT / "toolkit/python/scrna_runtime.py")
RUNTIME = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(RUNTIME)


class OutputLayoutTests(unittest.TestCase):
    def run_case(self, directory, exit_code, execute=True, explicit_manifest=False):
        root = Path(directory)
        output = root / "results"
        config = root / "config.json"
        config.write_text(json.dumps({
            "project": {"id": "layout_test"}, "output_dir": str(output),
            "executor": {"argv": [sys.executable, "-c", f"print('executor diagnostic'); raise SystemExit({exit_code})"]},
        }))
        argv = ["run.py", "--config", str(config)]
        if execute:
            argv.append("--execute")
        manifest = root / "explicit.json" if explicit_manifest else root / "_provenance/07-scrna-find-cluster-markers.manifest.json"
        if explicit_manifest:
            argv.extend(["--manifest", str(manifest)])
        # Input/scientific validation is covered by analysis-contract and E2E tests.
        with patch.object(sys, "argv", argv), patch.object(RUNTIME, "validate", return_value=([], [])):
            result = RUNTIME.main("07-scrna-find-cluster-markers")
        self.assertTrue(manifest.is_file())
        self.assertFalse((root / "07-scrna-find-cluster-markers.manifest.json").exists())
        if execute:
            self.assertEqual(result, exit_code)
            self.assertFalse((output / "run.log").exists())
            log = (output / "_provenance/run.log").read_text()
            self.assertIn("executor diagnostic", log)
            self.assertIn(f"EXIT {exit_code}", log)
        else:
            self.assertEqual(result, 0)
            self.assertFalse(output.exists())

    def test_success_and_failure_logs(self):
        for code in (0, 7):
            with self.subTest(code=code), tempfile.TemporaryDirectory() as directory:
                self.run_case(directory, code)

    def test_dry_run_does_not_create_results(self):
        with tempfile.TemporaryDirectory() as directory:
            self.run_case(directory, 0, execute=False)

    def test_explicit_manifest_path_is_respected(self):
        with tempfile.TemporaryDirectory() as directory:
            self.run_case(directory, 0, explicit_manifest=True)

    def test_qc_runners_preserve_failure_records(self):
        for skill in ("02-scrna-calculate-qc-metrics", "03-scrna-review-qc"):
            with self.subTest(skill=skill), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                (root / "pixi.toml").touch()
                (root / "input.rds").touch()
                output = root / "results"
                config = root / "config.json"
                config.write_text(json.dumps({
                    "pixi": {"project": str(root), "executable": sys.executable},
                    "input": {"object": str(root / "input.rds")},
                    "metadata": {"sample": "sample"}, "output_dir": str(output),
                }))
                spec = importlib.util.spec_from_file_location("qc_runner", ROOT / "skills" / skill / "scripts/run.py")
                module = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(module)
                # Python deliberately cannot execute a pixi command. No R/data processing occurs.
                with patch.object(sys, "argv", ["run.py", "--config", str(config), "--execute"]):
                    with self.assertRaisesRegex(SystemExit, "_provenance/run.log"):
                        module.main()
                record = json.loads((output / "_provenance/run_manifest.json").read_text())
                self.assertNotEqual(record["exit_status"], 0)
                self.assertTrue((output / "_provenance/run.log").read_text())
                self.assertEqual([p.name for p in output.iterdir()], ["_provenance"])


if __name__ == "__main__":
    unittest.main()
