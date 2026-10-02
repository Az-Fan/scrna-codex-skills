#!/usr/bin/env python3
import importlib.util
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "toolkit/python"))
import scrna_runtime as runtime
import check_dependencies as checker
spec = importlib.util.spec_from_file_location("manage_environments", ROOT / "scripts/manage_environments.py")
manager = importlib.util.module_from_spec(spec)
spec.loader.exec_module(manager)


class EnvironmentManagement(unittest.TestCase):
    def test_python_and_r_choose_the_same_root_and_overrides(self):
        rscript = Path.home() / "projects/scrna_envs/02-annotation/.pixi/envs/default/bin/Rscript"
        configs = [{}, {"runtime": {"pixi_root": "/tmp/config-env"}}, {"runtime": {"pixi_root": "${TEST_ENV_ROOT}/env"}}, {"benchmark": {"python_argv_prefix": ["~/test-python", "-u"]}, "runtime": {"sccoda_python": "~/missing-python"}}, {"benchmark": {"python_argv_prefix": ["custom-python", "-u"]}, "runtime": {"sccoda_python": "/tmp/missing-python"}}]
        for config in configs:
            with mock.patch.dict(os.environ, {"SCRNA_PIXI_ROOT": "/tmp/env-root", "TEST_ENV_ROOT": "/tmp/test-root"}):
                expected = {"integration": runtime.integration_python_prefix(config), "sccoda": runtime.sccoda_python(config)}
                expression = 'source("toolkit/R/runtime.R"); cfg<-jsonlite::fromJSON(Sys.getenv("TEST_ENV_CONFIG"),simplifyVector=FALSE); cat(jsonlite::toJSON(list(integration=as.list(scrna_python_prefix(cfg)),sccoda=scrna_python_prefix(cfg,"sccoda")),auto_unbox=TRUE))'
                env = {**os.environ, "TEST_ENV_CONFIG": json.dumps(config)}
                actual = subprocess.run([str(rscript), "-e", expression], cwd=ROOT, env=env, text=True, capture_output=True, check=True)
                self.assertEqual(json.loads(actual.stdout), expected)
        self.assertIs(checker.integration_python_prefix, runtime.integration_python_prefix)

    def test_bundle_matches_runtime_profiles_and_checksums(self):
        bundle = json.loads((manager.BUNDLE / "bundle.json").read_text())
        self.assertEqual(set(bundle["profiles"]), set(runtime.ENV_PROFILES.values()))
        manager.audit(bundle, Path("/tmp/absent-env"), list(bundle["profiles"]))

    def test_deployment_preserves_existing_config_and_binary_environment(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source_root = root / "source"
            source_root.mkdir()
            source = source_root / "pixi.toml"
            source.write_text("new")
            bundle = {"profiles": {"qc": {"environments": ["default"], "files": {"pixi.toml": {"source": "pixi.toml", "sha256": manager.digest(source)}}}}}
            target = root / "target"
            old = target / "qc/pixi.toml"
            old.parent.mkdir(parents=True)
            old.write_text("previous")
            binary = target / "qc/.pixi/envs/default/bin/Rscript"
            binary.parent.mkdir(parents=True)
            binary.write_text("binary")
            with mock.patch.object(manager, "BUNDLE", source_root):
                self.assertFalse(manager.audit(bundle, target, ["qc"])["qc"]["files"]["pixi.toml"]["matches"])
                self.assertEqual(old.read_text(), "previous")
                with self.assertRaisesRegex(ValueError, "--force"):
                    manager.deploy(bundle, target, ["qc"])
                backup = manager.deploy(bundle, target, ["qc"], force=True)
                self.assertEqual((Path(backup) / "qc/pixi.toml").read_text(), "previous")
                self.assertEqual(binary.read_text(), "binary")
                copy2 = manager.shutil.copy2
                def fail_during_copy(src, dst):
                    if Path(src) == source:
                        Path(dst).write_text("partial")
                        raise OSError("injected copy failure")
                    return copy2(src, dst)
                with mock.patch.object(manager.shutil, "copy2", side_effect=fail_during_copy):
                    with self.assertRaisesRegex(OSError, "injected copy failure"):
                        manager.deploy(bundle, target, ["qc"], force=True)
                self.assertEqual(old.read_text(), "new")
                source.write_text("corrupted")
                with self.assertRaisesRegex(ValueError, "checksum mismatch"):
                    manager.deploy(bundle, target, ["qc"], force=True)
                self.assertEqual(old.read_text(), "new")

    def test_supplemental_install_rejects_modified_source_before_installing(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            helper = root / "install_supplemental.R"
            helper.write_bytes((manager.BUNDLE / helper.name).read_bytes())
            spec = {"package": "NonexistentScrnaTest", "version": "0.1", "archive": "source.tar.gz", "urls": [], "checksum_algorithm": "sha256", "checksum": "0" * 64}
            (root / "supplemental-r.json").write_text(json.dumps({"packages": [spec]}))
            (root / "source.tar.gz").write_text("corrupted")
            rscript = Path.home() / "projects/scrna_envs/07-cell-abundance/.pixi/envs/default/bin/Rscript"
            result = subprocess.run([str(rscript), str(helper), "--install", "--library=" + str(root / "library"), "--cache=" + str(root)], text=True, capture_output=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("Cached source checksum mismatch", result.stderr)
            self.assertEqual(list((root / "library").iterdir()), [])

    def test_qc_commands_freeze_lock_and_disable_installation(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            project = root / "env"
            binary = project / ".pixi/envs/default/bin/Rscript"
            binary.parent.mkdir(parents=True)
            binary.touch()
            (project / "pixi.toml").write_text("[workspace]")
            (project / "pixi.lock").write_text("locked")
            pixi = root / "pixi"
            pixi.touch()
            obj = root / "input.rds"
            obj.touch()
            config = {"pixi": {"project": str(project), "executable": str(pixi)}, "input": {"object": str(obj)}, "metadata": {"sample": "sample"}, "output_dir": str(root / "out")}
            path = root / "config.json"
            path.write_text(json.dumps(config))
            for skill in ["02-scrna-calculate-qc-metrics", "03-scrna-review-qc"]:
                command = [sys.executable, str(ROOT / "skills" / skill / "scripts/run.py"), "--config", str(path)]
                plan = json.loads(subprocess.run(command, text=True, capture_output=True, check=True).stdout)
                self.assertIn("--frozen", plan["command"])
                self.assertIn("--no-install", plan["command"])
                self.assertEqual(plan["pixi_lock_sha256"], manager.digest(project / "pixi.lock"))
                binary.rename(binary.with_suffix(".backup"))
                failed = subprocess.run(command, text=True, capture_output=True)
                self.assertNotEqual(failed.returncode, 0)
                self.assertIn("deploy the environment", failed.stderr)
                binary.with_suffix(".backup").rename(binary)


if __name__ == "__main__":
    unittest.main()
