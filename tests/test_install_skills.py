#!/usr/bin/env python3
"""Verify that updates preserve old installations and migrate the retired alias."""
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import install_skills as installer


class InstallationTests(unittest.TestCase):
    def test_update_preserves_custom_files_and_migrates_legacy(self):
        with tempfile.TemporaryDirectory() as temporary:
            target = Path(temporary) / "skills"
            active = target / "04-scrna-apply-qc-filter"
            active.mkdir(parents=True)
            (active / "custom.txt").write_text("custom data")
            legacy = target / "03.1-scrna-apply-qc-filter"
            legacy.mkdir()
            (legacy / "custom.txt").write_text("legacy data")
            unrelated = target / "other-skill"
            unrelated.mkdir()
            with mock.patch.object(installer, "SKILLS", (active.name,)):
                backup = installer.install_skills(target, force=True)
            self.assertEqual((backup / active.name / "custom.txt").read_text(), "custom data")
            self.assertEqual((backup / legacy.name / "custom.txt").read_text(), "legacy data")
            self.assertFalse(legacy.exists())
            self.assertTrue((active / "scripts/runtime.R").is_file())
            self.assertFalse(list(active.rglob("*.pyc")))
            self.assertTrue(unrelated.is_dir())
            self.assertNotIn(target, backup.parents)

    def test_refusal_and_build_failure_leave_old_installation_intact(self):
        with tempfile.TemporaryDirectory() as temporary:
            target = Path(temporary) / "skills"
            active = target / "04-scrna-apply-qc-filter"
            active.mkdir(parents=True)
            (active / "custom.txt").write_text("unchanged")
            with mock.patch.object(installer, "SKILLS", (active.name,)):
                with self.assertRaises(SystemExit):
                    installer.install_skills(target)
                with mock.patch.object(installer, "build_skill", side_effect=RuntimeError("build failed")):
                    with self.assertRaises(RuntimeError):
                        installer.install_skills(target, force=True)
            self.assertEqual((active / "custom.txt").read_text(), "unchanged")

    def test_canonical_source_cannot_be_an_install_target(self):
        with self.assertRaises(SystemExit):
            installer.install_skills(installer.REPO / "skills", force=True)

    def test_failed_replacement_restores_all_previous_directories(self):
        with tempfile.TemporaryDirectory() as temporary:
            target = Path(temporary) / "skills"
            names = ("04-scrna-apply-qc-filter", "05-scrna-benchmark-integration")
            for name in names + tuple(installer.LEGACY_SKILLS):
                old = target / name
                old.mkdir(parents=True)
                (old / "custom.txt").write_text(name)
            original_rename = Path.rename

            def fail_second_publish(path, destination):
                if path.parent.name.startswith(".scrna-install-") and path.name == names[1]:
                    raise OSError("replacement failed")
                return original_rename(path, destination)

            with mock.patch.object(installer, "SKILLS", names), mock.patch.object(Path, "rename", fail_second_publish):
                with self.assertRaisesRegex(OSError, "replacement failed"):
                    installer.install_skills(target, force=True)
            for name in names + tuple(installer.LEGACY_SKILLS):
                self.assertEqual((target / name / "custom.txt").read_text(), name)
                self.assertFalse((target / name / "scripts").exists())


if __name__ == "__main__":
    unittest.main()
