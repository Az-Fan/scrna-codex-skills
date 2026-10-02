#!/usr/bin/env python3
"""Validate selected-branch requirements and the observed scIB/pandas incompatibility."""
import json
import subprocess
import sys
import unittest
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'toolkit/python'))
import check_dependencies as checker


class DependencyChecks(unittest.TestCase):
    def test_every_built_skill_has_a_dependency_profile(self):
        released = json.loads((ROOT / 'release/runtime-manifest.json').read_text())['skills']
        self.assertEqual(set(released), set(checker.PROFILES))

    def test_qs_is_required_only_when_selected_for_automatic_filter_output(self):
        skill = '04-scrna-apply-qc-filter'
        required, optional = checker.package_requirements(skill, {})
        self.assertNotIn('qs', required)
        self.assertIn('qs', optional)
        for config in [{'input': {'object': 'project.qs'}}, {'output': {'object_name': 'filtered.qs'}}, {'output': {'object_format': 'qs'}}]:
            required, optional = checker.package_requirements(skill, config)
            self.assertIn('qs', required)
            self.assertNotIn('qs', optional)
        for skill in ['06-scrna-preprocess-and-cluster', '10-scrna-score-programs']:
            self.assertIn('qs', checker.package_requirements(skill, {})[0])
            self.assertNotIn('qs', checker.package_requirements(skill, {'output': {'object_format': 'rds'}})[0])

    def test_selected_method_promotes_its_dependencies(self):
        config = {'tasks': [{'method': 'ucell', 'gene_sets': {'source': 'msigdb'}}]}
        required, optional = checker.package_requirements('10-scrna-score-programs', config)
        self.assertTrue({'UCell', 'msigdbr', 'ggplot2'}.issubset(required))
        self.assertNotIn('VISION', required)
        required, _ = checker.package_requirements('13-scrna-test-cell-abundance', {'analysis': {'methods': ['propeller']}})
        self.assertIn('speckle', required)
        self.assertNotIn('sccomp', required)
        self.assertNotIn('cmdstanr', required)

    def test_present_modules_do_not_hide_the_result_table_incompatibility(self):
        payload = {'modules': {'scib_metrics': True, 'pandas': True}, 'versions': {'scib_metrics': '0.5.7', 'pandas': '3.0.6'}}
        completed = subprocess.CompletedProcess([], 0, json.dumps(payload), '')
        with mock.patch.object(checker.subprocess, 'run', return_value=completed):
            present, error, versions = checker.probe_python(['python'], list(payload['modules']))
        self.assertTrue(all(present.values()))
        self.assertIn('pandas < 3', error)
        payload['versions']['pandas'] = '2.3.3'
        completed.stdout = json.dumps(payload)
        with mock.patch.object(checker.subprocess, 'run', return_value=completed):
            _, error, _ = checker.probe_python(['python'], list(payload['modules']))
        self.assertIsNone(error)


if __name__ == '__main__':
    unittest.main()
