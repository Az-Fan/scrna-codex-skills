import copy
import importlib.util
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("runtime", ROOT / "toolkit/python/scrna_runtime.py")
RUNTIME = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(RUNTIME)


class DEPreflightTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        source = self.root / "input.rds"
        source.write_text("fixture")
        self.config = {"project": {"id": "test"}, "input": {"object": str(source)},
                       "metadata": {"sample": "sample", "condition": "condition", "covariates": ["sex", "age"]},
                       "analysis": {"method": "pseudobulk_deseq2", "counts_source": {"kind": "raw_umi"}},
                       "comparison": {"numerator": "B", "denominator": "A"},
                       "output_dir": str(self.root / "results")}

    def errors(self, config):
        return RUNTIME.validate("11-scrna-run-differential-analysis", config, self.root / "config.json")[0]

    def test_numeric_boundaries_and_scalar_types(self):
        malformed = copy.deepcopy(self.config)
        malformed["analysis"] = None
        self.assertIn("analysis must be an object", self.errors(malformed))
        for field, minimum in {"min_samples_per_group": 2, "min_cells_per_sample_population": 1,
                               "min_total_count": 0, "min_count_per_sample": 0, "min_samples_expressed": 1}.items():
            for bad in (minimum - 1, 1.5, True, "2", None, [], float("nan"), float("inf")):
                with self.subTest(field=field, bad=bad):
                    config = copy.deepcopy(self.config)
                    config["analysis"][field] = bad
                    self.assertTrue(any(field in error for error in self.errors(config)))
            config = copy.deepcopy(self.config)
            config["analysis"][field] = minimum
            self.assertFalse(any(field in error for error in self.errors(config)))
        for field in ("lfc_threshold", "padj_threshold"):
            for bad in (-1, True, "0.05", None, [], float("nan"), float("inf")):
                config = copy.deepcopy(self.config)
                config["analysis"][field] = bad
                self.assertTrue(any(field in error for error in self.errors(config)))

    def test_unsupported_designs_and_quoted_column_names(self):
        for design in ("~ sex * condition", "~ sex + condition + sex:condition", "~ (sex + condition)^2",
                       "~ 0 + condition", "~ condition - 1", "~ log(age) + condition", "~ .", "~ sex"):
            config = copy.deepcopy(self.config)
            config["analysis"]["design"] = design
            self.assertTrue(any("design" in error for error in self.errors(config)), design)
        config = copy.deepcopy(self.config)
        config["metadata"].update(condition="condition+type", covariates=["batch*id"])
        config["analysis"]["design"] = "~ `batch*id` + `condition+type`"
        self.assertFalse(any("design" in error for error in self.errors(config)))


if __name__ == "__main__":
    unittest.main()
