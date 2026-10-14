"""Check the public entry point and the resources needed outside a source checkout."""

import json
from pathlib import Path
import subprocess
import sys
import unittest
import tempfile

from mqc.config import Config


class InstallationTests(unittest.TestCase):
    def test_reference_paths_are_independent_of_working_directory(self):
        with tempfile.TemporaryDirectory() as tmp:
            result = subprocess.run(
                [
                    sys.executable,
                    "-c",
                    "from mqc.config import Config; import json; "
                    'c=Config("output"); '
                    "print(json.dumps([str(c.general_library_path), str(c.general_library_path_meta), "
                    "str(c.metacyc_rules), str(c.modelseed_rxn)]))",
                ],
                cwd=tmp,
                capture_output=True,
                text=True,
                check=True,
            )
            files = json.loads(result.stdout)
            self.assertTrue(all(Path(p).is_file() for p in files))
            self.assertTrue((Path(tmp) / "output").is_dir())

    def test_cli_exposes_both_workflows_outside_checkout(self):
        with tempfile.TemporaryDirectory() as tmp:
            result = subprocess.run(
                [sys.executable, "-m", "mqc", "--help"],
                cwd=tmp,
                capture_output=True,
                text=True,
            )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("{check,correct}", result.stdout)

    def test_missing_input_fails_without_creating_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            result = subprocess.run(
                [sys.executable, "-m", "mqc", "missing.xml", "output"],
                cwd=tmp,
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("Model file does not exist", result.stderr)
            self.assertFalse((Path(tmp) / "output").exists())

    def test_all_configured_reference_files_are_installed(self):
        with tempfile.TemporaryDirectory() as tmp:
            config = Config(tmp)
            names = [
                "general_library_path",
                "general_library_path_meta",
                "metacyc_rules",
                "modelseed_met",
                "meta_met",
                "kegg_met",
                "bigg_metacyc",
                "modelseed_rxn",
                "virtual_met",
            ]
            for name in names:
                with self.subTest(resource=name):
                    self.assertTrue(Path(getattr(config, name)).is_file())


if __name__ == "__main__":
    unittest.main()
