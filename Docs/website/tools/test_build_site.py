"""Run with python3 -m unittest discover -s Docs/website/tools."""

from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import build_site


class ReferenceExtractionTests(unittest.TestCase):
    def test_comments_are_not_enum_values(self):
        with tempfile.TemporaryDirectory() as folder:
            include = Path(folder) / "Common/include"
            include.mkdir(parents=True)
            (include / "options.hpp").write_text(
                """static const MapType<std::string, int> Example_Map = {
  MakePair("ACTIVE", 1)
  // MakePair("DISABLED_LINE", 2)
  /* MakePair("DISABLED_BLOCK", 3) */
  MakePair("https://example.test/a/*b*/", 4)
};"""
            )
            self.assertEqual(
                build_site.read_maps(folder)["Example_Map"],
                [("ACTIVE", "1"), ("https://example.test/a/*b*/", "4")],
            )

    def test_current_disabled_model_values(self):
        su2 = Path(__file__).resolve().parents[3]
        maps = build_site.read_maps(su2)
        sst = dict(maps["SST_Options_Map"])
        self.assertIn("V1994m", sst)
        self.assertIn("V2003m", sst)
        self.assertNotIn("V1994", sst)
        self.assertNotIn("V2003", sst)
        self.assertNotIn("GENERALIZED_ALPHA", dict(maps["Time_Int_Map_FEA"]))

    def test_array_storage_is_not_a_default(self):
        su2 = Path(__file__).resolve().parents[3]
        arrays = [
            o
            for o in build_site.read_options(su2)
            if o["kind"] in ("DoubleArray", "UShortArray", "ULongArray")
        ]
        self.assertTrue(arrays)
        for option in arrays:
            with self.subTest(option=option["name"]):
                self.assertEqual(build_site.describe_type(option, {})[2], "")

    def test_parser_renames_accept_single_and_or_conditions(self):
        with tempfile.TemporaryDirectory() as folder:
            source = Path(folder) / "Common/src/CConfig.cpp"
            source.parent.mkdir(parents=True)
            source.write_text(
                r"""
                if (!option_name.compare("OLD"))
                    newString.append("Use NEW instead.\n\n");
                else if ( !option_name.compare("FIRST") ||
                          !option_name.compare("SECOND") )
                    newString.append("Both names were removed.\n");
                if (!option_name.compare("BAD_FIRST") !option_name.compare("BAD_SECOND"))
                    newString.append("No OR operator.");
            """
            )
            self.assertEqual(
                build_site.read_parser_renames(folder),
                {
                    "OLD": "Use NEW instead.",
                    "FIRST": "Both names were removed.",
                    "SECOND": "Both names were removed.",
                },
            )

    def test_parser_renames_reject_long_malformed_conditions_promptly(self):
        with tempfile.TemporaryDirectory() as folder:
            source = Path(folder) / "Common/src/CConfig.cpp"
            source.parent.mkdir(parents=True)
            source.write_text("if(" + ' !option_name.compare("OLD")' * 48)
            result = subprocess.run(
                [
                    sys.executable,
                    "-c",
                    "import build_site, sys; print(build_site.read_parser_renames(sys.argv[1]))",
                    folder,
                ],
                cwd=Path(__file__).resolve().parent,
                capture_output=True,
                text=True,
                timeout=5,
                check=True,
            )
            self.assertEqual(result.stdout.strip(), "{}")


if __name__ == "__main__":
    unittest.main()
