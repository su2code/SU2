"""Run with python3 -m unittest discover -s Docs/website/tools."""

from pathlib import Path
import tempfile
import unittest

import build_site


class ReferenceExtractionTests(unittest.TestCase):
    def test_comments_are_not_enum_values(self):
        with tempfile.TemporaryDirectory() as folder:
            include = Path(folder) / "Common/include"
            include.mkdir(parents=True)
            (include / "options.hpp").write_text('''static const MapType<std::string, int> Example_Map = {
  MakePair("ACTIVE", 1)
  // MakePair("DISABLED_LINE", 2)
  /* MakePair("DISABLED_BLOCK", 3) */
  MakePair("https://example.test/a/*b*/", 4)
};''')
            self.assertEqual(build_site.read_maps(folder)["Example_Map"],
                             [("ACTIVE", "1"), ("https://example.test/a/*b*/", "4")])

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
        arrays = [o for o in build_site.read_options(su2)
                  if o["kind"] in ("DoubleArray", "UShortArray", "ULongArray")]
        self.assertTrue(arrays)
        for option in arrays:
            with self.subTest(option=option["name"]):
                self.assertEqual(build_site.describe_type(option, {})[2], "")


if __name__ == "__main__":
    unittest.main()
