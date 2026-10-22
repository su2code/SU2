"""Regression checks for build isolation, links and searchable output fields."""

import json
from pathlib import Path
import tempfile
import sys
import unittest
from unittest.mock import patch

import build_docs
import check_site
import rebase_links


class DocsChecksTests(unittest.TestCase):
    def test_links_under_baseurl_with_escaped_anchor_and_relative_asset(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            page = root / "docs_v7/example/index.html"
            page.parent.mkdir(parents=True)
            (root / "asset.png").touch()
            page.write_text('<h2 id="a&amp;b">Title</h2><a href="#a%26b">ok</a>'
                            '<img src="../../asset.png"><a href="/SU2/missing/">broken</a>'
                            '<a href="#missing">broken</a><a href="https://example.test/">external</a>')
            problems, checked = check_site.check_links(root, "/SU2")
            self.assertEqual(checked, 4)
            self.assertEqual(len(problems), 2)
            self.assertTrue(any("missing anchor" in p for p in problems))
            self.assertTrue(any("missing file" in p for p in problems))

    def test_search_detects_missing_page_and_output_field(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "source"
            (source / "_docs_v7").mkdir(parents=True)
            (source / "_docs_v7/example.md").write_text(
                "---\ntitle: Example\npermalink: /docs_v7/example/\n---\nContent")
            (root / "search.json").write_text(json.dumps([]))
            with patch.object(check_site.build_site, "read_outputs", return_value={
                    "Flow": [{"name": "RMS_TEST", "literal": True}]}):
                problems, count = check_site.check_search(root, source, root, "/SU2")
            self.assertEqual(count, 0)
            self.assertEqual(len(problems), 2)

    def test_source_checks_reject_new_unknown_names_and_undocumented_options(self):
        with patch.object(check_site.build_site, "read_options", return_value=[
                {"name": "NEW_OPTION", "comment": None, "kind": "Double"}]), \
             patch.object(check_site.build_site, "read_maps", return_value={}), \
             patch.object(check_site.build_site, "read_code_names", return_value=set()), \
             patch.object(check_site.build_site, "check_docs", return_value=[("page", 1, "TYPO", "")]), \
             patch.object(check_site.build_site, "check_template", return_value=[(1, "OLD"), (2, "NEW")]):
            problems = check_site.check_source(Path("."), Path("."), {
                "template_options": ["OLD"], "missing_descriptions": [], "enum_values": []})
        self.assertEqual(len(problems), 3)

    def test_config_enum_values_ignore_source_code_and_check_disabled_variants(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "_docs_v7").mkdir()
            (root / "_docs_v7/example.md").write_text(
                "```\nSST_OPTIONS= V2003m\n```\n```\nSST_OPTIONS= V2003\n```\n"
                "```cpp\nSST_OPTIONS= SOMETHING_ELSE\n```\n")
            options = [{"name": "SST_OPTIONS", "kind": "EnumList", "args": ["", "", "", "SST_Map"]}]
            problems = check_site.check_enum_values(root, options, {"SST_Map": [("V2003m", "value")]})
            self.assertEqual(len(problems), 1)
            self.assertIn("unknown keyword V2003", next(iter(problems)))

    def test_rebase_links_handles_root_and_subpath_without_double_prefix(self):
        with tempfile.TemporaryDirectory() as directory:
            page = Path(directory) / "index.html"
            original = '<a href="/docs_v7/home/">docs</a><a href="/documents/example.pdf">pdf</a>'
            for base in ("", "/SU2"):
                page.write_text(original)
                with patch.object(sys, "argv", ["rebase_links.py", directory, base]):
                    self.assertEqual(rebase_links.main(), 0)
                    expected = page.read_text()
                    self.assertIn('href="' + base + '/docs_v7/home/"', expected)
                    self.assertIn('href="https://su2code.github.io/documents/example.pdf"', expected)
                    rebase_links.main()
                    self.assertEqual(page.read_text(), expected)

    def test_link_to_check_report_exists_during_its_own_validation(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            page = root / "docs_v7/example/index.html"
            page.parent.mkdir(parents=True)
            page.write_text('<a href="/SU2/docs-check.json">Check results</a>')
            with patch.object(sys, "argv", ["check_site.py", "--su2", directory,
                                          "--source", directory, "--site", directory]), \
                 patch.object(check_site, "check_source", return_value=[]), \
                 patch.object(check_site, "check_search", return_value=([], 0)):
                self.assertEqual(check_site.main(), 0)
            report = json.loads((root / "docs-check.json").read_text())
            self.assertEqual(report["new_problems"], [])
            self.assertEqual(report["references_checked"], 1)

    def test_build_refuses_sources_and_unmarked_nonempty_directory(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source, work = root / "source", root / "work"
            source.mkdir()
            work.mkdir()
            (work / "keep.txt").write_text("keep")
            for destination in (source, source / "build", root, work):
                with self.subTest(destination=destination), self.assertRaises(ValueError):
                    build_docs.prepare_directory(destination, (source,))
            self.assertEqual((work / "keep.txt").read_text(), "keep")


if __name__ == "__main__":
    unittest.main()
