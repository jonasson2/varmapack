"""Tests for standalone exports of the shared mathematical documentation."""

import importlib.util
import pathlib
import unittest


SCRIPT = pathlib.Path(__file__).with_name("process-markdown.py")
SPEC = importlib.util.spec_from_file_location("process_markdown", SCRIPT)
markdown = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(markdown)
ROOT = SCRIPT.parent.parent.resolve()
MATH = ROOT / "python/docs/math.md"
DRAFT_LINK = "[Draft](math-v0.2-draft.md)"
ONLINE_LINK = f"[Draft]({markdown.DOCS_URL}math-v0.2-draft.html)"


class MarkdownExportTests(unittest.TestCase):
    def test_sibling_document(self):
        self.assertEqual(markdown.export_doc_links(DRAFT_LINK, MATH, ROOT),
                         ONLINE_LINK)

    def test_fragment_and_query(self):
        link = "[Draft](math-v0.2-draft.md?view=full#recovering-innovation-covariance)"
        expected = ONLINE_LINK[:-1] + "?view=full#recovering-innovation-covariance)"
        self.assertEqual(markdown.export_doc_links(link, MATH, ROOT), expected)

    def test_unrelated_links_unchanged(self):
        links = [
            "[External](https://example.org/note.md)",
            "[Missing](not-a-document.md)",
            "[Section](#autocovariances)",
            "[README](../../README.md)",
            "![Image](math-v0.2-draft.md)",
        ]
        for link in links:
            with self.subTest(link=link):
                self.assertEqual(markdown.export_doc_links(link, MATH, ROOT), link)

    def test_include_resolves_relative_to_source(self):
        text = markdown.read_document(ROOT / "README.in", ROOT, [], 0, 1, False)
        self.assertIn(markdown.DOCS_URL + "math-v0.2-draft.html", text)
        self.assertNotIn("](math-v0.2-draft.md)", text)

    def test_vignette_export(self):
        text = markdown.read_document(MATH, ROOT, [], 0, 0, True)
        self.assertIn(markdown.DOCS_URL + "math-v0.2-draft.html", text)
        self.assertNotIn("# Mathematical Description", text)


if __name__ == "__main__":
    unittest.main()
