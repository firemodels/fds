"""Catalog publication checks without TeX, translation, or a browser."""
import json
from pathlib import Path
import tempfile
import unittest
from lxml import html
from build_catalog import catalog_entries, refresh_catalog


class CatalogTests(unittest.TestCase):
    def edition(self, root, output, key, language='en', complete=True):
        folder = root / output / 'html' / language
        folder.mkdir(parents=True)
        (folder / 'build-manifest.json').write_text(json.dumps({
            'manual': key, 'language': language, 'complete': complete,
            'pages': ['chapter.html']}))
        (folder / 'chapter.html').write_text('<html lang="en"><head><title>Chapter</title></head><body>Chapter</body></html>')
        (folder / 'index.html').write_text('Chapter')
        return folder

    def test_only_complete_correct_editions_with_present_pages(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            self.edition(root, 'tech-guide', 'tech')
            spanish = self.edition(root, 'tech-guide', 'tech', 'es')
            (spanish / 'chapter.html').unlink()
            self.edition(root, 'user-guide', 'user', complete=False)
            self.edition(root, 'config-guide', 'user')
            entries = catalog_entries(root)
            self.assertEqual([entry['key'] for entry in entries], ['tech'])
            self.assertEqual([e['code'] for e in entries[0]['editions']], ['en'])

    def test_generated_home_links_resolve_and_future_builds_join(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            self.edition(root, 'user-guide', 'user')
            refresh_catalog(root)
            doc = html.parse(str(root / 'index.html'))
            self.assertEqual(len(doc.xpath('//main/h1')), 1)
            self.assertEqual(doc.xpath('//section/h3/a/text()'), ['User Guide'])
            self.assertEqual(doc.xpath('//select[@id="catalog-language"]/option/@value'), ['en'])
            self.assertFalse(doc.xpath('//section/nav'))
            for href in doc.xpath('//section//a/@href'):
                self.assertTrue((root / href).is_file())
            self.edition(root, 'validation-guide', 'validation')
            refresh_catalog(root)
            data = json.loads((root / 'manuals.json').read_text())
            self.assertEqual([m['key'] for m in data['manuals']], ['validation', 'user'])
            self.assertTrue((root / 'catalog.js').is_file())

    def test_unrelated_index_is_preserved(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / 'index.html').write_text('Unrelated document')
            with self.assertRaisesRegex(RuntimeError, 'unrelated'):
                refresh_catalog(root)
            self.assertEqual((root / 'index.html').read_text(), 'Unrelated document')


if __name__ == '__main__':
    unittest.main()
