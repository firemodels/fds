"""Full-guide orchestration tests; no TeX conversion or real site publication."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from lxml import html
import full_guide as full
from translate_preview import translation_units


class FullGuideTests(unittest.TestCase):
    def test_source_inventory_includes_full_guide_and_external_figures(self):
        inv = full.inventory()
        self.assertEqual(inv['chapters'][0], 'Introduction_Chapter')
        self.assertEqual(inv['chapters'][-1], 'Wildland_Chapter')
        self.assertIn('Complex_Geometry_Chapter', inv['chapters'])
        self.assertEqual(inv['appendices'], 18)
        self.assertFalse(inv['missing'])
        self.assertTrue(any('Geometry_Figures' in str(p) for p in inv['figures'].values()))
        self.assertIn(r'\frontmatter', inv['sources'][full.MAIN])
        self.assertNotIn(r'\headerA', inv['sources'][full.MAIN])

    def test_geometry_disabled_and_comment_removal(self):
        text = full.uncomment('% hidden\n\\ifcompgeom yes \\else no \\fi 20\\%')
        self.assertNotIn('yes', full.active_geometry(text, False))
        self.assertIn(r'20\%', text)

    def test_comments_do_not_create_paragraphs_inside_align(self):
        source = '\\begin{align}\nx &= y \\\\\n% disabled row\n% another row\n&= z\n\\end{align}'
        clean = full.uncomment(source)
        self.assertNotIn('\n\n', clean)
        self.assertIn('\\\\\n&= z', clean)
        self.assertEqual(full.uncomment('a% comment\nb'), 'ab')

    def test_pgfplots_graphics_are_in_inventory(self):
        inv = full.inventory()
        for name in ['transport_vs_mixing_0461_a', 'transport_vs_mixing_0461_b']:
            self.assertIn('FIGURES/' + name, inv['figures'])
            self.assertIn('FIGURES/' + name, inv['plot_graphics'])

    def test_group_font_and_text_scripts_with_nested_braces(self):
        source = r'y_\mathrm{CO} + G_\mathrm{D}^\circ + x_\textrm{Reactant {dependency}} + z_{\mathrm{ok}}'
        expected = r'y_{\mathrm{CO}} + G_{\mathrm{D}}^\circ + x_{\textrm{Reactant {dependency}}} + z_{\mathrm{ok}}'
        self.assertEqual(full.reader.group_math_scripts(source), expected)
        self.assertEqual(full.reader.group_math_scripts(expected), expected)
        nested = r'x_\text{at $y_\mathrm{c}$}'
        self.assertEqual(full.reader.group_math_scripts(nested), r'x_{\text{at $y_{\mathrm{c}}$}}')

    def test_nested_primes_preserve_species_and_heat_release_scripts(self):
        source = r"\dot{q}_{\alpha'}'''=h_{\alpha'}(T) + x_{i'}''"
        expected = r"\dot{q}_{\alpha'}^{\prime\prime\prime}=h_{\alpha'}(T) + x_{i'}^{\prime\prime}"
        normalize = full.reader.group_nested_primes
        self.assertEqual(normalize(source), expected)
        self.assertEqual(normalize(expected), expected)
        unchanged = r"The user's text, q''' + q_{a}''' + x_{i'}'^n + x_{i'}' ^{n}"
        self.assertEqual(normalize(unchanged), unchanged)
        self.assertEqual(normalize(r"x_{a_{i'}}''"), r"x_{a_{i'}}^{\prime\prime}")

    def test_discovery_uses_rendered_order_not_filenames(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            fixtures = [('previewch9.html', 7, 'Chapter 1', 'Introduction'),
                        ('previewap1.html', 12, 'Appendix A', 'Nomenclature'),
                        ('previewli3.html', 2, '', 'Preface')]
            for name, order, marker, title in fixtures:
                (root / name).write_text(f'<html><body><h2><span class="titlemark">{marker}</span>'
                    f'<a id="x{order}-1000"></a>{title}</h2></body></html>')
            (root / 'preview.html').write_text('<html><body><div class="tableofcontents">'
                '<a href="previewch9.html" id="QQ2-7-1">Introduction</a></div></body></html>')
            with patch.object(full, 'expected_pages', return_value=['preface', 'introduction', 'appendix-a']):
                pages = full.discover_pages(root, {'chapters': ['Introduction_Chapter']})
            self.assertEqual([p['file'] for p in pages], ['manual-index.html', 'preface.html', 'introduction.html', 'appendix-a.html'])
            self.assertEqual(pages[-1]['number'], 'A')

    def test_discovery_reads_utf8_without_charset_metadata(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / 'previewap9.html'
            path.write_text('<html><body><h2><span class="titlemark">Appendix\u00a0I</span>'
                            '<a id="x32-288000I"></a>Lilly’s Analysis</h2></body></html>', encoding='utf-8')
            with patch.object(full, 'expected_pages', return_value=['appendix-i']):
                pages = full.discover_pages(Path(temp), {'chapters': []})
            self.assertEqual(pages[0]['file'], 'appendix-i.html')

    def test_citations_use_actual_anchor_destination(self):
        mapping = {'previewap1.html': 'appendix-a.html', 'previewli1.html': 'references.html'}
        ids = {'X0-Pope_2000': ['references.html']}
        link = full.reader.reader_link
        self.assertEqual(link('previewap1.html#X0-Pope_2000', 'overview.html', mapping, ids),
                         'references.html#X0-Pope_2000')
        self.assertEqual(link('#X0-Pope_2000', 'references.html', mapping, ids), '#X0-Pope_2000')
        external = 'https://example.org/previewap1.html#X0-Pope_2000'
        self.assertEqual(link(external, 'overview.html', mapping, ids), external)

    def edition(self, root, code, complete=True, include_page=True):
        folder = root / code
        folder.mkdir(parents=True)
        (folder / 'build-manifest.json').write_text(json.dumps({'language': code, 'complete': complete,
            'fds_revision': 'abc123', 'pages': ['overview.html']}))
        if include_page:
            (folder / 'overview.html').write_text('complete page')
        return folder

    def test_registry_lists_only_complete_present_editions_and_refreshes(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            self.edition(root, 'en')
            self.edition(root, 'es', include_page=False)
            full.refresh_registry(root)
            data = json.loads((root / 'languages.json').read_text())
            self.assertEqual([d['code'] for d in data['languages']], ['en'])
            (root / 'es/overview.html').write_text('página completa')
            full.refresh_registry(root)
            data = json.loads((root / 'languages.json').read_text())
            self.assertEqual([d['code'] for d in data['languages']], ['en', 'es'])
            (root / 'en/overview.html').unlink()
            full.refresh_registry(root)
            self.assertEqual([d['code'] for d in json.loads((root / 'languages.json').read_text())['languages']], ['es'])

    def test_failed_publish_restores_previous_edition(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            self.edition(root, 'en')
            with self.assertRaises(FileNotFoundError):
                full.publish(root / 'missing-candidate', root, 'en')
            self.assertEqual((root / 'en/overview.html').read_text(), 'complete page')
            self.assertFalse((root / '.previous-en').exists())
            self.assertFalse((root / 'languages.json').exists())

    def test_full_translation_includes_table_prose_not_math(self):
        doc = html.document_fromstring('<article><table><tr><td>Specific heat</td>'
                             '<td><math><mi>c</mi></math></td></tr></table></article>')
        units = translation_units(doc, include_tables=True)
        self.assertEqual([u['source'] for u in units], ['Specific heat'])


if __name__ == '__main__':
    unittest.main()
