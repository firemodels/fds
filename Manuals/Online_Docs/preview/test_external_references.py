"""Cross-manual PDF references and figure-cache regression checks."""
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import full_guide as full
from lxml import html
from manuals import MANUALS


class ExternalReferenceTests(unittest.TestCase):
    def test_pdf_destination_parameter_survives_tex4ht_id_filter(self):
        result = full.reader.reader_link('external-manuals/FDS_User_Guide.pdf#nameddest_section.17.1',
                                         'experiment.html', {}, {})
        self.assertEqual(result, 'external-manuals/FDS_User_Guide.pdf#nameddest=section.17.1')

    def test_omitted_html_label_uses_pdf_number_and_destination(self):
        link = html.fromstring('<a href="#missing">A.1<!-- tex4ht:ref: conservation --></a>')
        labels = {'conservation': {'number': 'C', 'destination': 'appendix.C'}}
        result = full.reader.missing_reference_link(link, '#missing', labels, {}, {})
        self.assertEqual(result, 'full-guide.pdf#nameddest=appendix.C')
        self.assertEqual(link.text_content(), 'C')
        linked = '<a href="#existing">7.2.9<!-- tex4ht:ref: conservation --></a>'
        link = html.fromstring(linked)
        self.assertEqual(full.reader.missing_reference_link(link, '#existing', labels,
                         {'existing': ['chapter.html']}, {}), '#existing')
        self.assertEqual(link.text_content(), '7.2.9')
        self.assertEqual(full.reader.missing_reference_link(html.fromstring('<a>broken</a>'),
                         '#missing', labels, {}, {}), '#missing')

    def inputs(self, root):
        guide = root / 'FDS_User_Guide'
        guide.mkdir()
        (guide / 'FDS_User_Guide.pdf').write_bytes(b'PDF fixture')
        (guide / 'FDS_User_Guide.aux').write_text(
            r'\newlabel{vegetation_model}{{17.1}{245}{Thermal Degradation}{section.17.1}{}}' + '\n')
        return {'sources': {'main.tex': r'\externaldocument[UG-]{FDS_User_Guide}'
                            r'Section \ref{UG-vegetation_model}; page \pageref{UG-vegetation_model}.'},
                'missing': []}

    def test_prefixed_labels_become_pdf_links_with_original_numbers(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            inv = self.inputs(root)
            with patch.object(full.reader, 'MANUALS', root):
                imported = full.external_references(inv)
            self.assertFalse(inv['missing'])
            source = inv['sources']['main.tex']
            output = full.resolve_external_references(source, imported['external_references'])
            self.assertNotIn(r'\externaldocument', output)
            self.assertIn(r'\href{external-manuals/FDS_User_Guide.pdf\#nameddest=section.17.1}{17.1}', output)
            self.assertIn(r'\href{external-manuals/FDS_User_Guide.pdf\#nameddest=page.245}{245}', output)
            self.assertEqual(full.resolve_external_references(r'\ref{local}', imported['external_references']), r'\ref{local}')
            self.assertEqual(len(imported['external_files']), 2)

    def test_unprefixed_import_does_not_override_local_label(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            inv = self.inputs(root)
            inv['sources']['main.tex'] = (r'\externaldocument{FDS_User_Guide}'
                                        r'\label{vegetation_model}\ref{vegetation_model}')
            with patch.object(full.reader, 'MANUALS', root):
                imported = full.external_references(inv)
            self.assertFalse(imported['external_references'])

    def test_missing_external_label_is_reported_before_conversion(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            inv = self.inputs(root)
            inv['sources']['main.tex'] += r'\ref{UG-missing}'
            with patch.object(full.reader, 'MANUALS', root):
                full.external_references(inv)
            self.assertIn('External reference label: UG-missing', inv['missing'])
            (root / 'FDS_User_Guide/FDS_User_Guide.aux').unlink()
            inv['missing'].clear()
            with patch.object(full.reader, 'MANUALS', root):
                full.external_references(inv)
            self.assertIn('External reference inputs:', inv['missing'][0])

    def test_cached_figures_avoid_poppler_and_keep_svg_mapping(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            manuals = root / 'Manuals'
            for folder in ['Bibliography', 'LaTeX_Style_Files']:
                (manuals / folder).mkdir(parents=True)
            (manuals / 'Bibliography/commoncommands.tex').write_text('')
            pdf = root / 'figure.pdf'
            pdf.write_bytes(b'PDF fixture')
            cache = root / 'output/figure-cache-v1'
            cache.mkdir(parents=True)
            stem = 'figure-' + full.reader.sha(pdf)[:12]
            for ext in ['.png', '.svg']:
                (cache / (stem + ext)).write_text('cached ' + ext)
            inv = {'manual': MANUALS['config'], 'sources': {MANUALS['config'].main: r'\includegraphics{FIGURES/figure}'},
                   'figures': {'FIGURES/figure': pdf}, 'plot_graphics': set()}
            with patch.object(full.reader, 'MANUALS', manuals), patch.object(full.reader, 'run') as run:
                tex = full.stage(inv, root / 'output/work')
                run.assert_not_called()
            self.assertIn(stem + '.png', (tex / 'preview.tex').read_text())
            self.assertIn(stem + '.svg', (tex / 'web-figures.json').read_text())
            self.assertEqual((tex / 'WEB_FIGURES' / (stem + '.svg')).read_text(), 'cached .svg')


if __name__ == '__main__':
    unittest.main()
