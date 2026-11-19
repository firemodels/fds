"""Manual selection and staging regressions; no external tools required."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from types import SimpleNamespace

import build_preview as reader
import full_guide as full
from accessibility import figure_labels
from lxml import html
from manual_sources import inventory
from manuals import MANUALS


class ManualTests(unittest.TestCase):
    def test_uncaptioned_histograms_keep_source_identity_and_review_flag(self):
        body = html.fromstring('<div><img src="assets/FDS_HGL_Depth_Histogram-0123456789ab.svg" alt="Technical diagram">'
                               '<img src="unknown.svg" alt="Technical diagram"></div>')
        figure_labels(body)
        images = body.xpath('.//img')
        self.assertEqual(images[0].get('alt'), 'Plot: FDS HGL Depth Histogram')
        self.assertEqual(images[0].get('data-description-review'), 'source-label-only')
        self.assertEqual(images[1].get('alt'), 'Technical diagram')
        figure_labels(body, 'es')
        self.assertEqual(images[0].get('alt'), 'Gráfico: FDS HGL Depth Histogram')

    def test_compound_math_grouping_preserves_operands_and_literal_examples(self):
        source = r'x^\frac{2}{5} + \sqrt \frac{\tau k}{\rho c} + y_\sqrt[3]{z}'
        expected = r'x^{\frac{2}{5}} + \sqrt{\frac{\tau k}{\rho c}} + y_{\sqrt[3]{z}}'
        self.assertEqual(reader.prepare_math(source), expected)
        self.assertEqual(reader.prepare_math(expected), expected)
        nested = r'x^\sqrt\frac{a}{b}'
        self.assertEqual(reader.prepare_math(nested), r'x^{\sqrt{\frac{a}{b}}}')
        literal = r'\begin{lstlisting}x^\frac{2}{5}\end{lstlisting} \verb|y_\mathrm{CO}|'
        self.assertEqual(reader.prepare_math(literal), literal)

    def test_recursive_sources_inline_chapters_appendices_and_graphics(self):
        manual = MANUALS['user']
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            guide = root / manual.directory
            guide.mkdir()
            (root / 'Bibliography').mkdir()
            (root / 'Bibliography/commoncommands.tex').write_text('')
            (guide / 'tables').mkdir()
            (guide / 'FIGURES').mkdir()
            (guide / 'FIGURES/image.png').write_bytes(b'fixture')
            (guide / 'FIGURES/overlay.png').write_bytes(b'fixture')
            (guide / manual.main).write_text(r'''
\includeonly{chapter}
\begin{document}\frontmatter\chapter*{Preface}
\mainmatter\part{Part One}\chapter{First}\input{tables/generated}
\include{chapter}\chapter{Third}\appendix\input{appendices}
\end{document}''')
            (guide / 'tables/generated.tex').write_text(r'\incgraphics{FIGURES/image}'
                r'\begin{picture}(100,100)\put(0,0){\includegraphics{FIGURES/overlay}}\end{picture}')
            (guide / 'chapter.tex').write_text(r'\chapter{Second}\addplot graphics {FIGURES/image};')
            (guide / 'appendices.tex').write_text(r'\chapter{A}\chapter{B}')
            inv = inventory(manual, root)
            self.assertEqual(inv['chapters'], ['chapter-1', 'chapter', 'chapter-3'])
            self.assertEqual(inv['part_sources'], [manual.main])
            self.assertEqual(inv['appendix_slugs'], ['appendix-a', 'appendix-b'])
            self.assertEqual(inv['appendix_sources'], ['appendices.tex'] * 2)
            self.assertIn('tables/generated.tex', inv['sources'])
            self.assertNotIn(r'\includeonly', inv['sources'][manual.main])
            self.assertEqual(inv['plot_graphics'], {'FIGURES/image', 'FIGURES/overlay'})
            self.assertFalse(inv['missing'])
            (guide / 'tables/generated.tex').unlink()
            self.assertIn('LaTeX input: tables/generated', inventory(manual, root)['missing'])

    def test_wrong_manual_output_is_rejected_before_conversion(self):
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp)
            edition = output / 'html/en'
            edition.mkdir(parents=True)
            (edition / 'build-manifest.json').write_text(json.dumps({'manual': 'tech'}))
            args = SimpleNamespace(manual='user', output=output)
            with patch.object(full, 'inventory') as inv:
                with self.assertRaisesRegex(RuntimeError, 'another manual'):
                    full.build_full(args)
                inv.assert_not_called()
            with self.assertRaisesRegex(RuntimeError, 'another manual'):
                full.english_edition(output, MANUALS['user'])

    def test_distinct_outputs_and_entry_pages(self):
        self.assertEqual(len({m.output for m in MANUALS.values()}), 5)
        self.assertEqual(MANUALS['tech'].entry, 'overview.html')
        self.assertEqual(MANUALS['user'].entry, 'chapter-1.html')

    def test_part_landing_pages_are_discovered_and_keep_source(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / 'previewpa1.html').write_text('<html><body><h1 class="partHead" id="basics">'
                '<span class="titlemark">Part I</span><a id="x10-9000I"></a>The Basics</h1></body></html>')
            inv = {'manual': MANUALS['user'], 'chapters': [],
                   'part_sources': [MANUALS['user'].main]}
            with patch.object(full, 'expected_pages', return_value=['part-1']):
                pages = full.discover_pages(root, inv)
            self.assertEqual(pages[0]['file'], 'part-1.html')
            self.assertEqual(pages[0]['eyebrow'], 'PART I')
            self.assertEqual(pages[0]['source'], MANUALS['user'].main)

    def test_multipage_tables_preserve_captions_headers_and_cells(self):
        source = r'''\begin{xltabular}{\textwidth}{|p{6cm}|X|l|}
\caption{Quantities}\label{quantities}\\
Name & Units & Value\\\endfirsthead
Name & Units & Value\\\endhead
\rowcolor{lavender}\ct{DENSITY} & kg/m$^3$ & 1\\
\end{xltabular}'''
        staged = full.prepare_tables(source)
        self.assertIn(r'\begin{longtable}{|p{6cm}|p{0.2\linewidth}|l|}', staged)
        self.assertIn(r'\end{longtable}', staged)
        self.assertNotIn('rowcolor', staged)
        self.assertNotIn('lavender', staged)
        for token in [r'\caption{Quantities}', r'\label{quantities}',
                      r'\endfirsthead', r'\endhead', r'\ct{DENSITY} & kg/m$^3$ & 1']:
            self.assertIn(token, staged)
        self.assertEqual(full.prepare_tables(staged), staged)

    def test_math_escapes_leave_code_mode_before_superscripts(self):
        source = r'\lstinline[mathescape=true]!ADA_Z$^1$!'
        self.assertEqual(full.prepare_tables(source), r'\lstinline!ADA_Z!$^1$')
        source = r'\lstinline[mathescape=true]!A$x^2$B!'
        self.assertEqual(full.prepare_tables(source), r'\lstinline!A!$x^2$\lstinline!B!')

    def test_print_rotation_and_boxed_code_keep_contents(self):
        source = r'\begin{tabular}[\textwidth]{c|c}A & B\end{tabular}'
        self.assertEqual(full.prepare_tables(source), r'\begin{tabular}{c|c}A & B\end{tabular}')
        aligned = r'\begin{tabular}[t]{c}A\end{tabular}'
        self.assertEqual(full.prepare_tables(aligned), aligned)
        source = r'\begin{landscape}\begin{sidewaystable}[p]content\end{sidewaystable}\end{landscape}'
        self.assertEqual(full.prepare_tables(source), r'\begin{table}[p]content\end{table}')
        source = r'\underbrace{x}_{\mbox{\ct{P_0}}} + \ct{P_0}'
        self.assertEqual(full.prepare_tables(source), r'\underbrace{x}_{\mbox{\texttt{P\_0}}} + \ct{P_0}')

    def test_smokeview_object_images_use_explicit_table_labels(self):
        body = html.fromstring('''<div><table><caption>Dual frame static objects</caption>
<tr><td><code>SMOKEVIEW_ID</code></td><td colspan="2">Image</td></tr>
<tr><td></td><td>inactive</td><td>active</td></tr>
<tr><td><code>nozzle</code></td><td><img src="off.png"></td><td><img src="on.png"></td></tr>
</table><table><tr><td>Unknown diagram</td><td><img src="unknown.png"></td></tr></table></div>''')
        figure_labels(body)
        images = body.xpath('.//img')
        self.assertEqual(images[0].get('alt'), 'Smokeview object nozzle, inactive')
        self.assertEqual(images[1].get('alt'), 'Smokeview object nozzle, active')
        self.assertEqual(images[0].get('data-description-review'), 'table-label-only')
        self.assertEqual(images[0].get('aria-describedby'), body.xpath('.//caption')[0].get('id'))
        self.assertIsNone(images[2].get('alt'))
        states = body.xpath('.//tr[@data-object-states="true"]/td')
        states[1].text, states[2].text = 'inactiva', 'activa'
        figure_labels(body, 'es')
        self.assertEqual(images[0].get('alt'), 'Objeto Smokeview nozzle, inactiva')
        self.assertEqual(images[1].get('alt'), 'Objeto Smokeview nozzle, activa')


if __name__ == '__main__':
    unittest.main()
