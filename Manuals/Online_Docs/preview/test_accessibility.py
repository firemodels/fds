import unittest
from lxml import etree, html

from accessibility import adapt_content, figure_labels, heading_levels, table_semantics
from audit_accessibility import inspect_document
from translate_preview import translation_units


class AccessibilityTests(unittest.TestCase):
    def test_headings_keep_anchors_and_translation_inputs(self):
        doc = html.document_fromstring('<article><h1 id="c">Chapter</h1>'
            '<h3 id="s">Section</h3><h4 id="ss">Subsection</h4></article>')
        before = [(u['id'], u['sha256']) for u in translation_units(doc)]
        heading_levels(doc)
        self.assertEqual(doc.xpath('//article/*/@id'), ['c', 's', 'ss'])
        self.assertEqual([n.tag for n in doc.xpath('//article/*')], ['h1', 'h2', 'h3'])
        self.assertEqual(before, [(u['id'], u['sha256']) for u in translation_units(doc)])

    def test_nested_caption_preserves_structured_math(self):
        doc = html.fromstring('<article><figure><div><img src="plot.svg">'
            '<figcaption>Figure 2.1: Flux <math><mfrac><mi>a</mi><mi>b</mi>'
            '</mfrac></math>.</figcaption></div></figure></article>')
        before = etree.tostring(doc.xpath('.//math')[0])
        figure_labels(doc)
        image = doc.xpath('.//img')[0]
        self.assertEqual(image.get('alt'), 'Figure 2.1')
        self.assertEqual(image.get('aria-describedby'), doc.xpath('.//figcaption')[0].get('id'))
        self.assertEqual(before, etree.tostring(doc.xpath('.//math')[0]))
        figure_labels(doc, 'es')  # Localization reuses the same caption anchor.
        self.assertEqual(image.get('aria-describedby'), 'a11y-caption-1')

    def test_grouped_headers_explicitly_cover_both_bounds(self):
        doc = html.fromstring('<article><table class="tabular" id="bands">'
            '<tr class="hline" id="rule"><td><hr></td></tr>'
            '<tr><td>Band #</td><td colspan="2">Bounds</td><td>Method</td></tr>'
            '<tr><td>1</td><td>500</td><td>880</td><td>modeled</td></tr>'
            '</table></article>')
        table_semantics(doc)
        self.assertEqual(doc.xpath('.//th/@scope'), ['col', 'colgroup', 'col'])
        self.assertEqual(doc.xpath('.//tr[last()]/td/@headers'),
                         ['bands-header-1', 'bands-header-2', 'bands-header-2', 'bands-header-3'])
        self.assertEqual(doc.xpath('.//*[@id="rule"]/@aria-hidden'), ['true'])

    def test_unknown_table_does_not_get_invented_headers(self):
        doc = html.fromstring('<article><table class="tabular"><tr><td>10</td><td>20</td></tr></table></article>')
        table_semantics(doc)
        self.assertFalse(doc.xpath('.//th'))
        self.assertEqual(doc.xpath('.//table/@data-table-review'), ['header-associations'])

    def test_layout_math_keeps_structure_and_single_scroll_region(self):
        doc = html.fromstring('<article><table class="equation"><tr><td>'
            '<math display="block"><mover><mi>q</mi><mo>˙</mo></mover></math>'
            '</td><td>(2.9)</td></tr></table><math display="block"><mi>x</mi></math></article>')
        before = [etree.tostring(n) for n in doc.xpath('.//math')]
        adapt_content(doc)
        self.assertEqual(before, [etree.tostring(n) for n in doc.xpath('.//math')])
        self.assertEqual(doc.xpath('.//table/@role'), ['presentation'])
        self.assertEqual(len(doc.xpath('.//div[@class="table-scroll"]')), 2)
        self.assertFalse(doc.xpath('.//div[@class="table-scroll"]//div[@class="table-scroll"]'))

    def test_audit_reports_real_errors_and_separates_human_review(self):
        doc = etree.ElementTree(html.document_fromstring('<html lang="en"><head><title>T</title></head>'
            '<body><a class="skip-link" href="#main-content">Skip</a><main id="main-content"><article>'
            '<h1>C</h1><h3>S</h3><img src="a.svg" alt="Technical diagram" aria-describedby="missing">'
            '<table class="tabular"><tr><td>Value</td></tr></table></article></main></body></html>'))
        errors, review = inspect_document(doc, 'page.html')
        self.assertEqual({e['rule'] for e in errors}, {'heading-order', 'image-alternative', 'label-target'})
        self.assertEqual({r['kind'] for r in review}, {'figure-description', 'table-headers'})


if __name__ == '__main__':
    unittest.main()
