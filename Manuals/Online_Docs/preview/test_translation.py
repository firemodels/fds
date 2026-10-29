"""Protect scientific markup while translating the surrounding prose."""
import unittest
from lxml import etree, html
from translate_preview import translation_units, apply_translation


class TranslationTests(unittest.TestCase):
    def unit(self):
        doc = html.fromstring('<main><article><p id="source">The density '
                             '<math><mi>ρ</mi></math> appears in '
                             '<a href="overview.html#eq-6">2.6</a>.</p></article></main>')
        return doc, translation_units(doc)[0]

    def test_preserves_math_and_reference(self):
        doc, unit = self.unit()
        before = etree.tostring(doc.xpath('.//math')[0])
        apply_translation(unit, 'La densidad ⟦0⟧ aparece en ⟦1⟧.')
        self.assertEqual(before.split(b'</math>')[0],
                         etree.tostring(doc.xpath('.//math')[0]).split(b'</math>')[0])
        self.assertEqual(doc.xpath('.//a/@href'), ['overview.html#eq-6'])
        self.assertEqual(doc.xpath('.//p/@id'), ['source'])

    def test_rejects_lost_or_reordered_math(self):
        for translation in ['La densidad aparece en ⟦1⟧.', '⟦1⟧ luego ⟦0⟧']:
            with self.subTest(translation=translation):
                _, unit = self.unit()
                with self.assertRaises(ValueError):
                    apply_translation(unit, translation)

    def test_changed_prose_changes_source_hash(self):
        doc, unit = self.unit()
        doc.xpath('.//p')[0].text = 'The pressure '
        self.assertNotEqual(unit['sha256'], translation_units(doc)[0]['sha256'])


if __name__ == '__main__':
    unittest.main()
