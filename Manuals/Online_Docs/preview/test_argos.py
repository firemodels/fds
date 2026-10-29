"""Offline integration checks using a deterministic translator double."""
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from lxml import etree, html
import full_guide as full
from translate_argos import generate_cache, translate_protected, digest
from translate_preview import translation_units, apply_translation


class Translator:
    provenance = {'provider': 'test', 'model_version': '1'}

    def __init__(self):
        self.calls = 0

    def translate(self, text):
        self.calls += 1
        return text.replace('The density', 'La densidad').replace('appears', 'aparece')


class ArgosTests(unittest.TestCase):
    def test_guarded_literals_and_math_survive_placeholder_loss(self):
        def dropping_model(text):
            return text.replace('⟦0⟧', '').replace('The density', 'La densidad')
        source = 'The density ⟦0⟧ in FDS is 10 kg.'
        value, fallback = translate_protected(source, dropping_model)
        self.assertTrue(fallback)
        self.assertEqual(value, 'La densidad ⟦0⟧ in FDS is 10 kg.')

    def test_reordered_placeholders_use_span_fallback(self):
        value, fallback = translate_protected('A value ⟦0⟧ and ⟦1⟧.',
            lambda text: text.replace('⟦0⟧ and ⟦1⟧', '⟦1⟧ and ⟦0⟧'))
        self.assertTrue(fallback)
        self.assertEqual(value, 'A value ⟦0⟧ and ⟦1⟧.')

    def test_cache_reuse_changed_text_and_model_invalidation(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / 'es/overview.json'
            translator = Translator()
            units = [{'id': 'b000', 'source': 'The density', 'sha256': digest('The density')}]
            generate_cache(units, path, translator)
            self.assertEqual(translator.calls, 1)
            units[0]['id'] = 'b020'
            blocks, stats = generate_cache(units, path, translator)
            self.assertEqual((translator.calls, blocks[0]['id'], stats['cached']), (1, 'b020', 1))
            units.append({'id': 'b021', 'source': 'appears', 'sha256': digest('appears')})
            generate_cache(units, path, translator)
            self.assertEqual(translator.calls, 2)
            translator.provenance = {'provider': 'test', 'model_version': '2'}
            generate_cache(units, path, translator)
            self.assertEqual(translator.calls, 4)

    def test_translation_preserves_markup_byte_for_byte(self):
        doc = html.fromstring('<article><p id="eq">The density '
                             '<math><mi>ρ</mi></math> appears in <a href="#eq">2.6</a>.</p></article>')
        before = [etree.tostring(n, with_tail=False) for n in doc.xpath('.//math | .//a')]
        for unit in translation_units(doc):
            value, _ = translate_protected(unit['source'], Translator().translate)
            apply_translation(unit, value)
        self.assertEqual(before, [etree.tostring(n, with_tail=False) for n in doc.xpath('.//math | .//a')])

    def test_non_english_never_starts_latex(self):
        with tempfile.TemporaryDirectory() as temp:
            args = SimpleNamespace(output=Path(temp), language='es', refresh_languages=False)
            with patch.object(full, 'inventory', side_effect=AssertionError('No TeX inventory')):
                with self.assertRaisesRegex(RuntimeError, 'Build English first'):
                    full.build_full(args)

    def test_failed_translation_leaves_published_spanish_and_registry(self):
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp)
            en = output / 'html/en'
            es = output / 'html/es'
            en.mkdir(parents=True)
            es.mkdir()
            (en / 'overview.html').write_text('<article><h1>Overview</h1></article>')
            (en / 'build-manifest.json').write_text(json.dumps({
                'complete': True, 'language': 'en', 'pages': ['overview.html'], 'fds_revision': 'abc'}))
            (es / 'overview.html').write_text('last good Spanish')
            (output / 'html/languages.json').write_text('last good registry')
            args = SimpleNamespace(language='es', check=False)
            with patch('translate_argos.ArgosTranslator', return_value=Translator()):
                with patch.object(full, 'translate_edition', side_effect=RuntimeError('model failed')):
                    with self.assertRaisesRegex(RuntimeError, 'model failed'):
                        full.build_translated(args, output)
            self.assertEqual((es / 'overview.html').read_text(), 'last good Spanish')
            self.assertEqual((output / 'html/languages.json').read_text(), 'last good registry')
            self.assertFalse(list(output.glob('edition-*')))

    def test_success_publishes_spanish_without_modifying_english(self):
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp)
            en = output / 'html/en'
            en.mkdir(parents=True)
            content = ('<html lang="en"><head><meta charset="utf-8"><title>Overview</title></head><body>'
                       '<article><h1 id="heading">The density</h1><p>appears '
                       '<math><mi>ρ</mi></math></p></article></body></html>')
            (en / 'overview.html').write_text(content)
            (en / 'build-manifest.json').write_text(json.dumps({
                'complete': True, 'language': 'en', 'pages': ['overview.html'], 'fds_revision': 'abc'}))
            args = SimpleNamespace(language='es', check=False)
            with patch('translate_argos.ArgosTranslator', return_value=Translator()):
                full.build_translated(args, output)
            self.assertEqual((en / 'overview.html').read_text(), content)
            spanish = html.parse(str(output / 'html/es/overview.html'))
            self.assertEqual(spanish.xpath('/html/@lang'), ['es'])
            self.assertEqual(spanish.xpath('//h1/text()'), ['La densidad'])
            self.assertEqual(spanish.xpath('//math/mi/text()'), ['ρ'])
            registry = json.loads((output / 'html/languages.json').read_text())
            self.assertEqual([entry['code'] for entry in registry['languages']], ['en', 'es'])


if __name__ == '__main__':
    unittest.main()
