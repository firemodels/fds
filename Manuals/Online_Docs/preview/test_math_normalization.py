"""Regression cases from the Tech Guide's converted MathML."""
import unittest
from lxml import etree, html
from build_preview import normalize_math


class MathTypographyTests(unittest.TestCase):
    def convert(self, content):
        body = html.fromstring('<article><math>' + content + '</math></article>')
        normalize_math(body)
        return body

    def test_divergence_and_dotted_mass_2_6(self):
        body = self.convert('<mi class="MathClass-op">∇</mi><mo>⁡</mo>'
                            '<mspace class="negthinspace" width="-0.17em"></mspace>'
                            '<mo>⋅</mo><mi>ρ</mi><mi>u</mi><mo>=</mo><mi>ṁ</mi>')
        self.assertEqual(body.xpath('string(.//mo[text()="⋅"]/@lspace)'), '0.0556em')
        self.assertFalse(body.xpath('.//mi[text()="ṁ"] | .//mspace'))
        self.assertTrue(body.xpath('.//mover/mi[text()="m"]'))
        self.assertTrue(body.xpath('.//mover[@accent="true"]/mo[text()="˙"][@accent="true"]'))
        self.assertFalse(body.xpath('.//mpadded'))
        once = etree.tostring(body)
        normalize_math(body)
        self.assertEqual(once, etree.tostring(body), 'Normalization must be idempotent')

    def test_gradients_2_9(self):
        body = self.convert('<mi>k</mi><mi class="MathClass-op">∇</mi><mo>⁡</mo><mi>T</mi>'
                            '<mo>+</mo><mi>∇</mi><mo>⁡</mo><msub><mi>Z</mi><mi>α</mi></msub>')
        self.assertEqual(len(body.xpath('.//mo[text()="∇"][@rspace="0em"]')), 2)
        self.assertFalse(body.xpath('.//mo[text()="⁡"]'))
        self.assertTrue(body.xpath('.//msub/mi[text()="α"]'))

    def test_laplacian_and_named_function(self):
        body = self.convert('<msup><mrow><mi>∇</mi><mo>⁡</mo></mrow><mn>2</mn></msup>'
                            '<mi>H</mi><mo>+</mo><mi>sin</mi><mo>⁡</mo><mi>x</mi>')
        self.assertTrue(body.xpath('.//msup/mrow/mo[text()="∇"]'))
        self.assertEqual(body.xpath('string(.//msup/mn)'), '2')
        self.assertEqual(len(body.xpath('.//mo[text()="⁡"]')), 1)

    def test_native_dot_does_not_move_other_accents(self):
        body = self.convert('<mover accent="true"><mi>q</mi><mo>˙</mo></mover>'
                            '<mover accent="true"><mi>p</mi><mo>¯</mo></mover>')
        self.assertFalse(body.xpath('.//mpadded'))
        self.assertTrue(body.xpath('.//mover/mo[text()="˙"][@stretchy="false"]'))
        self.assertTrue(body.xpath('.//mover/mo[text()="¯"]'))

    def test_removes_old_dot_lift_with_primes_and_subscripts_3_6(self):
        body = self.convert('<msubsup><mover accent="true"><mrow><mi>m</mi></mrow>'
                            '<mpadded class="dot-accent" voffset="0.03em" height="+0.03em">'
                            '<mo>˙</mo></mpadded></mover><mi>α</mi><mo>‴</mo></msubsup>')
        self.assertFalse(body.xpath('.//mpadded'))
        self.assertTrue(body.xpath('.//msubsup/mover/mi[text()="m"]'))
        self.assertEqual(body.xpath('string(.//msubsup/mi)'), 'α')
        self.assertEqual(body.xpath('string(.//msubsup/mo)'), '‴')

    def test_mixture_weight_overline_in_superscript(self):
        body = self.convert('<msup><mi>T</mi><mrow><mover accent="false" class="mml-overline">'
                            '<mrow><mi>W</mi></mrow><mo>¯</mo></mover><mi>g</mi></mrow></msup>'
                            '<mover accent="true"><mi>p</mi><mo>¯</mo></mover>')
        self.assertTrue(body.xpath('.//msup//mrow[@class="mixture-overline"]/mrow/mi[text()="W"]'))
        self.assertFalse(body.xpath('.//msup//mover | .//msup//mo[text()="¯" or text()="‾"]'))
        self.assertTrue(body.xpath('.//mover[mi="p"]/mo[text()="¯"][not(@stretchy)]'))
        once = etree.tostring(body)
        normalize_math(body)
        self.assertEqual(once, etree.tostring(body))

    def test_sp811_roman_differentials_not_diameter(self):
        body = self.convert('<mfrac><mrow><mstyle mathvariant="normal"><mi>d</mi></mstyle><mi>m</mi></mrow>'
                            '<mrow><mstyle mathvariant="normal"><mi>d</mi></mstyle><mi>t</mi></mrow></mfrac>'
                            '<mi>d</mi><mstyle mathvariant="normal"><mi>D</mi>'
                            '<mstyle mathvariant="italic"><mi>x</mi></mstyle></mstyle>')
        self.assertEqual(len(body.xpath('.//mi[text()="d"][@mathvariant="normal"]')), 2)
        self.assertEqual(len(body.xpath('.//mi[text()="d"][not(@mathvariant)]')), 1)
        self.assertTrue(body.xpath('.//mi[text()="D"][@mathvariant="normal"]'))
        self.assertFalse(body.xpath('.//mi[text()="m" or text()="t" or text()="x"][@mathvariant="normal"]'))

    def test_general_overline_preserves_group_scripts_and_nested_lines(self):
        body = self.convert('<msup><mover class="mml-overline" id="flux">'
                            '<mrow><mi>ρ</mi><msub><mi>Z</mi><mi>α</mi></msub></mrow>'
                            '<mo>¯</mo></mover><mi>FL</mi></msup>'
                            '<mover class="mml-overline"><mrow><mfrac><mi>ρ</mi>'
                            '<mover class="mml-overline"><mrow><mi>W</mi></mrow><mo>¯</mo>'
                            '</mover></mfrac></mrow><mo>¯</mo></mover>'
                            '<mover class="mml-overline"><mi>p</mi><mo>¯</mo></mover>'
                            '<mover accent="true"><mi>p</mi><mo>¯</mo></mover>')
        self.assertFalse(body.xpath('.//mover[@class="mml-overline"]'))
        self.assertTrue(body.xpath('.//msup/mrow[@id="flux"]//msub[mi="Z"][mi="α"]'))
        self.assertEqual(body.xpath('string(.//msup/mi)'), 'FL')
        self.assertTrue(body.xpath('.//mrow[@class="latex-overline"]//mfrac/mrow[@class="mixture-overline"]'))
        self.assertTrue(body.xpath('.//mrow[@class="latex-overline"]/mi[text()="p"]'))
        self.assertTrue(body.xpath('.//mover[@accent="true"]/mo[text()="¯"]'))
        once = etree.tostring(body)
        normalize_math(body)
        self.assertEqual(once, etree.tostring(body))

    def test_sp811_integral_measure(self):
        body = self.convert('<mrow><mo>∫</mo><mi>I</mi><mspace width="0.17em"/>'
                            '<mi>d</mi><msup><mi>s</mi><mo>′</mo></msup></mrow>')
        self.assertTrue(body.xpath('.//mi[text()="d"][@mathvariant="normal"]'))


if __name__ == '__main__':
    unittest.main()
