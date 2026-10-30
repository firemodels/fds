#!/usr/bin/env python3
"""Assemble the Spanish pilot from a cached, source-checked translation.

MathML and structural tags are protected placeholders. The cache is an external
translation artifact, not a second authored manual. Changed source blocks fail
closed until their translations have been refreshed.
"""
import argparse
from collections import Counter
import hashlib
import html as escape_html
import json
from pathlib import Path
import re
from urllib.parse import urlsplit

from lxml import etree, html

HERE = Path(__file__).resolve().parent
SITE = HERE.parent / 'build/tech-preview/html'
CACHE = HERE / 'translations/es-overview.json'
TOKEN = re.compile(r'⟦(\d+)⟧')
BLOCKS = {'p', 'li', 'h1', 'h2', 'h3', 'h4', 'h5', 'figcaption'}


def text(element):
    return ' '.join(''.join(element.itertext()).split())


def translation_units(doc, include_tables=False):
    article = doc.xpath('//article')[0]
    units = []
    blocks = BLOCKS | ({'td', 'th'} if include_tables else set())
    for node in article.iterdescendants():
        if node.tag not in blocks or any(a.tag in blocks | {'math'} | (set() if include_tables else {'table'}) for a in node.iterancestors()):
            continue
        protected = []

        def protect(markup):
            token = f'⟦{len(protected)}⟧'
            protected.append(markup)
            return token

        def contents(parent):
            result = parent.text or ''
            for child in parent:
                if not isinstance(child.tag, str):
                    result += child.tail or ''
                    continue
                classes = child.get('class', '')
                if child.tag in {'math', 'a', 'img', 'table', 'div', 'code', 'br'} or 'pcr' in classes or 'titlemark' in classes:
                    result += protect(etree.tostring(child, encoding='unicode', method='html', with_tail=False))
                else:
                    opening = '<' + child.tag + ''.join(' ' + k + '="' + escape_html.escape(v, quote=True) + '"'
                                                       for k, v in child.attrib.items()) + '>'
                    result += protect(opening) + contents(child) + protect('</' + child.tag + '>')
                result += child.tail or ''
            return result

        source = ' '.join(contents(node).split())
        if not re.search(r'[A-Za-z]', TOKEN.sub('', source)):
            continue
        # The source digest deliberately excludes serialized math/formatting:
        # immutable placeholders are always restored from the current build.
        units.append({'id': f'b{len(units):03}', 'node': node, 'source': source,
                      'sha256': hashlib.sha256(source.encode()).hexdigest(), 'protected': protected})
    return units


def apply_translation(unit, translated):
    if Counter(TOKEN.findall(unit['source'])) != Counter(TOKEN.findall(translated)):
        raise ValueError(f"Protected placeholder mismatch: {unit['id']}")
    # This pilot preserves protected-item order as well as exact markup.
    if TOKEN.findall(unit['source']) != TOKEN.findall(translated):
        raise ValueError(f"Protected placeholder order changed: {unit['id']}")
    markup = TOKEN.sub(lambda m: unit['protected'][int(m[1])], escape_html.escape(translated))
    replacement = html.fragment_fromstring(markup, create_parent=unit['node'].tag)
    node = unit['node']
    for child in list(node):
        node.remove(child)
    node.text = replacement.text
    node.extend(list(replacement))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--extract', action='store_true')
    parser.add_argument('--site', type=Path, default=SITE)
    args = parser.parse_args()
    site = args.site.resolve()
    doc = html.parse(str(site / 'overview.html'))
    units = translation_units(doc)
    if args.extract:
        print(json.dumps([{k: u[k] for k in ['id', 'sha256', 'source']} for u in units], ensure_ascii=False, indent=2))
        return
    cache = json.loads(CACHE.read_text())
    entries = {entry['id']: entry for entry in cache['blocks']}
    if set(entries) != {u['id'] for u in units}:
        raise ValueError('Translation inventory changed; refresh the cached translation.')
    original_math = [etree.tostring(n, with_tail=False) for n in doc.xpath('//math')]
    original_ids = doc.xpath('//@id')
    for unit in units:
        entry = entries[unit['id']]
        if entry['source_sha256'] != unit['sha256']:
            raise ValueError(f"Translation needs refresh: {unit['id']}")
        apply_translation(unit, entry['translation'])
    localize_reader(doc)
    # Captions provide the translated descriptions; image bytes stay shared.
    for figure in doc.xpath('//article//figure'):
        captions = figure.xpath('.//figcaption')
        if captions:
            for image in figure.xpath('.//img'):
                image.set('alt', text(captions[0]))
    assert original_math == [etree.tostring(n, with_tail=False) for n in doc.xpath('//math')], 'Math changed'
    assert original_ids == doc.xpath('//@id'), 'Anchors changed'
    target = site / 'es'
    target.mkdir(exist_ok=True)
    (target / 'overview.html').write_text('<!doctype html>\n' + etree.tostring(doc, encoding='unicode', method='html'))
    make_search(doc, target)
    for page in site.glob('*.html'):
        if html.parse(str(page)).xpath('//div[@class="header-links"]'):
            add_language_selector(page, spanish=False)
    add_language_selector(target / 'overview.html', spanish=True)
    validate_translation(target / 'overview.html')
    (target / 'translation-manifest.json').write_text(json.dumps({
        'source_page': '../overview.html', 'language': 'es',
        'cache_sha256': hashlib.sha256(CACHE.read_bytes()).hexdigest(),
        'provenance': cache['provenance'], 'translated_blocks': len(units),
        'unchanged_math_expressions': len(original_math),
    }, ensure_ascii=False, indent=2) + '\n')
    print(f'Spanish pilot: {len(units)} translated blocks; {len(original_math)} math expressions unchanged.')


def localize_reader(doc):
    """Localize interface text separately from the cached chapter content."""
    replacements = {
        '//a[@class="skip-link"]': 'Saltar al contenido',
        '//span[@class="brand-name"]/span': 'DOCUMENTACIÓN TÉCNICA',
        '//span[@class="edition"]': 'SEXTA EDICIÓN',
        '//div[@class="header-links"]/a': 'Manual PDF (inglés) ↗',
        '//div[@class="manual-name"]': 'Guía de referencia técnica',
        '//div[@class="volume"]': 'Volumen 1 · Modelo matemático',
        '//button[@class="search-trigger"]/span': '⌕  Buscar en este capítulo',
        '//div[@class="nav-heading"]': 'EN ESTA VISTA PREVIA',
        '//div[@class="sidebar-footer"]/p': 'Traducción automática de prueba. Solo este capítulo está disponible en español; las otras páginas y los textos de las figuras siguen en inglés.',
        '//div[@class="sidebar-footer"]/a': 'Descargar el PDF en inglés ↓',
        '//div[@class="breadcrumb"]/a': 'Guía de referencia técnica',
        '//div[@class="breadcrumb"]/span[last()]': 'Modelo matemático',
        '//div[@class="eyebrow"]': 'CAPÍTULO 2 · ESPAÑOL',
        '//aside[@class="page-toc"]/div': 'EN ESTA PÁGINA',
        '//a[@class="back-top"]': 'Volver arriba ↑',
        '//label[@id="search-title"]': 'Buscar en el manual',
        '//div[@id="search-status"]': 'Buscar en el capítulo traducido.',
        '//nav[@class="page-navigation"]/a/small': 'Siguiente · en inglés',
        '//footer[@class="page-footer"]/span[1]': 'FDS · Guía de referencia técnica',
        '//footer[@class="page-footer"]/a': 'Ver el código LaTeX original ↗',
    }
    for xpath, value in replacements.items():
        for node in doc.xpath(xpath):
            for child in list(node):
                node.remove(child)
            node.text = value
    doc.getroot().set('lang', 'es')
    doc.find('head/title').text = 'Descripción general del modelo FDS — Guía de referencia técnica'
    for node in doc.xpath('//div[@class="sidebar-footer"]/span'):
        node.tail = ' Traducción de prueba'
    footer = doc.xpath('//footer[@class="page-footer"]/span[2]')[0]
    footer.text = footer.text.replace('Source revision', 'Revisión del original')
    headings = {node.get('id'): text(node) for node in doc.xpath('//article//h2[@id] | //article//h3[@id]')}
    for link in doc.xpath('//aside[@class="page-toc"]//nav/a'):
        link.text = headings[link.get('href')[1:]]
    names = ['Descripción general del modelo FDS', 'Transporte de masa, especies y entalpía (inglés)',
             'Partículas lagrangianas (inglés)', 'Referencias (inglés)']
    for link, name in zip(doc.xpath('//nav[@class="chapter-nav"]/a'), names):
        link.findall('span')[-1].text = name
    pager = doc.xpath('//nav[@class="page-navigation"]/a')[0]
    pager.find('small').tail = 'Transporte de masa, especies y entalpía'
    labels = {'Open chapter navigation': 'Abrir la navegación por capítulos',
              'Chapter navigation': 'Navegación por capítulos', 'On this page': 'En esta página',
              'Previous and next pages': 'Páginas anterior y siguiente', 'Close search': 'Cerrar búsqueda'}
    for node in doc.xpath('//*[@aria-label]'):
        node.set('aria-label', labels.get(node.get('aria-label'), node.get('aria-label')))
    doc.xpath('//input[@id="search-input"]')[0].set('placeholder', 'Pruebe «malla escalonada» o «combustión»…')
    for node in doc.xpath('//a[@title]'):
        node.set('title', 'Abrir esta referencia en el manual PDF en inglés')
    for node in doc.xpath('//*[@href] | //*[@src]'):
        attr = 'href' if node.get('href') is not None else 'src'
        value = node.get(attr)
        url = urlsplit(value)
        if not url.scheme and not url.netloc and url.path and url.path != 'overview.html':
            node.set(attr, '../' + value)


def make_search(doc, target):
    entries, title, anchor, texts = [], '', '', []
    for node in doc.xpath('//article')[0]:
        if node.tag in {'h1', 'h2', 'h3'}:
            if texts:
                entries.append({'title': title, 'chapter': 'Descripción general del modelo FDS',
                                'url': 'overview.html#' + anchor, 'text': ' '.join(texts)})
            title, anchor, texts = text(node), node.get('id', ''), []
        elif isinstance(node.tag, str):
            texts.append(text(node))
    if texts:
        entries.append({'title': title, 'chapter': 'Descripción general del modelo FDS',
                        'url': 'overview.html#' + anchor, 'text': ' '.join(texts)})
    (target / 'search-index.json').write_text(json.dumps(entries, ensure_ascii=False))


def add_language_selector(path, spanish):
    doc = html.parse(str(path))
    for old in doc.xpath('//*[@class="language-picker"]'):
        old.getparent().remove(old)
    picker = etree.Element('label', {'class': 'language-picker'})
    label = etree.SubElement(picker, 'span')
    label.text = 'Idioma' if spanish else 'Language'
    select = etree.SubElement(picker, 'select', {'aria-label': 'Idioma' if spanish else 'Language'})
    overview = path.name == 'overview.html'
    for code, name, url in [('en', 'English', '../overview.html' if spanish else path.name),
                            ('es', 'Español' if overview else 'Español (Overview)',
                             'overview.html' if spanish else 'es/overview.html')]:
        option = etree.SubElement(select, 'option', value=url)
        option.set('data-keep-section', 'true' if overview or code == 'en' else 'false')
        option.text = name
        if (code == 'es') == spanish:
            option.set('selected', 'selected')
    doc.xpath('//div[@class="header-links"]')[0].insert(1, picker)
    path.write_text('<!doctype html>\n' + etree.tostring(doc, encoding='unicode', method='html'))


def validate_translation(path):
    doc = html.parse(str(path))
    for node in doc.xpath('//*[@href] | //*[@src] | //option[@value]'):
        value = node.get('href') or node.get('src') or node.get('value')
        url = urlsplit(value)
        if url.scheme or url.netloc:
            continue
        target = (path.parent / url.path).resolve() if url.path else path
        if not target.is_file():
            raise ValueError(f'Missing translated page asset/link: {value}')
        if url.fragment and target.suffix == '.html':
            assert url.fragment in html.parse(str(target)).xpath('//@id'), value


if __name__ == '__main__':
    main()
