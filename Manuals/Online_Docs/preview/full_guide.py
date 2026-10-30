"""Full FDS manual orchestration; --check never runs TeX or writes build files."""
from __future__ import annotations
import json
import os
from pathlib import Path
import re
import shutil
import tempfile
from urllib.parse import urlsplit

from lxml import etree, html
import build_preview as reader

from manuals import MANUALS, TECH
from manual_sources import GRAPHIC, uncomment, active_geometry, slug
from manual_sources import inventory as source_inventory

MAIN = TECH.main
LANGUAGES = {'en': 'English', 'es': 'Español'}
EXTERNAL_DOCUMENT = re.compile(r'\\externaldocument(?:\[([^\]]*)\])?\{([^{}]+)\}')
REFERENCE = re.compile(r'\\(ref|eqref|pageref)\{([^{}]+)\}')


def inventory(manual=TECH):
    inv = source_inventory(manual, reader.MANUALS)
    inv.update(external_references(inv))
    return inv


def external_references(inv):
    """Read the source's declared cross-manual labels from its PDF build."""
    references, files = {}, {}
    text = '\n'.join(inv['sources'].values())
    keys = {m[2] for m in REFERENCE.finditer(text)}
    local_labels = set(re.findall(r'\\label\{([^{}]+)\}', text))
    for match in EXTERNAL_DOCUMENT.finditer(text):
        prefix, name = match[1] or '', match[2]
        other = next((m for m in MANUALS.values() if m.directory == name), None)
        if other is None:
            inv['missing'].append('Unsupported external manual: ' + name)
            continue
        guide = reader.MANUALS / other.directory
        aux, pdf = guide / (name + '.aux'), guide / other.pdf
        if not aux.is_file() or not pdf.is_file():
            inv['missing'].append(f'External reference inputs: build {other.pdf} and its .aux files first')
            continue
        labels = reader.read_labels(guide)
        for path in guide.glob('*.aux'):
            files[str(path)] = path
        files[str(pdf)] = pdf
        for key in keys:
            if not key.startswith(prefix) or key in local_labels:
                continue
            label = labels.get(key[len(prefix):])
            if label is None:
                # With no prefix, only labels actually imported are external.
                if prefix:
                    inv['missing'].append('External reference label: ' + key)
                continue
            references[key] = {**label, 'pdf': pdf}
    return {'external_references': references, 'external_files': files}


def resolve_external_references(text, references):
    text = EXTERNAL_DOCUMENT.sub('', text)
    def replace(match):
        command, key = match.groups()
        if key not in references:
            return match[0]
        target = references[key]
        number = target['page'] if command == 'pageref' else target['number']
        # Printed page numbers differ from physical PDF indices after front matter.
        destination = 'page.' + target['page'] if command == 'pageref' else target['destination']
        fragment = 'nameddest=' + destination
        link = (r'\href{external-manuals/' + target['pdf'].name + r'\#' + fragment + '}{' + number + '}')
        return '(' + link + ')' if command == 'eqref' else link
    return REFERENCE.sub(replace, text)


def expected_pages(inv):
    return (['fds-developers', 'preface', 'disclaimer', 'acknowledgments', 'references']
            + [f'part-{i + 1}' for i in range(len(inv.get('part_sources', [])))]
            + [slug(name) for name in inv['chapters']]
            + inv.get('appendix_slugs', ['appendix-' + chr(97 + i) for i in range(inv['appendices'])]))


def cache_path(language, page):
    path = reader.HERE / 'translations' / language / (page + '.json')
    if not path.exists() and language == 'es' and page == 'overview':
        return reader.HERE / 'translations/es-overview.json'
    return path


def preflight(inv, language):
    manual = inv.get('manual', TECH)
    errors = list(inv['missing'])
    env_path = '/Library/TeX/texbin:/opt/homebrew/bin:' + os.environ.get('PATH', '')
    for tool in ['make4ht', 'biber', 'dvisvgm', 'kpsewhich', 'pdftocairo', 'pdftops']:
        if not shutil.which(tool, path=env_path):
            errors.append('Missing tool: ' + tool)
    if not (reader.MANUALS / manual.directory / manual.pdf).is_file():
        errors.append(f'Missing PDF manual: {reader.MANUALS / manual.directory / manual.pdf}')
    print(f"{manual.title}: {len(inv['chapters'])} chapters, {inv['appendices']} appendices, "
          f"{len(inv['figures'])} external figures; complex geometry {'enabled' if inv['complex_geometry'] else 'disabled'}.", flush=True)
    if errors:
        raise RuntimeError('Input checks failed:\n' + '\n'.join(errors))


def stage(inv, work):
    manual = inv.get('manual', TECH)
    tex = work / 'Manuals' / manual.directory
    tex.mkdir(parents=True)
    (work / '.fds-docs-work').write_text('Generated full-guide staging files.\n')
    for name in ['Bibliography', 'LaTeX_Style_Files']:
        dst = tex.parent / name
        dst.mkdir()
        for path in (reader.MANUALS / name).iterdir():
            if path.is_file() and path.suffix in {'.tex', '.bib', '.sty', '.cls', '.bst'}:
                shutil.copy2(path, dst / path.name)
    common = tex.parent / 'Bibliography/commoncommands.tex'
    common.write_text(common.read_text().replace(r'\usepackage[pdftex]{graphicx}', r'\usepackage{graphicx}')
                      .replace(r'\usepackage[pdftex,', r'\usepackage[')
                      .replace(r'\usepackage{colortbl}', '% Table colors supplied by reader CSS.')
                      .replace(r'\usepackage{xltabular}', '% Staged as longtable for TeX4ht.')
                      .replace('backref=true', 'backref=false'))
    folder = tex / 'WEB_FIGURES'
    folder.mkdir()
    cache = work.parent / 'figure-cache-v1'
    cache.mkdir(exist_ok=True)
    replacements, web_figures = {}, {}
    for i, (reference, source) in enumerate(sorted(inv['figures'].items()), 1):
        stem = re.sub(r'[^A-Za-z0-9_-]', '-', source.stem) + '-' + reader.sha(source)[:12]
        extensions = (['.eps'] if reference in inv['plot_graphics'] else
                      ['.png', '.svg'] if source.suffix.lower() == '.pdf' else [source.suffix.lower()])
        cached = all((cache / (stem + ext)).is_file() for ext in extensions)
        print(f"Preparing figure {i}/{len(inv['figures'])}: {reference}" + (' (cached)' if cached else ''), flush=True)
        if cached:
            for ext in extensions:
                shutil.copy2(cache / (stem + ext), folder / (stem + ext))
        if reference in inv['plot_graphics']:
            # PGFPlots backgrounds need real DVI graphics, not TeX4ht's BMP
            # placeholder specials (which silently disappear from the SVG).
            eps = folder / (stem + '.eps')
            if cached:
                pass
            elif source.suffix.lower() == '.pdf':
                reader.run(['pdftops', '-eps', str(source), str(eps)], tex)
            else:
                raster_to_eps(source, eps)
            if not cached:
                cache_figure(eps, cache / eps.name)
            replacements[reference] = str(eps.relative_to(tex))
            continue
        if source.suffix.lower() == '.pdf':
            raster = folder / (stem + '.png')
            vector = folder / (stem + '.svg')
            if not cached:
                reader.run(['pdftocairo', '-png', '-singlefile', '-r', '160', str(source), str(folder / stem)], tex)
                reader.run(['pdftocairo', '-svg', str(source), str(vector)], tex)
            web_figures[str(raster.relative_to(tex))] = str(vector.relative_to(tex))
        else:
            raster = folder / (stem + source.suffix.lower())
            shutil.copy2(source, raster)
        replacements[reference] = str(raster.relative_to(tex))
        if not cached:
            for ext in extensions:
                cache_figure(folder / (stem + ext), cache / (stem + ext))
    if inv.get('external_references'):
        external = tex / 'external-manuals'
        external.mkdir()
        for pdf in {target['pdf'] for target in inv['external_references'].values()}:
            shutil.copy2(pdf, external / pdf.name)
    for name, text in inv['sources'].items():
        text = resolve_external_references(text, inv.get('external_references', {}))
        text = GRAPHIC.sub(lambda m: m[0].replace('{' + m[1] + '}', '{' + replacements.get(m[1], m[1]) + '}'), text)
        text = prepare_tables(reader.prepare_math(text))
        target = tex / ('preview.tex' if name == manual.main else name)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(text)
    (tex / 'web-figures.json').write_text(json.dumps(web_figures, indent=2))
    (tex / 'preview.mk4').write_text('Make:htlatex {}\nMake:biber {}\nMake:htlatex {}\nMake:htlatex {}\n')
    return tex


def cache_figure(source, destination):
    # An interrupted copy must never become a usable cache entry.
    with tempfile.NamedTemporaryFile(dir=destination.parent, delete=False) as stream:
        temporary = Path(stream.name)
    try:
        shutil.copy2(source, temporary)
        temporary.replace(destination)
    finally:
        temporary.unlink(missing_ok=True)


def raster_to_eps(source, destination):
    """Give DVI/dvisvgm actual raster artwork inside a PGFPlots background."""
    from PIL import Image
    with Image.open(source) as original:
        rgba = original.convert('RGBA')
        background = Image.new('RGBA', rgba.size, 'white')
        background.alpha_composite(rgba)
        background.convert('RGB').save(destination, format='EPS')


def prepare_tables(text):
    """Use native longtable handling instead of TeX4ht's tabularx alias.

    The alias loses multipage captions and breaks listing cells. HTML controls
    table width; keep the cells, captions, labels, and repeated header structure.
    """
    # Print-page rotation produces invalid inline wrappers around block tables.
    # The reader presents these as upright, horizontally scrollable tables.
    text = re.sub(r'\\(?:begin|end)\{landscape\}', '', text)
    text = text.replace('{sidewaystable}', '{table}')
    # A legacy table supplies a width where tabular accepts only t/b/c alignment.
    # PDF LaTeX tolerates it, but TeX4ht treats it as an invalid control sequence.
    text = re.sub(r'(\\begin\{tabular\})\s*\[\\textwidth\]', r'\1', text)
    # Listings inside math text boxes emit HTML code tags across MathML nodes.
    # Keep literal monospaced labels as native text inside the existing box.
    def boxed_code(match):
        escapes = {'_': r'\_', '#': r'\#', '%': r'\%', '&': r'\&', '$': r'\$'}
        literal = ''.join(escapes.get(char, char) for char in match[1])
        return r'\mbox{\texttt{' + literal + '}}'
    text = re.sub(r'\\mbox\{\\ct\{([^{}\\]*)\}\}', boxed_code, text)
    pattern = re.compile(r'\\begin\{xltabular\}(?:\[[^\]]*\])?\s*(?=\{)')
    parts, pos = [], 0
    while match := pattern.search(text, pos):
        _, end = reader.group(text, match.end())
        while text[end].isspace():
            end += 1
        columns, end = reader.group(text, end)
        output, depth = [], 0
        for char in columns:
            output.append(r'p{0.2\linewidth}' if char == 'X' and depth == 0 else char)
            if char == '{':
                depth += 1
            elif char == '}':
                depth -= 1
        parts.extend([text[pos:match.start()], r'\begin{longtable}{', ''.join(output), '}'])
        pos = end
    parts.append(text[pos:])
    text = ''.join(parts).replace(r'\end{xltabular}', r'\end{longtable}')
    text = re.sub(r'\\(?:rowcolor|cellcolor)(?:\[[^\]]*\])?\{[^{}]*\}(?:\[[^\]]*\]){0,2}', '', text)
    def inline_math(match):
        delimiter, content = match.groups()
        pieces = re.split(r'(?<!\\)\$', content)
        if len(pieces) % 2 == 0:
            raise RuntimeError('Unbalanced math escape in inline listing: ' + content)
        # TeX4ht must close typewriter mode before entering native math.
        return ''.join(('$' + part + '$') if i % 2 else
                       (r'\lstinline' + delimiter + part + delimiter if part else '')
                       for i, part in enumerate(pieces))
    return re.sub(r'\\lstinline\[mathescape=true\]([^\w\s{])(.*?)\1', inline_math, text)


def discover_pages(tex, inv):
    manual = inv.get('manual', TECH)
    pages = []
    for path in tex.glob('preview*.html'):
        doc = html.parse(str(path), parser=html.HTMLParser(encoding='utf-8'))
        headings = doc.xpath('//body/h1 | //body/h2')
        if not headings:
            if path.name == 'preview.html' and doc.xpath('//div[@class="tableofcontents"]'):
                pages.append({'raw': path.name, 'file': 'manual-index.html', 'source': manual.main,
                              'number': None, 'eyebrow': 'CONTENTS', 'order': 0,
                              'synthetic_title': 'Manual contents'})
            continue
        if len(headings) != 1:
            raise RuntimeError(f'Expected one chapter per HTML page: {path}')
        heading = headings[0]
        marker = ' '.join(heading.xpath('./span[@class="titlemark"]//text()'))
        title = reader.plain(heading)
        number = re.search(r'Chapter\s+(\d+)', marker)
        appendix = re.search(r'Appendix\s+([A-Z]+)', marker)
        part = re.search(r'Part\s+([IVXLCDM]+)', marker)
        if part:
            values = {'I': 1, 'V': 5, 'X': 10, 'L': 50, 'C': 100, 'D': 500, 'M': 1000}
            digits = [values[c] for c in part[1]]
            index = sum(-v if i + 1 < len(digits) and v < digits[i + 1] else v for i, v in enumerate(digits))
            source = inv['part_sources'][index - 1]
            name, eyebrow, number = f'part-{index}', 'PART ' + part[1], None
        elif number:
            number = int(number[1])
            source = inv.get('chapter_sources', [name + '.tex' for name in inv['chapters']])[number - 1]
            name, eyebrow = slug(inv['chapters'][number - 1]), 'CHAPTER ' + str(number)
        elif appendix:
            number = appendix[1]
            name, eyebrow = 'appendix-' + number.lower(), 'APPENDIX ' + number
            source = dict(zip(inv.get('appendix_slugs', []), inv.get('appendix_sources', []))).get(name, 'Appendices.tex')
        else:
            number, source, name, eyebrow = None, manual.main, slug(title), 'FRONT MATTER'
            if name == 'references':
                source, eyebrow = None, 'REFERENCES'
        ids = heading.xpath('.//@id') + heading.xpath('./@id')
        order = [int(m[1]) for value in ids if (m := re.match(r'x(\d+)-', value))]
        if not order:
            raise RuntimeError(f'Cannot determine page order for {path}')
        pages.append({'raw': path.name, 'file': name + '.html', 'source': source,
                      'number': number, 'eyebrow': eyebrow, 'order': min(order)})
    pages.sort(key=lambda page: page['order'])
    filenames = [p['file'] for p in pages]
    if len(set(filenames)) != len(filenames):
        raise RuntimeError('Duplicate generated page names')
    missing = set(expected_pages(inv)) - {Path(p).stem for p in filenames}
    if missing:
        raise RuntimeError('Conversion omitted expected pages: ' + ', '.join(sorted(missing)))
    return pages


def translate_edition(site, pages, language, translator=None, cache_root=None, manual=TECH):
    """Translate prose offline, retaining all scientific markup and anchors."""
    from translate_preview import translation_units, apply_translation
    documents, titles, stats, heading_labels = {}, {}, {}, {}
    for position, page in enumerate(pages, 1):
        file = site / page['file']
        doc = html.parse(str(file), parser=html.HTMLParser(encoding='utf-8'))
        units = translation_units(doc, include_tables=True)
        if translator is not None:
            from translate_argos import generate_cache
            print(f'Translating page {position}/{len(pages)}: {file.stem}', flush=True)
            cache = cache_root / language / (file.stem + '.json')
            blocks, stats[file.stem] = generate_cache(units, cache, translator)
            entries = {entry['id']: entry for entry in blocks}
        else:
            cache = cache_path(language, file.stem)
            if not cache.is_file():
                raise RuntimeError(f'Missing translation cache: {cache}')
            entries = {entry['id']: entry for entry in json.loads(cache.read_text())['blocks']}
        if set(entries) != {unit['id'] for unit in units}:
            raise RuntimeError(f'Translation inventory changed: {cache}')
        math_before = [etree.tostring(n, with_tail=False) for n in doc.xpath('//math')]
        ids_before = doc.xpath('//@id')
        for unit in units:
            entry = entries[unit['id']]
            if entry['source_sha256'] != unit['sha256']:
                raise RuntimeError(f'Stale translation: {cache}, {unit["id"]}')
            apply_translation(unit, entry['translation'])
        assert ids_before == doc.xpath('//@id'), 'Translation changed anchors'
        assert math_before == [etree.tostring(n, with_tail=False) for n in doc.xpath('//math')], 'Translation changed math'
        titles[page['file']] = reader.plain(doc.xpath('//article/h1')[0])
        for heading in doc.xpath('//article//*[self::h1 or self::h2 or self::h3 or self::h4 or self::h5]'):
            label = reader.plain(heading)
            marker = heading.find('span')
            if marker is not None and marker.get('class') == 'titlemark':
                label = label.removeprefix(reader.plain(marker)).strip()
            for anchor in heading.xpath('.//@id | ./@id'):
                heading_labels[(page['file'], anchor)] = label
        documents[page['file']] = doc
    search = []
    for name, doc in documents.items():
        doc.getroot().set('lang', language)
        doc.find('head/title').text = titles[name] + ' — FDS'
        # Translate the reader separately from protected chapter markup.
        localize_full_interface(doc, manual)
        for link in doc.xpath('//article//div[@class="tableofcontents"]//a[@href]'):
            target = urlsplit(link.get('href'))
            label = heading_labels.get((target.path or name, target.fragment))
            if label is not None and not len(link):
                link.text = label
        for link in doc.xpath('//nav[@class="chapter-nav"]/a'):
            link.findall('span')[-1].text = titles.get(link.get('href'), link.findall('span')[-1].text)
        for link in doc.xpath('//nav[@class="page-navigation"]/a'):
            link.find('small').tail = titles.get(link.get('href'), link.find('small').tail)
        headings = {n.get('id'): reader.plain(n) for n in doc.xpath('//article//h2[@id] | //article//h3[@id]')}
        for link in doc.xpath('//aside[@class="page-toc"]//nav/a'):
            link.text = headings[link.get('href')[1:]]
        from accessibility import figure_labels
        figure_labels(doc.getroot(), language)
        for picker in doc.xpath('//label[@class="language-picker"]'):
            picker.find('span').text = 'Idioma'
            select = picker.find('select')
            select.set('aria-label', 'Idioma')
            for child in list(select):
                select.remove(child)
            etree.SubElement(select, 'option', value=name, selected='selected').text = LANGUAGES[language]
        (site / name).write_text('<!doctype html>\n' + etree.tostring(doc, encoding='unicode', method='html'))
        title, anchor, texts = titles[name], '', []
        for node in doc.xpath('//article')[0]:
            if node.tag in {'h1', 'h2', 'h3'}:
                if texts:
                    search.append({'title': title, 'chapter': titles[name], 'url': name + ('#' + anchor if anchor else ''), 'text': ' '.join(texts)})
                title, anchor, texts = reader.plain(node), node.get('id', ''), []
            elif isinstance(node.tag, str):
                texts.append(reader.plain(node))
        if texts:
            search.append({'title': title, 'chapter': titles[name], 'url': name + ('#' + anchor if anchor else ''), 'text': ' '.join(texts)})
    (site / 'search-index.json').write_text(json.dumps(search, ensure_ascii=False))
    return stats


def localize_full_interface(doc, manual=TECH):
    labels = {
        '//a[@class="skip-link"]': 'Saltar al contenido',
        '//span[@class="brand-name"]/span': 'DOCUMENTACIÓN TÉCNICA',
        '//span[@class="edition"]': 'SEXTA EDICIÓN',
        '//div[@class="header-links"]/a': 'Manual PDF (inglés) ↗',
        '//div[@class="manual-name"]': manual.title_es,
        '//div[@class="volume"]': manual.subtitle_es,
        '//button[@class="search-trigger"]/span': '⌕ Buscar en este manual',
        '//div[@class="nav-heading"]': 'EN ESTE MANUAL',
        '//div[@class="sidebar-footer"]/p': 'Traducción automática. Las figuras y el PDF conservan el contenido original en inglés.',
        '//div[@class="sidebar-footer"]/a': 'Descargar el PDF en inglés ↓',
        '//div[@class="breadcrumb"]/a': manual.title_es,
        '//div[@class="breadcrumb"]/span[last()]': manual.subtitle_es,
        '//aside[@class="page-toc"]/div': 'EN ESTA PÁGINA',
        '//a[@class="back-top"]': 'Volver arriba ↑',
        '//label[@id="search-title"]': 'Buscar en el manual',
        '//div[@id="search-status"]': 'Buscar en este manual.',
        '//footer[@class="page-footer"]/span[1]': 'FDS · ' + manual.title_es,
        '//footer[@class="page-footer"]/a': 'Ver el código LaTeX original ↗',
    }
    for xpath, label in labels.items():
        for node in doc.xpath(xpath):
            for child in list(node):
                node.remove(child)
            node.text = label
    for node in doc.xpath('//div[@class="eyebrow"]'):
        node.text = node.text.replace('CHAPTER', 'CAPÍTULO').replace('APPENDIX', 'APÉNDICE').replace('FRONT MATTER', 'PRELIMINARES').replace('REFERENCES', 'REFERENCIAS')
    for node in doc.xpath('//nav[@class="page-navigation"]/a/small'):
        node.text = 'Anterior' if node.text == 'Previous' else 'Siguiente'
    for node in doc.xpath('//input[@id="search-input"]'):
        node.set('placeholder', 'Buscar términos o temas…')
    for node in doc.xpath('//div[@class="sidebar-footer"]/span'):
        node.tail = ' Edición traducida'
    for node in doc.xpath('//footer[@class="page-footer"]/span[2]'):
        node.text = node.text.replace('Source revision', 'Revisión del original')
    aria = {'Open chapter navigation': 'Abrir la navegación por capítulos',
            'Chapter navigation': 'Navegación por capítulos', 'On this page': 'En esta página',
            'Previous and next pages': 'Páginas anterior y siguiente', 'Close search': 'Cerrar búsqueda',
            'Scrollable equation or table': 'Ecuación o tabla desplazable'}
    for node in doc.xpath('//*[@aria-label]'):
        node.set('aria-label', aria.get(node.get('aria-label'), node.get('aria-label')))


def publish(candidate, root, language):
    """Call only after validation. Registry advertises completed editions only."""
    root.mkdir(parents=True, exist_ok=True)
    destination, backup = root / language, root / ('.previous-' + language)
    if backup.exists():
        raise RuntimeError(f'Recover or remove the prior interrupted publication first: {backup}')
    old = destination.exists()
    if old:
        destination.rename(backup)
    try:
        candidate.rename(destination)
    except Exception:
        if old:
            backup.rename(destination)
        raise
    if old:
        shutil.rmtree(backup)
    refresh_registry(root)


def refresh_registry(root):
    # Reconstruct rather than merge: missing/failed editions cannot survive in
    # the registry. Serialize final publication in a multi-job Actions workflow.
    root.mkdir(parents=True, exist_ok=True)
    editions = []
    for code, label in LANGUAGES.items():
        proof = root / code / 'build-manifest.json'
        if not proof.is_file():
            continue
        data = json.loads(proof.read_text())
        pages = data.get('pages', [])
        if (data.get('complete') and data.get('language') == code and pages
                and all(Path(p).name == p and (root / code / p).is_file() for p in pages)):
            editions.append({'code': code, 'label': label, 'pages': pages, 'revision': data['fds_revision']})
    tmp = root / '.languages.json.tmp'
    tmp.write_text(json.dumps({'languages': editions}, ensure_ascii=False, indent=2) + '\n')
    tmp.replace(root / 'languages.json')
    (root / 'index.html').write_text('<!doctype html><meta charset="utf-8"><title>FDS Online Manual</title>'
        + ''.join(f'<p><a href="{entry["code"]}/index.html">{entry["label"]}</a></p>' for entry in editions))
    if root.name == 'html' and root.parent.name in {m.output for m in MANUALS.values()}:
        from build_catalog import refresh_catalog
        refresh_catalog(root.parent.parent)


def input_signature(inv):
    signature = {name: reader.hashlib.sha256(text.encode()).hexdigest() for name, text in inv['sources'].items()}
    # Conversion rules affect staged TeX even when manual sources are unchanged.
    signature['conversion-rules'] = 3
    signature.update({str(path): reader.sha(path) for path in inv.get('external_files', {}).values()})
    signature.update({str(path): reader.sha(path) for path in inv['figures'].values()})
    for directory in ['Bibliography', 'LaTeX_Style_Files']:
        signature.update({str(path): reader.sha(path) for path in (reader.MANUALS / directory).iterdir()
                          if path.is_file() and path.suffix in {'.tex', '.bib', '.sty', '.cls', '.bst'}})
    return signature


def english_edition(output, manual=TECH):
    source = output / 'html/en'
    manifest = source / 'build-manifest.json'
    message = f'Build English first: ./preview/build_manual.sh --manual {manual.key} --language en'
    if not manifest.is_file():
        raise RuntimeError('No completed English edition found. ' + message)
    data = json.loads(manifest.read_text())
    if data.get('manual', 'tech') != manual.key:
        raise RuntimeError(f'English edition belongs to another manual. {message}')
    names = data.get('pages', [])
    if (not data.get('complete') or data.get('language') != 'en' or not names
            or any(Path(name).name != name or not (source / name).is_file() for name in names)):
        raise RuntimeError('The English edition is incomplete. ' + message)
    return source, data, [{'file': name} for name in names]


def build_translated(args, output):
    from translate_argos import ArgosTranslator
    manual = MANUALS[getattr(args, 'manual', 'tech')]
    source, english, pages = english_edition(output, manual)
    reader.validate(source, pages)
    translator = ArgosTranslator(args.language)
    print(f'Translating completed English edition ({english["fds_revision"]}) '
          f'to {LANGUAGES[args.language]}; no LaTeX conversion.', flush=True)
    if args.check:
        print('English edition and installed Argos model checks passed. No translation performed.')
        return
    source_hashes = {p['file']: reader.sha(source / p['file']) for p in pages}
    candidate = Path(tempfile.mkdtemp(prefix='edition-', dir=output))
    try:
        shutil.copytree(source, candidate, dirs_exist_ok=True)
        (candidate / 'build-manifest.json').unlink()
        stats = translate_edition(candidate, pages, args.language, translator, output / 'translation-cache', manual)
        counts = reader.validate(candidate, pages)
        # Fail safely if another job replaced English while translation ran.
        if source_hashes != {p['file']: reader.sha(source / p['file']) for p in pages}:
            raise RuntimeError('English changed during translation; rerun this language build.')
        (candidate / 'index.html').write_text('<!doctype html><html lang="' + args.language + '">'
            f'<meta charset="utf-8"><meta http-equiv="refresh" content="0;url={manual.entry}">'
            f'<title>FDS</title><a href="{manual.entry}">' + LANGUAGES[args.language] + '</a></html>')
        manifest = {**english, 'complete': True, 'language': args.language, 'validation': counts,
                    'translation': translator.provenance, 'translation_stats': stats,
                    'english_pages_sha256': source_hashes,
                    'limitations': english.get('limitations', []) + [
                        'Machine translation; technical terminology requires human review.',
                        'Figure artwork, PDF, and protected scientific notation remain as in English.']}
        (candidate / 'build-manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2))
        publish(candidate, output / 'html', args.language)
    finally:
        if candidate.exists():
            shutil.rmtree(candidate)
    print(f'Ready: {output / "html" / args.language / "index.html"}')


def build_full(args):
    manual = MANUALS[getattr(args, 'manual', 'tech')]
    output = (args.output or reader.DOCS / 'build' / manual.output).resolve()
    for manifest in (output / 'html').glob('*/build-manifest.json'):
        if json.loads(manifest.read_text()).get('manual', 'tech') != manual.key:
            raise RuntimeError(f'Output belongs to another manual: {output}. Choose a separate --output directory.')
    if args.refresh_languages:
        refresh_registry(output / 'html')
        print(f'Refreshed available languages: {output / "html/languages.json"}')
        return
    if args.language != 'en':
        build_translated(args, output)
        return
    print('Checking English guide sources, figures, and tools…', flush=True)
    inv = inventory(manual)
    preflight(inv, args.language)
    if args.check:
        print('Input checks passed. No files converted or built.')
        return
    output.mkdir(parents=True, exist_ok=True)
    signature = input_signature(inv)
    work = output / 'work'
    record = output / 'conversion-inputs.json'
    if args.reuse_conversion:
        if not record.is_file() or json.loads(record.read_text()) != signature:
            raise RuntimeError('Full-guide conversion is missing or stale; rerun without --reuse-conversion.')
        tex = work / 'Manuals' / manual.directory
        print('Reusing the completed LaTeX conversion.', flush=True)
    else:
        if work.exists():
            if not (work / '.fds-docs-work').is_file():
                raise RuntimeError(f'Refusing to replace an unrecognized staging directory: {work}')
            # Import successful pre-cache builds. Their figure format and
            # content-hashed names are identical to figure-cache-v1.
            if record.is_file() and json.loads(record.read_text()).get('conversion-rules') in {2, 3}:
                old_figures = work / 'Manuals' / manual.directory / 'WEB_FIGURES'
                cache = output / 'figure-cache-v1'
                cache.mkdir(exist_ok=True)
                for figure in old_figures.glob('*'):
                    if figure.is_file() and not (cache / figure.name).exists():
                        figure.replace(cache / figure.name)
                print('Reusing saved figure conversions when source hashes match.', flush=True)
            shutil.rmtree(work)
        record.unlink(missing_ok=True)
        tex = stage(inv, work)
        print(f'Converting the full guide. Log: {tex / "conversion.log"}', flush=True)
        try:
            reader.run(['make4ht', '-u', '-f', 'html5-common_domfilters',
                        '-e', str(reader.HERE / 'full-guide.mk4'),
                        'preview.tex', 'mathml,svg,fn-in,2'], tex, tex / 'conversion.log', progress=True)
        except Exception as error:
            log = (tex / 'conversion.log').read_text()
            lines = re.sub(r'\x1b\[[0-9;]*m', '', log).splitlines()
            detail = []
            for i, line in enumerate(lines):
                if '[ERROR]' in line and '.tex' in line:
                    detail.append(line + ('\n' + lines[i + 1] if i + 1 < len(lines) else ''))
                    break
            raise RuntimeError('Conversion failed.\n' + '\n'.join(detail)
                               + f'\nFull log: {tex / "conversion.log"}') from error
        log = (tex / 'conversion.log').read_text()
        if any(marker in log for marker in ['[ERROR]', 'Unbalanced Tag', 'specials ignored']):
            details = '\n'.join(line for line in re.sub(r'\x1b\[[0-9;]*m', '', log).splitlines()
                                if any(marker in line for marker in ['[ERROR]', '[WARNING]', 'specials ignored']))
            raise RuntimeError(f'Conversion diagnostics require attention:\n{details}\nFull log: {tex / "conversion.log"}')
        record.write_text(json.dumps(signature, indent=2))
        print('LaTeX conversion completed; checking generated page inventory.', flush=True)
    pages = discover_pages(tex, inv)
    print(f'Assembling English reader: {len(pages)} pages.', flush=True)
    candidate = Path(tempfile.mkdtemp(prefix='edition-', dir=output))
    try:
        pages, revision = reader.make_reader(tex, candidate, {}, pages=pages, full=True, manual=manual)
        from translate_preview import translation_units
        translation_inputs = output / 'translation-inputs'
        translation_inputs.mkdir(exist_ok=True)
        for position, page in enumerate(pages, 1):
            doc = html.parse(str(candidate / page['file']))
            units = translation_units(doc, include_tables=True)
            print(f'Preparing translation input {position}/{len(pages)}: '
                  f'{page["file"]} ({len(units)} prose blocks)', flush=True)
            (translation_inputs / (Path(page['file']).stem + '.json')).write_text(json.dumps([
                {'id': unit['id'], 'source_sha256': unit['sha256'], 'source': unit['source']}
                for unit in units], ensure_ascii=False, indent=2))
        counts = reader.validate(candidate, pages, progress=True)
        print(f'Validation passed: {len(pages)} pages, '
              f'{sum(p["math"] for p in counts.values())} math expressions.', flush=True)
        manifest = {'complete': True, 'manual': manual.key, 'entry': manual.entry, 'language': args.language, 'fds_revision': revision,
                    'pages': [page['file'] for page in pages], 'validation': counts,
                    'source_sha256': signature, 'pdf_sha256': reader.sha(reader.MANUALS / manual.directory / manual.pdf),
                    'limitations': ['Uses local figure inputs; release archive download is not yet integrated.',
                                    'Print covers remain in the PDF; generated HTML includes front matter and all content chapters.',
                                    'Full-guide visual review is required after the first build.']}
        (candidate / 'build-manifest.json').write_text(json.dumps(manifest, indent=2))
        # A static selector is useful even without JS; JS refreshes availability
        # from the shared registry after other language builds complete.
        for page in pages:
            path = candidate / page['file']
            doc = html.parse(str(path))
            picker = etree.Element('label', {'class': 'language-picker'})
            etree.SubElement(picker, 'span').text = 'Idioma' if args.language == 'es' else 'Language'
            select = etree.SubElement(picker, 'select', {'aria-label': 'Idioma' if args.language == 'es' else 'Language'})
            etree.SubElement(select, 'option', value=page['file'], selected='selected').text = LANGUAGES[args.language]
            doc.xpath('//div[@class="header-links"]')[0].insert(1, picker)
            path.write_text('<!doctype html>\n' + etree.tostring(doc, encoding='unicode', method='html'))
        print('Publishing English edition and updating the language dropdown…', flush=True)
        publish(candidate, output / 'html', args.language)
    finally:
        if candidate.exists():
            shutil.rmtree(candidate)
    print(f'Ready: {output / "html" / args.language / "index.html"}')
    print(f'Serve with: python3 -m http.server 8766 --bind 127.0.0.1 --directory {output / "html"}')
