"""Full Volume 1 orchestration; --check never runs TeX or writes build files."""
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

MAIN = 'FDS_Technical_Reference_Guide.tex'
LANGUAGES = {'en': 'English', 'es': 'Español'}
GRAPHIC = re.compile(r'(?:\\includegraphics\*?|\\addplot\+?\s+graphics)(?:\s*\[[^\]]*\])?\s*\{([^{}]+)\}', re.S)
SLUGS = {'Equation_Chapter': 'overview', 'Mass_Chapter': 'mass-transport',
         'Particle_Chapter': 'particles', 'Introduction_Chapter': 'introduction'}


def uncomment(text):
    # TeX comments consume the newline as well. Leaving it behind creates
    # paragraph breaks where commented-out rows occur inside align/tabular.
    return re.sub(r'(?<!\\)%[^\n]*(?:\n|$)', '', text)


def active_geometry(text, enabled):
    def select(match):
        branches = match[1].split(r'\else', 1)
        return branches[0] if enabled else (branches[1] if len(branches) == 2 else '')
    return re.sub(r'\\ifcompgeom\b(.*?)\\fi\b', select, text, flags=re.S)


def slug(name):
    return SLUGS.get(name, re.sub(r'[^a-z0-9]+', '-', name.replace('_Chapter', '').lower()).strip('-'))


def inventory():
    main = uncomment((reader.GUIDE / MAIN).read_text())
    common = uncomment((reader.MANUALS / 'Bibliography/commoncommands.tex').read_text())
    flags = re.findall(r'\\compgeom(true|false)\b', common)
    enabled = bool(flags and flags[-1] == 'true')
    main = active_geometry(main, enabled)
    names = re.findall(r'\\include\{([^{}]+)\}', main.split(r'\mainmatter', 1)[1])
    chapters = [name for name in names if name != 'Appendices']
    # Preserve the original front matter and bibliography/appendix order, but
    # leave print-only covers and logo pages in the linked PDF.
    preamble = main.split(r'\begin{document}', 1)[0]
    master = preamble + '\\begin{document}\n' + main[main.index(r'\frontmatter'):]
    sources = {MAIN: master}
    sources.update({name + '.tex': uncomment((reader.GUIDE / (name + '.tex')).read_text()) for name in names})
    appendices = re.findall(r'\\chapter\{', sources['Appendices.tex'])
    figures, missing, plot_graphics = {}, [], set()
    for source in sources.values():
        for match in GRAPHIC.finditer(source):
            reference = match[1]
            if match[0].startswith(r'\addplot'):
                plot_graphics.add(reference)
            if reference in figures:
                continue
            path = (reader.GUIDE / reference).resolve()
            candidates = [path] if path.suffix else [path.with_suffix(ext) for ext in ['.pdf', '.png', '.jpg', '.jpeg', '.eps']]
            found = next((p for p in candidates if p.is_file()), None)
            if found is None:
                missing.append(reference)
            elif found.suffix.lower() == '.eps':
                missing.append(reference + ' (EPS requires a PDF or PNG counterpart)')
            else:
                figures[reference] = found
    for reference in plot_graphics:
        if reference in figures and figures[reference].suffix.lower() != '.pdf':
            missing.append(reference + ' (PGFPlots background currently requires a PDF source)')
    return {'chapters': chapters, 'appendices': len(appendices), 'sources': sources,
            'figures': figures, 'plot_graphics': plot_graphics,
            'missing': sorted(set(missing)), 'complex_geometry': enabled}


def expected_pages(inv):
    return (['fds-developers', 'preface', 'disclaimer', 'acknowledgments', 'references']
            + [slug(name) for name in inv['chapters']]
            + ['appendix-' + chr(97 + i) for i in range(inv['appendices'])])


def cache_path(language, page):
    path = reader.HERE / 'translations' / language / (page + '.json')
    if not path.exists() and language == 'es' and page == 'overview':
        return reader.HERE / 'translations/es-overview.json'
    return path


def preflight(inv, language):
    errors = ['Missing figure: ' + name for name in inv['missing']]
    env_path = '/Library/TeX/texbin:/opt/homebrew/bin:' + os.environ.get('PATH', '')
    for tool in ['make4ht', 'biber', 'dvisvgm', 'kpsewhich', 'pdftocairo', 'pdftops']:
        if not shutil.which(tool, path=env_path):
            errors.append('Missing tool: ' + tool)
    if not (reader.GUIDE / 'FDS_Technical_Reference_Guide.pdf').is_file():
        errors.append('Build the local Technical Guide PDF for the PDF manual link.')
    print(f"Volume 1: {len(inv['chapters'])} chapters, {inv['appendices']} appendices, "
          f"{len(inv['figures'])} external figures; complex geometry {'enabled' if inv['complex_geometry'] else 'disabled'}.", flush=True)
    if errors:
        raise RuntimeError('Input checks failed:\n' + '\n'.join(errors))


def stage(inv, work):
    tex = work / 'Manuals' / reader.GUIDE.name
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
                      .replace('backref=true', 'backref=false'))
    folder = tex / 'WEB_FIGURES'
    folder.mkdir()
    replacements, web_figures = {}, {}
    for i, (reference, source) in enumerate(sorted(inv['figures'].items()), 1):
        print(f"Preparing figure {i}/{len(inv['figures'])}: {reference}", flush=True)
        stem = re.sub(r'[^A-Za-z0-9_-]', '-', source.stem) + '-' + reader.sha(source)[:12]
        if reference in inv['plot_graphics']:
            # PGFPlots backgrounds need real DVI graphics, not TeX4ht's BMP
            # placeholder specials (which silently disappear from the SVG).
            eps = folder / (stem + '.eps')
            reader.run(['pdftops', '-eps', str(source), str(eps)], tex)
            replacements[reference] = str(eps.relative_to(tex))
            continue
        if source.suffix.lower() == '.pdf':
            raster = folder / (stem + '.png')
            vector = folder / (stem + '.svg')
            reader.run(['pdftocairo', '-png', '-singlefile', '-r', '160', str(source), str(folder / stem)], tex)
            reader.run(['pdftocairo', '-svg', str(source), str(vector)], tex)
            web_figures[str(raster.relative_to(tex))] = str(vector.relative_to(tex))
        else:
            raster = folder / (stem + source.suffix.lower())
            shutil.copy2(source, raster)
        replacements[reference] = str(raster.relative_to(tex))
    for name, text in inv['sources'].items():
        text = GRAPHIC.sub(lambda m: m[0].replace('{' + m[1] + '}', '{' + replacements[m[1]] + '}'), text)
        text = reader.group_nested_primes(reader.group_math_scripts(text))
        (tex / ('preview.tex' if name == MAIN else name)).write_text(text)
    (tex / 'web-figures.json').write_text(json.dumps(web_figures, indent=2))
    (tex / 'preview.mk4').write_text('Make:htlatex {}\nMake:biber {}\nMake:htlatex {}\nMake:htlatex {}\n')
    return tex


def discover_pages(tex, inv):
    pages = []
    for path in tex.glob('preview*.html'):
        doc = html.parse(str(path), parser=html.HTMLParser(encoding='utf-8'))
        headings = doc.xpath('//body/h2')
        if not headings:
            if path.name == 'preview.html' and doc.xpath('//div[@class="tableofcontents"]'):
                pages.append({'raw': path.name, 'file': 'manual-index.html', 'source': MAIN,
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
        if number:
            number = int(number[1])
            source = inv['chapters'][number - 1] + '.tex'
            name, eyebrow = slug(Path(source).stem), 'CHAPTER ' + str(number)
        elif appendix:
            number = appendix[1]
            source, name, eyebrow = 'Appendices.tex', 'appendix-' + number.lower(), 'APPENDIX ' + number
        else:
            number, source, name, eyebrow = None, MAIN, slug(title), 'FRONT MATTER'
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


def translate_edition(site, pages, language, translator=None, cache_root=None):
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
        for heading in doc.xpath('//article//*[self::h1 or self::h3 or self::h4 or self::h5]'):
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
        localize_full_interface(doc)
        for link in doc.xpath('//article//div[@class="tableofcontents"]//a[@href]'):
            target = urlsplit(link.get('href'))
            label = heading_labels.get((target.path or name, target.fragment))
            if label is not None and not len(link):
                link.text = label
        for link in doc.xpath('//nav[@class="chapter-nav"]/a'):
            link.findall('span')[-1].text = titles.get(link.get('href'), link.findall('span')[-1].text)
        for link in doc.xpath('//nav[@class="page-navigation"]/a'):
            link.find('small').tail = titles.get(link.get('href'), link.find('small').tail)
        headings = {n.get('id'): reader.plain(n) for n in doc.xpath('//article//h3[@id] | //article//h4[@id]')}
        for link in doc.xpath('//aside[@class="page-toc"]//nav/a'):
            link.text = headings[link.get('href')[1:]]
        for figure in doc.xpath('//article//figure'):
            caption = figure.find('figcaption')
            if caption is not None:
                for image in figure.xpath('.//img'):
                    image.set('alt', reader.plain(caption))
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
            if node.tag in {'h1', 'h3', 'h4'}:
                if texts:
                    search.append({'title': title, 'chapter': titles[name], 'url': name + ('#' + anchor if anchor else ''), 'text': ' '.join(texts)})
                title, anchor, texts = reader.plain(node), node.get('id', ''), []
            elif isinstance(node.tag, str):
                texts.append(reader.plain(node))
        if texts:
            search.append({'title': title, 'chapter': titles[name], 'url': name + ('#' + anchor if anchor else ''), 'text': ' '.join(texts)})
    (site / 'search-index.json').write_text(json.dumps(search, ensure_ascii=False))
    return stats


def localize_full_interface(doc):
    labels = {
        '//a[@class="skip-link"]': 'Saltar al contenido',
        '//span[@class="brand-name"]/span': 'DOCUMENTACIÓN TÉCNICA',
        '//span[@class="edition"]': 'SEXTA EDICIÓN',
        '//div[@class="header-links"]/a': 'Manual PDF (inglés) ↗',
        '//div[@class="manual-name"]': 'Guía de referencia técnica',
        '//div[@class="volume"]': 'Volumen 1 · Modelo matemático',
        '//button[@class="search-trigger"]/span': '⌕ Buscar en este manual',
        '//div[@class="nav-heading"]': 'EN ESTE MANUAL',
        '//div[@class="sidebar-footer"]/p': 'Traducción automática. Las figuras y el PDF conservan el contenido original en inglés.',
        '//div[@class="sidebar-footer"]/a': 'Descargar el PDF en inglés ↓',
        '//div[@class="breadcrumb"]/a': 'Guía de referencia técnica',
        '//div[@class="breadcrumb"]/span[last()]': 'Modelo matemático',
        '//aside[@class="page-toc"]/div': 'EN ESTA PÁGINA',
        '//a[@class="back-top"]': 'Volver arriba ↑',
        '//label[@id="search-title"]': 'Buscar en el manual',
        '//div[@id="search-status"]': 'Buscar en este manual.',
        '//footer[@class="page-footer"]/span[1]': 'FDS · Guía de referencia técnica',
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
            'Previous and next pages': 'Páginas anterior y siguiente', 'Close search': 'Cerrar búsqueda'}
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
    (root / 'index.html').write_text('<!doctype html><meta charset="utf-8"><title>FDS Technical Guide</title>'
        + ''.join(f'<p><a href="{entry["code"]}/index.html">{entry["label"]}</a></p>' for entry in editions))


def input_signature(inv):
    signature = {name: reader.hashlib.sha256(text.encode()).hexdigest() for name, text in inv['sources'].items()}
    signature.update({str(path): reader.sha(path) for path in inv['figures'].values()})
    for directory in ['Bibliography', 'LaTeX_Style_Files']:
        signature.update({str(path): reader.sha(path) for path in (reader.MANUALS / directory).iterdir()
                          if path.is_file() and path.suffix in {'.tex', '.bib', '.sty', '.cls', '.bst'}})
    return signature


def english_edition(output):
    source = output / 'html/en'
    manifest = source / 'build-manifest.json'
    message = 'Build English first: ./preview/build_tech_guide.sh --language en'
    if not manifest.is_file():
        raise RuntimeError('No completed English edition found. ' + message)
    data = json.loads(manifest.read_text())
    names = data.get('pages', [])
    if (not data.get('complete') or data.get('language') != 'en' or not names
            or any(Path(name).name != name or not (source / name).is_file() for name in names)):
        raise RuntimeError('The English edition is incomplete. ' + message)
    return source, data, [{'file': name} for name in names]


def build_translated(args, output):
    from translate_argos import ArgosTranslator
    source, english, pages = english_edition(output)
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
        stats = translate_edition(candidate, pages, args.language, translator, output / 'translation-cache')
        counts = reader.validate(candidate, pages)
        # Fail safely if another job replaced English while translation ran.
        if source_hashes != {p['file']: reader.sha(source / p['file']) for p in pages}:
            raise RuntimeError('English changed during translation; rerun this language build.')
        (candidate / 'index.html').write_text('<!doctype html><html lang="' + args.language + '">'
            '<meta charset="utf-8"><meta http-equiv="refresh" content="0;url=overview.html">'
            '<title>FDS</title><a href="overview.html">' + LANGUAGES[args.language] + '</a></html>')
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
    output = (args.output or reader.DOCS / 'build/tech-guide').resolve()
    if args.refresh_languages:
        refresh_registry(output / 'html')
        print(f'Refreshed available languages: {output / "html/languages.json"}')
        return
    if args.language != 'en':
        build_translated(args, output)
        return
    print('Checking English guide sources, figures, and tools…', flush=True)
    inv = inventory()
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
        tex = work / 'Manuals' / reader.GUIDE.name
        print('Reusing the completed LaTeX conversion.', flush=True)
    else:
        if work.exists():
            if not (work / '.fds-docs-work').is_file():
                raise RuntimeError(f'Refusing to replace an unrecognized staging directory: {work}')
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
        pages, revision = reader.make_reader(tex, candidate, {}, pages=pages, full=True)
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
        manifest = {'complete': True, 'language': args.language, 'fds_revision': revision,
                    'pages': [page['file'] for page in pages], 'validation': counts,
                    'source_sha256': signature, 'pdf_sha256': reader.sha(reader.GUIDE / 'FDS_Technical_Reference_Guide.pdf'),
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
