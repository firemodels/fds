#!/usr/bin/env python3
"""Convert FDS LaTeX with make4ht into a reader; use --full for all of Volume 1.

This pilot requires an already built local Tech Guide (PDF and .aux files) for
links to chapters outside the preview. It does not fetch release assets.
"""
from __future__ import annotations

import argparse
import hashlib
import html
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time
from urllib.parse import quote, unquote, urlsplit

from accessibility import heading_levels, adapt_content

from lxml import etree
from lxml import html as dom

HERE = Path(__file__).resolve().parent
DOCS = HERE.parent
MANUALS = DOCS.parent
GUIDE = MANUALS / 'FDS_Technical_Reference_Guide'
CHAPTERS = [('Equation_Chapter', 2, 'overview'),
            ('Mass_Chapter', 3, 'mass-transport'),
            ('Particle_Chapter', 9, 'particles')]


def run(args, cwd, log=None, progress=False):
    env = os.environ.copy()
    env['PATH'] = '/Library/TeX/texbin:/opt/homebrew/bin:' + env['PATH']
    env['TEXINPUTS'] = '.:../LaTeX_Style_Files:'
    env['XDG_CACHE_HOME'] = str(Path(cwd) / '.cache')
    # MacTeX does not discover Homebrew's Ghostscript library automatically.
    if 'LIBGS' not in env and Path('/opt/homebrew/lib/libgs.dylib').exists():
        env['LIBGS'] = '/opt/homebrew/lib/libgs.dylib'
    if log:
        with log.open('w') as stream:
            if not progress:
                subprocess.run(args, cwd=cwd, env=env, stdout=stream,
                               stderr=subprocess.STDOUT, check=True)
            else:
                started = last_update = time.monotonic()
                # Keep the complete log on disk while showing selected live
                # milestones. File polling also works when TeX is quiet.
                with log.open() as updates, subprocess.Popen(
                        args, cwd=cwd, env=env, stdout=stream,
                        stderr=subprocess.STDOUT) as process:
                    while True:
                        finished = process.poll() is not None
                        for line in updates:
                            if '[FDS-DOCS]' in line:
                                print(line.split('[FDS-DOCS]', 1)[1].strip(), flush=True)
                                last_update = time.monotonic()
                        now = time.monotonic()
                        if finished:
                            break
                        if now - last_update >= 30:
                            print(f'  Conversion running: {int(now - started)} seconds elapsed; '
                                  'details in conversion.log', flush=True)
                            last_update = now
                        time.sleep(0.25)
                    if process.returncode:
                        raise subprocess.CalledProcessError(process.returncode, args)
    else:
        return subprocess.check_output(args, cwd=cwd, env=env, text=True).strip()


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def group(text, pos):
    """Read a brace-balanced TeX argument (used only for .aux metadata)."""
    if text[pos] != '{':
        raise ValueError(f'Expected brace at {pos}')
    start, depth = pos + 1, 1
    for i in range(start, len(text)):
        if text[i] in '{}' and (i == 0 or text[i - 1] != '\\'):
            depth += 1 if text[i] == '{' else -1
            if depth == 0:
                return text[start:i], i + 1
    raise ValueError('Unclosed TeX argument')


def group_math_scripts(text):
    """Give TeX4ht an explicit group around font/text commands in scripts."""
    pattern = re.compile(r'([_^])\s*\\(?:text(?:rm|sf|tt|bf|it|normal)?|math(?:rm|sf|tt|bf|it|cal|bb|normal)|mbox|hbox)\s*(?=\{)')
    parts, pos = [], 0
    while match := pattern.search(text, pos):
        _, end = group(text, match.end())
        parts.append(text[pos:match.start()])
        parts.append(match[1] + '{' + group_math_scripts(text[match.start() + 1:end].strip()) + '}')
        pos = end
    parts.append(text[pos:])
    return ''.join(parts)


def group_nested_primes(text):
    """Avoid TeX4ht recursion when primes follow a primed, grouped subscript."""
    parts, pos = [], 0
    for match in re.finditer(r'(?<!\\)_\s*(?=\{)', text):
        if match.start() < pos:
            continue
        argument, end = group(text, match.end())
        primes = re.match(r"'+", text[end:])
        if "'" not in argument or not primes:
            continue
        stop = end + len(primes[0])
        # TeX combines x'^n into one superscript; leave that distinct case
        # alone rather than accidentally creating a second superscript.
        if text[stop:].lstrip().startswith('^'):
            continue
        parts.append(text[pos:end])
        parts.append('^{' + r'\prime' * len(primes[0]) + '}')
        pos = stop
    parts.append(text[pos:])
    return ''.join(parts)


def read_labels():
    result = {}
    for path in GUIDE.glob('*.aux'):
        for line in path.read_text().splitlines():
            if not line.startswith(r'\newlabel{'):
                continue
            key, end = group(line, len(r'\newlabel'))
            value, _ = group(line, end)
            fields, pos = [], 0
            while pos < len(value) and value[pos] == '{':
                field, pos = group(value, pos)
                fields.append(field)
            if len(fields) >= 4:
                result[key] = {'number': fields[0], 'destination': fields[3]}
    return result


def stage_inputs(work):
    labels = read_labels()
    full_pdf = GUIDE / 'FDS_Technical_Reference_Guide.pdf'
    if not labels or not full_pdf.exists():
        raise RuntimeError('Build the local Tech Guide PDF first; this pilot uses its .aux reference metadata.')
    sources = {name: (GUIDE / f'{name}.tex').read_text() for name, _, _ in CHAPTERS}
    local_labels = set(re.findall(r'\\label\{([^}]+)\}', '\n'.join(sources.values())))
    fallbacks = {}
    staged_manuals = work / 'Manuals'
    for directory in ['Bibliography', 'LaTeX_Style_Files']:
        dst = staged_manuals / directory
        dst.mkdir(parents=True, exist_ok=True)
        for path in (MANUALS / directory).iterdir():
            if path.is_file() and path.suffix in {'.tex', '.bib', '.sty', '.cls', '.bst'}:
                shutil.copy2(path, dst / path.name)
    common = staged_manuals / 'Bibliography/commoncommands.tex'
    # TeX4ht must select its own drivers. Only the staging copy is adapted.
    common.write_text(common.read_text().replace(r'\usepackage[pdftex]{graphicx}',
                                                 r'\usepackage{graphicx}')
                      .replace(r'\usepackage[pdftex,', r'\usepackage[')
                      .replace(r'\usepackage{colortbl}', '% Table colors supplied by reader CSS.')
                      .replace('backref=true', 'backref=false'))
    tex = staged_manuals / GUIDE.name
    tex.mkdir(parents=True, exist_ok=True)

    def reference(match):
        command, key = match.groups()
        if key in local_labels:
            return match.group(0)
        if key not in labels:
            raise RuntimeError(f'Unknown reference: {key}')
        target = labels[key]
        fallbacks[key] = target
        link = (r'\href{full-guide.pdf\#nameddest=' + target['destination']
                + '}{' + target['number'] + '}')
        return '(' + link + ')' if command == 'eqref' else link

    for name, text in sources.items():
        text = re.sub(r'\\(eqref|ref)\{([^}]+)\}', reference, text)
        # TeX4ht's MathML text handler needs explicit subscript grouping.
        text = group_nested_primes(group_math_scripts(text))
        (tex / f'{name}.tex').write_text(text)

    figures = tex / 'SCRIPT_FIGURES'
    figures.mkdir(exist_ok=True)
    plot = GUIDE / 'SCRIPT_FIGURES/particle_size_distribution.pdf'
    if not plot.exists():
        raise RuntimeError(f'Missing generated figure: {plot}')
    shutil.copy2(plot, figures / plot.name)
    # PNG supplies TeX dimensions; the final reader uses the vector SVG.
    run(['pdftocairo', '-png', '-singlefile', '-r', '160', str(plot),
         str(figures / plot.stem)], tex)
    run(['pdftocairo', '-svg', str(plot), str(figures / (plot.stem + '.svg'))], tex)
    particle = tex / 'Particle_Chapter.tex'
    particle.write_text(particle.read_text().replace(
        '{SCRIPT_FIGURES/particle_size_distribution}',
        '{SCRIPT_FIGURES/particle_size_distribution.png}'))
    body = '\n'.join(r'\setcounter{chapter}{' + str(number - 1) + '}\n'
                     + r'\input{' + name + '}' for name, number, _ in CHAPTERS)
    (tex / 'preview.tex').write_text(r'''\documentclass[11pt]{book}
\input{../Bibliography/commoncommands}
\addbibresource{../Bibliography/FDS_refs.bib}
\addbibresource{../Bibliography/FDS_general.bib}
\addbibresource{../Bibliography/FDS_mathcomp.bib}
\begin{document}
''' + body + '\n' + r'''\printbibliography[heading=bibintoc,title={References}]
\end{document}
''')
    (tex / 'preview.mk4').write_text(
        'Make:htlatex {}\nMake:biber {}\nMake:htlatex {}\nMake:htlatex {}\n')
    return tex, fallbacks


def plain(element):
    return ' '.join(''.join(element.itertext()).split())


def normalize_math(body):
    """Preserve differential operators and clearly separated math accents."""
    # MathML Core does not inherit legacy mstyle mathvariant onto single-letter
    # mi tokens. Preserve explicit LaTeX roman type on the token itself (SP 811,
    # 10.2.3), including differential d/D, units, and descriptive subscripts.
    for symbol in body.xpath('.//math//mi[not(@mathvariant)]'):
        styles = symbol.xpath('ancestor::*[@mathvariant]')
        if styles and styles[-1].get('mathvariant') == 'normal':
            symbol.set('mathvariant', 'normal')

    # Audited unmarked integral measures in 2.15 and 3.29: d followed by a
    # primed integration variable at the end of the integral's row. Never
    # globally romanize d: it also denotes nozzle diameter in Chapter 9.
    for symbol in body.xpath('.//math//mi[text()="d"][not(@mathvariant)]'):
        previous = symbol.getprevious()
        following = symbol.getnext()
        if (previous is not None and previous.tag == 'mspace'
                and following is not None and following.tag == 'msup'
                and following.getnext() is None
                and any('∫' in plain(n) for n in symbol.itersiblings(preceding=True))):
            symbol.set('mathvariant', 'normal')

    for nabla in body.xpath('.//math//*[self::mi or self::mo][text()="∇"]'):
        following = [node for node in nabla.itersiblings() if isinstance(node.tag, str)]
        application, negative_spaces, operator = [], [], None
        for node in following:
            if node.tag == 'mo' and plain(node) == '\u2061':
                # Nabla is not a named function such as sin. The inserted
                # application token also damages gradients and Laplacians.
                application.append(node)
            elif node.tag == 'mspace' and node.get('class') == 'negthinspace':
                negative_spaces.append(node)
            else:
                if node.tag == 'mo' and plain(node) in {'⋅', '×'}:
                    operator = node
                break
        for node in application:
            node.getparent().remove(node)
        nabla.tag = 'mo'
        nabla.attrib.pop('class', None)
        nabla.set('lspace', '0em')
        nabla.set('rspace', '0em')
        if operator is not None:
            if negative_spaces:
                for node in negative_spaces:
                    node.getparent().remove(node)
                # A TeX binary operator gets 4mu; \! removes 3mu on its left.
                operator.set('lspace', '0.0556em')
            elif operator.get('lspace') is None:
                operator.set('lspace', '0.2222em')
            operator.set('rspace', '0.2222em')

    for symbol in body.xpath('.//math//mi[text()="ṁ"]'):
        # A precomposed Latin extended character may use an upright fallback
        # glyph. A native accent over the mathematical m preserves italics.
        accent = etree.Element('mover', accent='true')
        base = etree.SubElement(accent, 'mi', dict(symbol.attrib))
        base.text = 'm'
        mark = etree.SubElement(accent, 'mo', accent='true', stretchy='false')
        mark.text = '˙'
        accent.tail = symbol.tail
        symbol.getparent().replace(symbol, accent)

    for accent in body.xpath('.//math//mover'):
        children = [node for node in accent if isinstance(node.tag, str)]
        if len(children) != 2:
            continue
        mark = children[1]
        # Remove the earlier lift, including when normalizing an existing tree.
        # Wrapping an accent in mpadded changes its layout as well as shifting
        # it, so reducing voffset alone did not fix the floating dots in 3.6.
        if mark.tag == 'mpadded' and mark.get('class') == 'dot-accent':
            wrapper = mark
            mark = wrapper[0]
            mark.tail = wrapper.tail
            accent.replace(wrapper, mark)
        if mark.tag != 'mo' or plain(mark) not in {'˙', '\u0307'}:
            continue
        # A direct accent token lets the font supply placement and attachment.
        accent.set('accent', 'true')
        mark.set('accent', 'true')
        mark.set('stretchy', 'false')
        base = children[0]
        if (base.tag == 'mrow' and not base.attrib and len(base) == 1
                and not (base.text or '').strip() and not (base[0].tail or '').strip()):
            symbol = base[0]
            symbol.tail = base.tail
            accent.replace(base, symbol)

    # A stretchy accent still rendered as a short glyph in the user's browser.
    # Draw every LaTeX \overline as a rule spanning its entire MathML base.
    # Process inner lines first so nested overlines retain both rules.
    for accent in reversed(body.xpath('.//math//mover')):
        children = [node for node in accent if isinstance(node.tag, str)]
        if (len(children) == 2
                and 'mml-overline' in accent.get('class', '').split()
                and children[1].tag == 'mo' and plain(children[1]) in {'¯', '‾', '\u0305'}):
            overline = etree.Element('mrow', dict(accent.attrib))
            overline.attrib.pop('accent', None)
            base_text = plain(children[0])
            overline.set('class', 'mixture-overline' if base_text == 'W' else 'latex-overline')
            overline.set('aria-label', base_text + ' with overline')
            overline.append(children[0])
            overline.tail = accent.tail
            accent.getparent().replace(accent, overline)


def reader_link(href, current, filename_map, id_files):
    """Resolve local anchors against their actual generated page, including citations."""
    url = urlsplit(href)
    if url.scheme or url.netloc:
        return href
    path, sep, fragment = href.partition('#')
    if path == 'preview.html' and fragment in filename_map:
        return filename_map[fragment]
    if not path or path in filename_map:
        destinations = id_files.get(fragment, [])
        if len(destinations) == 1:
            return (destinations[0] if destinations[0] != current else '') + '#' + fragment
        if path in filename_map:
            return filename_map[path] + (sep + fragment if sep else '')
    return href


def make_reader(tex, site, fallbacks, pages=None, full=False):
    site.mkdir(parents=True, exist_ok=True)
    assets = site / 'assets'
    assets.mkdir(exist_ok=True)
    for name in ['reader.css', 'reader.js']:
        shutil.copy2(HERE / name, assets / name)
    shutil.copy2(tex / 'preview.css', assets / 'tex4ht.css')
    font = Path(run(['kpsewhich', 'latinmodern-math.otf'], tex))
    if font.is_file():
        shutil.copy2(font, assets / 'latinmodern-math.otf')
    shutil.copy2(GUIDE / 'FDS_Technical_Reference_Guide.pdf', site / 'full-guide.pdf')
    filename_map = {'preview.html': 'index.html'}
    if pages is None:
        pages = [{'raw': f'previewch{i}.html', 'file': slug + '.html',
                  'source': name + '.tex', 'number': number}
                 for i, (name, number, slug) in enumerate(CHAPTERS, 1)]
        pages.append({'raw': 'previewli1.html', 'file': 'references.html', 'source': None, 'number': None})
    filename_map.update({page['raw']: page['file'] for page in pages})
    figure_map_path = tex / 'web-figures.json'
    figure_map = json.loads(figure_map_path.read_text()) if figure_map_path.exists() else {}

    id_files = {}
    for position, page in enumerate(pages, 1):
        if full:
            print(f'Reading page {position}/{len(pages)}: {page["file"]}', flush=True)
        tree = dom.parse(str(tex / page['raw']), parser=dom.HTMLParser(encoding='utf-8'))
        body = tree.find('body')
        for element in body.xpath('.//div[@class="crosslinks"] | .//div[@class="crosslinks-bottom"]'):
            element.drop_tree()
        headings = body.xpath('.//h2')
        if not headings and page.get('synthetic_title'):
            heading = etree.Element('h2', id='manual-index')
            heading.text = page['synthetic_title']
            body.insert(0, heading)
        else:
            heading = headings[0]
        marker = heading.find('span')
        if marker is not None:
            marker.drop_tree()
        for br in heading.findall('br'):
            br.drop_tree()
        heading.tag = 'h1'
        heading_levels(body)
        page['title'] = plain(heading)
        page['body'] = body
        page['toc'] = [(element.get('id'), plain(element), element.tag == 'h3')
                       for element in body.xpath('.//h2[@id] | .//h3[@id]')]
        for element in body.xpath('.//*[@id]'):
            # Shared empty bibliography helper anchors are not destinations.
            if element.get('id') == 'X0-':
                element.attrib.pop('id')
            else:
                id_files.setdefault(element.get('id'), []).append(page['file'])

    for page in pages:
        body = page['body']
        normalize_math(body)
        for link in body.iter():
            key = next((key for key in ['href', 'xlink:href', '{http://www.w3.org/1999/xlink}href']
                        if link.get(key) is not None), None)
            if key is None:
                continue
            href = link.get(key)
            link.set(key, reader_link(href, page['file'], filename_map, id_files))
            if 'full-guide.pdf' in href:
                link.set('title', 'Open this reference in the full PDF manual')
                link.set('class', 'pdf-reference')
        for image in body.xpath('.//img'):
            original = tex / image.get('src')
            if image.get('src') in figure_map:
                original = tex / figure_map[image.get('src')]
            elif original.stem == 'particle_size_distribution':
                original = original.with_suffix('.svg')
            if not original.is_file():
                raise RuntimeError(f'Missing converted image: {original}')
            shutil.copy2(original, assets / original.name)
            image.set('src', 'assets/' + original.name)
            image.set('loading', 'lazy')
            image.set('alt', 'Technical diagram')
        adapt_content(body)

    revision = run(['git', 'rev-parse', '--short', 'HEAD'], MANUALS)
    template = (HERE / 'reader.html').read_text()
    if full:
        template = template.replace('<body>', '<body data-language-manifest="../languages.json">')
        template = template.replace('IN THIS PREVIEW', 'IN THIS MANUAL').replace(
            'Selected chapters from the<br>FDS Technical Reference Guide.',
            'Complete Volume 1: Mathematical Model.')
        template = template.replace('Search the chapters included in this preview.', 'Search this manual.')
    # Refresh presentation assets automatically when their contents change.
    for name in ['reader.css', 'reader.js']:
        template = template.replace(f'assets/{name}', f'assets/{name}?v={sha(HERE / name)[:12]}')
    search = []
    for position, page in enumerate(pages):
        if full:
            print(f'Writing page {position + 1}/{len(pages)}: {page["file"]}', flush=True)
        nav = ''.join(
            f'<a class="chapter-link {"active" if p is page else ""}" '
            f'{"aria-current=\"page\"" if p is page else ""} href="{p["file"]}">'
            f'<span class="chapter-number">{p["number"] or ""}</span>'
            f'<span>{html.escape(p["title"])}</span></a>' for p in pages)
        toc = ''.join(f'<a class="{"subsection" if sub else ""}" href="#{anchor}">{html.escape(title)}</a>'
                      for anchor, title, sub in page['toc'])
        content = ''.join(etree.tostring(c, method='html', encoding='unicode') for c in page['body'])
        prev = pages[position - 1] if position > 0 else None
        nxt = pages[position + 1] if position + 1 < len(pages) else None
        pager = ''.join(f'<a href="{p["file"]}" class="{direction}"><small>{label}</small>{html.escape(p["title"])}</a>'
                        for p, direction, label in [(prev, 'previous', 'Previous'), (nxt, 'next', 'Next')]
                        if p is not None)
        if page['source']:
            source = assets / Path(page['source']).name
            shutil.copy2(GUIDE / page['source'], source)
            source_link = f'<a href="assets/{source.name}">View LaTeX source ↗</a>'
        else:
            source_link = ''
        values = {'TITLE': html.escape(page['title']), 'NAV': nav, 'TOC': toc,
                  'CONTENT': content, 'PAGER': pager, 'REVISION': revision,
                  'EYEBROW': page.get('eyebrow', 'CHAPTER ' + str(page['number']) if page['number'] else 'BIBLIOGRAPHY'),
                  'SOURCE': source_link}
        output = template
        for key, value in values.items():
            output = output.replace('{{' + key + '}}', value)
        (site / page['file']).write_text(output)

        title, anchor, texts = page['title'], '', []
        for element in page['body']:
            if element.tag in {'h1', 'h2', 'h3'}:
                if texts:
                    search.append({'title': title, 'chapter': page['title'],
                                   'url': page['file'] + ('#' + anchor if anchor else ''),
                                   'text': ' '.join(texts)})
                title, anchor, texts = plain(element), element.get('id', ''), []
            elif isinstance(element.tag, str):
                texts.append(plain(element))
        if texts:
            search.append({'title': title, 'chapter': page['title'],
                           'url': page['file'] + ('#' + anchor if anchor else ''), 'text': ' '.join(texts)})
    (site / 'index.html').write_text('<!doctype html><html lang="en"><meta charset="utf-8">'
        '<meta http-equiv="refresh" content="0;url=overview.html"><title>FDS Technical Guide</title>'
        '<a href="overview.html">Open the FDS Technical Guide</a></html>')
    (site / 'search-index.json').write_text(json.dumps(search, ensure_ascii=False))
    return pages, revision


def validate(site, pages, progress=False):
    docs = {p['file']: dom.parse(str(site / p['file'])) for p in pages}
    ids = {name: set(doc.xpath('//@id')) for name, doc in docs.items()}
    errors = []
    for position, (name, doc) in enumerate(docs.items(), 1):
        if progress:
            print(f'Validating page {position}/{len(docs)}: {name}', flush=True)
        # Validate before HTML's error recovery can conceal malformed math.
        for math in re.findall(r'<math\b.*?</math>', (site / name).read_text(), re.S):
            try:
                node = etree.fromstring(math.encode())
                if node.xpath('.//math'):
                    errors.append(f'{name}: nested math elements')
            except etree.XMLSyntaxError as error:
                errors.append(f'{name}: malformed MathML: {error}')
        for src in doc.xpath('//img/@src | //script/@src | //link[@rel="stylesheet"]/@href'):
            if not (site / urlsplit(src).path).is_file():
                errors.append(f'{name}: missing asset {src}')
        for href in doc.xpath('//a/@href'):
            url = urlsplit(href)
            if url.scheme or url.netloc:
                continue
            target = unquote(url.path) or name
            if not (site / target).is_file():
                errors.append(f'{name}: missing page {href}')
            if url.fragment and target in ids and unquote(url.fragment) not in ids[target]:
                errors.append(f'{name}: missing anchor {href}')
        if '??' in ''.join(doc.xpath('//article//text()')):
            errors.append(f'{name}: unresolved LaTeX reference')
    if errors:
        raise RuntimeError('\n'.join(errors))
    return {name: {'math': len(doc.xpath('//math')), 'figures': len(doc.xpath('//article//figure')),
                   'equations': len(doc.xpath('//td[@class="eq-no"]'))} for name, doc in docs.items()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, help='Build root; defaults depend on --full')
    parser.add_argument('--full', action='store_true', help='Build all of Volume 1, front matter, references, and appendices')
    parser.add_argument('--language', default='en', choices=['en', 'es'], help='Edition to build (default: en)')
    parser.add_argument('--check', action='store_true', help='Check full-guide inputs without converting or building')
    parser.add_argument('--refresh-languages', action='store_true', help='Refresh the full-site language registry after combining edition artifacts')
    parser.add_argument('--reuse-conversion', action='store_true', help='Restyle previously converted HTML')
    args = parser.parse_args()
    if args.full:
        from full_guide import build_full
        build_full(args)
        return
    if args.check or args.refresh_languages:
        parser.error('--check and --refresh-languages require --full')
    output = (args.output or DOCS / 'build/tech-preview').resolve()
    if args.reuse_conversion:
        tex = output / 'work/Manuals' / GUIDE.name
        fallbacks = json.loads((output / 'pdf-references.json').read_text())
    else:
        tex, fallbacks = stage_inputs(output / 'work')
        (output / 'pdf-references.json').write_text(json.dumps(fallbacks, indent=2))
        print('Converting original LaTeX chapters with make4ht…', flush=True)
        run(['make4ht', '-u', '-f', 'html5', 'preview.tex', 'mathml,svg,fn-in,2'],
            tex, tex / 'conversion.log')
        log = (tex / 'conversion.log').read_text()
        if any(marker in log for marker in ['[ERROR]', 'Unbalanced Tag', 'specials ignored']):
            raise RuntimeError(f'Conversion errors: inspect {tex / "conversion.log"}')
    pages, revision = make_reader(tex, output / 'html', fallbacks)
    counts = validate(output / 'html', pages)
    manifest = {'fds_revision': revision, 'converter': run(['make4ht', '--version'], tex),
                'chapters': [x[0] for x in CHAPTERS], 'validation': counts,
                'source_sha256': {f'{name}.tex': sha(GUIDE / f'{name}.tex') for name, _, _ in CHAPTERS},
                'figure_sha256': sha(GUIDE / 'SCRIPT_FIGURES/particle_size_distribution.pdf'),
                'pdf_sha256': sha(GUIDE / 'FDS_Technical_Reference_Guide.pdf'),
                'external_references': fallbacks,
                'limitations': ['Uses locally populated generated figures, not a freshly downloaded release.',
                                'References outside these chapters open the existing local PDF.',
                                'PDF auxiliary metadata must correspond to the local source.',
                                'Bibliography numbering follows citation order within this preview.',
                                'Sphinx/Pandoc comparison pending dependency installation.']}
    (output / 'html/build-manifest.json').write_text(json.dumps(manifest, indent=2))
    if args.language == 'es':
        print(run([sys.executable, str(HERE / 'translate_preview.py'),
                   '--site', str(output / 'html')], HERE))
    print(json.dumps(counts, indent=2))
    print(f'Preview ready: {output / "html"}')


if __name__ == '__main__':
    try:
        main()
    except (RuntimeError, subprocess.CalledProcessError) as error:
        print(f'Build failed: {error}', file=sys.stderr)
        sys.exit(1)
