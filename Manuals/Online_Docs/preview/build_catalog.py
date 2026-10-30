#!/usr/bin/env python3
"""Generate the shared manuals home page from completed editions only."""
import argparse
import hashlib
import html
import json
from pathlib import Path
import shutil
import tempfile

from manuals import MANUALS

HERE = Path(__file__).resolve().parent
LANGUAGES = {'en': 'English', 'es': 'Español'}
CATALOG_ORDER = ['tech', 'verification', 'validation', 'config', 'user']


def write_atomic(path, content):
    with tempfile.NamedTemporaryFile(mode='w', encoding='utf-8', dir=path.parent, delete=False) as stream:
        stream.write(content)
        temporary = Path(stream.name)
    temporary.replace(path)


def catalog_entries(root):
    entries = []
    for key in CATALOG_ORDER:
        manual = MANUALS[key]
        editions = []
        for code, label in LANGUAGES.items():
            folder = root / manual.output / 'html' / code
            try:
                proof = json.loads((folder / 'build-manifest.json').read_text())
            except (OSError, ValueError):
                continue
            pages = proof.get('pages', [])
            if (proof.get('complete') is True and proof.get('language') == code
                    and proof.get('manual', 'tech') == manual.key and pages
                    and (folder / 'index.html').is_file()
                    and all(Path(p).name == p and (folder / p).is_file() for p in pages)):
                editions.append({'code': code, 'label': label,
                                 'url': f'{manual.output}/html/{code}/index.html'})
        if editions:
            entries.append({'key': manual.key, 'title': manual.title,
                            'subtitle': manual.subtitle, 'title_es': manual.title_es,
                            'subtitle_es': manual.subtitle_es, 'editions': editions})
    return entries


def refresh_catalog(root, update_readers=False):
    root = Path(root).resolve()
    root.mkdir(parents=True, exist_ok=True)
    destination = root / 'index.html'
    marker = '<!-- FDS manuals catalog -->'
    if destination.exists() and marker not in destination.read_text():
        raise RuntimeError(f'Refusing to replace an unrelated home page: {destination}')
    entries = catalog_entries(root)
    groups = {'technical': [], 'user': []}
    for entry in entries:
        edition = next((e for e in entry['editions'] if e['code'] == 'en'), entry['editions'][0])
        link = f'<a href="{html.escape(edition["url"], quote=True)}">{html.escape(entry["title"])}</a>'
        group = 'user' if entry['key'] == 'user' else 'technical'
        groups[group].append(f'<section class="manual-card"><p class="volume">{html.escape(entry["subtitle"])}</p>'
                             f'<h3>{link}</h3></section>')
    cards = ''.join(f'<section class="publication"><h2>{heading}</h2><div class="manual-grid">{"".join(groups[group])}</div></section>'
                    for group, heading in [('technical', 'NIST SP 1018 · Technical Reference Guide'),
                                           ('user', 'NIST SP 1019 · User Guide')] if groups[group])
    content = (HERE / 'catalog.html').read_text().replace('{{CARDS}}', cards)
    available = {e['code'] for entry in entries for e in entry['editions']}
    options = ''.join(f'<option value="{code}">{label}</option>' for code, label in LANGUAGES.items() if code in available)
    content = content.replace('{{LANGUAGES}}', options or '<option value="en">English</option>')
    content = content.replace('{{EMPTY}}', '' if entries else '<p>No manuals have completed building yet.</p>')
    for name in ['catalog.css', 'catalog.js']:
        shutil.copy2(HERE / name, root / name)
        digest = hashlib.sha256((HERE / name).read_bytes()).hexdigest()[:12]
        content = content.replace(name + '"', name + '?v=' + digest + '"')
    # Known locations also discover converters that started before this feature.
    registry = {'kind': 'fds-manual-catalog', 'manuals': entries,
                'locations': [{'key': m.key, 'title': m.title, 'subtitle': m.subtitle,
                               'title_es': m.title_es, 'subtitle_es': m.subtitle_es,
                               'output': m.output} for m in [MANUALS[key] for key in CATALOG_ORDER]],
                'languages': LANGUAGES}
    write_atomic(root / 'manuals.json', json.dumps(registry, ensure_ascii=False, indent=2))
    write_atomic(destination, content)
    if update_readers:
        # Upgrade presentation assets without running TeX or translation.
        assets = {name: hashlib.sha256((HERE / name).read_bytes()).hexdigest()[:12]
                  for name in ['reader.js', 'reader.css']}
        from lxml import html as dom, etree
        for entry in entries:
            for edition in entry['editions']:
                folder = (root / edition['url']).parent
                for name in assets:
                    shutil.copy2(HERE / name, folder / 'assets' / name)
                for page in folder.glob('*.html'):
                    doc = dom.parse(str(page))
                    scripts = doc.xpath('//script[starts-with(@src,"assets/reader.js")]')
                    if not scripts:
                        continue
                    for script in scripts:
                        script.set('src', 'assets/reader.js?v=' + assets['reader.js'])
                    for style in doc.xpath('//link[starts-with(@href,"assets/reader.css")]'):
                        style.set('href', 'assets/reader.css?v=' + assets['reader.css'])
                    write_atomic(page, '<!doctype html>\n' + etree.tostring(doc, encoding='unicode', method='html'))
    return entries


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=HERE.parent / 'build')
    args = parser.parse_args()
    entries = refresh_catalog(args.root, update_readers=True)
    print(f'Manuals home: {args.root.resolve() / "index.html"} ({len(entries)} completed manuals)')
