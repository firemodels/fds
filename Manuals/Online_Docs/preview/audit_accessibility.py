#!/usr/bin/env python3
"""Audit generated editions; write reports under ignored build/.

Static checks always run. --browser adds axe, keyboard and reflow checks using
an isolated loopback server and Chromium. Neither mode certifies conformance.
"""
import argparse
from contextlib import contextmanager
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
import json
from pathlib import Path
from threading import Thread
from urllib.parse import quote

from lxml import html


def text(node):
    return ' '.join(''.join(node.itertext()).split())


def inspect_document(doc, filename):
    errors, review = [], []

    def error(rule, detail):
        errors.append({'page': filename, 'rule': rule, 'detail': detail})

    if not doc.getroot().get('lang'):
        error('language', 'Missing document language')
    if not doc.xpath('//head/title[normalize-space()]'):
        error('title', 'Missing page title')
    if not doc.xpath('//article'):
        return errors, review  # Redirect landing pages are checked by axe.
    if len(doc.xpath('//main')) != 1 or len(doc.xpath('//article//h1')) != 1:
        error('landmarks', 'Expected one main landmark and one chapter heading')
    if not doc.xpath('//a[@class="skip-link"][@href="#main-content"]'):
        error('skip-link', 'Missing skip-to-main link')
    level = 0
    for heading in doc.xpath('//article//*[self::h1 or self::h2 or self::h3 or self::h4 or self::h5 or self::h6]'):
        current = int(heading.tag[1])
        if current > level + 1:
            error('heading-order', f'{heading.tag}: {text(heading)[:100]}')
        level = current
    for node in doc.xpath('//*[@aria-describedby or @aria-labelledby or @headers]'):
        for attr in ('aria-describedby', 'aria-labelledby', 'headers'):
            for target in node.get(attr, '').split():
                if len(doc.xpath('//*[@id=$target]', target=target)) != 1:
                    error('label-target', f'{attr} references missing or ambiguous {target}')
    for image in doc.xpath('//article//img'):
        alt = image.get('alt', '').strip()
        if not alt or alt == 'Technical diagram':
            error('image-alternative', image.get('src'))
        review.append({'page': filename, 'kind': 'figure-description',
                       'target': image.get('src'), 'alt': alt,
                       'description': ' '.join(text(n) for target in image.get('aria-describedby', '').split()
                                               for n in doc.xpath('//*[@id=$target]', target=target)),
                       'task': ('This uncaptioned plot has only a source-name label. Supply a description of its essential information for the online edition.'
                                if image.get('data-description-review') == 'source-label-only' else
                                'Check that the caption conveys the figure’s essential information; record any supplemental description needed for the online edition.')})
    for table in doc.xpath('//article//table'):
        if table.get('role') == 'presentation':
            continue
        if table.get('data-table-review') or not table.xpath('.//th'):
            review.append({'page': filename, 'kind': 'table-headers', 'target': table.get('id'),
                           'task': 'Review complex header relationships; do not infer them from visual alignment.'})
        elif not table.xpath('.//th[@scope or @headers]'):
            error('table-headers', table.get('id'))
    overlines = doc.xpath('//math//*[@class="mixture-overline" or @class="latex-overline"]')
    if overlines:
        review.append({'page': filename, 'kind': 'math-speech', 'count': len(overlines),
                       'task': 'Check native MathML navigation and spoken overlines, scripts and accents with a screen reader. CSS-drawn overlines require particular attention.'})
    return errors, review


class QuietHandler(SimpleHTTPRequestHandler):
    def log_message(self, *args):
        pass


@contextmanager
def serve(root):
    server = ThreadingHTTPServer(('127.0.0.1', 0), partial(QuietHandler, directory=str(root)))
    thread = Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        yield f'http://127.0.0.1:{server.server_port}'
    finally:
        server.shutdown()
        server.server_close()
        thread.join()


def browser_checks(root, pages):
    from playwright.sync_api import sync_playwright
    from axe_playwright_python.sync_playwright import Axe
    results = []
    with serve(root) as origin, sync_playwright() as playwright:
        browser = playwright.chromium.launch()
        page = browser.new_page(viewport={'width': 1440, 'height': 1000})
        axe = Axe()
        for position, file in enumerate(pages, 1):
            path = file.relative_to(root).as_posix()
            print(f'Browser audit {position}/{len(pages)}: {path}', flush=True)
            page.goto(origin + '/' + quote(path), wait_until='networkidle')
            page.evaluate('document.fonts.ready')
            # Run the bundled engine's default rules, not a guessed WCAG tag
            # list; retain the engine version and all needs-review findings.
            result = axe.run(page).response
            results.append({'page': path, 'viewport': 'desktop', 'axe': result})
            if not page.locator('article').count():
                continue
            failures = []
            page.keyboard.press('Tab')
            if not page.locator('.skip-link').evaluate('(e) => e === document.activeElement'):
                failures.append('First Tab did not reach skip link')
            page.keyboard.press('Enter')
            if not page.locator('main').evaluate('(e) => e === document.activeElement'):
                failures.append('Skip link did not focus main')
            # Keyboard activation and native modal focus containment.
            page.locator('.search-trigger').focus()
            page.keyboard.press('Enter')
            if not page.locator('#search-input').evaluate('(e) => e === document.activeElement'):
                failures.append('Search did not focus input')
            for _ in range(4):
                page.keyboard.press('Tab')
                if not page.locator('#search-dialog').evaluate('(e) => e.contains(document.activeElement)'):
                    failures.append('Focus escaped search dialog')
            page.keyboard.press('Escape')
            if not page.locator('.search-trigger').evaluate('(e) => e === document.activeElement'):
                failures.append('Search did not restore focus')
            # 320 CSS px approximates the reflow width at 400% on 1280px.
            # This is not a substitute for actual browser/text zoom testing.
            page.set_viewport_size({'width': 320, 'height': 800})
            page.wait_for_timeout(250)
            if not page.locator('.sidebar').evaluate('(e) => e.inert'):
                failures.append('Closed mobile navigation is not inert')
            page.locator('.menu-toggle').focus()
            page.keyboard.press('Enter')
            if page.locator('.sidebar').evaluate('(e) => e.inert'):
                failures.append('Open mobile navigation remains inert')
            page.keyboard.press('Tab')
            # Header links precede the sidebar, so explicitly focus a sidebar
            # link before checking Escape returns to its toggle.
            page.locator('.search-trigger').focus()
            page.keyboard.press('Escape')
            if not page.locator('.menu-toggle').evaluate('(e) => e === document.activeElement'):
                failures.append('Mobile Escape did not restore focus')
            if page.evaluate('document.documentElement.scrollWidth > innerWidth + 1'):
                failures.append('Page overflows horizontally at 320 CSS px')
            results.append({'page': path, 'viewport': 'mobile', 'axe': axe.run(page).response,
                            'interaction_failures': failures})
            page.set_viewport_size({'width': 1440, 'height': 1000})
        browser.close()
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--site', type=Path, default=Path(__file__).resolve().parents[1] / 'build/tech-guide/html/en')
    parser.add_argument('--report', type=Path, default=Path(__file__).resolve().parents[1] / 'build/accessibility/report.json')
    parser.add_argument('--browser', action='store_true', help='Also run Chromium/axe and keyboard/reflow checks')
    args = parser.parse_args()
    root = args.site.resolve()
    pages = sorted(root.rglob('*.html'))
    if not pages:
        parser.error(f'No built HTML in {root}')
    report = {'site': str(root), 'pages': len(pages), 'errors': [], 'manual_review': [],
              'browser_status': 'not run', 'browser_results': [],
              'limitations': ['Automated checks do not establish WCAG or Section 508 conformance.',
                              'Manually test VoiceOver/Safari, NVDA/Firefox, keyboard, text zoom and 400% browser zoom.',
                              'Review translated descriptions and language pronunciation. The linked PDF is a separate accessibility deliverable.']}
    for index, page in enumerate(pages, 1):
        name = page.relative_to(root).as_posix()
        print(f'Static audit {index}/{len(pages)}: {name}', flush=True)
        doc = html.parse(str(page), parser=html.HTMLParser(encoding='utf-8'))
        errors, review = inspect_document(doc, name)
        report['errors'].extend(errors)
        report['manual_review'].extend(review)
    if args.browser:
        try:
            report['browser_results'] = browser_checks(root, pages)
            report['browser_status'] = 'completed'
        except Exception as error:
            report['browser_status'] = 'failed'
            report['errors'].append({'rule': 'browser-run', 'detail': str(error)})
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, ensure_ascii=False, indent=2) + '\n')
    violations = sum(len(r['axe'].get('violations', [])) + len(r.get('interaction_failures', []))
                     for r in report['browser_results'])
    print(f'{len(report["errors"])} static/run errors; '
          f'{len(report["manual_review"])} manual review items. Browser: {report["browser_status"]}.')
    if report["browser_status"] == "completed":
        print(f"Browser violations and interaction failures: {violations}")
    print(f'Report: {args.report.resolve()}')
    return 1 if report['errors'] or violations else 0


if __name__ == '__main__':
    raise SystemExit(main())
