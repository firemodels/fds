"""Inventory the manuals' actual input graphs before starting a conversion."""
import re
from pathlib import Path


GRAPHIC = re.compile(r'(?:\\includegraphics\*?|\\incgraphics|\\addplot\+?\s+graphics)(?:\s*\[[^\]]*\])?\s*\{([^{}]+)\}', re.S)
INPUT = re.compile(r'\\(include|input)\s*\{([^{}]+)\}')
SLUGS = {'Equation_Chapter': 'overview', 'Mass_Chapter': 'mass-transport',
         'Particle_Chapter': 'particles', 'Introduction_Chapter': 'introduction'}


def uncomment(text):
    return re.sub(r'(?<!\\)%[^\n]*(?:\n|$)', '', text)


def active_geometry(text, enabled):
    def select(match):
        branches = match[1].split(r'\else', 1)
        return branches[0] if enabled else (branches[1] if len(branches) == 2 else '')
    return re.sub(r'\\ifcompgeom\b(.*?)\\fi\b', select, text, flags=re.S)


def slug(name):
    return SLUGS.get(name, re.sub(r'[^a-z0-9]+', '-', name.replace('_Chapter', '').lower()).strip('-'))


def alphabetic(index):
    result = ''
    while index:
        index, digit = divmod(index - 1, 26)
        result = chr(97 + digit) + result
    return result


def inventory(manual, manuals_root):
    guide = manuals_root / manual.directory
    common = uncomment((manuals_root / 'Bibliography/commoncommands.tex').read_text())
    flags = re.findall(r'\\compgeom(true|false)\b', common)
    enabled = bool(flags and flags[-1] == 'true')
    clean = lambda text: active_geometry(uncomment(text), enabled)
    main = clean((guide / manual.main).read_text())
    if r'\frontmatter' not in main:
        raise RuntimeError(f'Cannot locate front matter in {manual.main}')
    preamble = main.split(r'\begin{document}', 1)[0]
    # A developer's temporary includeonly selection must not truncate a full
    # online edition. Disabled/commented chapters and geometry flags still apply.
    preamble = re.sub(r'\\includeonly\s*\{[^{}]*\}', '', preamble)
    master = preamble + '\\begin{document}\n' + main[main.index(r'\frontmatter'):]
    sources, missing, visiting = {manual.main: master}, [], set()

    def collect(name):
        if name in visiting:
            raise RuntimeError(f'Cyclic LaTeX input: {name}')
        visiting.add(name)
        for match in INPUT.finditer(sources[name]):
            reference = match[2]
            if '#' in reference or '\\' in reference:
                raise RuntimeError(f'Unsupported dynamic LaTeX input in {name}: {reference}')
            path = (guide / reference).resolve()
            if not path.suffix:
                path = path.with_suffix('.tex')
            # TeX resolves inputs from its working directory, not the directory
            # of the currently included file. Preserve that layout in staging.
            if path.parent in {(manuals_root / 'Bibliography').resolve(),
                               (manuals_root / 'LaTeX_Style_Files').resolve()}:
                if not path.is_file() and path.name != 'gitrevision.tex':
                    missing.append('LaTeX input: ' + reference)
                continue
            try:
                key = str(path.relative_to(guide.resolve()))
            except ValueError as error:
                raise RuntimeError(f'LaTeX input outside {manual.directory}: {reference}') from error
            if key in visiting:
                raise RuntimeError(f'Cyclic LaTeX input: {reference}')
            if key not in sources:
                if not path.is_file():
                    missing.append('LaTeX input: ' + reference)
                    continue
                sources[key] = clean(path.read_text())
                collect(key)
        visiting.remove(name)

    collect(manual.main)
    chapters, appendices, chapter_sources, appendix_sources, part_sources = [], [], [], [], []
    state = 'front'

    def headings(name):
        nonlocal state
        # Headings in literal examples are not document structure.
        text = re.sub(r'\\begin\{(lstlisting|verbatim)\}.*?\\end\{\1\}', '', sources[name], flags=re.S)
        tokens = re.compile(r'\\(mainmatter|appendix)\b|\\(?:chapter|part)\s*\{|\\(include|input)\s*\{([^{}]+)\}')
        for token in tokens.finditer(text):
            if token[1]:
                state = token[1]
            elif token[2]:
                path = (guide / token[3]).resolve()
                if not path.suffix:
                    path = path.with_suffix('.tex')
                try:
                    key = str(path.relative_to(guide.resolve()))
                except ValueError:
                    continue
                if key in sources:
                    headings(key)
            elif token[0].startswith(r'\part'):
                part_sources.append(name)
            elif state == 'mainmatter':
                number = len(chapters) + 1
                identifier = Path(name).stem if name != manual.main else f'chapter-{number}'
                # Chapter files normally have one chapter. Give any additional
                # chapters distinct, stable filenames rather than overwriting.
                if identifier in chapters:
                    identifier = f'chapter-{number}'
                chapters.append(identifier)
                chapter_sources.append(name)
            elif state == 'appendix':
                appendices.append('appendix-' + alphabetic(len(appendices) + 1))
                appendix_sources.append(name)

    headings(manual.main)
    figures, plot_graphics = {}, set()
    for source in sources.values():
        # TeX picture overlays need real DVI artwork, like PGFPlots backgrounds.
        for picture in re.finditer(r'\\begin\{picture\}.*?\\end\{picture\}', source, re.S):
            plot_graphics.update(m[1] for m in GRAPHIC.finditer(picture[0]))
        for match in GRAPHIC.finditer(source):
            reference = match[1]
            # The User Guide's incgraphics wrapper forwards its argument to
            # includegraphics. Its concrete invocations are inventoried above.
            if reference == '#1':
                continue
            if '#' in reference or '\\' in reference:
                missing.append('Unresolved graphic argument: ' + reference)
                continue
            if match[0].startswith(r'\addplot'):
                plot_graphics.add(reference)
            if reference in figures:
                continue
            path = (guide / reference).resolve()
            candidates = [path] if path.suffix else [path.with_suffix(ext) for ext in ('.pdf', '.png', '.jpg', '.jpeg', '.eps')]
            found = next((candidate for candidate in candidates if candidate.is_file()), None)
            if found is None:
                missing.append('Figure: ' + reference)
            elif found.suffix.lower() == '.eps':
                missing.append('Figure requires PDF or PNG counterpart: ' + reference)
            else:
                figures[reference] = found
    return {'manual': manual, 'chapters': chapters, 'chapter_sources': chapter_sources,
            'part_sources': part_sources,
            'appendices': len(appendices), 'appendix_slugs': appendices, 'appendix_sources': appendix_sources,
            'sources': sources, 'figures': figures, 'plot_graphics': plot_graphics,
            'missing': sorted(set(missing)), 'complex_geometry': enabled}
