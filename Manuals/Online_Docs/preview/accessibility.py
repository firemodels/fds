"""Shared, conservative accessibility adaptations for TeX4ht output.

Do not infer scientific descriptions or complex table headers from appearance.
The separate audit records those cases for human review.
"""
import re
import unicodedata
from lxml import etree


def text(node):
    return ' '.join(''.join(node.itertext()).split())


def heading_levels(body):
    # make_reader already promotes the chapter heading to h1. Compress absent
    # levels while retaining the source hierarchy (some appendices begin with
    # subsubsections). Keep classes, numbering, anchors and visual styling.
    stack = [(1, 1)]
    for node in body.xpath('.//*[self::h3 or self::h4 or self::h5 or self::h6]'):
        original = int(node.tag[1])
        while stack[-1][0] >= original:
            stack.pop()
        level = stack[-1][1] + 1
        node.tag = 'h' + str(level)
        stack.append((original, level))


def figure_labels(body, language='en'):
    """Associate the actual caption, retaining its navigable MathML.

    Captions are a baseline alternative, not proof of an adequate description
    of a plot. Keep that distinction machine readable for the audit report.
    """
    for index, figure in enumerate(body.xpath('.//figure'), 1):
        captions = figure.xpath('.//figcaption')
        if not captions:
            continue
        caption = captions[0]
        if not caption.get('id'):
            caption.set('id', f'a11y-caption-{index}')
        for part, image in enumerate(figure.xpath('.//img'), 1):
            label = re.match(r'(?:Figure|Figura)\s+[^:]+', text(caption))
            name = label[0] if label else ('Figura' if language == 'es' else 'Figure')
            if len(figure.xpath('.//img')) > 1:
                name += f', {"imagen" if language == "es" else "image"} {part}'
            image.set('alt', name)
            image.set('aria-describedby', caption.get('id'))
            image.set('data-description-review', 'caption-only')
    object_image_labels(body, language)


def object_image_labels(body, language='en'):
    """Use inspected Smokeview object-table labels, not inferred artwork meaning."""
    canonical = lambda node: unicodedata.normalize('NFKC', text(node))
    for index, table in enumerate(body.xpath('.//table'), 1):
        rows = table.xpath('./tr | ./tbody/tr')
        if not any(canonical(row).startswith('SMOKEVIEW_ID') for row in rows):
            continue
        captions = table.xpath('./caption')
        if not captions:
            continue
        caption = captions[0]
        if not caption.get('id'):
            caption.set('id', f'a11y-object-table-{index}')
        # The inspected dual-frame table explicitly labels these columns.
        states = None
        for row in rows:
            cells = row.xpath('./td | ./th')
            labels = [canonical(cell) for cell in cells]
            if row.get('data-object-states') == 'true' or labels == ['', 'inactive', 'active']:
                row.set('data-object-states', 'true')
                states = labels
        for row in rows:
            cells = row.xpath('./td | ./th')
            if not cells or not cells[0].xpath('.//code') or not row.xpath('.//img'):
                continue
            identifier = canonical(cells[0])
            for column, cell in enumerate(cells):
                for image in cell.xpath('.//img'):
                    label = ('Objeto Smokeview ' if language == 'es' else 'Smokeview object ') + identifier
                    if states and column < len(states) and states[column]:
                        label += ', ' + states[column]
                    image.set('alt', label)
                    image.set('aria-describedby', caption.get('id'))
                    image.set('data-description-review', 'table-label-only')


def table_semantics(body):
    for table in body.xpath('.//table'):
        classes = table.get('class', '').split()
        if any(c in classes for c in ('equation', 'equation-star', 'tabbing')):
            table.set('role', 'presentation')
            continue
        if 'tabular' not in classes:
            continue
        rows = table.xpath('./tr | ./tbody/tr')
        data = []
        for row in rows:
            if row.get('class') in {'hline', 'cline'} and not text(row):
                row.set('aria-hidden', 'true')
            else:
                data.append(row)
        if not data:
            continue
        cells = data[0].xpath('./td | ./th')
        # These header signatures were inspected in the Technical Guide.
        # Leave spectral-band schematics and future unknown tables for review.
        labels = [text(cell) for cell in cells]
        simple = (labels[:1] in [['Band #'], ['Flux Limiter'], ['Geometry'], ['Liquid']]
                  or labels == ['', 'Air', 'Fuel', 'Products'])
        if not simple:
            table.set('data-table-review', 'header-associations')
            continue
        columns = []
        for position, cell in enumerate(cells, 1):
            if not text(cell):
                columns.append(None)
                continue
            cell.tag = 'th'
            cell.set('scope', 'colgroup' if int(cell.get('colspan', '1')) > 1 else 'col')
            if not cell.get('id'):
                cell.set('id', f'{table.get("id", "table")}-header-{position}')
            columns.extend([cell.get('id')] * int(cell.get('colspan', '1')))
        for row in data[1:]:
            column = 0
            for cell in row.xpath('./td | ./th'):
                width = int(cell.get('colspan', '1'))
                headers = list(dict.fromkeys(h for h in columns[column:column + width] if h))
                if headers:
                    cell.set('headers', ' '.join(headers))
                column += width


def scroll_regions(body):
    # Do not create nested scroll areas, or hundreds of unnecessary tab stops.
    # JS makes only overflowing regions keyboard focusable; tabindex is the
    # fallback when scripts are disabled.
    for node in list(body.xpath('.//table | .//math[@display="block"]')):
        if node.xpath('ancestor::table | ancestor::math'):
            continue
        wrapper = etree.Element('div', {'class': 'table-scroll', 'tabindex': '0',
                                       'role': 'region', 'aria-label': 'Scrollable equation or table'})
        node.addprevious(wrapper)
        wrapper.append(node)


def adapt_content(body):
    figure_labels(body)
    table_semantics(body)
    scroll_regions(body)
