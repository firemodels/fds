"""Offline Argos adapter and resumable, model-specific translation cache."""
from __future__ import annotations

import hashlib
from importlib.metadata import version
import json
from pathlib import Path
import re
import sys

TOKEN = re.compile(r'⟦(\d+)⟧')

# These prose tokens are scientific data, not words for the model to rewrite.
# MathML, code, links, and markup already have protected placeholders.
LITERALS = re.compile(
    r'https?://[^\s]+|[\w.+-]+@[\w.-]+\.[A-Za-z]+|'
    r'\b[A-Z][A-Z0-9_]{1,}\b|'
    r'(?<!\w)[+-]?\d+(?:[.,]\d+)*(?:[eE][+-]?\d+)?'
    r'(?:\s*(?:kg|mg|km|cm|mm|kW|MW|kJ|MJ|Pa|kPa|MPa|Hz|mol|m|s|K|W|J|g|°C|°F|%)(?!\w))?'
)
POLICY = 'protected-prose-v1'


def digest(text):
    return hashlib.sha256(text.encode()).hexdigest()


def translate_protected(source, translate):
    """Try whole sentences; fall back to spans if the model alters placeholders."""
    originals = {}
    next_id = max([int(x) for x in TOKEN.findall(source)] + [-1]) + 1

    def protect(match):
        nonlocal next_id
        token = f'⟦{next_id}⟧'
        next_id += 1
        originals[token] = match[0]
        return token

    # Never match the digits inside existing markup placeholders.
    chunks = re.split(r'(⟦\d+⟧)', source)
    prepared = ''.join(part if TOKEN.fullmatch(part) else LITERALS.sub(protect, part) for part in chunks)
    expected = TOKEN.findall(prepared)
    result = translate(prepared) if re.search('[A-Za-z]', TOKEN.sub('', prepared)) else prepared
    fallback = TOKEN.findall(result) != expected or ('⟦' in TOKEN.sub('', result) or '⟧' in TOKEN.sub('', result))
    if fallback:
        translated = []
        for part in re.split(r'(⟦\d+⟧)', prepared):
            if TOKEN.fullmatch(part) or not re.search('[A-Za-z]', part):
                translated.append(part)
            else:
                leading = part[:len(part) - len(part.lstrip())]
                trailing = part[len(part.rstrip()):]
                translated.append(leading + translate(part.strip()).strip() + trailing)
        result = ''.join(translated)
    if TOKEN.findall(result) != expected:
        raise RuntimeError('Argos changed protected math/code placeholders; translation stopped.')
    for token, literal in originals.items():
        result = result.replace(token, literal)
    if not result.strip():
        raise RuntimeError('Argos returned an empty translation.')
    return result, fallback


class ArgosTranslator:
    def __init__(self, language):
        try:
            import argostranslate.package as package
            self.engine_version = version('argostranslate')
        except ImportError as error:
            raise RuntimeError(f'Argos is not installed for {sys.executable}. '
                               'Activate fds_python_env and install argostranslate.') from error
        packages = [p for p in package.get_installed_packages()
                    if p.from_code == 'en' and p.to_code == language]
        if len(packages) != 1:
            raise RuntimeError(f'Install exactly one direct Argos en → {language} model first; '
                               f'found {len(packages)}. The build does not download models.')
        self.package = packages[0]
        if not (Path(self.package.package_path) / 'stanza/resources.json').is_file():
            raise RuntimeError('This offline adapter requires a model with its bundled Stanza tokenizer. '
                               'No tokenizer will be downloaded during the build.')
        self.language = language
        self.engine = None
        # Hash actual model data, not just a mutable package version label.
        checksum = hashlib.sha256()
        for path in sorted(Path(self.package.package_path).rglob('*')):
            if path.is_file():
                checksum.update(str(path.relative_to(self.package.package_path)).encode())
                with path.open('rb') as stream:
                    for chunk in iter(lambda: stream.read(1024 * 1024), b''):
                        checksum.update(chunk)
        self.provenance = {
            'provider': 'argos', 'argos_version': self.engine_version,
            'model_version': str(self.package.package_version), 'model_sha256': checksum.hexdigest(),
            'source_language': 'en', 'target_language': language, 'policy': POLICY,
            'runtime': {name: version(name) for name in ['ctranslate2', 'sentencepiece', 'stanza']},
            'device': 'cpu',
        }

    def translate(self, text):
        if self.engine is None:
            import argostranslate.translate as translation
            import argostranslate.settings as settings
            settings.device = 'cpu'
            self.engine = translation.PackageTranslation(
                translation.Language('en', 'English'),
                translation.Language(self.language, self.language), self.package)
            # Older packages include Stanza models. Its default initialization
            # downloads metadata on every start; explicitly use the bundled data.
            if hasattr(self.engine.sentencizer, 'stanza_pipeline'):
                import stanza
                self.engine.sentencizer.stanza_pipeline = stanza.Pipeline(
                    lang=self.engine.sentencizer.stanza_lang_code,
                    dir=str(self.package.package_path / 'stanza'), processors='tokenize',
                    use_gpu=False, download_method=None, logging_level='WARNING')
        return self.engine.translate(text)


def generate_cache(units, path, translator):
    """Checkpoint each completed block; match text hashes across changed block IDs."""
    path = Path(path)
    cached = json.loads(path.read_text()) if path.is_file() else {}
    reusable = {}
    if cached.get('provenance') == translator.provenance:
        reusable = {b['source_sha256']: b for b in cached.get('blocks', [])}
    blocks, translated, reused, fallbacks = [], 0, 0, 0
    path.parent.mkdir(parents=True, exist_ok=True)
    for index, unit in enumerate(units, 1):
        previous = reusable.get(unit['sha256'])
        if (previous and previous.get('source') == unit['source']
                and TOKEN.findall(previous['translation']) == TOKEN.findall(unit['source'])):
            value, fallback = previous['translation'], previous.get('span_fallback', False)
            reused += 1
        else:
            value, fallback = translate_protected(unit['source'], translator.translate)
            translated += 1
        fallbacks += int(fallback)
        block = {'id': unit['id'], 'source_sha256': unit['sha256'], 'source': unit['source'],
                 'translation': value, 'span_fallback': fallback}
        blocks.append(block)
        reusable[unit['sha256']] = block
        # Retain later old blocks if interrupted partway through a refreshed page.
        checkpoint = {**reusable, **{b['source_sha256']: b for b in blocks}}
        data = {'provenance': translator.provenance, 'blocks': list(checkpoint.values())}
        temporary = path.with_suffix('.json.tmp')
        temporary.write_text(json.dumps(data, ensure_ascii=False, indent=2) + '\n')
        temporary.replace(path)
        if index % 10 == 0 or index == len(units):
            print(f'  {path.stem}: {index}/{len(units)} blocks ({translated} new, {reused} cached)', flush=True)
    return blocks, {'translated': translated, 'cached': reused, 'span_fallbacks': fallbacks}
