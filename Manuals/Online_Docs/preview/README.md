# Local Technical Guide preview

For accessibility scope, audit commands, and remaining human review, see
[Accessibility](ACCESSIBILITY.md). Accessibility is under development; these
manuals are not currently represented as fully Section 508 compliant.

## Build the complete Technical Guide (Volume 1)

From `Manuals/Online_Docs`, run:

```sh
./preview/build_tech_guide.sh --language en
python3 -m http.server 8766 --bind 127.0.0.1 --directory build/tech-guide/html
```

Open **http://127.0.0.1:8766/en/overview.html**. Port 8766 keeps the original
three-chapter preview at 8765 available for comparison. The launcher selects a
Python interpreter with lxml; `FDS_DOCS_PYTHON` can override it. The equivalent
Python command is `python3 preview/build_preview.py --full --language en`.

The full mode follows the source book's active chapter order and geometry flag,
including all 13 current chapters, all 18 appendices, developers, preface,
disclaimer, acknowledgments, and references. Print-only cover/logo pages stay
in the PDF. It stages all 18 referenced external figures, including the images
under `fig/fds/Geometry_Figures` and the two PGFPlots background images, converts PDF figures for browser use, and
retains the original LaTeX inputs. Generated figures must already exist locally;
release-archive downloading is still future work.
PGFPlots backgrounds use `pdftops -eps` during staging so their artwork survives
the DVI-to-SVG conversion. Other PDF figures use `pdftocairo` as before.

Useful commands:

```sh
# Check English tools, sources, and figure availability only:
./preview/build_tech_guide.sh --language en --check
# Restyle an unchanged, successfully converted full guide:
./preview/build_tech_guide.sh --language en --reuse-conversion
```

The English guide has been built and structurally validated. Visual review of
the full manual is ongoing.
Conversion diagnostics are written to
`build/tech-guide/work/Manuals/FDS_Technical_Reference_Guide/conversion.log`.
The English build reports each LaTeX pass and bibliography stage, emits an
elapsed-time update every 30 seconds during quiet conversion steps, then reports
page counts for reader assembly, translation input extraction, and validation.
The complete converter output remains in the log.
No successful edition is replaced until page inventory, references, and assets
pass validation. A failure leaves the existing served edition intact.

## Languages and Actions

`--language en` is the default. `--language es` selects a Spanish edition;
each successful full build is published under `html/<language>/`. The dropdown
fetches `html/languages.json` and offers only successfully built editions that
contain the current page. Switching preserves the section anchor. The current
edition remains usable if that registry cannot load.

Every non-English build **requires a completed English edition first**. It
copies `build/tech-guide/html/en`, translates the prose, and validates the new
edition. It never runs LaTeX, regenerates figures, or changes the English edition.
The translation uses that English edition's recorded source revision.

Activate your FDS Python environment, which must contain both `lxml` and
`argostranslate`. Install the direct English-to-Spanish Argos model once:

```sh
(cd ../../.github && python -m pip install -r requirements.txt)
python - <<'PY'
import argostranslate.package as package
package.update_package_index()
model = next(p for p in package.get_available_packages()
             if p.from_code == 'en' and p.to_code == 'es')
package.install_from_path(model.download())
PY
```

Then, from `Manuals/Online_Docs`:

```sh
./preview/build_tech_guide.sh --language es --check
./preview/build_tech_guide.sh --language es
```

The existing server will serve `http://127.0.0.1:8766/es/overview.html` after
`Ready:` appears. Refresh the English page to see Español in the dropdown.
The active virtual environment takes priority; `FDS_DOCS_PYTHON` remains an
explicit interpreter override. Do not reinstall from a particular source folder:
Python packages belong to the environment and Argos models to its package store.

No API key is needed. The build uses installed model data offline, including
Stanza's bundled sentence tokenizer. Translation caches live under
`build/tech-guide/translation-cache/es/`. They record model checksums, engine
versions, and source hashes. Blocks are checkpointed so interrupted builds can
resume; changed prose/model versions invalidate the relevant cache entries.
Math, code, link targets, numbers, and uppercase FDS keywords are protected.
If Argos changes placeholders, the adapter translates the intervening prose
spans separately and restores the protected items. This can reduce fluency;
technical terminology and grammar require review. Figure artwork and the PDF
remain English. The old assistant-generated Overview cache is not used here.

For Actions, build English once and upload the complete `html/en/` directory.
Each language job downloads it into `build/tech-guide/html/en/`, installs a
pinned Argos runtime and model, restores its translation cache, and runs the same
`--language` command. No TeX installation is needed in translation jobs. Pin the
model archive checksum rather than selecting the latest model each cron run.
Save caches per language and retain an immutable copy for reproducible releases;
Actions caches can be evicted. Currently `en` and `es` are implemented; another
language also needs translated reader interface labels.

Upload each successful `html/<language>/` edition and combine the artifacts into
one `build/tech-guide/html` directory. Then, in a single assembly/deployment job:

```sh
./preview/build_tech_guide.sh --refresh-languages
```

This regenerates the registry from the editions actually present. It does not
convert the manual or require translation caches. Do not run multiple builds
against the same work/output directory concurrently. Deployment and cron setup
are not created by these scripts.

## Three-chapter preview

This pilot converts the original Equation, Mass, and Particle chapters with
make4ht/TeX4ht, then adds a styled reader with chapter navigation, section links,
and a local search index. All chapter text is generated from the LaTeX source.

Requirements: Python 3.12+ with `lxml` (see [the FDS requirements](../../../.github/requirements.txt)), TeX Live with
make4ht, Biber, dvisvgm, Ghostscript, and Poppler's `pdftocairo`. A locally built
Technical Guide PDF, its `.aux` files, and populated `SCRIPT_FIGURES` are required.
The PDF and auxiliary files must correspond to the local source.

From `Manuals/Online_Docs`:

```sh
python3 preview/build_preview.py
python3 -m http.server 8765 --bind 127.0.0.1 --directory build/tech-preview/html
```

Open http://127.0.0.1:8765/overview.html. Use `--reuse-conversion` to rebuild the
reader after changing its CSS, HTML template, or JavaScript without rerunning TeX.

The initial run used the bundled Python interpreter at
`/Users/rmcdermo/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3`.

The build adapts only staging copies: it removes PDF-specific drivers and the
colortbl package, braces text subscripts for MathML, and links references outside
the preview to the full PDF. Table styling is supplied by CSS. On this Mac,
`LIBGS` points dvisvgm to the installed Homebrew Ghostscript library so legacy
picture diagrams retain their geometry. The particle plot is converted from its
existing generated PDF into SVG. Math uses native MathML and a local math font.

The converter rejects reported TeX errors and missing diagram content, and the
reader build checks local links, fragment targets, and assets. Details and input
hashes are recorded in `build/tech-preview/html/build-manifest.json`.

With `--language es`, the preview build also assembles a Spanish Overview pilot at `es/overview.html`, with
an English/Español selector that preserves the current section on Overview.
The selector appears on every preview page; from other chapters, the Spanish
option is labeled “Español (Overview)” and opens the translated Overview.
Navigation,
captions, and the chapter search index are translated. Other chapters, figure
artwork, bibliography entries, and the PDF remain English and are labeled.
The Spanish cache in `translations/es-overview.json` was generated by the
assistant for this trial; no external translation API is connected. It records
provenance, a terminology glossary, and source hashes for all 91 prose blocks.
The build rejects changed prose or missing/reordered protected placeholders.
MathML, numbers in equations, and anchor IDs must remain unchanged. The cache
is a generated translation artifact, not another manually authored manual.

To prepare a future translation refresh, run
`python3 preview/translate_preview.py --extract` and send the extracted prose
to the chosen translator with instructions to preserve every placeholder.
Review the updated translations, store their source hashes and provenance in
the cache, then rebuild. This legacy pilot is independent of the full-guide
Argos integration above. The deployment workflow itself has not been added.

Scientific typography follows NIST SP 811 (2008), especially §10.2.3.
The reader explicitly preserves Roman differential d/D in native MathML,
while keeping quantity variables italic. Run the regression checks with
`python3 -m unittest discover -s preview -p 'test_*.py'`.

Limitations: this is a three-chapter preview, uses existing local generated
figures rather than downloading a release archive, and renumbers bibliography
entries by first use in the preview. Browser visual/interaction checks were not
performed because the browser approval service failed. The Sphinx/Pandoc
comparison is pending dependency installation; its approval service also failed.
