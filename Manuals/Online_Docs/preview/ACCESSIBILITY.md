# Online manual accessibility

The engineering target is WCAG 2.2 AA. The Revised Section 508 Standards
incorporate WCAG 2.0 A/AA for public-facing electronic content; targeting a newer
WCAG version is not a compliance finding. See the Configuration Management
Plan's “Online Manuals and Accessibility” section for scope and current status.

## Build and audit

From `Manuals/Online_Docs`, with the FDS Python environment active:

```sh
# Install requirements from their one authoritative location.
(cd ../../.github && python -m pip install -r requirements.txt)
python -m playwright install chromium

./preview/build_tech_guide.sh --language en
./preview/build_tech_guide.sh --language es

# Static checks need only lxml; --browser also runs axe and keyboard/reflow checks.
python preview/audit_accessibility.py --site build/tech-guide/html/en --browser --report build/accessibility/en.json
python preview/audit_accessibility.py --site build/tech-guide/html/es --browser --report build/accessibility/es.json
```

The browser audit starts and stops its own loopback HTTP server. It does not
need the preview server or an API key. On an Ubuntu Actions runner, install
Chromium with `python -m playwright install --with-deps chromium`. Run the same
audit commands after building/downloading the edition artifacts. Upload the
JSON reports even on failure; do not put generated reports under version control.
No deployment or scheduled workflow is added by this accessibility pass.

The audit exits nonzero for structural errors, axe violations, keyboard/reflow
failures, or a failed browser run. A run without `--browser` explicitly records
`browser_status: "not run"`. Human-review items and axe `incomplete` findings
are retained in the report and must be resolved during conformance evaluation.
The axe engine version and its full results are recorded; its rule coverage is
not assumed to cover every criterion of WCAG 2.2 or Section 508.

## Implemented and remaining work

- Headings retain source hierarchy without missing levels. Anchors, numbering,
  and translation inputs are preserved.
- The reader provides named navigation, a keyboard skip link, modal search,
  visible focus, and hidden mobile navigation that cannot receive focus.
- Secondary text is darker; prose links are underlined. Wide equations and
  tables get keyboard-focusable scroll areas only when overflowing with JS
  enabled, with a focusable fallback when JS is disabled.
- Figure alternatives identify the figure and associate its full caption,
  including structured MathML. This is a baseline, not a reviewed description:
  every figure appears in the manual-review report. Improve insufficient
  descriptions through online accessibility metadata tied to source figures.
  This work does not retrofit the original LaTeX or PDF.
- Recognized simple data tables have headers and explicit cell associations.
  Unknown structures, including the nine spectral-band tables in Volume 1,
  remain flagged; no headers are guessed for them. Equation and tabbing layout
  tables use presentation semantics without hiding their content.
- Native MathML is retained. CSS-drawn overlines retain the approved visual
  treatment but need particular screen-reader evaluation; flattened spoken
  labels may not convey complex expressions faithfully.

Before claiming conformance, review all flagged figures/tables and exercise
VoiceOver with Safari and NVDA with Firefox, documenting actual versions.
Check equations 2.6, 2.9, 3.4–3.6, nested fractions, scripts and grouped overlines.
Test Tab/Shift+Tab, Enter, Escape, search, the language selector, mobile navigation,
and horizontally scrolling equations. Test 200% text enlargement and 400% browser
zoom, visible/unobscured focus, text spacing, and contrast in every interaction
state. The automated 320-CSS-pixel check is only a reflow approximation.

Check Spanish pronunciation and translated caption quality separately. Figure
artwork and the PDF remain English. PDF accessibility is outside this work; no LaTeX/PDF accessibility
retrofit is being made. Keep a dated assessment tied to
the exact source revision and edition; zero automated errors are not certification.

## Regression checks

```sh
python -m unittest discover -s preview -p 'test_*.py'
node preview/test_language_menu.cjs
node preview/test_accessibility.cjs
```

Initial verification on October 29, 2026 uses a separate reader fixture made from
the saved LaTeX conversion. The normal reuse guard detects a changed
`Bibliography/authors.tex`, so a full English rebuild is required before
publication. Live browser verification was blocked by the environment's approval
service; browser and assistive-technology results are not claimed as passed.
