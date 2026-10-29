# Online manual conversion rules

- Keep all pip requirements in the repository's `.github/requirements.txt`.
  Do not create a separate requirements file under the online documentation.
- Never track generated documentation, translation caches, converted figures,
  or other build products. Track source inputs, scripts, templates, and tests.
  Keep generated files ignored and exclude them from every staging command,
  including the legacy `preview/translations/es-overview.json` cache.

- Always use NIST SP 811 (2008) as the authority for scientific typography and
  SI notation in this project, especially sections 10.2.1–10.2.4 (printed
  pp. 34–35). Reference: https://doi.org/10.6028/NIST.SP.811e2008;
  user-supplied copy: /Users/rmcdermo/Papers/nistspecialpublication811e2008.pdf.
  Differential operators d (and the material-derivative D) are upright Roman;
  quantity symbols such as mass m and diameter d are italic. Units and
  descriptive subscripts are upright; variable subscripts are italic. Preserve
  these semantic distinctions during conversion and translation, not through
  global letter replacements. In native MathML, put mathvariant="normal" on
  the mi token: a legacy mstyle wrapper alone does not reliably enforce it.

- Keep scientific content in the original LaTeX. Fix conversion or presentation
  code, not generated HTML or a second maintained copy of a chapter.
- TeX4ht staging must brace font/text commands used directly after `_` or `^`
  (for example, `y_\mathrm{CO}` becomes `y_{\mathrm{CO}}`). Preserve nested
  arguments and avoid changing already grouped scripts. TeX comments consume
  their terminating newline; removing only their text introduces invalid blank
  paragraphs inside `align`. Figure discovery must cover PGFPlots
  `\addplot graphics` as well as `\includegraphics`.
  For PGFPlots PDF backgrounds, stage EPS via `pdftops -eps` so dvisvgm
  embeds the artwork. Raster placeholders can compile with bounding boxes
  yet silently disappear from the generated SVG; verify embedded images.
  Primes following a primed subscript can exhaust TeX4ht's grouping stack
  (Appendix P: `\dot{q}_{\alpha'}'''`). In staging, write the outer primes
  as `^{\prime\prime\prime}`; preserve the inner species prime and both
  script levels. Do not replace prose apostrophes or drop superscript structure.
  Repair known raw HTML serialization defects before make4ht DOM filters:
  misplaced `</mi>` after `</math>` and nonbreaking spaces in page-anchor tags.
  Do not accept HTML fallback parsing: it can corrupt surrounding inequalities
  and MathML. Read generated UTF-8 explicitly and resolve citation links against
  the actual unique destination anchor rather than trusting stale TeX4ht filenames.
- Full-guide builds take an explicit `--language` (English by default) and
  publish separate editions. The language dropdown reads the shared registry
  of successfully built editions and offers only editions containing the current
  page. Never advertise planned languages or publish incomplete translations.
  Check all inputs before conversion, retain the last successful edition on
  conversion/validation failure, and leave the original preview independent.
  Every non-English build requires a completed English HTML edition first;
  translate that edition without invoking TeX or regenerating figures. Use the
  active Python environment for both lxml and Argos. Translation builds must
  not download models or tokenizer metadata. Cache by source text, model bytes,
  and runtime version; checkpoint progress and preserve the last good published
  edition on failure. Validate protected placeholders before restoring markup.
- Plan multilingual editions as generated translations of the English source,
  with a language dropdown. Preserve math, symbols, code, FDS keywords, numbers,
  units, citations, and stable anchors. Cache translation outputs for repeatable
  builds; never require hand-maintained copies of translated manual chapters.
- When a mathematical typography defect is reported, audit every occurrence of
  that symbol/operator family across all converted chapters, including inline
  math, rather than patching only the reported equation.
- Specifically inspect every nabla after conversion: gradients (`\nabla T`,
  `\nabla Z_\alpha`), divergence (`\nabla\!\cdot` and `\Div`), curl,
  and powers such as `\nabla^2`. TeX4ht may emit nabla as a function-like `mi`
  followed by U+2061, introducing unwanted gaps. Normalize nabla and remove its
  spurious function-application marker; preserve actual function application
  for sin, cos, etc. Preserve explicit TeX spacing and script/operand structure.
- Accented variables must retain mathematical italics (for example, dotted m).
  Dot accents must remain legible above their base, including bold q and
  expressions with primes/subscripts. Apply the same spacing treatment to all
  single dot accents; do not move hats, bars, or unrelated accents.
  Both 0.12em and 0.03em mpadded lifts produced dots that were too high,
  including dotted m in equation 3.6. Use direct mover/mo accent attachment
  with the math font's native placement; do not reintroduce mpadded lifts.
  Judge spacing from rendered output, not merely passing structural tests.
- Every LaTeX `\overline{...}` must render as a drawn line spanning its full
  argument, whether a letter, product (such as rho Z), parenthesized group,
  fraction, or nested overline. Recognize TeX4ht's `mml-overline` marker and
  use a MathML base box with a CSS top border; do not substitute another accent
  glyph or assume stretchy attributes establish correct rendering. Preserve
  the full argument, subscripts/superscripts, anchors, and nested lines.
  Keep `\bar{...}` distinct as a short accent. This applies to all quantities,
  including pressure when its source specifically uses `\overline{p}`.
  Audit inline, display, and superscript occurrences, especially 3.4–3.6.
  The `\bW` macro already expands to `\overline{W}` and uses `mixture-overline`.
  The user preferred the drawn line, then requested a rightward adjustment:
  current alignment shifts the line 0.10em right and offsets its inner W back
  by the same amount, preserving the letter's position and equation spacing.
- Keep regression coverage for equations 2.6 (divergence and dotted m) and 2.9
  (gradients and dotted q), plus powered nabla and real named functions.
  Rebuild the reader and audit all generated occurrences after changes.
- Verify rendering visually when browser access is available. If it is blocked,
  report that limit and use the user's screenshots for joint review; do not
  claim that structural checks establish visual correctness.
