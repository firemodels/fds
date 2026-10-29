# Online FDS manuals from the existing LaTeX

Planning draft, 2026-10-28. Start with the Technical Reference Guide, Volume 1: Mathematical Model. This document proposes work; no converter or deployment has been implemented as part of this investigation.

## Objective and recommendation

Maintain technical content only in the existing LaTeX, bibliography, and figure sources. Generate both PDF and a searchable, navigable HTML manual automatically. Commit conversion code, configuration, and tests; never require manual edits to generated chapters, equations, captions, or references.

Target a **Sphinx-style documentation experience**. The user reviewed the FHE textbook example below and rejected its layout/navigation as a model for FDS. The proposed interface has a persistent chapter sidebar, search, a local section outline, breadcrumbs, previous/next navigation, and readable chapter pages. Evaluate a maintained Sphinx theme such as PyData or Read the Docs using an actual Tech Guide chapter.

Professional presentation is an explicit requirement. Start with a maintained theme and restrained FDS branding rather than a bespoke visual framework. Use consistent typography, comfortable reading width and spacing, crisp graphics, clear caption hierarchy, accessible contrast, and well-formatted equations and tables. Show the manual title/version and a PDF download clearly. The pilot should be polished enough to judge the intended finished product, not merely prove that HTML can be generated.

Keep the converter decision open. Compare Pandoc/Sphinx with TeX4ht through make4ht on a representative sample; consider LaTeXML if either needs extensive workarounds. make4ht remains a conversion candidate, not the proposed user interface. Use Python to collect inputs, run tools, validate output, and package the site; do not build a general LaTeX interpreter in Python.

| Candidate | Reason to consider it | What the pilot must establish |
| --- | --- | --- |
| make4ht / TeX4ht | Uses TeX processing; supports HTML, math, and graphics conversion, including TikZ-to-SVG. Likely a good fit for the existing manual. | Compatibility with FDS packages, macros, biblatex/Biber, diagrams, and equation numbering; acceptable HTML structure. |
| Pandoc + Sphinx | Existing starting point; good documentation navigation, search, and cross-reference infrastructure. | No silently dropped TeX; manageable translation of labels, multiline equations, tables, and figures into Sphinx nodes. |
| LaTeXML | Produces structured HTML/MathML and processes cross-references and graphics. | Package and macro coverage, bibliography fidelity, and amount of custom binding work. |

Require both content fidelity and the intended reading experience before choosing the stack. Prefer native Sphinx documents if conversion preserves the scientific content reliably. If a TeX-aware converter produces better results, measure the work needed to integrate its structured output into the documentation interface. Simply embedding raw HTML in Sphinx is insufficient: search, section navigation, and cross-references must understand the generated content. A separate HTML interface is an option only if its navigation and maintenance burden are demonstrated in the pilot.

## Findings from the current checkout

Inspected FDS revision: `b2bd16808fa3c76466b4a5eb83335575b611abf1`. Counts below are a static source inventory, not a successful conversion report.

- The 2023 prototype runs Pandoc on the User Guide and applies two `sed` substitutions for figure paths. Its Sphinx index points to a generated `test.rst`; there is no online-manual workflow in the current `.github/workflows` directory.
- The Tech Guide entry point includes 13 main chapters and `Appendices.tex`. With the main file, that is about 8,600 lines. The current common commands enable the Complex Geometry chapter with `\compgeomtrue`.
- The included sources contain 28 figure environments, including 9 `tikzpicture` and 4 `picture` environments, and 23 table environments. Several hundred displayed equations use custom `\be` / `\ee` macros; other forms include `align`, `eqnarray`, `gather`, and `multline`.
- The guide uses shared macros in `Bibliography/commoncommands.tex` and three bibliography databases with biblatex/Biber. Labels and references must survive across chapters and appendices.
- Five image inclusions in Complex Geometry refer to the separate `firemodels/fig` repository. They use explicit cropping options. Local inspection confirmed the referenced PNG files exist.
- `bot/Firebot/fds_guides2GH.sh` creates each figure archive from the contents of that manual's `SCRIPT_FIGURES` directory. Thus `FDS_TG_figures.tar.gz` supplies generated figures, not every figure in the guide. The guide currently references `particle_size_distribution` from that directory.
- The uploader replaces assets with `--clobber` and uploads revision information after the guides and archives. The release is updated through multiple operations, not one atomic snapshot.
- Live release metadata and archive contents were not retrieved in this investigation. Archive layout conclusions come from the local producer script; the first pilot must verify an actual downloaded bundle.

## Figure handling

| Source | Acquisition | Web representation |
| --- | --- | --- |
| Manual `FIGURES` and shared bibliography graphics | Pinned FDS checkout | Copy PNG/JPEG; convert PDF to SVG where suitable, with explicit high-resolution PNG fallback. |
| Manual `SCRIPT_FIGURES` | `FDS_TG_figures.tar.gz` from `firemodels/test_bundles`, release `FDS_TEST` | Extract into the staged manual's `SCRIPT_FIGURES`, then convert as above. |
| External `../../../fig/...` paths | Pinned `firemodels/fig` checkout, retaining the expected sibling layout | Copy/convert the referenced assets and preserve `trim`, `clip`, dimensions, and panel layout. |
| Inline TikZ, PGFPlots, and `picture` diagrams | The original LaTeX environment and its surrounding configuration | TeX-generated SVG; if necessary, automatically compile the figure body separately. Keep captions and labels as HTML. |

Resolve extensionless graphics with documented TeX-compatible precedence; fail if a referenced asset is missing. Preserve source-relative paths or use collision-resistant names. Cache conversions using source bytes, crop/scale options, converter configuration, and tool versions. Verify PDF page choice and font rendering. Do not assume an HTML image tag can directly display a PDF figure.

If a complex panel layout needs rendering as a single image, render its body automatically while retaining a selectable caption and a stable figure anchor. Keep original figure files available for download when useful. Derive initial accessible descriptions from captions; improve descriptions in the maintained source if needed.

## Proposed implementation sequence

1. **Build a representative conversion and navigation pilot.** Present the Equation chapter in a Sphinx-style interface, including persistent chapter navigation, search, and a local section outline. Use the original shared preamble and additional source-selected examples: Mass for `picture` and tables; Particle for the release figure; Complex Geometry for external images and cropping; Radiation/Appendices for multiline math, PGFPlots, tables, and appendix numbering. Generate sample inputs from the maintained sources rather than keeping another edited manual. Compare make4ht with a pinned Pandoc/Sphinx baseline, including navigation/search integration effort. Test LaTeXML if neither is satisfactory. Deliver browser-viewable samples and a short defect/maintenance comparison. A bare converter-output page does not satisfy this milestone.
2. **Convert the complete Tech Guide locally.** Follow the actual entry point, includes, conditionals, and bibliography order. Exclude unrelated notes such as `samr_notes.tex` and unused chapter files. Preserve preface, acknowledgments, disclaimer, scientific content, and references; replace print-only page layout with web presentation. Prefer external converter configuration for compatibility adaptations. Handle the hard-coded `pdftex` package drivers explicitly. Build a PDF from the same staged inputs for comparison.
3. **Provide one reproducible build command.** Implement a Python driver under this directory with separate acquisition, conversion, validation, and packaging stages. A proposed interface is `python build_online.py --manual tech --lock inputs.lock.json --output build/tech`. This command does not exist yet. Support an offline build from previously acquired inputs. Pin the TeX distribution/container, Python dependencies, converter, image tools, fonts, and browser math assets.
4. **Add publication through GitHub Actions.** Start with manual dispatch and preview artifacts. After the full-guide acceptance checks pass, add a daily scheduled build on the deployment repository's default branch. Download inputs, generate and validate the site, then publish the successful artifact through GitHub Pages. Keep the previous published site if acquisition, conversion, or validation fails. Record revisions visibly and link the matching PDF. PR checks can build previews without deployment credentials. A future notification from a completed Firebot upload can replace or supplement the schedule.
5. **Extend the same pipeline to other manuals.** Add per-manual entry points, bibliography dependencies, and release-asset names. Reuse the converter and validation rules. Test the User Guide next, then Verification and Validation, where figure volume and generated tables will need separate performance and dependency checks.

## Reproducibility and release consistency

Separate discovering the latest inputs from building an already selected set of inputs. The discovery stage writes a lock manifest containing the FDS commit, `fig` commit, source options, release/asset identifiers, archive SHA-256, revision-info contents, and toolchain/configuration identities. Cache or retain the actual archive bytes: a checksum alone cannot recover an overwritten release asset.

For the initial preview, record the exact source and figure versions and expose any unverified compatibility. An unchanged asset listing before and after download helps detect concurrent replacement but does not prove that all assets belong to one Firebot run. `FDS_INFO` is useful evidence, not an atomic manifest.

For dependable nightly publication, propose a producer-side manifest uploaded last, containing a generation ID, full source revisions, and hashes for each guide/archive. The consumer verifies every hash against that manifest and retries or skips incomplete generations. Prefer immutable, versioned bundle names or snapshots for retained builds. Until provenance is verified, do not describe a mixed source/bundle build as matching the released PDF.

Extract archives into a fresh staging directory, reject absolute/traversing paths and unsafe links, and avoid relying on files left by earlier builds. Normalize build dates using recorded source time, sort inventories, and keep run timestamps outside reproducible content. Rebuilding with the same lock and toolchain should produce identical normalized HTML/SVG/assets; test that property explicitly.

## Multilingual manuals and language selection

Add a language dropdown to the reader, initially offering English, Español,
Français, Deutsch, and Italiano. English LaTeX remains the maintained source;
other languages are generated translations, not separately authored manuals.

Prefer translation during the build, followed by static publication. The
dropdown selects a generated edition under `/en/`, `/es/`, `/fr/`, `/de/`, or
`/it/`. This keeps the reader fast, permits language-specific search, and avoids
requiring a translator API call or credentials in the reader's browser. The
translation provider remains a decision for a small pilot; compare scientific
terminology, glossary support, quality, cost, and API capabilities before
selecting a service. Do not provision or call a paid service as part of planning.

Proposed sequence:

1. Finish the English conversion and give paragraphs, sections, captions, and
   other translatable blocks stable identities derived from source locations or
   LaTeX labels. Keep URLs and anchors independent of translated heading text.
2. Extract prose, headings, navigation labels, captions, accessible descriptions,
   and textual table cells from the generated document. Translate complete
   contextual blocks, using protected placeholders for inline mathematics and
   other immutable content. Restore placeholders deterministically and reject
   responses that omit, duplicate, reorder semantically ordered protected items,
   or alter them. Do not send raw MathML or complete SVG markup for translation.
3. Preserve equations, variable names, mathematical symbols, code, FDS keywords,
   numeric values, units, reference numbers, citation keys, URLs, and label IDs.
   Keep bibliographic titles as published. Initially retain text embedded inside
   figure images/paths in English while translating captions and descriptions;
   translating diagram labels needs a separate source-aware process.
4. Maintain a small terminology glossary for FDS-specific language (for example,
   lumped species, heat release rate, and large-eddy simulation). Cache approved
   translation results by source content, context, target language, provider/model
   configuration, and glossary version. Reuse unchanged translations and request
   only changed or invalidated blocks. Human corrections belong in a translation
   memory keyed to their exact source, rather than manual edits to generated pages.
5. Generate a complete edition with localized navigation, its own search index,
   correct HTML language metadata, and shared figures/math assets. Switching
   languages should retain the current chapter and section. Remember an explicit
   language choice locally and keep English available. Clearly identify automatic
   translations and keep the PDF link labeled as an English PDF when applicable.
6. Pilot Spanish on the Overview chapter and obtain a bilingual technical review
   before expanding to French, German, and Italian. Check desktop/mobile layout,
   longer translated headings, scientific terminology, protected-token identity,
   search, and section-preserving language switching.

A fresh translator response is an external input, not guaranteed deterministic
even with fixed settings. Store its exact output and provenance so rebuilding
from the same source and translation snapshot is deterministic. Include the
translation snapshot hash in the build manifest. Keep credentials in the CI
secret store and set a translation budget before enabling scheduled requests.
If translation fails, retain the last complete published edition and clearly
identify its source revision; do not silently mix current English with stale
translated paragraphs. The dropdown should offer only editions that are ready.

## Acceptance criteria before publication

- All active chapters, appendices, figures, tables, captions, and cited bibliography entries are accounted for. No unsupported construct disappears silently.
- Internal links and label targets resolve. Equation, figure, table, chapter, and appendix numbers agree with the PDF built from the same inputs, including separately numbered lines and unnumbered equations. Do not let MathJax renumber independently without verifying the result.
- Custom math macros render correctly. Compare TeX-to-MathML and TeX-to-MathJax modes on the pilot; preserve readable, selectable mathematical content rather than rasterizing all equations.
- Every figure loads; cropping, scaling, panels, embedded fonts, and inline diagrams match the intended PDF content. Inspect all 28 current figures at desktop and narrow widths.
- Tables remain legible, captions can contain math, persistent chapter navigation and the local section outline work, and search finds representative prose and symbols where supported. Search results link to the relevant chapter/section. The mobile layout has usable navigation. Long equations and tables may scroll horizontally.
- The pilot demonstrates a professional visual design with coherent FDS branding, typography, spacing, captions, and navigation. Validate at desktop and mobile widths; the user can assess the actual styled chapter before full-guide rollout.
- A browser check waits for math rendering and detects missing images, math errors, and broken anchors. A clean build reports missing assets, unresolved references/citations, and unexpected converter warnings as failures; any unavoidable exception is explicit and reviewed.
- Two clean builds of the same locked inputs agree. Generated documents require no manual fixes.

## References and examples

- [make4ht documentation](https://www.kodymirus.cz/make4ht/make4ht-doc.html): build stages, filters, image conversion, and output options.
- [TeX4ht configuration documentation](https://www.kodymirus.cz/tex4ht-doc/Configurations.html): external configuration, graphics, TikZ, math, and bibliographies.
- [LaTeXML postprocessing documentation](https://math.nist.gov/~BMiller/LaTeXML/manual/usage/post/): references, bibliography, math, and graphics processing.
- [Pandoc manual](https://pandoc.org/MANUAL.html): supported formats, macro expansion, and raw TeX handling.
- [The Beginner's Textbook for Fully Homomorphic Encryption](https://fhetextbook.github.io/): evidence of make4ht use only; its [preface](https://fhetextbook.github.io/Preface.html) credits make4ht. The user rejected its layout/navigation as the desired FDS experience.
- [Learning Deep Representations of Data Distributions](https://ma-lab-berkeley.github.io/deep-representation-learning-book/): another mathematical book with an online edition. Its [repository build instructions](https://github.com/Ma-Lab-Berkeley/deep-representation-learning-book#building-the-website) describe make4ht, generated MathJax macros, and HTML postprocessing.

The first decision point is a pilot demonstrating both conversion fidelity and the proposed documentation interface, before full-guide rollout or scheduled deployment. Hosting location and the daily schedule can be chosen after that pilot; neither blocks evaluating the pilot now.
