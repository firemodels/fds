# FDS Repository Guidance

Follow the [Developer Commit Guidelines](https://github.com/firemodels/fds/wiki/Developer-Commit-Guidelines). Key requirements:

- Fortran: use `IMPLICIT NONE`, uppercase code, three-space indentation without tabs, and lines no longer than 132 characters. Follow existing spacing and `IF (...) THEN`, `ENDIF`, and `ENDDO` conventions.
- Use descriptive names, named constants, explicit real precision (`_EB`/`_FB`), and tolerances for real comparisons. Name long blocks and qualify their `CYCLE`/`EXIT` statements; alphabetize namelist inputs.
- Add verification tests for code changes and document new cases and inputs in the guides. Build affected source with repository makefiles in debug mode, without errors or warnings, using GNU and Intel when available.
- Update relevant manuals and release notes; mark beta features. Soft-wrap prose, follow NIST SP 811, and compile edited guides with `make_guide.sh`. Keep generated binaries and non-release PDFs out of commits.
- Keep changes and PRs small and focused, with no merge commits. Use descriptive commit prefixes such as `FDS Source`, `FDS User Guide`, or `Python`.

## Python dependency compatibility across Firemodels

- Keep Python dependencies compatible across FDS, CFAST, and SMV on Linux, macOS, and Windows; scope platform- or Python-specific pins with environment markers.
- Validate dependency changes against the supported Python versions and affected CI/bundle builds; explicitly report any untested configurations.
